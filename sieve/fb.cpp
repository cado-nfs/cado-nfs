#include "cado.h"
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cctype>
#include <gmp.h>
#include <pthread.h>
#ifdef HAVE_SYS_MMAN_H
#include <sys/mman.h>
#else
#define MAP_FAILED ((void *) -1)
#endif
#include "fb.hpp"
#include "mod_ul.h"
#include "verbose.h"
#include "getprime.h"
#include "gmp_aux.h"
#include "gzip.h"
#include "threadpool.hpp"


static unsigned int fb_log_2 (fbprime_t);
static bool fb_linear_root (fbroot_t *, const mpz_t *, fbprime_t);

/* strtoul(), but with const char ** for second argument.
   Otherwise it's not possible to do, e.g., strtoul(p, &p, 10) when p is
   of type const char *
*/
static inline unsigned long int
strtoul_const(const char *nptr, const char **endptr, const int base)
{
  char *end;
  unsigned long r;
  r = strtoul(nptr, &end, base);
  *endptr = end;
  return r;
}

static inline unsigned long long int
strtoull_const(const char *nptr, const char **endptr, const int base)
{
  char *end;
  unsigned long long r;
  r = strtoull(nptr, &end, base);
  *endptr = end;
  return r;
}

// Adapted from utils/ularith.h
// TODO: this function should go somewhere else...
static inline uint64_t
uint64_invmod(const uint64_t n)
{
    uint64_t r;
    ASSERT (n % UINT64_C(2) != UINT64_C(0));
    r = (UINT64_C(3) * n) ^ UINT64_C(2);
    r = UINT64_C(2) * r - (uint32_t) r * (uint32_t) r * (uint32_t) n;
    r = UINT64_C(2) * r - (uint32_t) r * (uint32_t) r * (uint32_t) n;
    r = UINT64_C(2) * r - (uint32_t) r * (uint32_t) r * (uint32_t) n;
    uint32_t k = (uint32_t)(r * n >> 32);
    k *= (uint32_t) r;
    r = r - ((uint64_t)k << 32);
    return r;
}

static inline redc_invp_t
compute_invq(fbprime_t q)
{
  if (q % 2 != 0) {
    if (sizeof(unsigned long) >= sizeof(redc_invp_t)) {
        return (redc_invp_t) (- ularith_invmod (q));
    } else {
        ASSERT(sizeof(redc_invp_t) == 8);
        return (redc_invp_t) (- uint64_invmod (q));
    }
  } else {
    return 0;
  }
}

/* If q = p^k, with k maximal and k > 1, return q.
   Otherwise return 0. If final_k is not NULL, write k there. */
fbprime_t
fb_is_power (fbprime_t q, unsigned long *final_k)
{
  unsigned long maxk, k;
  uint32_t p;

  maxk = fb_log_2(q);
  for (k = maxk; k >= 2; k--)
    {
      double dp = pow ((double) q, 1.0 / (double) k);
      double rdp = trunc(dp + 0.5);
      if (fabs(dp - rdp) < 0.001) {
        p = (uint32_t) rdp ;
        if (q % p == 0) {
          // ASSERT (fb_pow (p, k) == q);
          if (final_k != NULL)
            *final_k = k;
          return p;
        }
      }
    }
  return 0;
}


/* Allow construction of a root from a linear polynomial and a prime (power) */
fb_general_root::fb_general_root (fbprime_t q, const mpz_t *poly,
  const unsigned char nexp, const unsigned char oldexp)
  : exp(nexp), oldexp(oldexp)
{
  proj = fb_linear_root (&r, poly, q);
}


/* Allow assignment-construction of general entries from simple entries */
template <int Nr_roots>
fb_general_entry::fb_general_entry (const fb_entry_x_roots<Nr_roots> &e) {
  p = q = e.p;
  k = 1;
  invq = e.invq;
  for (int i = 0; i != Nr_roots; i++) {
    /* Use simple constructor for root */
    roots[i] = e.roots[i];
  }
}


/* Return whether this is a simple factor base prime.
   It is simple if it is a prime (not a power) and all its roots are simple. */
bool
fb_general_entry::is_simple() const
{
  bool is_simple = (k == 1);
  for (unsigned char i = 0; i != nr_roots; i++) {
    is_simple &= roots[i].is_simple();
  }
  return is_simple;
}


/* Read roots from a factor base file line and store them in roots.
   line must point at the first character of the first root on the line.
   linenr is used only for printing error messages in case of parsing error.
   Returns the number of roots read. */
void
fb_general_entry::read_roots (const char *lineptr, const unsigned char nexp,
                              const unsigned char oldexp,
                              const unsigned long linenr)
{
    unsigned long long last_t = 0;

    nr_roots = 0;
    while (*lineptr != '\0')
    {
        if (nr_roots == MAX_DEGREE) {
            verbose_output_print (1, 0,
                    "# Error, too many roots for prime (power) %" FBPRIME_FORMAT
                    " in factor base line %lu\n", q, linenr);
            exit(EXIT_FAILURE);
        }
        /* Projective roots r, i.e., ar == b (mod q), are stored as q + r in
           the factor base file; since q can be a 32-bit value, we read the
           root as a 64-bit integer first and subtract q if necessary. */
        const unsigned long long t = strtoull_const (lineptr, &lineptr, 10);
        if (nr_roots > 0 && t <= last_t) {
            verbose_output_print (1, 0,
                "# Error, roots must be sorted in the fb file, line %lu\n",
                linenr);
            exit(EXIT_FAILURE);
        }
        last_t = t;

        roots[nr_roots++] = fb_general_root(t, q, nexp, oldexp);
        if (*lineptr != '\0' && *lineptr != ',') {
            verbose_output_print (1, 0,
                    "# Incorrect format in factor base file line %lu\n",
                    linenr);
            exit(EXIT_FAILURE);
        }
        if (*lineptr == ',')
            lineptr++;
    }

    if (nr_roots == 0) {
        verbose_output_print (1, 0, "# Error, no root for prime (power) %"
                FBPRIME_FORMAT " in factor base line %lu\n", q, linenr - 1);
        exit(EXIT_FAILURE);
    }
}

/* Parse a factor base line.
   Return 1 if the line could be parsed and was a "short version", i.e.,
   without explicit old and new exponent.
   Return 2 if the line could be parsed and was a "long version".
   Otherwise return 0. */
void
fb_general_entry::parse_line (const char * lineptr, const unsigned long linenr)
{
    q = strtoul_const (lineptr, &lineptr, 10);
    if (q == 0) {
        verbose_output_print (1, 0, "# fb_read: prime is not an integer on line %lu\n",
                              linenr);
        exit (EXIT_FAILURE);
    } else if (*lineptr != ':') {
        verbose_output_print (1, 0,
                "# fb_read: prime is not followed by colon on line %lu",
                linenr);
        exit (EXIT_FAILURE);
    }

    lineptr++; /* Skip colon after q */
    const bool longversion = (strchr(lineptr, ':') != NULL);

    /* NB: a short version is not permitted for a prime power, so we
     * do the test for prime powers only for long version */
    p = q;
    unsigned char nexp = 1, oldexp = 0;
    k = 1;
    if (longversion) {
        unsigned long k_ul;
        const fbprime_t base = fb_is_power (q, &k_ul);
        ASSERT(ulong_isprime(base != 0 ? base : q));
        /* If q is not a power, then base==0, and we use p = q */
        if (base != 0) {
            p = base;
            k = static_cast<unsigned char>(k_ul);
        }

        /* read the multiple of logp, if any */
        /* this must be of the form  q:nlogp,oldlogp: ... */
        /* if the information is not present, it means q:1,0: ... */
        nexp = strtoul_const (lineptr, &lineptr, 10);

        if (nexp == 0) {
            verbose_output_print (1, 0, "# Error in fb_read: could not parse "
                "the integer after the colon of prime %" FBPRIME_FORMAT "\n",
                q);
            exit (EXIT_FAILURE);
        }
        if (*lineptr != ',') {
            verbose_output_print (1, 0,
		    "# fb_read: exp is not followed by comma on line %lu",
		    linenr);
            exit (EXIT_FAILURE);
        }
        lineptr++; /* skip comma */
        oldexp = strtoul_const (lineptr, &lineptr, 10);
        if (*lineptr != ':') {
            verbose_output_print (1, 0,
		    "# fb_read: oldlogp is not followed by colon on line %lu",
		    linenr);
            exit (EXIT_FAILURE);
        }
        ASSERT (nexp > oldexp);
        lineptr++; /* skip colon */
    }

    read_roots(lineptr, nexp, oldexp, linenr);

    /* exp and oldexp are a property of a root, not of a prime (power).
       The factor base file should specify them per root, but specifies
       them per prime instead - a bit of a design bug.
       For long version lines, we thus use the exp and oldexp values for all
       roots specified in that line. */
    for (unsigned char i = 1; i < nr_roots; i++) {
      roots[i].exp = roots[0].exp;
      roots[i].oldexp = roots[0].oldexp;
    }
}

void
fb_general_entry::fprint(FILE *out) const
{
  fprintf(out, "%" FBPRIME_FORMAT ": ", q);
  for (unsigned char i_root = 0; i_root != nr_roots; i_root++) {
    roots[i_root].fprint(out, q);
    if (i_root + 1 < nr_roots)
      fprintf(out, ",");
  }
  fprintf(out, "\n");
}

void
fb_general_entry::merge (const fb_general_entry &other)
{
  ASSERT_ALWAYS(p == other.p && q == other.q && k == other.k);
  for (unsigned char i_root = 0; i_root < other.nr_roots; i_root++) {
    ASSERT_ALWAYS(nr_roots < MAX_DEGREE);
    roots[nr_roots++] = other.roots[i_root];
  }
}

void
fb_general_entry::transform_roots(fb_general_entry::transformed_entry_t &result,
                                  const qlattice_basis &basis) const
{
  result.p = p;
  result.q = q;
  result.invq = invq;
  result.k = k;
  result.nr_roots = nr_roots;
  /* TODO: Use batch-inversion here */
  for (unsigned char i_root = 0; i_root != nr_roots; i_root++)
    roots[i_root].transform(result.roots[i_root], q, invq, basis);
}


template <int Nr_roots>
void
fb_entry_x_roots<Nr_roots>::transform_roots(fb_entry_x_roots<Nr_roots>::transformed_entry_t &result, const qlattice_basis &basis) const
{
  result.p = p;
  /* TODO: Use batch-inversion here */
  for (unsigned char i_root = 0; i_root != nr_roots; i_root++) {
    const unsigned long long t = fb_root_in_qlattice(p, roots[i_root], invq, basis);
    result.proj[i_root] = (t >= p);
    result.roots[i_root] = (t < p) ? t : (t - p);
  }
}

// FIXME: why do I have to make those instances explicit???
// If someone knows how to avoid that...

template void
fb_entry_x_roots<0>::transform_roots(fb_transformed_entry_x_roots<0> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<1>::transform_roots(fb_transformed_entry_x_roots<1> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<2>::transform_roots(fb_transformed_entry_x_roots<2> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<3>::transform_roots(fb_transformed_entry_x_roots<3> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<4>::transform_roots(fb_transformed_entry_x_roots<4> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<5>::transform_roots(fb_transformed_entry_x_roots<5> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<6>::transform_roots(fb_transformed_entry_x_roots<6> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<7>::transform_roots(fb_transformed_entry_x_roots<7> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<8>::transform_roots(fb_transformed_entry_x_roots<8> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<9>::transform_roots(fb_transformed_entry_x_roots<9> &, qlattice_basis const&) const; 
template void 
fb_entry_x_roots<10>::transform_roots(fb_transformed_entry_x_roots<10> &, qlattice_basis const&) const; 



template <int Nr_roots>
void
fb_entry_x_roots<Nr_roots>::fprint(FILE *out) const
{
  fprintf(out, "%" FBPRIME_FORMAT ": ", p);
  for (int i = 0; i != Nr_roots; i++) {
    fprintf(out, "%" FBROOT_FORMAT "%s", roots[i],
	    (i + 1 < Nr_roots) ? "," : "");
  }
  fprintf(out, "\n");
}

/* workaround for compiler bug. gcc-4.2.1 on openbsd 5.3 does not seem to
 * accept the access to foo.nr_roots when nr_roots is in fact a static
 * const member. That's a pity, really. I think the code below is
 * innocuous performance-wise. */
template <class FB_ENTRY_TYPE>
struct get_nroots {
    FB_ENTRY_TYPE const & r;
    get_nroots(FB_ENTRY_TYPE const & r) : r(r) {}
    inline unsigned char operator()() const { return r.nr_roots; }
};

template <int n>
struct get_nroots<fb_entry_x_roots<n> > {
    get_nroots(fb_entry_x_roots<n> const &) {}
    inline unsigned char operator()() const { return fb_entry_x_roots<n>::nr_roots; }
};


template <class FB_ENTRY_TYPE>
fb_vector<FB_ENTRY_TYPE>::~fb_vector()
{
  /* Free all the slices */
  typename std::map<const double, const slices_t *>::const_iterator it;
  for (it = cached_slices.begin(); it != cached_slices.end(); it++) {
    delete it->second;
  }
  clear();
}


template <class FB_ENTRY_TYPE>
void
fb_vector<FB_ENTRY_TYPE>::extract_bycost(std::vector<unsigned long> &p,
    fbprime_t pmax, fbprime_t td_thresh) const
{
  for (size_t i = 0; i < size; i++)
    data[i].extract_bycost(p, pmax, td_thresh);
}


template <class FB_ENTRY_TYPE>
plattices_vector_t *
fb_slice<FB_ENTRY_TYPE>::make_lattice_bases(const qlattice_basis &basis,
    const int logI, const sublat_t &sublat) const
{
  typename FB_ENTRY_TYPE::transformed_entry_t transformed;
  /* Create a transformed vector and store the index of the slice we currently
     transform */
  const unsigned long special_q = mpz_fits_ulong_p(basis.q) ? mpz_get_ui(basis.q) : 0;

  plattices_vector_t *result = new plattices_vector_t(get_index());
  slice_offset_t i_entry = 0;
  for (const FB_ENTRY_TYPE *it = begin(); it != end(); it++, i_entry++) {
    if (it->p == special_q) /* Assumes it->p != 0 */
      continue;
    it->transform_roots(transformed, basis);
    for (unsigned char i_root = 0; i_root != transformed.nr_roots; i_root++) {
      const fbroot_t r = transformed.get_r(i_root);
      const bool proj = transformed.get_proj(i_root);
      /* If proj and r > 0, then r == 1/p (mod p^2), so all hits would be in
         locations with p | gcd(i,j). */
      if (LIKELY(!proj || r == 0)) {
        plattice_info_t pli = plattice_info_t(transformed.get_q(), r, proj, logI);
        plattice_enumerate_t ple = plattice_enumerate_t(pli, i_entry, logI, sublat);
        // Skip (0,0) unless we have sublattices.
        if (!sublat.m)
          ple.next();
        if (LIKELY(pli.a0 != 0)) {
          result->push_back(ple);
        }
      }
    }
  }
  return result;
}


template <class FB_ENTRY_TYPE>
void
fb_slice<FB_ENTRY_TYPE>::fprint(FILE *out) const
{
  for (const FB_ENTRY_TYPE *it = begin(); it != end(); it++) {
    it->fprint(out);
  }
}


template <class FB_ENTRY_TYPE>
void
fb_vector<FB_ENTRY_TYPE>::_count_entries(size_t *nprimes, size_t *nroots, double *weight) const
{
  if (nprimes != NULL)
    *nprimes += size;
  double w = 0.;
  size_t nr = 0;
  for (size_t i = 0; i < size; i++) {
#if 0
      /* see compiler bug section above... */
      unsigned char nr_i = data[i].nr_roots;
#else
      /* use workaround */
      unsigned char nr_i = get_nroots<FB_ENTRY_TYPE>(data[i])();
#endif
      nr += nr_i;
      w += data[i].weight();
  }
  if (nroots != NULL)
    *nroots += nr;
  if (weight != NULL)
    *weight += w;
}

/* Compute an upper bound on this slice's weight by assuming all primes in the
   slice are equal on the first (and thus smallest) one. This relies on the
   vector being sorted. */
template <class FB_ENTRY_TYPE>
double
fb_vector<FB_ENTRY_TYPE>::est_weight_max(const size_t start, const size_t end) const
{
  return data[start].weight() * (end - start);
}

/* Estimate weight by the average of the weights of the two endpoints, i.e.,
   by trapezoidal rule. */
template <class FB_ENTRY_TYPE>
double
fb_vector<FB_ENTRY_TYPE>::est_weight_avg(const size_t start, const size_t end) const
{
  return (data[start].weight() + data[end-1].weight()) / 2. * (end - start);
}

/* Estimate weight by Simpson's rule on the weight of the two endpoints and
   of the midpoint */
template <class FB_ENTRY_TYPE>
double
fb_vector<FB_ENTRY_TYPE>::est_weight_simpson(const size_t start, const size_t end) const
{
  const size_t midpoint = start + (end - start) / 2;
  return (data[start].weight() + 4.*data[midpoint].weight() + data[end-1].weight()) / 6. * (end - start);
}

/* Estimate weight by using Merten's rule on the primes at the two endpoints */
template <class FB_ENTRY_TYPE>
double
fb_vector<FB_ENTRY_TYPE>::est_weight_mertens(const size_t start, const size_t end) const
{
  return log(log(data[end-1].get_q())) - log(log(data[start].get_q()));
}

/* Compute weight exactly with a sum over all entries */
template <class FB_ENTRY_TYPE>
double
fb_vector<FB_ENTRY_TYPE>::est_weight_sum(const size_t start, const size_t end) const
{
  double sum = 0.;
  for (size_t i = start; i < end; i++)
    sum += data[i].weight();
  return sum;
}

class wurst {
  double worst_err, worst_est, worst_val, mse;
  unsigned long nr;
public:
  wurst() : worst_err(0.), worst_est(0.), worst_val(0.), mse(0.), nr(0){}
  void update(const double est, const double val) {
    const double err = est / val - 1.;
    mse += err * err;
    nr++;
    if (fabs(err) > fabs(worst_err)) {
      worst_err = err;
      worst_est = est;
      worst_val = val;
    }
    verbose_output_print(0, 4, "%0.3g (%.3g)", val, err);
  }
  void print() {
    if (nr > 0)
      verbose_output_print(0, 4, "%0.3g vs. %0.3g (rel err. %.3g, MSE: %.3g)",
                           worst_est, worst_val, worst_err, mse / nr);
  }
};

static wurst worst_max, worst_avg, worst_simpson, worst_mertens;

void print_worst_weight_errors()
{
    verbose_output_start_batch();
    verbose_output_print(0, 4, "# Worst weight errors: max = ");
    worst_max.print();
    verbose_output_print(0, 4, ", avg = ");
    worst_avg.print();
    verbose_output_print(0, 4, ", simpson = ");
    worst_simpson.print();
    verbose_output_print(0, 4, ", mertens = ");
    worst_mertens.print();
    verbose_output_print(0, 4, "\n");
    verbose_output_end_batch();
}

template <class FB_ENTRY_TYPE>
double
fb_vector<FB_ENTRY_TYPE>::est_weight_compare(const size_t start, const size_t end) const
{
  const double _max = est_weight_max(start, end),
               avg = est_weight_avg(start, end),
               simpson = est_weight_simpson(start, end),
               mertens = est_weight_mertens(start, end),
               _sum = est_weight_sum(start, end);
  verbose_output_start_batch();
  verbose_output_print(0, 4, "# Slice [%zu, %zu] weight: max = ", start, end);
  worst_max.update(_max, _sum);
  verbose_output_print(0, 4, ", avg = ");
  worst_avg.update(avg, _sum);
  verbose_output_print(0, 4, ", simpson = ");
  worst_simpson.update(simpson, _sum);
  verbose_output_print(0, 4, ", mertens = ");
  worst_mertens.update(mertens, _sum);
  verbose_output_print(0, 4, ", sum = %0.3g\n", _sum);
  verbose_output_end_batch();
  return _sum;
}

template <class FB_ENTRY_TYPE>
double
fb_vector<FB_ENTRY_TYPE>::est_weight(const size_t start, const size_t end) const
{
  if (verbose_output_get(0, 4, 0) != NULL)
    return est_weight_compare(start, end);
  return est_weight_max(start, end);
}

/* For general vectors, we compute the weight the hard way, via a sum over
   all entries. General vectors are quite small so this should not take long */
template <>
double
fb_vector<fb_general_entry>::est_weight(const size_t start, const size_t end) const
{
  return est_weight_sum(start, end);
}

/* scale is the quantity by which we multiply the log_2 of the primes we
 * consider */
template <class FB_ENTRY_TYPE>
void
fb_vector<FB_ENTRY_TYPE>::make_slices(const double scale, const double max_weight,
                                      slice_index_t &next_index)
{
  /* Entries in data must be sorted in order of non-decreasing round(log(p));
     since the logarithm base varies, this implies order of non-decreasing p.
  */

  if (size == 0) {
    /* Nothing to do. Leave slices == NULL */
    ASSERT_ALWAYS(slices == NULL);
    return;
  }

  /* If we already have a slicing for this scale in the cache, use that one */
  if (!cached_slices.empty() &&
      cached_slices.find(scale) != cached_slices.end()) {
    slices = cached_slices[scale];
    return;
  }

  /* Otherwise create a new slicing */
  verbose_output_print (0, 3, "# Creating new slicing of ");
  if (FB_ENTRY_TYPE::is_general_type) {
    verbose_output_print (0, 3, "general vector");
  } else {
    verbose_output_print (0, 3, "vector with %d root%s",
                          FB_ENTRY_TYPE::fixed_nr_roots,
                          FB_ENTRY_TYPE::fixed_nr_roots != 1 ? "s" : "");
  }
  verbose_output_print (0, 3, " for scale=%1.2f\n", scale);

  slices_t *new_slices = new slices_t;

  size_t cur_slice_start = 0;

  while (cur_slice_start < size) {
    const unsigned char cur_logp = fb_log (data[cur_slice_start].p, scale, 0.);

    size_t next_slice_start = std::min(cur_slice_start + max_slice_len, size);
    ASSERT_ALWAYS(cur_slice_start < next_slice_start);
    ASSERT_ALWAYS(data[cur_slice_start].p <= data[next_slice_start - 1].p);

    /* See if last element of the current slice has a greater log(p) than
       the first element */
    if (cur_logp < fb_log (data[next_slice_start - 1].p, scale, 0.)) {
      /* Find out the first place where log(p) changes occurs, and make that
         next_slice_start */
      for (next_slice_start = cur_slice_start;
           cur_logp == fb_log (data[next_slice_start].p, scale, 0.);
           next_slice_start++);
      /* We know that a log(p) change occurs so the slice len cannot exceed
         the max allowed */
      ASSERT_ALWAYS(next_slice_start <= cur_slice_start + max_slice_len);
    }

    /* Maybe the slice's weight is greater than max_weight and we have to make
       the slice smaller */
    double weight;
    while ((weight = est_weight(cur_slice_start, next_slice_start)) > max_weight) {
      const size_t old = next_slice_start;
      next_slice_start = cur_slice_start + (next_slice_start - cur_slice_start) / 2;
      verbose_output_print (0, 3, "# Slice %u starting at offset %zu has too "
        "great weight %.3f > %.3f, shortening end from %zu to %zu\n",
        (unsigned int) next_index, cur_slice_start, weight, max_weight, old,
        next_slice_start);
      ASSERT_ALWAYS(next_slice_start > cur_slice_start);
    }

    verbose_output_print (0, 4, "# Slice %u starts at offset %zu (p = %"
        FBPRIME_FORMAT ", log(p) = %u) and ends at offset %zu (p = %"
        FBPRIME_FORMAT ", log(p) = %u), weight = %.3f\n",
           (unsigned int) next_index, cur_slice_start, data[cur_slice_start].p,
           (unsigned int) fb_log (data[cur_slice_start].p, scale, 0.), 
           next_slice_start - 1, data[next_slice_start - 1].p,
           (unsigned int) fb_log (data[next_slice_start - 1].p, scale, 0.),
           weight);
    fb_slice<FB_ENTRY_TYPE> s(data + cur_slice_start,
                              data + next_slice_start, cur_logp,
                              next_index++, weight);
    new_slices->push_back(s);

    cur_slice_start = next_slice_start;
  }

  /* Add new slicing to the cache */
  cached_slices[scale] = new_slices;
  /* Make the new slicing the currently used one */
  slices = new_slices;
}

template <class FB_ENTRY_TYPE>
void
fb_vector<FB_ENTRY_TYPE>::sort()
{
  std::sort(data, data+size);
}


template <class FB_ENTRY_TYPE>
void
fb_vector<FB_ENTRY_TYPE>::fprint(FILE *out) const
{
  /* If we have sliced the entries, we print them one slice at a time */
  if (!slices->empty()) {
    for (typename std::vector<fb_slice<FB_ENTRY_TYPE> >::const_iterator it = slices->begin(); it != slices->end(); it++) {
      fprintf(out, "#    Slice index = %u, logp = %hhu:\n", (unsigned int) it->get_index(), it->get_logp());
      it->fprint(out);
    }
  } else {
    /* Otherwise we print the whole vector */
    fprintf(out, "#    Not sliced\n");
    for (size_t i = 0; i < size; i++)
      data[i].fprint(out);
  }
}

// Return NULL on error.
template <class FB_ENTRY_TYPE>
void *
fb_vector<FB_ENTRY_TYPE>::mmap_fbc(void *p, void * const end)
{
  struct vec_t {
    size_t        size;
    FB_ENTRY_TYPE data[0];
  };
  vec_t *v = static_cast<vec_t *>(p);

  if (v->data > end || v->data + v->size > end) {
    verbose_output_print(1, 0, "# Could not read memory image of factor base: truncated file?\n");
    return NULL;
  }

  data      = v->data;
  size      = v->size;
  alloc     = 0;
  read_only = true;
  mmapped   = true;

  return static_cast<void *>(data + size);
}

template <class FB_ENTRY_TYPE>
bool
fb_vector<FB_ENTRY_TYPE>::dump_fbc(FILE *f) const
{
  ASSERT_ALWAYS(read_only);
  if (fwrite(&size, sizeof(size),          1,    f) != 1 ||
      fwrite( data, sizeof(FB_ENTRY_TYPE), size, f) != size) {
    verbose_output_print(1, 0, "# Could not write memory image of factor base to file\n");
    return false;
  }
  return true;
}

template <class FB_ENTRY_TYPE>
void
fb_vector<FB_ENTRY_TYPE>::clear()
{
  if (!mmapped)
    free(data);
  init();
}

template <class FB_ENTRY_TYPE>
void
fb_vector<FB_ENTRY_TYPE>::init()
{
  data      = NULL;
  size      = 0;
  alloc     = 0;
  read_only = false;
  mmapped   = false;
}


/* Append a factor base entry given in fb_cur to the
   correct vector, as determined by the number of roots.
   "Special" primes are added to the general-primes vector.
   This function is a de-multiplexer: Using "switch" to turn a run-time value
   (fb_cur->nr_roots) into a compile-time constant which can be used as a
   template specifier. */
void
fb_part::append(const fb_general_entry &fb_cur)
{
  /* Non-simple ones go in the general vector, or, if this is an only_general
     part, then all entries do */
  if (only_general || !fb_cur.is_simple()) {
    /* Append powers only up to this part's powlim */
    if (fb_cur.k == 1 || fb_cur.q < powlim)
      general_vector.append(fb_cur);
    return;
  }

  /* Simple ones go in the simple vector with the corresponding number of
     roots */
  get_slices(fb_cur.nr_roots)->append(fb_cur);
}

void
fb_part::finalize()
{
  if (!only_general) {
    for (int i_roots = 0; i_roots <= MAX_DEGREE; i_roots++) {
      get_slices(i_roots)->finalize();
    }
  }
  general_vector.finalize();
}

void
fb_part::fprint(FILE *out) const
{
  if (!only_general) {
    for (int i_roots = 0; i_roots <= MAX_DEGREE; i_roots++) {
      fprintf(out, "#   Entries with %d roots:\n", i_roots);
      get_slices(i_roots)->fprint(out);
    }
  }

  fprintf(out, "#   General entries (%s):\n",
	  only_general ? "contains all entries" :
	  "powers, ramified primes or primes with projective roots");
  general_vector.fprint(out);
}

void
fb_part::_count_entries(size_t *nprimes, size_t *nroots, double *weight) const
{
  if (!only_general) {
    for (int i_roots = 0; i_roots <= MAX_DEGREE; i_roots++)
      get_slices(i_roots)->_count_entries(nprimes, nroots, weight);
  }
  general_vector._count_entries(nprimes, nroots, weight);
}


void
fb_part::extract_bycost(std::vector<unsigned long> &p, fbprime_t pmax, fbprime_t td_thresh) const
{
  if (!only_general) {
    for (int i_roots = 0; i_roots <= MAX_DEGREE; i_roots++)
      get_slices(i_roots)->extract_bycost(p, pmax, td_thresh);
  }

  general_vector.extract_bycost(p, pmax, td_thresh);
}


/* scale is the quantity by which we multiply the log_2 of the primes we
 * consider */
void
fb_part::make_slices(const double scale, const double max_weight,
                     slice_index_t &next_index)
{
  if (!only_general) {
    for (int i_roots = 0; i_roots <= MAX_DEGREE; i_roots++)
      get_slices(i_roots)->make_slices(scale, max_weight, next_index);
    /* If we store all entries as general entries, we don't slice them,
       as those are the small primes in part 0 which get line sieved */
    /* FIXME: slicing the general vector is somewhat silly even if it gets
       bucket-sieved, as there are only few entries, scattered over the
       interval [2, powlim] (and perhaps a bigger prime or two with
       projective/ramified roots). Thus most slices contain only a few
       primes or prime powers. It would be better to let fill_in_buckets()
       write prime-hint updates instead when sieving the general vector,
       so that we don't have to write a pointer array for nearly-empty
       slices. */
    general_vector.make_slices(scale, max_weight, next_index);
  }
}

// Return NULL on error.
void *
fb_part::mmap_fbc(void *p, void * const end)
{
  if (!only_general) {
    for (int i_roots = 0; p != NULL && i_roots <= MAX_DEGREE; ++i_roots)
      p = get_slices(i_roots)->mmap_fbc(p, end);
  }
  if (p != NULL)
    p = general_vector.mmap_fbc(p, end);
  return p;
}

bool
fb_part::dump_fbc(FILE *f) const
{
  bool rc = true;
  if (!only_general) {
    for (int i_roots = 0; rc && i_roots <= MAX_DEGREE; ++i_roots)
      rc = get_slices(i_roots)->dump_fbc(f);
  }
  if (rc)
    rc = general_vector.dump_fbc(f);
  return rc;
}

void
fb_part::clear()
{
  if (!only_general) {
    for (int i_roots = 0; i_roots <= MAX_DEGREE; ++i_roots)
      get_slices(i_roots)->clear();
  }
  general_vector.clear();
}

/* powlim could well be one value per part, like thresholds */
fb_factorbase::fb_factorbase(const fbprime_t *thresholds,
                             const fbprime_t powlim,
			     const bool *only_general)
{
  for (size_t i = 0; i < FB_MAX_PARTS; i++) {
    ASSERT_ALWAYS(i == 0 || thresholds[i - 1] <= thresholds[i]);
    this->thresholds[i] = thresholds[i];
    // By default, only_general is true for part 0, and false for all others
    const bool og = (only_general == NULL) ? (i == 0) : only_general[i];
    parts[i] = new fb_part(powlim, og);
  }
  toplevel = 0;
}

fb_factorbase::~fb_factorbase()
{
  for (size_t i = 0; i < FB_MAX_PARTS; i++) {
    delete parts[i];
  }
}

/* Append a factor base entry to the factor base.
   The new entry is inserted into the correct part, as determined by the
   size of the prime p, and within that part, into the correct slice, as
   determined by the number of roots. */
void
fb_factorbase::append(const fb_general_entry &fb_cur)
{
  int i;
  static bool printed_too_large_prime_warning = false;

  /* Find the smallest threshold t such that t >= q */
  for (i = 0; i < FB_MAX_PARTS && fb_cur.q > thresholds[i]; i++);
  /* No prime > largest threshold should ever be added */
  if (i == FB_MAX_PARTS) {
    if (!printed_too_large_prime_warning) {
      verbose_output_print(1, 0, "# Factor base entry %" FBPRIME_FORMAT 
                           " is above factor base bound, skipping it "
                           "(and all other too large entries)\n", fb_cur.q);
      printed_too_large_prime_warning = true;
    }
    return; /* silently skip this entry */
  }
  if (i > toplevel) {
      toplevel = i;
  }
  parts[i]->append(fb_cur);
}

/* Remove newline, comment, and trailing space from a line. Write a
   '\0' character to the line at the position where removed part began (i.e.,
   line gets truncated).
   Return length in characters or remaining line, without trailing '\0'
   character.
*/
size_t
read_strip_comment (char *const line)
{
    size_t linelen, i;

    linelen = strlen (line);
    if (linelen > 0 && line[linelen - 1] == '\n')
        linelen--; /* Remove newline */
    for (i = 0; i < linelen; i++) /* Skip comments */
        if (line[i] == '#') {
            linelen = i;
            break;
        }
    while (linelen > 0 && isspace((int)(unsigned char)line[linelen - 1]))
        linelen--; /* Skip whitespace at end of line */
    line[linelen] = '\0';

    return linelen;
}

/* Read a factor base file, splitting it into pieces.
   
   Primes and prime powers up to smalllim go into fb_small. If smalllim is 0,
   all primes go into fb_small, and nothing is written to fb_pieces.
   
   If smalllim is not 0, then nr_pieces separate factor bases are made for
   primes/powers > smalllim; factor base entries from the file are written to 
   these pieces in round-robin manner.

   Pointers to the allocated memory of the factor bases are written to fb_small 
   and, if smalllim > 0, to fb_pieces[0, ..., nr_pieces-1].

   Returns 1 if everything worked, and 0 if not (i.e., if the file could not be 
   opened, or memory allocation failed)
*/

int
fb_factorbase::read(const char * const filename)
{
  fb_general_entry fb_cur, fb_last;
  bool had_entry = false;
  FILE *fbfile;
  // too small linesize led to a problem with rsa768;
  // it would probably be a good idea to get rid of fgets
  const size_t linesize = 1000;
  char line[linesize];
  unsigned long linenr = 0;
  fbprime_t maxprime = 0;
  unsigned long nr_primes = 0;
  
  fbfile = fopen_maybe_compressed (filename, "r");
  if (fbfile == NULL) {
    verbose_output_print (1, 0, "# Could not open file %s for reading\n",
        filename);
    return 0;
  }
  
  while (!feof(fbfile)) {
    /* Sadly, the size parameter of fgets() is of type int */
    if (fgets (line, static_cast<int>(linesize), fbfile) == NULL)
      break;
    linenr++;
    if (read_strip_comment(line) == (size_t) 0) {
      /* Skip empty/comment lines */
      continue;
    }
    
    fb_cur.parse_line (line, linenr);
    fb_cur.invq = compute_invq(fb_cur.q);
    if (!had_entry) {
      fb_last = fb_cur;
      had_entry = true;
    } else if (fb_cur.q == fb_last.q) {
      fb_last.merge(fb_cur);
    } else {
      append(fb_last);
      fb_last = fb_cur;
    }
    
    /* fb_fprint_entry (stdout, fb_cur); */
    if (fb_cur.p > maxprime)
      maxprime = fb_cur.p;
    nr_primes++;
  }

  if (had_entry)
    append(fb_last);
  
  verbose_output_print (0, 2, "# Factor base successfully read, %lu primes, "
			"largest was %" FBPRIME_FORMAT "\n",
			nr_primes, maxprime);
  
  fclose_maybe_compressed (fbfile, filename);
  
  finalize();
  return 1;
}

/* Return p^e. Trivial exponentiation for small e, no check for overflow */
fbprime_t
fb_pow (const fbprime_t p, const unsigned long e)
{
    fbprime_t r = 1;

    for (unsigned long i = 0; i < e; i++)
      r *= p;
    return r;
}

/* Returns floor(log_2(n)) for n > 0, and 0 for n == 0 */
static unsigned int
fb_log_2 (fbprime_t n)
{
  unsigned int k;
  for (k = 0; n > 1; n /= 2, k++);
  return k;
}

/* the idea here is that we want to return the rounding of log(p^k) minus
 * the rounding of log(p^(k-1)). So it's close, in spirit, to log(p)...
 *
 * Note that the old fb_log function was taking log_scale as something
 * relative to natural logarithm. So while we had a precomputed scale so
 * that the-log-we-want(p) = log2(p) * scale, here we were feedint
 * scale/log(2), which is pretty awkward.
 *
 * (log2 is C99. I do recall that some bsd libm lacks it. It's a bug,
 * period).
 */
unsigned char
fb_log (double n, double log_scale, double offset)
{
  const long l = floor (log2 (n) * log_scale + offset + 0.5);
  return static_cast<unsigned char>(l);
}

unsigned char
fb_log_delta (const fbprime_t p, const unsigned long newexp,
              const unsigned long oldexp, const double log_scale)
{
  return fb_log (fb_pow(p, newexp), log_scale, 0.)
         - fb_log (fb_pow(p, oldexp), log_scale, 0.);
}



/* Make one factor base entry for a linear polynomial poly[1] * x + poly[0]
   and the prime (power) q. We assume that poly[0] and poly[1] are coprime.
   Non-projective roots a/b such that poly[1] * a + poly[0] * b == 0 (mod q)
   with gcd(poly[1], q) = 1 are stored as a/b mod q.
   If do_projective != 0, also stores projective roots with gcd(q, f_1) > 1,
   but stores the reciprocal root.
   Returns true if the roots was projective, and false otherwise. */

static bool
fb_linear_root (fbroot_t *root, const mpz_t *poly, const fbprime_t q)
{
  modulusul_t m;
  residueul_t r0, r1;
  bool is_projective;

  modul_initmod_ul (m, q);
  modul_init_noset0 (r0, m);
  modul_init_noset0 (r1, m);

  modul_set_ul_reduced (r0, mpz_fdiv_ui (poly[0], q), m);
  modul_set_ul_reduced (r1, mpz_fdiv_ui (poly[1], q), m);

  /* We want poly[1] * a + poly[0] * b == 0 <=>
     a/b == - poly[0] / poly[1] */
  is_projective = (modul_inv (r1, r1, m) == 0); /* r1 = 1 / poly[1] */

  if (is_projective)
    {
      ASSERT_ALWAYS(mpz_gcd_ui(NULL, poly[1], q) > 1);
      /* Set r1 = poly[0] % q, r0 = poly[1] (mod q) */
      modul_set (r1, r0, m);
      modul_set_ul_reduced (r0, mpz_fdiv_ui (poly[1], q), m);
      int rc = modul_inv (r1, r1, m);
      ASSERT_ALWAYS(rc != 0);
    }

  modul_mul (r1, r0, r1, m); /* r1 = poly[0] / poly[1] */
  modul_neg (r1, r1, m); /* r1 = - poly[0] / poly[1] */

  *root = modul_get_ul (r1, m);

  modul_clear (r0, m);
  modul_clear (r1, m);
  modul_clearmod (m);

  return is_projective;
}

class fb_powers {
  struct fb_power_t {
    fbprime_t p, q;
    unsigned char k;
  };
  std::vector<fb_power_t> *powers;

  static int
  cmp_powers(fb_power_t a, fb_power_t b)
  {
    return a.q < b.q;
  }

  public:
  fb_powers(fbprime_t);
  ~fb_powers(){delete this->powers;};
  fb_power_t &operator[] (const size_t i) {
    return (*this->powers)[i];
  }
  const fb_power_t &operator[] (const size_t i) const {
    return (*this->powers)[i];
  }
  size_t size() const {return powers->size();}
};

/* Create a list of prime powers (with exponent >1) up to lim */
fb_powers::fb_powers (const fbprime_t lim)
{
  fbprime_t p;
  static prime_info pi;
  powers = new std::vector<fb_power_t>;
  
  prime_info_init(pi);

  for (p = 2; p <= lim / p; p = getprime_mt(pi)) {
    fbprime_t q = p;
    unsigned char k = 1;
    do {
      q *= p;
      k++;
      fb_power_t new_entry = {p, q, k};
      powers->push_back(new_entry);
    } while (q <= lim / p);
  }
  prime_info_clear(pi);
  
  std::sort (powers->begin(), powers->end(), cmp_powers);
}

/* Generate a factor base with primes <= bound and prime powers <= powbound
   for a linear polynomial. If projective != 0, adds projective roots
   (for primes that divide leading coefficient).
*/

void
fb_factorbase::make_linear (const mpz_t *poly)
			    
{
  fbprime_t next_prime;
  fb_general_entry fb_cur;
  fbprime_t powlim = 0;
  
  /* Find out the largest powlim among all parts */
  for (size_t i = 0; i < FB_MAX_PARTS; i++)
    powlim = MAX(powlim, parts[i]->powlim);

  /* Prepare for computing powers up to that limit */
  fb_powers *powers = new fb_powers(powlim);
  size_t next_pow = 0;

  prime_info(pi);

  verbose_output_vfprint(0, 1, gmp_vfprintf,
               "# Making factor base for polynomial g(x) = %Zd*x%s%Zd,\n"
               "# including primes up to %" FBPRIME_FORMAT
               " and prime powers up to %" FBPRIME_FORMAT ".\n",
               poly[1], (mpz_cmp_ui (poly[0], 0) >= 0) ? "+" : "",
               poly[0], thresholds[FB_MAX_PARTS-1], powlim);

  for (next_prime = 2; next_prime <= thresholds[FB_MAX_PARTS-1]; ) {
    /* Handle any prime powers that are smaller than next_prime */
    if (next_pow < powers->size() && (*powers)[next_pow].q <= next_prime) {
      /* The list of powers must not include primes */
      ASSERT_ALWAYS((*powers)[next_pow].q < next_prime);
      fb_cur.q = (*powers)[next_pow].q;
      fb_cur.p = (*powers)[next_pow].p;
      fb_cur.k = (*powers)[next_pow].k;
      next_pow++;
    } else {
      fb_cur.q = fb_cur.p = next_prime;
      fb_cur.k = 1;
      next_prime = getprime_mt(pi);
    }
    fb_cur.nr_roots = 1;
    fb_cur.roots[0].exp = fb_cur.k;
    fb_cur.roots[0].oldexp = fb_cur.k - 1U;
    fb_cur.roots[0].proj = fb_linear_root (&fb_cur.roots[0].r, poly, fb_cur.q);
    fb_cur.invq = compute_invq(fb_cur.q);
    append(fb_cur);
  }

  prime_info_clear(pi);

  delete (powers);
  finalize();
}

/*
  Parallel version, using thread pool.
*/

#define GROUP 1024

// A task will handle 1024 entries before returning the result to the
// master thread.
// This task_info structure contains:
//   - general info (poly, number of valid entries)
//   - input for the computation
//   - output of the computation
typedef struct {
  const mpz_t *poly;
  unsigned int n;

  fbprime_t p[GROUP];
  fbprime_t q[GROUP];
  unsigned char k[GROUP];

  fbroot_t r[GROUP];
  bool proj[GROUP];
  redc_invp_t invq[GROUP];
} task_info_t;


class make_linear_thread_param: public task_parameters {
public:
  task_info_t *T;
  make_linear_thread_param(task_info_t *_T) : T(_T) {}
  make_linear_thread_param() {}
};

class make_linear_thread_result: public task_result {
public:
  task_info_t *T;
  const make_linear_thread_param *orig_param;
  make_linear_thread_result(task_info_t *_T, const make_linear_thread_param *_p)
    : T(_T), orig_param(_p) {
      ASSERT_ALWAYS(T == orig_param->T);
    }
};

static task_result *
process_one_task(const worker_thread * worker MAYBE_UNUSED, const task_parameters *_param)
{
  const make_linear_thread_param *param =
    static_cast<const make_linear_thread_param *>(_param);
  task_info_t *T = param->T;
  for (unsigned int i = 0; i < T->n; ++i) {
    T->proj[i] = fb_linear_root (&T->r[i], T->poly, T->q[i]);
    T->invq[i] = compute_invq(T->q[i]);
  }
  return new make_linear_thread_result(T, param);
}


// Prepare a new task. Return 0 if there are no new task to schedule.
// Otherwise, return the number of ideals put in the task.
static int
get_new_task(task_info_t &T, fbprime_t &next_prime, prime_info& pi, const fbprime_t maxp,
    size_t &next_pow, const fb_powers &powers)
{
  unsigned int i;
  for (i = 0; i < GROUP && next_prime <= maxp; ++i) {
    if (next_pow < powers.size() && powers[next_pow].q <= next_prime) {
      ASSERT_ALWAYS(powers[next_pow].q < next_prime);
      T.q[i] = powers[next_pow].q;
      T.p[i] = powers[next_pow].p;
      T.k[i] = powers[next_pow].k;
      next_pow++;
    } else {
      T.q[i] = T.p[i] = next_prime;
      T.k[i] = 1;
      next_prime = getprime_mt(pi);
    }
  }
  T.n = i;
  return i;
}

static void
store_task_result(fb_factorbase *fb, task_info_t *T)
{
  fb_general_entry fb_cur;
  for (unsigned int j = 0; j < T->n; ++j) {
    fb_cur.q = T->q[j];
    fb_cur.p = T->p[j];
    fb_cur.k = T->k[j];
    fb_cur.nr_roots = 1;
    fb_cur.roots[0].exp = fb_cur.k;
    fb_cur.roots[0].oldexp = fb_cur.k - 1U;
    fb_cur.roots[0].proj = T->proj[j];
    fb_cur.roots[0].r = T->r[j];
    fb_cur.invq = T->invq[j];
    fb->append(fb_cur);
  }
}

void
fb_factorbase::make_linear_threadpool (const mpz_t *poly,
    const unsigned int nb_threads)
{
  /* Find out the largest powlim among all parts */
  fbprime_t powlim = 0;
  for (size_t i = 0; i < FB_MAX_PARTS; i++)
    powlim = MAX(powlim, parts[i]->powlim);
  /* Prepare for computing powers up to that limit */
  fb_powers powers(powlim);
  size_t next_pow = 0;

  verbose_output_vfprint(0, 1, gmp_vfprintf,
               "# Making factor base for polynomial g(x) = %Zd*x%s%Zd,\n"
               "# including primes up to %" FBPRIME_FORMAT
               " and prime powers up to %" FBPRIME_FORMAT
               " using threadpool of %u threads.\n",
               poly[1], (mpz_cmp_ui (poly[0], 0) >= 0) ? "+" : "",
               poly[0], thresholds[FB_MAX_PARTS-1], powlim, nb_threads);

#define MARGIN 3
  // Prepare more tasks, so that threads keep being busy.
  unsigned int nb_tab = nb_threads + MARGIN;
  task_info_t * T = new task_info_t[nb_tab];
  make_linear_thread_param * params = new make_linear_thread_param[nb_tab];
  for (unsigned int i = 0; i < nb_tab; ++i) {
    T[i].poly = poly;
    params[i].T = &T[i];
  }
  
  fbprime_t maxp = thresholds[FB_MAX_PARTS-1];
  fbprime_t next_prime = 2;
  prime_info pi;

  prime_info_init(pi);

  thread_pool pool(nb_threads);

  // Stage 0: prepare tasks
  unsigned int active_task = 0;
  for (unsigned int i = 0; i < nb_tab; ++i) {
    int ret;
    ret = get_new_task(T[i], next_prime, pi, maxp, next_pow, powers);
    if (!ret)
      break;
    pool.add_task(process_one_task, &params[i], 0);
    active_task++;
  }

  // Stage 1: while there are still primes, wait for a result and
  // schedule a new task.
  for(int cont = 1 ; cont && active_task ; ) {
    task_result *result = pool.get_result();
    make_linear_thread_result *res =
      static_cast<make_linear_thread_result *>(result);
    active_task--;
    task_info_t * curr_T = res->T;
    store_task_result(this, curr_T);
    cont = get_new_task(*curr_T, next_prime, pi, maxp, next_pow, powers);
    if (cont) {
      active_task++;
      pool.add_task(process_one_task, res->orig_param, 0);
    }
    delete result;
  }

  // Stage 2: purge last tasks
  for (unsigned int i = 0; i < active_task; ++i) {
    task_result *result = pool.get_result();
    make_linear_thread_result *res =
      static_cast<make_linear_thread_result *>(result);
    task_info_t * curr_T = res->T;
    store_task_result(this, curr_T);
    delete result;
  }

  delete [] T;
  delete [] params;
  prime_info_clear(pi);
  finalize();
}


void
fb_factorbase::finalize()
{
  for (size_t part = 0; part < FB_MAX_PARTS; part++) {
    parts[part]->finalize();
  }
}

void
fb_factorbase::fprint(FILE *out) const
{
  for (size_t part = 0; part < FB_MAX_PARTS; part++) {
    fprintf(out, "#  Factor base entries up to %" FBPRIME_FORMAT "\n", thresholds[part]);
    parts[part]->fprint(out);
  }
}

void
fb_factorbase::_count_entries(size_t *nprimes, size_t *nroots, double *weight) const
{
  for (size_t part = 0; part < FB_MAX_PARTS; part++) {
    parts[part]->_count_entries(nprimes, nroots, weight);
  }
}

void
fb_factorbase::extract_bycost(std::vector<unsigned long> &extracted, fbprime_t pmax, fbprime_t td_thresh) const
{
  for (size_t part = 0; part < FB_MAX_PARTS; part++) {
    parts[part]->extract_bycost(extracted, pmax, td_thresh);
  }
}

/* scale is the quantity by which we multiply the log_2 of the primes we
 * consider */
void
fb_factorbase::make_slices(const double scale, const double max_weight[FB_MAX_PARTS])
{
  slice_index_t next_index = 0;
  for (size_t part = 0; part < FB_MAX_PARTS; part++) {
    parts[part]->make_slices(scale, max_weight[part], next_index);
  }
}

// Return NULL on error.
void *
fb_factorbase::mmap_fbc(void *p, void * const end)
{
  for (size_t part = 0; p != NULL && part < FB_MAX_PARTS; ++part) {
    p = parts[part]->mmap_fbc(p, end);
    size_t sz = 0;
    parts[part]->_count_entries(&sz, NULL, NULL);
    if (sz)
      toplevel = part;
  }
  return p;
}

bool
fb_factorbase::dump_fbc(FILE *f) const
{
  bool rc = true;
  for (size_t part = 0; rc && part < FB_MAX_PARTS; ++part)
    rc = parts[part]->dump_fbc(f);
  return rc;
}

void
fb_factorbase::clear()
{
  for (size_t part = 0; part < FB_MAX_PARTS; ++part)
    parts[part]->clear();
  toplevel = 0;
}

bool
fb_mmap_fbc(fb_factorbase *fb[2] MAYBE_UNUSED, const char *filename MAYBE_UNUSED)
{
#ifdef HAVE_MMAP
  FILE *f = fopen(filename, "rb");
  if (f == NULL) {
    verbose_output_print(0, 1, "# Could not open file %s for reading\n",
                         filename);
    return false;
  }

  bool  rc  = true;
  long  sz  = 0;
  void *ptr = NULL;

  if (fseek(f, 0, SEEK_END) != 0) {
    verbose_output_print(1, 0, "# Could not seek to end of file %s\n",
                         filename);
    rc = false;
    goto cleanup;
  }

  sz = ftell(f);
  if (sz < 0) {
    verbose_output_print(1, 0, "# Could not get size of file %s\n",
                         filename);
    rc = false;
    goto cleanup;
  }

  ptr = mmap(NULL, sz, PROT_READ, MAP_PRIVATE, fileno(f), 0);
  if (ptr == MAP_FAILED) {
    verbose_output_print(1, 0, "# Could not map file %s to memory\n",
                         filename);
    rc  = false;
    ptr = NULL;
    goto cleanup;
  }

  {
    void *p   = ptr;
    void *end = static_cast<char *>(ptr) + sz;
    for (size_t side = 0; rc && side < 2; ++side) {
        if (fb[side])
            rc = (p = fb[side]->mmap_fbc(p, end)) != NULL;
    }
  }

cleanup:
  if (fclose(f) != 0) {
    verbose_output_print(1, 0, "# Could not close file %s\n",
                         filename);
    rc = false;
  }

  if (!rc && ptr != NULL) {
    if (munmap(ptr, sz) != 0) {
      verbose_output_print(1, 0, "# Could not unmap file %s\n",
                           filename);
    }
    for (size_t side = 0; side < 2; ++side)
      if (fb[side]) fb[side]->clear();
  }

  return rc;
#else
  return false;
#endif
}

void
fb_dump_fbc(fb_factorbase *fb[2], const char *filename)
{
  FILE *f = fopen(filename, "wb");
  if (f == NULL) {
    verbose_output_print(1, 0, "# Could not open file %s for writing\n",
                         filename);
    return;
  }

  bool rc = true;
  for (size_t side = 0; rc && side < 2; ++side)
    if (fb[side]) rc = fb[side]->dump_fbc(f);

  if (fclose(f) != 0) {
    verbose_output_print(1, 0, "# Could not close file %s\n",
                         filename);
    rc = false;
  }

  if (!rc)
    unlink(filename);
}


#ifdef TESTDRIVE

void output(fb_factorbase *fb, const char *name)
{
  size_t n_primes = 0, n_roots = 0;
  double weight = 0.;
  fb->_count_entries(&n_primes, &n_roots, &weight);
  verbose_output_print (0, 1,
	   "# Factor base %s (%zu primes, %zu roots, %f weight):\n",
	   name, n_primes, n_roots, weight);
  fb->fprint(stdout);
}

int main(int argc, char **argv)
{
  fbprime_t thresholds[4] = {200, 1000, 1000, 1000};
  fbprime_t powbound = 100;
  fb_factorbase *fb1 = new fb_factorbase(thresholds, powbound),
    *fb2 = new fb_factorbase(thresholds, powbound);

  mpz_t poly[2];

  mpz_init(poly[0]);
  mpz_init(poly[1]);

  mpz_set_ui(poly[0], 727);
  mpz_set_ui(poly[1], 210); /* Bunch of projective primes */

  fb1->make_linear(poly);
  fb1->make_slices(2.0 * M_LN2, 0.5);
  output(fb1, "from linear polynomial");

  // This line does (and should) fail to compile, as fb_factorbase is
  // NonCopyable:
  // fb_factorbase fb3(*fb1);

  if (argc > 1) {
    fb2->read(argv[1]);
    fb2->make_slices(2.0 * M_LN2, 0.5);
    output(fb2, "from file");
  }

  delete fb1;
  delete fb2;
  bool only_general[4] = {false, false, false, false};
  fb1 = new fb_factorbase(thresholds, powbound, only_general);
  fb1->make_linear(poly);
  fb1->make_slices(2.0 * M_LN2, 0.5);
  output(fb1, "from linear polynomial, only_general = false");

  printf("Trialdiv primes:\n");
  std::vector<unsigned long> *extracted = new std::vector<unsigned long>;
  fb1->extract_bycost(*extracted, 100, 200);
  for (std::vector<unsigned long>::iterator it = extracted->begin(); it != extracted->end(); it++) {
    printf("%lu ", *it);
  }
  printf("\n");

  delete fb1;
  mpz_clear(poly[0]);
  mpz_clear(poly[1]);
  exit(EXIT_SUCCESS);
}
#endif
