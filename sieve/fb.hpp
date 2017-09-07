/*****************************************************************
*                Functions for the factor base                  *
*****************************************************************/

#ifndef FB_H
#define FB_H

#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include <gmp.h>
#include <vector>
#include <map>
#include "las-config.h"
#include "cado_poly.h" // for MAX_DEGREE
#include "fb-types.h"
#include "las-qlattice.hpp"
#include "las-plattice.hpp"
#include "las-base.hpp"

/* Forward declaration so fb_general_entry can use it in constructors */
template <int Nr_roots>
class fb_entry_x_roots;

/* A root modulo a prime power q. q is specified externally */
struct fb_general_root {
  /* exp and oldexp are maximal such that:
     If not projective and a == br (mod p^k), then p^exp | F(a,b)
     -"-                   a == br (mod p^(k-1)), then p^oldexp | F(a,b)
     If projective and ar == b  -"- */
  fbroot_t r;
  bool proj;
  unsigned char exp, oldexp;

  fb_general_root (){}
  fb_general_root (const fbroot_t r, const unsigned char nexp=1,
                   const unsigned char oldexp=0, const bool proj=false) :
                   r(r), proj(proj), exp(nexp), oldexp(oldexp) {}
  /* Create a root from a linear polynomial */
  fb_general_root (fbprime_t q, const mpz_t *poly, const unsigned char nexp=1,
                   const unsigned char oldexp=0);

  /* Constructor from the old format of storing projective roots, which has q
     added to the root if the root is projective */
  fb_general_root (const unsigned long long old_r, const fbprime_t q,
                   const unsigned char nexp=1, const unsigned char oldexp=0) :
                   r((old_r >= q) ? (old_r - q) : old_r),
                   proj(old_r >= q), exp(nexp), oldexp(oldexp) {}

  /* A root is simple if it is not projective and the exp goes from 0 to 1 */
  bool is_simple() const {return exp == 1 && oldexp == 0 && !proj;}

  /* Convert a root to the old format of storing projective roots with q added */
  unsigned long long to_old_format(const fbprime_t q) const {
    return (unsigned long long) r + (proj ? q : 0);
  }

  /* Print one root. Projective roots are printed as r+q */
  void fprint(FILE *out, const fbprime_t q) const {
    fprintf(out, "%llu", to_old_format(q));
    if (oldexp != 0 || this->exp != 1)
      fprintf(out, ":%hhu:%hhu", oldexp, this->exp);
  }

  void transform(fb_general_root &result, const fbprime_t q,
                 const redc_invp_t invq,
                 const qlattice_basis &basis) const {
    unsigned long long t = to_old_format(q);
    t = fb_root_in_qlattice(q, t, invq, basis);
    result = fb_general_root(t, q, exp, oldexp);
  }
};

/* General entries are anything that needs auxiliary information:
   Prime powers, projective roots, ramified primes where exp != oldexp + 1,
   etc. They could, of course, also store the simple cases, but for those we
   use the simple struct to conserve memory and to decide algorithms (batch
   inversion, etc.) statically. */
class fb_general_entry : public _padded_pod<fb_general_entry> {
  void read_roots (const char *, unsigned char, unsigned char, unsigned long);
public:
  typedef fb_general_entry transformed_entry_t;
  fbprime_t q, p; /* q = p^k */
  redc_invp_t invq; /* invq = -1/q (mod 2^32), or (mod 2^64), depending on
		       the size of redc_invp_t */
  fb_general_root roots[MAX_DEGREE];
  unsigned char k, nr_roots;
  /* Static class members to allow fb_vector<> to distinguish between and
     operate on both kind of entries */
  static const bool is_general_type = true;
  static const unsigned char fixed_nr_roots = 0;

  fb_general_entry() {}
  template <int Nr_roots>
  fb_general_entry (const fb_entry_x_roots<Nr_roots> &e);
  fbprime_t get_q() const {return q;}
  fbroot_t get_r(const size_t i) const {return roots[i].r;};
  fbroot_t get_proj(const size_t i) const {return roots[i].proj;};
  void parse_line (const char *line, unsigned long linenr);
  void merge (const fb_general_entry &);
  void fprint(FILE *out) const;
  bool is_simple() const;
  void transform_roots(transformed_entry_t &, const qlattice_basis &) const;
  double weight() const {return 1./q * nr_roots;}
  void extract_bycost(std::vector<unsigned long> &extracted, fbprime_t pmax, fbprime_t td_thresh) const {
    if (k == 1 && p <= pmax && p <= nr_roots * td_thresh) {
      // printf("Extracting p = %" FBPRIME_FORMAT "\n", p);
      extracted.push_back(static_cast<unsigned long>(p));
    }
  }
  /* Allow sorting by p */
  bool operator<(const fb_general_entry &other) const {return this->p < other.p;}
  bool operator>(const fb_general_entry &other) const {return this->p > other.p;}
};

template <int Nr_roots>
class fb_transformed_entry_x_roots {
public:
  fbprime_t p;
  fbroot_t roots[Nr_roots];
  bool proj[Nr_roots];
  static const unsigned char k = 1, nr_roots = Nr_roots;
  /* Static class members to allow fb_vector<> to distinguish between and
     operate on both kind of entries */
  static const bool is_general_type = false;
  static const unsigned char fixed_nr_roots = Nr_roots;
  fbprime_t get_q() const {return p;}
  fbroot_t get_r(const size_t i) const {return roots[i];};
  fbroot_t get_proj(const size_t i) const {return proj[i];};
};

/* "Simple" factor base entries. We imply q=p, k=1, oldexp=0, exp=1,
   and projective=false for all roots. */
template <int Nr_roots>
class fb_entry_x_roots : public _padded_pod<fb_entry_x_roots<Nr_roots> > {
public:
  typedef fb_transformed_entry_x_roots<Nr_roots> transformed_entry_t;
  fbprime_t p;
  redc_invp_t invq; /* invq = -1/q (mod 2^32), or (mod 2^64), depending on
		       the size of redc_invp_t */
  fbroot_t roots[Nr_roots];
  /* Static class members to allow fb_vector<> to distinguish between and
     operate on both kind of entries */
  static const unsigned char k = 1, nr_roots = Nr_roots;
  static const bool is_general_type = false;
  static const unsigned char fixed_nr_roots = Nr_roots;
  fb_entry_x_roots() {};
  /* Allow assignment-construction from general entries */
  fb_entry_x_roots(const fb_general_entry &e) : p(e.p), invq(e.invq) {
    ASSERT_ALWAYS(Nr_roots == e.nr_roots);
    ASSERT(e.is_simple());
    for (int i = 0; i < Nr_roots; i++)
      roots[i] = e.roots[i].r;
  }
  fbprime_t get_q() const {return p;}
  fbroot_t get_r(const size_t i) const {return roots[i];};
  fbroot_t get_proj(const size_t i MAYBE_UNUSED) const {return false;};
  double weight() const {return 1./p * Nr_roots;}
  /* Allow sorting by p */
  bool operator<(const fb_entry_x_roots<Nr_roots> &other) const {return this->p < other.p;}
  bool operator>(const fb_entry_x_roots<Nr_roots> &other) const {return this->p > other.p;}
  void fprint(FILE *) const;
  void transform_roots(transformed_entry_t &, const qlattice_basis &) const;
  void extract_bycost(std::vector<unsigned long> &extracted, fbprime_t pmax, fbprime_t td_thresh) const {
    if (p <= pmax && p <= Nr_roots * td_thresh)
      extracted.push_back(static_cast<unsigned long>(p));
  }
};


class fb_countable_entries {
  friend class fb_part;
  virtual void _count_entries(size_t *nprimes, size_t *nroots, double *weight)
    const = 0;
  public:
  virtual ~fb_countable_entries(){}
  void count_entries(size_t *nprimes, size_t *nroots, double *weight) const {
    if (nprimes != NULL) *nprimes = 0;
    if (nroots != NULL) *nroots = 0;
    if (weight != NULL) *weight = 0;
    _count_entries(nprimes, nroots, weight);
  }
};

class fb_interface : public fb_countable_entries {
public:
  virtual ~fb_interface(){}
  virtual void append(const fb_general_entry &) = 0;
  virtual void fprint(FILE *) const = 0;
};


class fb_transformed_vector: 
      public std::vector<plattice_sieve_entry>, private NonCopyable {
  slice_index_t index;
public:
  fb_transformed_vector(const slice_index_t index) : index(index) {}
  slice_index_t get_index() const {return index;};
};


/* The work queue of slices that need to be sieved must return a common base
   class of fb_slice<FB_ENTRY_TYPE>; this base class is fb_slice_interface.
   It allows determining which kind of slice it is, so that the correct
   siever implementation can be called. It also allows generating a vector
   of general entries with the transformed roots. */
class fb_slice_interface {
  public:
  fb_slice_interface(){}
  virtual ~fb_slice_interface(){}
  virtual int get_nr_roots() const = 0;
  virtual bool is_general() const = 0;
  virtual plattices_vector_t * make_lattice_bases(const qlattice_basis &, int, const sublat_t &) const = 0;
  virtual unsigned char get_logp() const = 0;
  virtual slice_index_t get_index() const = 0;
  virtual fbprime_t get_prime(slice_offset_t offset) const = 0;
  virtual unsigned char get_k(slice_offset_t offset) const = 0;
  virtual double get_weight() const = 0;
};

class fb_vector_interface: public fb_interface {
  virtual void *mmap_fbc(void *p, void * const end) = 0;
  virtual bool  dump_fbc(FILE *f) const = 0;
  virtual void  clear() = 0;
  public:
  virtual ~fb_vector_interface(){}
  virtual const fb_slice_interface *get_slice(size_t) const = 0;
  virtual int get_nr_roots() const = 0;
  virtual bool is_general() const = 0;
  virtual void extract_bycost(std::vector<unsigned long> &extracted, fbprime_t pmax, fbprime_t td_thresh) const = 0;
  virtual void finalize() = 0;
  /* Create slices so that one slice contains at most max_slice_len entries,
   * and all entries in a slice have the same round(fb_log(p, scale)).
   * scale is the quantity by which we multiply the log_2 of the primes we
   * consider.
   */
  virtual void make_slices(double scale, double max_weight,
                           slice_index_t &next_index) = 0;
  virtual const fb_slice_interface *get_first_slice() const = 0;
  friend class fb_part;
};


/* A slice is a range of consecutive entries in an fb_vector. 
   The pointers of slice into the vector's data are "weak", i.e., they do not
   determine the pointed-to data's lifetime. In fact, slices are re-computed
   each time the logarithm base changes, and the slice pointers move around
   accordingly.
   Each slice in a factor base has a unique index. This index will be used to
   reconstruct factor base primes from the hint: each update stores a hint
   which is an offset within a slice, and the slice index is stored for ranges
   of updates, allowing the factor base entry to be looked up via (roughly)
   *(factorbase.get_slice(index) + offset), except the concrete type of the
   pointed-to data needs to be determined to be able to calculate the correct
   memory address in the "+ offset" operation.  */
template <class FB_ENTRY_TYPE>
class fb_slice : public fb_slice_interface {
  const FB_ENTRY_TYPE *_begin, *_end;
  unsigned char logp;
  slice_index_t index;
  double weight;
  public:
  typedef typename FB_ENTRY_TYPE::transformed_entry_t transformed_entry_t;
  fb_slice(const FB_ENTRY_TYPE *begin,
           const FB_ENTRY_TYPE *end,
           const unsigned char logp,
           const slice_index_t index,
           const double weight)
           : _begin(begin), _end(end), logp(logp), index(index),
             weight(weight) {}
  /* Iterators */
  const FB_ENTRY_TYPE *begin() const {return _begin;}
  const FB_ENTRY_TYPE *end() const {return _end;}
  /* Implement the fb_slice_interface */
  int get_nr_roots() const {return FB_ENTRY_TYPE::fixed_nr_roots;}
  bool is_general() const {return FB_ENTRY_TYPE::is_general_type;}
  unsigned char get_logp() const {return logp;};
  slice_index_t get_index() const {return index;}
  double get_weight() const {return weight;}
  void fprint(FILE *out) const;
  fbprime_t get_prime(const slice_offset_t offset) const {
    ASSERT_ALWAYS(_begin + offset < _end);
    return _begin[offset].p;
  }
  unsigned char get_k(const slice_offset_t offset) const {
    ASSERT_ALWAYS(_begin + offset < _end);
    return _begin[offset].k;
  }
  plattices_vector_t * make_lattice_bases(const qlattice_basis &, int, const sublat_t &) const;
};

/* A predicate to test whether an fb_slice has weight at most x, for use with
   STL search functions etc. */
class fb_slice_pred_weight_le {
  const double max_weight;
public:
  fb_slice_pred_weight_le (const double x) : max_weight(x) {}
  bool operator() (const fb_slice_interface *slice) const {return slice->get_weight() <= max_weight;}
};


template <class FB_ENTRY_TYPE>
class fb_vector : public fb_vector_interface, private NonCopyable {
  FB_ENTRY_TYPE *data;
  size_t         size, alloc;
  bool           read_only, mmapped;

  typedef std::vector<fb_slice<FB_ENTRY_TYPE> > slices_t;
  const slices_t *slices; /* If no slicing exists, e.g., because the fb_part
                             is only_general (and thus meant for line-sieving)
                             or because the vector is empty, then
                             slices == NULL */
  std::map<const double, const slices_t *> cached_slices;
  void sort();
  double est_weight_max(size_t, size_t) const;
  double est_weight_avg(size_t, size_t) const;
  double est_weight_simpson(size_t, size_t) const;
  double est_weight_mertens(size_t, size_t) const;
  double est_weight_sum(size_t, size_t) const;
  double est_weight_compare(size_t, size_t) const;
  double est_weight(size_t, size_t) const;
  void *mmap_fbc(void *p, void * const end);
  bool  dump_fbc(FILE *f) const;
  void  clear();
  void  init();
 public:
  static const size_t max_slice_len = 65536;

  fb_vector() : slices(NULL) { init(); }
  ~fb_vector();

  void make_slices(double scale, double max_weight, slice_index_t &next_index);

  typedef const FB_ENTRY_TYPE *const_iterator;

  const_iterator begin() const { return data;      }
  const_iterator end()   const { return data+size; }

  /* Implement fb_vector interface */
  int get_nr_roots() const {return FB_ENTRY_TYPE::fixed_nr_roots;}
  bool is_general() const {return FB_ENTRY_TYPE::is_general_type;}

  void append(const fb_general_entry &new_entry) {
    ASSERT_ALWAYS(!read_only);
    if (size >= alloc) {
      alloc += std::max(alloc, static_cast<size_t>(16));
      data   = static_cast<FB_ENTRY_TYPE *>(realloc(data, alloc * sizeof(FB_ENTRY_TYPE)));
      if (data == NULL)
        throw std::bad_alloc();
    }
    memset(data+size, 0, sizeof(FB_ENTRY_TYPE));
    data[size++] = new_entry;
  }
  
  void _count_entries(size_t *nprimes, size_t *nroots, double *weight) const;
  void fprint(FILE *out) const;
  void extract_bycost(std::vector<unsigned long> &extracted, fbprime_t pmax,
                      fbprime_t td_thresh) const;

  /* Finalizing requires only sorting the vector */
  void finalize() {
    ASSERT_ALWAYS(!read_only);
    alloc = size;
    data  = static_cast<FB_ENTRY_TYPE *>(realloc(data, alloc * sizeof(FB_ENTRY_TYPE)));
    sort();
    read_only = true;
  }

  const fb_slice<FB_ENTRY_TYPE> *get_first_slice() const {
    if (slices == NULL || slices->empty()) {
      return NULL;
    }
    return &(*slices)[0];
  }
  const fb_slice<FB_ENTRY_TYPE> *get_slice(const size_t n) const {
    if (slices == NULL || slices->empty())
      return NULL;
    const slice_index_t first = slices->front().get_index();
    const slice_index_t last = slices->back().get_index();
    const fb_slice<FB_ENTRY_TYPE> *slice = NULL;
    if (first <= n && n <= last) {
      slice = &(*slices)[n - first];
      ASSERT_ALWAYS(n == slice->get_index());
    }
    return slice;
  }
  friend class fb_part;
};

/* A "part" is the set of all factor base primes that get sieved over a given
   bucket region size.
   E.g., when we have only 1 level of bucket sorting, then the factor base has
   2 parts: the line-sieved primes, and the bucket-sieved primes. */
class fb_part: public fb_interface, private NonCopyable {
  
  /* How do we identify the number of slices in each array?
     Should we have
     size_t nr_entries[9];
     here, or let each array end with a "magic value", such as a slice with
     0 entries, or a NULL pointer? */

  const fbprime_t powlim;
  /* If true, all entries go in the general vector */
  const bool only_general;

  /* These vectors are filled when we read or generate the factor base.
     The slices point into the vectors' storage. */
  fb_vector<fb_entry_x_roots<0> > fb0_slices; /* From 0 to MAX_DEGREE */
  fb_vector<fb_entry_x_roots<1> > fb1_slices;
  fb_vector<fb_entry_x_roots<2> > fb2_slices;
  fb_vector<fb_entry_x_roots<3> > fb3_slices;
  fb_vector<fb_entry_x_roots<4> > fb4_slices;
  fb_vector<fb_entry_x_roots<5> > fb5_slices;
  fb_vector<fb_entry_x_roots<6> > fb6_slices;
  fb_vector<fb_entry_x_roots<7> > fb7_slices;
  fb_vector<fb_entry_x_roots<8> > fb8_slices;
  fb_vector<fb_entry_x_roots<9> > fb9_slices;
  fb_vector<fb_entry_x_roots<10> > fb10_slices;
  fb_vector<fb_general_entry> general_vector;

  fb_vector_interface *get_slices(const unsigned int n) {
    ASSERT_ALWAYS(n <= MAX_DEGREE);
    
    if (only_general)
     return NULL;
    
    switch (n) {
      case 0: return &fb0_slices;
      case 1: return &fb1_slices;
      case 2: return &fb2_slices;
      case 3: return &fb3_slices;
      case 4: return &fb4_slices;
      case 5: return &fb5_slices;
      case 6: return &fb6_slices;
      case 7: return &fb7_slices;
      case 8: return &fb8_slices;
      case 9: return &fb9_slices;
      case 10: return &fb10_slices;
      default: abort();
    }
  }
  const fb_vector_interface *get_slices(const unsigned int n) const {
    ASSERT_ALWAYS(n <= MAX_DEGREE);
    
    if (only_general)
     return NULL;
    
    switch (n) {
      case 0: return &fb0_slices;
      case 1: return &fb1_slices;
      case 2: return &fb2_slices;
      case 3: return &fb3_slices;
      case 4: return &fb4_slices;
      case 5: return &fb5_slices;
      case 6: return &fb6_slices;
      case 7: return &fb7_slices;
      case 8: return &fb8_slices;
      case 9: return &fb9_slices;
      case 10: return &fb10_slices;
      default: abort();
    }
  }
  void make_slices(double scale, double max_weight, slice_index_t &next_index);
  void *mmap_fbc(void *p, void * const end);
  bool  dump_fbc(FILE *f) const;
  void  clear();
public:
  fb_part(const fbprime_t _powlim, const bool only_general=false)
    : powlim(_powlim), only_general(only_general){}
  ~fb_part(){}
  void append(const fb_general_entry &);
  void fprint(FILE *) const;
  void _count_entries(size_t *nprimes, size_t *nroots, double *weight) const;
  void finalize();
  bool is_only_general() const {return only_general;}
  const fb_vector<fb_general_entry> *get_general_vector(){ return &general_vector; }
  const fb_slice_interface *get_first_slice() const {
    /* Find the first non-empty slices entry and return a pointer to it,
       or return NULL if all are empty */
    const fb_slice_interface *slice = NULL;
    if (!only_general) {
      for (unsigned int nr_roots = 0; slice == NULL && nr_roots <= MAX_DEGREE; nr_roots++) {
        const fb_vector_interface *slices = get_slices(nr_roots);
        if (slices != NULL)
          slice = slices->get_first_slice();
      }
    }
    if (slice == NULL)
      slice = general_vector.get_first_slice();

    return slice;
  }
  /* Returns the index of the first non-empty slice in this part.
     If there are no slices, i.e., this part is empty, returns 0.  */
  slice_index_t get_first_slice_index() const {
    const fb_slice_interface *slice = get_first_slice();
    return (slice == NULL) ? 0 : slice->get_index();
  }
  const fb_slice_interface *get_slice(const slice_index_t slice_idx) const {
    for (unsigned int nr_roots = 0; nr_roots <= MAX_DEGREE; nr_roots++) {
      const fb_vector_interface *slices = get_slices(nr_roots);
      if (slices != NULL) {
        const fb_slice_interface *slice;
        if ((slice = slices->get_slice(slice_idx)) != NULL) return slice;
      }
    }
    return general_vector.get_slice(slice_idx);
  }
  void extract_bycost(std::vector<unsigned long> &p, fbprime_t pmax, fbprime_t td_thresh) const;
  friend class fb_factorbase;
};

/* Splits the factor base for a polynomial into disjoint parts which are
   sieved over different sieve region sizes.
   For example,
   parts[0] will contain very small primes that are line-sieved,
   parts[1] contains bucket-sieved primes with 1 level of bucket sorting
   (i.e., hits get sorted into bucket regions of size 2^16)
*/
class fb_factorbase: public fb_interface, private NonCopyable {
  fb_part *parts[FB_MAX_PARTS];
  fbprime_t thresholds[FB_MAX_PARTS];
  void finalize();
  int toplevel;
  void *mmap_fbc(void *p, void * const end);
  bool  dump_fbc(FILE *f) const;
  void  clear();
 public:
  fb_factorbase(const fbprime_t *thresholds,
                fbprime_t powlim,
		const bool *only_general=NULL);
  ~fb_factorbase();
  int read(const char * const filename);
  void make_linear (const mpz_t *poly);
  void make_linear_parallel (const mpz_t *poly, const unsigned int nb_threads);
  void make_linear_threadpool(const mpz_t *poly, const unsigned int nb_threads);
  size_t size() const {return 0;}
  void fprint(FILE *) const;
  void append(const fb_general_entry &);
  void _count_entries(size_t *nprimes, size_t *nroots, double *weight) const;
  fb_part *get_part(const size_t n) const {ASSERT_ALWAYS(n < FB_MAX_PARTS); return parts[n];}
  void extract_bycost(std::vector<unsigned long> &extracted, fbprime_t pmax, fbprime_t td_thresh) const;
  void make_slices(double scale, const double max_weight[FB_MAX_PARTS]);
  const fb_slice_interface *get_slice(const slice_index_t slice_idx) const {
    for (unsigned int i_part = 0; i_part <= FB_MAX_PARTS; i_part++) {
      const fb_slice_interface *slice;
      if ((slice = parts[i_part]->get_slice(slice_idx)) != NULL)
        return slice;
    }
    return NULL;
  }
  int get_toplevel() const { return toplevel; }

  friend bool fb_mmap_fbc(fb_factorbase *[2], const char *);
  friend void fb_dump_fbc(fb_factorbase *[2], const char *);
};


/* round(x*y-z) */
unsigned char	fb_log (double x, double y, double z);
fbprime_t       fb_pow (fbprime_t, unsigned long);
fbprime_t       fb_is_power (fbprime_t, unsigned long *);
void print_worst_weight_errors();

bool fb_mmap_fbc(fb_factorbase *fb[2], const char *filename);
void fb_dump_fbc(fb_factorbase *fb[2], const char *filename);

#endif
