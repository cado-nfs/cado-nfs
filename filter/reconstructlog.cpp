#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cinttypes>    // for PRIu64, SCNu64
#include <cstdint>      // for uint64_t, int32_t, int64_t, uint8_t
#include <pthread.h>     // for pthread_mutex_lock, pthread_mutex_unlock
#ifdef HAVE_MINGW
#include <fcntl.h>   /* for _O_BINARY */
#endif
#include <gmp.h>
#include "bit_vector.h"  // for bit_vector_set, bit_vector, bit_vector_getbit
#include "cado_poly.h"  // cado_poly
#include "cxx_mpz.hpp"
#include "filter_io.h"  // earlyparsed_relation_ptr
#include "gmp_aux.h"     // for nbits, mpz_addmul_si
#include "gzip.h"       // fopen_maybe_compressed
#include "macros.h"
#include "memalloc.h"             // my_malloc_free_all
#include "mpz_poly.h"    // for mpz_poly_setcoeff_int64, mpz_poly, mpz_poly_...
#include "params.h"
#include "purgedfile.h" // purgedfile_read_firstline
#include "renumber.hpp"
#include "sm_utils.h"   // sm_side_info_clear
#include "stats.h"                     // stats_data_t
#include "timing.h"      // for wct_seconds
#include "typedefs.h"    // for ideal_merge_t, index_t, weight_t, PRid, prime_t
#include "verbose.h"    // verbose_decl_usage

#define DEBUG 0

stats_data_t stats; /* struct for printing progress */

/*********************** mutex for multi threaded version ********************/
/* used as mutual exclusion lock for reading the status of logarithms */
pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;

/**** Relations structure used for computing the logarithms from the rels ****/
typedef struct
{
  weight_t nb_unknown;
  ideal_merge_t *unknown;
  mpz_t log_known_part;
} log_rel_t;

/* Init a table of log_rel_t of size nrels (malloc + all mpz_init) */
static log_rel_t *
log_rel_init (uint64_t nrels)
{
  log_rel_t *rels;
  uint64_t i;
  rels = (log_rel_t *) malloc (nrels * sizeof (log_rel_t));
  FATAL_ERROR_CHECK(rels == NULL, "Cannot allocate memory");
  memset(rels, 0, nrels * sizeof(log_rel_t));
  for (i = 0; i < nrels; i++)
    mpz_init(rels[i].log_known_part);
  return rels;
}

/* Free what is allocated by log_rel_init */
static void
log_rel_free (log_rel_t *rels, uint64_t nrels)
{
  uint64_t i;
  my_malloc_free_all();
  for (i = 0; i < nrels; i++)
    mpz_clear (rels[i].log_known_part);
  free(rels);
}

/***** Light relations structure (for constructing the dependency graph) *****/
typedef struct
{
  weight_t len;
  index_t *needed;
} light_rel;

typedef light_rel *light_rels_t;

/* Init a light_rels_t of size nrels */
static light_rels_t
light_rels_init (uint64_t nrels)
{
  light_rels_t rels;
  rels = (light_rels_t) malloc (nrels * sizeof (light_rel));
  FATAL_ERROR_CHECK(rels == NULL, "Cannot allocate memory");
  memset(rels, 0, nrels * sizeof(light_rel));
  return rels;
}

/* Free what is allocated by light_rel_init (and during reading for the
 * "unknown" array in the light_rel structure */
static void
light_rels_free (light_rels_t rels)
{
  my_malloc_free_all();
  free(rels);
}

/****************** Struct for the tab of log (mpz_t *) *********************/
struct logtab
{
    uint64_t nprimes;
    uint64_t nknown;
    int nbsm;
    cxx_mpz ell;
    cado_poly_srcptr cpoly;
    sm_side_info *sm_info;
    private:
    mp_limb_t * data;
    public:
    struct mpz_ro_accessor {
        /* This is a ***NON-OWNING*** mpz-like accessor for an integer
         * stored somewhere.
         *
         * A typical construct might be to use this as in:
         *
         * mpz_addmul(foo, log[i], bar);
         *
         * here log[i] would create an mpz_ro_accessor object, and then
         * call its mpz_srcptr converter (member function below) to
         * return an mpz_srcptr that is appropriate to pass to
         * mpz_addmul.
         *
         * The fine point is whether the dtor for the temporary
         * mpz_ro_accessor object is called before or after entering the
         * mpz_addmul function.
         *
         * The C++ standard says: after. [C++11 § 12.2.3] Namely:
         *
         *      When an implementation introduces a temporary object of a
         *      class that has a non-trivial constructor (12.1, 12.8), it
         *      shall ensure that a constructor is called for the
         *      temporary object. Similarly, the destructor shall be
         *      called for a temporary with a non-trivial destructor
         *      (12.4). Temporary objects are destroyed as the last step
         *      in evaluating the full-expression (1.9) that (lexically)
         *      contains the point where they were created.
         *
         * Bottom line: the construct mpz_addmul(foo, log[i], bar) is
         * safe. (However mpz_srcptr z = log[i]; followed by
         * mpz_addmul(foo, z, bar) is not !)
         */
        private:
        static inline int mpz_normalized_size(mp_limb_t * p, mp_size_t n) {
            for( ; n && p[n-1] == 0 ; n--);
            return n;
        }
        protected:
        mpz_t ugly;
        public:
        /* The following construct would be valid with gmp-6+
         * MPZ_ROINIT_N function, but unfortunately MPZ_ROINIT_N sets the
         * _mp_alloc field to zero, which gets in the way of our usage in
         * the rw_accessor below.
        mpz_ro_accessor(logtab const & log, mp_limb_t * place) : ugly MPZ_ROINIT_N(place, mpz_normalized_size(place, mpz_size(log.ell))) {}
         */
        mpz_ro_accessor(logtab const & log, mp_limb_t * place) {
            ugly->_mp_d = place;
            ugly->_mp_alloc = mpz_size(log.ell);
            ugly->_mp_size = mpz_normalized_size(place, mpz_size(log.ell));
        }
        operator mpz_srcptr() const { return ugly; }
    };
    struct mpz_rw_accessor : public mpz_ro_accessor {
        private:
        /* This is ***NON-OWNING*** as well, but has operator= which
         * carries over to the parent structure */
        uint64_t h;
        logtab & log;
        public:
        bool is_known() const { return log.is_known(h); }
        mpz_rw_accessor(logtab & log, uint64_t h, mp_limb_t * place) : mpz_ro_accessor(log, place), h(h), log(log) {
        }
        mpz_rw_accessor& operator=(mpz_srcptr v) {
            if (is_known()) {
                gmp_fprintf (stderr, "ERROR, inconsistent log for h = %" PRIu64 " ; we previously had %Zd in the database, now we want to store %Zd\n", h,
                        (mpz_srcptr) *this, v);
                ASSERT_ALWAYS (mpz_cmp ((mpz_srcptr) *this, v) == 0);
                return *this;
            }
            if (mpz_cmp_ui (v, 0) < 0)
            {
                fprintf (stderr, "Warning, log is negative for h = %" PRIu64 "\n", h);
                cxx_mpz vv;
                mpz_mod(vv, v, log.ell);
                (*this) = vv;
            } else if (mpz_cmp (v, log.ell) >= 0) {
                fprintf (stderr, "Warning, log >= ell for h = %" PRIu64 "\n", h);
                cxx_mpz vv;
                mpz_mod(vv, v, log.ell);
                (*this) = vv;
            } else {
                ASSERT_ALWAYS(mpz_size(v) <= mpz_size(log.ell));
                mpn_zero(ugly->_mp_d, mpz_size(log.ell));
                if (mpz_cmp_ui (v, 0) == 0) {
                    fprintf (stderr, "Warning, log is zero for h = %" PRIu64 "\n", h);
                } else {
                    mpn_copyi(ugly->_mp_d, v->_mp_d, mpz_size(v));
                }
                // log of SM columns are not taken into account
                if (h < log.nprimes)
                    log.nknown++;
            }
            return *this;
        }
    };
    private:
    uint64_t smlog_index(int side, int idx_sm) const {
        uint64_t h = nprimes;
        for(int i = 0 ; i < side ; i++) {
            h += sm_info[i]->nsm;
        }
        return h + idx_sm;
    }
    public:
    bool is_zero(uint64_t h) const {
        mpz_ro_accessor z = (*this)[h];
        return mpz_cmp_ui ((mpz_srcptr) z, 0) == 0;
    }
    bool is_known(uint64_t h) const {
        mpz_ro_accessor z = (*this)[h];
        return mpz_cmp ((mpz_srcptr) z, ell) < 0;
    }
    mpz_ro_accessor operator[](uint64_t h) const {
        return mpz_ro_accessor(*this, data + h * mpz_size(ell));
    }
    mpz_rw_accessor operator[](uint64_t h) {
        return mpz_rw_accessor(*this, h, data + h * mpz_size(ell));
    }
    mpz_ro_accessor smlog(int side, int idx_sm) const {
        return (*this)[smlog_index(side, idx_sm)];
    }
    mpz_rw_accessor smlog(int side, int idx_sm) {
        return (*this)[smlog_index(side, idx_sm)];
    }
    void force_set(uint64_t h, mpz_srcptr v)
    {
        mp_limb_t * p = data + h * mpz_size(ell);
        ASSERT_ALWAYS(mpz_size(v) <= mpz_size(ell));
        mpn_zero(p, mpz_size(ell));
        mpn_copyi(p, v->_mp_d, mpz_size(v));
    }

    logtab(cado_poly_srcptr cpoly, sm_side_info * sm_info, uint64_t nprimes, mpz_srcptr ell)
        : nprimes(nprimes)
          , cpoly(cpoly)
          , sm_info(sm_info)
    {
        mpz_set(this->ell, ell);
        nknown = 0;
        nbsm = 0;
        for(int side = 0 ; side < cpoly->nb_polys ; side++) {
            nbsm += sm_info[side]->nsm;
        }
        data = (mp_limb_t *) malloc((nprimes + nbsm) * mpz_size(ell) * sizeof(mp_limb_t));
        FATAL_ERROR_CHECK(data == NULL, "Cannot allocate memory");
        /* set everything to the max value */
        memset(data, -1, (nprimes + nbsm) * mpz_size(ell) * sizeof(mp_limb_t));
    }
    logtab(logtab const &) = delete;
    ~logtab() { free(data); }
};

/************ Struct used for reading rels files with process_rels ***********/
struct read_data
{
  log_rel_t *rels;
  logtab & log;
  uint64_t nrels;
  cado_poly_ptr cpoly;
  sm_side_info *sm_info;
  renumber_t & renum_tab;
  read_data(logtab & log, uint64_t nrels,
          cado_poly_ptr cpoly, sm_side_info *sm_info,
          renumber_t & renum_tab)
      : log(log)
      , nrels(nrels)
      , cpoly(cpoly)
      , sm_info(sm_info)
      , renum_tab(renum_tab)
  {
      rels = log_rel_init (nrels);
      fflush(stdout);
  }
  read_data(read_data const &) = delete;
  ~read_data() {
      log_rel_free (rels, nrels);
  }
};

/************************ Handling of the SMs *******************************/
/* number of SM that must be used. */

/* Callback function called by filter_rels in compute_log_from_rels */
void *
thread_sm (void * context_data, earlyparsed_relation_ptr rel)
{
    read_data & data = * (read_data *) context_data;
    log_rel_t *lrel = &(data.rels[rel->num]);

    mpz_ptr l = lrel->log_known_part;
    int64_t a = rel->a;
    uint64_t b = rel->b;

    uint64_t nonvoidside = 0; /* bit vector of which sides appear in the rel */
    if (data.cpoly->nb_polys > 2) {
        for (weight_t i = 0; i < rel->nb; i++) {
          index_t h = rel->primes[i].h;
          int side = data.renum_tab.p_r_from_index(h).side;
          nonvoidside |= ((uint64_t) 1) << side;
        }
        /* nonvoidside must *not* be a power of two. If it is, then we
         * have a nasty problem similar to bug 21707: in a sense, we have
         * true gem of a relation that yields a trivial norm on one side,
         * but it's really too bad that we have no effective way to check
         * for it. */
        ASSERT_ALWAYS(nonvoidside & (nonvoidside - 1));
        /* one thing we might do at this point is recompute the norm from
         * a, b, and data.cpoly->pols[side], and see if we get \pm1.
         */
    } else {
        nonvoidside = 3;
    }

    if (rel->sm_size) {
        /* use the SM values which are already present in the input file,
         * because some goodwill computed them for us.
         */
        int c = 0;
        for(int side = 0 ; side < data.cpoly->nb_polys ; side++) {
            sm_side_info_srcptr S = data.sm_info[side];
            if (S->nsm > 0 && (nonvoidside & (((uint64_t) 1) << side))) {
#define xxxDOUBLECHECK_SM
#ifdef DOUBLECHECK_SM
                /* I doubt that this is really compatible with our
                 * changes in the SM mode.
                 */
                mpz_poly u;
                mpz_poly_init(u, MAX(1, S->f->deg-1));
                mpz_poly_setcoeff_int64(u, 0, a);
                mpz_poly_setcoeff_int64(u, 1, -b);
                compute_sm_piecewise(u, u, S);
                ASSERT_ALWAYS(u->deg < S->f->deg);
                ASSERT_ALWAYS(u->deg == S->f->deg - 1);
                for(int i = 0 ; i < S->nsm; i++) {
                    if (S->mode == SM_MODE_LEGACY_PRE2018)
                        ASSERT_ALWAYS(mpz_cmp(u->coeff[S->f->deg-1-i],rel->sm[c + i
]) == 0);
                    else
                        ASSERT_ALWAYS(mpz_cmp(u->coeff[i],rel->sm[c + i]) == 0);
                }

#endif
                ASSERT_ALWAYS(c + S->nsm <= rel->sm_size);
                for(int i = 0 ; i < S->nsm ; i++, c++) {
                    mpz_addmul(l, data.log.smlog(side, i), rel->sm[c]);
                }
                mpz_mod(l, l, data.log.ell);
#ifdef DOUBLECHECK_SM
                mpz_poly_clear(u);
#endif
            }
        }
    } else {
        mpz_srcptr ell = data.log.ell;
        for(int side = 0 ; side < data.cpoly->nb_polys ; side++) {
            sm_side_info_srcptr S = data.sm_info[side];
            if (S->nsm > 0 && (nonvoidside & (((uint64_t) 1) << side))) {
                mpz_poly u;
                mpz_poly_init(u, MAX(1, S->f->deg-1));
                mpz_poly_setcoeff_int64(u, 0, a);
                mpz_poly_setcoeff_int64(u, 1, -b);
                compute_sm_piecewise(u, u, S);
                ASSERT_ALWAYS(u->deg < S->f->deg);
                if (S->mode == SM_MODE_LEGACY_PRE2018) {
                    for(int i = S->f->deg-1-u->deg; i < S->nsm; i++) {
                        mpz_addmul (l, data.log.smlog(side, i), u->coeff[S->f->deg-1-i]);
                    }
                } else {
                    for(int i = 0; i < S->nsm; i++) {
                        mpz_addmul (l, data.log.smlog(side, i), u->coeff[i]);
                    }
                }
                mpz_mod(l, l, ell);
                mpz_poly_clear(u);
            }
        }
    }

    return NULL;
}

/****************** Computation of missing logarithms ************************/
/* Callback function called by filter_rels in compute_log_from_rels */
void *
thread_insert (void * context_data, earlyparsed_relation_ptr rel)
{
  read_data & data = * (read_data *) context_data;
  log_rel_t *lrel = &(data.rels[rel->num]);
  unsigned int next = 0;
  ideal_merge_t buf[REL_MAX_SIZE];
  int c = 0;
  for (unsigned int i = 0; i < rel->nb; i++)
  {
    index_t h = rel->primes[i].h;
    exponent_t e = rel->primes[i].e;

    if (data.log.is_known(h)) {
      mpz_addmul_si (lrel->log_known_part, data.log[h], e);
      c++;
    } else
      buf[next++] = (ideal_merge_t) {.id = h, .e = e};
  }
  if (c) mpz_mod(lrel->log_known_part, lrel->log_known_part, data.log.ell);

  lrel->unknown = idealmerge_my_malloc (next);
  lrel->nb_unknown = next;
  memcpy(lrel->unknown, buf, next * sizeof(ideal_merge_t));

  return NULL;
}

/* Return the number of unknown logarithms in the relation.
 * rels[i].nb_unknown may not be up-to-date (can only be greater than the actual
 * value) */
static inline weight_t
nb_unknown_log (read_data & data, uint64_t i)
{
  log_rel_t * lrel = &(data.rels[i]);
  weight_t j, k, len = lrel->nb_unknown;
  ideal_merge_t *p = lrel->unknown;
  int c = 0;
  for (j = 0, k = 0; k < len; k++)
  {
    pthread_mutex_lock (&lock);
    bool known = data.log.is_known(p[k].id);
    pthread_mutex_unlock (&lock);

    if (!known) {
      if (j != k)
        p[j] = p[k];
      j++;
    } else { // We know this log, add it to log_know_part
      mpz_addmul_si(lrel->log_known_part, data.log[p[k].id], p[k].e);
      c++;
    }
  }
  if (c) mpz_mod(lrel->log_known_part, lrel->log_known_part, data.log.ell);

  lrel->nb_unknown = j;
  return j;
}

/* In a relation with 1 missing logarithm of exponent e, compute its values,
 * i.e. compute   dest <- (-vlog / e) mod ell
 * Return 0 if dest was already known (i.e. computed between the call to
 * nb_unknown_log and the call to this functions), return 1 otherwise.
 */
static inline unsigned int
compute_missing_log (logtab::mpz_rw_accessor dest, mpz_t vlog, int32_t e, mpz_t ell)
{
  unsigned int ret;
  mpz_t tmp;
  mpz_init_set_si (tmp, e);
  mpz_invert (tmp, tmp, ell);
  mpz_neg (vlog, vlog);
  mpz_mul (vlog, vlog, tmp);
  mpz_mod (tmp, vlog, ell);

  pthread_mutex_lock (&lock);
  if (dest.is_known()) {
      // log was already computed by another thread
      ret = 0;
  } else {
      dest = tmp;
      ret = 1;
  }
  pthread_mutex_unlock (&lock);

  mpz_clear (tmp);
  return ret;
}

/* Compute all missing logarithms for relations in [start,end[.
 * Return the number of computed logarithms */
static uint64_t
log_do_one_part_of_iter (read_data & data, bit_vector not_used, uint64_t start,
                         uint64_t end)
{
  uint64_t i, computed = 0;

  for (i = start; i < end; i++)
  {
    if (bit_vector_getbit(not_used, (size_t) i))
    {
      weight_t nb = nb_unknown_log (data, i);
      if (nb <= 1)
      {
        bit_vector_clearbit(not_used, (size_t) i);
        mpz_ptr vlog = data.rels[i].log_known_part;
        if (nb == 0 && mpz_cmp_ui (vlog, 0) != 0)
        {
          gmp_fprintf (stderr, "Error, no unknown log in rel %" PRIu64 " and sum"
                       " of log is not zero, sum is: %Zd\n", i, vlog);
          exit (EXIT_FAILURE);
        }
        else if (nb == 1)
        {
          ideal_merge_t ideal = data.rels[i].unknown[0];
          computed +=
             compute_missing_log (data.log[ideal.id], vlog, ideal.e,
                                  data.log.ell);
        }
      }
    }
  }

  return computed;
}

#define log_do_one_iter_mono(d, bv, n) log_do_one_part_of_iter (d, bv, 0, n)

/************************** Dependency graph *********************************/
typedef struct
{
  uint8_t state;
  uint64_t i;
} node_dep_t;

typedef struct
{
  uint64_t size;
  node_dep_t *tab;
} graph_dep_t;

/* Macro for state of node_dep_t. UNKNOWN must be 0. */
#define NODE_DEP_LOG_UNKNOWN 0
#define NODE_DEP_LOG_KNOWN_FROM_LOGFILE 1
#define NODE_DEP_LOG_RECONSTRUCTED 2

#define GRAPH_DEP_IS_LOG_UNKNOWN(G, h) (G.tab[h].state == NODE_DEP_LOG_UNKNOWN)

graph_dep_t
graph_dep_init (uint64_t size)
{
  node_dep_t *tab = NULL;
  tab = (node_dep_t *) malloc (sizeof(node_dep_t) * size);
  FATAL_ERROR_CHECK(tab == NULL, "Cannot allocate memory");
  memset (tab, 0, sizeof(node_dep_t) * size);
  return (graph_dep_t) {.size = size, .tab = tab};
}

void
graph_dep_clear (graph_dep_t G)
{
  free(G.tab);
}

/* Set G[h].state accordingly to log[h] values */
void
graph_dep_set_log_already_known (graph_dep_t G, logtab const & log)
{
  for (uint64_t h = 0; h < log.nprimes; h++)
  {
    if (log.is_known(h))
      G.tab[h].state = NODE_DEP_LOG_KNOWN_FROM_LOGFILE;
  }
}

uint64_t
graph_dep_needed_rels_from_index (graph_dep_t G, index_t h, light_rels_t rels,
                                  bit_vector needed_rels)
{
  if (G.tab[h].state == NODE_DEP_LOG_UNKNOWN)
  {
    fprintf (stderr, "Error: logarithms of %" PRid" cannot be reconstructed "
                     "from this set of relations. Abort!\n", h);
    abort();
  }
  else if (G.tab[h].state == NODE_DEP_LOG_KNOWN_FROM_LOGFILE)
  {
    /* We know the wanted logarithm from linear algebra, no new relation is
     * necessary */
#if DEBUG >= 1
    fprintf (stderr, "DEBUG: h = %" PRid " is known from logfile\n", h);
#endif
    return 0;
  }
  else
  {
    uint64_t relnum = G.tab[h].i;
    bit_vector_setbit (needed_rels, relnum);
    uint64_t nadded = 1;
    weight_t nb_needed = rels[relnum].len;

#if DEBUG >= 1
    fprintf (stderr, "DEBUG: h = %" PRid " can be reconstructed\n", h);
    fprintf (stderr, "DEBUG:     relation %" PRIu64 " added\n", relnum);
    fprintf (stderr, "DEBUG:     depends of %u others logs\n", nb_needed-1);
#endif

    for (weight_t j = 0; j < nb_needed; j++)
    {
      index_t hh = rels[relnum].needed[j];
      if (!bit_vector_getbit (needed_rels, G.tab[hh].i))
      {
        nadded += graph_dep_needed_rels_from_index (G, hh, rels, needed_rels);
      }
    }
    return nadded;
  }
}

/* Trivial structure containing the data necessary for reading rels for dep. graph */
struct dep_read_data
{
  light_rels_t rels;
  graph_dep_t G;
  dep_read_data(light_rels_t rels, graph_dep_t G) : rels(rels), G(G) {}
};

/* Callback function called by filter_rels in compute_needed_rels */
void *
dep_thread_insert (void * context_data, earlyparsed_relation_ptr rel)
{
  dep_read_data & data = * (dep_read_data *) context_data;
  light_rels_t lrel = &(data.rels[rel->num]);
  unsigned int next = 0;
  index_t buf[REL_MAX_SIZE];

  for (unsigned int i = 0; i < rel->nb; i++)
  {
    index_t h = rel->primes[i].h;
    if (GRAPH_DEP_IS_LOG_UNKNOWN(data.G, h))
      buf[next++] = h;
  }

  lrel->needed = index_my_malloc (next);
  lrel->len = next;
  memcpy(lrel->needed, buf, next * sizeof(ideal_merge_t));

  return NULL;
}

/* Return the number of unknown logarithms in the relation and put in h the last
 * unknown logarithm of the relations (useful when the number of unknown
 * logarithms is 1) */
static inline weight_t
dep_nb_unknown_log (dep_read_data & data, uint64_t i, index_t *h)
{
  weight_t k, nb;
  index_t *p = data.rels[i].needed;

  for (nb = 0, k = 0; k < data.rels[i].len; k++)
  {
    pthread_mutex_lock (&lock);
    int unknow = GRAPH_DEP_IS_LOG_UNKNOWN (data.G, p[k]);
    pthread_mutex_unlock (&lock);

    if (unknow) // we do not know the log if this ideal
    {
      nb++;
      *h = p[k];
    }
  }
  return nb;
}

/* Compute all dependencies for relations in [start,end[.
 * Return the number of dependencies found */
static uint64_t
dep_do_one_part_of_iter (dep_read_data & data, bit_vector not_used,
                         uint64_t start, uint64_t end)
{
  uint64_t i, computed = 0;

  for (i = start; i < end; i++)
  {
    if (bit_vector_getbit(not_used, (size_t) i))
    {
      index_t h = 0; // Placate gcc
      weight_t nb = dep_nb_unknown_log(data, i, &h);
      if (nb <= 1)
      {
        bit_vector_clearbit(not_used, (size_t) i);
        if (nb == 1)
        {
          data.G.tab[h].state = NODE_DEP_LOG_RECONSTRUCTED;
          data.G.tab[h].i = i;
          computed++;
        }
      }
    }
  }

  return computed;
}

#define dep_do_one_iter_mono(d, bv, n) dep_do_one_part_of_iter (d, bv, 0, n)

/******************** Code for multi thread version **************************/
#define SIZE_BLOCK 1024

typedef struct {
  void *data;
  bit_vector_ptr not_used;
  uint64_t offset;
  uint64_t nb;
  uint64_t computed;
  int version; /* 0 means call dep_* , 1 means call log_* */
} thread_info;

void * thread_start(void *arg)
{
  thread_info *ti = (thread_info *) arg;
  bit_vector_ptr not_used = ti->not_used;
  uint64_t start = ti->offset;
  uint64_t end = start + ti->nb;

  if (ti->version == 0)
  {
    dep_read_data & data = * (dep_read_data *) ti->data;
    ti->computed = dep_do_one_part_of_iter (data, not_used, start, end);
  }
  else
  {
    read_data & data = * (read_data *) ti->data;
    ti->computed = log_do_one_part_of_iter (data, not_used, start, end);
  }
  return NULL;
}

static uint64_t
do_one_iter_mt (void *data, bit_vector not_used, int nt, uint64_t nrels,
                int version)
{
  // We'll use a rotating buffer of thread id.
  pthread_t *threads;
  threads = (pthread_t *) malloc( nt * sizeof(pthread_t));
  int active_threads = 0;  // number of running threads
  int threads_head = 0;    // next thread to wait / restart.

  // Prepare the main loop
  uint64_t i = 0; // counter of relation.
  uint64_t computed = 0;

  // Arguments for threads
  thread_info *tis;
  tis = (thread_info *) malloc( nt * sizeof(thread_info));
  for (int i = 0; i < nt; ++i)
  {
    tis[i].data = data;
    tis[i].not_used = not_used;
    tis[i].version = version;
    // offset and nb must be adjusted.
  }

  // Main loop
  while ((i < nrels) || (active_threads > 0))
  {
    // Start / restart as many threads as allowed
    if ((active_threads < nt) && (i < nrels))
    {
      tis[threads_head].offset = i;
      tis[threads_head].nb = MIN(SIZE_BLOCK, nrels-i);
      pthread_create(&threads[threads_head], NULL,
          &thread_start, (void *)(&tis[threads_head]));
      i += SIZE_BLOCK;
      active_threads++;
      threads_head++;
      if (threads_head == nt)
        threads_head = 0;
      continue;
    }
    // Wait for the next thread to finish in order to print result.
    pthread_join(threads[threads_head], NULL);
    active_threads--;
    computed += tis[threads_head].computed;

    // If we are at the end, no job will be restarted, but head still
    // must be incremented.
    if (i >= nrels)
    {
      threads_head++;
      if (threads_head == nt)
        threads_head = 0;
    }
  }

  free(tis);
  free(threads);
  return computed;
}

/* Compute all missing logarithms possible. Run through all the relations once.
 * Multi thread version
 * Return the number of computed logarithms */
static inline uint64_t
log_do_one_iter_mt (read_data & d, bit_vector bv, int nt, uint64_t nrels)
{
  return do_one_iter_mt ((void *) & d, bv, nt, nrels, 1);
}

/* Compute all missing logarithms possible. Run through all the relations once.
 * Multi thread version
 * Return the number of computed logarithms */
static inline uint64_t
dep_do_one_iter_mt (dep_read_data & d, bit_vector bv, int nt, uint64_t nrels)
{
  return do_one_iter_mt ((void *) & d, bv, nt, nrels, 0);
}

/***************** Important functions called by main ************************/
/* Read the logarithms computed by the linear algebra */
static void
read_log_format_LA (logtab & log, const char *logfile, const char *idealsfile,
                    sm_side_info *sm_info, int nb_polys)
{
  uint64_t i, ncols, col;
  index_t h;
  mpz_t tmp_log;
  FILE *flog = NULL, *fid = NULL;

  printf ("# Reading logarithms in LA format from %s\n", logfile);
  printf ("# Reading links between matrix columns and ideals from %s\n",
                                                                idealsfile);
  fflush(stdout);
  flog = fopen_maybe_compressed (logfile, "r");
  FATAL_ERROR_CHECK(flog == NULL, "Cannot open file for reading logarithms");
  fid = fopen_maybe_compressed (idealsfile, "r");
  FATAL_ERROR_CHECK(fid == NULL, "Cannot open ideals file");

  if (fscanf (fid, "# %" SCNu64 "\n", &ncols) != 1)
  {
    fprintf(stderr, "Error while reading first line of %s\n", idealsfile);
    abort();
  }

  mpz_init (tmp_log);
  i = 0;
  stats_init (stats, stdout, &i, nbits(ncols)-5, "Read", "logarithms", "", "logs");
  while (fscanf (fid, "%" SCNu64 " %" PRid "\n", &col, &h) == 2)
  {
    FATAL_ERROR_CHECK (col >= ncols, "Too big value of column number");
    FATAL_ERROR_CHECK (h >= log.nprimes, "Too big value of index");

    int ret = gmp_fscanf (flog, "%Zd\n", tmp_log);
    FATAL_ERROR_CHECK (ret != 1, "Error in file containing logarithms values");

    ASSERT_ALWAYS (col == i);
    log[h] = tmp_log;
    i++;
    if (stats_test_progress (stats))
      stats_print_progress (stats, i, 0, 0, 0);
  }
  stats_print_progress (stats, i, 0, 0, 1);
  ASSERT_ALWAYS (feof(fid));
  ASSERT_ALWAYS (i == ncols);

  for(int side = 0; side < nb_polys; side++)
  {
      for (int ism = 0; ism < sm_info[side]->nsm; ism++)
      {
        int ret = gmp_fscanf (flog, "%Zd\n", tmp_log);
        FATAL_ERROR_CHECK (ret != 1, "Error in file containing logarithms values");
        log.smlog(side, ism) = tmp_log;
      }
  }
  /* If we are not at the end of the file, it means that it remains some values
   * and we do not know to what "ideals" they correspond. Probably an error
   * somewhere, it is better to abort. */
  ASSERT_ALWAYS (feof(flog));

  for(int side = 0; side < nb_polys; side++)
  {
    if (sm_info[side]->nsm)
      printf ("# Logarithms for %d SM columns on side %d were also read\n",
              sm_info[side]->nsm, side);
  }
  mpz_clear (tmp_log);
  fclose_maybe_compressed (flog, logfile);
  fclose_maybe_compressed (fid, idealsfile);
}

/* Read the logarithms in output format of reconstructlog */
static void
read_log_format_reconstruct (logtab & log, MAYBE_UNUSED renumber_t const & renumb,
                             const char *filename)
{
  uint64_t nread = 0;
  index_t h;
  mpz_t tmp_log;
  FILE *f = NULL;
  int ret;

  printf ("# Reading logarithms in reconstruct format from %s\n", filename);
  fflush(stdout);
  f = fopen_maybe_compressed (filename, "r");
  FATAL_ERROR_CHECK(f == NULL, "Cannot open file for reading");

  mpz_init (tmp_log);
  stats_init (stats, stdout, &nread, nbits(renumb.get_size())-5, "Read", "logarithms", "",
              "logs");
  for (index_t i = 0; i < renumb.number_of_additional_columns(); i++) {
      ret = gmp_fscanf (f, "%" SCNid " added column %Zd\n", &h, tmp_log);
      ASSERT_ALWAYS (ret == 2);
      ASSERT_ALWAYS (renumb.is_additional_column(h));
      nread++;
      log[h] = tmp_log;
  }
  for (index_t i = 0; i < renumb.number_of_bad_ideals(); i++) {
      ret = gmp_fscanf (f, "%" SCNid " bad ideals %Zd\n", &h, tmp_log);
      ASSERT_ALWAYS (ret == 2);
      ASSERT_ALWAYS (renumb.is_bad(h));
      nread++;
      log[h] = tmp_log;
  }
  while (gmp_fscanf (f, "%" SCNid " %*" SCNpr " %*d %*s %Zd\n", &h, tmp_log)
          == 2)
  {
    nread++;
    log[h] = tmp_log;
    if (stats_test_progress (stats))
      stats_print_progress (stats, nread, 0, 0, 0);
  }
  stats_print_progress (stats, nread, 0, 0, 1);

  for (int nsm = 0; nsm < log.nbsm; nsm++)
  {
    unsigned int n, side;
    if (nsm == 0) /* h was already read by previous gmp_fscanf */
    {
      ret = gmp_fscanf (f, "SM %u %u %Zd\n", &side, &n, tmp_log);
      ASSERT_ALWAYS (ret == 3);
    }
    else
    {
      ret = gmp_fscanf (f, "%" SCNid " SM %u %u %Zd\n", &h, &side, &n, tmp_log);
      ASSERT_ALWAYS (ret == 4);
    }
    //    ASSERT_ALWAYS (n == nsm); // obsolete with new coding
    ASSERT_ALWAYS (h == (index_t) nsm + log.nprimes);
    log[h] = tmp_log;
  }
  ASSERT_ALWAYS (feof(f));

  mpz_clear (tmp_log);
  fclose_maybe_compressed (f, filename);
}

/* Write values of the known logarithms. */
/* TODO: use the fact that we now have a way to iterate through the
 * renumber table!!! */
static void
write_log (const char *filename, logtab & log, renumber_t const & tab, 
	   sm_side_info *sm_info)
{
  uint64_t i;
  FILE *f = NULL;

  printf ("# Opening %s for writing logarithms\n", filename);
  fflush(stdout);
  f = fopen_maybe_compressed (filename, "w");
  FATAL_ERROR_CHECK(f == NULL, "Cannot open file for writing");

  /* Divide all known logs by 'base' so that the first known non-zero logarithm
   * is equal to 1.
   * TODO: make a command line argument to choose this 'base'.
   */
  int base_already_set = 0;
  cxx_mpz base, scaled;
  for (i = 0; i < log.nprimes + log.nbsm; i++)
  {
    if (!log.is_known(i)) continue;
    if (log.is_zero(i)) continue;

      if (!base_already_set)
      {
        base_already_set = 1;
        /* base = 1/log[i] mod ell */
        int ret = mpz_invert (base, log[i], log.ell);
        ASSERT_ALWAYS (ret != 0);
        mpz_set_ui(scaled, 1);
        log.force_set(i, scaled);
      }
      else
      {
        mpz_mul (scaled, log[i], base);
        mpz_mod (scaled, scaled, log.ell);
        log.force_set(i, scaled);
      }
  }

  uint64_t nknown = 0;
  stats_init (stats, stdout, &nknown, nbits(tab.get_size())-5, "Wrote",
              "known logarithms", "ideals", "logs");
  i = 0;
  for( ; tab.is_additional_column(i) ; i++) {
      if (!log.is_known(i)) continue;
      nknown++;
      gmp_fprintf (f, "%" PRid " added column %Zd\n", 
              (index_t) i, (mpz_srcptr) log[i]);
  }
  for( ; tab.is_bad(i) ; i++) {
      if (!log.is_known(i)) continue;
      nknown++;
      gmp_fprintf (f, "%" PRid " bad ideals %Zd\n",
              (index_t) i, (mpz_srcptr) log[i]);
  }
  for( ; i < tab.get_size(); i++) {
      if (!log.is_known(i)) continue;
      nknown++;
      // TODO forward iterator
      renumber_t::p_r_side x = tab.p_r_from_index(i);
      if (x.side != tab.get_rational_side())
          gmp_fprintf (f, "%" PRid " %" PRpr " %d %" PRpr " %Zd\n",
                  (index_t) i, x.p, x.side, x.r, (mpz_srcptr) log[i]);
      else
          gmp_fprintf (f, "%" PRid " %" PRpr " %d rat %Zd\n",
                  (index_t) i, x.p, x.side, (mpz_srcptr) log[i]);
      if (stats_test_progress (stats))
          stats_print_progress (stats, nknown, i+1, 0, 0);
  }
  stats_print_progress (stats, nknown, tab.get_size(), 0, 1);
  for (int nsm = 0, i = tab.get_size(); nsm < log.nbsm; nsm++)
  {
    // compute side
    int side, nsm_tot = sm_info[0]->nsm, jnsm = nsm;
    for(side = 0; ((int)nsm) >= nsm_tot; side++){
	nsm_tot += sm_info[side+1]->nsm;
	jnsm -= sm_info[side]->nsm;
    }
    ASSERT_ALWAYS ((jnsm >= 0) && (jnsm < sm_info[side]->nsm));
    if (log.is_zero(i+nsm)) {
        printf("# Note: on side %d, log of SM number %d is zero\n", side, jnsm);
    } else {
        ASSERT_ALWAYS (log.is_known(i+nsm));
    }
    gmp_fprintf (f, "%" PRid " SM %d %d %Zd\n",
            (index_t) (i+nsm), side, jnsm,
            (mpz_srcptr) log[i+nsm]);
  }

  uint64_t missing = tab.get_size() - nknown;
  printf ("# factor base contains %" PRIu64 " elements\n"
          "# logarithms of %" PRIu64 " elements are known (%.1f%%)\n"
          "# logarithms of %" PRIu64 " elements are missing (%.1f%%)\n",
          tab.get_size(), nknown, 100.0 * nknown / (double) tab.get_size(),
          missing, 100.0 * missing / (double) tab.get_size());
  fclose_maybe_compressed (f, filename);
  ASSERT_ALWAYS (log.nknown == nknown);
}

/* Given a filename, compute all the possible logarithms of ideals appearing in
 * the file. Return the number of computed logarithms.
 * Modify its first argument bit_vector needed_rels */
static uint64_t
compute_log_from_rels (bit_vector needed_rels,
                       const char *relspfilename, uint64_t nrels_purged,
                       const char *relsdfilename, uint64_t nrels_del,
                       uint64_t nrels_needed, int nt,
                       read_data & data)
{
  double wct_tt0, wct_tt;
  uint64_t total_computed = 0, iter = 0, computed;
  uint64_t nrels = nrels_purged + nrels_del;
  ASSERT_ALWAYS (nrels_needed > 0);

  /* Reading all relations */
  printf ("# Reading relations from %s and %s\n", relspfilename, relsdfilename);
  if (nrels_needed != nrels)
    printf ("# Parsing only %" PRIu64 " needed relations out of %" PRIu64 "\n",
            nrels_needed, nrels);
#if DEBUG >= 1
  printf ("# DEBUG: Using %d thread(s) for thread_sm\n", nt);
#endif
  fflush(stdout);
  char *fic[3] = {(char *) relspfilename, (char *) relsdfilename, NULL};

  /* When purged.gz and relsdel.gz both have SM info included, we may
   * have an advantage in having more threads for thread_insert. Note
   * though that we'll most probably be limited by gzip throughput */

  int ni = 1;
  int ns = nt;
  
  if (!filename_matches_one_compression_format(relspfilename) && !filename_matches_one_compression_format(relsdfilename)) {
      printf("# Files %s and %s are uncompressed, limiting consumer threads\n",
              relspfilename,
              relsdfilename);
      ns = 4;
  }
  struct filter_rels_description desc[3] = {
                   { .f = thread_insert, .arg=(void*) &data, .n=ni},
                   { .f = thread_sm,     .arg=(void*) &data, .n=ns},
                   { .f = NULL,          .arg=0,    .n=0}
      };
  filter_rels2 (fic, desc,
          EARLYPARSE_NEED_AB_HEXA |
          EARLYPARSE_NEED_INDEX |
          EARLYPARSE_NEED_SM, /* It's fine (albeit slow) if we recompute them */
          needed_rels, NULL);

  /* computing missing log */
  printf ("# Starting to compute missing logarithms from rels\n");

  /* adjust the number of threads based on the number of needed relations */
  double ntm = ceil((nrels_needed + 0.0)/SIZE_BLOCK);
  if (nt > ntm)
    nt = (int) ntm;

  if (nt > 1)
    printf("# Using multithread version with %d threads\n", nt);
  else
    printf("# Using monothread version\n");

  wct_tt0 = wct_seconds();
  do
  {
    printf ("# Iteration %" PRIu64 ": starting... [%" PRIu64 " known logs]\n", iter, data.log.nknown);
    fflush(stdout);
    wct_tt = wct_seconds();

    if (nt > 1)
      computed = log_do_one_iter_mt (data, needed_rels, nt, nrels);
    else
      computed = log_do_one_iter_mono (data, needed_rels, nrels);
    total_computed += computed;

    printf ("# Iteration %" PRIu64 ": %" PRIu64 " new logarithms computed\n",
            iter, computed);
    printf ("# Iteration %" PRIu64 " took %.1fs (wall-clock time).\n",
            iter, wct_seconds() - wct_tt);

    iter++;
  } while (computed);

  printf ("# Computing %" PRIu64 " new logarithms took %.1fs (wall-clock "
          "time)\n", total_computed, wct_seconds() - wct_tt0);

  size_t c = bit_vector_popcount(needed_rels);
  if (c != 0)
    fprintf(stderr, "### Warning, %zu relations were not used\n", c);

  return total_computed;
}

/* Given a filename, compute all the relations needed to compute the logarithms
 * appearing in the file.
 * needed_rels should be initialized before calling this function. Its size
 * must be nrels_purged + nrels_del.
 * Output:
 *    bit_vector needed_rels, where bits of needed rels are set.
 *    Return the number of needed_rels.*/
static uint64_t
compute_needed_rels (bit_vector needed_rels,
                     const char *relspfilename, uint64_t nrels_purged,
                     const char *relsdfilename, uint64_t nrels_del,
                     logtab & log, const char *wanted_filename, int nt)
{
  double wct_tt0, wct_tt;
  // uint64_t total_computed = 0;
  uint64_t iter = 0, computed;
  uint64_t nrels = nrels_purged + nrels_del;
  graph_dep_t dep_graph = graph_dep_init (log.nprimes);
  light_rels_t rels = light_rels_init (nrels);

  graph_dep_set_log_already_known (dep_graph, log);

  dep_read_data data(rels, dep_graph);

  /* Init bit_vector to remember which relations were already used */
  bit_vector_set (needed_rels, 1);

  /* Reading all relations */
  printf ("# Reading relations from %s and %s\n", relspfilename, relsdfilename);
  fflush(stdout);
  char *fic[3] = {(char *) relspfilename, (char *) relsdfilename, NULL};
  filter_rels (fic, (filter_rels_callback_t) &dep_thread_insert, (void *) &data,
               EARLYPARSE_NEED_INDEX, NULL, NULL);

  /* computing dependencies */
  printf ("# Starting to compute dependencies from rels\n");

  /* adjust the number of threads based on the number of relations */
  double ntm = ceil((nrels + 0.0)/SIZE_BLOCK);
  if (nt > ntm)
    nt = (int) ntm;

  if (nt > 1)
    printf("# Using multithread version with %d threads\n", nt);
  else
    printf("# Using monothread version\n");

  wct_tt0 = wct_seconds();
  do
  {
    printf ("# Iteration %" PRIu64 ": starting...\n", iter);
    fflush(stdout);
    wct_tt = wct_seconds();

    if (nt > 1)
      computed = dep_do_one_iter_mt (data, needed_rels, nt, nrels);
    else
      computed = dep_do_one_iter_mono (data, needed_rels, nrels);
    // total_computed += computed;

    printf ("# Iteration %" PRIu64 ": %" PRIu64 " new dependencies computed\n",
            iter, computed);
    printf ("# Iteration %" PRIu64 " took %.1fs (wall-clock time).\n",
            iter, wct_seconds() - wct_tt);

    iter++;
  } while (computed);

  printf ("# Computing dependencies took %.1fs (wall-clock time)\n",
          wct_seconds() - wct_tt0);

  FILE *f = NULL;
  printf ("# Reading wanted logarithms from %s\n", wanted_filename);
  fflush(stdout);
  f = fopen_maybe_compressed (wanted_filename, "r");
  FATAL_ERROR_CHECK(f == NULL, "Cannot open file for reading");

  bit_vector_set (needed_rels, 0);
  index_t h;
  uint64_t nadded, nrels_necessary = 0, nwanted_log = 0;
  wct_tt = wct_seconds();
  while (fscanf (f, "%" SCNid "\n", &h) == 1)
  {
    FATAL_ERROR_CHECK (h >= log.nprimes, "Too big value of index");
    printf ("# Computing rels necessary for wanted log %" PRid "\n", h);
    fflush(stdout);
    nadded = graph_dep_needed_rels_from_index (dep_graph, h, rels, needed_rels);
    nrels_necessary += nadded;
    printf ("-> %" PRIu64 " needed relations were added (%" PRIu64 " so far)\n",
            nadded, nrels_necessary);
    nwanted_log++;
  }

  fclose_maybe_compressed (f, wanted_filename);
  printf ("# Reading %" PRIu64 " wanted logarithms took %.1fs\n", nwanted_log,
          wct_seconds() - wct_tt);
  printf ("# %" PRIu64 " relations are needed to compute these logarithms\n",
          nrels_necessary);
  ASSERT_ALWAYS (nrels_necessary == bit_vector_popcount (needed_rels));

  light_rels_free (rels);
  graph_dep_clear (dep_graph);
  return nrels_necessary;
}

/********************* usage functions and main ******************************/
static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "log", "input file containing known logarithms");
  param_list_decl_usage(pl, "logformat", "format of input log file: 'LA' or "
                                         "'reconstruct' (default is 'LA')");
  param_list_decl_usage(pl, "ell", "group order (see sm -ell parameter)");
  param_list_decl_usage(pl, "out", "output file for logarithms");
  param_list_decl_usage(pl, "renumber", "input file for renumbering table");
  param_list_decl_usage(pl, "poly", "input polynomial file");
  param_list_decl_usage(pl, "ideals", "link between matrix cols and ideals "
                                      "(see replay -ideals parameter)");
  param_list_decl_usage(pl, "purged", "file with purged relations "
                                      "(see purge -out parameter)");
  param_list_decl_usage(pl, "relsdel", "file with relations deleted by purge "
                                      "(see purge -outdel parameter)");
  param_list_decl_usage(pl, "nrels", "number of relations (same as purge "
                                     "-nrels parameter)");
  param_list_decl_usage(pl, "partial", "do not reconstruct everything "
                                       "that can be reconstructed");
  param_list_decl_usage(pl, "sm-mode", "SM mode (see sm-portability.h)");
  param_list_decl_usage(pl, "nsm", "number of SM's to add on side 0,1,...");
  param_list_decl_usage(pl, "mt", "number of threads (default 1)");
  param_list_decl_usage(pl, "wanted", "file containing list of wanted logs");
  param_list_decl_usage(pl, "force-posix-threads", "force the use of posix threads, do not rely on platform memory semantics");
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
  verbose_decl_usage(pl);
}

static void
usage (param_list pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}


// coverity[root_function]
int
main(int argc, char *argv[])
{
  char *argv0 = argv[0];

  uint64_t nrels_tot = 0, nrels_purged, nrels_del, nrels_needed;
  uint64_t nprimes;
  int mt = 1;
  int partial = 0;

  mpz_t ell;
  cxx_cado_poly cpoly;

  param_list pl;
  param_list_init(pl);
  declare_usage(pl);
  argv++,argc--;

  param_list_configure_switch(pl, "partial", &partial);
  param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);

#ifdef HAVE_MINGW
  _fmode = _O_BINARY;     /* Binary open for all files */
#endif

  if (argc == 0)
    usage(pl, argv0);

  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
    fprintf (stderr, "Unknown option: %s\n", argv[0]);
    usage(pl, argv0);
  }
  /* print command-line arguments */
  verbose_interpret_parameters(pl);
  param_list_print_command_line (stdout, pl);
  fflush(stdout);

  mpz_init (ell);
  const char * logfilename = param_list_lookup_string(pl, "log");
  const char * logformat = param_list_lookup_string(pl, "logformat");
  const char * idealsfilename = param_list_lookup_string(pl, "ideals");
  const char * relsdfilename = param_list_lookup_string(pl, "relsdel");
  const char * relspfilename = param_list_lookup_string(pl, "purged");
  const char * outfilename = param_list_lookup_string(pl, "out");
  const char * renumberfilename = param_list_lookup_string(pl, "renumber");
  const char * polyfilename = param_list_lookup_string(pl, "poly");
  const char * wantedfilename = param_list_lookup_string(pl, "wanted");
  param_list_parse_int(pl, "mt", &mt);
  const char *path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");

  /* Some checks on command line arguments */
  if (!param_list_parse_mpz(pl, "ell", ell) || mpz_cmp_ui (ell, 0) <= 0)
  {
    fprintf(stderr, "Error, missing -ell command line argument "
                    "(or ell <= 0)\n");
    usage (pl, argv0);
  }
  if (!param_list_parse_uint64(pl, "nrels", &nrels_tot) || nrels_tot == 0)
  {
    fprintf(stderr, "Error, missing -nrels command line argument "
                    "(or nrels = 0)\n");
    usage (pl, argv0);
  }
  if (logfilename == NULL)
  {
    fprintf(stderr, "Error, missing -log command line argument\n");
    usage (pl, argv0);
  }
  if (relspfilename == NULL)
  {
    fprintf(stderr, "Error, missing -purged command line argument\n");
    usage (pl, argv0);
  }
  if (relsdfilename == NULL)
  {
    fprintf(stderr, "Error, missing -relsdel command line argument\n");
    usage (pl, argv0);
  }
  if (outfilename == NULL)
  {
    fprintf(stderr, "Error, missing -out command line argument\n");
    usage (pl, argv0);
  }
  if (renumberfilename == NULL)
  {
    fprintf(stderr, "Error, missing -renumber command line argument\n");
    usage (pl, argv0);
  }
  if (polyfilename == NULL)
  {
    fprintf(stderr, "Error, missing -poly command line argument\n");
    usage (pl, argv0);
  }
  if (mt < 1)
  {
    fprintf(stderr, "Error: parameter mt must be at least 1\n");
    usage (pl, argv0);
  }

  if (logformat != NULL)
  {
    if (strcmp(logformat, "LA") != 0 && strcmp(logformat, "reconstruct") != 0)
    {
      fprintf(stderr, "Error, unknown -formatlog argument. Must be 'LA' or "
                      "'reconstruct'\n");
      usage (pl, argv0);
    }
  }
  if ((logformat == NULL || strcmp(logformat, "LA") == 0) &&
      idealsfilename == NULL)
  {
    fprintf(stderr, "Error, missing -ideals command line argument\n");
    usage (pl, argv0);
  }

  if (wantedfilename != NULL && !partial)
  {
    fprintf(stderr, "Warning, -wanted command line argument is ignored if "
                    "-partial is not set\n");
  }

  if (!cado_poly_read (cpoly, polyfilename))
  {
    fprintf (stderr, "Error reading polynomial file\n");
    exit (EXIT_FAILURE);
  }

  /* Read number of sm to be printed from command line */
  std::vector<int> nsm_arg(cpoly->nb_polys, -1);
  param_list_parse_int_args_per_side(pl, "nsm",
          nsm_arg.data(), cpoly->nb_polys,
          ARGS_PER_SIDE_DEFAULT_AS_IS);

  for(int side = 0; side < cpoly->nb_polys; side++) {
      if (nsm_arg[side] < 0)
          continue;
      if (nsm_arg[side] > cpoly->pols[side]->deg) {
          fprintf(stderr, "Error: nsm%d=%d can not exceed the degree=%d\n",
                  side, nsm_arg[side], cpoly->pols[side]->deg);
          exit (EXIT_FAILURE);
      }
  }

  const char * sm_mode_string = param_list_lookup_string(pl, "sm-mode");

  if (param_list_warn_unused(pl))
  {
    fprintf(stderr, "Error, unused parameters are given\n");
    usage(pl, argv0);
  }

  set_antebuffer_path (argv0, path_antebuffer);

  /* Init data for computation of the SMs. */
  sm_side_info * sm_info = new sm_side_info[cpoly->nb_polys];
  for (int side = 0; side < cpoly->nb_polys; side++)
  {
    sm_side_info_init(sm_info[side], cpoly->pols[side], ell);
    sm_side_info_set_mode(sm_info[side], sm_mode_string);
    fprintf(stdout, "\n# Polynomial on side %d:\n# F[%d] = ", side, side);
    mpz_poly_fprintf(stdout, cpoly->pols[side]);
    printf("# SM info on side %d:\n", side);
    sm_side_info_print(stdout, sm_info[side]);
    if (nsm_arg[side] >= 0)
      sm_info[side]->nsm = nsm_arg[side]; /* command line wins */
    printf("# Will use %d SMs on side %d\n", sm_info[side]->nsm, side);

    /* do some consistency checks */
    if (sm_info[side]->unit_rank != sm_info[side]->nsm)
    {
      fprintf(stderr, "# On side %d, unit rank is %d, computing %d SMs ; "
                      "weird.\n", side, sm_info[side]->unit_rank,
                      sm_info[side]->nsm);
      /* for the 0 case, we haven't computed anything: prevent the
       * user from asking SM data anyway */
      ASSERT_ALWAYS(sm_info[side]->unit_rank != 0);
    }
  }
  fflush(stdout);


  /* Reading renumber file */
  /* XXX legacy format insists on getting the badidealinfo file */
  printf ("\n###### Reading renumber file ######\n");
  renumber_t renumber_table(cpoly);
  renumber_table.read_from_file(renumberfilename, 1);
  nprimes = renumber_table.get_size();

  /* Read number of rows and cols on first line of purged file */
  {
      uint64_t nideals_purged;
      purgedfile_read_firstline (relspfilename, &nrels_purged, &nideals_purged);
      nrels_del = nrels_tot - nrels_purged;
  }

  /* Malloc'ing log tab and reading values of log */
  printf ("\n###### Reading known logarithms ######\n");
  fflush(stdout);

  logtab log(cpoly, sm_info, nprimes, ell);

  if (logformat == NULL || strcmp(logformat, "LA") == 0)
      read_log_format_LA (log, logfilename, idealsfilename, sm_info,
                                                            cpoly->nb_polys);
  else
    read_log_format_reconstruct (log, renumber_table, logfilename);

  /* Init bit_vector of rels that must be process by compute_log_from_rels */
  bit_vector rels_to_process;
  bit_vector_init (rels_to_process, nrels_tot);

  if (partial)
  {
    if (wantedfilename == NULL)
    {
      bit_vector_set (rels_to_process, 0);
      nrels_needed = 0;
    }
    else /* We compute needed rels for logarithms in wantedfilename */
    {
      printf ("\n###### Computing needed rels ######\n");
      nrels_needed =
        compute_needed_rels (rels_to_process, relspfilename, nrels_purged,
                             relsdfilename, nrels_del, log, wantedfilename, mt);
    }
  }
  else
  {
    bit_vector_set (rels_to_process, 1);
    nrels_needed = nrels_tot;
  }

  /* Computing logs using rels in purged file */
  printf ("\n###### Computing logarithms using rels ######\n");
  if (nrels_needed > 0)
  {
      read_data data(log, nrels_purged + nrels_del, cpoly, sm_info,
                      renumber_table);

      compute_log_from_rels (rels_to_process, relspfilename,
              nrels_purged, relsdfilename,
              nrels_del, nrels_needed, mt,
              data);
      printf ("# %" PRIu64 " logarithms are known.\n", log.nknown);
      extern double m_seconds;
      fprintf(stderr, "# %.2f\n", m_seconds);
  }
  else
    printf ("# All wanted logarithms are already known, skipping this step\n");
  fflush(stdout);

  /* Writing all the logs in outfile */
  printf ("\n###### Writing logarithms in a file ######\n");
  write_log (outfilename, log, renumber_table, sm_info);

  /* freeing and closing */
  mpz_clear(ell);

  for (int side = 0 ; side < cpoly->nb_polys ; side++)
    sm_side_info_clear (sm_info[side]);
  delete[] sm_info;

  bit_vector_clear(rels_to_process);
  param_list_clear (pl);
  return EXIT_SUCCESS;
}
