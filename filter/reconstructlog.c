#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>   /* for _O_BINARY */

#include "filter_config.h"

#include "utils_with_io.h"
#ifdef FOR_FFS
#include "utils-ffs.h"
#endif
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

/***** Light relations structure (for constructing the dependancy graph) *****/
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
struct logtab_struct
{
  uint64_t nprimes;
  uint64_t nknown;
  unsigned int nbsm;
  mpz_ptr ell;
  mpz_t *tab;
};

typedef struct logtab_struct logtab_t[1];
typedef struct logtab_struct * logtab_ptr;

static void
logtab_init (logtab_ptr log, uint64_t nprimes, int nbsm, mpz_t ell)
{
  log->tab = NULL;
  uint64_t size = nprimes + (uint64_t)nbsm;
  size_t ell_nbits = mpz_size(ell) * GMP_LIMB_BITS;

  log->nprimes = nprimes;
  log->nknown = 0;
  log->nbsm = nbsm;
  log->ell = ell;
  log->tab = (mpz_t *) malloc (size * sizeof(mpz_t));
  FATAL_ERROR_CHECK(log->tab == NULL, "Cannot allocate memory");
  for (uint64_t i = 0; i < size; i++)
  {
    mpz_init2 (log->tab[i], ell_nbits);
    mpz_set_si(log->tab[i], -1);
  }
}

static void
logtab_insert (logtab_ptr log, index_t h, mpz_t logvalue)
{
  if (mpz_cmp_ui (logvalue, 0) < 0)
  {
    fprintf (stderr, "Warning, log is negative for h = %" PRid "\n", h);
    mpz_mod (logvalue, logvalue, log->ell);
  }
  else if (mpz_cmp (logvalue, log->ell) >= 0)
  {
    fprintf (stderr, "Warning, log >= ell for h = %" PRid "\n", h);
    mpz_mod (logvalue, logvalue, log->ell);
  }
  if (mpz_cmp_ui (logvalue, 0) == 0)
    fprintf (stderr, "Warning, log is zero for h = %" PRid "\n", h);

  if (mpz_sgn (log->tab[h]) >= 0) //already known
  {
    ASSERT_ALWAYS (mpz_cmp (log->tab[h], logvalue) == 0);
  }
  else
  {
    mpz_set (log->tab[h], logvalue);
    if (h < log->nprimes) // log of SM columns are not taken into account
      log->nknown++;
  }
}

static void
logtab_clear (logtab_ptr log)
{
  uint64_t size = log->nprimes + log->nbsm;
  for (uint64_t i = 0; i < size; i++)
    mpz_clear (log->tab[i]);
  free (log->tab);
}

/************ Struct used for reading rels files with process_rels ***********/
struct read_data_s
{
  log_rel_t *rels;
  uint64_t nrels;
  logtab_ptr log;
  sm_side_info *sm_info;
  cado_poly_ptr poly;
  renumber_ptr renum_tab;
  const char * abunitsdirname[NB_POLYS_MAX];
  mpz_t * smlogs[NB_POLYS_MAX];    /* the known logarithms of the SMs */
};

typedef struct read_data_s read_data_t[1];
typedef struct read_data_s * read_data_ptr;

/* Init a read_data_t structure for nrels rels and nprimes primes. Assume the
 * table of log is already allocated, but not the table of log_rel_t     */
static void
read_data_init (read_data_ptr data, logtab_ptr log, uint64_t nrels,
                cado_poly_ptr poly, sm_side_info *sm_info,
                renumber_ptr renum_tab, const char *abunits0dirname,
                const char *abunits1dirname)
{
  data->rels = log_rel_init (nrels);
  data->log = log;
  data->abunitsdirname[0] = abunits0dirname;
  data->abunitsdirname[1] = abunits1dirname;
  data->poly = poly;
  data->renum_tab = renum_tab;
  data->nrels = nrels;
  data->sm_info = sm_info;
  for (int side = 0; side < poly->nb_polys; side++)
  {
    if (side == 0)
      data->smlogs[0] = &(log->tab[log->nprimes]);
    else
      data->smlogs[side] = data->smlogs[side-1] + data->sm_info[side-1]->nsm;
  }
  fflush(stdout);
}

static void
read_data_free (read_data_ptr data)
{
  log_rel_free (data->rels, data->nrels);
}

/************************ Handling of the SMs *******************************/
/* number of SM that must be used for side 1. Must be 0 for FFS */

/* Very naive code for finding the valuations of the units. */
static int
get_units_for_ab_from_file(int* u, int64_t a, uint64_t b, const char * abunitsfilename, unsigned int nbsm)
{
    FILE *in = fopen(abunitsfilename, "r");
    int64_t aa;
    uint64_t bb;
    int32_t ok = 0;
    while(fscanf(in, "%" PRId64 " %" PRIu64, &aa, &bb) != EOF){
	for(int i = 0; i < (int)nbsm; i++) {
            int ret = fscanf(in, "%d", u+i);
            ASSERT_ALWAYS(ret == 1);
        }
	if(aa == a && bb == b){
	    //	    printf("# GOTCHA a and b\n");
	    ok = 1;
	    break;
	}
    }
    fclose(in);
    if(ok == 0){
	fprintf(stderr, "GASP: no valuation(s) found ");
	fprintf(stderr, "for a=%" PRId64 " b=%" PRIu64 "\n", a, b);
	exit (EXIT_FAILURE);
    }
    return ok;
}

/* Semi-naive code for finding the valuations of the units. Read an index 
   first, then use the right file. We assume that directory dirname contains
   a file index, followed by files like 0, 1, 2, 3, ... so that their union
   contains the list of (a, b, u) by increasing values of a. We assume that
   all triples with the same a are in the same subfile.
*/
static int
get_units_for_ab(int* u, int64_t a, uint64_t b, const char * abunitsdirname,
		 unsigned int nbsm)
{
    size_t len = strlen(abunitsdirname);
    char *index = malloc((len+6+1) * sizeof(char));
    strncpy(index, abunitsdirname, len);
    strncpy(index+len, "/index", 6);
    index[len+6] = '\0';
    //printf("index=%s\n", index);
    /* read index file */
    FILE *in = fopen(index, "r");
    if(in == NULL){
	perror(index);
	return 0;
    }
    int n, ok = 0;
    int64_t aa;
    while(fscanf(in, "%" PRId64 " %d", &aa, &n) != EOF){
	if(a <= aa){
	    ok = 1;
	    break;
	}
    }
    fclose(in);
    free(index);
    if(ok == 0){
	fprintf(stderr, "Could not find abunit file for a=%" PRId64 "\n", a);
	exit (EXIT_FAILURE);
    }
    // TODO: use size of n
    char *fic = malloc((len+10) * sizeof(char));
    strncpy(fic, abunitsdirname, len);
    sprintf(fic+len, "/%d", n);
    //    printf("fic(%" PRId64 ")=%d => %s\n", a, n, fic);
    ok = get_units_for_ab_from_file(u, a, b, fic, nbsm);
    free(fic);
    return ok;
}

/* 
 * given S ==  data->sm_info[i], the range data->smlogs[i][0..S->nsm[
 * contains the logs of the SMs on side i.
 *
 */
/* Callback function called by filter_rels in compute_log_from_rels */
void *
thread_sm (void * context_data, earlyparsed_relation_ptr rel)
{
    read_data_ptr data = (read_data_ptr) context_data;
    log_rel_t *lrel = &(data->rels[rel->num]);

    mpz_ptr l = lrel->log_known_part;
    int64_t a = rel->a;
    uint64_t b = rel->b;

    uint64_t nonvoidside = 0; /* bit vector of which sides appear in the rel */
    for (weight_t i = 0; i < rel->nb; i++)
    {
      index_t h = rel->primes[i].h;
      int side = renumber_get_side_from_index (data->renum_tab, h, data->poly);
      nonvoidside |= ((uint64_t) 1) << side;
    }

    mpz_srcptr ell = data->log->ell;
    for(int side = 0 ; side < data->poly->nb_polys ; side++)
    {
      sm_side_info_srcptr S = data->sm_info[side];
      if (S->nsm > 0 && (nonvoidside & (((uint64_t) 1) << side)))
      {
        /* TODO: allow units for MNFS ??? */
        if (side < 2 && data->abunitsdirname[side]) {
            int e[S->nsm];
            /* This function returns an int, but an error is fatal, and
             * in fact trapped well before. If it derps, then there's no
             * point in trying to do anything, since we're dead already
             */
            get_units_for_ab(e, a, b, data->abunitsdirname[side], S->nsm);
            for(int i = 0; i < S->nsm; i++)
                mpz_addmul_si (l, data->smlogs[side][i], e[i]);
            mpz_mod(l, l, ell);
        } else {
            mpz_poly_t u;
            mpz_poly_init(u, MAX(1, S->f->deg-1));
            mpz_poly_setcoeff_int64(u, 0, a);
            mpz_poly_setcoeff_int64(u, 1, -b);
            compute_sm_piecewise(u, u, S);
            ASSERT_ALWAYS(u->deg < S->f->deg);
            for(int i = S->f->deg-1-u->deg; i < S->nsm; i++)
              mpz_addmul (l, data->smlogs[side][i], u->coeff[S->f->deg-1-i]);
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
  read_data_ptr data = (read_data_ptr) context_data;
  log_rel_t *lrel = &(data->rels[rel->num]);
  unsigned int next = 0;
  ideal_merge_t buf[REL_MAX_SIZE];
  int c = 0;
  for (unsigned int i = 0; i < rel->nb; i++)
  {
    index_t h = rel->primes[i].h;
    exponent_t e = rel->primes[i].e;

    if (mpz_sgn(data->log->tab[h]) >= 0) {
      mpz_addmul_si (lrel->log_known_part, data->log->tab[h], e);
      c++;
    } else
      buf[next++] = (ideal_merge_t) {.id = h, .e = e};
  }
  if (c) mpz_mod(lrel->log_known_part, lrel->log_known_part, data->log->ell);

  lrel->unknown = idealmerge_my_malloc (next);
  lrel->nb_unknown = next;
  memcpy(lrel->unknown, buf, next * sizeof(ideal_merge_t));

  return NULL;
}

/* Return the number of unknown logarithms in the relation.
 * rels[i].nb_unknown may not be up-to-date (can only be greater than the actual
 * value) */
static inline weight_t
nb_unknown_log (read_data_ptr data, uint64_t i)
{
  log_rel_t * lrel = &(data->rels[i]);
  weight_t j, k, len = lrel->nb_unknown;
  ideal_merge_t *p = lrel->unknown;
  int c = 0;
  for (j = 0, k = 0; k < len; k++)
  {
    pthread_mutex_lock (&lock);
    int sgn = mpz_sgn(data->log->tab[p[k].id]);
    pthread_mutex_unlock (&lock);

    if (sgn < 0) // we do not know the log if this ideal
    {
      if (j != k)
        p[j] = p[k];
      j++;
    }
    else { // We know this log, add it to log_know_part
      mpz_addmul_si(lrel->log_known_part, data->log->tab[p[k].id],
                          p[k].e);
        c++;
    }
  }
  if (c) mpz_mod(lrel->log_known_part, lrel->log_known_part, data->log->ell);

  lrel->nb_unknown = j;
  return j;
}

/* In a relation with 1 missing logarithm of exponent e, compute its values,
 * i.e. compute   dest <- (-vlog / e) mod ell
 * Return 0 if dest was already known (i.e. computed between the call to
 * nb_unknown_log and the call to this functions), return 1 otherwise.
 */
static inline unsigned int
compute_missing_log (mpz_t dest, mpz_t vlog, int32_t e, mpz_t ell)
{
  unsigned int ret;
  mpz_t tmp;
  mpz_init_set_si (tmp, e);
  mpz_invert (tmp, tmp, ell);
  mpz_neg (vlog, vlog);
  mpz_mul (vlog, vlog, tmp);
  mpz_mod (tmp, vlog, ell);

  pthread_mutex_lock (&lock);

  if (mpz_sgn(dest) < 0) // we do not already know the log
  {
    mpz_set (dest, tmp);
    ret = 1;
  }
  else // log was already computed by another thread
  {
    ret = 0;
  }
  pthread_mutex_unlock (&lock);

  mpz_clear (tmp);
  return ret;
}

/* Compute all missing logarithms for relations in [start,end[.
 * Return the number of computed logarithms */
static uint64_t
log_do_one_part_of_iter (read_data_ptr data, bit_vector not_used, uint64_t start,
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
        mpz_ptr vlog = data->rels[i].log_known_part;
        if (nb == 0 && mpz_cmp_ui (vlog, 0) != 0)
        {
          gmp_fprintf (stderr, "Error, no unknown log in rel %" PRIu64 " and sum"
                       " of log is not zero, sum is: %Zd\n", i, vlog);
          exit (EXIT_FAILURE);
        }
        else if (nb == 1)
        {
          ideal_merge_t ideal = data->rels[i].unknown[0];
          computed +=
             compute_missing_log (data->log->tab[ideal.id], vlog, ideal.e,
                                  data->log->ell);
        }
      }
    }
  }

  return computed;
}

#define log_do_one_iter_mono(d, bv, n) log_do_one_part_of_iter (d, bv, 0, n)

/************************** Dependancy graph *********************************/
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
graph_dep_set_log_already_known (graph_dep_t G, logtab_ptr log)
{
  for (uint64_t h = 0; h < log->nprimes; h++)
  {
    if (mpz_sgn(log->tab[h]) >= 0)
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

/* Structure containing the data necessary for reading rels for dep. graph */
struct dep_read_data_s
{
  light_rels_t rels;
  graph_dep_t G;
};

typedef struct dep_read_data_s dep_read_data_t[1];
typedef struct dep_read_data_s * dep_read_data_ptr;

/* Callback function called by filter_rels in compute_needed_rels */
void *
dep_thread_insert (void * context_data, earlyparsed_relation_ptr rel)
{
  dep_read_data_ptr data = (dep_read_data_ptr) context_data;
  light_rels_t lrel = &(data->rels[rel->num]);
  unsigned int next = 0;
  index_t buf[REL_MAX_SIZE];

  for (unsigned int i = 0; i < rel->nb; i++)
  {
    index_t h = rel->primes[i].h;
    if (GRAPH_DEP_IS_LOG_UNKNOWN(data->G, h))
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
dep_nb_unknown_log (dep_read_data_ptr data, uint64_t i, index_t *h)
{
  weight_t k, nb;
  index_t *p = data->rels[i].needed;

  for (nb = 0, k = 0; k < data->rels[i].len; k++)
  {
    pthread_mutex_lock (&lock);
    int unknow = GRAPH_DEP_IS_LOG_UNKNOWN (data->G, p[k]);
    pthread_mutex_unlock (&lock);

    if (unknow) // we do not know the log if this ideal
    {
      nb++;
      *h = p[k];
    }
  }
  return nb;
}

/* Compute all dependancies for relations in [start,end[.
 * Return the number of dependancies found */
static uint64_t
dep_do_one_part_of_iter (dep_read_data_ptr data, bit_vector not_used,
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
          data->G.tab[h].state = NODE_DEP_LOG_RECONSTRUCTED;
          data->G.tab[h].i = i;
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
    dep_read_data_ptr data = (dep_read_data_ptr) ti->data;
    ti->computed = dep_do_one_part_of_iter (data, not_used, start, end);
  }
  else
  {
    read_data_ptr data = (read_data_ptr) ti->data;
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
log_do_one_iter_mt (read_data_ptr d, bit_vector bv, int nt, uint64_t nrels)
{
  return do_one_iter_mt ((void *) d, bv, nt, nrels, 1);
}

/* Compute all missing logarithms possible. Run through all the relations once.
 * Multi thread version
 * Return the number of computed logarithms */
static inline uint64_t
dep_do_one_iter_mt (dep_read_data_ptr d, bit_vector bv, int nt, uint64_t nrels)
{
  return do_one_iter_mt ((void *) d, bv, nt, nrels, 0);
}

/***************** Important functions called by main ************************/
/* Read the logarithms computed by the linear algebra */
static void
read_log_format_LA (logtab_ptr log, const char *logfile, const char *idealsfile,
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
    FATAL_ERROR_CHECK (h >= log->nprimes, "Too big value of index");

    int ret = gmp_fscanf (flog, "%Zd\n", tmp_log);
    FATAL_ERROR_CHECK (ret != 1, "Error in file containing logarithms values");

    ASSERT_ALWAYS (col == i);
    logtab_insert (log, h, tmp_log);
    i++;
    if (stats_test_progress (stats))
      stats_print_progress (stats, i, 0, 0, 0);
  }
  stats_print_progress (stats, i, 0, 0, 1);
  ASSERT_ALWAYS (feof(fid));
  ASSERT_ALWAYS (i == ncols);

  unsigned int index = 0;
  for(int side = 0; side < nb_polys; side++)
  {
      for (int ism = 0; ism < sm_info[side]->nsm; ism++)
      {
        int ret = gmp_fscanf (flog, "%Zd\n", tmp_log);
        FATAL_ERROR_CHECK (ret != 1, "Error in file containing logarithms values");
        logtab_insert (log, log->nprimes+index, tmp_log);
        index++;
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
read_log_format_reconstruct (logtab_ptr log, MAYBE_UNUSED renumber_t renumb,
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
  stats_init (stats, stdout, &nread, nbits(renumb->size)-5, "Read", "logarithms", "",
              "logs");
  for (index_t i = 0; i < renumb->naddcols; i++)
  {
    ret = gmp_fscanf (f, "%" SCNid " added column %Zd\n", &h, tmp_log);
    ASSERT_ALWAYS (ret == 2 && h == i);
    nread++;
    logtab_insert (log, h, tmp_log);
  }
  for (int i = 0; i < renumb->bad_ideals.n; i++)
  {
    for (int k = 0; k < renumb->bad_ideals.nb[i]; k++)
    {
      ret = gmp_fscanf (f, "%" SCNid " bad ideals %Zd\n", &h, tmp_log);
      ASSERT_ALWAYS (ret == 2);
      nread++;
      logtab_insert (log, h, tmp_log);
    }
  }
  while (gmp_fscanf (f, "%" SCNid " %*" SCNpr " %*d %*s %Zd\n", &h, tmp_log)
          == 2)
  {
    nread++;
    logtab_insert (log, h, tmp_log);
    if (stats_test_progress (stats))
      stats_print_progress (stats, nread, 0, 0, 0);
  }
  stats_print_progress (stats, nread, 0, 0, 1);

  for (unsigned int nsm = 0; nsm < log->nbsm; nsm++)
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
    ASSERT_ALWAYS (h == (index_t) nsm + log->nprimes);
    logtab_insert (log, h, tmp_log);
  }
  ASSERT_ALWAYS (feof(f));

  mpz_clear (tmp_log);
  fclose_maybe_compressed (f, filename);
}

/* Write values of the known logarithms. */
static void
write_log (const char *filename, logtab_ptr log, renumber_t tab, 
	   cado_poly poly, sm_side_info *sm_info)
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
  mpz_t base;
  mpz_init (base);
  for (i = 0; i < log->nprimes + log->nbsm; i++)
  {
    if (mpz_sgn(log->tab[i]) > 0) /* The log is known and non-zero */
    {
      if (!base_already_set)
      {
        base_already_set = 1;
        /* base = 1/log->tab[i] mod ell */
        int ret = mpz_invert (base, log->tab[i], log->ell);
        ASSERT_ALWAYS (ret != 0);
        mpz_set_ui (log->tab[i], 1);
      }
      else
      {
        mpz_mul (log->tab[i], log->tab[i], base);
        mpz_mod (log->tab[i], log->tab[i], log->ell);
      }
    }
  }
  mpz_clear (base);

  uint64_t nknown = 0;
  stats_init (stats, stdout, &nknown, nbits(tab->size)-5, "Wrote",
              "known logarithms", "ideals", "logs");
  for (i = 0; i < tab->size; i++)
  {
	  if (mpz_sgn(log->tab[i]) >= 0) // we know the log of this ideal
    {
      nknown++;
      if (tab->table[i] == RENUMBER_SPECIAL_VALUE)
      {
        ASSERT_ALWAYS (mpz_cmp (log->tab[i], log->ell) < 0);
        if (renumber_is_additional_column (tab, i))
          gmp_fprintf (f, "%" PRid " added column %Zd\n", i, log->tab[i]);
        else
          gmp_fprintf (f, "%" PRid " bad ideals %Zd\n", i, log->tab[i]);
      }
      else
      {
        ASSERT_ALWAYS (mpz_cmp (log->tab[i], log->ell) < 0);
        p_r_values_t p, r;
        int side;
        renumber_get_p_r_from_index (tab, &p, &r, &side, i, poly);
        if (side != tab->rat)
          gmp_fprintf (f, "%" PRid " %" PRpr " %d %" PRpr " %Zd\n", i, p, side,
                                                               r, log->tab[i]);
        else
          gmp_fprintf (f, "%" PRid " %" PRpr " %d rat %Zd\n", i, p, side,
                                                                  log->tab[i]);
      }
      if (stats_test_progress (stats))
        stats_print_progress (stats, nknown, i+1, 0, 0);
    }
  }
  stats_print_progress (stats, nknown, tab->size, 0, 1);
  for (unsigned int nsm = 0, i = tab->size; nsm < log->nbsm; nsm++)
  {
    ASSERT_ALWAYS (mpz_sgn (log->tab[i+nsm]) >= 0);
    ASSERT_ALWAYS (mpz_cmp (log->tab[i+nsm], log->ell) < 0);
    // compute side
    int side, nsm_tot = sm_info[0]->nsm, jnsm = nsm;
    for(side = 0; ((int)nsm) >= nsm_tot; side++){
	nsm_tot += sm_info[side+1]->nsm;
	jnsm -= sm_info[side]->nsm;
    }
    ASSERT_ALWAYS ((jnsm >= 0) && (jnsm < sm_info[side]->nsm));
    gmp_fprintf (f, "%" PRid " SM %d %d %Zd\n", i+nsm, side, jnsm, log->tab[i+nsm]);
  }

  uint64_t missing = tab->size - nknown;
  printf ("# factor base contains %" PRIu64 " elements\n"
          "# logarithms of %" PRIu64 " elements are known (%.1f%%)\n"
          "# logarithms of %" PRIu64 " elements are missing (%.1f%%)\n",
          tab->size, nknown, 100.0 * nknown / (double) tab->size,
          missing, 100.0 * missing / (double) tab->size);
  fclose_maybe_compressed (f, filename);
  ASSERT_ALWAYS (log->nknown == nknown);
}

/* Given a filename, compute all the possible logarithms of ideals appearing in
 * the file. Return the number of computed logarithms.
 * Modify its first argument bit_vector needed_rels */
static uint64_t
compute_log_from_rels (bit_vector needed_rels,
                       const char *relspfilename, uint64_t nrels_purged,
                       const char *relsdfilename, uint64_t nrels_del,
                       uint64_t nrels_needed, int nt,
                       read_data_ptr data)
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
#if ! defined (FOR_FFS) && DEBUG >= 1
  printf ("# DEBUG: Using %d thread(s) for thread_sm\n", nt);
#endif
  fflush(stdout);
  char *fic[3] = {(char *) relspfilename, (char *) relsdfilename, NULL};
  struct filter_rels_description desc[3] = {
                   { .f = thread_insert, .arg=data, .n=1},
#ifndef FOR_FFS
                   { .f = thread_sm,     .arg=data, .n=nt},
#endif
                   { .f = NULL,          .arg=0,     .n=0}
      };
  filter_rels2 (fic, desc, EARLYPARSE_NEED_AB_HEXA | EARLYPARSE_NEED_INDEX,
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
    printf ("# Iteration %" PRIu64 ": starting...\n", iter);
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
                     logtab_ptr log, const char *wanted_filename, int nt)
{
  double wct_tt0, wct_tt;
  uint64_t total_computed = 0, iter = 0, computed;
  uint64_t nrels = nrels_purged + nrels_del;
  graph_dep_t dep_graph = graph_dep_init (log->nprimes);
  light_rels_t rels = light_rels_init (nrels);

  graph_dep_set_log_already_known (dep_graph, log);

  dep_read_data_t data;
  data->rels = rels;
  data->G = dep_graph;

  /* Init bit_vector to remember which relations were already used */
  bit_vector_set (needed_rels, 1);

  /* Reading all relations */
  printf ("# Reading relations from %s and %s\n", relspfilename, relsdfilename);
  fflush(stdout);
  char *fic[3] = {(char *) relspfilename, (char *) relsdfilename, NULL};
  filter_rels (fic, (filter_rels_callback_t) &dep_thread_insert, (void *) &data,
               EARLYPARSE_NEED_INDEX, NULL, NULL);

  /* computing dependencies */
  printf ("# Starting to compute dependancies from rels\n");

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
    total_computed += computed;

    printf ("# Iteration %" PRIu64 ": %" PRIu64 " new dependancies computed\n",
            iter, computed);
    printf ("# Iteration %" PRIu64 " took %.1fs (wall-clock time).\n",
            iter, wct_seconds() - wct_tt);

    iter++;
  } while (computed);

  printf ("# Computing dependancies took %.1fs (wall-clock time)\n",
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
    FATAL_ERROR_CHECK (h >= log->nprimes, "Too big value of index");
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
  param_list_decl_usage(pl, "gorder", "group order (see sm -gorder parameter)");
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
  param_list_decl_usage(pl, "partial", "(switch) do not reconstruct everything "
                                       "that can be reconstructed");
#ifndef FOR_FFS
  param_list_decl_usage(pl, "nsm", "number of SM's to add on side 0,1,...");
  param_list_decl_usage(pl, "abunits0", "units for all (a, b) pairs from purged and relsdels on side 0");
  param_list_decl_usage(pl, "abunits1", "units for all (a, b) pairs from purged and relsdels on side 1");
  param_list_decl_usage(pl, "explicit_units0", "use units for all (a, b) pairs from purged and relsdels on side 0");
  param_list_decl_usage(pl, "explicit_units1", "use units for all (a, b) pairs from purged and relsdels on side 1");
#endif
  param_list_decl_usage(pl, "mt", "number of threads (default 1)");
  param_list_decl_usage(pl, "wanted", "file containing list of wanted logs");
  param_list_decl_usage(pl, "force-posix-threads", "(switch)");
  param_list_decl_usage(pl, "path_antebuffer", "path to antebuffer program");
  verbose_decl_usage(pl);
}

static void
usage (param_list pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}


int
main(int argc, char *argv[])
{
  char *argv0 = argv[0];

  renumber_t renumber_table;
  uint64_t nrels_tot = 0, nrels_purged, nideals_purged, nrels_del, nrels_needed;
  uint64_t nprimes;
  int mt = 1;
  int partial = 0;
  int nsm_arg[NB_POLYS_MAX], nsm_tot;

  /* negative value means that the value that will be used is the value
   * computed later by sm_side_info_init */
  for (int side = 0; side < NB_POLYS_MAX; side++)
    nsm_arg[side] = -1;

  mpz_t ell;
  logtab_t log;
  cado_poly poly;

  param_list pl;
  param_list_init(pl);
  declare_usage(pl);
  argv++,argc--;

  param_list_configure_switch(pl, "partial", &partial);
  param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);

#ifndef FOR_FFS
  int units0 = 0;
  param_list_configure_switch(pl, "explicit_units0", &units0);
  int units1 = 0;
  param_list_configure_switch(pl, "explicit_units1", &units1);
  if(units0 != 0 && units1 != 0){
      fprintf (stderr, "units0 and units1 cannot be both set\n");
      exit (EXIT_FAILURE);
  }
#endif

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
  param_list_parse_uint64(pl, "nrels", &nrels_tot);
  param_list_parse_mpz(pl, "gorder", ell);
  param_list_parse_int(pl, "mt", &mt);
  const char *path_antebuffer = param_list_lookup_string(pl, "path_antebuffer");

  /* Some checks on command line arguments */
  if (mpz_cmp_ui (ell, 0) <= 0)
  {
    fprintf(stderr, "Error, missing -gorder command line argument "
                    "(or gorder <= 0)\n");
    usage (pl, argv0);
  }
  if (nrels_tot == 0)
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

  cado_poly_init (poly);
#ifndef FOR_FFS
  if (!cado_poly_read (poly, polyfilename))
#else
  if (!ffs_poly_read (poly, polyfilename))
#endif
  {
    fprintf (stderr, "Error reading polynomial file\n");
    exit (EXIT_FAILURE);
  }

#ifndef FOR_FFS
  /* Read number of sm to be printed from command line */
  param_list_parse_int_list (pl, "nsm", nsm_arg, poly->nb_polys, ",");
  for(int side = 0; side < poly->nb_polys; side++)
  {
    if (nsm_arg[side] > poly->pols[side]->deg)
    {
      fprintf(stderr, "Error: nsm%d=%d can not exceed the degree=%d\n",
                      side, nsm_arg[side], poly->pols[side]->deg);
      exit (EXIT_FAILURE);
    }
  }

  const char * abunits0dirname = NULL;
  abunits0dirname = param_list_lookup_string(pl, "abunits0");
  if(units0 == 0)
      abunits0dirname = NULL;

  const char * abunits1dirname = NULL;
  abunits1dirname = param_list_lookup_string(pl, "abunits1");
  if(units1 == 0)
      abunits1dirname = NULL;
#endif

  if (param_list_warn_unused(pl))
  {
    fprintf(stderr, "Error, unused parameters are given\n");
    usage(pl, argv0);
  }

  set_antebuffer_path (argv0, path_antebuffer);

  /* Init data for computation of the SMs. Need that now for nsm_tot */
  sm_side_info sm_info[NB_POLYS_MAX];
  nsm_tot = 0;
  for (int side = 0; side < poly->nb_polys; side++)
  {
#ifndef FOR_FFS
    sm_side_info_init(sm_info[side], poly->pols[side], ell);
    fprintf(stdout, "\n# Polynomial on side %d:\n# F[%d] = ", side, side);
    mpz_poly_fprintf(stdout, poly->pols[side]);
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
    nsm_tot += sm_info[side]->nsm;
#else
    sm_info[side]->nsm = 0;
#endif
  }
  fflush(stdout);


  /* Reading renumber file */
  printf ("\n###### Reading renumber file ######\n");
  renumber_init_for_reading (renumber_table);
  renumber_read_table (renumber_table, renumberfilename);
  nprimes = renumber_table->size;

  /* Read number of rows and cols on first line of purged file */
  purgedfile_read_firstline (relspfilename, &nrels_purged, &nideals_purged);
  nrels_del = nrels_tot - nrels_purged;

  /* Malloc'ing log tab and reading values of log */
  printf ("\n###### Reading known logarithms ######\n");
  fflush(stdout);
  logtab_init (log, nprimes, nsm_tot, ell);
  if (logformat == NULL || strcmp(logformat, "LA") == 0)
      read_log_format_LA (log, logfilename, idealsfilename, sm_info,
                                                            poly->nb_polys);
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
      read_data_t data;
      read_data_init (data, log, nrels_purged + nrels_del, poly, sm_info,
                      renumber_table, abunits0dirname, abunits1dirname);

      log->nknown += compute_log_from_rels (rels_to_process, relspfilename,
              nrels_purged, relsdfilename,
              nrels_del, nrels_needed, mt,
              data);
      printf ("# %" PRIu64 " logarithms are known.\n", log->nknown);
      read_data_free (data);
      extern double m_seconds;
      fprintf(stderr, "# %.2f\n", m_seconds);
  }
  else
    printf ("# All wanted logarithms are already known, skipping this step\n");
  fflush(stdout);

  /* Writing all the logs in outfile */
  printf ("\n###### Writing logarithms in a file ######\n");
  write_log (outfilename, log, renumber_table, poly, sm_info);

  /* freeing and closing */
  logtab_clear (log);
  mpz_clear(ell);

#ifndef FOR_FFS
  for (int side = 0 ; side < poly->nb_polys ; side++)
    sm_side_info_clear (sm_info[side]);
#endif

  renumber_clear (renumber_table);
  bit_vector_clear(rels_to_process);
  cado_poly_clear (poly);
  param_list_clear (pl);
  return EXIT_SUCCESS;
}
