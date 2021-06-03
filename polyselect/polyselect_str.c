#define TOKENPASTE(x, y) x ## y
#define TOKENPASTE2(x, y) TOKENPASTE(x, y)
#define LABEL_UNIQUE TOKENPASTE2(Label, __LINE__)

/* Data struct used for polyselect */
#include "cado.h" // IWYU pragma: keep
#include <float.h> // DBL_MAX
#include <math.h> // log
#include <string.h> // memset
#include <pthread.h> // pthread_mutex_lock
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "polyselect_str.h"
#include "cado_poly.h"
#include "getprime.h"   // getprime
#include "misc.h"   // nprimes_interval
#include "macros.h"

void match (unsigned long p1, unsigned long p2, int64_t i, mpz_srcptr m0,
            mpz_srcptr ad, unsigned long d, mpz_srcptr N, unsigned long q,
            mpz_srcptr rq);

void gmp_match (uint32_t p1, uint32_t p2, int64_t i, mpz_srcptr m0,
		mpz_srcptr ad, unsigned long d, mpz_srcptr N, uint64_t q,
		mpz_srcptr rq);

/* The following are primes used as factors in the special-q.
   Warning: if you add larger primes, you should ensure that any product
   still fits in an "unsigned long" (cf routine collision_on_sq).
   Here on a 32-bit machine the largest product of 4 primes is
   241*251*257*263 < 2^32, and on a 64-bit machine the largest product of 8
   primes is 233*239*241*251*257*263*269*271 < 2^64.
   LEN_SPECIAL_Q is defined in the header.
*/
#if ULONG_MAX == 4294967295UL
const unsigned int SPECIAL_Q[LEN_SPECIAL_Q] = {
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
  31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
  73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
  127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
  179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
  233, 239, 241, 251, 257, 263, 0 };
#else /* 64-bit */
const unsigned int SPECIAL_Q[LEN_SPECIAL_Q] = {
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
  31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
  73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
  127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
  179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
  233, 239, 241, 251, 257, 263, 269, 271, 0 };
#endif

static pthread_mutex_t lock=PTHREAD_MUTEX_INITIALIZER; 

//#define LESS_P

size_t expected_memory_usage_for_primes(unsigned long P)
{
  unsigned long Pmax = 2*P;
#ifdef LESS_P // if impatient for root finding
  Pmax = P + P/2;
#endif
  unsigned long maxprimes = nprimes_interval(P, Pmax);
  return maxprimes * sizeof (uint32_t);
}
    
/* init prime array */
unsigned long
initPrimes ( unsigned long P,
             uint32_t **primes )
{
  unsigned long p, nprimes = 0;
  unsigned long Pmax = 2*P;
#ifdef LESS_P // if impatient for root finding
  Pmax = P + P/2;
#endif
  unsigned long maxprimes = nprimes_interval(P, Pmax);

  *primes = (uint32_t*) malloc (maxprimes * sizeof (uint32_t));
  if ( (*primes) == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in initPrimes\n");
    exit (1);
  }

  prime_info pi;
  prime_info_init (pi);
  for (p = 2; p < P; p = getprime_mt (pi));

  while (p <= Pmax) {
    if (nprimes + 1 >= maxprimes) {
      maxprimes += maxprimes / 10;
      *primes = (uint32_t*) realloc (*primes, maxprimes * sizeof (uint32_t));
      if ( (*primes) == NULL) {
        fprintf (stderr, "Error, cannot reallocate memory in initPrimes\n");
        exit (1);
      }
    }
    (*primes)[nprimes++] = p;
    p = getprime_mt (pi);
  }

  prime_info_clear (pi);

  *primes = (uint32_t*) realloc (*primes, (nprimes) * sizeof (uint32_t));
  if ( (*primes) == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in initPrimes\n");
    exit (1);
  }

  return nprimes;
}


/* clear prime array */
void
printPrimes ( uint32_t *primes,
              unsigned long size )
{
  unsigned long i;
  for (i = 0; i < size; i++) {
    fprintf (stderr, "(%lu, %" PRIu32 ") ", i, primes[i]);
    if ((i+1) % 5 == 0)
      fprintf (stderr, "\n");
  }
  fprintf (stderr, "\n");
}


/* clear prime array */
void
clearPrimes ( uint32_t **primes )
{
  free (*primes);
}


/* init the header struct */
void
polyselect_poly_header_init (polyselect_poly_header_ptr header,
              mpz_ptr N,
              unsigned long d,
              mpz_ptr ad )
{
  /* compute Ntilde, m0 */
  mpz_init_set (header->N, N);
  mpz_init (header->Ntilde);
  mpz_init (header->m0);
  header->d = d;
  mpz_init_set (header->ad, ad);

  /* compute Ntilde, ... from N, ... */
  mpz_set (header->Ntilde, header->ad);
  mpz_mul_ui (header->Ntilde, header->Ntilde, header->d);
  mpz_pow_ui (header->Ntilde, header->Ntilde, header->d - 1);
  mpz_mul_ui (header->Ntilde, header->Ntilde, header->d);
  mpz_mul (header->Ntilde, header->Ntilde, header->N); /* d^d * ad^(d-1) * N */
  mpz_root (header->m0, header->Ntilde, header->d);
}


/* clear header struct */
void
polyselect_poly_header_clear (polyselect_poly_header_ptr header )
{
  mpz_clear (header->m0);
  mpz_clear (header->Ntilde);
  mpz_clear (header->N);
  mpz_clear (header->ad);
}

int
polyselect_poly_header_skip (polyselect_poly_header_srcptr header, unsigned long p)
{
  return header->d % p == 0 || mpz_divisible_ui_p (header->ad, p);
}

/* init polyselect_proots_t */
void
polyselect_proots_init (polyselect_proots_ptr R,
              unsigned long size )
{
  R->size = size;

  /* length of nr&roots is known now. lengths of roots[i] are TBD. */
  /* +1 for R->nr for end guard in collision_on_each_sq */
  R->nr = (uint8_t *) malloc ((size + 1) * sizeof (*(R->nr)));
  R->roots = (uint64_t **) malloc (size * sizeof (*(R->roots)));

  if (R->nr == NULL || R->roots == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in polyselect_proots_init().\n");
    exit (1);
  }
}


/* add a root to polyselect_proots_t */
void
polyselect_proots_add ( polyselect_proots_ptr R,
             unsigned long nr,
             uint64_t *roots,
             unsigned long index )
{
  unsigned int i;
  R->nr[index] = nr;

  if (nr != 0) {
    R->roots[index] = (uint64_t *) malloc (nr * sizeof (uint64_t));
    if (R->roots[index] == NULL) {
      fprintf (stderr, "Error, cannot allocate memory in polyselect_proots_add\n");
      exit (1);
    }

    for (i = 0; i < nr; i++)
      R->roots[index][i] = roots[i];
  }
  else
    R->roots[index] = NULL;
}


/* print roots */
void
polyselect_proots_print ( polyselect_proots_srcptr R,
               unsigned long size )
{
  unsigned int i, j;
  for (i = 0; i < size; i++) {
    if (R->nr[i] == 0) {
      fprintf (stderr, "NULL\n");
    }
    else {
      for (j = 0; j < R->nr[i]; j ++)
        fprintf (stderr, "%" PRIu64 " ", R->roots[i][j]);
      fprintf (stderr, "\n");
    }
  }
}


/* clear roots */
void
polyselect_proots_clear ( polyselect_proots_ptr R,
               unsigned long size )
{
  unsigned int i;

  free (R->nr);
  for (i = 0; i < size; i++)
    free (R->roots[i]);
  free (R->roots);
}


void
polyselect_qroots_init (polyselect_qroots_ptr R)
{
  R->alloc = 0;
  R->size = 0;
  R->q = NULL;
  R->nr = NULL;
  R->roots = NULL;
}

void
polyselect_qroots_realloc (polyselect_qroots_ptr R, unsigned long newalloc)
{
  ASSERT (newalloc >= R->size);
  R->alloc = newalloc;
  R->q = realloc (R->q, newalloc * sizeof (unsigned int));
  if (R->q == NULL)
  {
    fprintf (stderr, "Error, cannot reallocate memory in roots_realloc\n");
    exit (1);
  }
  R->nr = realloc (R->nr, newalloc * sizeof (unsigned int));
  if (R->nr == NULL)
  {
    fprintf (stderr, "Error, cannot reallocate memory in roots_realloc\n");
    exit (1);
  }
  R->roots = realloc (R->roots, newalloc * sizeof (uint64_t*));
  if (R->roots == NULL)
  {
    fprintf (stderr, "Error, cannot reallocate memory in roots_realloc\n");
    exit (1);
  }
}

/* reorder by decreasing number of roots (nr) */
void
polyselect_qroots_rearrange (polyselect_qroots_ptr R)
{
  if (R->size > 1) {
    unsigned int i, j, k, max, tmpq, tmpnr;
    uint64_t *tmpr = malloc (MAX_DEGREE * sizeof (uint64_t));

    for (i = 0; i < R->size; i ++) {
      max = i;
      for (j = i+1; j < R->size; j++) {
        if (R->nr[j] > R->nr[max]) {
          max = j;
        }
      }

      tmpq = R->q[i];
      tmpnr = R->nr[i];
      for (k = 0; k < MAX_DEGREE; k ++)
        tmpr[k] = R->roots[i][k];

      R->q[i] = R->q[max];
      R->nr[i] = R->nr[max];
      for (k = 0; k < MAX_DEGREE; k ++)
        R->roots[i][k] = R->roots[max][k];

      R->q[max] = tmpq;
      R->nr[max] = tmpnr;
      for (k = 0; k < MAX_DEGREE; k ++)
        R->roots[max][k] = tmpr[k];
    }
    free (tmpr);
  }
}

void
polyselect_qroots_add (polyselect_qroots_ptr R, unsigned int q, unsigned int nr, uint64_t *roots)
{
  unsigned int i;

  if (nr == 0)
    return;
  if (R->size == R->alloc)
    polyselect_qroots_realloc (R, R->alloc + R->alloc / 2 + 1);
  R->q[R->size] = q;
  R->nr[R->size] = nr;
  R->roots[R->size] = malloc (MAX_DEGREE * sizeof (uint64_t));
  if (R->roots[R->size] == NULL)
  {
    fprintf (stderr, "Error, cannot allocate memory in roots_add\n");
    exit (1);
  }
  for (i = 0; i < nr; i++)
    R->roots[R->size][i] = roots[i];
  R->size ++;
}

void
polyselect_qroots_print (polyselect_qroots_srcptr R)
{
  unsigned int i, j;
  for (i = 0; i < R->size; i++) {
    fprintf (stderr, "q: %u, r: ", R->q[i]);
    for (j = 0; j < R->nr[i]; j ++)
      fprintf (stderr, "%" PRIu64 " ", R->roots[i][j]);
    fprintf (stderr, "\n");
  }
}

void
polyselect_qroots_clear (polyselect_qroots_ptr R)
{
  unsigned int i;

  free (R->q);
  free (R->nr);
  for (i = 0; i < R->size; i++)
    free (R->roots[i]);
  free (R->roots);
}


/* init hash table */
void
polyselect_hash_init (polyselect_hash_ptr H, unsigned int init_size)
{
  H->alloc = init_size;
  H->slot = (struct polyselect_hash_slot_s *) malloc (H->alloc * sizeof (struct polyselect_hash_slot_s));
  if (H->slot == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in polyselect_hash_init\n");
    exit (1);
  }

  for (unsigned int j = 0; j < H->alloc; j ++) {
    H->slot[j].i = 0;
    H->slot[j].p = 0;
  }
  H->size = 0;
#ifdef DEBUG_HASH_TABLE
  H->coll = 0;
  H->coll_all = 0;
#endif
}

/* init_size is an approximation of the number of entries */
void
polyselect_shash_init (polyselect_shash_ptr H, unsigned int init_size)
{
  unsigned int init_size0 = init_size;

  /* round up to multiple of polyselect_SHASH_NBUCKETS */
  init_size = 1 + (init_size - 1) / polyselect_SHASH_NBUCKETS;
  init_size += init_size / 8 + 128; /* use 12.5% margin */
  if (init_size > init_size0)
    init_size = init_size0;
  H->alloc = init_size * (polyselect_SHASH_NBUCKETS + 1) + 8;
  /* + init_size for guard for the last buckets to avoid seg fault */
  /* + 8 for extreme guard (ASM X86 needs 8, C needs 5 when init_size is too small */
  H->mem = (uint64_t*) malloc (H->alloc * sizeof (uint64_t));
  if (!H->mem)
    {
      fprintf (stderr, "Error, cannot allocate memory in polyselect_shash_init\n");
      exit (1);
    }
  H->balloc = init_size;
}

void
polyselect_shash_reset (polyselect_shash_ptr H)
{
  H->base[0] = H->current[0] = H->mem;
  for (int j = 1; j <= polyselect_SHASH_NBUCKETS; j++)
    H->base[j] = H->current[j] = H->base[j-1] + H->balloc;
  /* Trick for prefetch T in polyselect_shash_find_collision after the end
     of the last bucket. */
  memset (H->base[polyselect_SHASH_NBUCKETS], 0, sizeof(**(H->base)) * 8);
}


/* rq is a root of N = (m0 + rq)^d mod (q^2) */
void
polyselect_hash_add (polyselect_hash_ptr H, unsigned long p, int64_t i, mpz_srcptr m0, mpz_srcptr ad,
          unsigned long d, mpz_srcptr N, unsigned long q, mpz_srcptr rq)
{
  uint32_t h;

  ASSERT(m0 != NULL);
  ASSERT(H->size < H->alloc);

  h = (uint32_t) i % H->alloc;

#ifdef DEBUG_HASH_TABLE
  if (H->slot[h].i != 0)
    H->coll ++;
#endif
  while (H->slot[h].i != 0)
  {
    if (H->slot[h].i == i) /* we cannot have H->slot[h].p = p, since for a
                       given prime p, all (p,i) values entered are different */
      match (H->slot[h].p, p, i, m0, ad, d, N, q, rq);
    if (UNLIKELY(++h == H->alloc))
      h = 0;
#ifdef DEBUG_HASH_TABLE
    H->coll_all ++;
#endif
  }
  H->slot[h].p = p;
  H->slot[h].i = i;
  H->size ++;
}

int
polyselect_shash_find_collision (polyselect_shash_srcptr H)
{
  static uint32_t size = 0, mask;
  uint64_t *Hj, *Hjm;
  uint32_t *T;
  uint32_t k;

#define polyselect_SHASH_RESEARCH(TH,I)				\
  do {							\
    key = ((I) >> 32) + (I);				\
    if (UNLIKELY(*TH)) do {				\
      if (UNLIKELY(*TH == key)) { free (T); return 1; }	\
    } while (*(++TH));					\
    *TH = key;						\
  } while (0)

#define polyselect_SHASH_TH_I(TH,I,IND)			\
  do {						\
    I = Hj[IND];				\
    TH = T + ((I >> LN2SHASH_NBUCKETS) & mask); \
    __builtin_prefetch(TH, 1, 3);		\
  } while (0)					\

  pthread_mutex_lock (&lock);
  if (!size) {
    size = H->balloc << 1;
    /* round up to power of 2 */
    size --;
    while (size & (size - 1))
      size &= size - 1;
    size <<= 2;
    ASSERT_ALWAYS((size & (size - 1)) == 0);
    mask = size - 1;
    size += 16; /* Guard to avoid to test the end of polyselect_hash_table when ++TH */
  }
  pthread_mutex_unlock (&lock);
  T = (uint32_t*) malloc (size * sizeof(*T));
  for (k = 0; k < polyselect_SHASH_NBUCKETS; k++) {
    Hj = H->base[k];
    Hjm = H->current[k];
    if (Hj == Hjm) continue;
    memset (T, 0, size * sizeof(*T));
    /* Here, a special guard at the end of polyselect_shash_init allows
       until Hjm[polyselect_SHASH_BUCKETS-1] + 5.
       So, it's not needed to test if Hj + 4 < Hjm to avoid prefetch problem. */
    uint32_t *Th0, *Th1, *Th2, *Th3, *Th4;
    uint64_t i0, i1, i2, i3, i4;
    unsigned int key;

    polyselect_SHASH_TH_I(Th0, i0, 0);
    polyselect_SHASH_TH_I(Th1, i1, 1);
    polyselect_SHASH_TH_I(Th2, i2, 2);
    polyselect_SHASH_TH_I(Th3, i3, 3);
    polyselect_SHASH_TH_I(Th4, i4, 4);
    Hj += 5;
    while (LIKELY(Hj < Hjm)) {
      __builtin_prefetch(((void *) Hj) + 0x280, 0, 3);
      polyselect_SHASH_RESEARCH(Th0, i0); polyselect_SHASH_TH_I(Th0, i0, 0);
      polyselect_SHASH_RESEARCH(Th1, i1); polyselect_SHASH_TH_I(Th1, i1, 1);
      polyselect_SHASH_RESEARCH(Th2, i2); polyselect_SHASH_TH_I(Th2, i2, 2);
      polyselect_SHASH_RESEARCH(Th3, i3); polyselect_SHASH_TH_I(Th3, i3, 3);
      polyselect_SHASH_RESEARCH(Th4, i4); polyselect_SHASH_TH_I(Th4, i4, 4);
      Hj += 5;
    }
    switch (Hj - Hjm) { /* no break: it's NOT an error! */
    case 0: polyselect_SHASH_RESEARCH(Th4, i4); no_break();
    case 1: polyselect_SHASH_RESEARCH(Th3, i3); no_break();
    case 2: polyselect_SHASH_RESEARCH(Th2, i2); no_break();
    case 3: polyselect_SHASH_RESEARCH(Th1, i1); no_break();
    case 4: polyselect_SHASH_RESEARCH(Th0, i0); // no_break();
    }
  }
  free (T);
  return 0;
}
#undef polyselect_SHASH_TH_I
#undef polyselect_SHASH_RESEARCH

/* return non-zero iff there is a collision */
#define PREFETCH 256
int
MAYBE_UNUSED polyselect_shash_find_collision_old (polyselect_shash_t H)
{
  static uint32_t size = 0, mask;
  uint64_t data[PREFETCH], *pdata, *edata, *ldata;
  uint32_t *Th;
  uint32_t key;
  uint64_t *Hj, *Hjm;
  uint32_t *T;
  uint64_t i;
  unsigned int j, k, l;
  
  if (!size) {
    size = H->balloc << 1;
    /* round up to power of 2 */
    size --;
    while (size & (size - 1))
      size &= size - 1;
    size <<= 2;
    ASSERT_ALWAYS((size & (size - 1)) == 0);
    mask = (size - 1);
    size += 16;
  }
  T = (uint32_t*) malloc (size * sizeof(*T));
  edata = data + PREFETCH;
  for (k = 0; k < polyselect_SHASH_NBUCKETS; k++) {
    Hj = H->base[k];
    Hjm = H->current[k];
    if (Hj == Hjm) continue;
    memset (T, 0, size * sizeof(*T));
    pdata = data;
    j = Hjm - Hj;
    for (l = 0; l < PREFETCH * 8; l += 64)
      __builtin_prefetch(((void *) data) + l, 1, 3);
    if (j > PREFETCH) j = PREFETCH;
    for (l = 0; l < j; l++) {
      i = Hj[l];
      key = (uint32_t) ((i >> 32) + i);
      i = (i >> LN2SHASH_NBUCKETS) & mask;
      __builtin_prefetch (T + i, 1, 3);
      i = (i << 32) + key;
      data[l] = i;
    }
    Hj += j;
    if (LIKELY(j == PREFETCH)) {
      while (LIKELY(Hj != Hjm)) {
	__builtin_prefetch(((void *) Hj) + 0x280, 0, 3);
	i = *pdata;
	key = (uint32_t) i;
	Th = T + (i >> 32);
	if (UNLIKELY(*Th)) do {
	  if (UNLIKELY(*Th == key)) { free (T); return 1; }
	} while (*(++Th));
	*Th = key;
	i = *Hj++;
	key = (uint32_t) ((i >> 32) + i);
	i = (i >> LN2SHASH_NBUCKETS) & mask;
	__builtin_prefetch (T + i, 1, 3);
	i = (i << 32) + key;
	*pdata++ = i;
	if (UNLIKELY (pdata == edata)) pdata = data;
      }
      ldata = pdata;
    }
    else
      ldata = data + l;
    do {
      i = *pdata++;
      key = (uint32_t) i;
      Th = T + (i >> 32);
      if (UNLIKELY (pdata == edata)) pdata = data;
      if (UNLIKELY(*Th)) do {
	if (UNLIKELY(*Th == key)) { free (T); return 1; }
      } while (*(++Th));
      *Th = key;
    } while (LIKELY(ldata != pdata));
  }
  free (T);
  return 0;
}

/* rq is a root of N = (m0 + rq)^d mod (q^2) */
void
polyselect_hash_add_gmp (polyselect_hash_ptr H, uint32_t p, int64_t i, mpz_srcptr m0, mpz_srcptr ad,
              unsigned long d, mpz_srcptr N, uint64_t q, mpz_srcptr rq)
{
  unsigned long h;

  if (H->size >= H->alloc)
    polyselect_hash_grow (H);
  if (i >= 0)
    h = ((int)i) % H->alloc;
  else
  {
    h = H->alloc - ( ((int)(-i)) % H->alloc);
    if (h == H->alloc)
      h = 0;
  }

  while (H->slot[h].p != 0)
  {
    if (m0 != NULL && H->slot[h].i== i && H->slot[h].p != p) {
      gmp_match (H->slot[h].p, p, i, m0, ad, d, N, q, rq);
    }
    if (++h == H->alloc)
      h = 0;
  }
  H->slot[h].p = p;
  H->slot[h].i = i;
  H->size ++;
}

void
polyselect_hash_clear (polyselect_hash_ptr H)
{
  free (H->slot);
}

void
polyselect_shash_clear (polyselect_shash_ptr H)
{
  free (H->mem);
}

void
polyselect_hash_grow (polyselect_hash_ptr H)
{
  unsigned long j, old_alloc;
  struct polyselect_hash_slot_s *old_slot;
  mpz_t tmp;

  mpz_init (tmp);
  mpz_set_ui (tmp, 0);

  old_alloc = H->alloc;
  old_slot = H->slot;
  H->alloc = 2 * old_alloc;
  H->slot = (struct polyselect_hash_slot_s *) malloc (H->alloc * sizeof (struct polyselect_hash_slot_s));
  if (H->slot == NULL) {
    fprintf (stderr, "Error, cannot allocate memory in polyselect_hash_init\n");
    exit (1);
  }
  memset (H->slot, 0, (sizeof(int64_t) + sizeof(uint32_t)) * H->alloc);
  H->size = 0;

  for (j = 0; j < old_alloc; j++)
    if (old_slot[j].p != 0)
      polyselect_hash_add (H, old_slot[j].p, old_slot[j].i, NULL, 0, 0, NULL, 0, tmp);

  free (old_slot);
  mpz_clear (tmp);
}

/****************************** polyselect_data_t functions ******************************/

void
polyselect_data_init (polyselect_data_ptr s)
{
  s->size = s->alloc = 0;
  s->x = NULL;
  s->sum = s->var = 0.0;
  s->min = DBL_MAX;
  s->max = -DBL_MAX;
}

void
polyselect_data_clear (polyselect_data_ptr s)
{
  free (s->x);
}

void
polyselect_data_add (polyselect_data_ptr s, double x)
{
  if (s->size == s->alloc)
    {
      s->alloc += 1 + s->alloc / 2;
      s->x = realloc (s->x, s->alloc * sizeof (double));
    }
  s->x[s->size++] = x;
  s->sum += x;
  s->var += x * x;
  if (x < s->min)
    s->min = x;
  if (x > s->max)
    s->max = x;
}

double
polyselect_data_mean (polyselect_data_srcptr s)
{
  return s->sum / (double) s->size;
}

double
polyselect_data_var (polyselect_data_srcptr s)
{
  double m = polyselect_data_mean (s);
  return s->var / (double) s->size - m * m;
}
