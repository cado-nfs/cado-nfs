/* purge --- remove singletons

Copyright 2008, 2009, 2010, 2011, 2012 Alain Filbois, Francois Morain,
                                       Paul Zimmermann

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* References:
 * On the Number Field Sieve Integer Factorisation Algorithm,
   Stefania Cavallar, PhD Thesis, University of Leiden, 2002.
*/

/*
  This program works in two passes over the relation files:
  - the first pass loads in memory only rational primes >= minpr and algebraic
    ideals >= minpa, but stores all ideals in the hash table, and keeps a count
    of the number of non-stored ideals for each relation.
    By default minpr and minpa are taken as rlim and alim respectively.
    Then simultaneously singleton removal is performed, and heavy relations
    are discarded, until the final excess is 'keep'.
  - the second pass goes through the relations again, and dumps the remaining
    ones in the format needed by 'merge'.

  This program uses the following data structures:
  rel_used[i]    - non-zero iff relation i is kept (so far)
  rel_compact[i] - list of 'h' indices in H table of considered (p,r) for row i
                   (terminated by a -1 sentinel)
  rel_weight[i]  - total weight of relation i.
  h = getHashAddr (H, p, r) - index of prime ideal (p, r) in hash table H
                              (rational primes use r = -2)
  GET_HASH_P(H,h) - prime corresponding to index h
  GET_HASH_R(H,h) - root  corresponding to index h (-2 for rational prime)
  H->hashcount[h] - number of occurrences of (p, r) in current relations
*/

/*
  HC_T : 8 bits.
  HT_T : (32 + 32 * need64) bits. Signed.
  HR_T : 32 if H.hm < 2^32-1, otherwise 64 bits.
  Note : sizeof(HT_T)>=sizeof(HR_T)
*/

#include "cado.h"
#include <gmp.h>
#include "mod_ul.c"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <errno.h>

#include "utils.h"
#ifdef FOR_FFS
#include "utils_ffs.h"
#endif

//#define STAT_FFS

//#define USE_CAVALLAR_WEIGHT_FUNCTION

#define MAX_FILES 1000000
#define DEFAULT_NPASS 50
#define DEFAULT_REQUIRED_EXCESS 0.1

typedef struct {
  volatile unsigned int ok;
  unsigned int num, end;
} fr_t;

/* Main variables */
static char rep_cado[4096];         /* directory of cado to find utils/antebuffer */

static hashtable_t H;
static HR_T **rel_compact  = NULL; /* see above */
static uint8_t *rel_weight = NULL; /* rel_weight[i] is the total weight of
                                      row i */
static char ** fic;
static char *pmin, *pminlessone;
static FILE *ofile;     /* For the principal file output. */
static bit_vector rel_used, Tbv;
static relation_stream rs;
static HR_T *sum; /* sum of row indices for primes with weight 2 */
static cado_poly pol;
static double wct0, W; /* total weight for last pass */
static size_t tot_alloc, tot_alloc0;
static HR_T nrel,
  nprimes = 0,
  nrelmax = 0,
  relsup,
  prisup,
  newnrel,
  newnprimes,
  Hsize,
  Hsizer,
  Hsizea,
  keep = 160;         /* default maximum final excess */
static HT_T minpr = UMAX(minpr);
static HT_T minpa = UMAX(minpa); /* negative values mean use minpr=rlim and
				    minpa=alim */
static int raw = 0, need64;
static unsigned int npt = 4;
static float w_ccc;
static uint8_t binfilerel;    /* True (1) if a rel_used relations file must be read */
static uint8_t boutfilerel;   /* True (1) if a rel_used relations file must be written */

#ifdef FOR_FFS
static FILE *ofile2;
static int pipe_2;
#ifdef STAT_FFS
uint64_t __stat_count[11] = {0,0,0,0,0,0,0,0,0,0,0};
uint64_t __stat_nonzero = 0;
#endif
#endif

static const unsigned char ugly[256] = {
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  0,   1,   2,   3,   4,   5,   6,   7,   8,   9, 255, 255, 255, 255, 255, 255,
  255,  10,  11,  12,  13,  14,  15, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255,  10,  11,  12,  13,  14,  15, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
  255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255 };

static const unsigned char nbbits[256] = {
  0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
  4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8 };

/*
  1<<NFR = Number of computation threads.
  The best is 2 (NFR = 1).
  if you have 2 cores or less, you could try NFR=1.
  if you have >4 cores, or hyperthreading, AND very slow cores (< 2 Ghz): NFR=2
*/
#define NFR (1)
/* 1<<NNFR = Number of sentences computed of (find root + hashkey of each term)
   computed in one pass. If the number is too big, the buffer between these
   threaded computations and insertRelation (which cannot be parallelized,
   because the hashkey is not bijective) is too big (memory waste) and
   the pipe-line is slow to start;
   If the number is too small, the computation threads are too often in
   nanosleep to keep CPU.
   NNFR=8 to 16, the greatest is the fastest for fast processors & fast
   memory; 14 seems to be the best choice (maximal speed for medium memory use).
*/
#define NNFR (14)

/* Size of relations buffer between parsing & insertRelation.
   CAREFUL! T_REL must be greater (at least double) than
   (1<<(NFR+NNFR+(NFR==0))).
   Stocks the sentences precomputed but not insered. About
   64K sentences for the optimal.
*/
#define T_REL MAX((1<<(NFR+NNFR+1+(NFR==0))),(1<<6))

/* The realistic minimal non-CPU waiting with nanosleep is about
   10 to 40 µs (1<<13 for nanosleep).
   But all the I/O between the threads have been buffered,
   and a thread do a nanosleep only if its buffer is empty.
   So I use here ~2ms (1<<21) to optimize CPU scheduler.
   Max pause is about 4 to 8ms (1<<22, 1<<23); after the program
   is slow down.
*/
static const struct timespec wait_classical = { 0, 1<<21 };

#define NB_RATP_OPT (12)
#define NB_ALGP_OPT (19)
#define NB_HK_OPT (NB_RATP_OPT + NB_ALGP_OPT - 4)
typedef struct {
  HR_T *hk;          /* The renumbers of primes in the rel. */
  relation_t rel;    /* Relation itself */
  unsigned int lhk;  /* Actual number of hk */
  unsigned int mhk;  /* Size max of the local hk; after free + malloc again */
  unsigned int ltmp; /* Size of renumbers for rel_compact */
  HR_T num;          /* Relation number */
} buf_rel_t;
static buf_rel_t *buf_rel;

typedef struct {
  rat_prime_t rp[NB_RATP_OPT];
  alg_prime_t ap[NB_ALGP_OPT];
  HR_T hk[NB_HK_OPT];
} buf_rel_data_t;
static buf_rel_data_t *buf_rel_data;

/* For the multithread sync */
static volatile unsigned long cpt_rel_a;
static volatile unsigned long cpt_rel_b;
static volatile unsigned int end_insertRelation = 0;

/* copied from utils/antebuffer.c */
#ifndef HAVE_NANOSLEEP
  int nanosleep(const struct timespec *req, struct timespec *rem) {
    if (rem == NULL) {
      /* Dummy to shut up the warning */
    }
#ifdef HAVE_USLEEP
    unsigned long usec = req->tv_sec * 1000000UL + req->tv_nsec / 1000UL;
    usleep(usec);
#else
    sleep(req->tv_sec);
#endif
    return 0;
  }
#endif

/* Be careful. 1<<13 = 8µs; in fact, about 30-50 µs at least
   with a very reactive machine. Use for debugging only.
*/
inline void attente_minimale_passive ()
{
  static const struct timespec wait_min = { 0, 1<<13 };
  nanosleep(&wait_min, NULL);
}

/* Don't use it, it's a CPU waste. */
inline void attente_minimale_active ()
{
  /* about 1 (for a nehalem 3 Ghz) to 5 µs */
  unsigned int i = (1<<6);
  while (i--)
    __asm__ volatile ( "\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
nop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\nnop\t\n\
" : : );
}

/********************** own memory allocation routines ***********************/

/* Rationale: calling one malloc() for each relation in insertNormalRelation()
   is expensive, since malloc() allocates some extra information to keep track
   of every memory blocks. Instead, we allocate memory in big blocks of size
   BLOCK_SIZE. */

#define BLOCK_SIZE (1<<20)  /* memory blocks are allocated of that # of HR_T's */
/* relcompact_list is a list of blocks, each one of BLOCK_SIZE HR_Ts */
static HR_T **relcompact_list = NULL;
static HR_T *myrelcompact;
static unsigned int relcompact_used = BLOCK_SIZE; /* usage of current block */
static unsigned int relcompact_size = 0;  /* minimal size of relcompact_list */
static int relcompact_current = -1; /* index of current block */

/* return a pointer to an array of n (HR_T) */
static HR_T *
my_malloc (unsigned int n)
{
  HR_T *ptr;

  if (relcompact_used + n > BLOCK_SIZE) {
    relcompact_used = 0;
    if (((unsigned int) (++relcompact_current)) == relcompact_size) {
      relcompact_size = relcompact_size ? (relcompact_size << 1) : (1<<16);
      if (!(relcompact_list = (HR_T **) realloc(relcompact_list, relcompact_size * sizeof(HR_T *)))) {
	  fprintf (stderr, "my_malloc_int: realloc error : %s\n", strerror (errno));
	  exit (1);
	}
      }
    SMALLOC(relcompact_list[relcompact_current], BLOCK_SIZE, "my_malloc_int 1");
    myrelcompact = relcompact_list[relcompact_current];
  }
  ptr = &(myrelcompact[relcompact_used]);
  relcompact_used += n;
  tot_alloc += n * sizeof (HR_T);
  return ptr;
}

static void
my_malloc_free_all (void)
{
  while (relcompact_current >= 0) {
    SFREE (relcompact_list[relcompact_current]);
    relcompact_current--;
  }
  SFREE (relcompact_list);
  relcompact_list = NULL;
  relcompact_used = BLOCK_SIZE;
  relcompact_size = 0;
}

/*****************************************************************************/

/* Print the relation 'rel' in matrix format, i.e., a line of the form:

   i a b k t_1 t_2 ... t_k

   i (decimal) is the row index from the nodup file (starting at 0)
   a (signed decimal) is a
   b (signed decimal) is b
   k (decimal) is the number of rational and algebraic primes in the relation
   t_1 ... t_k (hexadecimal) are the indices of the primes (starting at 0)

   Return the weight of the relation.

   Assumes the remaining primes in 'rel' are those with an odd exponent,
   and are all different.

   WARNING: the primes in the input relation are not necessarily sorted.

   WARNING2: the keys in *phk are not correct. There are only
   (HC1*A+HC2*B)%M. To compute the real keys, all keys must be computed
   by REALKEY.
*/

#define REALKEY(A,B)							\
  for (h = *phk++; H.ht[h].p != (HT_T) (A) || H.ht[h].r != (HT_T) (B); ) \
    if (++h >= H.hm) h = 0
#define FFSCOPYDATA(E) \
  t = p - op;								\
  for (i = (unsigned int) ((E) - 1); i--; p += t) memcpy (p, op, t)
#define WRITEP					\
  *p++ = ' ';					\
  p = u64toa16(p, (uint64_t) (H.hr[h] -1))

static int
fprint_rel_row (FILE *file, buf_rel_t *my_buf_rel)
{
  char buf[1<<12], *p;
  unsigned int nb_coeff;
  HR_T *phk = my_buf_rel->hk, h;
  rat_prime_t *prp;
  alg_prime_t *pap;

  p = u64toa10(buf, (uint64_t) my_buf_rel->num);   *p++ = ' ';
  p = d64toa10(p,   (int64_t)  my_buf_rel->rel.a); *p++ = ' ';
  p = u64toa10(p,   (uint64_t) my_buf_rel->rel.b); *p++ = ' ';

#ifndef FOR_FFS
  nb_coeff = my_buf_rel->rel.nb_ap + (my_buf_rel->rel.b ? my_buf_rel->rel.nb_rp : 1);
#else
  if (my_buf_rel->rel.b) {
    nb_coeff = 0;
    for (prp =  my_buf_rel->rel.rp;
         prp != &(my_buf_rel->rel.rp[my_buf_rel->rel.nb_rp]);
         nb_coeff += (prp++)->e);
    for (pap =  my_buf_rel->rel.ap;
         pap != &(my_buf_rel->rel.ap[my_buf_rel->rel.nb_ap]);
         nb_coeff += (pap++)->e);
  }
  else {
    nb_coeff = 1;
    for (pap =  my_buf_rel->rel.ap;
         pap != &(my_buf_rel->rel.ap[my_buf_rel->rel.nb_ap]);
         nb_coeff += (pap++)->e);
  }
#endif
  p = u64toa10(p, (uint64_t) nb_coeff);

#ifndef FOR_FFS
  if (my_buf_rel->rel.b) {
    for (prp =  my_buf_rel->rel.rp;
         prp != &(my_buf_rel->rel.rp[my_buf_rel->rel.nb_rp]);
         prp++) {
      REALKEY(prp->p, prp->p + 1);
      WRITEP;
    }
    for (pap =  my_buf_rel->rel.ap;
         pap != &(my_buf_rel->rel.ap[my_buf_rel->rel.nb_ap]);
         pap++) {
      REALKEY(pap->p, pap->r);
      WRITEP;
    }
  }
  else {
    REALKEY(my_buf_rel->rel.a, my_buf_rel->rel.a + 1);
    WRITEP;
    for (pap =  my_buf_rel->rel.ap;
         pap != &(my_buf_rel->rel.ap[my_buf_rel->rel.nb_ap]);
         pap++) {
      REALKEY(my_buf_rel->rel.a, pap->p);
      WRITEP;
    }
  }
#else
  if (my_buf_rel->rel.b) {
    char *op;
    size_t t;
    unsigned int i;

    for (prp =  my_buf_rel->rel.rp;
         prp != &(my_buf_rel->rel.rp[my_buf_rel->rel.nb_rp]);
         prp++)
      if (prp->e > 0) {
        op = p;
        REALKEY(prp->p, prp->p + 1);
        WRITEP;
        FFSCOPYDATA(prp->e);
      }
    for (pap =  my_buf_rel->rel.ap;
         pap != &(my_buf_rel->rel.ap[my_buf_rel->rel.nb_ap]);
         pap++)
      if (pap->e > 0) {
        op = p;
        REALKEY(pap->p, pap->r);
        WRITEP;
        FFSCOPYDATA(pap->e);
      }
  }
  else {
    char *op;
    size_t t;
    unsigned int i;

    REALKEY(my_buf_rel->rel.a, my_buf_rel->rel.a + 1);
    WRITEP;
    for (pap =  my_buf_rel->rel.ap;
         pap != &(my_buf_rel->rel.ap[my_buf_rel->rel.nb_ap]);
         pap++)
      if (pap->e > 0) {
        op = p;
        REALKEY(my_buf_rel->rel.a, pap->p);
        WRITEP;
        FFSCOPYDATA(pap->e);
      }
  }
#endif
  *p = '\n'; p[1] = 0;
  fputs(buf, file);

#ifndef FOR_FFS
  return nb_coeff;
#else
  return weight_rel_ffs (my_buf_rel->rel);
#endif
}

/* First we count the number of large primes; then we store all primes in
   the hash table, but not in the relations. This might end up with singletons
   here and there, but we don't care, since they will be dealt with in
   merge.

   Meaning of the different parameters:
   minpr, minpa: only ideals > minpr (resp. minpa) are considered on the
                 rational (resp. algebraic) side. This means that the output
                 might contain ideals <= minpr or minpa appearing only once.
*/

static inline void
insertNormalRelation (unsigned int j)
{
  ht_t *my_ht;
  HC_T *my_hc, *end_hc;
  HR_T *phk, *my_tmp;
  buf_rel_t *my_br;
  rat_prime_t *my_rp;
  alg_prime_t *my_ap;
  ht_t pr;
  HR_T h;
  unsigned int i, itmp;

  end_hc = &(H.hc[H.hm]);
  itmp = 0;
  my_br = &(buf_rel[j]);
  my_tmp = boutfilerel ? NULL : my_malloc(my_br->ltmp);
  phk = my_br->hk;
  i = my_br->rel.nb_rp;
  my_rp = my_br->rel.rp;
  while (i--)
    {
      /* we insert all ideals (even small ones) in the hash table since we
         need to renumber them in the second pass */

      /* The next code is exactly :
	 h = hashInsertWithKey(&H, (HT_T) buf_rel[j].rel.rp[i].p,
	 (HT_T) buf_rel[j].rel.rp[i].p + 1, *phk++, &np);
	 nprimes += np;
	 but it's so critical in perfs I prefer inlining and uses dirty
	 (really dirty!) tricks to help the compiler to speed up the code at the max.
      */
      pr = (ht_t) { (HT_T) my_rp->p, (HT_T) (my_rp->p + 1) };
      my_rp++;
      if (boutfilerel && pr.p < minpr) continue;
      h = *phk++;
      my_ht = &(H.ht[h]);
      my_hc = &(H.hc[h]);
      /* Highest critical section */
    t1:
      if (my_ht->p == pr.p && my_ht->r == pr.r) goto t11;
      if (!*my_hc) goto t12;
      my_ht++;
      my_hc++;
      if (my_hc != end_hc) goto t1;
      my_ht = H.ht;
      my_hc = H.hc;
      goto t1;
    t11:
      if (*my_hc != UMAX(*my_hc)) (*my_hc)++;
      if (!boutfilerel && pr.p >= minpr) my_tmp[itmp++] = my_hc - H.hc;
      continue;
    t12:
      *my_ht = pr;
      *my_hc = 1;
      nprimes++;
      if (!boutfilerel && pr.p >= minpr) my_tmp[itmp++] = my_hc - H.hc;
    }
  i = my_br->rel.nb_ap;
  my_ap = my_br->rel.ap;
  while (i--)
    {
      /*
	h = hashInsertWithKey (&H, (HT_T) buf_rel[j].rel.ap[i].p, (HT_T) buf_rel[j].rel.ap[i].r,
	*phk++, &np);
        nprimes += np;
      */
      pr = (ht_t) { (HT_T) my_ap->p, (HT_T) my_ap->r };
      my_ap++;
      if (boutfilerel && pr.p < minpa) continue;
      h = *phk++;
      my_ht = &(H.ht[h]);
      my_hc = &(H.hc[h]);
    t2:
      if (my_ht->p == pr.p && my_ht->r == pr.r) goto t21;
      if (!*my_hc) goto t22;
      my_ht++;
      my_hc++;
      if (my_hc != end_hc) goto t2;
      my_ht = H.ht;
      my_hc = H.hc;
      goto t2;
    t21:
      if (*my_hc != UMAX(*my_hc)) (*my_hc)++;
      if (!boutfilerel && pr.p >= minpa) my_tmp[itmp++] = my_hc - H.hc;
      continue;
    t22:
      *my_ht = pr;
      *my_hc = 1;
      nprimes++;
      if (!boutfilerel && pr.p >= minpa) my_tmp[itmp++] = my_hc - H.hc;
    }
  /* ASSERT_ALWAYS(++itmp == my_br->ltmp); */
  /* total relation weight */
  my_br->lhk = my_br->hk - phk;
  if (!boutfilerel) {
    my_tmp[itmp] = UMAX(*my_tmp); /* sentinel */
#ifndef FOR_FFS
    rel_weight[my_br->num] = my_br->rel.nb_rp + my_br->rel.nb_ap;
#else
    rel_weight[my_br->num] = weight_rel_ffs (my_br->rel);
#endif
    rel_compact[my_br->num] = my_tmp;
  }
}

/* The information is stored in the ap[].p part, which is odd, but convenient.
   rel->ap.p[0..deg[ contains the deg roots of f(x) mod p, leading to deg
   ideals for the factorization of (p). */
static inline void
insertFreeRelation (unsigned int j)
{
  HR_T *my_tmp, *phk, h;
  unsigned int i, itmp, np;

  /* the prime on the rational side is rel->a
     the prime ideal on the algebraic side are (rel->a, rel->ap[i].p) */

  itmp = 0;
  h = 0;
  my_tmp = boutfilerel ? NULL : my_malloc(buf_rel[j].ltmp);
  phk = buf_rel[j].hk;
  /* insert all ideals */
  if (!boutfilerel || (HT_T) buf_rel[j].rel.a >= minpr) {
    h = hashInsertWithKey(&H, (HT_T) buf_rel[j].rel.a, (HT_T) (buf_rel[j].rel.a + 1), *phk++, &np);
    nprimes += np;
    if (!boutfilerel && (HT_T) buf_rel[j].rel.a >= minpr) my_tmp[itmp++] = h;
  }

  if (!boutfilerel || (HT_T) buf_rel[j].rel.a >= minpa)
    for (i = 0; i < buf_rel[j].rel.nb_ap; i++) {
      h = hashInsertWithKey(&H, (HT_T) buf_rel[j].rel.a, (HT_T) buf_rel[j].rel.ap[i].p, *phk++, &np);
      nprimes += np; /* (H->hc[h] == 1); new ideal */
      if (!boutfilerel && (HT_T) buf_rel[j].rel.a >= minpa) my_tmp[itmp++] = h;
    }
  /* ASSERT_ALWAYS(++itmp == buf_rel[j].ltmp); */
  /* total relation weight */
  buf_rel[j].lhk = buf_rel[j].hk - phk;
  if (!boutfilerel) {
    my_tmp[itmp] = UMAX(*my_tmp);  /* sentinel */
#ifndef FOR_FFS
    rel_weight[buf_rel[j].num] = 1 + buf_rel[j].rel.nb_ap;
#else
    rel_weight[buf_rel[j].num] = weight_rel_ffs (buf_rel[j].rel);
#endif
    rel_compact[buf_rel[j].num] = my_tmp;
  }
}

/* Delete a relation: set rel_used[i] to 0, update the count of primes
   in that relation, and set rel_compact[i] to NULL.
   Warning: we only update the count of primes that we consider, i.e.,
   rational primes >= minpr and algebraic primes >= minpa.
*/
static void
delete_relation (HR_T i)
{
  HR_T *tab;
  HC_T *o;

  for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++) {
    o = &(H.hc[*tab]);
    ASSERT(*o);
    if (*o < UMAX(*o) && !(--(*o))) newnprimes--;
  }
  /* rel_compact[i] = NULL; */
  bit_vector_clearbit(rel_used, (size_t) i);
}

/* New pruning code, which optimizes the decrease of N*W where N is the number
   of rows, and W is the total weight. We consider the connected
   components of the relation R(i1,i2) iff i1 and i2 share a prime
   of weight 2. If we remove one component of n rows and total weight w,
   then we save w*N+n*W (neglecting 2nd order terms), thus we remove
   first the components with the largest value of n/N + w/W. */

typedef struct {
  float w;
  HR_T i;
} comp_t;

static int
compare (const void *v1, const void *v2)
{
  comp_t *w1 = (comp_t*) v1;
  comp_t *w2 = (comp_t*) v2;

  return (w1->w >= w2->w) ? -1 : 1;
}

float
weight_function_clique (HC_T w)
{
#ifdef USE_CAVALLAR_WEIGHT_FUNCTION
  if (w >= 3)
    return ldexpf (1, -(w-1));
  else if (w == 2)
    return 0.25;
  else
      return 0.0;
#else
    if (w >= 3)
      /* we use the multiplier 5 here, so that the average weight (assumed to
         be 1) is in the middle of the Count[10] array */
      return (float) 5.0 / (float) w;
    else
      return 0.0;
#endif
}

/* Compute connected component of row i for the relation R(i1,i2) if rows
   i1 and i2 share a prime of weight 2.
   Return number of rows of components, and put in w the total weight. */
static HR_T
compute_connected_component (HR_T i)
{
  HR_T *myrel_compact = rel_compact[i], h, k, n = 1;

  bit_vector_setbit(Tbv, (size_t) i); /* mark row as visited */
  while ((h = *myrel_compact++) != UMAX(h))
    {
      if (H.hc[h] == 2)
        {
          k = sum[h] - i; /* other row where prime of index h appears */
          if (!bit_vector_getbit(Tbv, (size_t) k)) /* row k not visited yet */
            n += compute_connected_component (k);
        }
      w_ccc += weight_function_clique (H.hc[h]);
    }
  return n;
}

/* Delete connected component of row i, assuming the bit-vector is set.
   Warning: we might have some H->hashcount[h] = 3, which is decreased
   to 2, but we don't want to treat that case. Thus we check in addition
   that sum[h] <> 0, which only occurs when H->hashcount[h] = 2 initially. */
static HR_T
delete_connected_component (HR_T i)
{
  HR_T *myrel_compact = rel_compact[i], h, k, w = 1;

  bit_vector_clearbit(Tbv, (size_t) i); /* mark row as visited */
  /* bit i of rel_used is cleared in delete_relation below */
  while ((h = *myrel_compact++) != UMAX(h)) {
    if (H.hc[h] == 2 && sum[h]) { /* first row that contains ideal of index h */
      k = sum[h] - i; /* other row where prime of index h appears */
      if (bit_vector_getbit(Tbv, (size_t) k) == 1) /* row k was not visited yet */
	w += delete_connected_component (k);
    }
  }
  delete_relation (i);
  return w;
}

static void
deleteHeavierRows (unsigned int npass)
/* int *nrel, int *nprimes, int nrelmax, int keep)
   &newnrel, &newnprimes, nrelmax, keep */
{
  static unsigned int count = 0;
  static HR_T chunk;
  comp_t *tmp = NULL; /* (weight, index) */
  HR_T *myrelcompact, i, h;
  double W = 0.0; /* total matrix weight */
  double N = 0.0; /* number of rows */
  unsigned int wceil, j, ltmp = 0, alloctmp = 0xff;
  long target;
#define MAX_WEIGHT 10
  unsigned int Count[MAX_WEIGHT], *pc; /* Count[w] is the # of components of weight >= w */

  if (newnrel - newnprimes <= keep)
    return;
  if (!count)
    chunk = (newnrel - newnprimes) / npass;

  /* first collect sums for primes with weight 2, and compute total weight */
  MEMSETZERO(sum, H.hm);
  for (i = 0; i < nrelmax; i++)
    if (bit_vector_getbit(rel_used, (size_t) i)) {
      for (myrelcompact = rel_compact[i]; (h = *myrelcompact++) != UMAX(h); )
	if (H.hc[h] == 2) sum[h] += i;
      N += 1.0;
      W += rel_weight[i]; /* row weight */
    }
  fprintf (stderr, "Step %u on %u: Matrix has %1.0f real (non null) rows "
                   "and weight %1.0f\n", count, npass, N, W);
  ASSERT_ALWAYS(N == (double) newnrel);

  /* now initialize bit table for relations used */
  bit_vector_init(Tbv, (size_t) nrelmax);
  bit_vector_neg (Tbv, rel_used);
  memset(Count, 0, sizeof(unsigned int) * MAX_WEIGHT);
  SMALLOC(tmp, alloctmp, "deleteHeavierRows 2");

  for (i = 0; i < nrelmax; i++)
    {
    if (!bit_vector_getbit(Tbv, (size_t) i)) {
      w_ccc = 0.0;
      h = compute_connected_component (i);
      wceil = (unsigned int) w_ccc + 1;
      /* don't consider weight-0 components */
      if ((w_ccc > 0.0) && ((wceil >= MAX_WEIGHT) || (Count[wceil] < chunk)))
	{
	  if (ltmp >= alloctmp)
            {
              alloctmp = 2 * alloctmp + 1;
              tmp = (comp_t *) realloc (tmp, alloctmp * sizeof (comp_t));
            }
	  tmp[ltmp++] = (comp_t) { w_ccc, i };
	  if (wceil > MAX_WEIGHT)
	    wceil = MAX_WEIGHT;
	  for (pc = &(Count[wceil]); pc-- != Count; (*pc)++);
	}
    }
    }
  qsort (tmp, ltmp, sizeof(comp_t), compare);

  /* remove heaviest components, assuming each one decreases the excess by 1;
     we remove only part of the excess at each call of deleteHeavierRows,
     hoping to get "better" components to remove at the next call. */
  if (++count < npass)
    {
      target = ((long) newnrel) - newnprimes - chunk;
      if (target < (long) keep) target = keep;
    }
  else
    target = keep; /* enough steps */

  for (j = 0; j < ltmp && newnrel > target + newnprimes; j++)
    newnrel -= delete_connected_component (tmp[j].i);
  fprintf (stderr, "deleted %u heavier connected components at %2.2lf\n",
                     j, seconds ());
  bit_vector_clear(Tbv);
  free (tmp);
}

#ifndef HAVE_SYNC_FETCH
static void
onepass_singleton_removal ()
{
  HR_T *tab, i;

  for (i = 0; i < nrelmax; i++)
    if (bit_vector_getbit(rel_used, (size_t) i))
      for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++)
	if (H.hc[*tab] == 1) {
	  delete_relation(i);
	  newnrel--;
	  break;
	}
}

#else /* ifndef HAVE_SYNC_FETCH */

typedef struct {
  unsigned int nb;
  pthread_t mt;
  HR_T begin, end, sup_nrel, sup_npri;
} ti_t;
static ti_t *ti;

/* Hightest criticality for performance. I inline all myself. */
static void
onepass_thread_singleton_removal (ti_t *mti)
{
  HR_T *tab, i;
  HC_T *o;
  bv_t j;

  mti->sup_nrel = mti->sup_npri = 0;
  for (i = mti->begin; i < mti->end; i++) {
    j = (((bv_t) 1) << (i & (BV_BITS - 1)));
    if (rel_used->p[i>>LN2_BV_BITS] & j)
      for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++)
	if (H.hc[*tab] == 1) {
	  for (tab = rel_compact[i]; *tab != UMAX(*tab); tab++) {
	    o = &(H.hc[*tab]);
	    ASSERT(*o);
	    if (*o < UMAX(*o) && !__sync_sub_and_fetch(o, 1))
	      (mti->sup_npri)++;
	  }
	  /* rel_compact[i] = NULL; */
	  rel_used->p[i>>LN2_BV_BITS] &= ~j;
 	  (mti->sup_nrel)++;
	  break;
	}
  }
  pthread_exit(NULL);
}

static void
onepass_singleton_parallel_removal (unsigned int nb_thread)
{
  pthread_attr_t attr;
  HR_T pas, incr;
  unsigned int i;
  int err;

  SMALLOC(ti, nb_thread, "onepass_singleton_parallel_removal :");
  ti[0].begin = 0;
  pas = (nrelmax / nb_thread) & ((HR_T) ~(BV_BITS -1));
  incr = 0;
  for (i = 0, incr = 0; i < nb_thread - 1; ) {
    incr += pas;
    ti[i].nb = i;
    ti[i++].end = incr;
    ti[i].begin = incr;
  }
  ti[i].nb = i;
  ti[i].end = nrelmax;
  pthread_attr_init (&attr);
  pthread_attr_setstacksize (&attr, 1<<16);
  pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
  for (i = 0; i < nb_thread; i++)
    if ((err = pthread_create (&ti[i].mt, &attr, (void *) onepass_thread_singleton_removal, &ti[i])))
    {
    fprintf (stderr, "onepass_singleton_parallel_removal : pthread_create error 1: %d. %s\n", err, strerror (errno));
    exit (1);
    }
  for (i = 0; i < nb_thread; i++) {
    pthread_join(ti[i].mt, NULL);
    newnrel -= ti[i].sup_nrel;
    newnprimes -= ti[i].sup_npri;
  }
  pthread_attr_destroy(&attr);
  SFREE(ti);
}
#endif /* ifdef HAVE_SYNC_FETCH */

static void
remove_singletons (unsigned int npass, double required_excess)
{
  HR_T oldnewnrel = 0, oldtmpnewnrel = 0;
#if HR == 32
  int32_t oldexcess = 0, excess;
#else
  int64_t oldexcess = 0, excess;
#endif
  int count = 0;

  SMALLOC(sum, H.hm, "remove_singletons 1");
  if (!binfilerel) newnrel = nrel;
  newnprimes = nprimes;
  excess = ((long) newnrel) - newnprimes;
  for ( ; newnrel != oldnewnrel || excess > (long) keep ; ) {
    /* delete heavy rows when we have reached a fixed point */
    if (newnrel == oldnewnrel) {
      /* check we have enough excess initially (at least required_excess) */
      if (count++ == 0 && (double) excess < required_excess*(double)newnprimes)
        {
          fprintf(stderr, "excess < %.2f * #primes. See -required_excess "
                          "argument.\n", required_excess);
          exit (1);
        }
      if (oldexcess > excess)
	fprintf (stderr, "   [each excess row deleted %2.2lf rows]\n",
		 (double) (oldtmpnewnrel - newnrel) / (double) (oldexcess - excess));
      oldexcess = excess;
      oldtmpnewnrel = newnrel;
      deleteHeavierRows (npass);
    }
    oldnewnrel = newnrel;
#ifdef HAVE_SYNC_FETCH
    ASSERT(npt);
    onepass_singleton_parallel_removal(npt);
#else
    onepass_singleton_removal();
#endif
    excess = ((long) newnrel) - newnprimes;
    fprintf (stderr, "   new_nrows=%lu new_ncols=%lu (%ld) at %2.2lf\n",
	     (unsigned long) newnrel, (unsigned long) newnprimes, (long) excess, seconds());
  }
  /* Warning: we might have empty rows, i.e., empty rel_compact[i] lists,
     if all primes in a relation are less than minpr or minpa. */
  nrel = newnrel;
  nprimes = newnprimes;
  SFREE(sum);
}

/* This function renumbers used primes (those with H->hashcount[i] > 1)
   and puts the corresponding index in H->hashcount[i].

   At return, nprimes is the number of used primes.

   We locate used primes and do not try to do fancy things as sorting w.r.t.
   weight, since this will probably be done later on.
   All rows will be 1 more that needed -> subtract 1 in fprint...! */
static void
renumber (const char *sos)
{
    FILE *fsos = NULL;
    HR_T i, nb = 1; /* we start at 1 here, but subtract 1 in fprint_rel_row */

    SMALLOC(H.hr, H.hm, "renumber 1");

    if (sos != NULL)
      {
	fprintf (stderr, "Output renumber table in file %s\n", sos);
	fsos = gzip_open (sos, "w");
        fprintf (fsos, "# each row contains 3 hexadecimal values: i p r\n"
		 "# i is the ideal index value (starting from 0)\n"
		 "# p is the corresponding prime\n"
		 "# r is the corresponding root (p+1 on the rational side)\n");
      }
    for (i = 0; i < H.hm; i++)
      if (H.hc[i]) {
	/* Since we consider only primes >= minpr or minpa,
	   smaller primes might appear only once here, thus we can't
	   assert H->hashcount[i] > 1, but H->hashcount[i] = 1 should
	   be rare if minpr/minpa are well chosen (not too large). */
	static int count = 0;
	if (H.hc[i] == 1 && (count ++ < 10))
	  {
	    if (GET_HASH_P(&H,i) == GET_HASH_R(&H,i) - 1)
	      fprintf (stderr, "Warning: singleton rational prime %lu\n",
		       (unsigned long) GET_HASH_P(&H,i));
	    else
	      fprintf (stderr, "Warning: singleton algebraic ideal (%lu,%lu)\n",
		       (unsigned long) GET_HASH_P(&H,i), (unsigned long) GET_HASH_R(&H,i));
	  }
	if (fsos)
	  fprintf(fsos, "%lx %lx %lx\n", (unsigned long) (nb - 1),
		  (unsigned long) GET_HASH_P(&H,i), (unsigned long) GET_HASH_R(&H,i));
	H.hr[i] = nb++;
      }
      else {
	H.hc[i] = UMAX(*(H.hc)); /* for getHashAddrAux */
	H.hr[i] = UMAX(*(H.hr));
      }
    if (fsos)
      gzip_close (fsos, sos);
    nb--;
    newnprimes = nb;
}

static void
prempt_load (prempt_t prempt_data) {
  char **p_files = prempt_data->files, *pprod;
  FILE *f;
  size_t try_load, load;
  char *p, *l, *pmax = &(prempt_data->buf[PREMPT_BUF]);

  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);

  while (*p_files) {
    if (!(f = popen (*p_files, "r"))) {
      fprintf (stderr, "prempt_load: popen error. %s\n", strerror (errno));
      exit (1);
    }
    p = strchr(*p_files, '/');
    if (p) {
      *p = 0;
      fprintf (stderr, "%s\n", *p_files);
      *p = '/';
      for ( ; ; ) {
	l = strchr(p, ' ');
	if (l) {
	  *l = 0;
	  fprintf(stderr, "   %-70s\n", p);
	  *l = ' ';
	  p = strchr(&(l[1]), '/');
	  if (!p) {
	    fprintf(stderr, "%s\n", &(l[1]));
	    break;
	  }
	}
	else {
	  fprintf(stderr, "   %-70s\n", p);
	  break;
	}
      }
    }
    pprod = (char *) prempt_data->pprod;
    for ( ; ; )
      if ((pprod != prempt_data->pcons) &&
	  (((((PREMPT_BUF + ((size_t) prempt_data->pcons))) - ((size_t) pprod)) &
	    (PREMPT_BUF - 1)) <= PREMPT_ONE_READ)) {
	nanosleep(&wait_classical, NULL);
      }
      else {
	try_load = MIN((PREMPT_BUF + ((size_t)pmin)) - ((size_t) pprod), PREMPT_ONE_READ);
	if ((load = fread (pprod, 1, try_load, f))) {
	  pprod += load;
	  if (pprod == pmax) pprod = pmin;
	  prempt_data->pprod = pprod;
	}
	else
	  if (feof(f)) {
	    pclose(f);
	    free(*p_files);
	    *p_files++ = NULL;
	    break;
	  }
	  else {
	    fprintf (stderr, "prempt_load: load error (%s) from\n%s\n", strerror (errno), *p_files);
	    exit (1);
	  }
      }
  }
  prempt_data->end = 1;
  pthread_exit (NULL);
}

static inline void
#ifndef FOR_FFS
relation_stream_get_fast (prempt_t prempt_data, unsigned int j)
#else
relation_stream_get_fast (prempt_t prempt_data, unsigned int j, int passtwo)
#endif
{
  buf_rel_t *mybufrel = &(buf_rel[j]);
  int64_t n;
  char *p;
  unsigned int k, i;
  unsigned long pr;
  unsigned char c, v;
  unsigned int ltmp;
#ifdef FOR_FFS
  unsigned int basis_ab = 16;
#else
  unsigned int basis_ab = 10;
#endif

#define LOAD_ONE(P) { c = *P; P = ((size_t) (P - pminlessone) & (PREMPT_BUF - 1)) + pmin; }

  p = (char *) prempt_data->pcons;

  LOAD_ONE(p);
  if (c == '-') {
    mybufrel->rel.a = -1;
    LOAD_ONE(p);
  }
  else
    mybufrel->rel.a = 1;
  for (n = 0 ; (v = ugly[c]) < basis_ab ; ) {
#ifdef FOR_FFS
    n = (n << 4) + v;
#else
    n = n * 10 + v;
#endif
    LOAD_ONE(p);
  }
  ASSERT_ALWAYS(c == ',');
  mybufrel->rel.a *= n;

  n = 0;
  LOAD_ONE(p);
  for ( ; (v = ugly[c]) < basis_ab ; ) {
#ifdef FOR_FFS
    n = (n << 4) + v;
#else
    n = n * 10 + v;
#endif
    LOAD_ONE(p);
  }
  ASSERT_ALWAYS(c == ':');
  mybufrel->rel.b = n;

  if (!rs->parse_only_ab) {
    /* Do something if we're also interested in primes */
    for ( k = 0, c = 0 ; ; ) {
    next_rat:
      if (c == ':') break;
      LOAD_ONE(p);
      for (pr = 0 ; (v = ugly[c]) < 16 ; ) {
	pr = (pr << 4) + v;
	LOAD_ONE(p);
      }
      ASSERT_ALWAYS(c == ',' || c == ':');
      /*
      if (pr)
        {
      */
	for (i = k; i--; )
	  {
	  if (mybufrel->rel.rp[i].p == pr)
	    {
#ifdef FOR_FFS
	      /* FOR_FFS, all .e != 0 must set to one in pass 1.
		 In passtwo, all .e must be keeped. */
	      if (passtwo)
#endif
		mybufrel->rel.rp[i].e++;
	      goto next_rat;
	    }
	  }
	if (mybufrel->rel.nb_rp_alloc == k) {
	  mybufrel->rel.nb_rp_alloc += mybufrel->rel.nb_rp_alloc >> 1;
	  if (k == NB_RATP_OPT) {
	    rat_prime_t *p = mybufrel->rel.rp;
	    SMALLOC(mybufrel->rel.rp, mybufrel->rel.nb_rp_alloc, "relation_stream_get_fast 1");
	    memcpy (mybufrel->rel.rp, p, NB_RATP_OPT * sizeof(rat_prime_t));
	  }
	  else
	    mybufrel->rel.rp = (rat_prime_t *)
	      realloc (mybufrel->rel.rp, mybufrel->rel.nb_rp_alloc * sizeof(rat_prime_t));
	}
	mybufrel->rel.rp[k++] = (rat_prime_t) { .p = pr, .e = 1};
	/*
      }
	*/
    }
    mybufrel->rel.nb_rp = k;

    for ( k = 0 ; ; ) {
    next_alg:
      if (c == '\n') break;
      LOAD_ONE(p);
      for (pr = 0 ; (v = ugly[c]) < 16 ; )
	{
	  pr = (pr << 4) + v;
	  LOAD_ONE(p);
	}
      ASSERT_ALWAYS(c == ',' || c == '\n');
      /*
      if (pr)
	{
      */
	for (i = k; i--; )
	  {
	    if (mybufrel->rel.ap[i].p == pr)
	      {
#ifdef FOR_FFS
		if (passtwo)
#endif
		  mybufrel->rel.ap[i].e++;
		
		goto next_alg;
	      }
	  }
	if (mybufrel->rel.nb_ap_alloc == k) {
	  mybufrel->rel.nb_ap_alloc += mybufrel->rel.nb_ap_alloc >> 1;
	  if (k == NB_ALGP_OPT) {
	    alg_prime_t *p = mybufrel->rel.ap;
	    SMALLOC(mybufrel->rel.ap, mybufrel->rel.nb_ap_alloc, "relation_stream_get_fast 2");
	    memcpy (mybufrel->rel.ap, p, NB_ALGP_OPT * sizeof(alg_prime_t));
	  }
	  else
	    mybufrel->rel.ap = (alg_prime_t *)
	      realloc (mybufrel->rel.ap, mybufrel->rel.nb_ap_alloc * sizeof(alg_prime_t));
	  }
	mybufrel->rel.ap[k++] = (alg_prime_t) { .p = pr,.r = -1, .e = 1};
	/*
	}
	*/
    }
    mybufrel->rel.nb_ap = k;

    if (mybufrel->rel.b > 0)
      {
	ltmp = 1;
#ifndef FOR_FFS
	for (k = 0, i = 0; i < mybufrel->rel.nb_rp; i++)
	  {
	    if (mybufrel->rel.rp[i].e & 1)
		{
		  mybufrel->rel.rp[k].p = mybufrel->rel.rp[i].p;
		  mybufrel->rel.rp[k].e = 1;
		  ltmp += ((HT_T) mybufrel->rel.rp[k++].p >= minpr);
		}
	  }
	mybufrel->rel.nb_rp = k;
	for (k = 0, i = 0; i < mybufrel->rel.nb_ap; i++)
	  {
	    if (mybufrel->rel.ap[i].e & 1)
		{
		  /* rel.ap[k].r will be computed later with another threads */
		  mybufrel->rel.ap[k].p = mybufrel->rel.ap[i].p;
		  mybufrel->rel.ap[k].e = 1;
		  ltmp += ((HT_T) mybufrel->rel.ap[k++].p >= minpa);
		}
	  }
	mybufrel->rel.nb_ap = k;
#else
	for (i = mybufrel->rel.nb_rp; i-- ;)
    {
	  ltmp += ((HT_T) mybufrel->rel.rp[i].p >= minpr);
#if defined FOR_FFS && defined STAT_FFS
	    if (passtwo && bit_vector_getbit(rel_used, (size_t) buf_rel[j].num))
      {
        if (abs(mybufrel->rel.rp[i].e) > 10)
          __stat_count[0]++;
        else
          __stat_count[abs(mybufrel->rel.rp[i].e)]++;
        __stat_nonzero++;
      }
#endif
	  }
  for (i = mybufrel->rel.nb_ap; i-- ;)
    {
	    ltmp += ((HT_T) mybufrel->rel.ap[i].p >= minpa);
#if defined FOR_FFS && defined STAT_FFS
	    if (passtwo && bit_vector_getbit(rel_used, (size_t) buf_rel[j].num))
      {
        if (abs(mybufrel->rel.ap[i].e) > 10)
          __stat_count[0]++;
        else
          __stat_count[abs(mybufrel->rel.ap[i].e)]++;
        __stat_nonzero++;
      }
#endif
    }
#endif
      }
    else
      {
	ltmp = ((HT_T) mybufrel->rel.a >= minpr) + 1;
	if ((HT_T) mybufrel->rel.a >= minpa)
	  ltmp += mybufrel->rel.nb_ap;
      }
    mybufrel->ltmp = ltmp;
  }
}

/* We don't use memory barrier nor (pre)processor orders for portability.
   So, if a ring buffer is empty, one problem exists if the producter
   produces one and the consumer takes it immediatly: the slot might be not
   complety written.
   Same problem exists when the buffer is full in the other sense.
   It's NOT a bug code, but the instructions reordonnancing of the
   optimiser compiler, which in the case of an empty buffer, may
   increase the producter counter BEFORE the end of the complete
   writing of the slot.
   The only << solution >> without synchronous barrier or (pre)processor
   order is a nanosleep after the end of the waiting loop.
   * empty buffer :
       if (A == B + 1) nanosleep (&wait_classical, NULL);
   * full buffer :
       if (A + 1 == B + SIZEBUF) nanosleep (&wait_classical, NULL);
   It's very dirty!
*/

void insertRelation() {
  unsigned int j;
  unsigned long cpy_cpt_rel_b;

  /*
  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);
  */
  cpy_cpt_rel_b = cpt_rel_b;
  for ( ; ; ) {
    while (cpt_rel_a == cpy_cpt_rel_b)
      if (!end_insertRelation) nanosleep (&wait_classical, NULL);
      else
	if (cpt_rel_a == cpy_cpt_rel_b) pthread_exit(NULL);
    j = (unsigned int) (cpy_cpt_rel_b & (T_REL - 1));
    if (cpt_rel_a == cpy_cpt_rel_b + 1) nanosleep (&wait_classical, NULL);
    if (bit_vector_getbit(rel_used, (size_t) buf_rel[j].num)) {
      if (buf_rel[j].rel.b)
	insertNormalRelation (j);
      else
	insertFreeRelation (j);
    }
    if (relation_stream_disp_progress_now_p(rs))
      fprintf(stderr, "read useful %lu relations in %.1fs"
	      " -- %.1f MB/s -- %.1f rels/s\n",
	      rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
    cpy_cpt_rel_b++;
    cpt_rel_b = cpy_cpt_rel_b;
  }
}

static int
no_singleton(buf_rel_t *br) {
  HR_T *phk, *end_phk;

  for (phk = br->hk, end_phk = &(phk[br->lhk]); phk != end_phk; phk++)
    if (H.hc[*phk] == 1) {
      relsup++;
      for (phk = br->hk, end_phk = &(phk[br->lhk]); phk != end_phk; phk++)
	if (H.hc[*phk] != UMAX(H.hc[*phk]) && !(--(H.hc[*phk]))) prisup++;
      return 0;
    }
  return 1;
}


static void
printrel() {
  unsigned int j, aff;
  unsigned long cpy_cpt_rel_b;

  /*
  pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
  pthread_setcanceltype(PTHREAD_CANCEL_ASYNCHRONOUS, NULL);
  */
  cpy_cpt_rel_b = cpt_rel_b;
  for ( ; ; ) {
    while (cpt_rel_a == cpy_cpt_rel_b)
      if (!end_insertRelation)
	nanosleep (&wait_classical, NULL);
      else
	if (cpt_rel_a == cpy_cpt_rel_b)
	  pthread_exit(NULL);
    j = (unsigned int) (cpy_cpt_rel_b & (T_REL - 1));
    if (cpt_rel_a == cpy_cpt_rel_b + 1) nanosleep (&wait_classical, NULL);
    aff = bit_vector_getbit(rel_used, (size_t) buf_rel[j].num);
    if (boutfilerel && aff) {
      aff = no_singleton(&(buf_rel[j]));
      if (!aff)
	bit_vector_clearbit(rel_used, (size_t) buf_rel[j].num);
    }
    if (aff) {
      if (raw)
	fprint_relation_raw(ofile, &(buf_rel[j].rel));
      else
	if (!boutfilerel)
	  W += (double) fprint_rel_row(ofile, &(buf_rel[j]));
    }
#ifdef FOR_FFS
    else
      if (!boutfilerel)
	fprint_relation_raw(ofile2, &(buf_rel[j].rel));
#endif
    if (relation_stream_disp_progress_now_p(rs))
      fprintf(stderr, "re-read & print useful %lu relations in %.1fs"
	      " -- %.1f MB/s -- %.1f rels/s\n",
	      rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
    cpy_cpt_rel_b++;
    cpt_rel_b = cpy_cpt_rel_b;
  }
}

static HR_T *buf_rel_new_hk(unsigned int j, unsigned int t) {
  if (buf_rel[j].mhk < t) {
    if (buf_rel[j].mhk != NB_HK_OPT) SFREE(buf_rel[j].hk);
    buf_rel[j].mhk = t;
    SMALLOC(buf_rel[j].hk, buf_rel[j].mhk, "buf_rel_new_hk");
  }
  return(buf_rel[j].hk);
}

static void threadfindroot(fr_t *mfr) {
  HR_T *phk;
  buf_rel_t *mybufrel;
  unsigned int i, j;

  for (;;)
    switch(mfr->ok) {
    case 0:
      nanosleep (&wait_classical, NULL);
      break;
    case 1 :
      for (j = mfr->num; j <= mfr->end; j++)
	{
	  mybufrel = &(buf_rel[j]);
	  if
#ifdef FOR_FFS
	    (!boutfilerel || bit_vector_getbit(rel_used, (size_t) mybufrel->num))
#else
	    (bit_vector_getbit(rel_used, (size_t) mybufrel->num))
#endif
	      {
		if (mybufrel->rel.b)
		  {
		    phk = buf_rel_new_hk(j, mybufrel->rel.nb_rp + mybufrel->rel.nb_ap);
		    for (i = 0; i < mybufrel->rel.nb_rp; i++)
		      if (!boutfilerel || mybufrel->rel.rp[i].p >= minpr)
			*phk++ = HKM(mybufrel->rel.rp[i].p, mybufrel->rel.rp[i].p + 1, H.hm);
		    for (i = 0; i < mybufrel->rel.nb_ap; i++)
		      if (!boutfilerel || mybufrel->rel.ap[i].p >= minpa) {
			mybufrel->rel.ap[i].r = (HT_T)
#ifdef FOR_FFS
			  findroot_ffs (mybufrel->rel.a, mybufrel->rel.b, mybufrel->rel.ap[i].p);
#else
			findroot (mybufrel->rel.a, mybufrel->rel.b, mybufrel->rel.ap[i].p);
#endif
			if (mybufrel->rel.ap[i].r == UMAX(mybufrel->rel.ap[i].r))
			  mybufrel->rel.ap[i].r = mybufrel->rel.ap[i].p;
			*phk++ = HKM(mybufrel->rel.ap[i].p, mybufrel->rel.ap[i].r, H.hm);
		      }
		  }
		else
		  {
		    phk = buf_rel_new_hk(j, mybufrel->rel.nb_ap + 1);
		    if (!boutfilerel || (uint64_t) mybufrel->rel.a >= minpr)
		      *phk++ = HKM(mybufrel->rel.a, mybufrel->rel.a + 1, H.hm);
		    if (!boutfilerel || (uint64_t) mybufrel->rel.a >= minpa)
		      for (i = 0; i < mybufrel->rel.nb_ap; i++)
			*phk++ = HKM(mybufrel->rel.a, mybufrel->rel.ap[i].p, H.hm);
		  }
		mybufrel->lhk = phk - mybufrel->hk;
	      }
	  else
	    mybufrel->lhk = 0;
	}
      mfr->ok = 0;
      break;
    case 2:
      mfr->ok = 3;
      pthread_exit(NULL);
    }
}

#define FIND_PR_IN_H							\
  h = HKM(pr.p, pr.r, H.hm);						\
  p = &(H.ht[h]);							\
  while (p->p != pr.p || p->r != pr.r) 					\
    if (++p == ep) p = H.ht;						\
  *phk++ = h = (HT_T) (p - H.ht)

static void
threadfindroot_exactphk(fr_t *mfr) {
  HR_T *phk;
  buf_rel_t *mybufrel;
  ht_t *p, *ep;
  ht_t pr;
  HT_T h;
  unsigned int i, j;

  ep = &(H.ht[H.hm]);
  for (;;)
    switch(mfr->ok) {
    case 0:
      nanosleep (&wait_classical, NULL);
      break;
    case 1 :
      for (j = mfr->num; j <= mfr->end; j++)
	{
	  mybufrel = &(buf_rel[j]);
	  if (bit_vector_getbit(rel_used, (size_t) mybufrel->num)) {
	    if (mybufrel->rel.b)
	      {
		phk = buf_rel_new_hk(j, mybufrel->rel.nb_rp + mybufrel->rel.nb_ap);
		for (i = 0; i < mybufrel->rel.nb_rp; i++)
		  if (mybufrel->rel.rp[i].p >= minpr) {
		    pr = (ht_t) { (HT_T) mybufrel->rel.rp[i].p, (HT_T) (mybufrel->rel.rp[i].p + 1) };
		    FIND_PR_IN_H;
		  }
		for (i = 0; i < mybufrel->rel.nb_ap; i++)
		  if (mybufrel->rel.ap[i].p >= minpa) {
#ifndef FOR_FFS
		    mybufrel->rel.ap[i].r = (HT_T)
		      findroot (mybufrel->rel.a, mybufrel->rel.b, mybufrel->rel.ap[i].p);
#else
		    mybufrel->rel.ap[i].r = (HT_T)
		      findroot_ffs (mybufrel->rel.a, mybufrel->rel.b, mybufrel->rel.ap[i].p);
#endif
		    if (mybufrel->rel.ap[i].r == UMAX(mybufrel->rel.ap[i].r))
		      mybufrel->rel.ap[i].r = mybufrel->rel.ap[i].p;
		    pr = (ht_t) { (HT_T) mybufrel->rel.ap[i].p, (HT_T) mybufrel->rel.ap[i].r };
		    FIND_PR_IN_H;
		  }
	      }
	    else
	      {
		phk = buf_rel_new_hk(j, mybufrel->rel.nb_ap + 1);
		if ((uint64_t) mybufrel->rel.a >= minpr) {
		  pr = (ht_t) { (HT_T) mybufrel->rel.a, (HT_T) (mybufrel->rel.a + 1) };
		  FIND_PR_IN_H;
		}
		if ((uint64_t) mybufrel->rel.a >= minpa)
		  for (i = 0; i < mybufrel->rel.nb_ap; i++) {
		    pr = (ht_t) { (HT_T) mybufrel->rel.a, (HT_T) mybufrel->rel.ap[i].p };
		    FIND_PR_IN_H;
		  }
	     }
	    mybufrel->lhk = phk - mybufrel->hk;
	  }
	  else
	    mybufrel->lhk = 0;
	}
      mfr->ok = 0;
      break;
    case 2:
      mfr->ok = 3;
      pthread_exit(NULL);
    }
}

/* Read all relations from file, and fills the rel_used and rel_compact arrays
   for each relation i:
   - rel_used[i] = 0 if the relation i is deleted
     rel_used[i] = 1 if the relation i is kept (so far)
   - rel_compact is an array, terminated by -1, of pointers to the entries
     in the hash table for the considered primes

     Trick: we only read relations for which rel_used[i]==1.
 */
static int
prempt_scan_relations_pass_one ()
{
  char *pcons, *pcons_old, *pcons_max, *p, **ff;
  pthread_attr_t attr;
  pthread_t thread_load, thread_relation, thread_fr[(1<<NFR)];
  fr_t fr[(1<<NFR)];
  prempt_t prempt_data;
  unsigned long cpy_cpt_rel_a;
  unsigned int length_line, i, k;
  int err;
  char c;

  memset (fr, 0, (1<<NFR) * sizeof(*fr));
  end_insertRelation = 0;

  SMALLOC (buf_rel, T_REL, "prempt_scan_relations_pass_one 1");
  MEMSETZERO(buf_rel, T_REL);
  SMALLOC (buf_rel_data, T_REL, "prempt_scan_relations_pass_one 2");
  for (i = T_REL; i--; ) {
    buf_rel[i].rel.rp = buf_rel_data[i].rp;
    buf_rel[i].rel.nb_rp_alloc = NB_RATP_OPT;
    buf_rel[i].rel.ap = buf_rel_data[i].ap;
    buf_rel[i].rel.nb_ap_alloc = NB_ALGP_OPT;
    buf_rel[i].hk = buf_rel_data[i].hk;
    buf_rel[i].mhk = NB_HK_OPT;
  }

  cpt_rel_a = cpt_rel_b = 0;
  cpy_cpt_rel_a = cpt_rel_a;
  nprimes = 0;
  relation_stream_init (rs);
  rs->pipe = 1;
  length_line = 0;

  prempt_data->files = prempt_open_compressed_rs (rep_cado, fic);

  SMALLOC (prempt_data->buf, PREMPT_BUF, "prempt_scan_relations_pass_one 3");
  pmin = prempt_data->buf;
  pminlessone = pmin - 1;
  prempt_data->pcons = pmin;
  prempt_data->pprod = pmin;
  pcons_max = &(prempt_data->buf[PREMPT_BUF]);
  prempt_data->end = 0;
  pthread_attr_init (&attr);
  pthread_attr_setstacksize (&attr, 1<<16);
  pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
  if ((err = pthread_create (&thread_load, &attr, (void *) prempt_load, prempt_data)))
    {
    fprintf (stderr, "prempt_scan_relations_pass_one: pthread_create error 1: %d. %s\n", err, strerror (errno));
    exit (1);
    }
  if ((err = pthread_create (&thread_relation, &attr, (void *) insertRelation, NULL)))
    {
    fprintf (stderr, "prempt_scan_relations_pass_one: pthread_create error 2: %d. %s\n", err, strerror (errno));
    exit (1);
    }
  for (i = 0; i < (1<<NFR); i++)
    if ((err = pthread_create (&(thread_fr[i]), &attr, (void *) threadfindroot, &(fr[i]))))
      {
    fprintf (stderr, "prempt_scan_relations_pass_one: pthread_create error 3: %d. %s\n", err, strerror (errno));
    exit (1);
      }

  pcons = (char *) prempt_data->pcons;
  for ( ; ; )
    {
      rs->pos += length_line;
      length_line = 0;
      prempt_data->pcons = pcons;

      while (pcons == prempt_data->pprod)
	if (!prempt_data->end)
	  nanosleep (&wait_classical, NULL);
	else
	  if (pcons == prempt_data->pprod)
	    goto end_of_files;
      if (pcons == prempt_data->pprod + sizeof(*pcons)) nanosleep (&wait_classical, NULL);

      rs->lnum++;
      if (*pcons != '#')
	{
	  while ((((PREMPT_BUF + ((size_t) prempt_data->pprod)) -
		   ((size_t) pcons)) & (PREMPT_BUF - 1)) <= ((unsigned int) RELATION_MAX_BYTES) &&
		 !prempt_data->end)
	    nanosleep(&wait_classical, NULL);

	  if (pcons > prempt_data->pprod)
	    {
	      p = &(pcons_max[-1]);
	      c = *p;
	      *p = '\n';
	      pcons_old = pcons;
	      while (*pcons++ != '\n');
	      *p = c;
	      length_line = (pcons - pcons_old);
	      if (pcons <= p)
		goto testendline;
	      pcons = pmin;
	      if (c == '\n')
		goto testendline;
	    }
	  p = &(((char *) prempt_data->pprod)[-1]);
	  c = *p;
	  *p = '\n';
	  pcons_old = pcons;
	  while (*pcons++ != '\n');
	  *p = c;
	  length_line += (pcons - pcons_old);
	testendline:
	  if (c != '\n' && pcons == prempt_data->pprod)
	    {
	      fprintf (stderr, "prempt_scan_relations_pass_one: "
		       "the last line has not a carriage return\n");
	      exit(1);
	    }
	  if (length_line > ((unsigned int) RELATION_MAX_BYTES))
	    {
	      fprintf (stderr, "prempt_scan_relations_pass_one: relation line size (%u) is "
		       "greater than RELATION_MAX_BYTES (%d)\n",
		       length_line, RELATION_MAX_BYTES);
	      exit(1);
	    }
	
	  while (cpy_cpt_rel_a == cpt_rel_b + T_REL) nanosleep(&wait_classical, NULL);
	  k = (unsigned int) (cpy_cpt_rel_a & (T_REL - 1));
	  if (cpy_cpt_rel_a + 1 == cpt_rel_b + T_REL) nanosleep(&wait_classical, NULL);
	  buf_rel[k].num = rs->nrels++;

#ifndef FOR_FFS
	  if (bit_vector_getbit(rel_used, (size_t) buf_rel[k].num))
	    relation_stream_get_fast (prempt_data, k);
#else
	  relation_stream_get_fast (prempt_data, k, 0);
#endif
	  /* Delayed find root computation by block of 1<<NNFR */
	  if (cpy_cpt_rel_a && !(k & ((1<<NNFR)-1)))
	    {
	    i = (k>>NNFR) & ((1<<NFR)-1);
	    while (fr[i].ok) nanosleep(&wait_classical, NULL);
	    if (k)
	      {
		fr[i].num = k - (1<<NNFR);
		fr[i].end = k - 1;
	      }
	    else
	      {
		fr[i].num = T_REL - (1<<NNFR);
		fr[i].end = T_REL - 1;
	      }
	    fr[i].ok = 1;
	    if (cpy_cpt_rel_a > (1<<(NFR+NNFR)))
	      cpt_rel_a = cpy_cpt_rel_a - (1<<(NFR+NNFR));
	    }
	  cpy_cpt_rel_a++;
	}
      else
	{
	  do
	    {
	      while (pcons == prempt_data->pprod)
		{
		  if (!prempt_data->end)
		    nanosleep (&wait_classical, NULL);
		  else
		    if (pcons == prempt_data->pprod)
		      {
			fprintf (stderr, "prempt_scan_relations_pass_one: at the end of"
				 " files, a line without \\n ?\n");
			exit (1);
		      }
		}
	      p = ((pcons <= prempt_data->pprod) ? (char *) prempt_data->pprod
		   : pcons_max) - 1;
	      c = *p;
	      *p = '\n';
	      pcons_old = pcons;
	      while (*pcons++ != '\n');
	      *p = c;
	      length_line += (pcons - pcons_old);
	      err = (pcons > p && c != '\n');
	      if (pcons == pcons_max) pcons = pmin;
	    }
	  while (err);
	}
    }

 end_of_files:
  if (cpy_cpt_rel_a) {
    k = (unsigned int) ((cpy_cpt_rel_a - 1) & (T_REL - 1));
    if (k & ((1<<NNFR)-1)) {
      i = ((k>>NNFR)+1) & ((1<<NFR)-1);
      while (fr[i].ok) nanosleep(&wait_classical, NULL);
      fr[i].num = k & ~((1<<NNFR)-1);
      fr[i].end = k;
      fr[i].ok = 1;
    }
  }
  for (i = 0; i < (1<<NFR); i++) {
    while (fr[i].ok) nanosleep(&wait_classical, NULL);
    fr[i].ok = 2;
    pthread_join(thread_fr[i], NULL);
  }
  cpt_rel_a = cpy_cpt_rel_a;
  while (cpy_cpt_rel_a != cpt_rel_b)
    nanosleep(&wait_classical, NULL);

  end_insertRelation = 1;
  pthread_join(thread_relation, NULL);
  /* if (pthread_tryjoin_np (thread_load, NULL)) */
  pthread_cancel(thread_load);
  pthread_join(thread_load, NULL);
  pthread_attr_destroy(&attr);

  free (prempt_data->buf);
  for (ff = prempt_data->files; *ff; free(*ff++));
  free (prempt_data->files);
  for (i = T_REL; i-- ; ) {
    if (buf_rel[i].rel.nb_rp_alloc != NB_RATP_OPT ) SFREE(buf_rel[i].rel.rp);
    if (buf_rel[i].rel.nb_ap_alloc != NB_ALGP_OPT)  SFREE(buf_rel[i].rel.ap);
    if (buf_rel[i].mhk             != NB_HK_OPT   ) SFREE(buf_rel[i].hk);
  }
  SFREE (buf_rel);
  SFREE (buf_rel_data);

  relation_stream_trigger_disp_progress(rs);
  fprintf (stderr, "End of read: %lu relations in %.1fs -- %.1f MB/s -- %.1f rels/s\n",
           rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
  relation_stream_clear(rs);

  if (rs->nrels != nrelmax) {
    fprintf (stderr, "Error, -nrels value should match the number of scanned relations\nexpected %lu relations, found %lu\n", (unsigned long) nrelmax, rs->nrels);
    exit (EXIT_FAILURE);
  }

  return 1;
}

/* ReRead all relations from files and output remaining relations
   (those with rel_used[i] <> 0) in ofile.

   For FFS ouput deleted relations (those with rel_used[i] == 0) in ofile2

   If raw is non-zero, output relations in CADO format
   (otherwise in format used by merge).

   Raw or not raw, a ring buffer is used to stock the relations.
   It maybe useless in some cases.
 */
static int
prempt_scan_relations_pass_two (const char *oname,
#ifdef FOR_FFS
				const char *oname2,
#else
				bit_vector_srcptr rel_used,
#endif
				HR_T nrows, HR_T ncols, int raw)
{
  char *pcons, *pcons_old, *pcons_max, *p, **f;
  pthread_attr_t attr;
  pthread_t thread_load, thread_printrel, thread_fr[(1<<NFR)];
  fr_t fr[(1<<NFR)];
  prempt_t prempt_data;
  unsigned long cpy_cpt_rel_a;
  unsigned int length_line, i, k;
  int err;
  char c;

  int pipe;

  ofile = fopen_compressed_w(oname, &pipe, NULL);
#ifdef FOR_FFS
  ofile2 = fopen_compressed_w(oname2, &pipe_2, NULL);
#endif
  if (!raw)
    fprintf (ofile, "%lu %lu\n", (unsigned long) nrows, (unsigned long) ncols);
  fprintf (stderr, "Final pass:\n");
  MEMSETZERO(fr, 1<<NFR);
  end_insertRelation = 0;

  SMALLOC(buf_rel, T_REL, "prempt_scan_relations_pass_two 1");
  MEMSETZERO(buf_rel, T_REL);
  SMALLOC (buf_rel_data, T_REL, "prempt_scan_relations_pass_one 2");
  for (i = T_REL; i--; ) {
    buf_rel[i].rel.rp = buf_rel_data[i].rp;
    buf_rel[i].rel.nb_rp_alloc = NB_RATP_OPT;
    buf_rel[i].rel.ap = buf_rel_data[i].ap;
    buf_rel[i].rel.nb_ap_alloc = NB_ALGP_OPT;
    buf_rel[i].hk = buf_rel_data[i].hk;
    buf_rel[i].mhk = NB_HK_OPT;
  }

  cpt_rel_a = cpt_rel_b = 0;
  cpy_cpt_rel_a = cpt_rel_a;
  relation_stream_init (rs);
  rs->pipe = 1;
  length_line = 0;

  prempt_data->files = prempt_open_compressed_rs (rep_cado, fic);

  SMALLOC(prempt_data->buf, PREMPT_BUF, "prempt_scan_relations_pass_two 3");
  pmin = prempt_data->buf;
  pminlessone = pmin - 1;
  prempt_data->pcons = pmin;
  prempt_data->pprod = pmin;
  pcons_max = &(prempt_data->buf[PREMPT_BUF]);
  prempt_data->end = 0;
  pthread_attr_init (&attr);
  pthread_attr_setstacksize (&attr, 1<<16);
  pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);

  if ((err = pthread_create (&thread_load, &attr, (void *) prempt_load, prempt_data)))
    {
    fprintf (stderr, "prempt_scan_relations_pass_two: pthread_create error 1: %d. %s\n", err, strerror (errno));
    exit (1);
    }
  W = 0.0;
  if ((err = pthread_create (&thread_printrel, &attr, (void *) printrel, NULL)))
    {
    fprintf (stderr, "prempt_scan_relations_pass_two: pthread_create error 2: %d. %s\n", err, strerror (errno));
    exit (1);
    }
  for (i = 0; i < (1<<NFR); i++)
    if ((err = pthread_create (&(thread_fr[i]), &attr, (void *) (boutfilerel ? threadfindroot_exactphk : threadfindroot), &(fr[i]))))
      {
    fprintf (stderr, "prempt_scan_relations_pass_two: pthread_create error 3: %d. %s\n", err, strerror (errno));
    exit (1);
      }
  pcons = (char *) prempt_data->pcons;
  for ( ; ; )
    {
      rs->pos += length_line;
      length_line = 0;
      prempt_data->pcons = pcons;

      while (pcons == prempt_data->pprod)
	if (!prempt_data->end)
	  nanosleep (&wait_classical, NULL);
	else
	  if (pcons == prempt_data->pprod)
	    goto end_of_files;
      if (pcons == prempt_data->pprod + sizeof(*pcons)) nanosleep (&wait_classical, NULL);

      rs->lnum++;
      if (*pcons != '#')
	{
	  while ((((PREMPT_BUF + ((size_t) prempt_data->pprod)) -
		   ((size_t) pcons)) & (PREMPT_BUF - 1)) <= ((unsigned int) RELATION_MAX_BYTES) &&
		 !prempt_data->end)
	      nanosleep(&wait_classical, NULL);

	  if (pcons > prempt_data->pprod)
	    {
	      p = &(pcons_max[-1]);
	      c = *p;
	      *p = '\n';
	      pcons_old = pcons;
	      while (*pcons++ != '\n');
	      *p = c;
	      length_line = (pcons - pcons_old);
	      if (pcons <= p)
		goto lineok;
	      pcons = pmin;
	      if (c == '\n')
		goto lineok;
	    }
	  p = &(((char *) prempt_data->pprod)[-1]);
	  c = *p;
	  *p = '\n';
	  pcons_old = pcons;
	  while (*pcons++ != '\n');
	  *p = c;
	  length_line += (pcons - pcons_old);
	lineok:
	  while (cpy_cpt_rel_a == cpt_rel_b + T_REL) nanosleep(&wait_classical, NULL);
	  k = (unsigned int) (cpy_cpt_rel_a & (T_REL - 1));
	  if (cpy_cpt_rel_a + 1 == cpt_rel_b + T_REL) nanosleep(&wait_classical, NULL);
	  buf_rel[k].num = rs->nrels++;

#ifndef	FOR_FFS
	  if (bit_vector_getbit(rel_used, (size_t) buf_rel[k].num))
	    relation_stream_get_fast (prempt_data, k);
#else		
	  relation_stream_get_fast (prempt_data, k, 1);
#endif
	  /* Delayed find root computation by block of 1<<NNFR */
	  if (cpy_cpt_rel_a && !(k & ((1<<NNFR)-1)))
	    {
	      i = (k>>NNFR) & ((1<<NFR)-1);
	      while (fr[i].ok) nanosleep(&wait_classical, NULL);
	      if (k)
		{
		  fr[i].num = k - (1<<NNFR);
		  fr[i].end = k - 1;
		}
	      else
		{
		  fr[i].num = T_REL - (1<<NNFR);
		  fr[i].end = T_REL - 1;
		}
	      fr[i].ok = 1;
	      if (cpy_cpt_rel_a > (1<<(NFR+NNFR)))
		cpt_rel_a = cpy_cpt_rel_a - (1<<(NFR+NNFR));
	    }
	  cpy_cpt_rel_a++;
	}
      else
	{
	  do
	    {
	      while (pcons == prempt_data->pprod)
		{
		  if (!prempt_data->end)
		    nanosleep (&wait_classical, NULL);
		}
	      p = ((pcons <= prempt_data->pprod) ? (char *) prempt_data->pprod
		   : pcons_max) - 1;
	      c = *p;
	      *p = '\n';
	      pcons_old = pcons;
	      while (*pcons++ != '\n');
	      *p = c;
	      length_line += (pcons - pcons_old);
	      err = (pcons > p && c != '\n');
	      if (pcons == pcons_max) pcons = pmin;
	    }
	  while (err);
	}
    }

 end_of_files:
  if (cpy_cpt_rel_a) {
    k = (unsigned int) ((cpy_cpt_rel_a - 1) & (T_REL - 1));
    if (k & ((1<<NNFR)-1)) {
      i = ((k>>NNFR)+1) & ((1<<NFR)-1);
      while (fr[i].ok) nanosleep(&wait_classical, NULL);
      fr[i].num = k & ~((1<<NNFR)-1);
      fr[i].end = k;
      fr[i].ok = 1;
    }
  }
  for (i = 0; i < (1<<NFR); i++) {
    while (fr[i].ok) nanosleep(&wait_classical, NULL);
    fr[i].ok = 2;
    pthread_join(thread_fr[i], NULL);
  }
  cpt_rel_a = cpy_cpt_rel_a;
  while (cpy_cpt_rel_a != cpt_rel_b)
    nanosleep(&wait_classical, NULL);

  end_insertRelation = 1;
  pthread_join(thread_printrel, NULL);
  /* if (pthread_tryjoin_np (thread_load, NULL)) */
  pthread_cancel(thread_load);
  pthread_join(thread_load, NULL);
  pthread_attr_destroy(&attr);

  free (prempt_data->buf);
  for (f = prempt_data->files; *f; free(*f++));
  free (prempt_data->files);
  for (i = T_REL; i-- ; ) {
    if (buf_rel[i].rel.nb_rp_alloc != NB_RATP_OPT ) SFREE(buf_rel[i].rel.rp);
    if (buf_rel[i].rel.nb_ap_alloc != NB_ALGP_OPT)  SFREE(buf_rel[i].rel.ap);
    if (buf_rel[i].mhk             != NB_HK_OPT   ) SFREE(buf_rel[i].hk);
  }
  SFREE (buf_rel);
  SFREE (buf_rel_data);

  relation_stream_trigger_disp_progress(rs);
  fprintf (stderr, "End of re-read: %lu relations in %.1fs -- %.1f MB/s -- %.1f rels/s\n",
           rs->nrels, rs->dt, rs->mb_s, rs->rels_s);

#if defined FOR_FFS && defined STAT_FFS
  fprintf (stderr, "# of non zero coeff: %lu\n", __stat_nonzero);
  for (int i = 1; i <= 10 ; i++)
    fprintf (stderr, "# of coeffs of abs value %d: %lu(%.2f%%)\n", i,
             __stat_count[i], 100 * (double) __stat_count[i]/__stat_nonzero);
  fprintf (stderr, "# of coeffs of abs value > 10: %lu(%.2f%%)\n",
           __stat_count[0], 100 * (double) __stat_count[0]/__stat_nonzero);
#endif

  /* write excess to stdout */
  if (!raw)
    printf ("NROWS:%lu WEIGHT:%1.0f WEIGHT*NROWS=%1.2e\n",
                 (unsigned long) nrows, W, W * (double) nrows);
  printf ("EXCESS: %lu\n", ((long) nrows) - ncols);
  fflush (stdout);

  relation_stream_clear(rs);

  if (pipe)  pclose(ofile);  else fclose(ofile);
#ifdef FOR_FFS
  if (pipe_2) pclose(ofile2); else fclose(ofile2);
#endif

  return 1;
}

static void
usage (void)
{
  fprintf (stderr, "Usage: purge [options] -poly polyfile -out purgedfile -nrels nnn [-basepath <path>] [-subdirlist <sl>] [-filelist <fl>] file1 ... filen\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "       -keep    nnn - prune if excess > nnn (default 160)\n");
  fprintf (stderr, "       -minpa   nnn - purge alg. primes >= nnn (default alim)\n");
  fprintf (stderr, "       -minpr   nnn - purge rat. primes >= nnn (default rlim)\n");
  fprintf (stderr, "       -nprimes nnn - expected number of prime ideals\n");
  fprintf (stderr, "       -sos sosfile - to keep track of the renumbering\n");
  fprintf (stderr, "       -raw         - output relations in CADO format\n");
  fprintf (stderr, "       -npthr   nnn - threads number for suppress singletons\n");
  fprintf (stderr, "       -inprel  file_rel_used : load active relations\n");
  fprintf (stderr, "       -outrel  file_rel_used : write active relations\n");
  fprintf (stderr, "       -npass   nnn - number of step of clique removal (default %d)\n", DEFAULT_NPASS);
  fprintf (stderr, "       -required_excess nnn - percentage of excess required at the end of the first singleton removal step (default %.2f)\n",
  DEFAULT_REQUIRED_EXCESS);
#ifdef FOR_FFS
  fprintf (stderr, "       -outdel file - output file for deleted relations\n");
#endif
  exit (1);
}

#ifndef FOR_FFS
/* estimate the number of primes <= B */
static uint64_t
approx_phi (long B)
{
  return (B <= 1) ? 0 : (uint64_t) ((double) B / log ((double) B));
}
#else
/* estimate the number of ideals of degree <= B */
static int
approx_ffs (int d)
{
#ifdef USE_F2
  ASSERT_ALWAYS(d <= 36); /* otherwise the result is > 2^31 */
  if (d <= 0)
    return 0;
  else
  /* for d >= 9, between 1 and 10% greater than the real value */
    return (int) ((double) 1.12 * pow (2.0, (double) (d + 1 - log(d)/log(2))));
#elif USE_F3
  ASSERT_ALWAYS(d <= 22); /* otherwise the result is > 2^31 */
  if (d <= 0)
    return 0;
  else
  /* for d >= 6, between 1 and 10% greater than the real value */
    return (int) ((double) 1.05 * pow (3.0, (double) (d+0.4-log(d)/log(3))));
#else
  ASSERT_ALWAYS(0);
#endif
}
#endif

static void
set_rep_cado (char *argv0) {
  char *p;

  strcpy(rep_cado, argv0);
  p = strrchr(rep_cado, '/');
  if (p)
    strcpy (&(p[1]), "../");
  else
    strcat(rep_cado, "../");
}

int
main (int argc, char **argv)
{
  int k;
  size_t mysize;
  param_list pl;
  unsigned int npass = DEFAULT_NPASS;
  double required_excess = DEFAULT_REQUIRED_EXCESS;

  set_rep_cado(argv[0]);
  wct0 = wct_seconds ();
  fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
  for (k = 1; k < argc; k++)
    fprintf (stderr, " %s", argv[k]);
  fprintf (stderr, "\n");

  param_list_init(pl);

  param_list_configure_switch(pl, "raw", &raw);

  argv++,argc--;

  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) continue;
    /* Since we accept file names freeform, we decide to never abort
     * on unrecognized options */
    if (!strcmp(*argv, "--help")) usage();
    break;
  }

#if HR == 32
  param_list_parse_uint(pl, "nrels", (unsigned int *) &nrelmax);
  param_list_parse_uint(pl, "nprimes", (unsigned int *) &nprimes);
  param_list_parse_uint(pl, "keep", (unsigned int *) &keep);
#else
  param_list_parse_uint64(pl, "nrels", (uint64_t *) &nrelmax);
  param_list_parse_uint64(pl, "nprimes", (uint64_t *) &nprimes);
  param_list_parse_uint64(pl, "keep", (uint64_t *) &keep);
#endif

#if HT == 32
  param_list_parse_uint(pl, "minpr", (unsigned int *) &minpr);
  param_list_parse_uint(pl, "minpa", (unsigned int *) &minpa);
#else
  param_list_parse_uint64(pl, "minpr", (uint64_t *) &minpr);
  param_list_parse_uint64(pl, "minpa", (uint64_t *) &minpa);
#endif

  /* param_list_parse_uint(pl, "npthr", (unsigned int *) &npt); */
  const char * snpt = param_list_lookup_string(pl, "npthr");
  if (snpt) {
    char *p;
    if ((p = strchr(snpt, 'x'))) {
      unsigned int x, y;
      *p = 0;
      if (sscanf(snpt, "%u", &x) && sscanf(&p[1], "%u", &y))
	npt = x * y;
      else
	usage();
    } else
      if (!sscanf(snpt, "%u", &npt))
	usage();
  }
  param_list_parse_uint(pl, "npass", &npass);
  param_list_parse_double(pl, "required_excess", &required_excess);
  const char * filelist = param_list_lookup_string(pl, "filelist");
  const char * basepath = param_list_lookup_string(pl, "basepath");
  const char * subdirlist = param_list_lookup_string(pl, "subdirlist");
  const char * purgedname = param_list_lookup_string(pl, "out");
  const char * sos = param_list_lookup_string(pl, "sos");
  const char * infilerel = param_list_lookup_string(pl, "inrel");
  const char * outfilerel = param_list_lookup_string(pl, "outrel");
#ifdef FOR_FFS
  const char * deletedname = param_list_lookup_string(pl, "outdel");
#endif
  binfilerel = (infilerel != 0);
  boutfilerel = (outfilerel != 0);

  cado_poly_init (pol);

  const char * tmp;

  ASSERT_ALWAYS((tmp = param_list_lookup_string(pl, "poly")) != NULL);
#ifndef FOR_FFS
  cado_poly_read(pol, tmp);
#else
  ffs_poly_read(pol, tmp);
#endif

  if (param_list_warn_unused(pl)) {
    usage();
  }

  if ((basepath || subdirlist) && !filelist) {
    fprintf(stderr, "-basepath / -subdirlist only valid with -filelist\n");
    usage();
  }

  if (nrelmax == 0)
    {
      fprintf (stderr, "Error, missing -nrels ... option (or nrels=0)\n");
      usage ();
    }

#ifdef FOR_FFS /* For FFS we need to remember the renumbering of primes*/
  if (sos == NULL)
    {
      fprintf (stderr, "Error, missing -sos option.\n");
      exit(1);
    }
  if (deletedname == NULL)
    {
      fprintf (stderr, "Error, missing -outdel option.\n");
      exit(1);
    }
#endif

  /* the current code assumes the "rational" side has degree 1 */
#ifndef FOR_FFS
  ASSERT_ALWAYS(pol->rat->degree == 1);
#endif

  /* On a 32-bit computer, even 1 << 32 would overflow. Well, we could set
     map[ra] = 2^32-1 in that case, but not sure we want to support 32-bit
     primes on a 32-bit computer... */
#ifdef FOR_FFS
#ifdef USE_F2
  need64 = (pol->rat->lpb >= 32) || (pol->alg->lpb >= 32);
#elif USE_F3
  need64 = (pol->rat->lpb >= 16) || (pol->alg->lpb >= 16);
#else
  ASSERT_ALWAYS(0);
#endif
#else
  need64 = (pol->rat->lpb > 32) || (pol->alg->lpb > 32);
#endif

  if (need64 && sizeof (HT_T) < 8)
    {
      fprintf (stderr, "Error, too large LPBs for a 32-bit program\n");
      exit(1);
    }

  if (minpr == UMAX(minpr)) minpr = pol->rat->lim;
  if (minpa == UMAX(minpa)) minpa = pol->alg->lim;

  fprintf (stderr, "Weight function used during clique removal:\n"
                   "  0     1     2     3     4     5     6     7\n"
                   "%0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f %0.3f\n",
                   weight_function_clique(0), weight_function_clique(1),
                   weight_function_clique(2), weight_function_clique(3),
                   weight_function_clique(4), weight_function_clique(5),
                   weight_function_clique(6), weight_function_clique(7));

  fprintf (stderr, "Number of relations is %lu\n", (unsigned long) nrelmax);
  if (nprimes > 0) Hsize = nprimes;
  else
    {
      /* Estimating the number of needed primes (remember that hashInit
         multiplies by a factor 1.5). */
#ifndef FOR_FFS
      Hsizer = approx_phi (1L << pol->rat->lpb);
      Hsizea = approx_phi (1L << pol->alg->lpb);
#else
      Hsizer = approx_ffs (pol->rat->lpb);
      Hsizea = approx_ffs (pol->alg->lpb);
#endif
      Hsize = Hsizer + Hsizea;
    }
  fprintf (stderr, "Estimated number of prime ideals: %lu\n", (unsigned long) Hsize);
  tot_alloc0 = H.hm * (sizeof(HC_T) + sizeof(ht_t));

  mysize = (nrelmax + 7) >> 3;
  bit_vector_init(rel_used, nrelmax);
  if (!rel_used->p) {
    fprintf (stderr, "Main: malloc error (%lu bytes).\n", (unsigned long) mysize);
    exit (1);
  }
  if (binfilerel) {
    FILE *in;
    void *p, *pl, *pf;
    int pipe;

    fprintf (stderr, "Loading rel_used file %s, %lu bytes\n",
	     infilerel, (unsigned long) mysize);
    if (!(in = fopen_compressed_r (infilerel, &pipe, NULL))) {
      fprintf (stderr, "Purge main: rel_used file %s cannot be read.\n", infilerel);
      exit (1);
    }
    newnrel = 0;
    for (p = (void *) rel_used->p, pf = p + mysize; (pf - p); ) {
      pl = p + fread(p, 1, pf - p, in);
      while (p != pl) newnrel += nbbits[*((uint8_t *) p++)];
      if ((p == pf) ^ !feof(in)) {
	fprintf (stderr, "Purge main: rel_used file length is incorrect versus nrels.\n");
	exit (1);
      }
    }
    fclose(in);
    fprintf (stderr, "Number of used relations is %lu\n", (unsigned long) newnrel);
  }
  else
    bit_vector_set(rel_used, 1);

  if (nrelmax & (BV_BITS - 1))
    rel_used->p[nrelmax>>LN2_BV_BITS] &= (((bv_t) 1)<<(nrelmax & (BV_BITS - 1))) - 1;

  tot_alloc0 += mysize;
  fprintf (stderr, "Allocated rel_used of %uMb (total %zuMb so far)\n",
	   nrelmax >> 20, tot_alloc0 >> 20);

  if (!boutfilerel) {
    SMALLOC(rel_compact, nrelmax, "main 1");
    SMALLOC(rel_weight, nrelmax, "main 2");
  tot_alloc0 += nrelmax * (sizeof (HR_T *) + sizeof (HC_T));
  /* %zu is the C99 modifier for size_t */
  fprintf (stderr, "Allocated rel_compact of %zu MB (total %zu MB so far)\n",
	   ((size_t) nrelmax * sizeof (HR_T *)) >> 20, tot_alloc0 >> 20);
  }
  /* Build the file list (ugly). It is the concatenation of all
   *  b s p
   * where:
   *    b is the basepath (empty if not given)
   *    s ranges over all subdirs listed in the subdirlist (empty if no
   *    such list)
   *    p ranges over all paths listed in the filelist.
   *
   * If files are provided directly on the command line, the basepath
   * and subdirlist arguments are ignored.
   */

  if (!filelist) {
    fic = argv;
  } else if (!subdirlist) {
    fic = filelist_from_file (basepath, filelist, 0);
  } else {
    /* count the number of files in the filelist */
    int nfiles = 0;
    int nsubdirs = 0;
    char ** fl = filelist_from_file (NULL, filelist, 0);
    for(char ** p = fl ; *p ; p++, nfiles++);

    char ** sl = filelist_from_file (basepath, subdirlist, 1);
    for(char ** p = sl ; *p ; p++, nsubdirs++);

    SMALLOC(fic, nsubdirs * nfiles + 1, "main 3");
    char ** full = fic;
    for(char ** f = fl ; *f ; f++) {
      for(char ** s = sl ; *s ; s++, full++) {
	int ret = asprintf(full, "%s/%s", *s, *f);
	ASSERT_ALWAYS(ret >= 0);
      }
    }
    *full=NULL;
    filelist_clear(fl);
    filelist_clear(sl);
  }

  nrel = nrelmax;

  /************************** first pass *************************************/

  tot_alloc = tot_alloc0;

#ifndef FOR_FFS
  fprintf (stderr, "Pass 1, filtering ideals >= %lu on rat. side and "
           "%lu on alg. side:\n", (unsigned long) minpr, (unsigned long) minpa);
#else
  fprintf (stderr, "Pass 1, filtering ideals of degree >= %lu on rat. side "
           "and %lu on alg. side:\n", (unsigned long) minpr, (unsigned long) minpa);
#ifdef USE_F2
  minpr = 1 << minpr;
  minpa = 1 << minpa;
#elif USE_F3
  minpr = 1 << (2*minpr);
  minpa = 1 << (2*minpa);
#else
  ASSERT_ALWAYS(0);
#endif
#endif

  hashInit (&H, Hsize, 1);

  prempt_scan_relations_pass_one ();

  if (!binfilerel)
    fprintf (stderr, "   nrels=%lu, nprimes=%lu; excess=%ld\n",
	     (unsigned long) nrel, (unsigned long) nprimes, ((long) nrel) - nprimes);
  else
    fprintf (stderr, "   nrels=%lu, nrels used=%lu, nprimes=%lu; excess=%ld\n",
	     (unsigned long) nrel, (unsigned long) newnrel, (unsigned long) nprimes, ((long) newnrel) - nprimes);

  if (!boutfilerel) {
    remove_singletons (npass, required_excess);
    fprintf (stderr, "   nrel=%lu, nprimes=%lu; excess=%ld\n",
	     (unsigned long) nrel, (unsigned long) nprimes, ((long) nrel) - nprimes);
    if (nrel <= nprimes) /* covers case nrel = nprimes = 0 */
      {
	fprintf(stderr, "number of relations <= number of ideals\n");
	exit (1);
      }
  }
  hashCheck (&H);

  if (!boutfilerel) {
    my_malloc_free_all ();

    fprintf (stderr, "Freeing rel_compact array...\n");
    /* we do not use it anymore */
    free (rel_compact);
    free (rel_weight);

    /*************************** second pass ***********************************/

    /* we renumber the primes in order of apparition in the hashtable */
    renumber (sos);

    /* reread the relation files and convert them to the new coding */
    fprintf (stderr, "Storing remaining relations...\n");
  }
  relsup = 0;
  prisup = 0;

#ifndef FOR_FFS
  /* reread (purgedname, fic, rel_used, nrel_new, nprimes_new, raw); */
  prempt_scan_relations_pass_two (purgedname, rel_used, nrel, nprimes, raw);
#else
  /* reread (purgedname, deletedname, fic, rel_used, nrel_new, nprimes_new,                                                                raw); */
  prempt_scan_relations_pass_two (purgedname, deletedname, nrel, nprimes, raw);
#endif

  if (boutfilerel) {
    FILE *out;
    int pipe;
    void *p, *pf;

    if (nrelmax & (BV_BITS - 1))
      rel_used->p[nrelmax>>LN2_BV_BITS] &= (((bv_t) 1)<<(nrelmax & (BV_BITS - 1))) - 1;

    fprintf (stderr, "Relations with at least one singleton found and suppress:%lu\nNumber of primes suppress : %lu\nWriting rel_used file %s, %lu bytes\n",
	     (unsigned long) relsup, (unsigned long) prisup, outfilerel, (unsigned long) mysize);
    if (!(out = fopen_compressed_w (outfilerel, &pipe, NULL))) {
      fprintf (stderr, "Purge main: rel_used file %s cannot be written.\n", outfilerel);
      exit (1);
    }
    for (p = (void *) rel_used->p, pf = p + mysize; p != pf;
	 p += fwrite(p, 1, pf - p, out));
    fclose (out);
  }

  hashFree (&H);
  bit_vector_clear(rel_used);
  cado_poly_clear (pol);

  if (filelist) filelist_clear(fic);

  param_list_clear(pl);

  print_timing_and_memory (wct0);

  return 0;
}
