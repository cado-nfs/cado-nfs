/* purge --- remove singletons

Copyright 2008 Francois Morain, Paul Zimmermann

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
  This program works in several passes over the relation files (currently 3):
  - the first pass considers only rational primes >= minpr and algebraic
    primes >= minpa. This means that only those primes are loaded in memory,
    and only relations with singletons in that range can be removed. By
    default minpr and minpa are taken as rlim and alim respectively.
    Usually more than half of the relations are discarded in that pass.
  - the second pass takes minpr=minpa=0, i.e., all primes are now considered
    (but only for relations surviving the first pass). It removes the remaining
    relations having singletons (typically very few).
    Also, if the excess exceeds 'keep', it prunes the heavier relations.
  - the final pass goes through the relations again, and dumps the remaining
    ones in the format needed by 'merge'.

  Remark: more passes can be implemented to decrease the memory usage in the
  first one. Simply increase the initial minpr/minpa, and decrease them at
  each pass. The only requirement is that the pass before the final one
  must have minpr=minpa=0, to ensure no singleton remains for 'merge'.

  This program uses the following data structures:
  rel_used[i]    - non-zero iff relation i is kept (so far)
  rel_compact[i] - list of 'h' indices in H table of considered (p,r) for row i
                   (terminated by a -1 sentinel)
  h = getHashAddr (H, p, r) - index of prime ideal (p, r) in hash table H
                              (rational primes use r = -2)
  GET_HASH_P(H,h) - prime corresponding to index h
  GET_HASH_R(H,h) - root  corresponding to index h (-2 for rational prime)
  H->hashcount[h] - number of occurrences of (p, r) in current relations
                   (in the last pass, H->hashcount[] stores the index of the
                    prime ideal).
*/

#include "cado.h"
#include <gmp.h>
#include "mod_ul.c"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "utils.h"

#include "hashpair.h"
#include "gzip.h"

#define MAX_FILES 1000000

/********************** own memory allocation routines ***********************/

/* Rationale: calling one malloc() for each relation in insertNormalRelation()
   is expensive, since malloc() allocates some extra information to keep track
   of every memory blocks. Instead, we allocate memory in big blocks of size
   BLOCK_SIZE. */

#define BLOCK_SIZE 1000000 /* memory blocks are allocated of that # of int's */

/* relcompact_list is a list of blocks, each one of BLOCK_SIZE ints */
int **relcompact_list = NULL;
int relcompact_current = -1; /* index of current block */
unsigned long relcompact_used = BLOCK_SIZE; /* usage of current block */

/* return a pointer to an array of n ints */
static int*
my_malloc_int (unsigned long n)
{
  int *ptr;

  if (relcompact_used + n > BLOCK_SIZE)
    {
#if 0
        /* NEVER EVER rely on the fact that a shrinking something by
         * realloc does not change the pointer !
         */
      /* first shrink current block */
      if (relcompact_current >= 0)
        relcompact_list[relcompact_current] = (int*) realloc (relcompact_list[relcompact_current], relcompact_used * sizeof (int));
#endif
      /* allocate a new block */
      relcompact_current ++;
      relcompact_list = (int**) realloc (relcompact_list, (relcompact_current + 1) * sizeof (int*));
      relcompact_list[relcompact_current] = (int*) malloc (BLOCK_SIZE * sizeof (int));
      relcompact_used = 0;
    }
  ptr = relcompact_list[relcompact_current] + relcompact_used;
  relcompact_used += n;
  return ptr;
}

static void
my_malloc_free_all (void)
{
  while (relcompact_current >= 0)
    free (relcompact_list[relcompact_current--]);
  free (relcompact_list);
  relcompact_list = NULL;
  relcompact_used = BLOCK_SIZE;
}

/*****************************************************************************/

/* dirty trick to distinguish rational primes: we store -2 for their root */
static unsigned long minus2;

/* Adds in table_ind[] the indices of the rational primes in 'rel'.
   All primes are assumed to be different (and to appear with an odd exponent).
   nb_coeff is the current index in table_ind[] (input and output).
 */
static void
fprint_rat (int *table_ind, int *nb_coeff, relation_t *rel, hashtable_t *H)
{
    int i, nbc = *nb_coeff;

    for (i = 0; i < rel->nb_rp; i++)
      table_ind[nbc++] = H->hashcount[getHashAddr (H, rel->rp[i].p, minus2)];
    *nb_coeff = nbc;
}

/* Adds in table_ind[] the indices of the algebraic primes in 'rel'.
   All primes are assumed to be different (and to appear with an odd exponent).
   nb_coeff is the current index in table_ind[] (input and output).
 */
static void
fprint_alg (int *table_ind, int *nb_coeff, relation_t rel, hashtable_t *H)
{
    int i, nbc = *nb_coeff;

    for (i = 0; i < rel.nb_ap; i++)
      table_ind[nbc++] = H->hashcount[getHashAddr (H, rel.ap[i].p,
                                                      rel.ap[i].r)];
    *nb_coeff = nbc;
}

/* Adds a free relation in table_ind[]: a is the corresponding prime */
static void
fprint_free (int *table_ind, int *nb_coeff, relation_t rel, hashtable_t *H)
{
    long p = rel.a;
    int i, nbc = *nb_coeff, index;

    index = getHashAddr (H, p, minus2);
    ASSERT(H->hashcount[index] >= 0);
    table_ind[nbc++] = H->hashcount[index];
    for(i = 0; i < rel.nb_ap; i++)
      {
	index = getHashAddr (H, p, rel.ap[i].p);
        ASSERT_ALWAYS(H->hashcount[index] >= 0);
        table_ind[nbc++] = H->hashcount[index];
    }
    *nb_coeff = nbc;
}

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
*/
static int
fprint_rel_row (FILE *file, int irel, relation_t rel, hashtable_t *H)
{
  int i;
  int *table_ind;
  int nb_coeff;

  table_ind = (int*) malloc ((rel.nb_rp + rel.nb_ap) * sizeof (int));

  nb_coeff = 0;

  if (rel.b == 0) /* free relation */
      fprint_free (table_ind, &nb_coeff, rel, H);
  else
    {
      /* adds rational primes in table_ind */
      fprint_rat (table_ind, &nb_coeff, &rel, H);

      /* adds algebraic primes in table_ind */
      fprint_alg (table_ind, &nb_coeff, rel, H);
    }

  fprintf (file, "%d %"PRId64" %"PRIu64" %d",
          irel, rel.a, rel.b, nb_coeff);
  for (i = 0; i < nb_coeff; ++i)
    /* due to the +1 in renumber */
    fprintf (file, " " PURGE_INT_FORMAT, table_ind[i] - 1);
  fprintf (file, "\n");

  free (table_ind);

  return nb_coeff;
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
static void
insertNormalRelation (int **rel_compact, int irel,
                      int *nprimes, hashtable_t *H, relation_t *rel,
                      unsigned long minpr, unsigned long minpa,
                      unsigned long *tot_alloc)
{
    int *tmp = NULL, ltmp = 0, i, j, h;

    reduce_exponents_mod2 (rel);

    /* first count number of "large" primes */
    for (j = 0; j < rel->nb_rp; j++)
      ltmp += (rel->rp[j].p >= minpr); /* only consider primes >= minpr */
    for (j = 0; j < rel->nb_ap; j++)
      ltmp += (rel->ap[j].p >= minpa); /* only consider primes >= minpa */

    /* ltmp is the number of considered primes in the relation.
       We might have ltmp=0 if all primes are less than minpr, maxpr. */
    tmp = my_malloc_int (ltmp + 1);
    *tot_alloc += (ltmp + 1) * sizeof (int);
    i = 0; /* number of entries in tmp */

    /* trick: we use the same hash table for rational and algebraic
       primes, but use a fake root -2 for rational primes, which
       ensures there is no collision with algebraic primes */
    for (j = 0; j < rel->nb_rp; j++)
      if (rel->rp[j].p >= minpr)
        {
          /* we only insert 'considered' primes, since we cannot get
             a correct count of the other ones when we delete a relation */
          h = hashInsert (H, rel->rp[j].p, minus2);
          *nprimes += (H->hashcount[h] == 1); /* new prime */
          tmp[i++] = h;
        }

    for (j = 0; j < rel->nb_ap; j++)
      if (rel->ap[j].p >= minpa)
        {
          /* FIXME: instead of computing r = a/b mod p, we could
             first search for a potential match (p,r) in the hash table,
             then check that a = r*b mod b, and only compute a/b mod p
             for new entries, thus trading a modular inverse by a modular
             multiplication for all but the first occurrences of that ideal */
          rel->ap[j].r = findroot (rel->a, rel->b, rel->ap[j].p);
          h = hashInsert (H, rel->ap[j].p, rel->ap[j].r);
          *nprimes += (H->hashcount[h] == 1); /* new prime */
          tmp[i++] = h;
        }

    tmp[i] = -1; /* sentinel */
    rel_compact[irel] = tmp;
}

/* The information is stored in the ap[].p part, which is odd, but convenient.
   rel->ap.p[0..deg[ contains the deg roots of f(x) mod p, leading to deg
   ideals for the factorization of (p). */
static void
insertFreeRelation (int **rel_compact, int irel, int *nprimes,
                    hashtable_t *H, relation_t *rel,
                    unsigned long minpr, unsigned long minpa,
                    unsigned long *tot_alloc)
{
    unsigned long p = rel->a; /* rel->b == 0 */
    int j, h, *tmp = NULL, itmp, ltmp;

    /* the prime on the rational side is rel->a
       the prime ideal on the algebraic side are (rel->a, rel->ap[j].p) */
    ltmp = 1 + (p >= minpr);
    if (p >= minpa)
      ltmp += rel->nb_ap;
    tmp = my_malloc_int (ltmp);
    *tot_alloc += ltmp * sizeof (int);
    itmp = 0;
    if (p >= minpr)
      {
        h = hashInsert (H, p, minus2);
        *nprimes += (H->hashcount[h] == 1); /* new prime */
        tmp[itmp++] = h;
      }
    for (j = 0; j < rel->nb_ap; j++)
      if (p >= minpa)
        {
          h = hashInsert(H, p, rel->ap[j].p);
          *nprimes += (H->hashcount[h] == 1); /* new ideal */
          tmp[itmp++] = h;
        }
    tmp[itmp++] = -1;
    rel_compact[irel] = tmp;
}

/* Read all relations from file, and fills the rel_used and rel_compact arrays
   for each relation i:
   - rel_used[i] = 0 if the relation i is deleted
     rel_used[i] = 1 if the relation i is kept (so far)
   - rel_compact is an array, terminated by -1, of pointers to the entries
     in the hash table for the considered primes

     Trick: we only read relations for which rel_used[i]=1.
 */
static int
scan_relations (char **ficname, int *nprimes, hashtable_t *H,
                bit_vector_srcptr rel_used, int nrelmax, int **rel_compact,
		long minpr, long minpa, unsigned long *tot_alloc)
{
    *nprimes = 0;

    ASSERT(rel_compact != NULL);

    relation_stream rs;
    relation_stream_init(rs);

    for( ; *ficname ; ficname++) {
	fprintf(stderr, "   %-70s\n", *ficname);
        relation_stream_openfile(rs, *ficname);
        for ( ; ; ) {
            int irel = rs->nrels;
            if (relation_stream_get(rs, NULL) < 0)
                break;
            ASSERT_ALWAYS(rs->nrels <= nrelmax);
            if (bit_vector_getbit(rel_used, irel)) {
                if (rs->rel.b > 0)
                    insertNormalRelation (rel_compact, irel, nprimes,
                            H, &(rs->rel), minpr, minpa, tot_alloc);
                else
                    insertFreeRelation (rel_compact, irel, nprimes,
                            H, &(rs->rel), minpr, minpa, tot_alloc);
            }
            if (!relation_stream_disp_progress_now_p(rs))
                continue;

            fprintf(stderr, 
                    "read %d relations in %.1fs"
                    " -- %.1f MB/s -- %.1f rels/s\n",
                    rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
        }
        relation_stream_closefile(rs);
    }
    fprintf(stderr, 
            "read %d relations in %.1fs -- %.1f MB/s -- %.1f rels/s\n",
            rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
    if (rs->nrels != nrelmax) {
        fprintf (stderr, "Error, -nrels value should match the number of scanned relations\n");
        fprintf (stderr, "expected %d relations, found %d\n", nrelmax,
                 rs->nrels);
        exit (EXIT_FAILURE);
    }

    relation_stream_clear(rs);

    return 1;
}

/* Return a non-zero value iff some prime (ideal) in the array tab[] is single
   (tab[j] is the index in the hash table of the corresponding prime ideal).
*/
static int
has_singleton (int *tab, hashtable_t *H)
{
    int j;

    for (j = 0; tab[j] != -1; j++)
      if (H->hashcount[tab[j]] == 1)
        return 1;
    return 0;
}

/* Delete a relation: set rel_used[i] to 0, update the count of primes
   in that relation, and set rel_compact[i] to NULL.
   Warning: we only update the count of primes that we consider, i.e.,
   rational primes >= minpr and algebraic primes >= minpa.
*/
static void
delete_relation (int i, int *nprimes, hashtable_t *H,
                 bit_vector_ptr rel_used, int **rel_compact)
{
  int j, *tab = rel_compact[i];

  for (j = 0; tab[j] != -1; j++)
    {
      H->hashcount[tab[j]]--; /* remove one occurrence of prime 'j' */
      *nprimes -= (H->hashcount[tab[j]] == 0);
    }
  rel_compact[i] = NULL;
  bit_vector_clearbit(rel_used, i);
}

/* New pruning code, which optimizes the decrease of N*W where N is the number
   of rows, and W is the total weight. We consider the connected
   components of the relation R(i1,i2) iff i1 and i2 share a prime
   of weight 2. If we remove one component of n rows and total weight w,
   then we save w*N+n*W (neglecting 2nd order terms), thus we remove
   first the components with the largest value of n/N + w/W. */

/* Define the weight of a connected component: those with the largest cost
   are removed first in the pruning; w (unsigned long) is the component weight,
   W (double) is the total matrix weight, n (int) is the number of rows of the
   component, and N (double) is the number of rows of the matrix. This macro
   must return a double. */
#if 0
/* optimize the decrease of W*N */
#define COST(w,W,n,N) ((double) (w) / (W) + (double) (n) / (N))
/* optimize the decrease of W */
#define COST(w,W,n,N) ((double) (w) / (W))
#endif

/* optimize the decrease of N */
#define COST(w,W,n,N) ((double) (n) / (N))

static int
compare (const void *v1, const void *v2)
{
    float w1 = *((float*) v1);
    float w2 = *((float*) v2);

    return (w1 >= w2) ? -1 : 1;
}

/* Compute connected component of row i for the relation R(i1,i2) if rows
   i1 and i2 share a prime of weight 2.
   Return number of rows of components, and put in w the total weight. */
static int
compute_connected_component (bit_vector_ptr T, uint32_t i, hashtable_t *H,
                             unsigned long *w, int **rel_compact,
                             uint32_t *sum)
{
  int j, h, n;
  uint32_t k;

  n = 1;         /* current row */
  bit_vector_setbit(T, i); /* mark row as visited */
  for (j = 0; (h = rel_compact[i][j]) != -1; j++)
    if (H->hashcount[h] == 2)
      {
        k = sum[h] - i; /* other row where prime of index h appears */
        if (!bit_vector_getbit(T, k)) /* row k was not visited yet */
          n += compute_connected_component (T, k, H, w, rel_compact, sum);
      }
  *w += j; /* add weight of current row */
  return n;
}

typedef struct {
  float w;
  uint32_t i;
} comp_t;

static void
deleteHeavierRows (hashtable_t *H, int *nrel, int *nprimes, bit_vector_ptr rel_used,
                   int **rel_compact, int nrelmax, int keep)
{
  uint32_t *sum; /* sum of row indices for primes with weight 2 */
  int i, j, h, n, ltmp = 0;
  unsigned long w;
  double W = 0.0; /* total matrix weight */
  double N = 0.0; /* numebr of rows */
  comp_t *tmp = NULL; /* (weight, index) */
  int target;

  if ((*nrel - *nprimes) <= keep)
    return;

  /* first collect sums for primes with weight 2, and compute total weight */
  sum = (uint32_t*) malloc (H->hashmod * sizeof (uint32_t));
  memset (sum, 0, H->hashmod * sizeof (uint32_t));
  for (i = 0; i < nrelmax; i++)
    if (bit_vector_getbit(rel_used, i)) {
        for (j = 0; (h = rel_compact[i][j]) != -1; j++)
          if (H->hashcount[h] == 2)
            sum[h] += i;
        N += 1.0;
        W += (double) j;
      }
  fprintf (stderr, "Matrix has %1.0f rows and weight %1.0f\n", N, W);
  ASSERT_ALWAYS(N == (double) *nrel);

  /* now initialize bit table for relations used */
  bit_vector T;
  bit_vector_init_set(T, nrelmax, 0);

  for (i = 0; i < nrelmax; i++)
    if (bit_vector_getbit(rel_used, i) && bit_vector_getbit(T,i) == 0)
      {
        w = 0;
        n = compute_connected_component (T, i, H, &w, rel_compact, sum);
        ltmp ++;
        tmp = (comp_t*) realloc (tmp, ltmp * sizeof (comp_t));
        tmp[ltmp - 1].w = (float) COST(w,W,n,N);
        tmp[ltmp - 1].i = i;
      }

  qsort (tmp, ltmp, sizeof(comp_t), compare);

  /* remove heaviest components, assuming each one decreases the excess by 1;
     we remove only half of the excess at each call of deleteHeavierRows,
     hoping to get "better" components to remove at the next call. */
  target = (*nrel - *nprimes + keep) / 2;
  for (i = 0; i < ltmp && (*nrel) - (*nprimes) > target; i ++)
    {
      /* it suffices to remove one relation per component, the other ones
         will be removed by onepass_singleton_removal */
      delete_relation (tmp[i].i, nprimes, H, rel_used, rel_compact);
      *nrel -= 1;
    }

  bit_vector_clear(T);
  free (sum);
  free (tmp);
}

static void
onepass_singleton_removal (int nrelmax, int *nrel, int *nprimes,
                           hashtable_t *H, bit_vector_ptr rel_used, int **rel_compact)
{
  int i;

  for (i = 0; i < nrelmax; i++)
    if (bit_vector_getbit(rel_used, i) && has_singleton (rel_compact[i], H))
      {
        delete_relation (i, nprimes, H, rel_used, rel_compact);
        *nrel -= 1;
      }
}

static void
remove_singletons (int *nrel, int nrelmax, int *nprimes, hashtable_t *H,
                   bit_vector_ptr rel_used, int **rel_compact, int keep, int final)
{
  int old, newnrel = *nrel, newnprimes = *nprimes;

    do{
	old = newnrel;
        if (final) /* otherwise the excess is wrong */
          deleteHeavierRows (H, &newnrel, &newnprimes, rel_used, rel_compact,
                             nrelmax, keep);
	if(newnrel != old)
	    fprintf (stderr, "deleted heavier relations: %d %d at %2.2lf\n",
                     newnrel, newnprimes, seconds ());
	onepass_singleton_removal (nrelmax, &newnrel, &newnprimes, H,
                                   rel_used, rel_compact);
	fprintf (stderr, "   new_nrows=%d new_ncols=%d (%d) at %2.2lf\n",
		newnrel, newnprimes, newnrel-newnprimes, seconds());
    } while(newnrel != old);

    /* Warning: we might have empty rows, i.e., empty rel_compact[i] lists,
       if all primes in a relation are less then minpr or minpa. */

    *nrel = newnrel;
    *nprimes = newnprimes;
}

/* This function renumbers used primes (those with H->hashcount[i] > 1)
   and puts the corresponding index in H->hashcount[i].

   At return, nprimes is the number of used primes.

   We locate used primes and do not try to do fancy things as sorting w.r.t.
   weight, since this will probably be done later on.
   All rows will be 1 more that needed -> subtract 1 in fprint...! */
static void
renumber (int *nprimes, hashtable_t *H, const char *sos)
{
    FILE *fsos = NULL;
    unsigned int i;
    int nb = 1; /* we start at 1 here, but subtract 1 in fprint_rel_row */

    if (sos != NULL)
      {
	fprintf (stderr, "Output renumber table in file %s\n", sos);
	fsos = gzip_open (sos, "w");
        fprintf (fsos, "# each row contains 3 hexadecimal values: i p r\n");
        fprintf (fsos, "# i is the ideal index value (starting from 0)\n");
        fprintf (fsos, "# p is the corresponding prime\n");
        fprintf (fsos, "# r is the corresponding root (-2=fffffffe on the rational side)\n");
      }
    for(i = 0; i < H->hashmod; i++)
      if (H->hashcount[i] == 0)
        H->hashcount[i] = -1;
      else
        {
          /* Since we consider only primes >= minpr or minpa,
             smaller primes might appear only once here, thus we can't
             assert H->hashcount[i] > 1, but H->hashcount[i] = 1 should
             be rare if minpr/minpa are well chosen (not too large). */
          static int count = 0;
          if (H->hashcount[i] == 1 && count ++ < 10)
            {
              if (GET_HASH_R(H,i) == minus2)
                fprintf (stderr, "Warning: singleton rational prime %"
                        PRIu64 "\n",
                         GET_HASH_P(H,i));
              else
                fprintf (stderr, "Warning: singleton algebraic ideal (%"
                        PRIu64",%"PRIu64")\n",
                         GET_HASH_P(H,i), GET_HASH_R(H,i));
            }
          H->hashcount[i] = nb++;
          if (fsos != NULL)
            fprintf(fsos, "%x %" PRIx64 " %" PRIx64 "\n",
                    H->hashcount[i] - 1, GET_HASH_P(H,i), GET_HASH_R(H,i));
	}
    if (fsos != NULL)
      gzip_close (fsos, sos);
    nb--;
    *nprimes = nb;
}

/* Read again the relation files ficname[0], ..., ficname[nbfic-1],
   and output remaining relations (those with rel_used[i] <> 0) in ofile.

   If raw is non-zero, output relations in CADO format
   (otherwise in format used by merge).
*/
static void
reread (const char *oname, char ** ficname,
        hashtable_t *H, bit_vector_srcptr rel_used, int nrows, int ncols, int raw)
{
  FILE *ofile;
  int ret, nr = 0;
  double W = 0.0; /* total weight */
  int pipe;

  ofile = fopen_compressed_w(oname, &pipe, NULL);
  if (raw == 0)
    fprintf (ofile, "%d %d\n", nrows, ncols);
  fprintf (stderr, "Final pass:\n");

  relation_stream rs;
  relation_stream_init(rs);

  for ( ; *ficname ; ficname++) {
      relation_stream_openfile(rs, *ficname);
      fprintf(stderr, "   %-70s\n", *ficname);
      for ( ; ; ) {
          int irel = rs->nrels;
          if (relation_stream_get(rs, NULL) < 0)
              break;
          // ASSERT_ALWAYS(rs->nrels <= nrelmax);
          if (bit_vector_getbit(rel_used, irel)) {
              if (raw == 0)
              {
                  if (rs->rel.b > 0)
                  {
                      reduce_exponents_mod2 (&rs->rel);
                      computeroots (&rs->rel);
                  }
                  W += (double) fprint_rel_row (ofile, irel, rs->rel, H);
              }
              else
                  fprint_relation_raw (ofile, rs->rel);
              nr++;
              if (nr >= nrows)
                  ret = 0; /* we are done */
          }
          if (!relation_stream_disp_progress_now_p(rs))
              continue;

          fprintf(stderr, 
                  "re-read %d relations in %.1fs"
                  " -- %.1f MB/s -- %.1f rels/s\n",
                  rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
      }
      relation_stream_closefile(rs);
  }
  relation_stream_trigger_disp_progress(rs);
  fprintf(stderr, 
          "re-read %d relations in %.1fs"
          " -- %.1f MB/s -- %.1f rels/s\n",
          rs->nrels, rs->dt, rs->mb_s, rs->rels_s);
  relation_stream_clear(rs);
  if (pipe) pclose(ofile); else fclose(ofile);

  /* write excess to stdout */
  if (raw == 0)
    printf ("WEIGHT: %1.0f WEIGHT*NROWS=%1.2e\n", W, W * (double) nr);
  printf ("EXCESS: %d\n", nrows - ncols);
}

static void
usage (void)
{
  fprintf (stderr, "Usage: purge [options] -poly polyfile -out purgedfile -nrels nnn [-basepath <path>] [-subdirlist <sl>] [-filelist <fl>] file1 ... filen\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, "       -excess  nnn - initial excess must be >= nnn (default 1)\n");
  fprintf (stderr, "       -keep    nnn - prune if excess > nnn (default 160)\n");
  fprintf (stderr, "       -minpa   nnn - purge alg. primes >= nnn (default alim)\n");
  fprintf (stderr, "       -minpr   nnn - purge rat. primes >= nnn (default rlim)\n");
  fprintf (stderr, "       -nprimes nnn - expected number of prime ideals\n");
  fprintf (stderr, "       -sos sosfile - to keep track of the renumbering\n");
  fprintf (stderr, "       -raw         - output relations in CADO format\n");
  exit (1);
}

/* Estimate the number of primes <= B. The 0.85 factor accounts for the
   fact that the combinatorial explosion happens before B/log(B). */
static int
approx_phi (long B)
{
  return (B <= 1) ? 0 : (int) (0.85 * (double) B / log ((double) B));
}

int
main (int argc, char **argv)
{
    hashtable_t H;
    int **rel_compact = NULL;
    int ret, k;
    int nrel, nprimes = 0;
    unsigned int nrelmax = 0;
    int nrel_new, nprimes_new, Hsize, Hsizer, Hsizea;
    long excess = 1;    /* minimum initial excess */
    long keep = 160;    /* maximum final excess */
    long minpr = -1, minpa = -1; /* -1 means use minpr=rlim, minpa=alim */
    cado_poly pol;
    unsigned long tot_alloc, tot_alloc0;
    int need64 = 0; /* non-zero if large primes are > 2^32 */
    int final, pass = 0;
    int raw = 0;
    char ** fic;

    fprintf (stderr, "%s.r%s", argv[0], CADO_REV);
    for (k = 1; k < argc; k++)
      fprintf (stderr, " %s", argv[k]);
    fprintf (stderr, "\n");

    param_list pl;
    param_list_init(pl);

    param_list_configure_knob(pl, "raw", &raw);

    argv++,argc--;

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Since we accept file names freeform, we decide to never abort
         * on unrecognized options */
        if (strcmp(*argv, "--help") == 0)
            usage();
        break;
        // fprintf (stderr, "Unknown option: %s\n", argv[0]);
        // abort(); 
    }

    param_list_parse_uint(pl, "nrels", &nrelmax);
    param_list_parse_int(pl, "nprimes", &nprimes);
    param_list_parse_long(pl, "minpr", &minpr);
    param_list_parse_long(pl, "minpa", &minpa);
    param_list_parse_long(pl, "excess", &excess);
    param_list_parse_long(pl, "keep", &keep);

    const char * filelist = param_list_lookup_string(pl, "filelist");
    const char * basepath = param_list_lookup_string(pl, "basepath");
    const char * subdirlist = param_list_lookup_string(pl, "subdirlist");
    const char * purgedname = param_list_lookup_string(pl, "out");
    const char * sos = param_list_lookup_string(pl, "sos");

    cado_poly_init (pol);

    const char * tmp;

    ASSERT_ALWAYS((tmp = param_list_lookup_string(pl, "poly")) != NULL);
    cado_poly_read(pol, tmp);

    if (param_list_warn_unused(pl)) {
        usage();
    }

    if ((basepath || subdirlist) && !filelist) {
        fprintf(stderr, "-basepath / -subdirlist only valid with -filelist\n");
        usage();
    }

    if (nrelmax == 0) {
        fprintf (stderr, "Error, missing -nrels ... option\n");
        usage();
    }

    /* On a 32-bit computer, even 1 << 32 would overflow. Well, we could set
       map[ra] = 2^32-1 in that case, but not sure we want to support 32-bit
       primes on a 32-bit computer... */
    need64 = (pol[0].lpbr >= 32) || (pol[0].lpba >= 32);

    if (need64 && sizeof (long) < 8)
      {
        fprintf (stderr, "Error, too large LPBs for a 32-bit computer\n");
        usage();
      }

    minus2 = (need64) ? 18446744073709551614UL : 4294967294UL;

    if (minpr == -1)
      minpr = pol[0].rlim;
    if (minpa == -1)
      minpa = pol[0].alim;

    fprintf (stderr, "Using minpr=%ld minpa=%ld\n", minpr, minpa);
    fprintf (stderr, "Number of relations is %u\n", nrelmax);
    if (nprimes > 0)
	Hsize = nprimes;
    else{
      /* Estimating the number of needed primes (remember that hashInit
         multiplies by a factor 1.5). */
      Hsizer = approx_phi (1L << pol[0].lpbr) - approx_phi (minpr);
      Hsizea = approx_phi (1L << pol[0].lpba) - approx_phi (minpa);
      Hsize = Hsizer + Hsizea;
    }
    fprintf (stderr, "initializing hash tables with Hsize=%d...\n", Hsize);
    hashInit (&H, Hsize, 1, need64);
    tot_alloc0 = H.hashmod * H.size;

    bit_vector rel_used;
    bit_vector_init_set(rel_used, nrelmax, 1);
    tot_alloc0 += nrelmax;
    fprintf (stderr, "Allocated rel_used of %uMb (total %luMb so far)\n",
             nrelmax >> 20,
             tot_alloc0 >> 20);

    /* FIXME: we could allocate rel_compact[] at each pass according to the
       number of *remaining* relations. We could also manage so that it takes
       only 32 bits per entry on a 64-bit machine. */
    rel_compact = (int **) malloc (nrelmax * sizeof (int *));
    tot_alloc0 += nrelmax * sizeof (int*);
    /* %zu is the C99 modifier for size_t */
    fprintf (stderr, "Allocated rel_compact of %zuMb (total %luMb so far)\n",
             (nrelmax * sizeof (int *)) >> 20,
             tot_alloc0 >> 20);

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
        fic = filelist_from_file(basepath, filelist);
    } else {
        /* count the number of files in the filelist */
        int nfiles = 0;
        int nsubdirs = 0;
        char ** fl = filelist_from_file(NULL, filelist);
        for(char ** p = fl ; *p ; p++, nfiles++);

        char ** sl = filelist_from_file(basepath, subdirlist);
        for(char ** p = sl ; *p ; p++, nsubdirs++);

        fic = malloc((nsubdirs * nfiles + 1) * sizeof(char*));
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
    while (1)
      {
        final = minpr <= 0 && minpa <= 0;
        tot_alloc = tot_alloc0;

        fprintf (stderr, "Pass %d:\n", ++pass);
        ret = scan_relations (fic, &nprimes, &H, rel_used, nrelmax,
                              rel_compact, minpr, minpa, &tot_alloc);
        ASSERT (ret);

        fprintf (stderr, "   nrels=%d, nprimes=%d; excess=%d\n",
                 nrel, nprimes, nrel - nprimes);

        fprintf (stderr, "   Starting singleton removal...\n");
        nrel_new = nrel;
        nprimes_new = nprimes;
        if (final && (nrel_new < nprimes_new + excess))
          {
            fprintf (stderr, "Initial excess is below requested %ld, stopping.\n",
                     excess);
            exit (1); /* initial excess is too small */
          }
        remove_singletons (&nrel_new, nrelmax, &nprimes_new, &H, rel_used,
                           rel_compact, keep, final);
        
        fprintf (stderr, "   nrel=%d, nprimes=%d; excess=%d\n",
                 nrel_new, nprimes_new, nrel_new - nprimes_new);

        if (nrel_new <= nprimes_new) /* covers case nrel = nprimes = 0 */
          exit (1); /* if the excess is <= here, it will be more negative
                       when we decrease minpr and minpa in the next pass */

        my_malloc_free_all ();
        if (final)
          break;
        else
          minpr = minpa = 0;
        hashClear (&H); /* Reset all tables to 0. FIXME: in the 2nd pass,
                           we might reuse some of the information computed
                           in the 1st pass (roots of primes, etc). */
        nrel = nrel_new;
      }

    fprintf (stderr, "Freeing rel_compact array...\n");
    /* we do not use it anymore */
    free (rel_compact);

    /* we renumber the primes in order of apparition in the hashtable */
    fprintf (stderr, "Renumbering primes...\n");
    renumber (&nprimes_new, &H, sos);

    /* reread the relation files and convert them to the new coding */
    fprintf (stderr, "Storing remaining relations...\n");
    reread (purgedname, fic, &H, rel_used, nrel_new, nprimes_new, raw);

    hashFree (&H);
    bit_vector_clear(rel_used);
    cado_poly_clear (pol);

    if (filelist) filelist_clear(fic);

    param_list_clear(pl);

    return 0;
}
