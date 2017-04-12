/* Characters
   
   Copyright 2009, 2010 Andreas Enge, Pierrick Gaudry, Fran\c{c}ois Morain, Emmanuel Thom\'e, Paul Zimmermann
   
   This file is part of CADO-NFS.
   
   CADO-NFS is free software; you can redistribute it and/or modify it
   under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation; either version 2.1 of the License, or
   (at your option) any later version.
   
   CADO-NFS is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more
   details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with CADO-NFS; see the file COPYING.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
   Boston, MA 02110-1301, USA.
*/

// Input:
//
// * [k] A bit matrix of (small_nrows) rows and (some multiple of 64)
//   cols, whose column span _contains_ the kernel of the matrix which
//   has been fed to the linear algebra.
// * A list of the (npurged) a,b pairs. This is obtained from the
//   purgedfile.
// * A matrix of (small_nrows) rows and (npurged) cols, which indicates
//   the contents of each relation-set. This is obtained from the
//   indexfile.
//
// Output:
//
// * A subspace of the column span of [k], whose vectors happen to also
//   cancel the characters used. By construction, it might contain zero
//   vectors, which are eventually discarded from the output.
//
// Algorithm.
//
// * [big character matrix, bcmat] First we build a matrix of size
//   npurged x nchars, where each row corresponds to the character values
//   at the corresponding pair (a,b)
// * [small character matrix, scmat] Then this matrix is multiplied by
//   the corresponding matrix taken from the indexfile.
//              scmat = index * bcmat
//   We thus obtain a matrix of size small_nrows x nchars, containing the
//   character values for each relation set.
// * [heavy block, h] According to the value of the ``skip'' parameter,
//   ``dense'' block which was taken out of the matrix during replay is
//   now read again. This is a matrix of size small_nrows x skip -- the
//   concatenation of scmat and h is thus a matrix of size (small_nrows)
//   x (nchars + skip), although this concatenation is not computed.
//
// * [t] At this point, if k were completely satisfactory, we would have:
//            transpose(k) * [scmat | h] == 0
//   this is a priori not the case. Therefore we compute the product:
//            t = transpose(k) * [scmat | h]
// * [kb] Now we compute the transpose of the left nullspace of t, namely
//   a matrix kb of size ncols(k) * ncols(k) (number of ``kernel
//   vectors'' out of bwc), and satisfying:
//            transpose(kb) * t == 0
//            transpose(k * kb) * [scmat | h] == 0
// * [nk] The ``new kernel'' k*kb is computed as nk.
// * [nkt] Its transpose is computed, so that each kernel vector can be
//   printed separately.
//
// Notes.
//
// Because [k] is only guaranteed to _contain_ the kernel in its column
// span, it is reasonable to first compute a basis of this column span.
// This is done on a heuristic basis, taking the first 4k coordinates
// only.

#include "cado.h"

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "mod_ul.c"

#include "cado-endian.h"
#include "portability.h"
#include "utils.h"
#include "utils_with_io.h"
#include "blockmatrix.h"
#include "gauss.h"
#include "worker-threads.h"


typedef struct {
    unsigned long p;      /* algebraic prime */
    unsigned long r;      /* corresponding root: r = a/b mod p */
    int e;                /* exponent (may want negative exponent in sqrt) */
} alg_prime_t;

/* Calculates a 64-bit word with the values of the characters chi(a,b), where
 * chi ranges from chars to chars+64
 */
uint64_t eval_64chars(int64_t a, int64_t b, alg_prime_t * chars, cado_poly_ptr pol)
{
    /* FIXME: do better. E.g. use 16-bit primes, and a look-up table. Could
     * beat this. */
    uint64_t v = 0;
    for(int i = 0 ; i < 64 ; i++) {
        alg_prime_t * ch = chars + i;
        int res;
        if (ch->p == 0) {
            // all special characters are identified by p==0
            if (ch->r == 0) {
                res = 0;        // trivial character
            } else if (ch->r == 1) {
                res = 1;        // parity character
            } else if (ch->r == 2) {
                /* Special: rational sign (sign of m1*a-m2*b) */
                mpz_t tmp1, tmp2;
                int ratside = cado_poly_get_ratside(pol);
		ASSERT_ALWAYS(ratside != -1);

                mpz_poly_ptr po = pol->pols[ratside];

                /* first perform a quick check: if m1*a and -m2*b have
                 * the same sign, we don't need to actually compute them.
                 */
                int sign_m1 = mpz_sgn(po->coeff[0]);
                int sign_m2 = mpz_sgn(po->coeff[1]);
                int m1a_negative = (a < 0) ? (sign_m1 < 0) : (sign_m1 > 0);
                int mm1b_negative = (b < 0) ? (sign_m2 > 0) : (sign_m2 < 0);
                if (m1a_negative == mm1b_negative) {
                    res = m1a_negative;
                } else {
                    mpz_init(tmp1);
                    mpz_mul_int64(tmp1, po->coeff[1], a);
                    mpz_init(tmp2);
                    mpz_mul_int64(tmp2, po->coeff[0], b);
                    mpz_sub(tmp1, tmp1, tmp2);
                    res = mpz_sgn(tmp1) < 0;
                    mpz_clear(tmp1);
                    mpz_clear(tmp2);
                }
            } else if (ch->r == 3) {
                res = (b==0);        // parity of the number of free relations
            } else {
                abort();
            }
        } else {
            /* Compute a+b*r (mod p) */
            residueul_t ra, rb, rr;
            modulusul_t mp;
            modul_initmod_ul(mp, ch->p);
            modul_init(ra, mp);
            modul_init(rb, mp);
            modul_init(rr, mp);
            if (a < 0) {
                modul_set_ul(ra, (unsigned long)(-a), mp);
                modul_neg(ra, ra, mp);
            } else {
                modul_set_ul(ra, (unsigned long)a, mp);
            }
            if (b < 0) {
                modul_set_ul(rb, (unsigned long)(-b), mp);
                modul_neg(rb, rb, mp);
            } else {
                modul_set_ul(rb, (unsigned long)b, mp);
            }
            modul_set_ul_reduced(rr, ch->r, mp);
            modul_mul(rr, rb, rr, mp);
            modul_add(rr, rr, ra, mp);
            res = modul_jacobi(rr, mp);
            modul_clear(ra, mp);
            modul_clear(rb, mp);
            modul_clear(rr, mp);
            modul_clearmod(mp);
            
            // If res is 0, it means that ch->p divides the norm, which
            // should not be, unless the special-q has been chosen to be
            // larger than lpb.
            // A fix, for that case that should not happen in real life
            // would be to artificially enlarge the large prime bound
            // if we had to sieve beyond it, so that subsequent program
            // (including characters) behaves correctly.
            if (res == 0) {
                fprintf (stderr, "Error, Jacobi symbol is 0 for a = %" PRId64 
                         ", b = %" PRIu64 ", p = %lu, r = %lu\n", 
                         a, b, ch->p, ch->r);
                ASSERT_ALWAYS(res != 0);
            }
            res = res < 0;   // -1->1, 1->0
        }
        v |= ((uint64_t) res) << i;
    }
    return v;
}

struct charbatch {
    uint64_t * W;
    int64_t * A;
    int64_t *B;
    unsigned int n;
    alg_prime_t * chars;
    cado_poly_ptr pol;
};
void eval_64chars_batch_thread(struct worker_threads_group * g, int tnum, void * t)
{
    struct charbatch * ss = (struct charbatch *) t;

    for(unsigned int z = tnum * ss->n / g->n ; z < (tnum + 1) * ss->n / g->n ; z++) {
        int64_t a = ss->A[z];
        int64_t b = ss->B[z];
        ss->W[z] = eval_64chars(a,b,ss->chars,ss->pol);
    }
    return;
}

static alg_prime_t * create_characters(int nchars[2],
        cado_poly pol, unsigned long *lpb)
{
    unsigned long p;
    int ret;
    mpz_t pp;
    unsigned long *roots;

    ASSERT_ALWAYS(nchars[0] + nchars[1] > 0);

    int nchars2 = iceildiv(nchars[0] + nchars[1], 64) * 64;

    mpz_init (pp);
    roots = malloc(MAX(pol->pols[0]->deg, pol->pols[1]->deg)
            * sizeof(unsigned long));

    alg_prime_t * chars = malloc(nchars2 * sizeof(alg_prime_t));

    int nspecchar = 2;
    /* force parity */
    chars[0] = (alg_prime_t) { .p = 0, .r = 1 };
    /* force parity of the free relations. This is really only because
     * we're lazy -- it's been asserted that it eases stuff at some
     * point, but nobody remembers the why and how. */
    chars[1] = (alg_prime_t) { .p = 0, .r = 3 };
    if (nchars[0] == 0 || nchars[1] == 0) {
        ASSERT_ALWAYS(MIN(pol->pols[0]->deg, pol->pols[1]->deg) == 1);
        /* force rational sign */
        chars[2] = (alg_prime_t) { .p = 0, .r = 2 };
        nspecchar++;
    }

    /* we might want to force evenness of the number of relations as well. Easy
     * to put this in chars[1] if needed (and add the appropriate stuff above
     * of course).  */

    /* Rational characters. Normally we have none. But the -nratchars
     * option inserts some */
    /* we want some prime beyond the (rational) large prime bound */

    int i = nspecchar;
    for (int side = 0; side < 2; ++side) {
        if (nchars[side] == 0)
            continue;
        mpz_set_ui (pp, 1UL << lpb[side]);

        int j = 0;
        do {
            mpz_nextprime(pp, pp);
            p = mpz_get_ui(pp);
            ret = mpz_poly_roots_ulong (roots, pol->pols[side], p);
            for (int k = 0; k < ret; ++k) {
                if (i == nchars[0] + nchars[1])
                    break;
                if (j == nchars[side])
                    break;
                chars[i].p = p;
                chars[i].r = roots[k];
                ++i;
                ++j;
            }
        } while(j < nchars[side] && (i < nchars[0] + nchars[1]));
    }

    /* pad with trivial characters */
    for(int i = nchars[0] + nchars[1] ; i < nchars2 ; i++) {
        chars[i] = (alg_prime_t) { .p = 0, .r = 0 };
    }
    if (nchars[0] + nchars[1] < nchars2) {
        fprintf(stderr, "Note: total %d characters, "
                "including %d trivial padding characters\n",
                nchars2, nchars2-(nchars[0] + nchars[1]));
    }

    free(roots);
    mpz_clear(pp);

    return chars;
}

typedef struct
{
  int64_t *a;
  int64_t *b;
} chars_data_t;

void *
thread_chars (void * context_data, earlyparsed_relation_ptr rel)
{
  chars_data_t *data = (chars_data_t *) context_data;
  int nab = EARLYPARSED_RELATION_MAX_AB;
  for( ; nab > 2 && !rel->ab[nab-1] ; nab--);
  ASSERT_ALWAYS(nab == 2 && rel->ab[0] >= 0);

  data->a[rel->num] = rel->ab[0];
  data->b[rel->num] = rel->ab[1];
  return NULL;
}

// The big character matrix has (number of purged rels) rows, and (number of
// characters) cols

static blockmatrix big_character_matrix(alg_prime_t * chars, unsigned int nchars2, const char * purgedname, cado_poly_ptr pol, struct worker_threads_group * g)
{
    uint64_t nrows, ncols;
    purgedfile_read_firstline (purgedname, &nrows, &ncols);

    int64_t *all_A = (int64_t *) malloc (nrows * sizeof(int64_t));
    int64_t *all_B = (int64_t *) malloc (nrows * sizeof(int64_t));
    ASSERT_ALWAYS(all_A != NULL && all_B != NULL);
    blockmatrix res = blockmatrix_alloc(nrows, nchars2);
    blockmatrix_set_zero(res);

    /* For each rel, read the a,b-pair and init the corresponding poly pairs[] */
    fprintf(stderr, "Reading %" PRIu64 " (a,b) pairs from %s\n", nrows,
                     purgedname);
    chars_data_t data = {.a = all_A, .b=all_B};
    char *fic[2] = {(char *) purgedname, NULL};
    filter_rels (fic, (filter_rels_callback_t) thread_chars, &data,
          EARLYPARSE_NEED_AB, NULL, NULL);

    fprintf(stderr, "Computing %u characters for %" PRIu64 " (a,b) pairs\n",
            nchars2, nrows);

    for(uint64_t i = 0 ; i < nrows; ) {
        static const int batchsize = 16384;
        int64_t A[batchsize];
        int64_t B[batchsize];
        uint64_t W[batchsize];
        int bs = 0;
        while (bs < batchsize && i+bs < nrows) {
            A[bs] = all_A[i+bs];
            B[bs] = all_B[i+bs];
            bs++;
        }
        struct charbatch ss = { .W=W,.A=A,.B=B,.pol=pol,.n=bs };
        for(unsigned int cg = 0 ; cg < nchars2 ; cg+=64) {
            ss.chars = chars + cg;
            worker_threads_do(g, eval_64chars_batch_thread, &ss);
            for(int z = 0 ; z < bs ; z++) {
                *blockmatrix_subrow_ptr(res, i+z, cg) = W[z];
            }
        }
        i += bs;
    }
    free (all_A);
    free (all_B);

    return res;
}

/* The small character matrix has only (number of relation-sets) rows -- its
 * number of cols is still (number of characters) */
static blockmatrix small_character_matrix(blockmatrix bcmat, const char * indexname)
{
    FILE * ix = fopen_maybe_compressed(indexname, "r");
    uint64_t small_nrows;
    int ret;

    ret = fscanf(ix, "%" SCNu64 "\n", &small_nrows);
    ASSERT(ret == 1);

    unsigned int nchars2 = bcmat->ncols;

    blockmatrix res = blockmatrix_alloc(small_nrows, nchars2);

    for(uint64_t i = 0 ; i < small_nrows ; i++) {
        int nc;
        ret = fscanf(ix, "%d", &nc); ASSERT_ALWAYS(ret == 1);
        for(unsigned int cg = 0 ; cg < nchars2 ; cg+=64) {
            res->mb[(i/64) + (cg/64) * res->nrblocks][i%64] = 0;
        }
        for(int k = 0 ; k < nc ; k++) {
            unsigned int col;
            ret = fscanf(ix, "%x", &col); ASSERT_ALWAYS(ret == 1);
            ASSERT_ALWAYS(col < bcmat->nrows);
            for(unsigned int cg = 0 ; cg < nchars2 ; cg+=64) {
                res->mb[(i/64) + (cg/64) * res->nrblocks][i%64] ^=
                    bcmat->mb[(col/64) + (cg/64) * bcmat->nrblocks][col%64];
            }
        }
    }
    fclose_maybe_compressed(ix, indexname);
    return res;
}

/* We support both ascii and binary format, which is close to a bug. */

static blockmatrix
read_heavyblock_matrix_binary(const char * heavyblockname)
{
    FILE * f = fopen(heavyblockname, "rb");

    if (f == NULL) {
        fprintf(stderr, "Warning: %s not found, assuming empty\n", heavyblockname);
        return blockmatrix_alloc(0,0);
    }

    unsigned int nrows, ncols;

    /* If we're binary, we insist on having the companion files as well,
     * which provide a quick hint at the matrix dimensions */

    {
        char * rwname = derived_filename(heavyblockname, "rw", ".bin");
        struct stat sbuf[1];
        int rc = stat(rwname, sbuf);
        if (rc < 0) { perror(rwname); exit(1); }
        nrows = sbuf->st_size / sizeof(uint32_t);
        free(rwname);
    }
    {
        char * cwname = derived_filename(heavyblockname, "cw", ".bin");
        struct stat sbuf[1];
        int rc = stat(cwname, sbuf);
        if (rc < 0) { perror(cwname); exit(1); }
        ncols = sbuf->st_size / sizeof(uint32_t);
        free(cwname);
    }

    blockmatrix res = blockmatrix_alloc(nrows, ncols);

    /* Sometimes the heavy block width is not a multiple of 64. Thus we pad
     * with zeros */
    blockmatrix_set_zero(res);

    for(unsigned int i = 0 ; i < nrows ; i++) {
        uint32_t len;
        int r = fread32_little(&len, 1, f);
        ASSERT_ALWAYS(r == 1);
        for( ; len-- ; ) {
            uint32_t v;
            r = fread32_little(&v, 1, f); ASSERT_ALWAYS(r == 1);
            res->mb[(i/64) + (v/64) * res->nrblocks][i%64] ^= ((uint64_t)1) << (v%64);
        }
    }
    fclose (f);
    return res;
}

static blockmatrix
read_heavyblock_matrix_ascii(const char * heavyblockname)
{
    FILE * f = fopen(heavyblockname, "r");

    if (f == NULL) {
        fprintf(stderr, "Warning: %s not found, assuming empty\n", heavyblockname);
        return NULL;
    }

    unsigned int nrows, ncols;
    int rc = fscanf(f,"%u %u", &nrows, &ncols);
    ASSERT_ALWAYS(rc == 2);

    blockmatrix res = blockmatrix_alloc(nrows, ncols);

    /* Sometimes the heavy block width is not a multiple of 64. Thus we pad
     * with zeros */
    blockmatrix_set_zero(res);

    for(unsigned int i = 0 ; i < nrows ; i++) {
        uint32_t len;
        int r = fscanf(f, "%" SCNu32, &len); ASSERT_ALWAYS(r == 1);
        for( ; len-- ; ) {
            uint32_t v;
            r = fscanf(f, "%" SCNu32, &v); ASSERT_ALWAYS(r == 1);
            res->mb[(i/64) + (v/64) * res->nrblocks][i%64] ^= ((uint64_t)1) << (v%64);
        }
    }
    fclose (f);
    return res;
}

static blockmatrix
read_heavyblock_matrix (const char * heavyblockname)
{
  if (heavyblockname == NULL)
    return blockmatrix_alloc (0, 0);
  else if (has_suffix(heavyblockname, ".bin"))
    return read_heavyblock_matrix_binary(heavyblockname);
  else
    return read_heavyblock_matrix_ascii(heavyblockname);
}

int compute_transpose_of_blockmatrix_kernel(blockmatrix kb, blockmatrix t)
{
    /* gauss.c's kernel() function takes its input with a different ordering.
     * It's tiny data anyway. */

    fprintf(stderr, "Computing left nullspace of %u x %u matrix\n",
            t->nrows, t->ncols);
    unsigned int tiny_nrows = t->nrows;
    unsigned int tiny_ncols = t->ncols;
    unsigned int tiny_limbs_per_row = iceildiv(tiny_ncols, 64);
    unsigned int tiny_limbs_per_col = iceildiv(tiny_nrows, 64);
    unsigned int tiny_chars = FLAT_BYTES_WITH_READAHEAD(t->nrows, t->ncols);
    unsigned int tiny_64bit_words = tiny_chars / sizeof(uint64_t);

    /* we need some readahead zones because of the block matrix structure */
    uint64_t * tiny = malloc (tiny_chars);
    memset(tiny, 0, tiny_chars);
    
    blockmatrix_copy_to_flat(tiny, tiny_limbs_per_row, 0, 0, t);

    /* The kernel matrix is essentially a square matrix of tiny_nrows rows and
     * columns (tiny_nrows is the same as the number of kernel vectors)
     */
    /* we need some readahead zones because of the block matrix structure */
    uint64_t * kerdata = malloc(FLAT_BYTES_WITH_READAHEAD(t->nrows, t->nrows));
    memset(kerdata, 0, FLAT_BYTES_WITH_READAHEAD(t->nrows, t->nrows));
    unsigned int kerdata_64bit_words = FLAT_BYTES_WITH_READAHEAD(t->nrows, t->nrows) / sizeof(uint64_t);


    uint64_t ** myker = (uint64_t **) malloc(tiny_nrows * sizeof(uint64_t *));
    ASSERT(myker != NULL);
    for (unsigned int i = 0; i < tiny_nrows; ++i)
        myker[i] = kerdata + i * tiny_limbs_per_col;
    /* gauss.c knows about mp_limb_t's only */
    ASSERT_ALWAYS(sizeof(uint64_t) % sizeof(mp_limb_t) == 0);
    swap_words_if_needed (tiny, tiny_64bit_words);
    int dim = kernel((mp_limb_t *) tiny,
            (mp_limb_t **) myker,
            tiny_nrows, tiny_ncols,
            sizeof(uint64_t) / sizeof(mp_limb_t) * tiny_limbs_per_row,
            sizeof(uint64_t) / sizeof(mp_limb_t) * tiny_limbs_per_col);
    swap_words_if_needed (tiny, tiny_64bit_words); /* FIXME: this is maybe not
                                                      needed since tiny is
                                                      destroyed, but keep it
                                                      for debugging */
    swap_words_if_needed(kerdata, kerdata_64bit_words);
    free(tiny);
    /* Now take back our kernel to block format, and multiply. Exciting. */
    if (kb)
      blockmatrix_copy_transpose_from_flat(kb, kerdata, tiny_limbs_per_col,
                                           0, 0);
    free(myker);
    free(kerdata);
    return dim;
}

/* This only builds a basis, not an echelonized basis */
blockmatrix blockmatrix_column_reduce(blockmatrix m, unsigned int max_rows_to_consider)
{
    blockmatrix t = blockmatrix_submatrix(m, 0, 0, MIN(max_rows_to_consider, m->nrows), m->ncols);

    unsigned int tiny_nrows = t->ncols;
    unsigned int tiny_ncols = t->nrows;
    unsigned int tiny_limbs_per_row = iceildiv(tiny_ncols, 64);
    unsigned int tiny_limbs_per_col = iceildiv(tiny_nrows, 64);

    unsigned int tiny_nlimbs = tiny_nrows * tiny_limbs_per_row;
    uint64_t * tiny = malloc(tiny_nlimbs * sizeof(uint64_t));
    memset(tiny, 0, tiny_nlimbs * sizeof(uint64_t));
    blockmatrix_copy_transpose_to_flat(tiny, tiny_limbs_per_row, 0, 0, t);

    uint64_t * sdata = (uint64_t *) malloc(tiny_nrows * tiny_limbs_per_col * sizeof(uint64_t));
    memset(sdata, 0, tiny_nrows * tiny_limbs_per_col * sizeof(uint64_t *));
    unsigned int sdata_64bit_words = tiny_nrows * tiny_limbs_per_col;
    free(t);

    swap_words_if_needed (tiny, tiny_nlimbs);
    int rank = spanned_basis(
            (mp_limb_t *) sdata,
            (mp_limb_t *) tiny,
            tiny_nrows,
            tiny_ncols,
            sizeof(uint64_t) / sizeof(mp_limb_t) * tiny_limbs_per_row,
            sizeof(uint64_t) / sizeof(mp_limb_t) * tiny_limbs_per_col,
            NULL
            );
    swap_words_if_needed (tiny, tiny_nlimbs);
    swap_words_if_needed (sdata, sdata_64bit_words);
    free(tiny);

    blockmatrix s = blockmatrix_alloc(m->ncols, rank);
    blockmatrix_copy_transpose_from_flat(s, sdata, tiny_limbs_per_col, 0, 0);
    blockmatrix k2 = blockmatrix_alloc(m->nrows, rank);
    blockmatrix_mul_smallb(k2, m, s);
    blockmatrix_free(s);
    free(sdata);
    return k2;
}

static void
declare_usage (param_list pl)
{
  param_list_decl_usage (pl, "purged", "output-from-purge file");
  param_list_decl_usage (pl, "index",  "index file");
  param_list_decl_usage (pl, "out",    "output file");
  param_list_decl_usage (pl, "heavyblock", "heavyblock output file");
  param_list_decl_usage (pl, "poly",   "polynomial file");
  param_list_decl_usage (pl, "nchar",  "number of characters");
  param_list_decl_usage (pl, "lpb0",   "large prime bound on side 0");
  param_list_decl_usage (pl, "lpb1",   "large prime bound on side 1");
  param_list_decl_usage (pl, "t",      "number of threads");
  param_list_decl_usage (pl, "ker",      "input kernel file");
  param_list_decl_usage(pl, "force-posix-threads", "(switch)");
}

int main(int argc, char **argv)
{
    const char * heavyblockname = NULL;
    int nchars;
    alg_prime_t *chars;
    cado_poly pol;
    const char *purgedname = NULL;
    const char *indexname = NULL;
    const char *outname = NULL;
    int nthreads = 1;
    unsigned long lpb[2] = {0,0};
    char *argv0 = argv[0];

    /* print the command line */
    fprintf (stderr, "%s.r%s", argv[0], cado_revision_string);
    for (int i = 1; i < argc; i++)
      fprintf (stderr, " %s", argv[i]);
    fprintf (stderr, "\n");

    param_list pl;
    param_list_init(pl);
    declare_usage(pl);

    argc--,argv++;
    const char *bw_kernel_file = NULL;

    param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;

        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        param_list_print_usage(pl, argv0, stderr);
        exit (EXIT_FAILURE);
    }
    purgedname = param_list_lookup_string(pl, "purged");
    indexname = param_list_lookup_string(pl, "index");
    outname = param_list_lookup_string(pl, "out");
    heavyblockname = param_list_lookup_string(pl, "heavyblock");
    bw_kernel_file = param_list_lookup_string(pl, "ker");

    cado_poly_init (pol);

    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "poly")) == NULL)
      {
        fprintf (stderr, "Error: parameter -poly is mandatory\n");
        param_list_print_usage (pl, argv0, stderr);
        exit (EXIT_FAILURE);
      }
    cado_poly_read(pol, tmp);

    if (param_list_parse_int(pl, "nchar", &nchars) == 0)
      {
        fprintf (stderr, "Error: parameter -nchar is mandatory\n");
        param_list_print_usage (pl, argv0, stderr);
        exit (EXIT_FAILURE);
      }
    if (param_list_parse_ulong(pl, "lpb0", &lpb[0]) == 0)
      {
        fprintf (stderr, "Error: parameter -lpb0 is mandatory\n");
        param_list_print_usage (pl, argv0, stderr);
        exit (EXIT_FAILURE);
      }
    if (param_list_parse_ulong(pl, "lpb1", &lpb[1]) == 0)
      {
        fprintf (stderr, "Error: parameter -lpb1 is mandatory\n");
        param_list_print_usage (pl, argv0, stderr);
        exit (EXIT_FAILURE);
      }

    param_list_parse_int(pl, "t", &nthreads);

    if (purgedname == NULL || indexname == NULL || outname == NULL) {
        fprintf (stderr,
                "Error: parameters -purged, -index and -out are mandatory\n");
        param_list_print_usage (pl, argv0, stderr);
        exit (EXIT_FAILURE);
    }

    /* Put characters on all algebraic sides */
    int nch[2] = {0, 0};
    for (int side = 0; side < 2; ++side) {
        if (pol->pols[side]->deg > 1)
            nch[side] = nchars;
    }

    struct worker_threads_group * g = worker_threads_init (nthreads);
    chars = create_characters (nch, pol, lpb);
    int nchars2 = iceildiv(nch[0] + nch[1], 64) * 64;
    double tt=wct_seconds();
    blockmatrix bcmat = big_character_matrix(chars, nchars2, purgedname, pol, g);
    free(chars);
    worker_threads_clear(g);

    fprintf(stderr, "done building big character matrix at %.1f\n", wct_seconds()-tt);

    blockmatrix scmat = small_character_matrix(bcmat, indexname);
    fprintf(stderr, "done building small character matrix at %.1f\n", wct_seconds()-tt);

    blockmatrix_free(bcmat);
    
    unsigned int small_nrows = scmat->nrows;

    /* It's ok if heavyblockname == 0. After all sufficiently many
     * characters should be enough */
    blockmatrix h = read_heavyblock_matrix(heavyblockname);
    if (h->ncols == 0) { h->nrows = small_nrows; }
    if (h->ncols)
          fprintf(stderr, "done reading heavy block of size %u x %u at %.1f\n",
                  h->nrows, h->ncols, wct_seconds()-tt);
    ASSERT_ALWAYS(h->nrows == small_nrows);

    /* Now do dot products of these matrices by the kernel vectors
     * supplied on input */

    /* First compute how many kernel vectors we have */

    unsigned int total_kernel_cols = 0;
    unsigned int kncols;
    struct stat sbuf[1];
    int rc = stat(bw_kernel_file, sbuf);
    if (rc < 0) { perror(bw_kernel_file); exit(1); }
    ASSERT_ALWAYS(sbuf->st_size % small_nrows == 0);
    unsigned int ncols = 8 * (sbuf->st_size / small_nrows);
    fprintf(stderr, "%s: %u vectors\n", bw_kernel_file, ncols);
    ASSERT_ALWAYS(ncols % 64 == 0);
    total_kernel_cols += kncols = ncols;

    fprintf(stderr, "Total: %u kernel vectors\n", total_kernel_cols);

    /* kmat is the join of all kernel vectors */
    blockmatrix k = blockmatrix_alloc(small_nrows, total_kernel_cols);
    int j0 = 0;
    blockmatrix_read_from_flat_file(k, 0, j0, bw_kernel_file, small_nrows, kncols);

    fprintf(stderr, "done reading %u kernel vectors at %.1f\n",
            total_kernel_cols, wct_seconds() - tt);

    blockmatrix k2 = blockmatrix_column_reduce(k, 4096);
    blockmatrix_free(k);
    k = k2;
    total_kernel_cols = k->ncols;
    fprintf(stderr, "Info: input kernel vectors reduced to dimension %u\n",
            total_kernel_cols);

    /* tmat is the product of the character matrices times the kernel vector */
    blockmatrix t = blockmatrix_alloc(total_kernel_cols, scmat->ncols + h->ncols);
    blockmatrix tc = blockmatrix_submatrix(t, 0, 0, total_kernel_cols, scmat->ncols);
    blockmatrix th = blockmatrix_submatrix(t, 0, scmat->ncols, total_kernel_cols, h->ncols);

    blockmatrix_mul_Ta_b(tc, k, scmat);
    blockmatrix_mul_Ta_b(th, k, h);

    fprintf(stderr, "done multiplying matrices at %.1f\n", wct_seconds() - tt);

    free(tc);
    free(th);
    blockmatrix_free(scmat);
    blockmatrix_free(h);


    blockmatrix kb = blockmatrix_alloc(k->ncols, k->ncols);
    blockmatrix_set_zero(kb);
    int dim = compute_transpose_of_blockmatrix_kernel(kb, t);
    blockmatrix_free(t);
    fprintf(stderr, "dim of ker = %d\n", dim);

    blockmatrix kbsub = blockmatrix_submatrix(kb, 0, 0, kb->nrows, dim);

    blockmatrix nk = blockmatrix_alloc(small_nrows, kbsub->ncols);
    blockmatrix_mul_smallb(nk, k, kbsub);
    free(kbsub);
    blockmatrix_free(k);
    blockmatrix_free(kb);

    /* Sanity check: count the number of zero dependencies */
    unsigned int nonzero_deps = 0;
    for(unsigned int j = 0 ; j < nk->ncblocks ; j++) {
        uint64_t sanity = 0 ;
        for(unsigned int i = 0 ; i < nk->nrows ; i+=64) {
            for(unsigned int ii = 0 ; ii < 64 && i + ii < nk->nrows ; ii++)
                sanity |= nk->mb[j*nk->stride + ii/64][i];
        }
        // do popcount.
        for( ; sanity ; sanity>>=1) nonzero_deps += sanity&1UL;
    }
    if (!nonzero_deps) {
        fprintf(stderr, "Error, all dependencies are zero !\n");
        exit(1);
    }
    blockmatrix_write_to_flat_file(outname, nk, 0, 0, nk->nrows, nk->ncols);

    fprintf(stderr, "Wrote %d non-zero dependencies to %s\n",
            nonzero_deps, outname);
    if (nonzero_deps < (unsigned int) dim || nk->ncols % 64) {
        fprintf(stderr, "This includes %u discarded zero dependencies, as well as %u padding zeros\n",
                dim - nonzero_deps,
                nk->ncblocks * 64 - dim);
    }
    blockmatrix_free(nk);

    cado_poly_clear(pol);
    param_list_clear(pl);

    return 0;
}
