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
//   cols, whose column span _contains_ a subspace of the left kernel of
//   the matrix which has been fed to the linear algebra.
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

#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <cinttypes>
#include <cstdint>

#include <algorithm>
#include <numeric>
#include <utility>
#include <vector>

#include <sys/stat.h>

#include <gmp.h>
#include "fmt/base.h"

#include "gmp_aux.h"
#include "bblas_gauss.h"
#include "blockmatrix.hpp"
#include "cado_poly.hpp"
#include "filter_io.hpp"
#include "fix-endianness.h"
#include "gzip.h"
#include "macros.h"
#include "misc.h"
#include "arith/mod_ul.h"
#include "mpz_poly.h"
#include "params.hpp"
#include "purgedfile.h"
#include "rootfinder.h"
#include "submatrix_range.hpp"
#include "timing.h"
#include "version_info.h"

typedef struct {
    unsigned long p;      /* algebraic prime */
    unsigned long r;      /* corresponding root: r = a/b mod p */
    // int e;                /* exponent (may want negative exponent in sqrt) */
} alg_prime_t;

/* Computes chi(a,b), returns -1, 0 or 1 */
static int
eval_one_non_special_char(alg_prime_t const & chi, int64_t a, uint64_t b)
{
    /* Compute b*r-a (mod p) */
    residueul_t ra, rb, rr;
    modulusul_t mp;
    modul_initmod_ul(mp, chi.p);
    modul_init(ra, mp);
    modul_init(rb, mp);
    modul_init(rr, mp);
    if (a < 0) {
        modul_set_ul(ra, (unsigned long)(-a), mp);
        modul_neg(ra, ra, mp);
    } else {
        modul_set_ul(ra, a, mp);
    }
    // coverity[result_independent_of_operands]
    ASSERT(b <= ULONG_MAX);
    modul_set_ul(rb, (unsigned long)b, mp);
    modul_set_ul_reduced(rr, chi.r, mp);
    modul_mul(rr, rb, rr, mp);
    modul_sub(rr, rr, ra, mp);
    int res = modul_jacobi(rr, mp);
    modul_clear(ra, mp);
    modul_clear(rb, mp);
    modul_clear(rr, mp);
    modul_clearmod(mp);
    return res;
}

static int
eval_one_non_special_char(alg_prime_t const & chi, mpz_srcptr a, mpz_srcptr b)
{
    int64_t sa = (int64_t) mpz_tdiv_ui(a, chi.p);
    if (mpz_sgn(a) < 0) {
        sa = -sa;
    }
    return eval_one_non_special_char(chi, sa, mpz_tdiv_ui(b, chi.p));
}

static inline cxx_mpz eval_helper(cxx_mpz_poly const & f, int64_t a, uint64_t b)
{
    cxx_mpz t;
    mpz_poly_homogeneous_eval_siui(t, f, a, b);
    return t;
}

static inline cxx_mpz eval_helper(cxx_mpz_poly const & f, cxx_mpz const & a, cxx_mpz const & b)
{
    cxx_mpz t;
    mpz_poly_homogeneous_eval(t, f, a, b);
    return t;
}

/* Calculates a 64-bit word with the values of the characters chi(a,b), where
 * chi ranges from chars to chars+64
 */
template<typename relation_type>
static uint64_t eval_64chars(
        relation_type const & rel,
        alg_prime_t const * chars,
        cxx_cado_poly const & cpoly)
{
    auto const & a = rel.a;
    auto const & b = rel.b;
    /* FIXME: do better. E.g. use 16-bit primes, and a look-up table. Could
     * beat this. */
    uint64_t v = 0;
    for(int i = 0 ; i < 64 ; i++) {
        alg_prime_t const & chi = *(chars + i);
        int res;
        if (chi.p == 0) {
            // all special characters are identified by p==0
            if (chi.r == 0) {
                res = 0;        // trivial character
            } else if (chi.r == 1) {
                res = 1;        // parity character
            } else if (chi.r == 2 || chi.r >> 10 == 1) {
                /* Special: sign */
                int side;
                if (chi.r == 2) {
                    side = cpoly.get_ratside();
                    ASSERT_ALWAYS(side != -1);
                } else {
                    side = chi.r - 1024;
                    ASSERT_ALWAYS(side < cpoly.nsides());
                }
                mpz_poly_srcptr po = cpoly[side];
                int const deg = po->deg;
                ASSERT_ALWAYS(deg > 0);
                bool a_is_pos = a > 0;
                /* first perform a quick check for deg=1 */
                res = a_is_pos ? mpz_sgn(mpz_poly_coeff_const(po, 1))
                               : -mpz_sgn(mpz_poly_coeff_const(po, 1));
                if (deg > 1 || mpz_sgn(mpz_poly_coeff_const(po, 0)) != res) {
                    cxx_mpz t = eval_helper(po, a, b);
                    res = mpz_sgn(t) < 0;
                } else {
                    res = res < 0;
                }
            } else if (chi.r == 3) {
                // parity of the number of free relations
                res = (b==0);
            } else {
                fmt::print(stderr, "Error, unrecognized special character: "
                                   "r={}\n", chi.r);
                abort();
            }
        } else {
            res = eval_one_non_special_char(chi, a, b);

            // If res is 0, it means that chi.p divides the norm, which
            // should not be, unless the special-q has been chosen to be
            // larger than lpb.
            // A fix, for that case that should not happen in real life
            // would be to artificially enlarge the large prime bound
            // if we had to sieve beyond it, so that subsequent program
            // (including characters) behaves correctly.
            if (res == 0) {
                fmt::print(stderr, "Error, Jacobi symbol is 0 for a = {}, "
                                   "b = {} , p = {}, r = {}\n",
                                   cxx_mpz(a), cxx_mpz(b), chi.p, chi.r);
                fprintf (stderr, "Please check lpb0/lpb1 are large enough\n");
                ASSERT_ALWAYS(res != 0);
            }
            res = res < 0;   // -1->1, 1->0
        }
        v |= ((uint64_t) res) << i;
    }
    return v;
}

static std::vector<alg_prime_t>
create_characters(
        std::vector<int> const & nchars,
        cxx_cado_poly const & cpoly,
        const unsigned long *lpb)
{
    unsigned long p;
    int ret;
    cxx_mpz pp;
    std::vector<unsigned long> roots;

    const int nchars_tot = std::accumulate(nchars.begin(), nchars.end(), 0);
    ASSERT_ALWAYS(nchars_tot > 0);

    int const nchars2 = iceildiv(nchars_tot, 64) * 64;

    for (int side = 0; side < cpoly.nsides(); ++side) {
        if (cpoly[side]->deg > (int) roots.size()) {
            roots.resize(cpoly[side]->deg);
        }
    }
    std::vector<alg_prime_t> chars(nchars2);

    int nspecchar = 2;
    /* force parity */
    chars[0] = (alg_prime_t) { .p = 0, .r = 1 };
    /* force parity of the free relations. This is really only because
     * we're lazy -- it's been asserted that it eases stuff at some
     * point, but nobody remembers the why and how. */
    chars[1] = (alg_prime_t) { .p = 0, .r = 3 };
    if (std::ranges::any_of(nchars, [](int nc) { return nc == 0; })) {
        ASSERT_ALWAYS(cpoly.get_ratside() != -1);
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

    cxx_gmp_randstate rstate;

    int i = nspecchar;
    for (int side = 0; side < cpoly.nsides(); ++side) {
        if (nchars[side] == 0)
            continue;
        mpz_set_ui (pp, 1UL << lpb[side]);

        int j = 0;
        do {
            mpz_nextprime(pp, pp);
            p = mpz_get_ui(pp);
            ret = mpz_poly_roots_ulong (roots.data(), cpoly[side], p,
                                        rstate);
            for (int k = 0; k < ret && j < nchars[side] && i < nchars_tot;
                                                                ++k, ++i, ++j) {
                chars[i].p = p;
                chars[i].r = roots[k];
            }
        } while(j < nchars[side] && i < nchars_tot);
    }

    /* pad with trivial characters */
    for(int i = nchars_tot; i < nchars2 ; i++) {
        chars[i] = (alg_prime_t) { .p = 0, .r = 0 };
    }
    if (nchars_tot < nchars2) {
        fmt::print(stderr, "Note: total {} characters, including {} trivial "
                           "padding characters\n", nchars2, nchars2-nchars_tot);
    }

    return chars;
}

static std::vector<alg_prime_t>
create_sign_characters_only(cxx_cado_poly const & cpoly)
{
    std::vector<alg_prime_t> chars(iceildiv(cpoly.nsides(), 64) * 64);

    for(unsigned long i = 0; i < chars.size() ; i++) {
        if ((int) i < cpoly.nsides()) { /* sign character for side i */
            chars[i] = (alg_prime_t) { .p = 0, .r = 1024 + i };
        } else { /* trivial character */
            chars[i] = (alg_prime_t) { .p = 0, .r = 0 };
        }
    }
    if ((size_t) cpoly.nsides() < chars.size()) {
        fmt::print(stderr, "Note: total {} characters, including {} trivial "
                           "padding characters\n", chars.size(),
                           chars.size()-cpoly.nsides());
    }

    return chars;
}

struct big_characters_data {
    cxx_cado_poly const & cpoly;
    std::string purgedname;
    std::vector<alg_prime_t> const & chars;
    blockmatrix res;

    static uint64_t purgedfile_nrows(const char * purgedname)
    {
        uint64_t nrows, ncols;
        purgedfile_read_firstline (purgedname, &nrows, &ncols);
        return nrows;
    }

    big_characters_data(cxx_cado_poly const & cpoly,
            const char * purgedname,
            std::vector<alg_prime_t> const & chars)
        : cpoly(cpoly)
        , purgedname(purgedname)
        , chars(chars)
        , res(purgedfile_nrows(purgedname), chars.size())
    {
        res.set_zero();
    }

    template<typename relation_type>
    void filter(int nthreads) {
        /* For each rel, read the a,b-pair, compute the characters and
         * store them in res.
         */
        fmt::print(stderr, "Reading {} pairs from {} and computing {} characters\n",
                           res.nrows(), purgedname, chars.size());

        filter_rels<relation_type>(purgedname, nullptr, nullptr,
            cado::filter_io_details::multithreaded_call(nthreads,
            [&](relation_type & rel) {
                for(size_t cg = 0 ; cg < chars.size() ; cg+=64) {
                    res[rel.num][cg] = eval_64chars(rel, chars.data() + cg, cpoly);
                }
            }));
    }
};

// The big character matrix has (number of purged rels) rows, and (number of
// characters) cols
static blockmatrix big_character_matrix(
        std::vector<alg_prime_t> const & chars,
        const char * purgedname,
        cxx_cado_poly const & cpoly,
        int nthreads,
        int largeab)
{
    big_characters_data B(cpoly, purgedname, chars);
    if (!largeab) {
        using relation_type = cado::relation_building_blocks::ab_block<uint64_t, 16>;
        B.filter<relation_type>(nthreads);
    } else {
        using relation_type = cado::relation_building_blocks::ab_block<cxx_mpz, 16>;
        B.filter<relation_type>(nthreads);
    }

    return std::move(B.res);
}

/* The small character matrix has only (number of relation-sets) rows -- its
 * number of cols is still (number of characters) */
static blockmatrix small_character_matrix(blockmatrix & bcmat, const char * indexname)
{
    FILE * ix = fopen_maybe_compressed(indexname, "r");
    uint64_t small_nrows;
    int ret;

    ret = fscanf(ix, "%" SCNu64 "\n", &small_nrows);
    ASSERT(ret == 1);

    unsigned int const nchars2 = bcmat.ncols();

    blockmatrix res(small_nrows, nchars2);

    for(uint64_t i = 0 ; i < small_nrows ; i++) {
        int nc;
        // coverity[tainted_argument]
        ret = fscanf(ix, "%d", &nc); ASSERT_ALWAYS(ret == 1);
        for(unsigned int cg = 0 ; cg < nchars2 ; cg+=64) {
            res[i][cg] = 0;
        }
        for(int k = 0 ; k < nc ; k++) {
            unsigned int col;
            ret = fscanf(ix, "%x", &col); ASSERT_ALWAYS(ret == 1);
            ASSERT_ALWAYS(col < bcmat.nrows());
            for(unsigned int cg = 0 ; cg < nchars2 ; cg+=64) {
                res[i][cg] ^= bcmat[col][cg];
            }
        }
    }
    fclose_maybe_compressed(ix, indexname);
    return res;
}

/* We support both ascii and binary format, which is close to a bug. */

static blockmatrix
read_heavyblock_matrix_binary(unsigned int exp_nrows, const char * heavyblockname)
{
    FILE * f = fopen(heavyblockname, "rb");

    if (f == NULL) {
        fprintf(stderr, "Warning: %s not found, assuming empty\n", heavyblockname);
        return blockmatrix(0,0);
    }

    unsigned int nrows, ncols;

    /* If we're binary, we insist on having the companion files as well,
     * which provide a quick hint at the matrix dimensions */

    {
        char * rwname = derived_filename(heavyblockname, "rw", ".bin");
        struct stat sbuf[1];
        int const rc = stat(rwname, sbuf);
        if (rc < 0) { perror(rwname); exit(1); }
        nrows = sbuf->st_size / sizeof(uint32_t);
        free(rwname);
    }
    {
        char * cwname = derived_filename(heavyblockname, "cw", ".bin");
        struct stat sbuf[1];
        int const rc = stat(cwname, sbuf);
        if (rc < 0) { perror(cwname); exit(1); }
        ncols = sbuf->st_size / sizeof(uint32_t);
        free(cwname);
    }

    ASSERT_ALWAYS(nrows == exp_nrows);

    blockmatrix res(nrows, ncols);

    /* Sometimes the heavy block width is not a multiple of 64. Thus we pad
     * with zeros */
    res.set_zero();

    for(unsigned int i = 0 ; i < nrows ; i++) {
        uint32_t len;
        int r = fread32_little(&len, 1, f);
        ASSERT_ALWAYS(r == 1);
        for( ; len-- ; ) {
            uint32_t v;
            r = fread32_little(&v, 1, f); ASSERT_ALWAYS(r == 1);
            res[i][v-(v%64)] ^= ((uint64_t)1) << (v%64);
        }
    }
    fclose (f);
    return res;
}

static blockmatrix
read_heavyblock_matrix_ascii(unsigned int exp_nrows, const char * heavyblockname)
{
    FILE * f = fopen(heavyblockname, "r");

    if (f == NULL) {
        fprintf(stderr, "Warning: %s not found, assuming empty\n", heavyblockname);
        return blockmatrix(0,0);
    }

    unsigned int nrows, ncols;
    int const rc = fscanf(f,"%u %u", &nrows, &ncols);
    ASSERT_ALWAYS(rc == 2);
    ASSERT_ALWAYS(nrows == exp_nrows);

    blockmatrix res(nrows, ncols);

    /* Sometimes the heavy block width is not a multiple of 64. Thus we pad
     * with zeros */
    res.set_zero();

    for(unsigned int i = 0 ; i < nrows ; i++) {
        uint32_t len;
        int r = fscanf(f, "%" SCNu32, &len); ASSERT_ALWAYS(r == 1);
        for( ; len-- ; ) {
            uint32_t v;
            r = fscanf(f, "%" SCNu32, &v); ASSERT_ALWAYS(r == 1);
            res[i][v] ^= ((uint64_t)1) << (v%64);
        }
    }
    fclose (f);
    return res;
}

static blockmatrix
read_heavyblock_matrix (unsigned int nrows, const char * heavyblockname)
{
  if (heavyblockname == NULL)
    return blockmatrix(nrows, 0);
  else if (has_suffix(heavyblockname, ".bin"))
    return read_heavyblock_matrix_binary(nrows, heavyblockname);
  else
    return read_heavyblock_matrix_ascii(nrows, heavyblockname);
}

static int compute_transpose_of_blockmatrix_kernel(blockmatrix & kb, blockmatrix & t)
{
    /* gauss.c's kernel() function takes its input with a different ordering.
     * It's tiny data anyway. */

    fprintf(stderr, "Computing left nullspace of %u x %u matrix\n",
            t.nrows(), t.ncols());
    unsigned int const tiny_nrows = t.nrows();
    unsigned int const tiny_ncols = t.ncols();
    unsigned int const tiny_limbs_per_row = iceildiv(tiny_ncols, 64);
    unsigned int const tiny_limbs_per_col = iceildiv(tiny_nrows, 64);
    unsigned int const tiny_chars = FLAT_BYTES_WITH_READAHEAD(t.nrows(), t.ncols());
    unsigned int const tiny_64bit_words = tiny_chars / sizeof(uint64_t);

    /* we need some readahead zones because of the block matrix structure */
    uint64_t * tiny = (uint64_t *) malloc (tiny_chars);
    memset(tiny, 0, tiny_chars);
    
    blockmatrix::copy_to_flat(tiny, tiny_limbs_per_row, t);

    /* The kernel matrix is essentially a square matrix of tiny_nrows rows and
     * columns (tiny_nrows is the same as the number of kernel vectors)
     */
    /* we need some readahead zones because of the block matrix structure */
    uint64_t * kerdata = (uint64_t *) malloc(FLAT_BYTES_WITH_READAHEAD(t.nrows(), t.nrows()));
    memset(kerdata, 0, FLAT_BYTES_WITH_READAHEAD(t.nrows(), t.nrows()));
    unsigned int const kerdata_64bit_words = FLAT_BYTES_WITH_READAHEAD(t.nrows(), t.nrows()) / sizeof(uint64_t);


    uint64_t ** myker = (uint64_t **) malloc(tiny_nrows * sizeof(uint64_t *));
    ASSERT(myker != NULL);
    for (unsigned int i = 0; i < tiny_nrows; ++i)
        myker[i] = kerdata + i * tiny_limbs_per_col;
    /* gauss.c knows about mp_limb_t's only */
    ASSERT_ALWAYS(sizeof(uint64_t) % sizeof(mp_limb_t) == 0);
    blockmatrix::swap_words_if_needed (tiny, tiny_64bit_words);
    int const dim = kernel((mp_limb_t *) tiny,
            (mp_limb_t **) myker,
            tiny_nrows, tiny_ncols,
            sizeof(uint64_t) / sizeof(mp_limb_t) * tiny_limbs_per_row,
            sizeof(uint64_t) / sizeof(mp_limb_t) * tiny_limbs_per_col);
    blockmatrix::swap_words_if_needed (tiny, tiny_64bit_words); /* FIXME: this is maybe not
                                                      needed since tiny is
                                                      destroyed, but keep it
                                                      for debugging */
    blockmatrix::swap_words_if_needed(kerdata, kerdata_64bit_words);
    free(tiny);
    /* Now take back our kernel to block format, and multiply. Exciting. */
    blockmatrix::copy_transpose_from_flat(kb, kerdata, tiny_limbs_per_col);
    free(myker);
    free(kerdata);
    return dim;
}

/* This only builds a basis, not an echelonized basis */
static void blockmatrix_column_reduce(blockmatrix & m, unsigned int max_rows_to_consider)
{
    auto t = m.view(submatrix_range(0, 0, MIN(max_rows_to_consider, m.nrows()), m.ncols()));

    unsigned int const tiny_nrows = t.ncols();
    unsigned int const tiny_ncols = t.nrows();
    unsigned int const tiny_limbs_per_row = iceildiv(tiny_ncols, 64);
    unsigned int const tiny_limbs_per_col = iceildiv(tiny_nrows, 64);

    unsigned int const tiny_nlimbs = tiny_nrows * tiny_limbs_per_row;
    uint64_t * tiny = (uint64_t *) malloc(tiny_nlimbs * sizeof(uint64_t));
    memset(tiny, 0, tiny_nlimbs * sizeof(uint64_t));
    blockmatrix::copy_transpose_to_flat(tiny, tiny_limbs_per_row, t);

    uint64_t * sdata = (uint64_t *) malloc(tiny_nrows * tiny_limbs_per_col * sizeof(uint64_t));
    memset(sdata, 0, tiny_nrows * tiny_limbs_per_col * sizeof(uint64_t));
    unsigned int const sdata_64bit_words = tiny_nrows * tiny_limbs_per_col;

    blockmatrix::swap_words_if_needed (tiny, tiny_nlimbs);
    int const rank = spanned_basis(
            (mp_limb_t *) sdata,
            (mp_limb_t *) tiny,
            tiny_nrows,
            tiny_ncols,
            sizeof(uint64_t) / sizeof(mp_limb_t) * tiny_limbs_per_row,
            sizeof(uint64_t) / sizeof(mp_limb_t) * tiny_limbs_per_col,
            NULL
            );
    blockmatrix::swap_words_if_needed (tiny, tiny_nlimbs);
    blockmatrix::swap_words_if_needed (sdata, sdata_64bit_words);
    free(tiny);

    blockmatrix s(m.ncols(), rank);
    blockmatrix::copy_transpose_from_flat(s, sdata, tiny_limbs_per_col);
    blockmatrix k2(m.nrows(), rank);
    blockmatrix::mul_smallb(k2, m, s);
    std::swap(k2, m);
    free(sdata);
}

static void
declare_usage (cxx_param_list & pl)
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
  param_list_decl_usage (pl, "ker",    "input kernel file");
  param_list_decl_usage (pl, "nratchars", "number of characters on rational "
                                          "side");
  param_list_decl_usage(pl, "large-ab", "enable support for a and b larger than"
                                        "64 bits");
  param_list_decl_usage (pl, "only-sign-chars", "use only the sign character "
                                                "on each side");
}

// coverity[root_function]
int main(int argc, char const * argv[])
{
    const char * heavyblockname = NULL;
    int nchars, nratchars = 0;
    cxx_cado_poly cpoly;
    const char *purgedname = NULL;
    const char *indexname = NULL;
    const char *outname = NULL;
    int nthreads = 1;
    unsigned long lpb[2] = {0,0};
    bool only_sign_chars = false;
    double const cpu0 = seconds ();
    double const wct0 = wct_seconds ();

    /* print the command line */
    fmt::print(stderr, "# ({}) {}\n",
            cado_revision_string,
            collect_command_line(argc, argv));

    cxx_param_list pl;
    declare_usage(pl);

    const char *bw_kernel_file = NULL;

    pl.configure_switch("large-ab");
    pl.configure_switch("only-sign-chars");

    cado::filter_io_details::configure(pl);

    param_list_process_command_line(pl, &argc, &argv, false);

    purgedname = param_list_lookup_string(pl, "purged");
    indexname = param_list_lookup_string(pl, "index");
    outname = param_list_lookup_string(pl, "out");
    heavyblockname = param_list_lookup_string(pl, "heavyblock");
    bw_kernel_file = param_list_lookup_string(pl, "ker");

    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "poly")) == NULL)
        pl.fail("Error: parameter -poly is mandatory\n");
    cpoly.read(tmp);
    if (cpoly.nsides() < 1 || cpoly.nsides() > 2)
        pl.fail("Error: number of polys should be 1 or 2, got {}\n",
                           cpoly.nsides());

    if (param_list_parse_int(pl, "nchar", &nchars) == 0)
        pl.fail("Error: parameter -nchar is mandatory\n");

    if (only_sign_chars && nchars != cpoly.nsides()) {
        fmt::print(stderr, "Error: with -only-sign-chars, -nchars should be "
                           "equal to the number of sides ({}), got {}\n",
                           cpoly.nsides(), nchars);
        exit (EXIT_FAILURE);
    }

    pl.parse("only_sign_chars", only_sign_chars);

    /* parse the optional -nratchars option */
    param_list_parse_int(pl, "nratchars", &nratchars);
    if (nratchars && cpoly.get_ratside() == -1) { /* no rat side */
        fmt::print(stderr, "Error: nratchars is non-zero ({}) but poly has no "
                           "rational side\n", nratchars);
        exit (EXIT_FAILURE);
    }
    /* parse lpb{side} for each side of cpoly */
    for (int side = 0; side < cpoly.nsides(); ++side) {
        std::string arg = fmt::format("lpb{}", side);
        if (param_list_parse_ulong(pl, arg.c_str(), &lpb[side]) == 0)
            pl.fail("Error: parameter {} is mandatory\n", arg);
    }

    param_list_parse_int(pl, "t", &nthreads);

    if (purgedname == NULL || indexname == NULL || outname == NULL)
        pl.fail("Error: parameters -purged, -index and -out are mandatory\n");

    /* Put nchars characters on all algebraic sides and nratchars on the
     * rational side (if it exists).
     */
    std::vector<alg_prime_t> chars;
    if (only_sign_chars) {
        chars = create_sign_characters_only(cpoly);
    } else {
        std::vector<int> nch(cpoly.nsides());
        for (int side = 0; side < cpoly.nsides(); ++side) {
            nch[side] = cpoly[side]->deg > 1 ? nchars : nratchars;
        }
        chars = create_characters (nch, cpoly, lpb);
    }

    cado::filter_io_details::interpret_parameters(pl);

    blockmatrix bcmat = big_character_matrix(chars, purgedname,
            cpoly, nthreads, pl.parse<bool>("large-ab"));

    fprintf(stderr, "done building big character matrix at wct=%.1fs\n", wct_seconds()-wct0);

    blockmatrix const scmat = small_character_matrix(bcmat, indexname);
    fprintf(stderr, "done building small character matrix at wct=%.1fs\n", wct_seconds()-wct0);

    unsigned int const small_nrows = scmat.nrows();

    /* It's ok if heavyblockname == 0. After all sufficiently many
     * characters should be enough */
    blockmatrix const h = read_heavyblock_matrix(small_nrows, heavyblockname);
    if (h.ncols())
          fprintf(stderr, "done reading heavy block of size %u x %u at wct=%.1fs\n",
                  h.nrows(), h.ncols(), wct_seconds()-wct0);
    ASSERT_ALWAYS(h.nrows() == small_nrows);

    /* Now do dot products of these matrices by the kernel vectors
     * supplied on input */

    /* First compute how many kernel vectors we have */

    unsigned int total_kernel_cols = 0;
    unsigned int kncols;
    struct stat sbuf[1];
    int const rc = stat(bw_kernel_file, sbuf);
    if (rc < 0) { perror(bw_kernel_file); exit(1); }
    ASSERT_ALWAYS(sbuf->st_size % small_nrows == 0);
    unsigned int const ncols = 8 * (sbuf->st_size / small_nrows);
    fprintf(stderr, "%s: %u vectors\n", bw_kernel_file, ncols);
    ASSERT_ALWAYS(ncols % 64 == 0);
    total_kernel_cols += kncols = ncols;

    fprintf(stderr, "Total: %u kernel vectors\n", total_kernel_cols);

    /* kmat is the join of all kernel vectors */
    blockmatrix k(small_nrows, total_kernel_cols);
    int const j0 = 0;
    k.read_from_flat_file(0, j0, bw_kernel_file, small_nrows, kncols);

    fprintf(stderr, "done reading %u kernel vectors at wct=%.1fs\n",
            total_kernel_cols, wct_seconds() - wct0);

    {
        blockmatrix_column_reduce(k, 4096);
    }
    total_kernel_cols = k.ncols();
    fprintf(stderr, "Info: input kernel vectors reduced to dimension %u\n",
            total_kernel_cols);

    /* tmat is the product of the character matrices times the kernel vector */
    blockmatrix t(total_kernel_cols, scmat.ncols() + h.ncols());
    blockmatrix::mul_Ta_b(t.view(submatrix_range(0, 0, total_kernel_cols, scmat.ncols())), k, scmat);
    blockmatrix::mul_Ta_b(t.view(submatrix_range(0, scmat.ncols(), total_kernel_cols, h.ncols())), k, h);
    fprintf(stderr, "done multiplying matrices at wct=%.1fs\n", wct_seconds() - wct0);

    // blockmatrix_free(scmat);
    // blockmatrix_free(h);


    blockmatrix kb(k.ncols(), k.ncols());
    kb.set_zero();
    int const dim = compute_transpose_of_blockmatrix_kernel(kb, t);
    // blockmatrix_free(t);
    
    fprintf(stderr, "dim of ker = %d\n", dim);

    auto kbsub = kb.view(submatrix_range(0, 0, kb.nrows(), dim));

    blockmatrix::mul_smallb(k, kbsub);

    // blockmatrix nk(small_nrows, kbsub.ncols());
    // blockmatrix::mul_smallb(nk, k, kbsub);
    // free(kbsub);
    // blockmatrix_free(k);
    // blockmatrix_free(kb);

    /* Sanity check: count the number of zero dependencies */
    unsigned int nonzero_deps = 0;
    for(unsigned int j = 0 ; j < k.ncols() ; j+=64) {
        uint64_t sanity = 0 ;
        for(unsigned int i = 0 ; i < k.nrows() ; i+=64) {
                sanity |= k[i][j];
        }
        // do popcount.
        for( ; sanity ; sanity>>=1) nonzero_deps += sanity&1UL;
    }
    if (!nonzero_deps) {
        fprintf(stderr, "Error, all dependencies are zero !\n");
        exit(1);
    }
    k.write_to_flat_file(outname, 0, 0, k.nrows(), k.ncols());

    fprintf(stderr, "Wrote %d non-zero dependencies to %s\n",
            nonzero_deps, outname);
    if (nonzero_deps < (unsigned int) dim || k.ncols() % 64) {
        fprintf(stderr, "This includes %u discarded zero dependencies, as well as %u padding zeros\n",
                dim - nonzero_deps,
                k.ncols_padded() - dim);
    }
    // blockmatrix_free(k);

    /* print total time and memory usage */
    print_timing_and_memory (stdout, cpu0, wct0);

    return 0;
}
