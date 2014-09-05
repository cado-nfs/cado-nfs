#include "cado.h"

#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <list>
#include <string.h>
#include <cstdio>
#include <gmp.h>
#include <errno.h>
#include <unistd.h>
#include <ctype.h>
#include <utility>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <functional>
#include <iomanip>
#include <sstream>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include "bwc_config.h"
#include "macros.h"
#include "utils.h"
#include "bw-common.h"
#include "filenames.h"

#include "gf2x-fft.h"
#include "lingen_mat_types.hpp"

/* Provide workalikes of usual interfaces for some ungifted systems */
#include "portability.h"

/* Name of the source a file */
char input_file[FILENAME_MAX];

/* threshold for the recursive algorithm */
unsigned int lingen_threshold = 0;

/* threshold for cantor fft algorithm */
unsigned int cantor_threshold = UINT_MAX;

/* {{{ macros used here only -- could be bumped to macros.h if there is
 * need.
 */
#define NEVER_HAPPENS(tst, action) do {					\
    if (UNLIKELY(tst)) {						\
        fprintf(stderr, "under attack by martians\n");	        	\
        action								\
    }									\
} while (0)

#define WARN_ERRNO_DIAG(tst, func, arg) do {				\
    if (UNLIKELY(tst)) {						\
        fprintf(stderr, "%s(%s): %s\n", func, arg, strerror(errno));	\
    }									\
} while (0)
/* }}} */

/* output: F0x is the candidate number x. All coordinates of the
 * candidate are grouped, and coefficients come in order (least
 * significant first), separated by a carriage return.
 */

// applies a permutation on the source indices/*{{{*/
template<typename T>
void permute(std::vector<T>& a, unsigned int p[])
{
    std::vector<T> b(a.size(), (unsigned int) -1);
    for(unsigned int i = 0 ; i < a.size() ; i++) {
        b[p[i]] = a[i];
    }
    a.swap(b);
}/*}}}*/


/* Needed operations on critical path.
 *
 * compute_ctaf: compute the coefficient of degree t of the product of E
 * (m*b matrix) by PI (a b*b matrix), to be stored into small-e (m*b
 * matrix)
 *
 * zeroing out a column of small-e (m*b)
 *
 * do a column gaussian elimination on small-e:
 *   - column after column,
 *     - find the first non-zero bit in the column.
 *     - add this column to other column.
 *
 * Then polynomial operations.
 *
 * => small-e is better stored column-wise.
 * => big-E as well, and PI also.
 *
 * Extra un-critical stuff includes
 * - initialization -> computation of F0
 *                  -> computation of initial E = A * F0
 * - termination -> computation of final F = F0 * PI
 *
 * (while uncritical, the computation of the initial E and the final F
 * may become a problem for debugging situations when the block size is
 * largish).
 *
 * all this stuff has linear complexity.
 */

    /*{{{*/ /* debug code for polmat */
void dbmat(bmat const * x)
{
    using namespace std;
    for(unsigned int i = 0 ; i < x->nrows ; i++) {
        for(unsigned int j = 0 ; j < x->ncols ; j++) {
            std::cout << x->coeff(i,j);
        }
        std::ostringstream w;
        unsigned int j = 0;
        for( ; j < x->ncols ; ) {
            unsigned long z = 0;
            for(unsigned int k = 0 ; j < x->ncols && k < ULONG_BITS ; k++,j++) {
                // z <<= 1;
                // z |= x->coeff(i,j);
                // This fancy bit positioning matches the syntax of
                // read_hexstring and write_hexstring
                z |= ((unsigned long) x->coeff(i,j)) << (ULONG_BITS-4-k+((k&3)<<1));
            }
            w << setw(ULONG_BITS/4) << setfill('0') << hex << z;
        }
        w.flush();
        std::cout << ":" << w.str() << "\n";
    }
}
void dpmat(polmat const * pa, unsigned int k)
{
    bmat x;
    pa->extract_coeff(x, k);
    dbmat(&x);
}/*}}}*/

namespace globals {
    unsigned int m,n;
    unsigned int t,t0,sequence_length;
    std::vector<unsigned int> delta;
    std::vector<unsigned int> chance_list;

    polmat E;
#ifndef NDEBUG
    polmat E_saved;
#endif
    bmat e0;

    // F0 is exactly the n x n identity matrix, plus the X^(s-exponent)e_{cnum}
    // vectors. Here we store the cnum,exponent pairs.
    std::vector<std::pair<unsigned int, unsigned int> > f0_data;

    double start_time;
}

// To multiply on the right an m x n matrix A by F0, we start by copying
// A into the first n columns. Since we're also dividing out by X^t0, the
// result has to be shifted t0 positions to the right.
// Afterwards, column n+j of the result is column cnum[j] of A, shifted
// exponent[j] positions to the right.
void compute_E_from_A(polmat const &a)/*{{{*/
{
    using namespace globals;
    polmat tmp_E(m, m + n, a.ncoef - t0);
    for(unsigned int j = 0 ; j < n ; j++) {
        tmp_E.import_col_shift(j, a, j, - (int) t0);
    }
    for(unsigned int j = 0 ; j < m ; j++) {
        unsigned int cnum = f0_data[j].first;
        unsigned int exponent = f0_data[j].second;
        tmp_E.import_col_shift(n + j, a, cnum, - (int) exponent);
    }
    E.swap(tmp_E);
    E.clear_highbits();
}/*}}}*/

#ifdef  DO_EXPENSIVE_CHECKS     /* {{{ */
/* revive this old piece of code from old commits. */
void multiply_slow(polmat& dst, polmat /* const */ &a, polmat /* const */ &b)
{
    BUG_ON(a.ncols != b.nrows);

    unsigned int Gaw = (a.ncoef+ULONG_BITS-1)/ULONG_BITS;
    unsigned int Gbw = (b.ncoef+ULONG_BITS-1)/ULONG_BITS;

    unsigned long tmp[Gaw + Gbw];

    polmat res(a.nrows, b.ncols, (Gaw + Gbw) * ULONG_BITS);
    for (uint j = 0; j < b.ncols; j++) {
        for (uint i = 0; i < a.nrows; i++) {
            for (uint k = 0; k < a.ncols; k++) {
                int da = a.deg(k);
                int db = b.deg(j);
                if (da == -1 || db == -1)
                    continue;
                /* const */ unsigned long * ax = a.poly(i,k);
                /* const */ unsigned long * bx = b.poly(k,j);
                unsigned int aw = 1 + da / ULONG_BITS;
                unsigned int bw = 1 + db / ULONG_BITS;
                unsigned int az = (da+1) & (ULONG_BITS-1);
                unsigned int bz = (db+1) & (ULONG_BITS-1);

                /* It's problematic if the high bits aren't cleared, as
                 * gf2x works at the word level */
                /* becuase import_col_shift isn't careful enough, we have
                 * to do the cleanup by ourself.
                 */
                if (az) { ax[aw-1] &= (1UL << az) - 1; }
                if (bz) { bx[bw-1] &= (1UL << bz) - 1; }
                // if (az) { ASSERT((ax[aw-1]>>az)==0); }
                // if (bz) { ASSERT((bx[bw-1]>>bz)==0); }
                mul_gf2x(tmp,ax,aw,bx,bw);
                
                for(unsigned int w = 0 ; w < aw+bw ; w++) {
                    res.poly(i,j)[w] ^= tmp[w];
                }
            }
        }
    }
    for (uint j = 0; j < b.ncols; j++) {
        res.setdeg(j);
    }
    dst.swap(res);
}

struct checker {
    polmat e;
    unsigned int t0;
    checker() {
        e.copy(globals::E_saved);
        t0=globals::t;
    }
    int check(polmat& pi) {
        using namespace globals;
        printf("Checking %u..%u\n",t0,t);
        multiply_slow(e,e,pi);
        unsigned int v = e.valuation();
        ASSERT(v == t-t0);
        e.xdiv_resize(t-t0, e.ncoef-(t-t0));
        return v == t-t0;
    }
};
#endif /* }}} */

/* {{{ utility */
template<typename iterator>
std::string intlist_to_string(iterator t0, iterator t1)
{
    std::ostringstream cset;
    for(iterator t = t0 ; t != t1 ; ) {
        iterator u;
        for(u = t ; u != t1 ; u++) {
            if ((typename std::iterator_traits<iterator>::difference_type) (*u - *t) != (u-t)) break;
        }
        if (t != t0) cset << ',';
        if (u-t == 1) {
            cset << *t;
        } else {
            cset << *t << "-" << u[-1];
        }
        t = u;
    }
    return cset.str();
}
/* }}} */

// F is in fact F0 * PI.
// To multiply on the *left* by F, we cannot work directly at the column
// level (we could work at the row level, but it does not fit well with
// the way data is organized). So it's merely a matter of adding
// polynomials.
void compute_final_F_from_PI(polmat& F, polmat const& pi)/*{{{*/
{
    printf("Computing final F from PI\n");
    using namespace globals;
    // We take t0 rows, so that we can do as few shifts as possible
    // tmpmat is used only within the inner loop.
    polmat tmpmat(t0,1,pi.ncoef);
    polmat tmp_F(n, m + n, globals::t0 + pi.ncoef);
    using namespace std;
    for(unsigned int i = 0 ; i < n ; i++) {
        // What contributes to entries in this row ?
        vector<pair<unsigned int, unsigned int> > l;
        set<unsigned int> sexps;
        // l.push_back(make_pair(i,0));
        for(unsigned int j = 0 ; j < m ; j++) {
            if (f0_data[j].first == i) {
                l.push_back(make_pair(n + j, f0_data[j].second));
                sexps.insert(f0_data[j].second);
            }
        }
        vector<unsigned int> exps(sexps.begin(), sexps.end());
        // So the i-th row of f0 has a 1 at position i, and then
        // X^(t0-f.second) at position n+j whenever f.first == i.
        //
        // Now fill in the row.
        for(unsigned int j = 0 ; j < m + n ; j++) {
            // tmpmat is a single column.
            tmpmat.zcol(0);
            for(unsigned int k = 0 ; k < l.size() ; k++) {
                tmpmat.addpoly(l[k].second, 0, pi.poly(l[k].first, j));
            }
            tmp_F.addpoly(i,j,pi.poly(i,j));
            for(unsigned int k = 0 ; k < exps.size() ; k++) {
                tmpmat.xmul_poly(exps[k],0,t0-exps[k]);
                tmp_F.addpoly(i,j,tmpmat.poly(exps[k],0));
            }
        }
    }
    for(unsigned int j = 0 ; j < m + n ; j++) {
        int pideg = pi.deg(j);
        tmp_F.deg(j) = t0 + pideg;
        if (pideg == -1) {
            tmp_F.deg(j) = -1;
        }
    }
    F.swap(tmp_F);
}/*}}}*/

/* {{{ main input ; read_data_for_series */
/* ondisk_length denotes something which might be larger than the
 * actually required sequence length, because the krylov computation is
 * allowed to run longer
 */
void read_data_for_series(polmat& A MAYBE_UNUSED, unsigned int ondisk_length)
{
    using namespace globals;

    /* It's slightly non-trivial. polmat entries are polynomials, so
     * we've got to cope with the non-trivial striding */
    
    { /* {{{ check file size */
        struct stat sbuf[1];
        int rc = stat(input_file, sbuf);
        DIE_ERRNO_DIAG(rc<0,"stat",input_file);
        ssize_t expected = m * n / CHAR_BIT * ondisk_length;

        if (sbuf->st_size != expected) {
            fprintf(stderr, "%s does not have expected size %zu\n",
                    input_file, expected);
            exit(1);
        }
    } /* }}} */

    /* We've got matrices of m times n bits. */
    FILE * f = fopen(input_file, "rb");
    DIE_ERRNO_DIAG(f == NULL, "fopen", input_file);


    ASSERT_ALWAYS(n % ULONG_BITS == 0);
    size_t ulongs_per_mat = m*n/ULONG_BITS;
    unsigned long * buf;
    size_t rz;

    buf = (unsigned long *) malloc(ulongs_per_mat * sizeof(unsigned long));

    printf("Using A(X) div X in order to consider Y as starting point\n");
    rz = fread(buf, sizeof(unsigned long), ulongs_per_mat, f);
    ASSERT_ALWAYS(rz == ulongs_per_mat);
    sequence_length--;

    polmat a(m, n, sequence_length);

    /* We hardly have any other option than reading the stuff bit by bit
     * :-(( */

    unsigned long mask = 1UL;
    unsigned int offset = 0;
    for(unsigned int k = 0 ; k < sequence_length ; k++) {
        rz = fread(buf, sizeof(unsigned long), ulongs_per_mat, f);
        ASSERT_ALWAYS(rz == ulongs_per_mat);

        unsigned long * v = buf;
        unsigned long lmask = 1UL;

        for(unsigned int i = 0 ; i < m ; i++) {
            for(unsigned int j = 0 ; j < n ; j++) {
                a.poly(i,j)[offset] |= mask & -((*v & lmask) != 0);
                lmask <<= 1;
                v += (lmask == 0);
                lmask += (lmask == 0);
            }
        }

        mask <<= 1;
        offset += mask == 0;
        mask += mask == 0;
    }
    for (unsigned int j = 0; j < n; j ++) {
        a.deg(j) = sequence_length - 1;
    }

    free(buf);
    fclose(f);

    A.swap(a);
}
/* }}} */


/* {{{ main output ; bw_commit_f */
void bw_commit_f(polmat& F)
{
    using namespace globals;
    using namespace std;

    /* Say n' columns are said interesting. We'll pad these to n */
    printf("Writing F files\n");

    unsigned int * pick = (unsigned int *) malloc(n * sizeof(unsigned int));
    memset(pick, 0, n * sizeof(unsigned int));
    unsigned int nres = 0;

    unsigned int ncoef = 0;
    for (unsigned int j = 0; j < m + n; j++) {
        if (chance_list[j]) {
            ASSERT_ALWAYS(F.deg(j) != -1);
            ncoef = std::max(ncoef, (unsigned int) (1 + F.deg(j)));
            ASSERT_ALWAYS(nres < n);
            pick[nres++] = j;
        }
    }
    string s = intlist_to_string(pick, pick + nres);
    printf("Picking columns %s\n", s.c_str());

    ASSERT_ALWAYS(n % ULONG_BITS == 0);
    size_t ulongs_per_mat = n*n/ULONG_BITS;
    unsigned long * buf;
    size_t rz;

    FILE * f = fopen(LINGEN_F_FILE, "wb");
    DIE_ERRNO_DIAG(f == NULL, "fopen", LINGEN_F_FILE);

    buf = (unsigned long *) malloc(ulongs_per_mat * sizeof(unsigned long));

    /* We're writing the reversal of the F polynomials in the output
     * file. The reversal is taken column-wise. Since not all delta_i's
     * are equal, we have to take into account the possibility that the
     * data written to disk corresponds to possibly unrelated terms in
     * the polynomial expansion of F.
     */

    unsigned long * fmasks;
    unsigned int * foffsets = 0;
    fmasks = (unsigned long *) malloc(nres  * sizeof(unsigned long));
    foffsets = (unsigned int *) malloc(nres  * sizeof(unsigned int));

    for(unsigned int i = 0 ; i < nres ; i++) {
        unsigned int ii = pick[i];
        unsigned int d = delta[ii];

        fmasks[i] = 1UL << (d % ULONG_BITS);
        foffsets[i] = d / ULONG_BITS;
    }

    for(unsigned int k = 0 ; k < ncoef ; k++) {
        memset(buf, 0, ulongs_per_mat * sizeof(unsigned long));

        unsigned long * v = buf;
        unsigned long lmask = 1UL;

        for(unsigned int i = 0 ; i < nres ; i++) {
            for(unsigned int j = 0 ; j < n ; j++) {
                if (foffsets[i] != UINT_MAX)
                    *v |= lmask & -((F.poly(j,pick[i])[foffsets[i]] & fmasks[i]) != 0);
                lmask <<= 1;
                v += (lmask == 0);
                lmask += (lmask == 0);
            }
        }

        rz = fwrite(buf, sizeof(unsigned long), ulongs_per_mat, f);
        ASSERT_ALWAYS(rz == ulongs_per_mat);

        for(unsigned int i = 0 ; i < nres ; i++) {
            if (foffsets[i] == UINT_MAX)
                continue;

            fmasks[i] >>= 1;
            foffsets[i] -= fmasks[i] == 0;
            fmasks[i] += ((unsigned long) (fmasks[i] == 0)) << (ULONG_BITS-1);
        }
    }
    fclose(f);
    free(fmasks);
    free(foffsets);
    free(pick);
    free(buf);
}
/*}}}*/


/* {{{ polmat I/O ; this concerns only the pi matrices */
/* This changes from the format used pre-bwc. */
void write_polmat(polmat const& P, const char * fn)/*{{{*/
{
    FILE * f = fopen(fn, "wb");
    if (f == NULL) {
        perror(fn);
        exit(1);
    }

    size_t rz;
    size_t nw = ( P.ncoef + ULONG_BITS - 1) / ULONG_BITS;

    fprintf(f, "%u %u %u\n", P.nrows, P.ncols, P.ncoef);
    for(unsigned int j = 0 ; j < P.ncols ; j++) {
        fprintf(f, "%s%d", j?" ":"", P.deg(j));
    }
    fprintf(f, "\n");

    fprintf(f, "%zu %zu\n", nw, sizeof(unsigned long));
    for(unsigned int i = 0 ; i < P.nrows ; i++) {
        for(unsigned int j = 0 ; j < P.ncols ; j++) {
            rz = fwrite(P.poly(i,j), sizeof(unsigned long), nw, f);
            ASSERT_ALWAYS(rz == nw);
        }
    }
    fclose(f);
}

/*}}}*/

void write_pi(polmat const& P, unsigned int t1, unsigned int t2)
{
    char * tmp;
    int rc = asprintf(&tmp, LINGEN_PI_PATTERN, t1, t2);
    ASSERT_ALWAYS(rc >= 0);
    write_polmat(P, tmp);
    // printf("written %s\n", tmp);
    free(tmp);
}

void unlink_pi(unsigned int t1, unsigned int t2)
{
    char * tmp;
    int rc = asprintf(&tmp, LINGEN_PI_PATTERN, t1, t2);
    ASSERT_ALWAYS(rc >= 0);
    rc = unlink(tmp);
    WARN_ERRNO_DIAG(rc < 0, "unlink", tmp);
    free(tmp);
}
/* }}} */

bool recover_f0_data()/*{{{*/
{
    std::ifstream f(LINGEN_BOOTSTRAP_FILE);

    using namespace globals;
    for(unsigned int i = 0 ; i < m ; i++) {
        unsigned int exponent,cnum;
        if (!(f >> cnum >> exponent))
            return false;
        f0_data.push_back(std::make_pair(cnum,exponent));
    }
    printf("recovered " LINGEN_BOOTSTRAP_FILE " from disk\n");
    return true;
}/*}}}*/

bool write_f0_data()/*{{{*/
{
    std::ofstream f(LINGEN_BOOTSTRAP_FILE);

    using namespace globals;
    for(unsigned int i = 0 ; i < m ; i++) {
        unsigned int cnum = f0_data[i].first;
        unsigned int exponent = f0_data[i].second;
        if (!(f << " " << cnum << " " << exponent))
            return false;
    }
    f << std::endl;
    printf("written " LINGEN_BOOTSTRAP_FILE " to disk\n");
    return true;
}/*}}}*/

// the computation of F0 says that the column number cnum of the
// coefficient X^exponent of A increases the rank.
//
// These columns are gathered into an invertible matrix at step t0, t0
// being strictly greater than the greatest exponent above. This is so
// that there is no trivial column dependency in F0.
void set_t0_delta_from_F0()/*{{{*/
{
    using namespace globals;
    t0 = f0_data.back().second + 1;     // see above for the +1
    for(unsigned int j = 0 ; j < m + n ; j++) {
        delta[j] = t0;
        // Even irrespective of the value of the exponent, we get t0 for
        // delta.
    }
}/*}}}*/

void compute_f_init(polmat& A)/*{{{*/
{
    using namespace globals;

    /* build F_INIT */
    printf("Computing t0\n");

    /* For each integer i between 0 and m-1, we have a column, picked
     * from column cnum[i] of coeff exponent[i] of A which, once reduced modulo
     * the other ones, has coefficient at row pivots[i] unequal to zero.
     */
    /* We can't use a VLA here because it's not in C++ for non-POD types,
     * and we don't want to open up the can of worms of allowing vector<>
     * here either, because that would mean allowing copy. We care about
     * forbidding copies entirely.
     */
    bcol * pcols = new bcol[m];
    unsigned int pivots[m], exponent[m], cnum[m];

    unsigned int r = 0;
    for(unsigned int k=0;r < m && k<A.ncoef;k++) {
        bmat amat;
        A.extract_coeff(amat, k);
        for(unsigned int j=0;r < m && j<n;j++) {
            /* copy column j, coeff k into acol */
            bcol acol;
            amat.extract_col(acol, j);
            /* kill as many coeffs as we can */
            for(unsigned int v=0;v<r;v++) {
                unsigned int u=pivots[v];
                /* the v-th column in the matrix reduced_rank is known to
                 * kill coefficient u (more exactly, to have a -1 as u-th
                 * coefficient, and zeroes for the other coefficients
                 * referenced in the pivots[0] to pivots[v-1] indices).
                 */

                acol.add(pcols[v],acol.coeff(u));
            }
            unsigned int u = acol.ffs();
            if (u == UINT_MAX) {
                printf("[X^%d] A, col %d does not increase rank (still %d)\n",
                        k,j,r);
                if (k * n > m + 40) {
                    printf("The choice of starting vectors was bad. "
                            "Cannot find %u independent cols within A\n",m);
                    exit(1);
                }
                continue;
            }

            /* Bingo, it's a new independent col. */
            pivots[r]=u;
            cnum[r]=j;
            exponent[r]=k;

            f0_data.push_back(std::make_pair(cnum[r], exponent[r]));

            /* TODO: For non-binary stuff, multiply the column so that
             * acol[u] becomes -1
             */
            pcols[r].swap(acol);
            r++;

            if (r == m)
                printf("[X^%d] A, col %d increases rank to %d (head row %d)\n",
                        k,j,r,u);
        }
    }
    delete[] pcols;

    t0 = exponent[r-1] + 1;
    printf("Found satisfying init data for t0=%d\n", t0);
                    
    /*
    printf("Init e0 matrix\n");
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            std::cout << A.coeff(i,j,t0);
        }
        for(unsigned int j = 0 ; j < m ; j++) {
            std::cout << A.coeff(i,cnum[j],exponent[j]);
        }
        std::cout << "\n";
    }
    */


    if (r!=m) {
        printf("This amount of data is insufficient. "
                "Cannot find %u independent cols within A\n",m);
        exit(1);
    }

}/*}}}*/

void print_deltas()/*{{{*/
{
    using namespace globals;
    printf("[t=%4u] delta =", t);
    unsigned int last = UINT_MAX;
    unsigned int nrep = 0;
    for(unsigned int i = 0 ; i < m + n ; i++) {
        unsigned int d = delta[i];
        if (d == last) {
            nrep++;
            continue;
        }
        // Flush the pending repeats
        if (last != UINT_MAX) {
            printf(" %u", last);
            if (nrep > 1)
                printf(" [%u]", nrep);
        }
        last = d;
        nrep = 1;
    }
    ASSERT_ALWAYS(last != UINT_MAX);
    printf(" %u", last);
    if (nrep > 1)
        printf(" [%u]", nrep);
    printf("\n");
}/*}}}*/

static void rearrange_ordering(polmat & PI, unsigned int piv[])/*{{{*/
{
    /* Sort the columns. It might seem merely cosmetic and useless to
     * sort w.r.t both the global and local nominal degrees. In fact, it
     * is crucial for the corectness of the computations. (Imagine a
     * 2-step increase, starting with uneven global deltas, and hitting
     * an even situation in the middle. One has to sort out the local
     * deltas to prevent trashing the whole picture).
     *
     * The positional sort, however, *is* cosmetic (makes debugging
     * easier).
     */
    using namespace std;
    using namespace globals;
    typedef pair<pair<unsigned int,unsigned int>, int> corresp_t;
    vector<corresp_t> corresp(m+n);
    for(unsigned int i = 0 ; i < m + n ; i++) {
        int pideg = PI.deg(i);
        if (pideg == -1) { /* overflowed ! */
            pideg = INT_MAX;
        }
        corresp[i] = make_pair(make_pair(delta[i],pideg), i);
    }
    sort(corresp.begin(), corresp.end(), less<corresp_t>());
    unsigned int p[m+n];
    for(unsigned int i = 0 ; i < m + n ; i++) {
        p[corresp[i].second] = i;
    }
    permute(delta, p);
    permute(chance_list, p);
    if (piv) {
        for(unsigned int i = 0 ; i < m ; i++) {
            piv[i] = p[piv[i]];
        }
    }
    PI.perm(p);
    e0.perm(p);
}/*}}}*/

static bool gauss(unsigned int piv[], polmat& PI)/*{{{*/
{
    /* Do one step of (column) gaussian elimination on e0. The columns
     * are assumed to be in the exact order that corresponds to the
     * ordering of the columns of PI (therefore relative to delta)
     */

    /* Note that we do *NOT* modify E. It seems to be the fastest option,
     * although the other deserves being investigated as well (see
     * old/lingen2.c, there are two versions of the quadratic algorithm).
     */
    unsigned int i,j,k;
    unsigned int rank;
    rank = 0 ;

    using namespace globals;

    /*
    std::cout << "Input matrix\n";
    dbmat(&e0);
    */

    std::vector<unsigned int> overflowed;

    for(j = 0 ; j < e0.ncols ; j++) {
        /* Find the pivot inside the column. */
        i = e0.ffs(j);
        if (i == UINT_MAX)
            continue;
        ASSERT(rank < e0.nrows && rank < e0.ncols);
        // std::cout << fmt("col % is the %-th pivot\n") % j % rank;
        piv[rank++] = j;
        /* Cancel this coeff in all other columns. */
        for(k = j + 1 ; k < e0.ncols ; k++) {
            /* TODO : Over the binary field, this branch avoiding trick
             * could most probably be deleted, I doubt it gains anything. */
            unsigned long c = e0.coeff(i,k);
            /* add c times column j to column k */
            e0.acol(k,j,i,c);
            // E.acol(k,j,c);
            // This one is tempting, but it's a wrong assert.
            // ASSERT(PI.deg(j) <= PI.deg(k));
            // ASSERT(delta[j] <= delta[k]);
            ASSERT(std::make_pair(delta[j],PI.deg(j)) <= std::make_pair(delta[k],PI.deg(k)));
            PI.acol(k,j,c);
        }
        PI.xmul_col(j);
        // E.xmul_col(j);
        delta[j]++;
        if (PI.deg(j) >= (int) PI.ncoef - 1) {
            overflowed.push_back(j);
            /* Then don't hesitate. Set the degree to -1 altoghether.
             * There's not much meaning remaining in this vector anyway,
             * so we'd better trash it */
        }
    }
    ASSERT_ALWAYS (rank == m);

    /* Normally our bound is set up so that we're assured to have at
     * least one generator found.
     */
    if (!overflowed.empty()) {
        std::string s = intlist_to_string(overflowed.begin(), overflowed.end());

        printf("%-8u** %zu cols (%s) exceed maxdeg=%d (normal at the end) **\n",
                t, overflowed.size(), s.c_str(), PI.ncoef - 1);
        unsigned ctot = 0;
        for (unsigned int j = 0; j < m + n; j++) {
            ctot += chance_list[j];
        }
        ASSERT_ALWAYS(ctot > 0);
    }

    return !overflowed.empty();
    /*
    std::cout << "Invertible e0 matrix\n";
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < m ; j++) {
            std::cout << e0.coeff(i,piv[j]);
        }
        std::cout << "\n";
    }
    std::cout << "[ rows";
    for(unsigned int j = 0 ; j < m ; j++) {
        std::cout << " " << piv[j];
    }
    std::cout << " ]\n";
    */

    /*
    std::cout << "Invertible e0 matrix\n";
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < m ; j++) {
            std::cout << e0.coeff(i,piv[j]);
        }
        std::cout << "\n";
    }
    std::cout << "[ rows";
    for(unsigned int j = 0 ; j < m ; j++) {
        std::cout << " " << piv[j];
    }
    std::cout << " ]\n";
    */
}/*}}}*/

#if 0/* {{{ */
const char *pi_meta_filename = "pi-%d-%d";

static int retrieve_pi_files(struct t_poly ** p_pi, int t_start)
{
    DIR * pi_dir;
    struct dirent * curr;
    int n_pi_files;
    struct couple {int s; int e;} * pi_files;
    const char *pattern;
    int i;
    struct t_poly * left = NULL, * right = NULL;
    ft_order_t order;
    int o_i;
    struct dft_bb * dft_left, * dft_right, * dft_prod;
    double tt;

    *p_pi=NULL;

    if ((pi_dir=opendir("."))==NULL) {
        perror(".");
        return t_start;
    }

    printf("Scanning directory %s for pi files\n", ".");

    pattern = strrchr(pi_meta_filename,'/');
    if (pattern == NULL) {
        pattern = pi_meta_filename;
    } else {
        pattern++;
    }

    for(n_pi_files=0;(curr=readdir(pi_dir))!=NULL;) {
        int s,e;
        if (sscanf(curr->d_name,pattern,&s,&e)==2) {
            printf("Found %s\n", curr->d_name);
            if (s>e) {
                printf("but that's a stupid one\n");
                continue;
            }
            n_pi_files++;
        }
    }

    if (n_pi_files==0) {
        printf("Found no pi files\n");
        return t_start;
    }

    pi_files=(struct couple *) malloc(n_pi_files*sizeof(struct couple));

    rewinddir(pi_dir);

    for(i=0;i<n_pi_files && (curr=readdir(pi_dir))!=NULL;) {
        if (sscanf(curr->d_name,
                    pattern,
                    &(pi_files[i].s),
                    &(pi_files[i].e))==2)
        {
            if (pi_files[i].s>pi_files[i].e)
                continue;
            i++;
        }
    }
    n_pi_files=i;
    closedir(pi_dir);

    /* The rule is: only look at the best candidate. It's not worth
     * bothering about more subtle cases */
    for(;;) {
        int t_max;
        int best=-1;
        FILE *f;
        char filename[FILENAME_LENGTH];

        printf("Scanning for data starting at t=%d\n",t_start);
        t_max=-1;
        for(i=0;i<n_pi_files;i++) {
            if (pi_files[i].s==t_start && pi_files[i].e != -1) {
                printf("candidate : ");
                printf(pattern,pi_files[i].s,pi_files[i].e);
                printf("\n");
                if (pi_files[i].e>t_max) {
                    t_max=pi_files[i].e;
                    best=i;
                }
            }
        }
        if (t_max==-1) {
            printf("Could not find such data\n");
            break;
        }

        sprintf(filename,pi_meta_filename,t_start,t_max);
        printf("trying %s\n", filename);
        f=fopen(filename,"rb");
        if (f==NULL) {
            perror(filename);
            pi_files[best].e=-1;
            continue;
        }
        /* Which degree can we expect for t_start..t_max ?
         */

        unsigned int pideg;
        pideg = iceildiv(m * (t_max - t_start), (m+n));
        pideg += 10;
        if (t_max > sequence_length) {
            pideg += t_max - sequence_length;
        }

        right=tp_read(f, pideg);
        fclose(f);

        if (right==NULL) {
            printf("%s : bad or nonexistent data\n",filename);
            pi_files[best].e=-1;
            continue;
        }

        if (left==NULL) {
            left=right;
            right=NULL;
            t_start=t_max;
            continue;
        }

        printf("Beginning multiplication\n");
        *p_pi=tp_comp_alloc(left,right);
        core_if_null(*p_pi,"*p_pi");

        order.set((*p_pi)->degree+1);
        o_i = order;

        dft_left=fft_tp_dft(left,order,&tt);
        printf("DFT(pi_left,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_left,"dft_left");

        dft_right=fft_tp_dft(right,order,&tt);
        printf("DFT(pi_right,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_right,"dft_right");

        dft_prod=fft_bbb_conv(dft_left,dft_right,&tt);
        printf("CONV(pi_left,pi_right,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_prod,"dft_prod");

        fft_tp_invdft(*p_pi,dft_prod,&tt);
        printf("IDFT(pi,%d) : %.2fs\n",o_i,tt);

        tp_free(left);
        tp_free(right);
        dft_bb_free(dft_left);
        dft_bb_free(dft_right);
        dft_bb_free(dft_prod);
        left=*p_pi;
        right=NULL;
        t_start=t_max;
    }
    free(pi_files);

    *p_pi=left;
    return t_start;
}

static void banner_traditional(int t, int deg, double inner, double * last)
{
    *last+=inner;
    if ((t+1)%print_min==0 || t==deg) {
        reclevel_prolog();
        printf("avg=%.1f	"
                "step:	%.2fs	last %d : %.2fs\n",
                ((double)global_sum_delta)/(m+n),inner,
                1+(t%print_min),*last);
        *last=0.0;
    }
}
#endif/*}}}*/

/* There are two choices for the quadratic algorithm. The first one
 * (default) seems to be a bit faster. Once I've got a more reasonable
 * implementation, it'll be time to do a serious comparison.
 *
 * The complexity curve is steeper with the first version.
 */

#if 0/*{{{*/
static void bw_traditional_algo_1(struct e_coeff * ec, int * delta,
        struct t_poly * pi, int check_chance)
{
    unsigned int * perm;
    int t;
    bw_mbmat e;
    unsigned int * pivlist;
    double inner_1,inner_2,last;
    int k;

    perm=(unsigned int *) malloc((m+n)*sizeof(unsigned int));
    pivlist=(unsigned int *) malloc(m*sizeof(unsigned int));

    mbmat_alloc(e);
    mbmat_zero(e);

    ASSERT(!ec_is_twisted(ec));
    last=0.0;

    for(t=0;t<=ec->degree;t++) {
        compute_ctaf(e,ec,pi,t,t?pivlist:NULL, &inner_1);
        column_order(perm,delta,pi);
        tp_apply_perm(pi,perm);
        if (check_chance)
            bw_check_chance(e,pi->clist);
        bw_gauss_onestep(e,ec,pi,delta,pivlist,	&inner_2);
        t_counter++;
        banner_traditional(t, ec->degree, inner_1 + inner_2, &last);
    }
    printf("DELTA : ( ");
    for(k=0;k< m + n ;k++) printf("%d ",delta[k]);
    printf(")\n");

    ec_advance(ec,ec->degree+1);	/* cosmetic */
    mbmat_free(e);
    free(pivlist);
    free(perm);
}

static void bw_traditional_algo_2(struct e_coeff * ec, int * delta,
        struct t_poly * pi, int check_chance)
{
    unsigned int * perm;
    int t;
    int deg;
    double inner,last;

    perm=(unsigned int *) malloc((m+n)*sizeof(unsigned int));

    ASSERT(!ec_is_twisted(ec));

    deg=ec->degree;
    last=0.0;

    for(t=0;t<=deg;t++) {
        if (check_chance)
            bw_check_chance(mbpoly_coeff(ec->p,0),ec->clist);
        column_order(perm,delta,pi);
        tp_apply_perm(pi,perm);
        ec_apply_perm(ec,perm);
        bw_gauss_onestep(mbpoly_coeff(ec->p,0),ec,pi,delta,NULL,&inner);
        ec_advance(ec,1);
        t_counter++;
        banner_traditional(t, deg, inner, &last);
    }

    free(perm);
}
#endif/*}}}*/

static unsigned long extract_coeff_degree_t(unsigned int t, unsigned long const * a, unsigned int da, unsigned long const * b, unsigned int db)/*{{{*/
{
    unsigned long c = 0;
    /*
    0 <= s <= t
    0 <= s <= na
    0 <= t-s <= nb
    t-nb <= s <= t
    */
    unsigned int high = std::min(t, da);
    unsigned int low = std::max(0u, t - db);
    for(unsigned int s = low ; s <= high ; s++) {
        unsigned int si = s / ULONG_BITS;
        unsigned int ss = s % ULONG_BITS;
        unsigned int ri = (t-s) / ULONG_BITS;
        unsigned int rs = (t-s) % ULONG_BITS;
        c ^= a[si] >> ss & b[ri] >> rs;
    }
    return c & 1UL;
}/*}}}*/


static void extract_coeff_degree_t(unsigned int tstart, unsigned int dt, unsigned int piv[], polmat const& PI)/*{{{*/
{
    using namespace std;
    using namespace globals;
    vector<bool> known(m+n,false);
    if (piv != NULL) {
        for(unsigned int i = 0 ; i < m ; i++)
            known[piv[i]] = true;
    }

    vector<unsigned int> z;

    for(unsigned int j = 0 ; j < m+n ; j++) {
        if (known[j]) continue;
        e0.zcol(j);
        unsigned long some_nonzero = 0;
        for(unsigned int i = 0 ; i < m ; i++) {
            unsigned long c = 0;
            for(unsigned int k = 0 ; k < m + n ; k++) {
                c ^= extract_coeff_degree_t(dt,
                        E.poly(i, k), E.ncoef + 1,
                        PI.poly(k, j), PI.deg(j));
            }
            e0.addcoeff(i,j,c);
            some_nonzero |= c;
        }
        if (!some_nonzero) {
            z.push_back(j);
        }
    }
    vector<unsigned int> ncha(m + n, 0);
    if (z.empty()) {
        chance_list.swap(ncha);
    } else {
        for(unsigned int i = 0 ; i < z.size() ; i++) {
            ++ncha[z[i]];
        }

        /* resets the global chance_list counter */
        for(unsigned int i = 0 ; i < m + n ; i++) {
            if (ncha[i] == 0) {
                chance_list[i] = 0;
            } else {
                chance_list[i] += ncha[i];
            }
        }

        printf("%-8u%zucols=0:", tstart + dt, z.size());

        vector<pair<unsigned int, unsigned int> > zz;
        for(unsigned int i = 0 ; i < z.size() ; i++) {
            zz.push_back(make_pair(z[i], chance_list[z[i]]));
        }

        // Now print this out more nicely.
        sort(zz.begin(),zz.end());
        for( ; zz.size() ; ) {
            unsigned int mi = UINT_MAX;
            for(unsigned int i = 0 ; i < zz.size() ; i++) {
                if (zz[i].second < mi) {
                    mi = zz[i].second;
                }
            }
            printf(" [");
            for(unsigned int i = 0 ; i < zz.size() ; ) {
                unsigned int j;
                for(j = i; j < zz.size() ; j++) {
                    if (zz[j].first-zz[i].first != j-i) break;
                }
                if (i) printf(",");
                if (zz[i].first == zz[j-1].first - 1) {
                    printf("%u,%u", zz[i].first, zz[j-1].first);
                } else if (zz[i].first < zz[j-1].first) {
                    printf("%u..%u", zz[i].first, zz[j-1].first);
                } else {
                    printf("%u", zz[i].first);
                }
                i = j;
            }
            printf("]");
            if (mi > 1)
                printf("*%u",mi);

            vector<pair<unsigned int, unsigned int> > zz2;
            for(unsigned int i = 0 ; i < zz.size() ; i++) {
                if (zz[i].second > mi) {
                    zz2.push_back(make_pair(zz[i].first,zz[i].second));
                }
            }
            zz.swap(zz2);
        }
        printf("\n");
    }
}/*}}}*/

static unsigned int pi_deg_bound(unsigned int d)/*{{{*/
{
    using namespace globals;
    /* How many coefficients should we allocate for the transition matrix pi
     * corresponding to a degree d increase ?
     *
     * The average value is
     *
     * \left\lceil d\frac{m}{m+n}\right\rceil
     *
     * With large probability, we expect the growth of the degree in the
     * different columns to be quite even. Only at the end of the computation
     * does it begin to go live.
     *
     * The + 11 excess here is sort of a safety net.
     */
    return iceildiv(d * m, (m + n)) + 11;
}/*}}}*/

/*
 * Rule : in the following, ec can be trashed at will.
 *
 * Compute pi_left, of degree (ec->degree + 1)*(m/(m+n)), such that ec * pi is
 * divisible by X^(ec->degree + 1) (that is, all coefficients up to
 * degree ec->degree in the product are forced to 0).
 *
 */

static bool go_quadratic(polmat& pi)/*{{{*/
{
    using namespace globals;

    unsigned int piv[m];
    unsigned int deg = E.ncoef - 1;

    polmat tmp_pi(m + n, m + n, pi_deg_bound(deg) + 1);
    for(unsigned int i = 0 ; i < m + n ; i++) {
        tmp_pi.addcoeff(i,i,0,1UL);
        tmp_pi.deg(i) = 0;
    }

    rearrange_ordering(tmp_pi, NULL);

    unsigned int tstart = t;
    bool finished = false;
    for (unsigned int dt = 0; !finished && dt <= deg ; dt++) {
#ifdef  VERBOSE
        double delta;
        delta = seconds() - start_time;
        double percent = (double) dt / (deg + 1);
        percent = percent * percent;
        double estim_final = delta / percent;
        percent *= 100.0;
        printf("%5.0f / est %-7.0f (%2.0f%%) ",
                delta, estim_final, percent);
        print_deltas();
#endif
	extract_coeff_degree_t(tstart, dt, dt ? piv : NULL, tmp_pi);
        finished = gauss(piv, tmp_pi);
        rearrange_ordering(tmp_pi, piv);
        t++;
        // if (t % 60 < 10 || deg-dt < 30)
        // write_pi(tmp_pi,tstart,t);
    }
    pi.swap(tmp_pi);
#ifdef  DO_EXPENSIVE_CHECKS
        write_pi(pi,tstart,t);
#else
    if (bw->checkpoints)
        write_pi(pi,tstart,t);
#endif

#ifdef  VERBOSE
    print_deltas();
#endif

#ifdef  DO_EXPENSIVE_CHECKS
    printf("Checking\n");
    multiply_slow(E_saved,E_saved,pi);
    unsigned int v = E_saved.valuation();
    ASSERT(v == t-tstart);
    E_saved.xdiv_resize(t-tstart, E_saved.ncoef-(t-tstart));
    for (uint j = 0; j < E_saved.ncols; j++) {
        E_saved.setdeg(j);
    }
#endif

    return finished;
}/*}}}*/

/* {{{ a tool for measuring timings and displaying progress of a
 * recursive algorithm
 */
struct recursive_tree_timer_t {
    struct spent_time {
        unsigned int step;
        double total;
        double proper;
        spent_time() : step(0), total(0), proper(0) {}
    };

    std::vector<spent_time> spent;

    std::vector<std::pair<double, double> > stack;

    void push() {
        unsigned int level = stack.size();

        if (spent.size() <= level) {
            assert(spent.size() == level);
            spent.push_back(spent_time());
        }

        ASSERT(spent.size() > level);

        double st = seconds();
        double children = 0;
        if (spent.size() > level + 1) {
            children = -spent[level+1].total;
        }

        stack.push_back(std::make_pair(st, children));
    }

    void pop(unsigned int t) {
        unsigned int level = stack.size() - 1;

        double st = stack.back().first;
        double children = stack.back().second;
        stack.pop_back();

        // unsigned int t1 = t;
        double dtime = seconds() - st;

        if (spent.size() > level + 1) {
            children += spent[level+1].total;
        }

        spent[level].step++;
        spent[level].total += dtime;

        double ptime = dtime - children;

        spent[level].proper += ptime;

        double pct_loc = spent[level].step / (double) (1 << level);

        /* make up some guess about the total time of all levels */
        unsigned int outermost = level;
        for(unsigned int back = 0 ; back <= level ; back++) {
            if (spent[level-back].step == 0)
                break;
            outermost = level-back;
        }
        unsigned int innermost = spent.size() - 1;
        double estim_tot = 0;
        double spent_tot = 0;
        for(unsigned int i = outermost ; i <= innermost ; i++) {
            estim_tot += spent[i].proper * (1 << i) / (double) spent[i].step;
            spent_tot += spent[i].proper;
        }
        double estim_above = 0;
        double spent_above = 0;
        for(unsigned int i = level ; i <= innermost ; i++) {
            estim_above += spent[i].proper * (1 << i) / (double) spent[i].step;
            spent_above += spent[i].proper;
        }

        bool leaf = level == spent.size() - 1;

        printf("%-8u", t);
        /*
        if (leaf) {
            printf(" %.2f", ptime);
        } else {
            printf(" %.2f+%.2f", ptime, children);
        }
        */
        // printf(" (total %.2f)\n", spent_tot);
        // printf("      est:");
        char * buf;
        int rc = 0;
        if (leaf) {
            rc = asprintf(&buf, "[%u]: %.1f/%.1f",
                    level, spent[level].proper, spent[level].proper / pct_loc);
            ASSERT_ALWAYS(rc >= 0);
        } else {
            rc = asprintf(&buf, "[%u,%u+]: %.1f/%.1f,%.1f",
                    level, level, spent[level].proper, spent[level].proper / pct_loc,
                    estim_above);
            ASSERT_ALWAYS(rc >= 0);
        }
        printf("%-36s", buf);
        free(buf);
        rc = asprintf(&buf, "[%u+]: %.1f/%.1f (%.0f%%)",
                outermost, spent_tot, estim_tot, 100.0 * spent_tot/estim_tot);
        ASSERT_ALWAYS(rc >= 0);
        printf("%-27s", buf);
            free(buf);
        printf("\n");
    }
    void final_info()
    {
        for(unsigned int i = 0 ; i < spent.size() ; i++) {
            printf("[%d] spent %.2f",
                    i, spent[i].proper);
            if (i < spent.size() - 1) {
                printf(" (%.2f counting children)",
                        spent[i].total);
            }
            printf("\n");
        }
        printf("Total computation took %.2f\n", spent[0].total);
    }
};      /* }}} */


static bool compute_lingen(polmat& pi, recursive_tree_timer_t&);

template<typename fft_type>
static bool go_recursive(polmat& pi, recursive_tree_timer_t& tim)
{
    using namespace globals;

    /* E is known up to O(X^E.ncoef), so we'll consider this is a problem
     * of degree E.ncoef -- this is exactly the number of increases we
     * have to make
     */
    unsigned int length_E = E.ncoef;
    unsigned int rlen = length_E / 2;
    unsigned int llen = length_E - rlen;

    ASSERT(llen && rlen && llen + rlen == length_E);

    unsigned int expected_pi_deg = pi_deg_bound(llen);
    unsigned int kill;
    unsigned int tstart = t;
    unsigned int tmiddle;
    bool finished_early;
#ifdef  DO_EXPENSIVE_CHECKS
    checker c;
#endif

#ifdef	HAS_CONVOLUTION_SPECIAL
    kill=llen;
#else
    kill=0;
#endif

    // std::cout << "Recursive call, degree " << length_E << std::endl;
    // takes lengths.
    fft_type o(length_E, expected_pi_deg + 1,
            /* length_E + expected_pi_deg - kill, */
            m + n);

    tpolmat<fft_type> E_hat;

    /* The transform() calls expect a number of coefficients, not a
     * degree. */
    transform(E_hat, E, o, length_E);

    /* ditto for this one */
    E.resize(llen);

    polmat pi_left;
    finished_early = compute_lingen(pi_left, tim);
    E.clear();
    int pi_l_deg = pi_left.maxdeg();

    if (t < t0 + llen) {
        ASSERT(finished_early);
    }
    if (pi_l_deg >= (int) expected_pi_deg) {
        printf("%-8u" "deg(pi_l) = %d >= %u ; escaping\n",
                t, pi_l_deg, expected_pi_deg);
        finished_early=1;
    }

    if (finished_early) {
        printf("%-8u" "deg(pi_l) = %d ; escaping\n",
                t, pi_l_deg);
        pi.swap(pi_left);
        return true;
    }


    tmiddle = t;

    // printf("deg(pi_l)=%d, bound is %d\n",pi_l_deg,expected_pi_deg);

    /* Since we break early now, this should no longer occur except in
     * VERY pathological cases. Two uneven increases might go unnoticed
     * if the safety margin within pi_deg_bound is loose enough. However
     * this will imply that the product has larger degree than expected
     * -- and will keep going on the steeper slope. Normally this would
     * go with ``finishing early'', unless there's something odd with
     * the rational form of the matrix -- typically  a matrix with
     * hu-u-uge kernel might trigger odd behaviour to this respect.
     *
     * For the record: (20080128):
     *
     * ./doit.pl msize=500 dimk=99 mn=8 vectoring=8 modulus=2 dens=4
     * multisols=1 tidy=0 seed=714318 dump=1
     *
     * Another case which failed (20101129) on a 64-bit Core 2:
     * cadofactor.pl params=params.c100 n=8629007704268343292699373415320999727349017259324300710654086838576679934757298577716056231336635221 bwmt=4x3
     * (increasing the safety net in pi_deg_bound from +10 to +11 made it work)
     */
    ASSERT_ALWAYS(pi_l_deg < (int) expected_pi_deg);
#if 0
    if (pi_l_deg >= (int) expected_pi_deg) {/*{{{*/
        printf("Warning : pi grows above its expected degree...\n");
        /*
        printf("order %d , while :\n"
                "deg=%d\n"
                "deg(pi_left)=%d\n"
                "llen-1=%d\n"
                "hence, %d is too big\n",
                so_i,
                deg,pi_left->degree,llen - 1,
                deg+pi_left->degree-kill + 1);
                */

        int n_exceptional = 0;
        for(unsigned int i = 0 ; i < m + n ; i++) {
            n_exceptional += chance_list[i] * m;
        }
        if (!n_exceptional) {
            die("This should only happen at the end of the computation\n",1);
        }
        pi.swap(pi_left);
        return;
    }/*}}}*/
#endif

    tpolmat<fft_type> pi_l_hat;
    /* The transform() calls expect a number of coefficients, not a
     * degree ! */
    transform(pi_l_hat, pi_left, o, pi_l_deg + 1);
    pi_left.clear();

    tpolmat<fft_type> E_middle_hat;

    /* TODO XXX Do special convolutions here */
    compose(E_middle_hat, E_hat, pi_l_hat, o);

    E_hat.clear();
    /* pi_l_hat is used later on ! */

    /* The transform() calls expect a number of coefficients, not a
     * degree ! */
    itransform(E, E_middle_hat, o, length_E + pi_l_deg - kill + 1);
    E_middle_hat.clear();

    /* Make sure that the first llen-kill coefficients of all entries of
     * E are zero. It's only a matter of verification, so this does not
     * have to be critical. Yet, it's an easy way of spotting bugs as
     * well... */

    ASSERT_ALWAYS(E.valuation() >= llen - kill);

    /* This chops off some data */
    E.xdiv_resize(llen - kill, rlen);

    polmat pi_right;
    finished_early = compute_lingen(pi_right, tim);
    int pi_r_deg = pi_right.maxdeg();
    E.clear();

    tpolmat<fft_type> pi_r_hat;
    /* The transform() calls expect a number of coefficients, not a
     * degree ! */
    transform(pi_r_hat, pi_right, o, pi_r_deg + 1);
    pi_right.clear();

    tpolmat<fft_type> pi_hat;
    compose(pi_hat, pi_l_hat, pi_r_hat, o);
    pi_l_hat.clear();
    pi_r_hat.clear();

    /* The transform() calls expect a number of coefficients, not a
     * degree ! */
    itransform(pi, pi_hat, o, pi_l_deg + pi_r_deg + 1);

    if (bw->checkpoints) {
        write_pi(pi,tstart,t);
        {
            unlink_pi(tstart,tmiddle);
            unlink_pi(tmiddle,t);
        }
    }

#ifdef  DO_EXPENSIVE_CHECKS
    c.check(pi);
#endif

    if (finished_early) {
        printf("%-8u" "deg(pi_r) = %d ; escaping\n",
                t, pi.maxdeg());
    }
    return finished_early;
}

static bool compute_lingen(polmat& pi, recursive_tree_timer_t & tim)
{
    /* reads the data in the global thing, E and delta. ;
     * compute the linear generator from this.
     */
    using namespace globals;
    unsigned int deg_E = E.ncoef - 1;

    bool b;

    // unsigned int t0 = t;
    
    tim.push();

    if (deg_E <= lingen_threshold) {
        b = go_quadratic(pi);
    } else if (deg_E < cantor_threshold) {
        /* The bound is such that deg + deg/4 is 64 words or less */
        b = go_recursive<fake_fft>(pi, tim);
    } else {
        /* Presently, c128 requires input polynomials that are large
         * enough.
         */
        b = go_recursive<c128_fft>(pi, tim);
    }

    tim.pop(t);

    return b;
}


void recycle_old_pi(polmat& pi MAYBE_UNUSED)
{
#if 0/*{{{*/
    if (pi_left!=NULL && !(!check_input && ec->degree<new_t-t_counter)) {
        int dg_kill;
        struct dft_bb * dft_pi_left;
        struct dft_mb * dft_e_left, * dft_e_middle;
        ft_order_t order;
        int o_i;

        printf("Beginning multiplication (e_left*pi_left)\n");

        /* That's a pity to bring the DFT so far, but we have to
         * check the input is correct */

        dg_kill=new_t-t_counter;

#ifdef	HAS_CONVOLUTION_SPECIAL
        if (check_input) {
            order.set(1+pi_left->degree+ec->degree);
        } else {
            order.set(1+pi_left->degree+ec->degree-dg_kill);
        }
#else
        order.set(1+pi_left->degree+ec->degree);
#endif
        o_i = order;

        dft_e_left=fft_ec_dft(ec,order,&tt);
        printf("DFT(e_left,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_e_left,"dft_e_left");

        dft_pi_left=fft_tp_dft(pi_left,order,&tt);
        printf("DFT(pi_left,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_pi_left,"dft_pi_left");

        ec_untwist(ec);
#ifdef  HAS_CONVOLUTION_SPECIAL
        if (check_input) {
            dft_e_middle=fft_mbb_conv_sp(dft_e_left,dft_pi_left,
                    0,&tt);
        } else {
            dft_e_middle=fft_mbb_conv_sp(dft_e_left,dft_pi_left,
                    dg_kill,&tt);
            ec_advance(ec,dg_kill);
        }
#else
        dft_e_middle=fft_mbb_conv(dft_e_left,dft_pi_left,&tt);
#endif

        printf("CONV(e_left,pi_left,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_e_middle,"dft_e_middle");

#ifdef  HAS_CONVOLUTION_SPECIAL
        if (check_input) {
            fft_mb_invdft(ec->p,dft_e_middle,
                    ec->degree,&tt);
            printf("IDFT(e,%d) : %.2fs\n",o_i,tt);
            printf("Verifying the product\n");
            check_zero_and_advance(ec, dg_kill);
            printf("Input OK\n");
        } else {
            fft_mb_invdft(ec->p,dft_e_middle,
                    ec->degree-dg_kill,&tt);
            printf("IDFT(e,%d) : %.2fs\n",o_i,tt);
        }
#else
        fft_mb_invdft(ec->p,dft_e_middle,
                ec->degree,&tt);
        printf("IDFT(e,%d) : %.2fs\n",o_i,tt);
        printf("Verifying the product\n");
        check_zero_and_advance(ec, dg_kill);
        printf("Input OK\n");
#endif
        tp_act_on_delta(pi_left,global_delta);

        dft_bb_free(dft_pi_left);
        dft_mb_free(dft_e_left);
        dft_mb_free(dft_e_middle);
        t_counter=new_t;

        if (check_pi)
            save_pi(pi_left,t_start,-1,t_counter);
    }
    if (pi_left!=NULL && !check_input && ec->degree<new_t-t_counter) {
        printf("We are not interested in the computation of e(X)\n");
        ec_advance(ec,new_t-t_counter);
        tp_act_on_delta(pi_left,global_delta);
        t_counter=new_t;
        if (check_pi)
            save_pi(pi_left,t_start,-1,t_counter);
    }
#endif/*}}}*/
}

extern "C" {
    void usage();
}

void usage()
{
    fprintf(stderr, "Usage: ./bw-lingen-binary --subdir <dir> -t <threshold>\n");
    exit(1);
}

int main(int argc, char *argv[])
{
    using namespace globals;

    setvbuf(stdout,NULL,_IONBF,0);
    setvbuf(stderr,NULL,_IONBF,0);

    param_list pl;
    param_list_init(pl);
    param_list_configure_alias(pl, "lingen-threshold", "lingen_threshold");
    param_list_configure_alias(pl, "cantor-threshold", "cantor_threshold");
    bw_common_init(bw, pl, &argc, &argv);
    param_list_parse_uint(pl, "lingen-threshold", &lingen_threshold);
    param_list_parse_uint(pl, "cantor-threshold", &cantor_threshold);
    if (param_list_warn_unused(pl)) usage();
    param_list_clear(pl);

    m = n = 0;

    if (bw->m == -1) { fprintf(stderr, "no m value set\n"); exit(1); } 

    if (lingen_threshold == 0) {
        fprintf(stderr, "no lingen_threshold value set\n");
        exit(1);
    }

    unsigned int n0, n1, j0, j1;

    { /* {{{ detect the input file -- there must be only one file. */
        DIR * dir = opendir(".");
        struct dirent * de;
        input_file[0]='\0';
        for( ; (de = readdir(dir)) != NULL ; ) {
            int len;
            int rc = sscanf(de->d_name, A_FILE_PATTERN "%n",
                    &n0, &n1, &j0, &j1, &len);
            /* rc is expected to be 4 or 5 depending on our reading of the
             * standard */
            if (rc < 4 || len != (int) strlen(de->d_name)) {
                continue;
            }
            if (input_file[0] != '\0') {
                fprintf(stderr, "Found two possible file names %s and %s\n",
                        input_file, de->d_name);
                exit(1);
            }
            size_t clen = std::max((size_t) len, sizeof(input_file));
            memcpy(input_file, de->d_name, clen);
        }
        closedir(dir);
    } /* }}} */

    if (bw->n == 0) {
        bw->n = n1 - n0;
    } else if (bw->n != (int) (n1 - n0)) {
        fprintf(stderr, "n value mismatch (config says %d, A file says %u)\n",
                bw->n, n1 - n0);
        exit(1);
    }

    m = bw->m;
    n = bw->n;

    if (bw->end == 0) {
        ASSERT_ALWAYS(bw->start == 0);
        bw->end = j1 - j0;
    } else if (bw->end - bw->start > (int) (j1 - j0)) {
        fprintf(stderr, "sequence file %s is too short\n", input_file);
        exit(1);
    }
        
    sequence_length = bw->end - bw->start;

/* bw_init
 *
 * fill in the structures with the data available on disk. a(X) is
 * fetched this way. f(X) is chosen in order to make [X^0]e(X)
 * nonsingular. Once this condition is satisfied, e(X) is computed. All
 * further computations are done with e(X).
 *
 * a(X) is thrown away as soon as e(X) is computed.
 *
 * At a given point of the algorithm, one can bet that most of the data
 * for e(X) is useless.
 *
 * XXX Therefore it is wise to shrink e(x) aggressively. The new
 * organisation of the data prevents e(x) from being properly swapped
 * out. XXX TODO !!!
 *
 */
    polmat A;

    printf("Reading scalar data in polynomial ``a''\n");
    read_data_for_series(A, j1 - j0);

    /* Data read stage completed. */
    /* TODO. Prepare the FFT engine for handling polynomial
     * multiplications up to ncmax coeffs */
    // unsigned int ncmax = (sequence_length<<1)+2;
    // ft_order_t::init(ncmax, ops_poly_args);

    delta.assign(m + n, (unsigned int) -1);
    chance_list.assign(m + n, 0);

    if (!recover_f0_data()) {
        // This is no longer useful
	// give_poly_rank_info(A, read_coeffs - 1);
	compute_f_init(A);
	write_f0_data();
    }
    set_t0_delta_from_F0();

    t = t0;
    printf("t0 = %d\n", t0);

    // A must be understood as A_computed + O(X^sequence_length)
    // therefore it does not make sense to compute coefficients of E at
    // degrees above sequence_length-1
    printf("Computing value of E(X)=A(X)F(X) (degree %d) [ +O(X^%d) ]\n",
            sequence_length-1, sequence_length);
    
    // multiply_slow(E, A, F0, t0, sequence_length-1);
    compute_E_from_A(A);

    printf("Throwing out a(X)\n");
    A.clear();


#if 0/*{{{*/
    if (check_pi) {
        if (retrieve_pi_files(&pi_left)) {
            recycle_old_pi(pi_left);
        }
    } else {
        printf("Not reading pi files due to --nopi option\n");
    }
#endif/*}}}*/

    using namespace globals;
    using namespace std;

    printf("E: %d coeffs, t=%u\n", E.ncoef, t);

    { bmat tmp_e0(m,m + n); e0.swap(tmp_e0); }

    for(unsigned int i = 0 ; i < m + n ; i++) {
        E.deg(i) = E.ncoef - 1;
    }

#ifndef NDEBUG
    E_saved.copy(E);
#endif

    start_time = seconds();

    // E.resize(deg + 1);
    polmat pi_left;
    recursive_tree_timer_t tim;
    compute_lingen(pi_left, tim);
    tim.final_info();

    int nresults = 0;
    for (unsigned int j = 0; j < m + n; j++) {
        nresults += chance_list[j];
    }
    if (nresults == 0) {
        fprintf(stderr, "No solution found\n");
        exit(1);
    }

    polmat F;
    compute_final_F_from_PI(F, pi_left);
    bw_commit_f(F);

#if 0/*{{{*/
    if (ec->degree>=0) {
        bw_lingen(ec,global_delta,&pi_right);
    } else {
        pi_right=pi_left;
        pi_left=NULL;
    }

    if (pi_left!=NULL) {
        struct dft_bb * dft_pi_left, * dft_pi_right, * dft_pi_prod;
        ft_order_t order;
        int o_i;

        printf("Beginning multiplication (pi_left*pi_right)\n");
        pi_prod=tp_comp_alloc(pi_left,pi_right);
        core_if_null(pi_prod,"pi_prod");

        order.set(pi_prod->degree+1);
        o_i = order;

        dft_pi_left=fft_tp_dft(pi_left,order,&tt);
        printf("DFT(pi_left,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_pi_left,"dft_pi_left");

        dft_pi_right=fft_tp_dft(pi_right,order,&tt);
        printf("DFT(pi_right,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_pi_right,"dft_pi_right");

        dft_pi_prod=fft_bbb_conv(dft_pi_left,dft_pi_right,&tt);
        printf("CONV(pi_left,pi_right,%d) : %.2fs\n",o_i,tt);
        core_if_null(dft_pi_prod,"dft_pi_prod");

        fft_tp_invdft(pi_prod,dft_pi_prod,&tt);
        printf("IDFT(pi,%d) : %.2fs\n",o_i,tt);

        tp_free(pi_left);
        tp_free(pi_right);
        dft_bb_free(dft_pi_left);
        dft_bb_free(dft_pi_right);
        dft_bb_free(dft_pi_prod);
        pi_left=NULL;
        pi_right=NULL;
        if (check_pi) {
            save_pi(pi_prod,t_start,new_t,t_counter);
        }
        /* new_t is the inner value that has to be discarded */
    } else {
        pi_prod=pi_right;
        /* Don't save here, since it's already been saved by
         * lingen (or it comes from the disk anyway).
         */
    }
#endif/*}}}*/
    // print_chance_list(sequence_length, chance_list);

    return 0;
}

/* vim: set sw=4 sta et: */
