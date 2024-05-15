#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <memory>
#include <cinttypes>
#include <cstdio>
#include <cstring>
#include <climits>
#include <cstdint>              // for uint32_t
#include <ctime>                // for time
#include <cstdlib>
#include <memory>

#include <vector>
#include <array>
#include <set>
#include <map>
#include <tuple>
#include <algorithm>
#include <utility>               // for pair, make_pair
#include <sstream>

#include <dirent.h>
#include <pthread.h>             // for pthread_mutex_lock, pthread_mutex_un...

#include <gmp.h>                 // for mpz_cmp_ui, gmp_randclear, gmp_randi...

#include "bw-common.h"
#include "cxx_mpz.hpp"
#include "matmul_top.hpp"
#include "matmul_top_comm.hpp"
#include "arith-generic.hpp"
#include "arith-cross.hpp"
#include "parallelizing_info.hpp"
#include "params.h"
#include "portability.h"
#include "select_mpi.h"
#include "verbose.h"    // verbose_enabled
#include "macros.h"
#include "mmt_vector_pair.hpp"
#include "utils_cxx.hpp"

struct sfile_info {/*{{{*/
    unsigned int s0,s1;
    unsigned int iter0;
    unsigned int iter1;
    char name[NAME_MAX];
    char name_pattern[NAME_MAX];
    static bool match(sfile_info& v, const char * name) {
        int k;
        int rc;
        unsigned int s0,s1;
        unsigned int iter0, iter1;
        rc = sscanf(name, "S.sols%u-%u.%u-%u%n", &s0, &s1, &iter0, &iter1, &k);
        if (rc < 4 || k != (int) strlen(name))
            return false;
        v.s0 = s0;
        v.s1 = s1;
        v.iter0 = iter0;
        v.iter1 = iter1;
        size_t res = strlcpy(v.name, name, NAME_MAX);
        ASSERT_ALWAYS(res < NAME_MAX);
        snprintf(v.name_pattern, sizeof(name_pattern), "S.sols%%u-%%u.%u-%u", iter0, iter1);
        return true;
    }
};/*}}}*/

static bool operator<(sfile_info const & a, sfile_info const & b)/*{{{*/
{
    if (a.iter1 < b.iter1) return true;
    if (a.iter1 > b.iter1) return false;
    if (a.iter0 < b.iter0) return true;
    if (a.iter0 > b.iter0) return false;
    if (a.s1 < b.s1) return true;
    if (a.s1 > b.s1) return false;
    if (a.s0 < b.s0) return true;
    if (a.s0 > b.s0) return false;
    return false;
}/*}}}*/

int exitcode = 0;

std::vector<sfile_info> prelude(parallelizing_info_ptr pi)/*{{{*/
{
    int leader = pi->m->jrank == 0 && pi->m->trank == 0;
    int char2 = mpz_cmp_ui(bw->p, 2) == 0;
    int splitwidth = char2 ? 64 : 1;
    std::vector<sfile_info> res;
    serialize_threads(pi->m);
    if (leader) {
        /* It's our job to collect the directory data.  */
        DIR * dir = opendir(".");
        struct dirent * de;
        for( ; (de = readdir(dir)) != NULL ; ) {
            sfile_info sf;
            if (!sfile_info::match(sf, de->d_name))
                continue;
            res.push_back(sf);
            if (bw->interval && sf.iter1 % bw->interval != 0) {
                fprintf(stderr,
                        "Warning: %s is not a checkpoint at a multiple of "
                        "the interval value %d -- this might indicate a "
                        "severe bug with leftover data, likely to corrupt "
                        "the final computation\n", sf.name, bw->interval);
            }
        }
        closedir(dir);
        sort(res.begin(), res.end());

        std::vector<sfile_info> res2;
        auto it = res.begin();
        for( ; it < res.end() ; ) {
            if (it->s0 < bw->solutions[0] || it->s0 >= bw->solutions[1]) {
                it++;
                continue;
            }
            /* try to read a complete set of files at this iteration */
            unsigned int x = bw->solutions[0];
            unsigned int z = bw->solutions[1];
            bool ok = true;
            auto current = *it;
            for(unsigned int y ; ok && x < z ; x = y, it++) {
                y = x + splitwidth;
                ok = it->s0 == x && it->s1 == y;
            }
            if (ok) {
                current.s0 = bw->solutions[0];
                current.s1 = bw->solutions[1];
                res2.push_back(current);
            }
        }
        std::swap(res, res2);

        unsigned int prev_iter = 0;
        for(size_t i = 0 ; i < res.size() ; i++) {
            sfile_info const & cur(res[i]);
            if (cur.iter0 && cur.iter0 != prev_iter) {
                fprintf(stderr, "Within the set of S files, file "
                        "%s seems to be the first "
                        "to come after iteration %u, therefore there is "
                        "a gap for the range %u..%u\n",
                        cur.name, prev_iter, prev_iter, cur.iter0);
                exit(EXIT_FAILURE);
            }
            prev_iter = cur.iter1;
        }
    }

    serialize(pi->m);
    unsigned long s = res.size();
    pi_bcast(&s, 1, BWC_PI_UNSIGNED_LONG, 0, 0, pi->m);
    if (!leader)
        res.assign(s, sfile_info());
    pi_bcast(res.data(), s * sizeof(sfile_info), BWC_PI_BYTE, 0, 0, pi->m);
    serialize(pi->m);
    return res;
}/*}}}*/

void fprint_signed(FILE * f, arith_generic * A, arith_generic::elt const & x)
{
    std::ostringstream os;
    arith_generic::elt * minus = A->alloc();
    A->neg(*minus, x);
    if (A->cmp(*minus, x) < 0) {
        os << "-";
        A->cxx_out(os, *minus);
    } else {
        A->cxx_out(os, x);
    }
    A->free(minus);
    fprintf(f, "%s", os.str().c_str());
}

void allgather(std::vector<unsigned int>& v, pi_comm_ptr wr)/*{{{*/
{
    /* want to collectively merge all vectors "v". boost mpi
     * would be great for that, really */

    /* first merge to one leader per node */

    std::vector<unsigned int> *mainv;
    mainv = &v;
    /* Yes, we're sharing a pointer, here */
    pi_thread_bcast(&mainv, sizeof(mainv), BWC_PI_BYTE, 0, wr);

    for(unsigned int j = 1 ; j < wr->ncores ; ++j) {
        if (wr->trank == j)
            mainv->insert(mainv->end(), v.begin(), v.end());
        serialize_threads(wr);
    }

    if (wr->trank == 0) {
        std::vector<unsigned int> allv;
        std::vector<int> sizes(wr->njobs, 0);
        std::vector<int> displs(wr->njobs, 0);
        sizes[wr->jrank] = v.size();
        int total = v.size();
        MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_INT, MPI_SUM,
                wr->pals);
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                sizes.data(), 1, MPI_INT, wr->pals);
        for(unsigned int i = 1 ; i < wr->njobs ; i++)
            displs[i] = displs[i-1] + sizes[i-1];
        allv.assign(total, 0);
        std::copy(v.begin(), v.end(), allv.begin() + displs[wr->jrank]);
        MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
                allv.data(), sizes.data(), displs.data(), MPI_UNSIGNED,
                wr->pals);
        std::swap(v, allv);
    }

    serialize_threads(wr);
    v = *mainv;
    serialize_threads(wr);
}/*}}}*/

void broadcast(std::vector<unsigned int>& v, pi_comm_ptr wr)/*{{{*/
{
    int total = v.size();
    pi_bcast(&total, 1, BWC_PI_INT, 0, 0, wr);
    if (wr->jrank || wr->trank) v.assign(total, 0);
    pi_bcast(v.data(), total, BWC_PI_UNSIGNED, 0, 0, wr);
    serialize(wr);
}/*}}}*/

std::vector<unsigned int> indices_of_zero_or_nonzero_values(mmt_vec & y, unsigned int maxidx, int want_nonzero)/*{{{*/
{
    arith_generic * A = y.abase;
    parallelizing_info_ptr pi = y.pi;

    std::vector<unsigned int> myz;

    if (pi->wr[y.d]->trank == 0 && pi->wr[y.d]->jrank == 0) {
        for(unsigned int i = 0 ; i < maxidx ; i++) {
            if (y.i0 <= i && i < y.i1) {
                if (!!want_nonzero == !A->is_zero(A->vec_item(y.v, i - y.i0))) {
                    myz.push_back(i);
                }
            }
        }

        /* in fact, a single gather at node 0 thread 0 would do */
        allgather(myz, pi->wr[!y.d]);
    }

    broadcast(myz, pi->m);    /* And broadcast that to everyone as well. */

    return myz;
}/*}}}*/

std::vector<unsigned int> indices_of_zero_values(mmt_vec & y, unsigned int maxidx)/*{{{*/
{
    return indices_of_zero_or_nonzero_values(y, maxidx, 0);
}/*}}}*/

std::vector<unsigned int> indices_of_nonzero_values(mmt_vec & y, unsigned int maxidx)/*{{{*/
{
    return indices_of_zero_or_nonzero_values(y, maxidx, 1);
}/*}}}*/

std::vector<unsigned int> get_possibly_wrong_columns(matmul_top_data & mmt)/*{{{*/
{
    parallelizing_info_ptr pi = mmt.pi;
    int tcan_print = bw->can_print && pi->m->trank == 0;

    std::vector<unsigned int> allz;

    if (mmt.n0[bw->dir] >= mmt.n0[!bw->dir]) return allz;

    /*
       if (tcan_print) {
       const char * name[2] = { "row", "column" };
       printf("// The original matrix has more %ss than %ss.\n// With nullspace=%s,this carries a slight danger that some parasite elements of a degree 2 nilpotent nullspace mask the real solution.\n// Trying to detect this situation.\n", name[!bw->dir], name[bw->dir], bw_dirtext[bw->dir]);
       }
        */

    /* Comments below assume that we're in the typical case
     * bw->dir==1, nrows > ncols */

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

    if (pi->m->trank == 0 && !bw->seed) {
        /* note that bw is shared between threads, thus only thread 0 should
         * test and update it here.
         * at pi->m->jrank > 0, we don't care about the seed anyway
         */
        bw->seed = time(NULL);
        MPI_Bcast(&bw->seed, 1, MPI_INT, 0, pi->m->pals);
    }
    serialize_threads(pi->m);

    gmp_randseed_ui(rstate, bw->seed);
    if (tcan_print)
        printf("// Random generator seeded with %d\n", bw->seed);

    /* Do that in the opposite direction compared to ymy */
    mmt_vector_pair zmz(mmt, !bw->dir);
    mmt_vec & z = zmz[0];
    mmt_vec & mz = zmz[zmz.size()-1];


    mmt_vec_set_random_inconsistent(z, rstate);
    mmt_vec_truncate_above_index(mmt, z, mmt.n0[bw->dir]);
    mmt_vec_apply_T(mmt, z);
    mmt_vec_twist(mmt, z);
    matmul_top_mul(mmt, zmz.vectors(), NULL);
    mmt_vec_untwist(mmt, z);
    mmt_apply_identity(mz, z);
    mmt_vec_allreduce(mz);
    mmt_vec_unapply_T(mmt, mz);
    serialize(pi->m);

    /* Check the indices of columns in the principal part. We're
     * interested in column indices which are still zero. So it's
     * really a loop until mmt.n0[bw->dir] */

    ASSERT_ALWAYS(mz.d == bw->dir);
    allz = indices_of_zero_values(mz, mmt.n0[bw->dir]);


    /* Do a second check with the *FULL* vector. Coordinates that are
     * still zero are no reason to be worried */

    mmt_vec_set_random_inconsistent(z, rstate);
    mmt_vec_truncate_above_index(mmt, z, mmt.n0[!bw->dir]);
    mmt_vec_apply_T(mmt, z);
    mmt_vec_twist(mmt, z);
    matmul_top_mul(mmt, zmz.vectors(), NULL);
    mmt_vec_untwist(mmt, z);
    mmt_apply_identity(mz, z);
    mmt_vec_allreduce(mz);
    mmt_vec_unapply_T(mmt, mz);
    serialize(pi->m);

    ASSERT_ALWAYS(mz.d == bw->dir);
    std::set<unsigned int> allz_set(allz.begin(), allz.end());
    for(auto j : indices_of_zero_values(mz, mmt.n0[bw->dir])) {
        allz_set.erase(j);
    }
    allz.assign(allz_set.begin(), allz_set.end());

    gmp_randclear(rstate);

    return allz;
}/*}}}*/

/* Take a vector, a priori sparse, and pick its coefficients at indices
 * marked by [rows], and put them, in order, in the j-th column of the
 * matrix pointed to by matrix, which has cblocks blocks of
 * A->simd_groupsize(A) entries. The column number j is thus made of elements
 * of the blocks whose index is congruent to
 * (j/groupsize)-th block mod cblocks ; A is my.abase.
 */
void compress_vector_to_sparse(arith_generic::elt * matrix, unsigned int j, unsigned int cblocks, mmt_vec const & my, std::vector<unsigned int> & rows)
{
    arith_generic * A = my.abase;

    unsigned int own_i0 = my.i0 + mmt_my_own_offset_in_items(my);
    unsigned int own_i1 = own_i0 + mmt_my_own_size_in_items(my);
    cxx_mpz v;
    unsigned int jq = j / A->simd_groupsize();
    // unsigned int jr = j % A->simd_groupsize(A);

    int char2 = mpz_cmp_ui(bw->p, 2) == 0;
    ASSERT_ALWAYS(char2 || A->simd_groupsize() == 1);

    for(unsigned int ii = 0 ; ii < rows.size() ; ii++) {
        unsigned int i = rows[ii];
        if (own_i0 <= i && i < own_i1) {
            arith_generic::elt const & src = A->vec_item(my.v, i - my.i0);
            arith_generic::elt & dst = A->vec_item(matrix, ii * cblocks + jq);
            A->set(dst, src);
        }
    }
}

struct abase_proxy {

    parallelizing_info_ptr pi;
    std::unique_ptr<arith_generic> A;
    pi_datatype_ptr A_pi;

    abase_proxy(parallelizing_info_ptr pi, int width)
        : pi(pi)
        , A(arith_generic::instance(bw->p, width))
    {
        A_pi = pi_alloc_arith_datatype(pi, A.get());
    }
    abase_proxy(abase_proxy&&) = default;
    abase_proxy& operator=(abase_proxy&&) = default;
    static abase_proxy most_natural(parallelizing_info_ptr pi) {
        return abase_proxy(pi, mpz_cmp_ui(bw->p, 2) == 0 ? 64 : 1);
    }
    std::map<arith_generic *, std::shared_ptr<arith_cross_generic>> tdict;
    arith_cross_generic * templates(arith_generic * A1) {
        auto it = tdict.find(A1);
        if (it == tdict.end())
            tdict[A1] = std::shared_ptr<arith_cross_generic>(arith_cross_generic::instance(A.get(), A1));
        return tdict[A1].get();
    }
    ~abase_proxy()
    {
        pi_free_arith_datatype(pi, A_pi);
    }
};

struct rhs /*{{{*/ {
    matmul_top_data & mmt;
    arith_generic * A;
    unsigned int nrhs;
    arith_generic::elt * rhscoeffs;
    abase_proxy natural;
    arith_generic * Av;

    rhs(rhs const&) = delete;

    rhs(matmul_top_data & mmt, const char * rhs_name, unsigned int solutions[2]) : mmt(mmt), A(mmt.abase), natural(abase_proxy::most_natural(mmt.pi)) /* {{{ */
    {
        nrhs = 0;
        rhscoeffs = NULL;
        Av = NULL;

        if (!rhs_name) return;

        parallelizing_info_ptr pi = mmt.pi;
        int tcan_print = bw->can_print && pi->m->trank == 0;
        int leader = pi->m->jrank == 0 && pi->m->trank == 0;

        /* This is just for a check -- in truth, it might be that the
         * code here works correctly for inhomogeneous characteristic 2,
         * but that would be pure chance, as it was never tested */
        int char2 = mpz_cmp_ui(bw->p, 2) == 0;
        ASSERT_ALWAYS(!char2);
        ASSERT_ALWAYS(A->simd_groupsize() == 1);

        Av = natural.A.get();

        if (leader)
            get_rhs_file_header(rhs_name, NULL, &nrhs, NULL);
        pi_bcast(&nrhs, 1, BWC_PI_UNSIGNED, 0, 0, pi->m);

        if (tcan_print) {
            printf("** Informational note about GF(p) inhomogeneous system:\n");
            printf("   Original matrix dimensions: %" PRIu32" %" PRIu32"\n", mmt.n0[0], mmt.n0[1]);
            printf("   We expect to obtain a vector of size %" PRIu32"\n",
                    mmt.n0[!bw->dir] + nrhs);
            printf("   which we hope will be a kernel vector for (M_square||RHS).\n"
                    "   We will discard the coefficients which correspond to padding columns\n"
                    "   This entails keeping coordinates in the intervals\n"
                    "   [0..%" PRIu32"[ and [%" PRIu32"..%" PRIu32"[\n"
                    "   in the result.\n"
                    "** end note.\n",
                    mmt.n0[bw->dir], mmt.n0[!bw->dir], mmt.n0[!bw->dir] + nrhs);
        }

        /* everyone allocates something */
        rhscoeffs = A->alloc(nrhs, ALIGNMENT_ON_ALL_BWC_VECTORS);

        if (leader) {
            // yeah, we asserted that we're GF(p) at this point anyway.
            // coverity[dead_error_line]
            unsigned int splitwidth = char2 ? 64 : 1;
            ASSERT_ALWAYS(Av->simd_groupsize() == splitwidth);

            if (char2 || solutions[1] != solutions[0] + splitwidth) {
                ASSERT_ALWAYS(0);/* never tested. I did attempt to code it right for the simd case though, but did not test. */
            }
            unsigned int Av_multiplex = (solutions[1] - solutions[0]) / splitwidth;
            for(unsigned int j = 0 ; j < nrhs ; j++) {
                for(unsigned int i = 0 ; i < Av_multiplex ; i++) {
                    char * tmp;
                    unsigned int s0 = solutions[0] + i * splitwidth;
                    unsigned int s1 = solutions[0] + (i + 1) * splitwidth;
                    int rc = asprintf(&tmp, "F.sols%u-%u.%u-%u.rhs", s0, s1, j, j+splitwidth);
                    if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_LOADING_MKSOL_FILES)) {
                        printf("loading %s\n", tmp);
                    }
                    ASSERT_ALWAYS(rc >= 0);
                    FILE * f = fopen(tmp, "rb");
                    ASSERT_ALWAYS(f);
                    rc = fread(
                            Av->vec_subvec(rhscoeffs, j * Av_multiplex + i),
                            Av->elt_stride(),
                            1, f);
                    ASSERT_ALWAYS(rc == 1);
                    if (Av->is_zero(Av->vec_item(rhscoeffs, j * Av_multiplex + i))) {
                        printf("Notice: coefficient for vector V%u-%u in file %s is zero\n", j, j+1, tmp);
                    }
                    fclose(f);
                    free(tmp);
                }
            }
        }
        pi_bcast(rhscoeffs, nrhs, mmt.pitype, 0, 0, pi->m);
    }/*}}}*/

    ~rhs() {/*{{{*/
        if (rhscoeffs)
            A->free(rhscoeffs);
    }/*}}}*/
    void fwrite_rhs_coeffs(FILE * f, unsigned int i=0) /* {{{ */
    {
        unsigned int Av_multiplex = A->simd_groupsize() / Av->simd_groupsize();
        for(unsigned int j = 0 ; j < nrhs ; j++) {
            int rc = fwrite(
                    Av->vec_subvec(rhscoeffs, j * Av_multiplex + i),
                    Av->elt_stride(),
                    1, f);
            ASSERT_ALWAYS(rc == 1);
        }
    }/*}}}*/
    void fprint_rhs_coeffs(FILE * f2) /* {{{ */
    {
        for(uint32_t i = 0 ; i < nrhs ; i++) {
            std::ostringstream os;
            A->cxx_out(os, A->vec_item(rhscoeffs, i));
            fprintf(f2, "%s\n", os.str().c_str());
        }
    }/*}}}*/
    operator bool() const { return nrhs; }

    void add_contribution(mmt_vec & y) const/*{{{*/
    {
        if (!nrhs) return;

        parallelizing_info_ptr pi = mmt.pi;
        arith_generic * A = mmt.abase;
        ASSERT_ALWAYS(y.abase == A);
        unsigned int unpadded = MAX(mmt.n0[0], mmt.n0[1]);
        size_t eblock = mmt_my_own_size_in_items(y);

        abase_proxy natural = abase_proxy::most_natural(pi);
        arith_generic * Av = natural.A.get();
        pi_datatype_ptr Av_pi = natural.A_pi;

        mmt_vec vi(mmt,Av,Av_pi, bw->dir, /* shared ! */ 1, mmt.n[bw->dir]);

        for(unsigned int j = 0 ; j < nrhs ; j++) {
            int ok = mmt_vec_load(vi, "V%u-%u.0", unpadded, j);
            ASSERT_ALWAYS(ok);

            natural.templates(A)->addmul_tiny(
                    mmt_my_own_subvec(y),
                    mmt_my_own_subvec(vi),
                    A->vec_subvec(rhscoeffs, j),
                    eblock);
        }

        /* addmul_tiny degrades consistency ! */
        y.consistency = 1;
        mmt_vec_broadcast(y);
    }/*}}}*/
};
/*}}}*/

std::tuple<int, int> check_zero_and_padding(mmt_vec & y, unsigned int maxidx)/*{{{*/
{

    /* Here, we want to make sure that we have something non-zero in
     * the **input coordinate space**. It matters little to us if we
     * found a solution of (M||zero-pad) * x = 0 with x having
     * non-zero coordinates only in the zero-pad part.
     *
     * Therefore, the following check is unfortunately not good
     * enough:

     int is_zero = A->vec_is_zero(A,
     mmt_my_own_subvec(y), mmt_my_own_size_in_items(y));

     * instead, we want to check only up to index mmt.n0[bw->dir]
     * (and bw->dir is y.d). (this is valid because at this point, y
     * is untwisted and has T unapplied).
     */
    size_t my_input_coordinates;
    size_t my_pad_coordinates;
    ASSERT_ALWAYS(y.d == bw->dir);
    size_t my_i0 = y.i0 + mmt_my_own_offset_in_items(y);
    size_t my_i1 = my_i0 + mmt_my_own_size_in_items(y);

    if (my_i0 >= maxidx) {
        my_input_coordinates = 0;
        my_pad_coordinates =  mmt_my_own_size_in_items(y);
    } else if (my_i1 >= maxidx) {
        my_input_coordinates = maxidx - my_i0;
        my_pad_coordinates =  my_i1 - maxidx;
    } else {
        my_input_coordinates = mmt_my_own_size_in_items(y);
        my_pad_coordinates = 0;
    }
    int input_is_zero = y.abase->vec_is_zero(
            mmt_my_own_subvec(y),
            my_input_coordinates);
    int pad_is_zero = y.abase->vec_is_zero(
            y.abase->vec_subvec(
                mmt_my_own_subvec(y), my_input_coordinates),
            my_pad_coordinates);

    pi_allreduce(NULL, &input_is_zero, 1, BWC_PI_INT, BWC_PI_MIN, y.pi->m);
    pi_allreduce(NULL, &pad_is_zero, 1, BWC_PI_INT, BWC_PI_MIN, y.pi->m);

    return std::make_tuple(input_is_zero, pad_is_zero);
}/*}}}*/

std::tuple<int, int, int> test_one_vector(matmul_top_data & mmt, mmt_vector_pair & ymy, rhs const & R)
{
    arith_generic * A = mmt.abase;
    parallelizing_info_ptr pi = mmt.pi;

    mmt_vec & y = ymy[0];

    int input_is_zero;
    int pad_is_zero;
    int hamming_out = -1;
    std::tie(input_is_zero, pad_is_zero) = check_zero_and_padding(y, mmt.n0[bw->dir]);
    if (!input_is_zero) {
        serialize(pi->m);
        mmt_vec_apply_T(mmt, y);
        serialize(pi->m);
        mmt_vec_twist(mmt, y);
        matmul_top_mul(mmt, ymy.vectors(), NULL);
        mmt_vec_untwist(mmt, y);
        serialize(pi->m);
        /* Add the contributions from the right-hand side vectors, to see
         * whether that makes the sum equal to zero */
        if (R) {
            R.add_contribution(y);
        } else {
            // unapply_T is only valid with respect to something we'll
            // multiply *again*. 
            mmt_vec_unapply_T(mmt, y);
        }

        /* This "is zero" check is also valid on the padded matrix of
         * course, so we don't heave the same headache as above */
        int is_zero = A->vec_is_zero(
                mmt_my_own_subvec(y), mmt_my_own_size_in_items(y));
        pi_allreduce(NULL, &is_zero, 1, BWC_PI_INT, BWC_PI_MIN, pi->m);

        hamming_out = is_zero ? 0 : mmt_vec_hamming_weight(y);
    }
    return std::make_tuple(input_is_zero, pad_is_zero, hamming_out);
}

/* Use y_saved as input (and leave it untouched). Store result in both y
 * and my */
std::tuple<int, int, int> expanded_test(matmul_top_data & mmt, mmt_vector_pair & ymy, mmt_vec const & y_saved, rhs const& R)
{
    parallelizing_info_ptr pi = mmt.pi;
    mmt_vec & y = ymy[0];
    mmt_vec & my = ymy[ymy.size()-1];
    mmt_full_vec_set(y, y_saved);
    auto res = test_one_vector(mmt, ymy, R);

    /* Need to get the indices with respect to !bw->dir...  */
    serialize(pi->m);
    mmt_apply_identity(my, y);
    mmt_vec_allreduce(my);
    mmt_vec_unapply_T(mmt, my);
    serialize(pi->m);
    return res;
}

/* This is for a specific class of parasites that may be caused by empty
 * columns in the matrix.
 */
class parasite_fixer {/*{{{*/
    matmul_top_data & mmt;
    arith_generic * A;
    parallelizing_info_ptr pi;

    // typedef std::map<std::pair<unsigned int, unsigned int>, cxx_mpz> pre_matrix_t;

    arith_generic::elt * matrix;
    public:
    bool attempt_to_fix;

    std::vector<unsigned int> cols;
    std::vector<unsigned int> rows;
    size_t nrows() const { return rows.size(); }
    std::vector<std::pair<std::array<unsigned int, 2>, int> > pivot_list;

    parasite_fixer(matmul_top_data & mmt) : mmt(mmt), A(mmt.abase), pi(mmt.pi) {/*{{{*/
        matrix = NULL;
        int tcan_print = bw->can_print && pi->m->trank == 0;

        cols = get_possibly_wrong_columns(mmt);

        attempt_to_fix = !cols.empty() && cols.size() < 64;

        if (!cols.empty() && tcan_print) {
            printf("# %zu possibly wrong coordinates detected"
                    " in solution vector because the matrix"
                    " has a non-trivial nilpotent space.\n",
                    cols.size());
            if (attempt_to_fix) {
                printf("# Will try to fix\n");
            } else {
                printf("# Deliberately avoiding attempt to fix because"
                        " this number of vectors is large."
                        " It is rather an indication that something is wrong."
                        " Proceeding anyway\n");
            }
        }
        compute_pivot_list();
    }/*}}}*/

    parasite_fixer(parasite_fixer const&) = delete;

    /*{{{ row_coordinates_of_nonzero_cols(matmul_top_data & mmt) */
    std::vector<unsigned int> row_coordinates_of_nonzero_cols(matmul_top_data & mmt, std::vector<unsigned int> const& cols)
    {
        // int tcan_print = bw->can_print && pi->m->trank == 0;

        arith_generic * A = mmt.abase;

        mmt_vector_pair ymy(mmt, bw->dir);
        mmt_vec & y = ymy[0];
        mmt_vec & my = ymy[ymy.size()-1];

        /* Now try to see which indices are potentially affected */
        unsigned int B = A->simd_groupsize();
        for(unsigned int jjq = 0 ; jjq < cols.size() ; jjq+=B) {
            mmt_full_vec_set_zero(y);
            for(unsigned int jjr = 0 ; jjr < B && (jjq + jjr < cols.size()) ; jjr++) {
                unsigned int j = cols[jjq + jjr];
                mmt_vec_add_basis_vector_at(y, jjr, j);
            }
            mmt_vec_apply_T(mmt, y);
            mmt_vec_twist(mmt, y);
            matmul_top_mul(mmt, ymy.vectors(), NULL);
            mmt_vec_untwist(mmt, y);
            /* Not entirely clear to me if I should unapply_T here or not
            */
            mmt_apply_identity(my, y);
            mmt_vec_allreduce(my);
            mmt_vec_unapply_T(mmt, my);
            serialize(pi->m);

            /* This reads the global list */
            std::vector<unsigned int> kk = indices_of_nonzero_values(my, mmt.n0[!bw->dir]);

#if 0
            if (tcan_print) {
                printf("# List of the %zu non-zero coordinates attached to input coordinate(s)", kk.size());
                for(unsigned int jjr = 0 ; jjr < B && (jjq + jjr < cols.size()) ; jjr++) {
                    unsigned int j = cols[jjq + jjr];
                    printf(" %u", j);
                }
                printf("\n");
                for(auto k : kk)
                    printf("#\t%u\n", k);
            }
#endif
            rows.insert(rows.end(), kk.begin(), kk.end());
            serialize(pi->m);
        }
        std::sort(rows.begin(), rows.end());
        auto it = std::unique(rows.begin(), rows.end());
        rows.resize(it - rows.begin());

        return rows;
    }/*}}}*/

#if 0
    /* this one could live outside the class */
    void pre_matrix_to_matrix(pre_matrix_t & pre_matrix)/*{{{*/
    {
        if (pi->wr[!bw->dir]->jrank == 0 && pi->wr[!bw->dir]->trank == 0) {
            unsigned int kk = 0;
            for(unsigned int ii = 0 ; ii < rows.size() ; ii++) {
                unsigned int i = rows[ii];
                for(unsigned int jj = 0 ; jj < cols.size() ; jj++) {
                    unsigned int j = cols[jj];
                    auto it = pre_matrix.find(std::make_pair(i,j));
                    if (it != pre_matrix.end()) {
                        A->set_mpz(A, A->vec_coeff_ptr(A, matrix, kk), (mpz_srcptr) it->second);
                    }
                    kk++;
                }
            }
        }

        pi_allreduce(NULL, matrix,
                rows.size() * cols.size(),
                mmt.pitype, BWC_PI_SUM, pi->m);
    }/*}}}*/
#endif

    void debug_print_local_matrix(arith_generic::elt * matrix,
            std::vector<unsigned int> const & rows,
            std::vector<unsigned int> const & cols,
            arith_generic::elt * nz = NULL)/*{{{*/
    {
        size_t nr = rows.size();
        size_t nc = cols.size();
        printf("# Dump of the full matrix%s as seen by J%uT%u\n",
                nz ? " (with coefficients of the vector encountered)" : "",
                pi->m->jrank, pi->m->trank);
        unsigned int kk = 0;
        unsigned int B = A->simd_groupsize();
        unsigned int cblocks = iceildiv(nc, B);
        printf("#\t\t");
        for(unsigned int jj = 0 ; jj < nc ; jj++) {
            printf("[%u] ", cols[jj]);
        }
        printf("\n");
        for(unsigned int ii = 0 ; ii < nr ; ii++) {
            printf("#\t%u\t", rows[ii]);
            for(unsigned int jj = 0 ; jj < cblocks ; jj++) {
                printf(" ");
                fprint_signed(stdout, A, A->vec_item(matrix, kk));
                kk++;
            }
            if (nz) {
                printf(" ");
                fprint_signed(stdout, A, A->vec_item(nz, ii));
            }
            printf("\n");
        }
    }/*}}}*/

#if 0
    void debug_print_all_local_matrices(arith_generic::elt * nz = NULL)/*{{{*/
    {
        for(unsigned int jr = 0 ; jr < pi->m->njobs ; jr++) {
            for(unsigned int tr = 0 ; tr < pi->m->ncores ; tr++) {
                if (pi->m->jrank == jr && pi->m->trank == tr)
                    debug_print_local_matrix(rows.size(), cols.size(), nz);
                serialize(pi->m);
            }
        }
    }/*}}}*/
#endif

    void compute_pivot_list() {/*{{{*/
        if (!attempt_to_fix) return;

        int tcan_print = bw->can_print && pi->m->trank == 0;
        int leader = pi->m->jrank == 0 && pi->m->trank == 0;

        rows = row_coordinates_of_nonzero_cols(mmt, cols);

        arith_generic * A = mmt.abase;
        arith_generic::elt * dummy = A->alloc();

        int char2 = mpz_cmp_ui(bw->p, 2) == 0;

        /* code is similar to row_coordinates_of_nonzero_cols() */
        mmt_vector_pair ymy(mmt, bw->dir);
        mmt_vec & y = ymy[0];
        mmt_vec & my = ymy[ymy.size()-1];

        /* 1, -1: coeff is 1 or -1.
         * 2: coeff is something else, and lookup is needed (char!=2
         * only)
         */
        std::map<unsigned int, std::pair<unsigned int, int> > pivots;
        std::set<unsigned int> scols(cols.begin(), cols.end());
        std::set<unsigned int> srows;
        for(unsigned int ii = 0 ; ii < rows.size() ; ii++)
            srows.insert(ii);

        unsigned int B = A->simd_groupsize();

        matrix = A->alloc(iceildiv(cols.size(), B) * rows.size(), ALIGNMENT_ON_ALL_BWC_VECTORS);
        A->vec_set_zero(matrix, iceildiv(cols.size(), B) * rows.size());

        for(unsigned int drop = UINT_MAX, spin=0 ; drop && !scols.empty() ; spin++) {
            drop = 0;

            std::vector<unsigned int> vcols(scols.begin(), scols.end());
            std::vector<unsigned int> vrows;
            for(auto x : srows) vrows.push_back(rows[x]);

            if (tcan_print)
                printf("# Pass %d: checking for error fixing pivots within a matrix of dimension %zu*%zu\n", spin, vrows.size(), vcols.size());

            arith_generic::elt * mat;
            arith_generic::elt ** pmat = &mat;
            unsigned int cblocks = iceildiv(vcols.size(), B);
            if (spin) {
                mat = A->alloc(cblocks * vrows.size(), ALIGNMENT_ON_ALL_BWC_VECTORS);
                A->vec_set_zero(*pmat, cblocks * vrows.size());
            } else {
                pmat = &matrix;
            }

            /* Now try to see which indices are potentially affected by these
             * columns. */
            for(unsigned int jjq = 0 ; jjq < cols.size() ; jjq += B) {
                mmt_full_vec_set_zero(y);
                for(unsigned int jjr = 0 ; jjr < B && (jjq + jjr < vcols.size()) ; jjr++) {
                    unsigned int j = vcols[jjq + jjr];
                    mmt_vec_add_basis_vector_at(y, jjr, j);
                }
                mmt_vec_apply_T(mmt, y);
                mmt_vec_twist(mmt, y);
                matmul_top_mul(mmt, ymy.vectors(), NULL);
                mmt_vec_untwist(mmt, y);
                /* Not entirely clear to me if I should unapply_T here or not
                */
                mmt_apply_identity(my, y);
                mmt_vec_allreduce(my);
                mmt_vec_unapply_T(mmt, my);
                serialize(pi->m);

                compress_vector_to_sparse(*pmat, jjq, cblocks, my, vrows);

                serialize(pi->m);
            }

            pi_allreduce(NULL, *pmat,
                    cblocks * vrows.size(),
                    mmt.pitype, BWC_PI_SUM, pi->m);

            if (leader) {
                printf("Print matrix of size %zu*%zu\n", srows.size(), scols.size());
                debug_print_local_matrix(matrix, vrows, vcols);
            }

            size_t ii = 0;
            for(auto xi : srows) {
                /* xi is the index of the row within the set of error
                 * rows.
                 * ii is the index within the set of the error rows that
                 * are being considered within this pass.
                 */
                const arith_generic::elt * row = A->vec_subvec(*pmat, ii * cblocks);
                int w = A->vec_simd_hamming_weight(row, cblocks);
                if (w == 1) {
                    int p = A->vec_simd_find_first_set(*dummy, row, cblocks);
                    if (char2) {
                        pivots[xi] = std::make_pair(vcols[p], 1);
                    } else {
                        /* Is the first non-zero value a potential pivot
                         * (+1 or -1 ?). And how do we do that with an
                         * interface that looks reasonable with
                         * characteristic two as well???
                         */
                        if (A->cmp(*dummy, 1) == 0) {
                            pivots[xi] = std::make_pair(vcols[p], 1);
                        } else if (A->neg(*dummy, *dummy), A->cmp(*dummy, 1) == 0) {
                            pivots[xi] = std::make_pair(vcols[p], -1);
                        } else if (pivots[xi].second == 0) {
                            pivots[xi] = std::make_pair(vcols[p], 2);
                        }
                    }
                }
                ii++;
            }

            for(auto const & pp : pivots) {
                unsigned int xi = pp.first;
                unsigned int j = pp.second.first;
                int v = pp.second.second;
                if (v == 2) continue;
                if (tcan_print)
                    printf("Found pivot for row %u:"
                            " column %u has coefficient %d\n",
                            rows[xi], j, v);
                srows.erase(xi);
                if (scols.erase(j)) {
                    drop++;
                    /* XXX colum j may have already been deleted */
                    std::array<unsigned int, 2> xij {{ xi, j }};
                    pivot_list.push_back({xij, v});
                }
            }
            if (spin)
                A->free(mat);

            if (tcan_print)
                printf("# Pass %d: number of cols has dropped by %u. We have %zu rows (at most) and %zu columns left\n", spin, drop, srows.size(), scols.size());
            /* use this marker to indicate synchronization */
            std::array<unsigned int, 2> xij {{ 0u, 0u }};
            pivot_list.push_back({xij, 0});
        }
        if (!scols.empty()) {
            /* we exited with drop==0 then, so srows correctly reflects
             * the last set of active rows */
            for(auto r : srows) {
                auto it = pivots.find(r);
                if (it != pivots.end()) {
                    printf("For row %u, only pivot found is non-trivial. Please implement that code\n", rows[r]);
                } else {
                    printf("For row %u, linear algebra reduction is needed. Please implement real code\n", rows[r]);
                }
            }
        }
        ASSERT_ALWAYS(scols.empty());

        A->free(dummy);

        serialize(pi->m);
    }/*}}}*/

    ~parasite_fixer() {/*{{{*/
        if (matrix)
            A->free(matrix);
    }/*}}}*/

    std::tuple<int, int, int> attempt(matmul_top_data & mmt, mmt_vector_pair & ymy, mmt_vec & y_saved, rhs const& R)/*{{{*/
    {
        mmt_vec & my = ymy[ymy.size()-1];
        int tcan_print = bw->can_print && pi->m->trank == 0;
        int leader = pi->m->jrank == 0 && pi->m->trank == 0;

        int input_is_zero;
        int pad_is_zero;
        int hamming_out;
        std::tuple<int, int, int> res;

        res = expanded_test(mmt, ymy, y_saved, R);
        std::tie(input_is_zero, pad_is_zero, hamming_out) = res;

        if (input_is_zero || hamming_out == 0) return res;

        std::vector<unsigned int> nz_pos;
        nz_pos = indices_of_nonzero_values(my, mmt.n0[!bw->dir]);

        if (tcan_print) {
            printf("# Input vector has %s input, %s padding\n",
                    input_is_zero ? "zero" : "non-zero",
                    pad_is_zero ? "zero" : "non-zero");
            printf("# Output has Hamming weight %d\n", hamming_out);
        }
        if (!std::includes(rows.begin(), rows.end(), nz_pos.begin(), nz_pos.end())) {
            if (tcan_print)
                printf("# Note: cannot attempting to fix %zu wrong coordinates, not included in the set of %zu known possibly wrong ones\n", nz_pos.size(), rows.size());
            return res;
        }

        if (tcan_print)
            printf("# Note: all the non-zero coordinates are included in the output of the \"possibly wrong\" columns, which is a good sign\n");

        arith_generic::elt * nz;


        ASSERT_ALWAYS(my.abase == mmt.abase);

        nz = A->alloc(rows.size(), ALIGNMENT_ON_ALL_BWC_VECTORS);
        A->vec_set_zero(nz, rows.size());
        compress_vector_to_sparse(nz, 0, 1, my, rows);
        pi_allreduce(NULL, nz, rows.size(), mmt.pitype, BWC_PI_SUM, pi->m); 

        if (leader) debug_print_local_matrix(matrix, rows, cols, nz);

        // debug_print_all_local_matrices(nz);

        for(auto pp : pivot_list) {
            unsigned int ii = pp.first[0];
            // unsigned int i = rows[ii];
            unsigned int j = pp.first[1];
            int v = pp.second;
            if (!v) {
                serialize(pi->m);
                res = expanded_test(mmt, ymy, y_saved, R);
                std::tie(input_is_zero, pad_is_zero, hamming_out) = res;
                A->vec_set_zero(nz, rows.size());
                compress_vector_to_sparse(nz, 0, 1, my, rows);
                pi_allreduce(NULL, nz, rows.size(), mmt.pitype, BWC_PI_SUM, pi->m); 

                if (leader) debug_print_local_matrix(matrix, rows, cols, nz);
                continue;
            }

            /* Everyone has this coefficient */
            arith_generic::elt const & error = A->vec_item(nz, ii);

            /* y_saved is shared across threads in the wiring direction,
             * so we only need to touch our very own data in there.
             *
             * (is this is ever meant to change, see mmt_full_vec_set for
             * instance)
             */
            ASSERT_ALWAYS(mmt_vec_is_shared(y_saved));

            size_t own_i0 = y_saved.i0 + mmt_my_own_offset_in_items(y_saved);
            size_t own_i1 = own_i0 + mmt_my_own_size_in_items(y_saved);

            if (own_i0 <= j && j < own_i1) {
                arith_generic::elt & source = A->vec_item(y_saved.v, j - y_saved.i0);
                printf("Row %u, coefficient is ", rows[ii]);
                fprint_signed(stdout, A, source);
                if (v == -1) {
                    printf(" ; fixing by adding to coordinate %u\n", j);
                    A->add_and_reduce(source, error);
                } else if (v == 1) {
                    printf(" ; fixing by subtracting from coordinate %u\n", j);
                    A->sub_and_reduce(source, error);
                } else {
                    ASSERT_ALWAYS(0);
                }
            }
            /* On the other hand, we're not shared across MPI nodes here
             * anyway, so we have some work to do ! */
            y_saved.consistency = 1;
            mmt_vec_broadcast(y_saved);
        }
        serialize(pi->m);

        if (tcan_print)
            printf("# After fix, Hamming weight is %d\n", hamming_out);

        // if (leader) debug_print_local_matrix(matrix, rows, cols, nz);

        A->free(nz);

        return res;
    }/*}}}*/
};/*}}}*/

void * gather_prog(parallelizing_info_ptr pi, param_list pl, void * arg MAYBE_UNUSED)
{
    ASSERT_ALWAYS(!pi->interleaved);

    int tcan_print = bw->can_print && pi->m->trank == 0;
    int leader = pi->m->jrank == 0 && pi->m->trank == 0;

    unsigned int solutions[2] = { bw->solutions[0], bw->solutions[1], };
    int char2 = mpz_cmp_ui(bw->p, 2) == 0;

    /* Define and initialize our arithmetic back-ends. More or less the
     * same deal as for mksol (for the "solutions" part). See comments
     * there. */

    /* {{{ main arithmetic backend: for the solutions we compute. */
    // unsigned int A_multiplex MAYBE_UNUSED = 1;
    unsigned int A_width = solutions[1]-solutions[0];
    if ((char2 && (A_width != 64 && A_width != 128 && A_width != 256))
            || (!char2 && A_width > 1))
    {
        fprintf(stderr,
                "We cannot support computing %u solutions at a time "
                "with one single Spmv operation, given the currently "
                "implemented code\n",
                A_width);
        exit(EXIT_FAILURE);
    }

    abase_proxy abase_solutions(pi, A_width);
    arith_generic * A = abase_solutions.A.get();

    /* }}} */

    matmul_top_data mmt(A, pi, pl, bw->dir);

    parasite_fixer pfixer(mmt);

    mmt_vector_pair ymy(mmt, bw->dir);
    mmt_vec & y = ymy[0];

    mmt_vec y_saved(mmt,0,0, bw->dir, mmt_vec_is_shared(y), mmt.n[bw->dir]);

    /* this is really a misnomer, because in the typical case, M is
     * rectangular, and then the square matrix does induce some padding.
     * This is however not the padding we are interested in. The padding
     * we're referring to in the naming of this variable is the one which
     * is related to the number of jobs and threads: internal dimensions
     * are arranged to be multiples */
    unsigned int unpadded = MAX(mmt.n0[0], mmt.n0[1]);

    std::vector<sfile_info> sl = prelude(pi);
    if (sl.empty()) {
        if (tcan_print) {
            fprintf(stderr, "Found zero S files for solution range %u..%u. "
                    "Problem with command line ?\n",
                    solutions[0], solutions[1]);
            pthread_mutex_lock(pi->m->th->m);
            exitcode=1;
            pthread_mutex_unlock(pi->m->th->m);
        }
        serialize_threads(pi->m);
        return NULL;
    }

    rhs R(mmt, param_list_lookup_string(pl, "rhs"), solutions);

    if (tcan_print)
        printf("Trying to build solutions %u..%u\n", solutions[0], solutions[1]);


    { /* {{{ Collect now the sum of the LHS contributions */
        mmt_vec svec(mmt,0,0, bw->dir, /* shared ! */ 1, mmt.n[bw->dir]);
        for(size_t i = 0 ; i < sl.size() ; i++) {
            if (tcan_print && verbose_enabled(CADO_VERBOSE_PRINT_BWC_LOADING_MKSOL_FILES)) {
                printf("loading %s\n", sl[i].name);
            }
            int ok = mmt_vec_load(svec, sl[i].name_pattern, unpadded, solutions[0]);
            ASSERT_ALWAYS(ok);

            A->vec_add_and_reduce(
                    mmt_my_own_subvec(y),
                    mmt_my_own_subvec(svec),
                    mmt_my_own_size_in_items(y));
        }
        y.consistency = 1;
        mmt_vec_broadcast(y);
    } /* }}} */

    printf("Hamming weight of sum: %lu\n", mmt_vec_hamming_weight(y));

    /* Note that for the inhomogeneous case, we'll do the loop only
     * once, since we end with a break. */
    int winning_iter = 0;

    /* This is used as an out-of-loop check for an exceptional condition.
     * This is slighly clumsy, I should do it better */
    int is_zero = 0; /* placate old compilers */

    /* The mmt_vec_unapply_T below is disturbing. Does that mean that the
     * vectors we feed are with T applied ? All of this is slightly
     * hairy...
     * 
     * I _think_ that this may be related to "apply_T" and "unapply_T"
     * being named in a way that reflects a direction which is not what
     * we have in mind all the time. "apply" for column vectors means
     * v<-T^-1*v !!
     *
     * We have, in order
     *
     * a vector read from disk
     * unapply_T
     *
     * apply_T
     * twist
     * matmul_top_mul
     * untwist
     * unapply_T
     *
     * which, for bw->dir==1, handles the following data (for a column
     * vector V on the left, for the same vector written as a row vector
     * on the right)
     *
     *  load :  T^-1 * V
     *  unapply : V
     *
     *  apply   : T^-1 * V
     *  twist   : Sc * T^-1 * V
     *  mul     : Sc*Mpad*T*Sc^-1 * T^-1 * V = Sc*Mpad*V
     *  untwist : Mpad*V
     *  unapply : T*Mpad*V
     *
     *  apply   : Mpad*V
     *  twist   : Sc*Mpad*V
     *  mul     : Sc*(Mpad*T)^2*T^-1*V
     *  untwist : (Mpad*T)^2*T^-1*V
     *  unapply : T*(Mpad*T)^2*T^-1*V
     *
     *  apply   : (Mpad*T)^2*T^-1*V
     *  (and so on)
     */
    mmt_vec_unapply_T(mmt, y);

    for(int i = 1 ; R ? (i<=1) : (i < 10) ; i++) {

        int input_is_zero;
        int pad_is_zero;
        int hamming_out;

        mmt_full_vec_set(y_saved, y);
        if (pfixer.attempt_to_fix) {
            /* uses y_saved as input */
            auto res = pfixer.attempt(mmt, ymy, y_saved, R);
            std::tie(input_is_zero, pad_is_zero, hamming_out) = res;
        } else {
            /* uses y as input */
            auto res = test_one_vector(mmt, ymy, R);
            std::tie(input_is_zero, pad_is_zero, hamming_out) = res;
        }
        is_zero = hamming_out == 0;

        /* {{{ save file. If input is zero, bail out */
        {
            char * tmp;
            int rc = asprintf(&tmp, "%sK.sols%%u-%%u.%d", input_is_zero ? "zero" : "",
                    i-1);
            ASSERT_ALWAYS(rc >= 0);
            mmt_vec_save(y_saved, tmp, unpadded, solutions[0]);
            if (input_is_zero) {
                if (tcan_print)
                    fprintf(stderr,
                            "Using %sM^%u%s as input: Found zero vector."
                            " (coordinates on the padding part are %s zero)."
                            "\n"
                            "No solution found."
                            " Most certainly a bug."
                            "\n",
                            bw->dir ? "" : "V * ",
                            i-1,
                            bw->dir ? " * V" : "",
                            pad_is_zero ? "also" : "NOT");
                if (tcan_print && !pad_is_zero) {
                    char * tmp2;
                    int rc = asprintf(&tmp2, tmp, solutions[0], solutions[1]);
                    ASSERT_ALWAYS(rc >= 0);
                    fprintf(stderr,
                            "For reference, this useless vector (non-zero out, zero in) is stored in %s."
                            "\n",
                            tmp2);
                    free(tmp2);
                }
                serialize(pi->m);
                pthread_mutex_lock(pi->m->th->m);
                exitcode=1;
                pthread_mutex_unlock(pi->m->th->m);
                free(tmp);
                return NULL;
            }
            free(tmp);
        }
        /* }}} */

        if (tcan_print) {
            char * hwinfo = NULL;
            if (hamming_out && tcan_print) {
                int rc = asprintf(&hwinfo, ", Hamming weight is %d", hamming_out);
                ASSERT_ALWAYS(rc >= 0);
            }
            if (R) {
                const char * strings[] = {
                    "V * M^%u + R is %s%s\n", /* unsupported anyway */
                    "M^%u * V + R is %s%s\n",};
                printf(strings[bw->dir], i, is_zero ? "zero" : "NOT zero",
                        hwinfo ? hwinfo : "");
            } else {
                const char * strings[] = {
                    "V * M^%u is %s%s [K.sols%u-%u.%u contains V * M^%u]!\n",
                    "M^%u * V is %s%s [K.sols%u-%u.%u contains M^%u * V]!\n",
                };
                printf(strings[bw->dir], i,
                        is_zero ? "zero" : "NOT zero",
                        hwinfo ? hwinfo : "",
                        solutions[0], solutions[1], i-1, i-1);
            }
            if (hamming_out && tcan_print)
                free(hwinfo);
        }

        if (is_zero) {
            winning_iter = i-1;
            break;
        }

    }

    /* we exit with an untwisted vector. */

    if (!is_zero) {
        int nb_nonzero_coeffs=mmt_vec_hamming_weight(y);
        if (tcan_print) {
            printf("Solution range %u..%u: no solution found [%d non zero coefficients in result], most probably a bug\n", solutions[0], solutions[1], nb_nonzero_coeffs);
            if (nb_nonzero_coeffs < bw->n) {
                printf("There is some likelihood that by combining %d different solutions (each entailing a separate mksol run), or even with iterates of a single solution (then with a single mksol run) a full solution can be recovered. Ask for support.\n", nb_nonzero_coeffs + 1);
                /* A case that we do see is with
                 * test_bwc_modp_homogeneous_minimal_mn4 and seed=8 (in
                 * 64-bit mode). Iterates are confined to a 1-dimensional
                 * space, which means that it is easy to combine S and
                 * M*S into a full solution. I'm a bit lazy, though.
                 */
            }
        }
        pthread_mutex_lock(pi->m->th->m);
        exitcode=1;
        pthread_mutex_unlock(pi->m->th->m);
        return NULL;
    }

    if (leader) {
        /* About the ASSERT below: The 0-characteristic subspace
         * nilpotent part is not really a topic of interest in the case
         * of inhomogenous linear systems.  Exiting the loop with i>1 in
         * the inhomogenous case is an error, it seems, and I believe
         * that the solution provided was never really a solution. On the
         * contrary, the i>1 case in homogenous situations is perfectly
         * valid, and we're handling it now instead of using
         * winning_iter==0 always.
         */
        ASSERT_ALWAYS(winning_iter == 0 || !R);
        int splitwidth = char2 ? 64 : 1;
        unsigned int Av_multiplex = (solutions[1] - solutions[0]) / splitwidth;
        for(unsigned int i = 0 ; i < Av_multiplex ; ++i) {
            unsigned int sol0 = solutions[0] + i * splitwidth;
            unsigned int sol1 = sol0 + splitwidth;
            char * tmp;
            int rc = asprintf(&tmp, "K.sols%u-%u.%u", sol0, sol1, winning_iter);
            ASSERT_ALWAYS(rc >= 0);
            /* {{{ append the RHS coefficients if relevant */
            if (R) {
                printf("Expanding %s so as to include the coefficients for the %d RHS columns\n", tmp, R.nrhs);
                FILE * f = fopen(tmp, "ab");
                ASSERT_ALWAYS(f);
                rc = fseek(f, 0, SEEK_END);
                ASSERT_ALWAYS(rc >= 0);
                R.fwrite_rhs_coeffs(f, i);
                fclose(f);
            }

            printf("%s is now a %s nullspace vector for %s"
                    " (for a square M of dimension %" PRIu32"x%" PRIu32").\n",
                    tmp,
                    bw_dirtext[bw->dir],
                    R ? "(M|RHS)" : "M",
                    unpadded, unpadded);
            /* }}} */

            char * tmp2;
            rc = asprintf(&tmp2, "K.sols%u-%u.%u.txt", sol0, sol0 + splitwidth, winning_iter);
            ASSERT_ALWAYS(rc >= 0);
            /* {{{ write an ascii version while we're at it */
            {
                FILE * f  = fopen(tmp, "rb");
                FILE * f2 = fopen(tmp2, "w");
                ASSERT_ALWAYS(f);
                ASSERT_ALWAYS(f2);
                arith_generic::elt * data = A->alloc(1);
                for(uint32_t i = 0 ; i < mmt.n0[bw->dir] ; i++) {
                    size_t rc = fread(data, A->elt_stride() / Av_multiplex, 1, f);
                    ASSERT_ALWAYS(rc == 1);
                    std::ostringstream os;
                    A->cxx_out(os, *data);
                    fprintf(f2, "%s\n", os.str().c_str());
                }
                A->free(data);
                R.fprint_rhs_coeffs(f2);
                fclose(f);
                fclose(f2);
            }

            printf("%s (in ascii) is now a %s nullspace vector for %s"
                    " (for the original M, of dimension %" PRIu32"x%" PRIu32").\n",
                    tmp2,
                    bw_dirtext[bw->dir],
                    R ? "(M|RHS)" : "M",
                    mmt.n0[0], mmt.n0[1]);
            /* }}} */
            free(tmp2);
            free(tmp);
        }
    }

    serialize(pi->m);

    return NULL;
}

// coverity[root_function]
int main(int argc, char * argv[])
{
    param_list pl;

    bw_common_init(bw, &argc, &argv);
    param_list_init(pl);
    parallelizing_info_init();

    bw_common_decl_usage(pl);
    parallelizing_info_decl_usage(pl);
    matmul_top_decl_usage(pl);
    /* declare local parameters and switches */
    param_list_decl_usage(pl, "rhs",
            "file with the right-hand side vectors for inhomogeneous systems mod p (only the header is read by this program, while the actual contents are recovered from the V*.0 files)");
    param_list_decl_usage(pl, "rhscoeffs",
            "for the solution vector(s), this corresponds to the contribution(s) on the columns concerned by the rhs");

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    param_list_remove_key(pl, "interleaving");

    bw_common_interpret_parameters(bw, pl);
    parallelizing_info_lookup_parameters(pl);
    matmul_top_lookup_parameters(pl);
    /* interpret our parameters */

    ASSERT_ALWAYS(!param_list_lookup_string(pl, "ys"));
    ASSERT_ALWAYS(param_list_lookup_string(pl, "solutions"));

    param_list_lookup_string(pl, "rhs");
    param_list_lookup_string(pl, "rhscoeffs");

    if (param_list_warn_unused(pl)) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank) param_list_print_usage(pl, bw->original_argv[0], stderr);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    pi_go(gather_prog, pl, 0);

    parallelizing_info_finish();
    param_list_clear(pl);
    bw_common_clear(bw);

    return exitcode;
}

