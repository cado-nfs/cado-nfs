#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstring>
#include <climits>
#include <cstdint>
#include <ctime>
#include <cstdlib>

#include <algorithm>
#include <array>
#include <stdexcept>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <pthread.h>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "gmp_aux.h"
#include "bw-common.h"
#include "cxx_mpz.hpp"
#include "matmul_top.hpp"
#include "matmul_top_comm.hpp"
#include "matmul_top_vec.hpp"
#include "arith-generic.hpp"
#include "arith-cross.hpp"
#include "abase_proxy.hpp"
#include "parallelizing_info.hpp"
#include "params.h"
#include "portability.h"
#include "select_mpi.h"
#include "verbose.h"
#include "macros.h"
#include "mmt_vector_pair.hpp"
#include "utils_cxx.hpp"
#include "bwc_filenames.hpp"

static int exitcode = 0;

static std::vector<bwc_iteration_range> prelude(parallelizing_info_ptr pi)/*{{{*/
{
    int const leader = pi->m->jrank == 0 && pi->m->trank == 0;
    std::vector<bwc_iteration_range> iteration_ranges;
    serialize_threads(pi->m);

    if (leader) {
        std::map<bwc_iteration_range, std::vector<bwc_S_file>> by_iter;
        for(auto const & S : bwc_file_base::ls<bwc_S_file>(".")) {
            /* filter out the files that do not correspond to the
             * solution range we're interested in! */
            if (S.s0 < bw->solutions[0] || S.s0 >= bw->solutions[1])
                continue;
            if (S.s1 <= bw->solutions[0] || S.s1 > bw->solutions[1])
                continue;

            /* sanity check */
            if (bw->interval && S.n1 % bw->interval != 0)
                fmt::print(stderr,
                        "Warning: {} is not a checkpoint at a multiple of "
                        "the interval value {} -- this might indicate a "
                        "severe bug with leftover data, likely to corrupt "
                        "the final computation\n", S, bw->interval);

            /* prepare for hole checking */
            by_iter[S].push_back(S);
        }

        unsigned int prev_iter = 0;
        for(auto const & cur : by_iter) {
            unsigned int n0 = cur.first[0];
            unsigned int n1 = cur.first[1];
            if (n0 != prev_iter)
                throw std::runtime_error(fmt::format(
                            "Within the set of S files, file(s) "
                            "{} seems come first "
                            "after iteration {}, therefore there is "
                            "a gap for the range {}..{}\n",
                            join(cur.second, " "), prev_iter, prev_iter, n0));
            prev_iter = n1;

            /* It also makes sense to verify that we have a collection of
             * files that do match our solution range. The mmt_vec_load
             * will collect files based on that.
             */
            unsigned int prev_s = bw->solutions[0];
            for(auto const & S : cur.second) {
                if (S.s0 != prev_s)
                    throw std::runtime_error(fmt::format(
                                "For iterations {}-{}, we don't have"
                                " the files for solution columns {}-{}",
                                n0, n1,
                                prev_s, S.s0));
                prev_s = S.s1;
            }
            if (prev_s != bw->solutions[1])
                throw std::runtime_error(fmt::format(
                            "For iterations {}-{}, we don't have"
                            " the files for solution columns {}-{}",
                            n0, n1,
                            prev_s, bw->solutions[1]));

            iteration_ranges.emplace_back(cur.first);
        }
    }

    serialize(pi->m);
    unsigned long s = iteration_ranges.size();
    pi_bcast(&s, 1, BWC_PI_UNSIGNED_LONG, 0, 0, pi->m);
    if (!leader)
        iteration_ranges.assign(s, {});
    pi_bcast(iteration_ranges.data(),
            s * sizeof(decltype(iteration_ranges)::value_type), BWC_PI_BYTE,
            0, 0, pi->m);
    serialize(pi->m);
    return iteration_ranges;
}/*}}}*/

static void fprint_signed(FILE * f, arith_generic * A, arith_generic::elt const & x)
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
    fmt::print(f, "{}", os.str());
}

static std::vector<unsigned int> indices_of_zero_or_nonzero_values(mmt_vec & y, unsigned int maxidx, int want_nonzero)/*{{{*/
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
        parallelizing_info_experimental::allgather(myz, pi->wr[!y.d]);
    }

    parallelizing_info_experimental::broadcast(myz, pi);    /* And broadcast that to everyone as well. */

    return myz;
}/*}}}*/

static std::vector<unsigned int> indices_of_zero_values(mmt_vec & y, unsigned int maxidx)/*{{{*/
{
    return indices_of_zero_or_nonzero_values(y, maxidx, 0);
}/*}}}*/

static std::vector<unsigned int> indices_of_nonzero_values(mmt_vec & y, unsigned int maxidx)/*{{{*/
{
    return indices_of_zero_or_nonzero_values(y, maxidx, 1);
}/*}}}*/

static std::vector<unsigned int> get_possibly_wrong_columns(matmul_top_data & mmt)/*{{{*/
{
    parallelizing_info_ptr pi = mmt.pi;
    int const tcan_print = bw->can_print && pi->m->trank == 0;

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

    cxx_gmp_randstate rstate;

    if (pi->m->trank == 0 && !bw->seed) {
        /* note that bw is shared between threads, thus only thread 0 should
         * test and update it here.
         * at pi->m->jrank > 0, we don't care about the seed anyway
         */
        bw->seed = (int) time(nullptr);
        MPI_Bcast(&bw->seed, 1, MPI_INT, 0, pi->m->pals);
    }
    serialize_threads(pi->m);

    gmp_randseed_ui(rstate, bw->seed);
    if (tcan_print)
        fmt::print("// Random generator seeded with {}\n", bw->seed);

    /* Do that in the opposite direction compared to ymy */
    mmt_vector_pair zmz(mmt, !bw->dir);
    mmt_vec & z = zmz[0];
    mmt_vec & mz = zmz[zmz.size()-1];


    mmt_vec_set_random_inconsistent(z, rstate);
    mmt_vec_truncate_above_index(mmt, z, mmt.n0[bw->dir]);
    mmt_vec_apply_T(mmt, z);
    mmt_vec_twist(mmt, z);
    matmul_top_mul(mmt, zmz.vectors(), nullptr);
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
    matmul_top_mul(mmt, zmz.vectors(), nullptr);
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

    return allz;
}/*}}}*/

/* Take a vector, a priori sparse, and pick its coefficients at indices
 * marked by [rows], and put them, in order, in the j-th column of the
 * matrix pointed to by matrix, which has cblocks blocks of
 * A->simd_groupsize(A) entries. The column number j is thus made of elements
 * of the blocks whose index is congruent to
 * (j/groupsize)-th block mod cblocks ; A is my.abase.
 */
static void compress_vector_to_sparse(arith_generic::elt * matrix, unsigned int j, unsigned int cblocks, mmt_vec const & my, std::vector<unsigned int> & rows)
{
    arith_generic * A = my.abase;

    unsigned int const own_i0 = my.i0 + mmt_my_own_offset_in_items(my);
    unsigned int const own_i1 = own_i0 + mmt_my_own_size_in_items(my);
    cxx_mpz const v;
    unsigned int const jq = j / A->simd_groupsize();
    // unsigned int jr = j % A->simd_groupsize(A);

    int const char2 = mpz_cmp_ui(bw->p, 2) == 0;
    ASSERT_ALWAYS(char2 || A->simd_groupsize() == 1);

    for(unsigned int ii = 0 ; ii < rows.size() ; ii++) {
        unsigned int const i = rows[ii];
        if (own_i0 <= i && i < own_i1) {
            arith_generic::elt const & src = A->vec_item(my.v, i - my.i0);
            arith_generic::elt & dst = A->vec_item(matrix, ii * cblocks + jq);
            A->set(dst, src);
        }
    }
}

struct rhs /*{{{*/ {
    matmul_top_data & mmt;
    arith_generic * A = nullptr;
    unsigned int nrhs = 0;
    arith_generic::elt * rhscoeffs = nullptr;
    abase_proxy natural;
    arith_generic * Av = nullptr;

    rhs(rhs const&) = delete;
    rhs(rhs &&) = delete;
    rhs& operator=(rhs const&) = delete;
    rhs& operator=(rhs &&) = delete;

    rhs(matmul_top_data & mmt,
            const char * rhs_name,
            bwc_solution_range const & solutions)
        : mmt(mmt)
        , A(mmt.abase)
        , natural(abase_proxy::most_natural(mmt.pi)) /* {{{ */
    {
        if (!rhs_name) return;

        parallelizing_info_ptr pi = mmt.pi;
        int const tcan_print = bw->can_print && pi->m->trank == 0;
        int const leader = pi->m->jrank == 0 && pi->m->trank == 0;

        /* This is just for a check -- in truth, it might be that the
         * code here works correctly for inhomogeneous characteristic 2,
         * but that would be pure chance, as it was never tested */
        int const char2 = mpz_cmp_ui(bw->p, 2) == 0;
        ASSERT_ALWAYS(!char2);
        ASSERT_ALWAYS(A->simd_groupsize() == 1);

        Av = natural.A.get();

        if (leader)
            get_rhs_file_header(rhs_name, nullptr, &nrhs, nullptr);
        pi_bcast(&nrhs, 1, BWC_PI_UNSIGNED, 0, 0, pi->m);

        if (tcan_print) {
            fmt::print("** Informational note about GF(p) inhomogeneous system:\n");
            fmt::print("   Original matrix dimensions: {} {}\n",
                    mmt.n0[0], mmt.n0[1]);
            fmt::print("   We expect to obtain a vector of size {}\n",
                    mmt.n0[!bw->dir] + nrhs);
            fmt::print(
                    "   which we hope will be a kernel vector for (M_square||RHS).\n"
                    "   We will discard the coefficients which correspond to padding columns\n"
                    "   This entails keeping coordinates in the intervals\n"
                    "   [0..{}[ and [{}..{}[\n"
                    "   in the result.\n"
                    "** end note.\n",
                    mmt.n0[bw->dir], mmt.n0[!bw->dir], mmt.n0[!bw->dir] + nrhs);
        }

        /* everyone allocates something */
        rhscoeffs = A->alloc(nrhs, ALIGNMENT_ON_ALL_BWC_VECTORS);

        if (leader) {
            // yeah, we asserted that we're GF(p) at this point anyway.
            // coverity[dead_error_line]
            unsigned int const splitwidth = char2 ? 64 : 1;
            ASSERT_ALWAYS(Av->simd_groupsize() == splitwidth);

            if (char2 || solutions[1] != solutions[0] + splitwidth) {
                ASSERT_ALWAYS(0);/* never tested. I did attempt to code it right for the simd case though, but did not test. */
            }
            unsigned int const Av_multiplex = (solutions[1] - solutions[0]) / splitwidth;
            for(unsigned int j = 0 ; j < nrhs ; j++) {
                for(unsigned int i = 0 ; i < Av_multiplex ; i++) {
                    unsigned int const s0 = solutions[0] + i * splitwidth;
                    unsigned int const s1 = solutions[0] + (i + 1) * splitwidth;
                    std::string tmp = fmt::format("F.sols{}-{}.{}-{}.rhs", s0, s1, j, j+splitwidth);
                    if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_LOADING_MKSOL_FILES)) {
                        fmt::print("loading {}\n", tmp);
                    }
                    const auto f = fopen_helper(tmp, "rb");
                    const size_t rc = fread(
                            Av->vec_subvec(rhscoeffs, j * Av_multiplex + i),
                            Av->elt_stride(),
                            1, f.get());
                    ASSERT_ALWAYS(rc == 1);
                    if (Av->is_zero(Av->vec_item(rhscoeffs, j * Av_multiplex + i))) {
                        fmt::print("Notice: coefficient for vector V{}-{} in file {} is zero\n", j, j+1, tmp);
                    }
                }
            }
        }
        pi_bcast(rhscoeffs, nrhs, mmt.pitype, 0, 0, pi->m);
    }/*}}}*/

    ~rhs() {/*{{{*/
        if (rhscoeffs)
            A->free(rhscoeffs);
    }/*}}}*/

    void fwrite_rhs_coeffs(FILE * f, unsigned int i=0) const /* {{{ */
    {
        unsigned int const Av_multiplex = A->simd_groupsize() / Av->simd_groupsize();
        for(unsigned int j = 0 ; j < nrhs ; j++) {
            size_t const rc = fwrite(
                    Av->vec_subvec(rhscoeffs, j * Av_multiplex + i),
                    Av->elt_stride(),
                    1, f);
            ASSERT_ALWAYS(rc == 1);
        }
    }/*}}}*/

    void fprint_rhs_coeffs(FILE * f2) const /* {{{ */
    {
        for(uint32_t i = 0 ; i < nrhs ; i++) {
            std::ostringstream os;
            A->cxx_out(os, A->vec_item(rhscoeffs, i));
            fmt::print(f2, "{}\n", os.str());
        }
    }/*}}}*/

    explicit operator bool() const { return nrhs; }

    void add_contribution(mmt_vec & y) const/*{{{*/
    {
        if (!nrhs) return;

        parallelizing_info_ptr pi = mmt.pi;
        arith_generic * A = mmt.abase;
        ASSERT_ALWAYS(y.abase == A);
        unsigned int const unpadded = MAX(mmt.n0[0], mmt.n0[1]);
        size_t const eblock = mmt_my_own_size_in_items(y);

        abase_proxy natural = abase_proxy::most_natural(pi);
        arith_generic * Av = natural.A.get();
        pi_datatype_ptr Av_pi = natural.A_pi;

        mmt_vec vi(mmt,Av,Av_pi, bw->dir, /* shared ! */ 1, mmt.n[bw->dir]);

        for(unsigned int j = 0 ; j < nrhs ; j++) {
            int const ok = mmt_vec_load(vi, bwc_V_file::pattern(0), unpadded, j);
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

static std::tuple<int, int> check_zero_and_padding(mmt_vec & y, unsigned int maxidx)/*{{{*/
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
    size_t const my_i0 = y.i0 + mmt_my_own_offset_in_items(y);
    size_t const my_i1 = my_i0 + mmt_my_own_size_in_items(y);

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

    pi_allreduce(nullptr, &input_is_zero, 1, BWC_PI_INT, BWC_PI_MIN, y.pi->m);
    pi_allreduce(nullptr, &pad_is_zero, 1, BWC_PI_INT, BWC_PI_MIN, y.pi->m);

    return std::make_tuple(input_is_zero, pad_is_zero);
}/*}}}*/

static std::tuple<int, int, int> test_one_vector(matmul_top_data & mmt, mmt_vector_pair & ymy, rhs const & R)
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
        matmul_top_mul(mmt, ymy.vectors(), nullptr);
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
        pi_allreduce(nullptr, &is_zero, 1, BWC_PI_INT, BWC_PI_MIN, pi->m);

        hamming_out = is_zero ? 0 : mmt_vec_hamming_weight(y);
    }
    return std::make_tuple(input_is_zero, pad_is_zero, hamming_out);
}

/* Use y_saved as input (and leave it untouched). Store result in both y
 * and my */
static std::tuple<int, int, int> expanded_test(matmul_top_data & mmt, mmt_vector_pair & ymy, mmt_vec const & y_saved, rhs const& R)
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

    // using pre_matrix_t = std::map<std::pair<unsigned int, unsigned int>, cxx_mpz>;

    arith_generic::elt * matrix = nullptr;
    public:
    bool attempt_to_fix;

    std::vector<unsigned int> cols;
    std::vector<unsigned int> rows;
    size_t nrows() const { return rows.size(); }
    std::vector<std::pair<std::array<unsigned int, 2>, int> > pivot_list;

    explicit parasite_fixer(matmul_top_data & mmt)
        : mmt(mmt)
        , A(mmt.abase)
        , pi(mmt.pi)
    {/*{{{*/
        int const tcan_print = bw->can_print && pi->m->trank == 0;

        cols = get_possibly_wrong_columns(mmt);

        attempt_to_fix = !cols.empty() && cols.size() < 64;

        if (!cols.empty() && tcan_print) {
            fmt::print("# {} possibly wrong coordinates detected"
                    " in solution vector because the matrix"
                    " has a non-trivial nilpotent space.\n",
                    cols.size());
            if (attempt_to_fix) {
                fmt::print("# Will try to fix\n");
            } else {
                fmt::print("# Deliberately avoiding attempt to fix because"
                        " this number of vectors is large."
                        " It is rather an indication that something is wrong."
                        " Proceeding anyway\n");
            }
        }
        compute_pivot_list();
    }/*}}}*/

    parasite_fixer(parasite_fixer const&) = delete;
    parasite_fixer(parasite_fixer &&) = delete;
    parasite_fixer& operator=(parasite_fixer const&) = delete;
    parasite_fixer& operator=(parasite_fixer &&) = delete;

    /*{{{ row_coordinates_of_nonzero_cols(matmul_top_data & mmt) */
    std::vector<unsigned int> row_coordinates_of_nonzero_cols(matmul_top_data & mmt, std::vector<unsigned int> const& cols)
    {
        // int tcan_print = bw->can_print && pi->m->trank == 0;

        arith_generic * A = mmt.abase;

        mmt_vector_pair ymy(mmt, bw->dir);
        mmt_vec & y = ymy[0];
        mmt_vec & my = ymy[ymy.size()-1];

        /* Now try to see which indices are potentially affected */
        unsigned int const B = A->simd_groupsize();
        for(unsigned int jjq = 0 ; jjq < cols.size() ; jjq+=B) {
            mmt_full_vec_set_zero(y);
            for(unsigned int jjr = 0 ; jjr < B && (jjq + jjr < cols.size()) ; jjr++) {
                unsigned int const j = cols[jjq + jjr];
                mmt_vec_add_basis_vector_at(y, jjr, j);
            }
            mmt_vec_apply_T(mmt, y);
            mmt_vec_twist(mmt, y);
            matmul_top_mul(mmt, ymy.vectors(), nullptr);
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
        std::ranges::sort(rows);
        auto [ last, end ] = std::ranges::unique(rows);
        rows.erase(last, end);

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

        pi_allreduce(nullptr, matrix,
                rows.size() * cols.size(),
                mmt.pitype, BWC_PI_SUM, pi->m);
    }/*}}}*/
#endif

    void debug_print_local_matrix(arith_generic::elt * matrix,
            std::vector<unsigned int> const & rows,
            std::vector<unsigned int> const & cols,
            arith_generic::elt * nz = nullptr)/*{{{*/
    {
        size_t const nr = rows.size();
        size_t const nc = cols.size();
        fmt::print("# Dump of the full matrix{} as seen by J{}T{}\n",
                nz ? " (with coefficients of the vector encountered)" : "",
                pi->m->jrank, pi->m->trank);
        unsigned int kk = 0;
        unsigned int const B = A->simd_groupsize();
        unsigned int const cblocks = iceildiv(nc, B);
        fmt::print("#\t\t");
        for(unsigned int jj = 0 ; jj < nc ; jj++) {
            fmt::print("[{}] ", cols[jj]);
        }
        fmt::print("\n");
        for(unsigned int ii = 0 ; ii < nr ; ii++) {
            fmt::print("#\t{}\t", rows[ii]);
            for(unsigned int jj = 0 ; jj < cblocks ; jj++) {
                fmt::print(" ");
                fprint_signed(stdout, A, A->vec_item(matrix, kk));
                kk++;
            }
            if (nz) {
                fmt::print(" ");
                fprint_signed(stdout, A, A->vec_item(nz, ii));
            }
            fmt::print("\n");
        }
    }/*}}}*/

#if 0
    void debug_print_all_local_matrices(arith_generic::elt * nz = nullptr)/*{{{*/
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

        int const tcan_print = bw->can_print && pi->m->trank == 0;
        int const leader = pi->m->jrank == 0 && pi->m->trank == 0;

        rows = row_coordinates_of_nonzero_cols(mmt, cols);

        arith_generic * A = mmt.abase;
        arith_generic::elt * dummy = A->alloc();

        int const char2 = mpz_cmp_ui(bw->p, 2) == 0;

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

        unsigned int const B = A->simd_groupsize();

        matrix = A->alloc(iceildiv(cols.size(), B) * rows.size(), ALIGNMENT_ON_ALL_BWC_VECTORS);
        A->vec_set_zero(matrix, iceildiv(cols.size(), B) * rows.size());

        for(unsigned int drop = UINT_MAX, spin=0 ; drop && !scols.empty() ; spin++) {
            drop = 0;

            std::vector<unsigned int> vcols(scols.begin(), scols.end());
            std::vector<unsigned int> vrows;
            vrows.reserve(srows.size());
            for(auto x : srows) vrows.push_back(rows[x]);

            if (tcan_print)
                fmt::print("# Pass {}: checking for error fixing pivots within a matrix of dimension {}*{}\n", spin, vrows.size(), vcols.size());

            arith_generic::elt * mat;
            arith_generic::elt ** pmat = &mat;
            unsigned int const cblocks = iceildiv(vcols.size(), B);
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
                    unsigned int const j = vcols[jjq + jjr];
                    mmt_vec_add_basis_vector_at(y, jjr, j);
                }
                mmt_vec_apply_T(mmt, y);
                mmt_vec_twist(mmt, y);
                matmul_top_mul(mmt, ymy.vectors(), nullptr);
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

            pi_allreduce(nullptr, *pmat,
                    cblocks * vrows.size(),
                    mmt.pitype, BWC_PI_SUM, pi->m);

            if (leader) {
                fmt::print("Print matrix of size {}*{}\n", srows.size(), scols.size());
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
                size_t const w = A->vec_simd_hamming_weight(row, cblocks);
                if (w == 1) {
                    size_t const p = A->vec_simd_find_first_set(*dummy, row, cblocks);
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
                const unsigned int xi = pp.first;
                const unsigned int j = pp.second.first;
                int const v = pp.second.second;
                if (v == 2) continue;
                if (tcan_print)
                    fmt::print("Found pivot for row {}:"
                            " column {} has coefficient {}\n",
                            rows[xi], j, v);
                srows.erase(xi);
                if (scols.erase(j)) {
                    drop++;
                    /* XXX colum j may have already been deleted */
                    std::array<unsigned int, 2> const xij { xi, j };
                    pivot_list.emplace_back(xij, v);
                }
            }
            if (spin)
                A->free(mat);

            if (tcan_print)
                fmt::print("# Pass {}: number of cols has dropped by {}."
                        " We have {} rows (at most) and {} columns left\n",
                        spin, drop, srows.size(), scols.size());
            /* use this marker to indicate synchronization */
            std::array<unsigned int, 2> const xij { 0U, 0U };
            pivot_list.emplace_back(xij, 0);
        }
        if (!scols.empty()) {
            /* we exited with drop==0 then, so srows correctly reflects
             * the last set of active rows */
            for(auto r : srows) {
                auto it = pivots.find(r);
                if (it != pivots.end()) {
                    fmt::print("For row {}, only pivot found is non-trivial."
                            " Please implement that code\n", rows[r]);
                } else {
                    fmt::print("For row {}, linear algebra reduction is needed."
                            " Please implement real code\n", rows[r]);
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
        int const tcan_print = bw->can_print && pi->m->trank == 0;
        int const leader = pi->m->jrank == 0 && pi->m->trank == 0;

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
            fmt::print("# Input vector has {} input, {} padding\n",
                    input_is_zero ? "zero" : "non-zero",
                    pad_is_zero ? "zero" : "non-zero");
            fmt::print("# Output has Hamming weight {}\n", hamming_out);
        }
        if (!std::includes(rows.begin(), rows.end(), nz_pos.begin(), nz_pos.end())) {
            if (tcan_print)
                fmt::print("# Note: cannot attempting to fix {} wrong coordinates, not included in the set of {} known possibly wrong ones\n", nz_pos.size(), rows.size());
            return res;
        }

        if (tcan_print)
            fmt::print("# Note: all the non-zero coordinates are included in the output of the \"possibly wrong\" columns, which is a good sign\n");

        arith_generic::elt * nz;


        ASSERT_ALWAYS(my.abase == mmt.abase);

        nz = A->alloc(rows.size(), ALIGNMENT_ON_ALL_BWC_VECTORS);
        A->vec_set_zero(nz, rows.size());
        compress_vector_to_sparse(nz, 0, 1, my, rows);
        pi_allreduce(nullptr, nz, rows.size(), mmt.pitype, BWC_PI_SUM, pi->m); 

        if (leader) debug_print_local_matrix(matrix, rows, cols, nz);

        // debug_print_all_local_matrices(nz);

        for(auto pp : pivot_list) {
            unsigned int const ii = pp.first[0];
            // unsigned int i = rows[ii];
            unsigned int const j = pp.first[1];
            int const v = pp.second;
            if (!v) {
                serialize(pi->m);
                res = expanded_test(mmt, ymy, y_saved, R);
                std::tie(input_is_zero, pad_is_zero, hamming_out) = res;
                A->vec_set_zero(nz, rows.size());
                compress_vector_to_sparse(nz, 0, 1, my, rows);
                pi_allreduce(nullptr, nz, rows.size(), mmt.pitype, BWC_PI_SUM, pi->m); 

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

            size_t const own_i0 = y_saved.i0 + mmt_my_own_offset_in_items(y_saved);
            size_t const own_i1 = own_i0 + mmt_my_own_size_in_items(y_saved);

            if (own_i0 <= j && j < own_i1) {
                arith_generic::elt & source = A->vec_item(y_saved.v, j - y_saved.i0);
                fmt::print("Row {}, coefficient is ", rows[ii]);
                fprint_signed(stdout, A, source);
                if (v == -1) {
                    fmt::print(" ; fixing by adding to coordinate {}\n", j);
                    A->add_and_reduce(source, error);
                } else if (v == 1) {
                    fmt::print(" ; fixing by subtracting from coordinate {}\n", j);
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
            fmt::print("# After fix, Hamming weight is {}\n", hamming_out);

        // if (leader) debug_print_local_matrix(matrix, rows, cols, nz);

        A->free(nz);

        return res;
    }/*}}}*/
};/*}}}*/

static void * gather_prog(parallelizing_info_ptr pi, cxx_param_list & pl, void * arg MAYBE_UNUSED)
{
    ASSERT_ALWAYS(!pi->interleaved);

    int const tcan_print = bw->can_print && pi->m->trank == 0;
    int const leader = pi->m->jrank == 0 && pi->m->trank == 0;

    bwc_solution_range sol_range { bw->solutions[0], bw->solutions[1], };
    int const char2 = mpz_cmp_ui(bw->p, 2) == 0;

    /* Define and initialize our arithmetic back-ends. More or less the
     * same deal as for mksol (for the "solutions" part). See comments
     * there. */

    /* {{{ main arithmetic backend: for the solutions we compute. */
    // unsigned int A_multiplex MAYBE_UNUSED = 1;
    unsigned int const A_width = sol_range[1] - sol_range[0];
    if ((char2 && (A_width != 64 && A_width != 128 && A_width != 256))
            || (!char2 && A_width > 1))
    {
        fmt::print(stderr,
                "We cannot support computing {} solutions at a time "
                "with one single Spmv operation, given the currently "
                "implemented code\n",
                A_width);
        exit(EXIT_FAILURE);
    }

    abase_proxy const abase_solutions(pi, A_width);
    arith_generic * A = abase_solutions.A.get();

    /* }}} */

    matmul_top_data mmt(A, pi, pl, bw->dir);

    parasite_fixer pfixer(mmt);

    mmt_vector_pair ymy(mmt, bw->dir);
    mmt_vec & y = ymy[0];

    mmt_vec y_saved(mmt, nullptr, nullptr, bw->dir, mmt_vec_is_shared(y), mmt.n[bw->dir]);

    /* this is really a misnomer, because in the typical case, M is
     * rectangular, and then the square matrix does induce some padding.
     * This is however not the padding we are interested in. The padding
     * we're referring to in the naming of this variable is the one which
     * is related to the number of jobs and threads: internal dimensions
     * are arranged to be multiples */
    unsigned int const unpadded = MAX(mmt.n0[0], mmt.n0[1]);

    auto sl = prelude(pi);
    if (sl.empty()) {
        if (tcan_print) {
            fmt::print(stderr, "Found zero S files for solution range {}..{}. "
                    "Problem with command line ?\n",
                    sol_range[0], sol_range[1]);
            pthread_mutex_lock(pi->m->th->m);
            exitcode=1;
            pthread_mutex_unlock(pi->m->th->m);
        }
        serialize_threads(pi->m);
        return nullptr;
    }

    rhs R(mmt, param_list_lookup_string(pl, "rhs"), sol_range);

    if (tcan_print)
        fmt::print("Trying to build solutions {}..{}\n", sol_range[0], sol_range[1]);

    serialize(mmt.pi->m);

    { /* {{{ Collect now the sum of the LHS contributions */
        mmt_vec svec(mmt, nullptr, nullptr, bw->dir, /* shared ! */ 1, mmt.n[bw->dir]);
        for(auto const & S : sl) {
            auto pat = bwc_S_file::pattern(S);
            /* sl is an iteration range */
            if (tcan_print && verbose_enabled(CADO_VERBOSE_PRINT_BWC_LOADING_MKSOL_FILES)) {
                fmt::print("loading contributions from iterations {}-{}\n",
                        S.n0, S.n1);
            }
            int const ok = mmt_vec_load(svec, pat, unpadded, sol_range[0]);
            ASSERT_ALWAYS(ok);

            A->vec_add_and_reduce(
                    mmt_my_own_subvec(y),
                    mmt_my_own_subvec(svec),
                    mmt_my_own_size_in_items(y));
        }
        y.consistency = 1;
        mmt_vec_broadcast(y);
    } /* }}} */

    serialize(mmt.pi->m);

    auto w = mmt_vec_hamming_weight(y);

    if (tcan_print)
        fmt::print("Hamming weight of sum: {}\n", w);

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
        auto Kpat = bwc_K_file::pattern(i-1);
        if (input_is_zero)
            Kpat = std::string("zero").append(Kpat);
        bwc_K_file Kfile(sol_range, i-1);
        {
            mmt_vec_save(y_saved, Kpat, unpadded, sol_range[0]);
            if (input_is_zero) {
                if (tcan_print)
                    fmt::print(stderr,
                            "Using {}M^{}{} as input: Found zero vector."
                            " (coordinates on the padding part are {} zero)."
                            "\n"
                            "No solution found."
                            " Most certainly a bug."
                            "\n",
                            bw->dir ? "" : "V * ",
                            i-1,
                            bw->dir ? " * V" : "",
                            pad_is_zero ? "also" : "NOT");
                if (tcan_print && !pad_is_zero) {
                    fmt::print(stderr,
                            "For reference,"
                            " this useless vector (non-zero out, zero in)"
                            " is stored in zero{}.\n",
                            Kfile);
                }
                serialize(pi->m);
                pthread_mutex_lock(pi->m->th->m);
                exitcode=1;
                pthread_mutex_unlock(pi->m->th->m);
                return nullptr;
            }
        }
        /* }}} */

        if (tcan_print) {
            std::string hwinfo;
            std::string zinfo = is_zero ? "zero" : "NOT zero";
            if (hamming_out)
                hwinfo = fmt::format(", Hamming weight is {}", hamming_out);

            if (R) {
                const char * strings[] = {
                    "V * M^{} + R is {}{}\n", /* unsupported anyway */
                    "M^{} * V + R is {}{}\n",};
                fmt::print(fmt::runtime(strings[bw->dir]), i, zinfo, hwinfo);
            } else {
                const char * strings[] = {
                    "V * M^{} is {}{} [{} contains V * M^{}]!\n",
                    "M^{} * V is {}{} [{} contains M^{} * V]!\n",
                };
                fmt::print(fmt::runtime(strings[bw->dir]), i,
                        zinfo,
                        hwinfo,
                        Kfile, i-1);
            }
        }

        if (is_zero) {
            winning_iter = i-1;
            break;
        }

    }

    /* we exit with an untwisted vector. */

    if (!is_zero) {
        unsigned long const nb_nonzero_coeffs = mmt_vec_hamming_weight(y);
        if (tcan_print) {
            fmt::print("Solution range {}..{}: no solution found"
                    " [{} non zero coefficients in result],"
                    " most probably a bug\n",
                    sol_range[0], sol_range[1], nb_nonzero_coeffs);
            if (nb_nonzero_coeffs < (unsigned long) bw->n) {
                fmt::print("There is some likelihood that by combining"
                        " {} different solutions"
                        " (each entailing a separate mksol run),"
                        " or even with iterates of a single solution"
                        " (then with a single mksol run)"
                        " a full solution can be recovered."
                        " Ask for support.\n", nb_nonzero_coeffs + 1);
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
        return nullptr;
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
        int const splitwidth = char2 ? 64 : 1;
        unsigned int const Av_multiplex = (sol_range[1] - sol_range[0]) / splitwidth;
        for(unsigned int i = 0 ; i < Av_multiplex ; ++i) {
            auto tmp = std::string(bwc_K_file(sol_range, winning_iter));
            /* {{{ append the RHS coefficients if relevant */
            if (R) {
                fmt::print("Expanding {} so as to include"
                        " the coefficients for the {} RHS columns\n",
                        tmp, R.nrhs);
                auto f = fopen_helper(tmp, "ab");
                const int rc = fseek(f.get(), 0, SEEK_END);
                ASSERT_ALWAYS(rc >= 0);
                R.fwrite_rhs_coeffs(f.get(), i);
            }

            fmt::print("{} is now a {} nullspace vector for {}"
                    " (for a square M of dimension {}x{}).\n",
                    tmp,
                    bw_dirtext[bw->dir],
                    R ? "(M|RHS)" : "M",
                    unpadded, unpadded);
            /* }}} */

            auto tmp2 = tmp + ".txt";
            /* {{{ write an ascii version while we're at it */
            {
                auto f = fopen_helper(tmp, "rb");
                auto f2 = fopen_helper(tmp2, "w");
                arith_generic::elt * data = A->alloc(1);
                for(uint32_t i = 0 ; i < mmt.n0[bw->dir] ; i++) {
                    size_t const rc = fread(data, A->elt_stride() / Av_multiplex, 1, f.get());
                    ASSERT_ALWAYS(rc == 1);
                    std::ostringstream os;
                    A->cxx_out(os, *data);
                    fmt::print(f2.get(), "{}\n", os.str());
                }
                A->free(data);
                R.fprint_rhs_coeffs(f2.get());
            }

            fmt::print("{} (in ascii) is now a {} nullspace vector for {}"
                    " (for the original M, of dimension {}x{}).\n",
                    tmp2,
                    bw_dirtext[bw->dir],
                    R ? "(M|RHS)" : "M",
                    mmt.n0[0], mmt.n0[1]);
            /* }}} */
        }
    }

    serialize(pi->m);

    return nullptr;
}

// coverity[root_function]
int main(int argc, char const * argv[])
{
    cxx_param_list pl;

    bw_common_init(bw, &argc, &argv);

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

    pi_go(gather_prog, pl, nullptr);

    parallelizing_info_finish();

    bw_common_clear(bw);

    return exitcode;
}

