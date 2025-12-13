#include "cado.h" // IWYU pragma: keep

/* This is the interface for doing products of polynomial matrices by
 * caching transforms, and transferring them over MPI. There are
 * currently two pitfalls.
 *
 * 1/ we have no streaming of computations and communications. This means
 * two things. First, communication/computation overlap can be
 * beneficial: while transforms of rows on the left term of a product are
 * computed (assuming transforms of the right term are ready), it is
 * possible to communicate the already computed row transforms. A second
 * point is that if we manage to get these computed soon, then there is
 * no need to keep them in memory. Such an approach would make it
 * possible to reduce the overall memory footprint. The place to fix this
 * is to modify the control flow in bigmatpoly_mul_caching_adj
 *
 * 2/ we are transferring full-length transforms, while it would 
 * be sufficient to transfer the truncated transforms of course.
 */

#include <climits>
#include <cstdint>
#include <cstdio>

#include <array>

#include "lingen_bigmatpoly_ft.hpp"
#include "lingen_fft_select.hpp"
#include "lingen_matpoly_bigmatpoly_ft_common.hpp"
#include "lingen_matpoly_ft.hpp"
#include "lingen_memory_pool.hpp"
#include "lingen_mul_substeps.hpp"
#include "macros.h"
#include "select_mpi.h"

template<bool is_binary>
class matpoly;
class tree_stats;
template <typename T, typename fft_type> struct OP_CTX;

template<typename fft_type> struct OP_CTX<bigmatpoly<is_binary_fft<fft_type>::value>, fft_type> : public OP_CTX_base<bigmatpoly<is_binary_fft<fft_type>::value>> {
    MPI_Datatype mpi_ft;
    tree_stats & stats;
    static constexpr bool is_binary = is_binary_fft<fft_type>::value;
    using T = bigmatpoly<is_binary>;
    using FFT = fft_type;
    template<typename... Args>
    OP_CTX(tree_stats & stats, fft_type const & fti, Args&&... args) : OP_CTX_base<T>(args...), stats(stats) {
        MPI_Type_contiguous(fti.get_alloc_sizes()[0], MPI_BYTE, &mpi_ft);
        MPI_Type_commit(&mpi_ft);
    }
    ~OP_CTX() {
        MPI_Type_free(&mpi_ft);
    }
    using OP_CTX_base<bigmatpoly<is_binary>>::a;
    using OP_CTX_base<bigmatpoly<is_binary>>::b;
    using OP_CTX_base<bigmatpoly<is_binary>>::c;

    int a_irank() const { return a.irank(); }
    int b_irank() const { return b.irank(); }
    int a_jrank() const { return a.jrank(); }
    int b_jrank() const { return b.jrank(); }
    int mesh_inner_size() const { return a.n1; }
    static const bool uses_mpi = true;
    void mesh_checks() const {
        ASSERT_ALWAYS(a.get_model().is_square());
        ASSERT_ALWAYS(a.get_model() == b.get_model());
        ASSERT_ALWAYS(a.irank() == b.irank());
        ASSERT_ALWAYS(a.jrank() == b.jrank());
    }
    void alloc_c_if_needed(size_t size) const {
        if (c.m != a.m || c.n != a.n || c.get_size() != size)
            c = T(a.ab, a.get_model(), a.m, b.n, size);
    }
    matpoly<is_binary> const & a_local() { return a.my_cell(); }
    matpoly<is_binary> const & b_local() { return b.my_cell(); }
    matpoly<is_binary> & c_local() { return c.my_cell(); }
    private:
    void do_allgather(void * p, MPI_Comm com, int n) {
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, 
                p, n,
                mpi_ft, com);
    }
    public:
    void a_allgather(void * p, int n) {
        do_allgather(p, a.get_model().com[1], n);
    }
    void b_allgather(void * p, int n) {
        do_allgather(p, b.get_model().com[2], n);
    }
    template<typename OP> void doit(OP & op, lingen_call_companion::mul_or_mp_times * M) {
        size_t ram = SIZE_MAX;
        if (M) {
            std::array<size_t, 3> alloc_sizes = op.fti.get_alloc_sizes();
            ram = M->ram(alloc_sizes);
            if (ram > M->ram()) {
                int rank;
                MPI_Comm_rank(a.get_model().com[0], &rank);
                if (rank == 0)
                    fprintf(stderr, "Transform size for %s with input operand sizes (%zu, %zu) is (%zu,%zu,%zu), which exceeds expected (%zu,%zu,%zu) (anticipated for operand sizes (%zu, %zu). Adjusting memory\n",
                            op.op_name(),
                            a.get_size(),
                            b.get_size(),
                            alloc_sizes[0],
                            alloc_sizes[1],
                            alloc_sizes[2],
                            M->fft_alloc_sizes[0],
                            M->fft_alloc_sizes[1],
                            M->fft_alloc_sizes[2],
                            M->asize,
                            M->bsize
                           );
                /* save it to the object. this means that we must pay
                 * attention to taking it by reference in the calling
                 * function.
                 */
                M->fft_alloc_sizes = alloc_sizes;
            }
        }
        try {
            typename matpoly_ft<fft_type>::memory_guard const dummy(ram);
            mp_or_mul<OP_CTX<bigmatpoly<is_binary>, fft_type>, OP>(*this, op, M)();
        } catch (memory_pool_exception const & e) {
            int rank;
            MPI_Comm_rank(a.get_model().com[0], &rank);
            fprintf(stderr, "Rank %d raised a memory pool exception: %s\n", rank, e.what());
            typename matpoly_ft<fft_type>::memory_guard const dummy(SIZE_MAX);
            mp_or_mul<OP_CTX<bigmatpoly<is_binary>, fft_type>, OP>(*this, op, M)();
        }
    }
};

template<typename fft_type>
auto bigmatpoly_ft<fft_type>::mp_caching(tree_stats & stats, bigmatpoly<is_binary> const & a, bigmatpoly<is_binary> const & b, lingen_call_companion::mul_or_mp_times * M) -> bigmatpoly<is_binary>
{
    op_mp<is_binary, fft_type> op(a, b, UINT_MAX);
    bigmatpoly<is_binary> c(a.get_model());
    OP_CTX<bigmatpoly<is_binary>, fft_type>(stats, op.fti, c, a, b).doit(op, M);
    return c;
}
template<typename fft_type>
auto bigmatpoly_ft<fft_type>::mul_caching(tree_stats & stats, bigmatpoly<is_binary> const & a, bigmatpoly<is_binary> const & b, lingen_call_companion::mul_or_mp_times * M) -> bigmatpoly<is_binary>
{
    op_mul<is_binary, fft_type> op(a, b, UINT_MAX);
    bigmatpoly<is_binary> c(a.get_model());
    OP_CTX<bigmatpoly<is_binary>, fft_type>(stats, op.fti, c, a, b).doit(op, M);
    return c;
}

#ifdef LINGEN_BINARY
template class bigmatpoly_ft<gf2x_fake_fft_info>;
template class bigmatpoly_ft<gf2x_cantor_fft_info>;
template class bigmatpoly_ft<gf2x_ternary_fft_info>;
#else
template class bigmatpoly_ft<fft_transform_info>;
#endif

