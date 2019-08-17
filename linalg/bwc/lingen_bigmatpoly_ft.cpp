#include "cado.h"
#include "mpfq_layer.h"
#include <cstdlib>
#include "portability.h"
#include "macros.h"
#include "lingen_matpoly_ft.hpp"
#include "lingen_bigmatpoly_ft.hpp"
#include "logline.h"
#include "fmt/format.h"

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

#include "lingen_matpoly_bigmatpoly_ft_common.hpp"

template<> struct OP_CTX<bigmatpoly> : public OP_CTX_base<bigmatpoly> {
    MPI_Datatype mpi_ft;
    tree_stats & stats;
    typedef bigmatpoly T;
    template<typename... Args>
    OP_CTX(tree_stats & stats, fft_transform_info const * fti, Args&&... args) : OP_CTX_base<T>(args...), stats(stats) {
        size_t fft_alloc_sizes[3];
        fft_get_transform_allocs(fft_alloc_sizes, fti);
        size_t tsize = fft_alloc_sizes[0];
        MPI_Type_contiguous(tsize, MPI_BYTE, &mpi_ft);
        MPI_Type_commit(&mpi_ft);
    }
    ~OP_CTX() {
        MPI_Type_free(&mpi_ft);
    }
    inline int a_irank() const { return a.irank(); }
    inline int b_irank() const { return b.irank(); }
    inline int a_jrank() const { return a.jrank(); }
    inline int b_jrank() const { return b.jrank(); }
    inline int mesh_size() const { return a.n1; }
    static const bool uses_mpi = true;
    inline void mesh_checks() const {
        ASSERT_ALWAYS(a.get_model().is_square());
        ASSERT_ALWAYS(a.get_model() == b.get_model());
        ASSERT_ALWAYS(a.irank() == b.irank());
        ASSERT_ALWAYS(a.jrank() == b.jrank());
        ASSERT_ALWAYS(a.get_model().com[1] == b.get_model().com[2]);
    }
    void alloc_c_if_needed(size_t size) const {
        if (c.m != a.m || c.n != a.n || c.size != size)
            c = T(a.ab, a.get_model(), a.m, b.n, size);
    }
    inline matpoly const & a_local() { return a.my_cell(); }
    inline matpoly const & b_local() { return b.my_cell(); }
    inline matpoly & c_local() { return c.my_cell(); }
    inline void do_allgather(void * p, int n) {
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, 
                p, n,
                mpi_ft, a.get_model().com[1]);
    }
    inline void begin_smallstep(std::string const & func, unsigned int n) {
        stats.begin_smallstep(func, n);
    }
    inline void end_smallstep() {
        stats.end_smallstep();
    }
    inline void skip_smallstep(std::string const & func, unsigned int n) {
        stats.skip_smallstep(func, n);
    }
    inline bool local_smallsteps_done() {
        return stats.local_smallsteps_done();
    }
    template<typename OP> void doit(OP & op, lingen_call_companion::mul_or_mp_times * M) {
        if (M && op.get_transform_ram() > M->per_transform_ram) {
            fprintf(stderr, "Transform size for %s with input operand sizes (%zu, %zu) is %zu, which exceeds expected %zu (anticipated for operand sizes (%zu, %zu). Updating\n",
                    OP::name,
                    a.size,
                    b.size,
                    op.get_transform_ram(),
                    M->per_transform_ram,
                    M->asize,
                    M->bsize
                   );
            size_t ntransforms = M->ram / M->per_transform_ram;
            ASSERT_ALWAYS(M->ram % M->per_transform_ram == 0);
            M->per_transform_ram = op.get_transform_ram();
            M->ram = ntransforms * M->per_transform_ram;
        }
        matpoly_ft::memory_pool_guard dummy(M ? M->ram : SIZE_MAX);
        mp_or_mul(*this, op, op.fti, M ? & M->S : NULL);
    }
};


void bigmatpoly_mp_caching_adj(tree_stats & stats, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, unsigned int adj, lingen_call_companion::mul_or_mp_times * M)
{
    op_mp OP(a, b, adj);
    OP_CTX<bigmatpoly>(stats, OP.fti, c, a, b).doit(OP, M);
}
void bigmatpoly_mul_caching_adj(tree_stats & stats, bigmatpoly & c, bigmatpoly const & a, bigmatpoly const & b, unsigned int adj, lingen_call_companion::mul_or_mp_times * M)
{
    op_mul OP(a, b, adj);
    OP_CTX<bigmatpoly>(stats, OP.fti, c, a, b).doit(OP, M);
}

