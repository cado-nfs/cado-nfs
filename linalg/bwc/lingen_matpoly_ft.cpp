#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cstdio>

#include "lingen_fft_select.hpp"
#include "lingen_mul_substeps.hpp"
#include "lingen_matpoly_ft.hpp"
#include "lingen_matpoly_bigmatpoly_ft_common.hpp"

template <typename T, typename fft_type> struct OP_CTX;

template<typename fft_type>
struct OP_CTX<matpoly<is_binary_fft<fft_type>::value>, fft_type>
    : public OP_CTX_base<matpoly<is_binary_fft<fft_type>::value>>
{
    tree_stats & stats;
    static constexpr bool is_binary = is_binary_fft<fft_type>::value;
    using T = matpoly<is_binary>;
    using FFT = fft_type;
    template<typename... Args>
    OP_CTX(tree_stats & stats, Args&&... args) : OP_CTX_base<T>(args...), stats(stats) {}
    int a_irank() const { return 0; }
    int b_irank() const { return 0; }
    int a_jrank() const { return 0; }
    int b_jrank() const { return 0; }
    int mesh_inner_size() const { return 1; }
    static const bool uses_mpi = false;
    void mesh_checks() const { }
    using OP_CTX_base<matpoly<is_binary>>::a;
    using OP_CTX_base<matpoly<is_binary>>::b;
    using OP_CTX_base<matpoly<is_binary>>::c;
    void alloc_c_if_needed(size_t size) {
        if (c.m != a.m || c.n != a.n || c.capacity() != size)
            c = T(a.ab, a.m, b.n, size);
    }
    T const & a_local()const  { return a; }
    T const & b_local() const { return b; }
    T & c_local() { return c; }
    void a_allgather(void *, int) const {}
    void b_allgather(void *, int) const {}
    template<typename OP> void doit(OP & op, lingen_call_companion::mul_or_mp_times * M) {
        size_t ram = SIZE_MAX;
        if (M) {
            std::array<size_t, 3> alloc_sizes = op.fti.get_alloc_sizes();
            ram = M->ram(alloc_sizes);
            if (ram > M->ram()) {
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
            mp_or_mul<OP_CTX<T, fft_type>, OP>(*this, op, M)();
        } catch (memory_pool_exception const & e) {
            fprintf(stderr, "Memory pool exception: %s\n", e.what());
            typename matpoly_ft<fft_type>::memory_guard const dummy(SIZE_MAX);
            mp_or_mul<OP_CTX<T, fft_type>, OP>(*this, op, M)();
        }
    }
};

template<typename fft_type> typename matpoly_ft<fft_type>::memory_pool_type matpoly_ft<fft_type>::memory;

template<typename fft_type>
auto matpoly_ft<fft_type>::mp_caching_adj(tree_stats & stats, matpoly<is_binary> const & a, matpoly<is_binary> const & b, unsigned int adj, lingen_call_companion::mul_or_mp_times * M) -> matpoly<is_binary>
    /*{{{*/
{
    op_mp<is_binary, fft_type> op(a, b, adj);
    matpoly<is_binary> c;
    OP_CTX<matpoly<is_binary>, fft_type>(stats, c, a, b).doit(op, M);
    return c;
} /* }}} */

template<typename fft_type>
auto matpoly_ft<fft_type>::mul_caching_adj(tree_stats & stats, matpoly<is_binary> const & a, matpoly<is_binary> const & b, unsigned int adj, lingen_call_companion::mul_or_mp_times * M) -> matpoly<is_binary>/*{{{*/
{
    op_mul<is_binary, fft_type> op(a, b, adj);
    matpoly<is_binary> c;
    OP_CTX<matpoly<is_binary>, fft_type>(stats, c, a, b).doit(op, M);
    return c;
} /* }}} */

#ifdef LINGEN_BINARY
template class matpoly_ft<gf2x_fake_fft_info>;
template class matpoly_ft<gf2x_cantor_fft_info>;
template class matpoly_ft<gf2x_ternary_fft_info>;
#else
template class matpoly_ft<fft_transform_info>;
#endif

