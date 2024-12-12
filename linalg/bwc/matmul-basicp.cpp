#include "cado.h" // IWYU pragma: keep

#include <cstdarg>         // for va_list, va_end, va_start
#include <cstddef>         // for ptrdiff_t
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <climits>
#include <gmp.h>
#include "matmul.hpp"
#include "matmul-common.hpp"
#include "arith-hard.hpp"

#include "matmul_facade.hpp"

#include "arith-modp.hpp"
#include "macros.h"
#include "params.h"

/* This extension is used to distinguish between several possible
 * implementations of the product */
#define MM_EXTENSION   "-basicp"
#define MM_MAGIC_FAMILY        0xb001UL
#define MM_MAGIC_VERSION       0x1001UL
#define MM_MAGIC (MM_MAGIC_FAMILY << 16 | MM_MAGIC_VERSION)

/* This selects the default behaviour as to which is our best code
 * for multiplying. If this flag is 1, then a multiplication matrix times
 * vector (direction==1) performs best if the in-memory structure
 * reflects the non-transposed matrix. Similarly, a vector times matrix
 * multiplication (direction==0) performs best if the in-memory structure
 * reflects the transposed matrix. When the flag is 1, the converse
 * happens.
 * This flag depends on the implementation, and possibly even on the cpu
 * type under certain circumstances.
 */
#define MM_DIR0_PREFERS_TRANSP_MULT   1

template<typename Arith>
struct matmul_basicp : public matmul_interface {
    /* repeat the fields from the public interface */
    /* now our private fields */
    arith_hard * xab;
    std::vector<uint32_t> q;

    void build_cache(matrix_u32 &&) override;
    int reload_cache_private() override;
    void save_cache_private() override;
    void mul(void *, const void *, int) override;
    ~matmul_basicp() override = default;

    matmul_basicp(matmul_public &&, arith_concrete_base *, cxx_param_list &, int);

    matmul_basicp(matmul_basicp const &) = delete;
    matmul_basicp& operator=(matmul_basicp const &) = delete;
    matmul_basicp(matmul_basicp &&) noexcept = default;
    matmul_basicp& operator=(matmul_basicp &&) noexcept = default;
};

template<typename Arith>
matmul_basicp<Arith>::matmul_basicp(
        matmul_public && P,
        arith_concrete_base * pxx,
        cxx_param_list & pl,
        int optimized_direction)
    : matmul_interface(std::move(P))
    , xab((arith_hard *) pxx) // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
{
    int const suggest = optimized_direction ^ MM_DIR0_PREFERS_TRANSP_MULT;
    store_transposed = suggest;
    param_list_parse(pl, "mm_store_transposed", 
            store_transposed);
    if (store_transposed != suggest) {
        fprintf(stderr, "Warning, mm_store_transposed"
                " overrides suggested matrix storage ordering\n");
    }           
}

template<typename Arith>
void matmul_basicp<Arith>::build_cache(matrix_u32 && m)
{
    q = std::move(m.p);
}

template<typename Arith>
int matmul_basicp<Arith>::reload_cache_private()
{
    auto f = matmul_common_reload_cache_fopen(sizeof(arith_hard::elt), *this, MM_MAGIC);
    if (!f) return 0;

    uint32_t datasize;
    MATMUL_COMMON_READ_ONE32(datasize, f.get());
    resize_and_check_meaningful(q, datasize, f.get());
    MATMUL_COMMON_READ_MANY32(q.data(), datasize, f.get());

    return 1;
}

template<typename Arith>
void matmul_basicp<Arith>::save_cache_private()
{
    auto f = matmul_common_save_cache_fopen(sizeof(arith_hard::elt), *this, MM_MAGIC);
    if (!f) return;

    MATMUL_COMMON_WRITE_ONE32(q.size(), f.get());
    MATMUL_COMMON_WRITE_MANY32(q.data(), q.size(), f.get());
}

template<typename Arith>
void matmul_basicp<Arith>::mul(void * xdst, void const * xsrc, int d)
{
    ASM_COMMENT("multiplication code");
    arith_hard * x = xab;
    auto const * src = (const typename Arith::elt *) xsrc;
    auto * dst = (typename Arith::elt *) xdst;

    /* d == 1: matrix times vector product */
    /* d == 0: vector times matrix product */

    /* However the matrix may be stored either row-major
     * (store_transposed == 0) or column-major (store_transposed == 1)
     */

    x->vec_set_zero(dst, dim[!d]);

    uint32_t const * qq = q.data();

    // NOLINTBEGIN(readability-static-accessed-through-instance)
    if (d == !store_transposed) {
        ARITH_MODP_TEMPORARY_ALLOC(x, elt_ur_for_add, rowsum);
        ASM_COMMENT("critical loop");
        for(unsigned int i = 0 ; i < dim[!d] ; i++) {
            uint32_t len = *qq++;
            unsigned int j = 0;
            x->set_zero(rowsum);
            for( ; len-- ; ) {
                j = *qq++;
                auto const c = int32_t(*qq++);
                ASSERT(j < dim[d]);
                if (c == 1) {
                    x->add(rowsum, x->vec_item(src, j));
                } else if (c == -1) {
                    x->sub(rowsum, x->vec_item(src, j));
                } else if (c > 0) {
                    x->addmul_ui(rowsum, x->vec_item(src, j), c);
                } else {
                    x->submul_ui(rowsum, x->vec_item(src, j), -c);
                }
            }
            x->reduce(x->vec_item(dst, i), rowsum);
        }
        ASM_COMMENT("end of critical loop");
    } else {
        auto * tdst = x->alloc<typename Arith::elt_ur_for_add>(dim[!d]);
        x->vec_set_zero(tdst, dim[!d]);
        if (iteration[d] == 10) {
            fprintf(stderr, "Warning: Doing many iterations with transposed code (not a huge problem for impl=basicp)\n");
        }
        ASM_COMMENT("critical loop (transposed mult)");
        for(unsigned int i = 0 ; i < dim[d] ; i++) {
            uint32_t len = *qq++;
            unsigned int j = 0;
            for( ; len-- ; ) {
                j = *qq++;
                auto const c = int32_t(*qq++);
                ASSERT(j < dim[!d]);
                if (c == 1) {
                    x->add(x->vec_item(tdst, j), x->vec_item(src, i));
                } else if (c == -1) {
                    x->sub(x->vec_item(tdst, j), x->vec_item(src, i));
                } else if (c > 0) {
                    x->addmul_ui(x->vec_item(tdst, j), x->vec_item(src, i), c);
                } else {
                    x->submul_ui(x->vec_item(tdst, j), x->vec_item(src, i), -c);
                }
            }
        }
        for(unsigned int j = 0 ; j < dim[!d] ; j++) {
            x->reduce(x->vec_item(dst, j), x->vec_item(tdst, j));
        }
        ASM_COMMENT("end of critical loop (transposed mult)");
        x->free<typename Arith::elt_ur_for_add>(tdst);
    }
    ASM_COMMENT("end of multiplication code");
    // NOLINTEND(readability-static-accessed-through-instance)

    iteration[d]++;
}

// NOLINTNEXTLINE(misc-use-internal-linkage)
matmul_interface * CADO_CONCATENATE4(new_matmul_, ARITH_LAYER, _, MM_IMPL)(
        matmul_public && P,
        arith_generic * arith,
        cxx_param_list & pl,
        int optimized_direction)
{
    return new matmul_basicp<arith_hard>(std::move(P), arith->concrete(), pl, optimized_direction);
}

/* vim: set sw=4: */
