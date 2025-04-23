#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>

#include <utility>
#include <vector>

#include "arith-generic.hpp"
#include "arith-hard.hpp"
#include "macros.h"
#include "matmul.hpp"
#include "matmul-common.hpp"
#include "matrix_u32.hpp"
#include "params.h"

/* This extension is used to distinguish between several possible
 * implementations of the product */
#define MM_EXTENSION   "-basic"
#define MM_MAGIC_FAMILY        0xa001UL
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
struct matmul_basic : public matmul_interface {
    /* now our private fields */
    Arith * xab;
    std::vector<uint32_t> q;
    void build_cache(matrix_u32 &&) override;
    int reload_cache_private() override;
    void save_cache_private() override;
    void mul(void *, const void *, int) override;
    ~matmul_basic() override = default;

    matmul_basic(matmul_public &&, arith_concrete_base *, cxx_param_list &, int);

    matmul_basic(matmul_basic const &) = delete;
    matmul_basic& operator=(matmul_basic const &) = delete;
    matmul_basic(matmul_basic &&) noexcept = default;
    matmul_basic& operator=(matmul_basic &&) noexcept = default;
};

template<typename Arith>
matmul_basic<Arith>::matmul_basic(matmul_public && P, arith_concrete_base * pxx, cxx_param_list & pl, int optimized_direction)
    : matmul_interface(std::move(P))
    , xab((Arith *) pxx) // NOLINT(cppcoreguidelines-pro-type-cstyle-cast)
{
    int const suggest = optimized_direction ^ MM_DIR0_PREFERS_TRANSP_MULT;
    store_transposed = suggest;
    param_list_parse(pl, "mm_store_transposed", store_transposed);
    if (store_transposed != suggest) {
        fprintf(stderr, "Warning, mm_store_transposed"
                " overrides suggested matrix storage ordering\n");
    }           
}

template<typename Arith>
void matmul_basic<Arith>::build_cache(matrix_u32 && m)
{
    q = std::move(m.p);
}

template<typename Arith>
int matmul_basic<Arith>::reload_cache_private()
{
    auto f = matmul_common_reload_cache_fopen(sizeof(typename Arith::elt), *this, MM_MAGIC);
    if (!f) return 0;

    uint32_t datasize;

    MATMUL_COMMON_READ_ONE32(datasize, f.get());
    resize_and_check_meaningful(q, datasize, f.get());
    MATMUL_COMMON_READ_MANY32(q.data(), datasize, f.get());

    return 1;
}

template<typename Arith>
void matmul_basic<Arith>::save_cache_private()
{
    auto f = matmul_common_save_cache_fopen(sizeof(typename Arith::elt), *this, MM_MAGIC);
    if (!f) return;

    MATMUL_COMMON_WRITE_ONE32(q.size(), f.get());
    MATMUL_COMMON_WRITE_MANY32(q.data(), q.size(), f.get());
}

template<typename Arith>
void matmul_basic<Arith>::mul(void * xdst, void const * xsrc, int d)
{
    ASM_COMMENT("multiplication code");
    uint32_t const * qq = q.data();
    Arith * x = xab;
    auto const * src = (typename Arith::elt const *) xsrc;
    auto * dst = (typename Arith::elt *) xdst;

    if (d == !store_transposed) {
        x->vec_set_zero(dst, dim[!d]);
        ASM_COMMENT("critical loop");
        for(unsigned int i = 0 ; i < dim[!d] ; i++) {
            uint32_t len = *qq++;
            unsigned int j = 0;
            for( ; len-- ; ) {
                j = *qq++;
                ASSERT(j < dim[d]);
                x->add(x->vec_item(dst, i), x->vec_item(src, j));
            }
        }
        ASM_COMMENT("end of critical loop");
    } else {
        if (iteration[d] == 10) {
            fprintf(stderr, "Warning: Doing many iterations with transposed code (not a huge problem for impl=basic)\n");
        }
        x->vec_set_zero(dst, dim[!d]);
        ASM_COMMENT("critical loop (transposed mult)");
        for(unsigned int i = 0 ; i < dim[d] ; i++) {
            uint32_t len = *qq++;
            unsigned int j = 0;
            for( ; len-- ; ) {
                j = *qq++;
                ASSERT(j < dim[!d]);
                x->add(x->vec_item(dst, j), x->vec_item(src, i));
            }
        }
        ASM_COMMENT("end of critical loop (transposed mult)");
    }
    ASM_COMMENT("end of multiplication code");
    iteration[d]++;
}

// NOLINTNEXTLINE(misc-use-internal-linkage)
matmul_interface * CADO_CONCATENATE4(new_matmul_, ARITH_LAYER, _, MM_IMPL)(
        matmul_public && P,
        arith_generic * arith,
        cxx_param_list & pl,
        int optimized_direction)
{
    return new matmul_basic<arith_hard>(std::move(P), arith->concrete(), pl, optimized_direction);
}

/* vim: set sw=4: */
