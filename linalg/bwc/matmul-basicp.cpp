#include "cado.h" // IWYU pragma: keep

#include <cstdarg>         // for va_list, va_end, va_start
#include <cstddef>         // for ptrdiff_t
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <climits>
#include <gmp.h>
#include "matmul.h"
#include "matmul-common.h"
#include "mpfq_layer.h"

#include "matmul_facade.h"

#include "arith-modp.hpp"
#include "macros.h"
#include "params.h"

/* define "gfp" as being our c++ type built from the number of words in
 * the underlying mpfq data type.
 *
 * It's a bit hacky, but mpfq should provide a ``number of words''
 * implementation info.
 */
template<typename F, unsigned int m> struct our_gfp_type {
    static const int mpfq_base_field_width = sizeof(F)/sizeof(unsigned long);
    typedef arith_modp::gfp<mpfq_base_field_width> type;
};

template<typename F> struct our_gfp_type<F, UINT_MAX> {
    static const int mpfq_base_field_width = 0;
    /* We *intentionally* do not provide a variable-width GF(p) type with
     * the C++ code. That wouldn't be totally impossible, but I can assure
     * that it would be a royal pain (something like a 200+-line patch of
     * barely parseable c++ hacks to arith-modp.hpp -- tried it out and
     * gave up...).
     */
};

/* If this line complains that ::type is not a type, then see above */
typedef our_gfp_type<abelt,abimpl_max_characteristic_bits()>::type gfp;



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

struct matmul_basicp_data_s {
    /* repeat the fields from the public interface */
    struct matmul_public_s public_[1];
    /* now our private fields */
    size_t datasize;
    abdst_field xab;
    uint32_t * q;
};

void MATMUL_NAME(clear)(matmul_ptr mm0)
{
    struct matmul_basicp_data_s * mm = (struct matmul_basicp_data_s *) mm0;
    matmul_common_clear(mm->public_);
    free(mm->q);
    free(mm);
}

matmul_ptr MATMUL_NAME(init)(void* xx, param_list pl, int optimized_direction)
{
    struct matmul_basicp_data_s * mm;
    mm = (struct matmul_basicp_data_s *) malloc(sizeof(struct matmul_basicp_data_s));
    memset(mm, 0, sizeof(struct matmul_basicp_data_s));
    mm->xab = (abdst_field) xx;

    int suggest = optimized_direction ^ MM_DIR0_PREFERS_TRANSP_MULT;
    mm->public_->store_transposed = suggest;
    if (pl) {
        param_list_parse_int(pl, "mm_store_transposed", 
                &mm->public_->store_transposed);
        if (mm->public_->store_transposed != suggest) {
            fprintf(stderr, "Warning, mm_store_transposed"
                    " overrides suggested matrix storage ordering\n");
        }           
    }   

    return (matmul_ptr) mm;
}

void MATMUL_NAME(build_cache)(matmul_ptr mm0, uint32_t * data, size_t size)
{
    ASSERT_ALWAYS(data);

    struct matmul_basicp_data_s * mm = (struct matmul_basicp_data_s *) mm0;
    unsigned int nrows_t = mm->public_->dim[ mm->public_->store_transposed];
    
    uint32_t * ptr = data;
    unsigned int i = 0;

    /* count coefficients */
    for( ; i < nrows_t ; i++) {
        unsigned int weight = 0;
        weight += *ptr;
        mm->public_->ncoeffs += weight;

        ptr++;
        ptr += weight;
        ptr += weight;
    }

    mm->q = data;

    mm->datasize = nrows_t + 2 * mm->public_->ncoeffs;

    ASSERT_ALWAYS(size == mm->datasize);
    ASSERT_ALWAYS(ptr - data == (ptrdiff_t) mm->datasize);
}

int MATMUL_NAME(reload_cache)(matmul_ptr mm0)
{
    FILE * f;
    struct matmul_basicp_data_s * mm = (struct matmul_basicp_data_s *) mm0;
    f = matmul_common_reload_cache_fopen(sizeof(abelt), mm->public_, MM_MAGIC);
    if (!f) return 0;

    MATMUL_COMMON_READ_ONE32(mm->datasize, f);
    mm->q = (uint32_t *) malloc(mm->datasize * sizeof(uint32_t));
    MATMUL_COMMON_READ_MANY32(mm->q, mm->datasize, f);
    fclose(f);

    return 1;
}

void MATMUL_NAME(save_cache)(matmul_ptr mm0)
{
    FILE * f;

    struct matmul_basicp_data_s * mm = (struct matmul_basicp_data_s *) mm0;
    f = matmul_common_save_cache_fopen(sizeof(abelt), mm->public_, MM_MAGIC);
    if (!f) return;

    MATMUL_COMMON_WRITE_ONE32(mm->datasize, f);
    MATMUL_COMMON_WRITE_MANY32(mm->q, mm->datasize, f);

    fclose(f);
}

void MATMUL_NAME(mul)(matmul_ptr mm0, void * xdst, void const * xsrc, int d)
{
    struct matmul_basicp_data_s * mm = (struct matmul_basicp_data_s *) mm0;
    ASM_COMMENT("multiplication code");
    uint32_t * q = mm->q;
    abdst_field x = mm->xab;
    const gfp::elt * src = (const gfp::elt *) xsrc;
    gfp::elt * dst = (gfp::elt *) xdst;

    gfp::preinv preinverse;
    gfp::elt prime;

    {
        mpz_t p;
        mpz_init(p);
        abfield_characteristic(x, p);
        prime = p;
        mpz_clear(p);
    }

    gfp::compute_preinv(preinverse, prime);

    /* d == 1: matrix times vector product */
    /* d == 0: vector times matrix product */

    /* However the matrix may be stored either row-major
     * (store_transposed == 0) or column-major (store_transposed == 1)
     */

    gfp::elt::zero(dst, mm->public_->dim[!d]);

    if (d == !mm->public_->store_transposed) {
        gfp::elt_ur rowsum;
        ASM_COMMENT("critical loop");
        for(unsigned int i = 0 ; i < mm->public_->dim[!d] ; i++) {
            uint32_t len = *q++;
            unsigned int j = 0;
            rowsum.zero();
            for( ; len-- ; ) {
                j = *q++;
                int32_t c = *(int32_t*)q++;
                ASSERT(j < mm->public_->dim[d]);
                if (c == 1) {
                    gfp::add(rowsum, src[j]);
                } else if (c == -1) {
                    gfp::sub(rowsum, src[j]);
                } else if (c > 0) {
                    gfp::addmul_ui(rowsum, src[j], c, prime, preinverse);
                } else {
                    gfp::submul_ui(rowsum, src[j], -c, prime, preinverse);
                }
            }
            gfp::reduce(dst[i], rowsum, prime, preinverse);
        }
        ASM_COMMENT("end of critical loop");
    } else {
        gfp::elt_ur * tdst = new gfp::elt_ur[mm->public_->dim[!d]];
        // gfp::elt::zero(tdst, mm->public_->dim[!d]);
        if (mm->public_->iteration[d] == 10) {
            fprintf(stderr, "Warning: Doing many iterations with transposed code (not a huge problem for impl=basicp)\n");
        }
        ASM_COMMENT("critical loop (transposed mult)");
        for(unsigned int i = 0 ; i < mm->public_->dim[d] ; i++) {
            uint32_t len = *q++;
            unsigned int j = 0;
            for( ; len-- ; ) {
                j = *q++;
                int32_t c = *(int32_t*)q++;
                ASSERT(j < mm->public_->dim[!d]);
                if (c == 1) {
                    gfp::add(tdst[j], src[i]);
                } else if (c == -1) {
                    gfp::sub(tdst[j], src[i]);
                } else if (c > 0) {
                    gfp::addmul_ui(tdst[j], src[i], c, prime, preinverse);
                } else {
                    gfp::submul_ui(tdst[j], src[i], -c, prime, preinverse);
                }
            }
        }
        for(unsigned int j = 0 ; j < mm->public_->dim[!d] ; j++) {
            gfp::reduce(dst[j], tdst[j], prime, preinverse);
        }
        ASM_COMMENT("end of critical loop (transposed mult)");
        delete[] tdst;
    }
    ASM_COMMENT("end of multiplication code");

    mm->public_->iteration[d]++;
}

void MATMUL_NAME(report)(matmul_ptr mm0 MAYBE_UNUSED, double scale MAYBE_UNUSED) {
}

void MATMUL_NAME(auxv)(matmul_ptr mm0 MAYBE_UNUSED, int op MAYBE_UNUSED, va_list ap MAYBE_UNUSED)
{
}

void MATMUL_NAME(aux)(matmul_ptr mm0, int op, ...)
{
    va_list ap;
    va_start(ap, op);
    MATMUL_NAME(auxv) (mm0, op, ap);
    va_end(ap);
}

/* vim: set sw=4: */
