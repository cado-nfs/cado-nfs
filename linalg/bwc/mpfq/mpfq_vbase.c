#include "cado.h"
/* MPFQ generated file -- do not edit */

#include <stdarg.h>
#include "mpfq_vbase.h"
#ifdef COMPILE_MPFQ_BINARY_FIELD_m128
#include "mpfq_m128.h"
#include "mpfq_m128_t.h"
#endif /* COMPILE_MPFQ_BINARY_FIELD_m128 */
#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k1
#include "mpfq_u64k1.h"
#include "mpfq_u64k1_t.h"
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k1 */
#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k2
#include "mpfq_u64k2.h"
#include "mpfq_u64k2_t.h"
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k2 */
#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k3
#include "mpfq_u64k3.h"
#include "mpfq_u64k3_t.h"
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k3 */
#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k4
#include "mpfq_u64k4.h"
#include "mpfq_u64k4_t.h"
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k4 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_1
#include "mpfq_p_1.h"
#include "mpfq_p_1_t.h"
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_1 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_10
#include "mpfq_p_10.h"
#include "mpfq_p_10_t.h"
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_10 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_11
#include "mpfq_p_11.h"
#include "mpfq_p_11_t.h"
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_11 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_12
#include "mpfq_p_12.h"
#include "mpfq_p_12_t.h"
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_12 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_13
#include "mpfq_p_13.h"
#include "mpfq_p_13_t.h"
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_13 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_14
#include "mpfq_p_14.h"
#include "mpfq_p_14_t.h"
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_14 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_15
#include "mpfq_p_15.h"
#include "mpfq_p_15_t.h"
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_15 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_2
#include "mpfq_p_2.h"
#include "mpfq_p_2_t.h"
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_2 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_3
#include "mpfq_p_3.h"
#include "mpfq_p_3_t.h"
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_3 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_4
#include "mpfq_p_4.h"
#include "mpfq_p_4_t.h"
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_4 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_5
#include "mpfq_p_5.h"
#include "mpfq_p_5_t.h"
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_5 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_6
#include "mpfq_p_6.h"
#include "mpfq_p_6_t.h"
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_6 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_7
#include "mpfq_p_7.h"
#include "mpfq_p_7_t.h"
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_7 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_8
#include "mpfq_p_8.h"
#include "mpfq_p_8_t.h"
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_8 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_9
#include "mpfq_p_9.h"
#include "mpfq_p_9_t.h"
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_9 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_pz
#include "mpfq_pz.h"
#include "mpfq_pz_t.h"
#endif /* COMPILE_MPFQ_PRIME_FIELD_pz */
void mpfq_vbase_oo_field_init_byfeatures(mpfq_vbase_ptr v, ...)
{
        va_list ap;
        va_start(ap, v);
        mpz_t p;
        mpz_init_set_ui(p, 2);
        int groupsize = 1;
        const char * mandatory_tag = NULL;
        for(int a ; (a = va_arg(ap, int)) != 0 ; ) {
            if (a == MPFQ_PRIME_MPZ) {
                mpz_set(p, va_arg(ap, mpz_srcptr));
            } else if (a == MPFQ_SIMD_GROUPSIZE) {
                groupsize = va_arg(ap, int);
            } else if (a == MPFQ_MANDATORY_TAG) {
                mandatory_tag = va_arg(ap, const char *);
            } else if (a == MPFQ_PRIME_MPN) {
                fprintf(stderr, "Feature code MPFQ_PRIME_MPN unsupported\n");
                exit(EXIT_FAILURE);
            } else {
                /* We do not support MPFQ_PRIME_MPN. Only MPFQ_PRIME_MPZ*/
                fprintf(stderr, "Feature code %d unsupported\n", a);
                exit(EXIT_FAILURE);
            }
        }
        va_end(ap);
        if (mandatory_tag) {
            if (0) {
#ifdef COMPILE_MPFQ_BINARY_FIELD_m128
            } else if (strcmp(mandatory_tag, "m128") == 0 && groupsize == 128 && mpz_cmp_ui(p, 2) == 0) {
                mpfq_m128_oo_field_init(v);
#endif /* COMPILE_MPFQ_BINARY_FIELD_m128 */
#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k1
            } else if (strcmp(mandatory_tag, "u64k1") == 0 && groupsize == 64 && mpz_cmp_ui(p, 2) == 0) {
                mpfq_u64k1_oo_field_init(v);
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k1 */
#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k2
            } else if (strcmp(mandatory_tag, "u64k2") == 0 && groupsize == 128 && mpz_cmp_ui(p, 2) == 0) {
                mpfq_u64k2_oo_field_init(v);
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k2 */
#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k3
            } else if (strcmp(mandatory_tag, "u64k3") == 0 && groupsize == 192 && mpz_cmp_ui(p, 2) == 0) {
                mpfq_u64k3_oo_field_init(v);
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k3 */
#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k4
            } else if (strcmp(mandatory_tag, "u64k4") == 0 && groupsize == 256 && mpz_cmp_ui(p, 2) == 0) {
                mpfq_u64k4_oo_field_init(v);
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k4 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_1
            } else if (strcmp(mandatory_tag, "p_1") == 0 && groupsize == 1 && mpz_size(p) == 1) {
                mpfq_p_1_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_1 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_10
            } else if (strcmp(mandatory_tag, "p_10") == 0 && groupsize == 1 && mpz_size(p) == 10) {
                mpfq_p_10_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_10 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_11
            } else if (strcmp(mandatory_tag, "p_11") == 0 && groupsize == 1 && mpz_size(p) == 11) {
                mpfq_p_11_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_11 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_12
            } else if (strcmp(mandatory_tag, "p_12") == 0 && groupsize == 1 && mpz_size(p) == 12) {
                mpfq_p_12_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_12 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_13
            } else if (strcmp(mandatory_tag, "p_13") == 0 && groupsize == 1 && mpz_size(p) == 13) {
                mpfq_p_13_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_13 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_14
            } else if (strcmp(mandatory_tag, "p_14") == 0 && groupsize == 1 && mpz_size(p) == 14) {
                mpfq_p_14_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_14 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_15
            } else if (strcmp(mandatory_tag, "p_15") == 0 && groupsize == 1 && mpz_size(p) == 15) {
                mpfq_p_15_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_15 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_2
            } else if (strcmp(mandatory_tag, "p_2") == 0 && groupsize == 1 && mpz_size(p) == 2) {
                mpfq_p_2_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_2 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_3
            } else if (strcmp(mandatory_tag, "p_3") == 0 && groupsize == 1 && mpz_size(p) == 3) {
                mpfq_p_3_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_3 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_4
            } else if (strcmp(mandatory_tag, "p_4") == 0 && groupsize == 1 && mpz_size(p) == 4) {
                mpfq_p_4_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_4 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_5
            } else if (strcmp(mandatory_tag, "p_5") == 0 && groupsize == 1 && mpz_size(p) == 5) {
                mpfq_p_5_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_5 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_6
            } else if (strcmp(mandatory_tag, "p_6") == 0 && groupsize == 1 && mpz_size(p) == 6) {
                mpfq_p_6_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_6 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_7
            } else if (strcmp(mandatory_tag, "p_7") == 0 && groupsize == 1 && mpz_size(p) == 7) {
                mpfq_p_7_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_7 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_8
            } else if (strcmp(mandatory_tag, "p_8") == 0 && groupsize == 1 && mpz_size(p) == 8) {
                mpfq_p_8_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_8 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_9
            } else if (strcmp(mandatory_tag, "p_9") == 0 && groupsize == 1 && mpz_size(p) == 9) {
                mpfq_p_9_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_9 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_pz
            } else if (strcmp(mandatory_tag, "pz") == 0 && groupsize == 1 && 1) {
                mpfq_pz_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_pz */
            } else {
                gmp_fprintf(stderr, "Unsupported combination: mandatory tag = %s, group size = %d, p = %Zd, %zu limbs\n", mandatory_tag, groupsize, p, mpz_size(p));
                exit(1);
            }
        /* now for the case where there is no mandatory tag */
#ifdef COMPILE_MPFQ_BINARY_FIELD_m128
        } else if (groupsize == 128 && mpz_cmp_ui(p, 2) == 0) {
            mpfq_m128_oo_field_init(v);
#endif /* COMPILE_MPFQ_BINARY_FIELD_m128 */
#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k1
        } else if (groupsize == 64 && mpz_cmp_ui(p, 2) == 0) {
            mpfq_u64k1_oo_field_init(v);
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k1 */
#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k2
        } else if (groupsize == 128 && mpz_cmp_ui(p, 2) == 0) {
            mpfq_u64k2_oo_field_init(v);
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k2 */
#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k3
        } else if (groupsize == 192 && mpz_cmp_ui(p, 2) == 0) {
            mpfq_u64k3_oo_field_init(v);
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k3 */
#ifdef COMPILE_MPFQ_BINARY_FIELD_u64k4
        } else if (groupsize == 256 && mpz_cmp_ui(p, 2) == 0) {
            mpfq_u64k4_oo_field_init(v);
#endif /* COMPILE_MPFQ_BINARY_FIELD_u64k4 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_1
        } else if (groupsize == 1 && mpz_size(p) == 1) {
            mpfq_p_1_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_1 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_10
        } else if (groupsize == 1 && mpz_size(p) == 10) {
            mpfq_p_10_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_10 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_11
        } else if (groupsize == 1 && mpz_size(p) == 11) {
            mpfq_p_11_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_11 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_12
        } else if (groupsize == 1 && mpz_size(p) == 12) {
            mpfq_p_12_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_12 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_13
        } else if (groupsize == 1 && mpz_size(p) == 13) {
            mpfq_p_13_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_13 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_14
        } else if (groupsize == 1 && mpz_size(p) == 14) {
            mpfq_p_14_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_14 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_15
        } else if (groupsize == 1 && mpz_size(p) == 15) {
            mpfq_p_15_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_15 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_2
        } else if (groupsize == 1 && mpz_size(p) == 2) {
            mpfq_p_2_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_2 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_3
        } else if (groupsize == 1 && mpz_size(p) == 3) {
            mpfq_p_3_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_3 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_4
        } else if (groupsize == 1 && mpz_size(p) == 4) {
            mpfq_p_4_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_4 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_5
        } else if (groupsize == 1 && mpz_size(p) == 5) {
            mpfq_p_5_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_5 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_6
        } else if (groupsize == 1 && mpz_size(p) == 6) {
            mpfq_p_6_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_6 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_7
        } else if (groupsize == 1 && mpz_size(p) == 7) {
            mpfq_p_7_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_7 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_8
        } else if (groupsize == 1 && mpz_size(p) == 8) {
            mpfq_p_8_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_8 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_p_9
        } else if (groupsize == 1 && mpz_size(p) == 9) {
            mpfq_p_9_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_p_9 */
#ifdef COMPILE_MPFQ_PRIME_FIELD_pz
        } else if (groupsize == 1 && 1) {
            mpfq_pz_oo_field_init(v);
#endif /* COMPILE_MPFQ_PRIME_FIELD_pz */
        } else {
            gmp_fprintf(stderr, "Unsupported combination: group size = %d, p = %Zd, %zu limbs\n", groupsize, p, mpz_size(p));
            exit(1);
        }
        v->field_specify(v, MPFQ_PRIME_MPZ, p);
        v->field_specify(v, MPFQ_SIMD_GROUPSIZE, &groupsize);
        mpz_clear(p);
}

void mpfq_vbase_oo_init_templates(mpfq_vbase_tmpl_ptr w, mpfq_vbase_ptr v0, mpfq_vbase_ptr v1)
{
    const char * s0 = v0->impl_name();
    const char * s1 = v1->impl_name();
    if (0) {
#if defined(COMPILE_MPFQ_BINARY_FIELD_m128)
    } else if (strcmp(s0, "m128") == 0 && strcmp(s1, "m128") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_m128_m128_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_m128_m128_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_m128_m128_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_m128) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_m128) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k1)
    } else if (strcmp(s0, "m128") == 0 && strcmp(s1, "u64k1") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_m128_u64k1_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_m128_u64k1_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_m128_u64k1_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_m128) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k1) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_m128) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k2)
    } else if (strcmp(s0, "m128") == 0 && strcmp(s1, "u64k2") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_m128_u64k2_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_m128_u64k2_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_m128_u64k2_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_m128) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k2) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_m128) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k3)
    } else if (strcmp(s0, "m128") == 0 && strcmp(s1, "u64k3") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_m128_u64k3_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_m128_u64k3_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_m128_u64k3_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_m128) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k3) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_m128) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k4)
    } else if (strcmp(s0, "m128") == 0 && strcmp(s1, "u64k4") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_m128_u64k4_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_m128_u64k4_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_m128_u64k4_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_m128) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k4) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_m128) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k1)
    } else if (strcmp(s0, "u64k1") == 0 && strcmp(s1, "m128") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k1_m128_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k1_m128_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k1_m128_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_m128) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k1) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_u64k1)
    } else if (strcmp(s0, "u64k1") == 0 && strcmp(s1, "u64k1") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k1_u64k1_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k1_u64k1_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k1_u64k1_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_u64k1) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_u64k1) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k2)
    } else if (strcmp(s0, "u64k1") == 0 && strcmp(s1, "u64k2") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k1_u64k2_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k1_u64k2_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k1_u64k2_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_u64k1) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k2) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_u64k1) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k3)
    } else if (strcmp(s0, "u64k1") == 0 && strcmp(s1, "u64k3") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k1_u64k3_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k1_u64k3_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k1_u64k3_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_u64k1) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k3) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_u64k1) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k4)
    } else if (strcmp(s0, "u64k1") == 0 && strcmp(s1, "u64k4") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k1_u64k4_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k1_u64k4_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k1_u64k4_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_u64k1) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k4) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_m128) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k2)
    } else if (strcmp(s0, "u64k2") == 0 && strcmp(s1, "m128") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k2_m128_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k2_m128_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k2_m128_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_m128) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k2) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_u64k1) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k2)
    } else if (strcmp(s0, "u64k2") == 0 && strcmp(s1, "u64k1") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k2_u64k1_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k2_u64k1_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k2_u64k1_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_u64k1) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k2) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_u64k2)
    } else if (strcmp(s0, "u64k2") == 0 && strcmp(s1, "u64k2") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k2_u64k2_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k2_u64k2_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k2_u64k2_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_u64k2) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_u64k2) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k3)
    } else if (strcmp(s0, "u64k2") == 0 && strcmp(s1, "u64k3") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k2_u64k3_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k2_u64k3_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k2_u64k3_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_u64k2) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k3) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_u64k2) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k4)
    } else if (strcmp(s0, "u64k2") == 0 && strcmp(s1, "u64k4") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k2_u64k4_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k2_u64k4_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k2_u64k4_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_u64k2) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k4) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_m128) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k3)
    } else if (strcmp(s0, "u64k3") == 0 && strcmp(s1, "m128") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k3_m128_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k3_m128_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k3_m128_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_m128) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k3) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_u64k1) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k3)
    } else if (strcmp(s0, "u64k3") == 0 && strcmp(s1, "u64k1") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k3_u64k1_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k3_u64k1_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k3_u64k1_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_u64k1) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k3) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_u64k2) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k3)
    } else if (strcmp(s0, "u64k3") == 0 && strcmp(s1, "u64k2") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k3_u64k2_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k3_u64k2_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k3_u64k2_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_u64k2) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k3) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_u64k3)
    } else if (strcmp(s0, "u64k3") == 0 && strcmp(s1, "u64k3") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k3_u64k3_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k3_u64k3_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k3_u64k3_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_u64k3) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_u64k3) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k4)
    } else if (strcmp(s0, "u64k3") == 0 && strcmp(s1, "u64k4") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k3_u64k4_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k3_u64k4_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k3_u64k4_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_u64k3) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k4) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_m128) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k4)
    } else if (strcmp(s0, "u64k4") == 0 && strcmp(s1, "m128") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k4_m128_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k4_m128_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k4_m128_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_m128) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k4) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_u64k1) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k4)
    } else if (strcmp(s0, "u64k4") == 0 && strcmp(s1, "u64k1") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k4_u64k1_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k4_u64k1_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k4_u64k1_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_u64k1) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k4) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_u64k2) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k4)
    } else if (strcmp(s0, "u64k4") == 0 && strcmp(s1, "u64k2") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k4_u64k2_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k4_u64k2_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k4_u64k2_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_u64k2) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k4) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_u64k3) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k4)
    } else if (strcmp(s0, "u64k4") == 0 && strcmp(s1, "u64k3") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k4_u64k3_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k4_u64k3_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k4_u64k3_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_u64k3) && defined(COMPILE_MPFQ_BINARY_FIELD_u64k4) */
#if defined(COMPILE_MPFQ_BINARY_FIELD_u64k4)
    } else if (strcmp(s0, "u64k4") == 0 && strcmp(s1, "u64k4") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k4_u64k4_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_u64k4_u64k4_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_u64k4_u64k4_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_BINARY_FIELD_u64k4) */
#if defined(COMPILE_MPFQ_PRIME_FIELD_p_1)
    } else if (strcmp(s0, "p_1") == 0 && strcmp(s1, "p_1") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_1_p_1_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_1_p_1_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_p_1_p_1_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_PRIME_FIELD_p_1) */
#if defined(COMPILE_MPFQ_PRIME_FIELD_p_10)
    } else if (strcmp(s0, "p_10") == 0 && strcmp(s1, "p_10") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_10_p_10_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_10_p_10_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_p_10_p_10_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_PRIME_FIELD_p_10) */
#if defined(COMPILE_MPFQ_PRIME_FIELD_p_11)
    } else if (strcmp(s0, "p_11") == 0 && strcmp(s1, "p_11") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_11_p_11_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_11_p_11_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_p_11_p_11_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_PRIME_FIELD_p_11) */
#if defined(COMPILE_MPFQ_PRIME_FIELD_p_12)
    } else if (strcmp(s0, "p_12") == 0 && strcmp(s1, "p_12") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_12_p_12_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_12_p_12_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_p_12_p_12_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_PRIME_FIELD_p_12) */
#if defined(COMPILE_MPFQ_PRIME_FIELD_p_13)
    } else if (strcmp(s0, "p_13") == 0 && strcmp(s1, "p_13") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_13_p_13_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_13_p_13_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_p_13_p_13_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_PRIME_FIELD_p_13) */
#if defined(COMPILE_MPFQ_PRIME_FIELD_p_14)
    } else if (strcmp(s0, "p_14") == 0 && strcmp(s1, "p_14") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_14_p_14_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_14_p_14_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_p_14_p_14_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_PRIME_FIELD_p_14) */
#if defined(COMPILE_MPFQ_PRIME_FIELD_p_15)
    } else if (strcmp(s0, "p_15") == 0 && strcmp(s1, "p_15") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_15_p_15_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_15_p_15_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_p_15_p_15_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_PRIME_FIELD_p_15) */
#if defined(COMPILE_MPFQ_PRIME_FIELD_p_2)
    } else if (strcmp(s0, "p_2") == 0 && strcmp(s1, "p_2") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_2_p_2_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_2_p_2_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_p_2_p_2_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_PRIME_FIELD_p_2) */
#if defined(COMPILE_MPFQ_PRIME_FIELD_p_3)
    } else if (strcmp(s0, "p_3") == 0 && strcmp(s1, "p_3") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_3_p_3_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_3_p_3_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_p_3_p_3_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_PRIME_FIELD_p_3) */
#if defined(COMPILE_MPFQ_PRIME_FIELD_p_4)
    } else if (strcmp(s0, "p_4") == 0 && strcmp(s1, "p_4") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_4_p_4_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_4_p_4_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_p_4_p_4_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_PRIME_FIELD_p_4) */
#if defined(COMPILE_MPFQ_PRIME_FIELD_p_5)
    } else if (strcmp(s0, "p_5") == 0 && strcmp(s1, "p_5") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_5_p_5_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_5_p_5_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_p_5_p_5_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_PRIME_FIELD_p_5) */
#if defined(COMPILE_MPFQ_PRIME_FIELD_p_6)
    } else if (strcmp(s0, "p_6") == 0 && strcmp(s1, "p_6") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_6_p_6_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_6_p_6_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_p_6_p_6_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_PRIME_FIELD_p_6) */
#if defined(COMPILE_MPFQ_PRIME_FIELD_p_7)
    } else if (strcmp(s0, "p_7") == 0 && strcmp(s1, "p_7") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_7_p_7_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_7_p_7_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_p_7_p_7_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_PRIME_FIELD_p_7) */
#if defined(COMPILE_MPFQ_PRIME_FIELD_p_8)
    } else if (strcmp(s0, "p_8") == 0 && strcmp(s1, "p_8") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_8_p_8_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_8_p_8_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_p_8_p_8_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_PRIME_FIELD_p_8) */
#if defined(COMPILE_MPFQ_PRIME_FIELD_p_9)
    } else if (strcmp(s0, "p_9") == 0 && strcmp(s1, "p_9") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_9_p_9_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_p_9_p_9_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_p_9_p_9_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_PRIME_FIELD_p_9) */
#if defined(COMPILE_MPFQ_PRIME_FIELD_pz)
    } else if (strcmp(s0, "pz") == 0 && strcmp(s1, "pz") == 0) {
        w->add_dotprod = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_pz_pz_wrapper_add_dotprod;
        w->addmul_tiny = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int)) mpfq_pz_pz_wrapper_addmul_tiny;
        w->transpose = (void (*) (mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *)) mpfq_pz_pz_wrapper_transpose;
#endif /* defined(COMPILE_MPFQ_PRIME_FIELD_pz) */
    } else {
        abort();
    }
}


/* vim:set ft=cpp: */
