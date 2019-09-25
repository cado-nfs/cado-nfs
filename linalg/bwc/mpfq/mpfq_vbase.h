#ifndef MPFQ_VBASE_H_
#define MPFQ_VBASE_H_

/* MPFQ generated file -- do not edit */

#include <stddef.h>
#include <stdio.h>
#include <gmp.h>

#ifdef __cplusplus
extern "C" {
#endif

struct mpfq_vbase_s;
typedef struct mpfq_vbase_s * mpfq_vbase_ptr;
typedef struct mpfq_vbase_s const * mpfq_vbase_srcptr;

struct mpfq_vbase_tmpl_s;
typedef struct mpfq_vbase_tmpl_s * mpfq_vbase_tmpl_ptr;
typedef struct mpfq_vbase_tmpl_s const * mpfq_vbase_tmpl_srcptr;

struct mpfq_vbase_s {
    void * obj; /* pointer to global implementation private fields */
    const char * (*impl_name)();
    unsigned long (*impl_max_characteristic_bits)();
    unsigned long (*impl_max_degree)();
    void (*field_characteristic)(mpfq_vbase_ptr, mpz_ptr);
    mpz_srcptr (*field_characteristic_srcptr)(mpfq_vbase_ptr);
    unsigned long (*field_characteristic_bits)(mpfq_vbase_ptr);
    int (*field_degree)(mpfq_vbase_ptr);
    void (*field_init)(mpfq_vbase_ptr);
    void (*field_clear)(mpfq_vbase_ptr);
    void (*field_specify)(mpfq_vbase_ptr, unsigned long, const void *);
    void (*field_setopt)(mpfq_vbase_ptr, unsigned long, void *);
    void (*init)(mpfq_vbase_ptr, void *);
    void (*clear)(mpfq_vbase_ptr, void *);
    ptrdiff_t (*elt_stride)(mpfq_vbase_ptr);
    void (*set)(mpfq_vbase_ptr, void *, const void *);
    void (*set_ui)(mpfq_vbase_ptr, void *, unsigned long);
    void (*set_zero)(mpfq_vbase_ptr, void *);
    unsigned long (*get_ui)(mpfq_vbase_ptr, const void *);
    void (*set_mpn)(mpfq_vbase_ptr, void *, mp_limb_t *, size_t);
    void (*set_mpz)(mpfq_vbase_ptr, void *, mpz_srcptr);
    void (*get_mpn)(mpfq_vbase_ptr, mp_limb_t *, const void *);
    void (*get_mpz)(mpfq_vbase_ptr, mpz_ptr, const void *);
    void (*random)(mpfq_vbase_ptr, void *, gmp_randstate_t);
    void (*random2)(mpfq_vbase_ptr, void *, gmp_randstate_t);
    void (*add)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*sub)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*neg)(mpfq_vbase_ptr, void *, const void *);
    void (*mul)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*sqr)(mpfq_vbase_ptr, void *, const void *);
    int (*is_sqr)(mpfq_vbase_ptr, const void *);
    int (*sqrt)(mpfq_vbase_ptr, void *, const void *);
    void (*pow)(mpfq_vbase_ptr, void *, const void *, unsigned long *, size_t);
    void (*powz)(mpfq_vbase_ptr, void *, const void *, mpz_srcptr);
    void (*frobenius)(mpfq_vbase_ptr, void *, const void *);
    void (*add_ui)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*sub_ui)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*mul_ui)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*normalize)(mpfq_vbase_ptr, void *);
    int (*inv)(mpfq_vbase_ptr, void *, const void *);
    void (*hadamard)(mpfq_vbase_ptr, void *, void *, void *, void *);
    void (*elt_ur_init)(mpfq_vbase_ptr, void *);
    void (*elt_ur_clear)(mpfq_vbase_ptr, void *);
    ptrdiff_t (*elt_ur_stride)(mpfq_vbase_ptr);
    void (*elt_ur_set)(mpfq_vbase_ptr, void *, const void *);
    void (*elt_ur_set_elt)(mpfq_vbase_ptr, void *, const void *);
    void (*elt_ur_set_zero)(mpfq_vbase_ptr, void *);
    void (*elt_ur_set_ui)(mpfq_vbase_ptr, void *, unsigned long);
    void (*elt_ur_add)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*elt_ur_neg)(mpfq_vbase_ptr, void *, const void *);
    void (*elt_ur_sub)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*mul_ur)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*sqr_ur)(mpfq_vbase_ptr, void *, const void *);
    void (*reduce)(mpfq_vbase_ptr, void *, void *);
    void (*addmul_si_ur)(mpfq_vbase_ptr, void *, const void *, long);
    int (*cmp)(mpfq_vbase_ptr, const void *, const void *);
    int (*cmp_ui)(mpfq_vbase_ptr, const void *, unsigned long);
    int (*is_zero)(mpfq_vbase_ptr, const void *);
    int (*asprint)(mpfq_vbase_ptr, char * *, const void *);
    int (*fprint)(mpfq_vbase_ptr, FILE *, const void *);
    int (*print)(mpfq_vbase_ptr, const void *);
    int (*sscan)(mpfq_vbase_ptr, void *, const char *);
    int (*fscan)(mpfq_vbase_ptr, FILE *, void *);
    int (*scan)(mpfq_vbase_ptr, void *);
    int (*read)(mpfq_vbase_ptr, FILE *, void *);
    int (*importdata)(mpfq_vbase_ptr, FILE *, void *, int, int);
    int (*write)(mpfq_vbase_ptr, FILE *, const void *);
    int (*exportdata)(mpfq_vbase_ptr, FILE *, const void *, int, int);
    void (*vec_init)(mpfq_vbase_ptr, void *, unsigned long);
    void (*vec_reinit)(mpfq_vbase_ptr, void *, unsigned long, unsigned long);
    void (*vec_clear)(mpfq_vbase_ptr, void *, unsigned long);
    void (*vec_set)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*vec_set_zero)(mpfq_vbase_ptr, void *, unsigned long);
    void (*vec_setcoeff)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*vec_setcoeff_ui)(mpfq_vbase_ptr, void *, unsigned long, unsigned long);
    void (*vec_getcoeff)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*vec_add)(mpfq_vbase_ptr, void *, const void *, const void *, unsigned long);
    void (*vec_neg)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*vec_rev)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*vec_sub)(mpfq_vbase_ptr, void *, const void *, const void *, unsigned long);
    void (*vec_scal_mul)(mpfq_vbase_ptr, void *, const void *, const void *, unsigned long);
    void (*vec_conv)(mpfq_vbase_ptr, void *, const void *, unsigned long, const void *, unsigned long);
    void (*vec_random)(mpfq_vbase_ptr, void *, unsigned long, gmp_randstate_t);
    void (*vec_random2)(mpfq_vbase_ptr, void *, unsigned long, gmp_randstate_t);
    long (*vec_cmp)(mpfq_vbase_ptr, const void *, const void *, unsigned long);
    long (*vec_is_zero)(mpfq_vbase_ptr, const void *, unsigned long);
    void * (*vec_subvec)(mpfq_vbase_ptr, void *, long);
    const void * (*vec_subvec_const)(mpfq_vbase_ptr, const void *, long);
    void * (*vec_coeff_ptr)(mpfq_vbase_ptr, void *, long);
    const void * (*vec_coeff_ptr_const)(mpfq_vbase_ptr, const void *, long);
    long (*vec_asprint)(mpfq_vbase_ptr, char * *, const void *, unsigned long);
    long (*vec_fprint)(mpfq_vbase_ptr, FILE *, const void *, unsigned long);
    long (*vec_print)(mpfq_vbase_ptr, const void *, unsigned long);
    long (*vec_sscan)(mpfq_vbase_ptr, void *, unsigned long *, const char *);
    long (*vec_fscan)(mpfq_vbase_ptr, FILE *, void *, unsigned long *);
    long (*vec_scan)(mpfq_vbase_ptr, void *, unsigned long *);
    long (*vec_read)(mpfq_vbase_ptr, FILE *, void *, unsigned long);
    long (*vec_write)(mpfq_vbase_ptr, FILE *, const void *, unsigned long);
    long (*vec_import)(mpfq_vbase_ptr, void *, unsigned long, long, long);
    long (*vec_export)(mpfq_vbase_ptr, const void *, unsigned long, long, long);
    long (*vec_hamming_weight)(mpfq_vbase_ptr, const void *, unsigned long);
    long (*vec_find_first_set)(mpfq_vbase_ptr, const void *, unsigned long);
    long (*vec_simd_hamming_weight)(mpfq_vbase_ptr, const void *, unsigned long);
    long (*vec_simd_find_first_set)(mpfq_vbase_ptr, const void *, unsigned long);
    void (*vec_ur_init)(mpfq_vbase_ptr, void *, unsigned long);
    void (*vec_ur_set_zero)(mpfq_vbase_ptr, void *, unsigned long);
    void (*vec_ur_set_vec)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*vec_ur_reinit)(mpfq_vbase_ptr, void *, unsigned long, unsigned long);
    void (*vec_ur_clear)(mpfq_vbase_ptr, void *, unsigned long);
    void (*vec_ur_set)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*vec_ur_setcoeff)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*vec_ur_getcoeff)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*vec_ur_add)(mpfq_vbase_ptr, void *, const void *, const void *, unsigned long);
    void (*vec_ur_sub)(mpfq_vbase_ptr, void *, const void *, const void *, unsigned long);
    void (*vec_ur_neg)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*vec_ur_rev)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*vec_scal_mul_ur)(mpfq_vbase_ptr, void *, const void *, const void *, unsigned long);
    void (*vec_conv_ur)(mpfq_vbase_ptr, void *, const void *, unsigned long, const void *, unsigned long);
    void (*vec_reduce)(mpfq_vbase_ptr, void *, void *, unsigned long);
    void * (*vec_ur_subvec)(mpfq_vbase_ptr, void *, long);
    const void * (*vec_ur_subvec_const)(mpfq_vbase_ptr, const void *, long);
    void * (*vec_ur_coeff_ptr)(mpfq_vbase_ptr, void *, long);
    const void * (*vec_ur_coeff_ptr_const)(mpfq_vbase_ptr, const void *, long);
    ptrdiff_t (*vec_elt_stride)(mpfq_vbase_ptr, long);
    ptrdiff_t (*vec_ur_elt_stride)(mpfq_vbase_ptr, long);
    void (*poly_init)(mpfq_vbase_ptr, void *, unsigned int);
    void (*poly_clear)(mpfq_vbase_ptr, void *);
    void (*poly_set)(mpfq_vbase_ptr, void *, const void *);
    void (*poly_setmonic)(mpfq_vbase_ptr, void *, const void *);
    void (*poly_setcoeff)(mpfq_vbase_ptr, void *, const void *, unsigned int);
    void (*poly_setcoeff_ui)(mpfq_vbase_ptr, void *, unsigned long, unsigned int);
    void (*poly_getcoeff)(mpfq_vbase_ptr, void *, const void *, unsigned int);
    int (*poly_deg)(mpfq_vbase_ptr, const void *);
    void (*poly_add)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*poly_sub)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*poly_set_ui)(mpfq_vbase_ptr, void *, unsigned long);
    void (*poly_add_ui)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*poly_sub_ui)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*poly_neg)(mpfq_vbase_ptr, void *, const void *);
    void (*poly_scal_mul)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*poly_mul)(mpfq_vbase_ptr, void *, const void *, const void *);
    int (*poly_divmod)(mpfq_vbase_ptr, void *, void *, const void *, const void *);
    void (*poly_precomp_mod)(mpfq_vbase_ptr, void *, const void *);
    void (*poly_mod_pre)(mpfq_vbase_ptr, void *, const void *, const void *, const void *);
    void (*poly_gcd)(mpfq_vbase_ptr, void *, const void *, const void *);
    void (*poly_xgcd)(mpfq_vbase_ptr, void *, void *, void *, const void *, const void *);
    void (*poly_random)(mpfq_vbase_ptr, void *, unsigned int, gmp_randstate_t);
    void (*poly_random2)(mpfq_vbase_ptr, void *, unsigned int, gmp_randstate_t);
    int (*poly_cmp)(mpfq_vbase_ptr, const void *, const void *);
    int (*poly_asprint)(mpfq_vbase_ptr, char * *, const void *);
    int (*poly_fprint)(mpfq_vbase_ptr, FILE *, const void *);
    int (*poly_print)(mpfq_vbase_ptr, const void *);
    int (*poly_sscan)(mpfq_vbase_ptr, void *, const char *);
    int (*poly_fscan)(mpfq_vbase_ptr, FILE *, void *);
    int (*poly_scan)(mpfq_vbase_ptr, void *);
    int (*simd_groupsize)(mpfq_vbase_ptr);
    int (*simd_hamming_weight)(mpfq_vbase_ptr, const void *);
    int (*simd_find_first_set)(mpfq_vbase_ptr, const void *);
    unsigned long (*simd_get_ui_at)(mpfq_vbase_ptr, const void *, int);
    void (*simd_set_ui_at)(mpfq_vbase_ptr, void *, int, unsigned long);
    void (*simd_add_ui_at)(mpfq_vbase_ptr, void *, const void *, int, unsigned long);
    void (*simd_set_ui_all)(mpfq_vbase_ptr, void *, unsigned long);
    void (*add_dotprod)(mpfq_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*mul_constant_ui)(mpfq_vbase_ptr, void *, const void *, unsigned long);
    void (*oo_field_init)(mpfq_vbase_ptr);
    void (*oo_field_clear)(mpfq_vbase_ptr);
};

struct mpfq_vbase_tmpl_s {
    void (*add_dotprod)(mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*addmul_tiny)(mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *, const void *, unsigned int);
    void (*transpose)(mpfq_vbase_ptr, mpfq_vbase_ptr, void *, const void *);
};
typedef struct mpfq_vbase_s mpfq_vbase[1];
typedef struct mpfq_vbase_tmpl_s mpfq_vbase_tmpl[1];

void mpfq_vbase_oo_field_init_byfeatures(mpfq_vbase_ptr, ...);
void mpfq_vbase_oo_init_templates(mpfq_vbase_tmpl_ptr, mpfq_vbase_ptr, mpfq_vbase_ptr);

#ifdef __cplusplus
}
#endif

#endif  /* MPFQ_VBASE_H_ */

/* vim:set ft=cpp: */
