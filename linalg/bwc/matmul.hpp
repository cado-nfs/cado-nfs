#ifndef MATMUL_HPP_
#define MATMUL_HPP_

#include <stddef.h>
#include <stdarg.h>
#include <stdint.h>

/* This header is common to the different matrix product implementations
 */
#include "params.h"
#include "arith-generic.hpp"
#include "raw_matrix_u32.h"
#include "bwc_config.h" // BUILD_DYNAMICALLY_LINKABLE_BWC // IWYU pragma: keep

struct matmul_public_s; // IWYU pragma: keep

typedef struct matmul_public_s  * matmul_t;
typedef struct matmul_public_s  * matmul_ptr;

typedef matmul_ptr (*matmul_init_t)(void *, param_list pl, int);
typedef void (*matmul_build_cache_t)(matmul_ptr, uint32_t *, size_t);
typedef int (*matmul_reload_cache_t)(matmul_ptr);
typedef void (*matmul_save_cache_t)(matmul_ptr);
typedef void (*matmul_mul_t)(matmul_ptr, void *, const void *, int);
typedef void (*matmul_report_t)(matmul_ptr, double scale);
typedef void (*matmul_clear_t)(matmul_ptr mm);
typedef void (*matmul_auxv_t)(matmul_ptr mm, int op, va_list ap);
typedef void (*matmul_aux_t)(matmul_ptr mm, int op, ...);

struct matmul_bindings_s {
	matmul_build_cache_t	build_cache;
	matmul_reload_cache_t	reload_cache;
	matmul_save_cache_t	save_cache;
	matmul_mul_t		mul;
	matmul_report_t		report;
	matmul_clear_t		clear;
	matmul_init_t		init;
	matmul_auxv_t		auxv;
	matmul_aux_t		aux;
        const char *            impl;
};

struct matmul_public_s {
    /* The fields here must be exposed by all implementations */
    unsigned int dim[2];        /* dim[0] is nrows, dim[1] is ncols. Really. */
                                /* However, data needs not be stored in
                                 * that order */
    uint64_t ncoeffs;

    int store_transposed;       /* For matrix building purposes, this
                                   indicates whether the implementation
                                   layer expects the flat matrix data in
                                   row major order (0) or transposed, in
                                   column major order (1). */
    int iteration[2];           /* [0]: number of vector times matrix products
                                 * [1]: number of matrix times vector products
                                 */

    const char * locfile;
    int no_save_cache;          /* if true, cache file is not saved */

    char * cachefile_name;
    char * local_cache_copy;

    /* Now the virtual method table */
    struct matmul_bindings_s bind[1];

#ifdef  BUILD_DYNAMICALLY_LINKABLE_BWC
    void * solib_handle;        /* passed to dlclose() eventually */
    /* The rest of the implementation-dependent storage comes right after
     * that, in memory. Therefore, only the pointer may be manipulated
     * freely. The datasize field merely indicates the __on-disk__ data
     * size. The in-memory data size is allowed to involve pointers, and
     * is not required to match the on-disk structure exactly.
     */
#endif  /* BUILD_DYNAMICALLY_LINKABLE_BWC */

    char * report_string;
    size_t report_string_size;
};


extern void matmul_decl_usage(param_list_ptr pl);

extern void matmul_lookup_parameters(param_list_ptr pl);

/* This is defined only in the low-level code (e.g. matmul-bucket.cpp) */
extern const char * matmul_mpfq_name();

/* This function is one of the main code multiplexing points. Its get
 * called with a pointer to an abstract virtual base, and the actual
 * instantiation is a derived class that also inherits from a concrete
 * class. Internally, the low-level implementation will only see the
 * concrete class (with some wild type-punning).
 *
 * arith_generic (the vbase) has a virtual void * concrete() const member
 * function that is defined in the derived class, and returns the
 * concrete instance.
 *
 */
extern matmul_ptr matmul_init(arith_generic *, unsigned int, unsigned int, const char * locfile, const char * impl, param_list pl, int);

/* These two functions are the means of initializing the mm layer. No
 * bare _init() function is publicly accessible */
extern void matmul_build_cache(matmul_ptr mm, matrix_u32_ptr m);
extern int matmul_reload_cache(matmul_ptr mm);

/* Exit point */
extern void matmul_clear(matmul_ptr mm);

/* Convenience. Once the data has been built, it may be saved */
extern void matmul_save_cache(matmul_ptr);

/* matmul_mul() is the main entry point for the cpu-intensive critical loop.
 * besides the obvious, there's a direction argument which tells whether
 * we're concerned with the linear operator represented by the matrix
 * (d==1) or its transpose (d==0)
 *
 * Note that not all implementations properly handle both values of d. On
 * paper, both are dual to each other, but doing so equally efficiently
 * with the same in-memory structure is not easy.
 */
extern void matmul_mul(matmul_ptr, void *, const void *, int);

/* matmul_aux is used on some occasions for obtaining private
 * informations on the matmul structure. It's really a virtual method
 * hiding its name.
 *
 * _aux() is only the top-level. auxv is the one that does something.
 */
#define MATMUL_AUX_GET_READAHEAD        10
#define MATMUL_AUX_ZERO_STATS        11
extern void matmul_auxv(matmul_ptr, int, va_list ap);
extern void matmul_aux(matmul_ptr, int, ...);


/* Used for statistics only */
extern void matmul_report(matmul_ptr, double scale);


#endif	/* MATMUL_HPP_ */
