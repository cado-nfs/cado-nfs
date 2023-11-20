#include "cado.h" // IWYU pragma: keep
#include <stdio.h>        // for fprintf, stderr, asprintf
#include <stdlib.h>       // for free, abort
#include <stdarg.h>
#include <errno.h>
#include <string.h>

#include <sys/stat.h>
#include <unistd.h>
#include <pthread.h> // IWYU pragma: keep
#ifdef HAVE_STATVFS_H
#include <sys/statvfs.h>
#endif

#include "bwc_config.h" // BUILD_DYNAMICALLY_LINKABLE_BWC // IWYU pragma: keep

#ifdef  BUILD_DYNAMICALLY_LINKABLE_BWC
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include "solib-naming.h" // IWYU pragma: keep
#endif

#include <gmp.h>          // for mp_bits_per_limb
#include "macros.h"       // for FATAL_ERROR_CHECK, iceildiv
#include "matmul.hpp"
#include "verbose.h"
#include "portability.h" // asprintf // IWYU pragma: keep
#include "params.h"

void matmul_decl_usage(param_list_ptr pl)
{
    param_list_decl_usage(pl, "mm_impl",
            "name of the lower layer matmul implementation");
    param_list_decl_usage(pl, "mm_store_transposed",
            "override the default setting for the matrix storage ordering");

    param_list_decl_usage(pl, "l1_cache_size",
            "internal, for mm_impl=bucket and mm_impl=sliced");
    param_list_decl_usage(pl, "l2_cache_size",
            "internal, for mm_impl=bucket");
    param_list_decl_usage(pl, "cache_line_size",
            "internal");
#if 0
    param_list_decl_usage(pl, "mm_threaded_nthreads",
            "internal, for mm_impl=threaded");
    param_list_decl_usage(pl, "mm_threaded_sgroup_size",
            "internal, for mm_impl=threaded");
    param_list_decl_usage(pl, "mm_threaded_offset1",
            "internal, for mm_impl=threaded");
    param_list_decl_usage(pl, "mm_threaded_offset2",
            "internal, for mm_impl=threaded");
    param_list_decl_usage(pl, "mm_threaded_offset3",
            "internal, for mm_impl=threaded");
    param_list_decl_usage(pl, "mm_threaded_densify_tolerance",
            "internal, for mm_impl=threaded");
#endif

    param_list_decl_usage(pl, "matmul_bucket_methods",
            "internal, for mm_impl=bucket");

    param_list_decl_usage(pl, "local_cache_copy_dir",
            "path to a local directory where a secondary copy of the cache will be saved");
    param_list_decl_usage(pl, "no_save_cache",
            "skip saving the cache file to disk");
}

void matmul_lookup_parameters(param_list_ptr pl)
{
    param_list_lookup_string(pl, "mm_impl");
    param_list_lookup_string(pl, "mm_store_transposed");

    param_list_lookup_string(pl, "l1_cache_size");
    param_list_lookup_string(pl, "l2_cache_size");
    param_list_lookup_string(pl, "cache_line_size");
#if 0
    param_list_lookup_string(pl, "mm_threaded_nthreads");
    param_list_lookup_string(pl, "mm_threaded_sgroup_size");
    param_list_lookup_string(pl, "mm_threaded_offset1");
    param_list_lookup_string(pl, "mm_threaded_offset2");
    param_list_lookup_string(pl, "mm_threaded_offset3");
    param_list_lookup_string(pl, "mm_threaded_densify_tolerance");
#endif
    param_list_lookup_string(pl, "matmul_bucket_methods");

    param_list_lookup_string(pl, "local_cache_copy_dir");
    param_list_lookup_string(pl, "no_save_cache");
}



#ifndef  BUILD_DYNAMICALLY_LINKABLE_BWC
#define MATMUL_NAME_INNER(abase, impl, func) matmul_ ## abase ## _ ## impl ## _ ## func
#define MATMUL_NAME_INNER0(abase, impl, func) MATMUL_NAME_INNER(abase, impl, func)
#define CONFIGURE_MATMUL_LIB(d_, i_)				\
    } else if (strcmp(impl, #i_) == 0 && strcmp(dimpl, #d_) == 0) {		\
        extern void MATMUL_NAME_INNER0(d_, i_, rebind)(matmul_ptr mm);	\
        return MATMUL_NAME_INNER0(d_, i_, rebind);


/* this is a function taking two strings, and returning a pointer to a
 * function taking a matmul_ptr and returning void */
void (*get_rebinder(const char * impl, const char * dimpl))(matmul_ptr mm)
{
    if (0) {
        return NULL;
#define DO(x, y) CONFIGURE_MATMUL_LIB(x, y)
    COOKED_BWC_BACKENDS;
    /* There's a cmake-defined macro in cado_config.h which exposes the
     * configuration settings. Its expansion is something like:
     *
    CONFIGURE_MATMUL_LIB(b64     , bucket)
    CONFIGURE_MATMUL_LIB(b128    , bucket)
    CONFIGURE_MATMUL_LIB(b64     , basic)
    CONFIGURE_MATMUL_LIB(b128    , basic)
    CONFIGURE_MATMUL_LIB(b64     , sliced)
    CONFIGURE_MATMUL_LIB(b128    , sliced)
    CONFIGURE_MATMUL_LIB(b64     , threaded)
    CONFIGURE_MATMUL_LIB(b128    , threaded)
    CONFIGURE_MATMUL_LIB(p1      , basicp)
    CONFIGURE_MATMUL_LIB(p2      , basicp)
    CONFIGURE_MATMUL_LIB(p3      , basicp)
    CONFIGURE_MATMUL_LIB(p4      , basicp)
    CONFIGURE_MATMUL_LIB(p8      , basicp)
    */
    } else {
        fprintf(stderr, "Cannot find the proper rebinder for data backend = %s and matmul backend = %s ; are the corresponding configuration lines present in local.sh or linalg/bwc/CMakeLists.txt ?\n", dimpl, impl);
        return NULL;
    }
}
#endif

matmul_ptr matmul_init(arith_generic * x, unsigned int nr, unsigned int nc, const char * locfile, const char * impl, param_list pl, int optimized_direction)
{
    ASSERT_ALWAYS(pl);

    struct matmul_public_s fake[1];
    memset(fake, 0, sizeof(fake));

    if (!impl) { impl = param_list_lookup_string(pl, "mm_impl"); }
    if (!impl) {
        const char * tmp = param_list_lookup_string(pl, "prime");
        if (tmp && strcmp(tmp, "2") != 0) {
#ifdef HAVE_CXX11
            impl = "zone";
#else   /* HAVE_CXX11 */
            impl = "basicp";
#endif  /* HAVE_CXX11 */
        } else {
            impl = "bucket";
        }
    }

    typedef void (*rebinder_t)(matmul_ptr mm);
    rebinder_t rebinder;

#ifdef  BUILD_DYNAMICALLY_LINKABLE_BWC
    char solib[256];
    snprintf(solib, sizeof(solib),
            SOLIB_PREFIX "matmul_%s_%s" SOLIB_SUFFIX,
            x->impl_name().c_str(), impl);

    static pthread_mutex_t pp = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&pp);
    void * handle = dlopen(solib, RTLD_NOW);
    if (handle == NULL) {
        fprintf(stderr, "loading %s: %s\n", solib, dlerror());
        abort();
    }
    rebinder = (rebinder_t) dlsym(handle, "matmul_solib_do_rebinding");
    if (rebinder == NULL) {
        fprintf(stderr, "loading %s: %s\n", solib, dlerror());
        abort();
    }
    pthread_mutex_unlock(&pp);
#else   /* BUILD_DYNAMICALLY_LINKABLE_BWC */
    std::string impl_name = x->impl_name();
    rebinder = get_rebinder(impl, impl_name.c_str());
#endif   /* BUILD_DYNAMICALLY_LINKABLE_BWC */

    (*rebinder)(fake);
    matmul_ptr mm = fake->bind->init(x->concrete(), pl, optimized_direction);
    if (mm == NULL) return NULL;
#ifdef  BUILD_DYNAMICALLY_LINKABLE_BWC
    mm->solib_handle = handle;
#endif
    (*rebinder)(mm);

    mm->dim[0] = nr;
    mm->dim[1] = nc;

    // this just copies the locfile field from the parent structure
    // matmul_top... Except that for benches, the matmul structure lives
    // outside of this.
    mm->locfile = locfile;
    param_list_parse_int(pl, "no_save_cache", &mm->no_save_cache);

    if (locfile) {
        int rc = asprintf(&mm->cachefile_name, "%s-%s%s.bin", locfile, mm->bind->impl, mm->store_transposed ? "T" : "");
        FATAL_ERROR_CHECK(rc < 0, "out of memory");
    } else {
        mm->cachefile_name = NULL;
    }

    mm->local_cache_copy = NULL;

    const char * local_cache_copy_dir = param_list_lookup_string(pl, "local_cache_copy_dir");
    if (local_cache_copy_dir && mm->cachefile_name) {
        char * basename;
        basename = strrchr(mm->cachefile_name, '/');
        if (basename == NULL) {
            basename = mm->cachefile_name;
        }
        struct stat sbuf[1];
        int rc = stat(local_cache_copy_dir, sbuf);
        if (rc < 0) {
            fprintf(stderr, "Warning: accessing %s is not possible: %s\n",
                    local_cache_copy_dir, strerror(errno));
            return mm;
        }
        rc = asprintf(&mm->local_cache_copy, "%s/%s", local_cache_copy_dir, basename);
        FATAL_ERROR_CHECK(rc < 0, "out of memory");
    }

    return mm;
}

void matmul_build_cache(matmul_ptr mm, matrix_u32_ptr m)
{
    /*
    if (m == NULL) {
        fprintf(stderr, "Called matmul_build_cache() with NULL as matrix_u32_ptr argument. A guaranteed abort() !\n");
    }
    */
    /* read from ->locfile if the raw_matrix_u32 structure has not yet
     * been created. */
    mm->bind->build_cache(mm, m ? m->p : NULL, m ? m->size : 0);
}

static void save_to_local_copy(matmul_ptr mm)
{
    if (mm->no_save_cache ||
        !mm->local_cache_copy ||
        !mm->cachefile_name)
        return;

    struct stat sbuf[1];
    int rc;

    rc = stat(mm->cachefile_name, sbuf);
    if (rc < 0) {
        fprintf(stderr, "stat(%s): %s\n", mm->cachefile_name, strerror(errno));
        return;
    }
    unsigned long fsize = sbuf->st_size;

#ifdef HAVE_STATVFS_H
    /* Check for remaining space on the filesystem */
    char * dirname = strdup(mm->local_cache_copy);
    char * last_slash = strrchr(dirname, '/');
    if (last_slash == NULL) {
        free(dirname);
        dirname = strdup(".");
    } else {
        *last_slash = 0;
    }


    struct statvfs sf[1];
    rc = statvfs(dirname, sf);
    if (rc >= 0) {
        unsigned long mb = sf->f_bsize * sf->f_bavail;

        if (fsize > mb * 0.5) {
            fprintf(stderr, "Copying %s to %s would occupy %lu MB out of %lu MB available, so more than 50%%. Skipping copy\n",
                    mm->cachefile_name, dirname, fsize >> 20, mb >> 20);
            free(dirname);
            return;
        }
        fprintf(stderr, "%lu MB available on %s\n", mb >> 20, dirname);
    } else {
        fprintf(stderr, "Cannot do statvfs on %s (skipping check for available disk space): %s\n", dirname, strerror(errno));
    }
    free(dirname);
#endif


    if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_MAJOR_INFO)) {
        fprintf(stderr, "Also saving cache data to %s (%lu MB)\n", mm->local_cache_copy, fsize >> 20);
    }

    char * normal_cachefile = mm->cachefile_name;
    mm->cachefile_name = mm->local_cache_copy;
    mm->bind->save_cache(mm);
    mm->cachefile_name = normal_cachefile;
}

int matmul_reload_cache(matmul_ptr mm)
{
    struct stat sbuf[2][1];
    int rc;
    if (!mm->cachefile_name)
        return 0;
    rc = stat(mm->cachefile_name, sbuf[0]);
    if (rc < 0) {
        return 0;
    }
    int local_is_ok = mm->local_cache_copy != NULL;
    if (local_is_ok) {
        /* coverity[toctou] */
        rc = stat(mm->local_cache_copy, sbuf[1]);
        local_is_ok = rc == 0;
    }

    if (local_is_ok && (sbuf[0]->st_size != sbuf[1]->st_size)) {
        fprintf(stderr, "%s and %s differ in size ; latter ignored\n",
                mm->cachefile_name,
                mm->local_cache_copy);
        unlink(mm->local_cache_copy);
        local_is_ok = 0;
    }

    if (local_is_ok && (sbuf[0]->st_mtime > sbuf[1]->st_mtime)) {
        fprintf(stderr, "%s is newer than %s ; latter ignored\n",
                mm->cachefile_name,
                mm->local_cache_copy);
        unlink(mm->local_cache_copy);
        local_is_ok = 0;
    }

    if (!local_is_ok) {
        // no local copy.
        rc = mm->bind->reload_cache(mm);
        if (rc == 0)
            return 0;
        // succeeded in loading data.
        save_to_local_copy(mm);
        return 1;
    } else {
        char * normal_cachefile = mm->cachefile_name;
        mm->cachefile_name = mm->local_cache_copy;
        rc = mm->bind->reload_cache(mm);
        mm->cachefile_name = normal_cachefile;
        return rc;
    }
}

void matmul_save_cache(matmul_ptr mm)
{
    if (mm->no_save_cache) return;
    mm->bind->save_cache(mm);
    save_to_local_copy(mm);
}

void matmul_mul(matmul_ptr mm, void * dst, void const * src, int d)
{
    /* d == 1: this is a matrix-times-vector product.
     * d == 0: this is a vector-times-matrix product.
     *
     * In terms of number of rows and columns of vectors, this has
     * implications of course.
     *
     * A matrix times vector product takes a src vector of mm->dim[1]
     * coordinates, and outputs a dst vector of mm->dim[0] coordinates.
     *
     * This external interface is _not_ concerned with the value of the
     * mm->store_transposed flag. That one only relates to how the
     * implementation layer expects the data to be presented when the
     * cache is being built.
     */
    mm->bind->mul(mm, dst, src, d);
}

void matmul_report(matmul_ptr mm, double scale)
{
    mm->bind->report(mm, scale);
}

void matmul_clear(matmul_ptr mm)
{
    if (mm->report_string_size) {
        free(mm->report_string);
        mm->report_string = NULL;
    }
    if (mm->cachefile_name != NULL) free(mm->cachefile_name);
    if (mm->local_cache_copy != NULL) free(mm->local_cache_copy);
    mm->locfile = NULL;
#ifdef BUILD_DYNAMICALLY_LINKABLE_BWC
    void * handle = mm->solib_handle;
    mm->bind->clear(mm);
    dlclose(handle);
#else
    mm->bind->clear(mm);
#endif
}

void matmul_auxv(matmul_ptr mm, int op, va_list ap)
{
    mm->bind->auxv (mm, op, ap);
}

void matmul_aux(matmul_ptr mm, int op, ...)
{
    va_list ap;
    va_start(ap, op);
    matmul_auxv (mm, op, ap);
    va_end(ap);
}
