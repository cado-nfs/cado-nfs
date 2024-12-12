#include "cado.h" // IWYU pragma: keep

#include <cstdio>        // for fprintf, stderr, asprintf
#include <cstdlib>       // for free, abort
#include <cstdarg>
#include <cerrno>
#include <cstring>

#include <sys/stat.h>
#include <unistd.h>
#include <pthread.h> // IWYU pragma: keep
#ifdef HAVE_STATVFS_H
#include <sys/statvfs.h>
#endif

#include "bwc_config.h" // BUILD_DYNAMICALLY_LINKABLE_BWC // IWYU pragma: keep

#ifdef BUILD_DYNAMICALLY_LINKABLE_BWC
#include <dlfcn.h>
#include "solib-naming.h" // IWYU pragma: keep
#endif

#include <gmp.h>          // for mp_bits_per_limb

#include "macros.h"       // for FATAL_ERROR_CHECK, iceildiv
#include "matmul.hpp"
#include "verbose.h"
#include "portability.h" // asprintf // IWYU pragma: keep
#include "params.h"

void matmul_decl_usage(cxx_param_list & pl)
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

void matmul_lookup_parameters(cxx_param_list & pl)
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
#define CONFIGURE_MATMUL_LIB(d_, i_)			\
    } else if (impl == #i_ && dimpl == #d_) {		\
        extern matmul_interface_ctor_t CADO_CONCATENATE4(new_matmul_, d_, _, i_);       \
        return &CADO_CONCATENATE4(new_matmul_, d_, _, i_);


/* this is a function taking two strings, and returning a pointer to a
 * function taking a matmul_ptr and returning void */
static matmul_interface_ctor_t * reach_ctor(std::string const & impl, std::string const & dimpl)
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
        fmt::print(stderr, "Cannot find the proper rebinder for data backend = {} and matmul backend = {} ; are the corresponding configuration lines present in local.sh or linalg/bwc/CMakeLists.txt ?\n", dimpl, impl);
        return NULL;
    }
}
#endif

std::shared_ptr<matmul_interface> matmul_interface::create(arith_generic * x,
        unsigned int nr,
        unsigned int nc,
        std::string const & locfile,
        std::string const & impl,
        cxx_param_list & pl,
        int optimized_direction)
{
    if (impl.empty()) {
        const char * tmp = param_list_lookup_string(pl, "mm_impl");

        if (tmp) return create(x, nr, nc, locfile, tmp, pl, optimized_direction);

        tmp = x->is_characteristic_two() ? "bucket" : "zone";

        return create(x, nr, nc, locfile, tmp, pl, optimized_direction);
    }

    matmul_public stem;

    // there's not much more that we can init, but it's already
    // something...
    stem.dim[0] = nr;
    stem.dim[1] = nc;
    stem.locfile = locfile;
    param_list_parse(pl, "no_save_cache", stem.no_save_cache);

    std::shared_ptr<matmul_interface> mm;

    matmul_interface_ctor_t * ctor;

#ifdef  BUILD_DYNAMICALLY_LINKABLE_BWC
    {
        void * handle;
        std::string solib = fmt::format(
                SOLIB_PREFIX "matmul_{}_{}" SOLIB_SUFFIX,
                x->impl_name(), impl);
        static std::mutex pp;
        std::lock_guard<std::mutex> const dummy(pp);
        handle = dlopen(solib.c_str(), RTLD_NOW);
        if (handle == nullptr) {
            fmt::print(stderr, "loading {}: {}\n", solib, dlerror());
            abort();
        }
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
        auto ctor_getter = (matmul_interface_ctor_t *(*)()) dlsym(handle, "matmul_solib_reach_ctor");
        ctor = (*ctor_getter)();
        // stem.solib_handle.reset(handle);
        if (ctor == nullptr) {
            fmt::print(stderr, "loading {}: {}\n", solib, dlerror());
            abort();
        }
        /* We're using the fact that shared_ptr virtualizes the dtor,
         * which is quite helpful, here.
         */
        mm = std::unique_ptr<matmul_interface, remove_shared_lib_after_object> {
                (*ctor)(std::move(stem), x, pl, optimized_direction),
                { handle },
        };

    }
#else   /* BUILD_DYNAMICALLY_LINKABLE_BWC */
    {
        /* the naive deleter works well in that case. We're still
         * interested in having a shared_ptr for other reasons
         * (interleaving), but beyond that no custom deletion is
         * necessary. */
        ctor = reach_ctor(impl, x->impl_name());
        ASSERT_ALWAYS(ctor != nullptr);
        mm = std::unique_ptr<matmul_interface> {
                (*ctor)(std::move(stem), x, pl, optimized_direction)
        };
    }
#endif   /* BUILD_DYNAMICALLY_LINKABLE_BWC */

    if (!mm) return mm;

    // NOTE that store_transposed is a priori set up by the
    // implementation.

    if (!locfile.empty()) {
        mm->cachefile_name = fmt::format("{}-{}{}.bin",
                locfile,
                impl,
                mm->store_transposed ? "T" : "");
    }

    const char * local_cache_copy_dir = param_list_lookup_string(pl, "local_cache_copy_dir");
    if (local_cache_copy_dir && !mm->cachefile_name.empty()) {
        struct stat sbuf[1];
        int const rc = stat(local_cache_copy_dir, sbuf);
        if (rc < 0) {
            fprintf(stderr, "Warning: accessing %s is not possible: %s\n",
                    local_cache_copy_dir, strerror(errno));
            return mm;
        }

        auto it = mm->cachefile_name.rfind('/');
        it = (it == std::string::npos) ? 0 : (it + 1);

        mm->local_cache_copy = fmt::format("{}/{}",
                local_cache_copy_dir,
                mm->cachefile_name.substr(it));
    }

    return mm;
}

void matmul_interface::save_to_local_copy()
{
    if (no_save_cache ||
        local_cache_copy.empty() ||
        cachefile_name.empty())
        return;

    struct stat sbuf[1];
    int rc;

    rc = stat(cachefile_name.c_str(), sbuf);
    if (rc < 0) {
        fmt::print("stat({}): {}\n", cachefile_name, strerror(errno));
        return;
    }
    size_t const fsize = sbuf->st_size;

#ifdef HAVE_STATVFS_H
    /* Check for remaining space on the filesystem */
        auto it = local_cache_copy.rfind('/');
        std::string dirname = (it == std::string::npos) ? "." : local_cache_copy.substr(0, it);

    struct statvfs sf[1];
    rc = statvfs(dirname.c_str(), sf);
    if (rc >= 0) {
        size_t const mb = sf->f_bsize * sf->f_bavail;

        if (fsize > size_t(0.5*double(mb))) {
            fmt::print(stderr, "Copying {} to {} would occupy {} MB out of {} MB available, so more than 50%. Skipping copy\n",
                    cachefile_name, dirname, fsize >> 20, mb >> 20);
            return;
        }
        fmt::print(stderr, "{} MB available on {}n", mb >> 20, dirname);
    } else {
        fmt::print(stderr, "Cannot do statvfs on {} (skipping check for available disk space): {}\n", dirname, strerror(errno));
    }
#endif


    if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_MAJOR_INFO)) {
        fmt::print(stderr, "Also saving cache data to {} ({} MB)\n",
                local_cache_copy,
                fsize >> 20);
    }

    std::string const normal_cachefile = cachefile_name;
    cachefile_name = local_cache_copy;
    save_cache();
    cachefile_name = normal_cachefile;
}

/* This is wicked. It's the main entry point of cache reloading, but it
 * dares play games with the cache name before the instance method gets
 * called. In turn, that instance method calls
 * matmul_common_reload_cache_fopen.
 *
 * It would make more sense if matmul_common_reload_cache_fopen chose
 * which file to load in its implementation, rather than here.
 */

int matmul_interface::reload_cache()
{
    struct stat sbuf[2][1];
    int rc;
    if (cachefile_name.empty())
        return 0;
    rc = stat(cachefile_name.c_str(), sbuf[0]);
    if (rc < 0) {
        return 0;
    }
    int local_is_ok = !local_cache_copy.empty();
    if (local_is_ok) {
        /* coverity[toctou] */
        rc = stat(local_cache_copy.c_str(), sbuf[1]);
        local_is_ok = rc == 0;
    }

    if (local_is_ok && (sbuf[0]->st_size != sbuf[1]->st_size)) {
        fmt::print(stderr, "{} and {} differ in size ; latter ignored\n",
                cachefile_name,
                local_cache_copy);
        unlink(local_cache_copy.c_str());
        local_is_ok = 0;
    }

    if (local_is_ok && (sbuf[0]->st_mtime > sbuf[1]->st_mtime)) {
        fmt::print(stderr, "{} is newer than {} ; latter ignored\n",
                cachefile_name,
                local_cache_copy);
        unlink(local_cache_copy.c_str());
        local_is_ok = 0;
    }

    if (!local_is_ok) {
        // no local copy.
        rc = reload_cache_private();
        if (rc == 0)
            return 0;
        // succeeded in loading data.
        save_to_local_copy();
        return 1;
    } else {
        std::string const normal_cachefile = cachefile_name;
        cachefile_name = local_cache_copy;
        rc = reload_cache_private();
        cachefile_name = normal_cachefile;
        return rc;
    }
}

void matmul_interface::save_cache()
{
    // nothing fancy here. Maybe more stuff could come, but it's probably
    // better to keep it clean. For symmetry, we have a front-end /
    // back-end mechanism for reload_cache, so we have one here, that's
    // it.
    save_cache_private();
}

#ifdef  BUILD_DYNAMICALLY_LINKABLE_BWC
void matmul_interface::remove_shared_lib_after_object::operator()(matmul_interface * mm) {
    delete mm;
    dlclose(solib_handle);
}
#endif
