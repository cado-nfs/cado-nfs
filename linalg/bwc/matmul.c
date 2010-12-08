#include "cado.h"
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/statvfs.h>


#include "bwc_config.h"
#include "matmul.h"

#define MATMUL_DEFAULT_IMPL bucket

/* Include all matmul implementations here */
#include "matmul-basic.h"
#include "matmul-sliced.h"
#include "matmul-threaded.h"
#include "matmul-bucket.h"


/* Now some cpp glue which sets up the different options */
#define MATMUL_NAME(kind,func) matmul_ ## kind ## _ ## func
#define MATMUL_T(func) matmul_ ## func ## _t

#define REBIND_F(mm, kind,func) \
        mm->bind->func = (MATMUL_T(func)) & MATMUL_NAME(kind, func)

#define SET_IMPL(mm, impl_) do { mm->bind->impl = # impl_ ; } while (0)

#define REBIND_ALL(mm, kind) do {					\
        REBIND_F(mm, kind, build_cache);				\
        REBIND_F(mm, kind, reload_cache);				\
        REBIND_F(mm, kind, save_cache);			        	\
        REBIND_F(mm, kind, mul);					\
        REBIND_F(mm, kind, report);					\
        REBIND_F(mm, kind, clear);					\
        REBIND_F(mm, kind, init);					\
        REBIND_F(mm, kind, auxv);					\
        REBIND_F(mm, kind, aux);					\
        SET_IMPL(mm, kind);                                             \
    } while (0)

/* these are just trampoline functions. A macro in the .h would serve the
 * same purpose, but could be misleading. Here at least, ctags bring the
 * reader quickly to the fact that we have macros.
 */

void do_rebinding(matmul_ptr mm, const char * impl)
{
#define CHECK_REBIND(mm, K) \
    if (strcmp(impl, #K) == 0) { REBIND_ALL(mm, K); } else

    if (impl == NULL) { REBIND_ALL(mm, MATMUL_DEFAULT_IMPL); } else
    CHECK_REBIND(mm, sliced)        // no semicolon !
    CHECK_REBIND(mm, threaded)      // no semicolon !
    CHECK_REBIND(mm, basic)         // no semicolon !
    CHECK_REBIND(mm, bucket)         // no semicolon !
    {   
        fprintf(stderr, "no implementation %s known (update %s ?)\n",
            impl, __FILE__);
        exit(1);
    }
}

matmul_ptr matmul_init(abobj_ptr x, unsigned int nr, unsigned int nc, const char * locfile, const char * impl, param_list pl, int optimized_direction)
{
    struct matmul_public_s fake[1];
    memset(fake, 0, sizeof(fake));
    do_rebinding(fake, impl);
    matmul_ptr mm = fake->bind->init(x, pl, optimized_direction);
    if (mm == NULL) return NULL;
    do_rebinding(mm, impl);

    mm->dim[0] = nr;
    mm->dim[1] = nc;

    // this just copies the locfile field from the parent structure
    // matmul_top... Except that for benches, the matmul structure lives
    // outside of this.
    mm->locfile = locfile;

    int rc = asprintf(&mm->cachefile_name, "%s-%s%s.bin", locfile, mm->bind->impl, mm->store_transposed ? "T" : "");
    FATAL_ERROR_CHECK(rc < 0, "out of memory");

    mm->local_cache_copy = NULL;

    if (!pl)
        return mm;

    const char * local_cache_copy_dir = param_list_lookup_string(pl, "local_cache_copy_dir");
    if (local_cache_copy_dir) {
        char * basename;
        basename = strrchr(mm->cachefile_name, '/');
        if (basename == NULL) {
            basename = mm->cachefile_name;
        }
        struct stat sbuf[1];
        rc = stat(local_cache_copy_dir, sbuf);
        if (rc < 0) {
            fprintf(stderr, "Warning: accessing %s is not possible: %s\n",
                    local_cache_copy_dir, strerror(errno));
            return mm;
        }
        int rc = asprintf(&mm->local_cache_copy, "%s/%s", local_cache_copy_dir, basename);
        FATAL_ERROR_CHECK(rc < 0, "out of memory");
    }

    return mm;
}

void matmul_build_cache(matmul_ptr mm, uint32_t * ptr)
{
    mm->bind->build_cache(mm, ptr);
}

static void save_to_local_copy(matmul_ptr mm)
{
    if (mm->local_cache_copy == NULL)
        return;

    struct stat sbuf[1];
    int rc;

    rc = stat(mm->cachefile_name, sbuf);
    if (rc < 0) {
        fprintf(stderr, "stat(%s): %s\n", mm->cachefile_name, strerror(errno));
        return;
    }
    unsigned long fsize = sbuf->st_size;

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
    if (rc < 0) {
        fprintf(stderr, "Cannot do statvfs on %s: %s\n", dirname, strerror(errno));
        free(dirname);
    }
    unsigned long mb = sf->f_bsize;
    mb *= sf->f_bavail;


    if (fsize > mb * 0.5) {
        fprintf(stderr, "Copying %s to %s would occupy %lu MB out of %lu MB available, so more than 50%%. Skipping copy\n",
                mm->cachefile_name, dirname, fsize >> 20, mb >> 20);
        free(dirname);
        return;
    }
    fprintf(stderr, "%lu MB available on %s\n", mb >> 20, dirname);

    free(dirname);

    fprintf(stderr, "Also saving cache data to %s (%lu MB)\n", mm->local_cache_copy, fsize >> 20);

    char * normal_cachefile = mm->cachefile_name;
    mm->cachefile_name = mm->local_cache_copy;
    mm->bind->save_cache(mm);
    mm->cachefile_name = normal_cachefile;
}

int matmul_reload_cache(matmul_ptr mm)
{
    struct stat sbuf[2][1];
    int rc;
    rc = stat(mm->cachefile_name, sbuf[0]);
    if (rc < 0) {
        return 0;
    }
    int local_is_ok = mm->local_cache_copy != NULL;
    if (local_is_ok) {
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
    mm->bind->save_cache(mm);
    save_to_local_copy(mm);
}

void matmul_mul(matmul_ptr mm, abt * dst, abt const * src, int d)
{
    mm->bind->mul(mm, dst, src, d);
}

void matmul_report(matmul_ptr mm, double scale)
{
    mm->bind->report(mm, scale);
}

void matmul_clear(matmul_ptr mm)
{
    if (mm->cachefile_name != NULL) free(mm->cachefile_name);
    if (mm->local_cache_copy != NULL) free(mm->local_cache_copy);
    mm->locfile = NULL;
    mm->bind->clear(mm);
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
