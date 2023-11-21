#include "cado.h" // IWYU pragma: keep
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>   // for abort
#include "matmul.hpp"   // for matmul_public_s
#include "matmul-common.hpp"
#include "verbose.h"

#define MM_COMMON_MAGIC 0xb0010003UL

const char * rowcol[2] = { "row", "col", };

/* Factor out some stuff which turns out to appear fairly often */

FILE * matmul_common_reload_cache_fopen(size_t stride, struct matmul_public_s * mm, uint32_t magic)
{
    if (!mm->cachefile_name) return NULL;
    FILE * f = fopen(mm->cachefile_name, "rb");
    if (f == NULL) return NULL;

    // mm->cachefile_name is a cache file for mm->locfile (which in
    // general never exists)

    uint32_t magic_check;
    MATMUL_COMMON_READ_ONE32(magic_check, f);
    
    if (magic_check != magic) {
        fprintf(stderr, "Wrong magic in cached matrix file\n");
        fclose(f);
        return NULL;
    }   
    
    MATMUL_COMMON_READ_ONE32(magic_check, f);
    if (magic_check != MM_COMMON_MAGIC) {
        fprintf(stderr, "Wrong magic in cached matrix file\n");
        fclose(f);
        return NULL;
    }   

    /* Four reserved bytes for alignment */
    MATMUL_COMMON_READ_ONE32(magic_check, f);

    uint32_t nbytes_check;
    MATMUL_COMMON_READ_ONE32(nbytes_check, f);
    /* It's not fatal. It only deserves a warning */
    if (nbytes_check != stride) {
        fprintf(stderr, "Warning: cached matrix file fits data with different stride\n");
    }

    MATMUL_COMMON_READ_ONE32(mm->dim[0], f);
    MATMUL_COMMON_READ_ONE32(mm->dim[1], f);

    MATMUL_COMMON_READ_ONE64(mm->ncoeffs, f);

    return f;
}

FILE * matmul_common_save_cache_fopen(size_t stride, struct matmul_public_s * mm, uint32_t magic)
{
    if (!mm->cachefile_name) return NULL;
    FILE * f = fopen(mm->cachefile_name, "wb");
    if (f == NULL) {
        fprintf(stderr, "Cannot open %s for writing: %s\n", mm->cachefile_name, strerror(errno));
        abort();
    }

    if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_MAJOR_INFO))
        printf("Saving %s to cache file %s\n", mm->locfile, mm->cachefile_name);

    MATMUL_COMMON_WRITE_ONE32(magic,f);
    MATMUL_COMMON_WRITE_ONE32(MM_COMMON_MAGIC,f);
    MATMUL_COMMON_WRITE_ONE32(0,f);
    MATMUL_COMMON_WRITE_ONE32(stride,f);
    MATMUL_COMMON_WRITE_ONE32(mm->dim[0],f);
    MATMUL_COMMON_WRITE_ONE32(mm->dim[1],f);
    MATMUL_COMMON_WRITE_ONE64(mm->ncoeffs,f);
    return f;
}

void matmul_common_clear(struct matmul_public_s * mm MAYBE_UNUSED)
{
}
