#include "cado.h"
#include <errno.h>
#include <string.h>
#include "bwc_config.h"
#include "matmul-common.h"
#include "params.h"
#include "readmat-easy.h"

#define MM_COMMON_MAGIC 0xb0010001UL

const char * rowcol[2] = { "row", "col", };

/* Factor out some stuff which turns out to appear fairly often */

FILE * matmul_common_reload_cache_fopen(size_t stride, struct matmul_public_s * mm, uint32_t magic)
{
    FILE * f = fopen(mm->cachefile_name, "r");
    if (f == NULL) return NULL;

    // mm->cachefile_name is a cache file for mm->locfile (which in
    // general never exists)

    uint32_t magic_check;
    MATMUL_COMMON_READ_ONE32(magic_check, f);
    
    if (magic_check != magic) {
        fprintf(stderr, "Wrong magic in cached matrix file\n");
        return NULL;
    }   
    
    MATMUL_COMMON_READ_ONE32(magic_check, f);
    if (magic_check != MM_COMMON_MAGIC) {
        fprintf(stderr, "Wrong magic in cached matrix file\n");
        return NULL;
    }   

    /* Four reserved bytes for alignment */
    MATMUL_COMMON_READ_ONE32(magic_check, f);

    uint32_t nbytes_check;
    MATMUL_COMMON_READ_ONE32(nbytes_check, f);
    /* It's not fatal. It only deserves a warning */
    if (nbytes_check != stride) {
        fprintf(stderr, "Warning: cached matrix file fits data with different striding\n");
    }

    MATMUL_COMMON_READ_ONE32(mm->dim[0], f);
    MATMUL_COMMON_READ_ONE32(mm->dim[1], f);

    MATMUL_COMMON_READ_ONE64(mm->ncoeffs, f);

    return f;
}

FILE * matmul_common_save_cache_fopen(size_t stride, struct matmul_public_s * mm, uint32_t magic)
{
    FILE * f = fopen(mm->cachefile_name, "w");
    if (f == NULL) {
        fprintf(stderr, "Cannot open %s for writing: %s\n", mm->cachefile_name, strerror(errno));
        abort();
    }

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

/* Okay, this matrix reading stage is really a memory hog. But it's
 * done only once, so we practically don't care. */
uint32_t * matmul_common_read_stupid_data(struct matmul_public_s * mm)
{
    uint32_t * data;
    unsigned int nr, nc;
    if (mm->store_transposed) {
        read_easy(mm->locfile, NULL, &data, &nr, &nc);
    } else {
        read_easy(mm->locfile, &data, NULL, &nr, &nc);
    }
    if (mm->dim[0] == 0 && mm->dim[1] == 0) {
        mm->dim[0] = nr;
        mm->dim[1] = nc;
    }
    ASSERT_ALWAYS(mm->dim[0] == nr);
    ASSERT_ALWAYS(mm->dim[1] == nc);
    return data;
}
