#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <stdint.h>
#include <string.h>
#include <stddef.h>

#include "macros.h"
#include "readmat.h"
#include "readmat-easy.h"

#define CHECK_REALLOC(ptr__,alloc__,size__,n__) do {    		\
        if (size__ + n__ >= alloc__) {					\
            alloc__ = size__ + n__;                                     \
            alloc__ += alloc__ / 2;                                     \
            ptr__ = realloc(ptr__, alloc__ * sizeof(uint32_t));		\
        }								\
    } while (0)

void read_easy_inner(const char * filename,
        uint32_t ** p_direct, uint32_t ** p_transposed,
        unsigned int * p_nr, unsigned int * p_nc,
        size_t * sz0, size_t * sz1)
{
    filter_matrix_t smat;
    FILE * f;

    filter_matrix_init(smat);
    f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "fopen(%s): %s\n", filename, strerror(errno));
        exit(1);
    }
    read_matrix_header(f, smat);

    *sz0 = *sz1 = 0;
    uint32_t * data;
    size_t alloc = 0;
    size_t size = 0;

    uint32_t * colweights;
    colweights = malloc(smat->ncols * sizeof(uint32_t));
    memset(colweights, 0, smat->ncols * sizeof(uint32_t));

    data = NULL;

    if (p_nr) *p_nr = smat->nrows;
    if (p_nc) *p_nc = smat->ncols;

    /* do it the bovine way. read just about everything in memory. */
    for(unsigned int i = 0 ; i < smat->nrows ; i++) {
        read_matrix_row(f,smat,smat->data,1);
        CHECK_REALLOC(data, alloc, size, 1 + smat->data[0]);
        data[size++] = smat->data[0];
        /* read_matrix_row now returns a sorted row */
        for(unsigned int j = 0 ; j < smat->data[0] ; j++) {
            data[size++] = smat->data[1+j];
            colweights[smat->data[1+j]]++;
        }
    }
    fclose(f);

    if (p_direct) {
        *p_direct = data;
    }
    *sz0 = size;

    if (p_transposed == NULL) {
        free(colweights);
        return;
    }

    uint32_t ** colptrs;
    colptrs = malloc(smat->ncols * sizeof(uint32_t *));

    if (smat->ncols > smat->nrows) {
        size += smat->ncols - smat->nrows;
    }

    uint32_t * res = malloc(size * sizeof(uint32_t));
    uint32_t * ptr = res;
    for(unsigned int j = 0 ; j < smat->ncols ; j++) {
        uint32_t w = colweights[j];
        ASSERT((size_t) (ptr - res) < size);
        *ptr++ = w;
        colptrs[j] = ptr;
        ptr += w;
    }

    uint32_t * v = data;
    for(unsigned int i = 0 ; i < smat->nrows ; i++) {
        uint32_t rlen = *v++;
        for(unsigned int j = 0 ; j < rlen ; j++) {
            uint32_t c = *v++;
            *(colptrs[c]++) = i;
        }
    }

    if (p_direct == NULL) {
        free(data);
    }

    ASSERT(p_transposed != NULL);
    *p_transposed = res;
    *sz1 = size;

    free(colweights);
    free(colptrs);
    filter_matrix_clear(smat);
}

void read_easy(const char * filename,
        uint32_t ** p_direct, uint32_t ** p_transposed,
        unsigned int * p_nr, unsigned int * p_nc)
{
    char * cname;
    int rc;
    rc = asprintf(&cname, "%s-easy.bin", filename);
    ASSERT_ALWAYS(rc >= 0);
    char * cnameT;
    rc = asprintf(&cnameT, "%s-easyT.bin", filename);
    ASSERT_ALWAYS(rc >= 0);

    FILE * f;
    if (p_direct) {
        f = fopen(cname, "r");
        if (f) {
            printf("Reusing intermediary cache file %s\n", cname);
            unsigned int nr;
            unsigned int nc;
            size_t fsize, rs;
            rs = fread(&nr, sizeof(unsigned int), 1, f);
            ASSERT_ALWAYS(rs == 1);
            rs = fread(&nc, sizeof(unsigned int), 1, f);
            ASSERT_ALWAYS(rs == 1);
            rs = fread(&fsize, sizeof(size_t), 1, f);
            ASSERT_ALWAYS(rs == 1);
            if (p_nr) { *p_nr = nr; p_nr = NULL; }
            if (p_nc) { *p_nc = nc; p_nc = NULL; }
            *p_direct = malloc(fsize * sizeof(uint32_t));
            rs = fread(*p_direct, sizeof(uint32_t), fsize, f);
            ASSERT_ALWAYS(rs == fsize);
            p_direct = NULL;
            fclose(f);
        }
    }

    if (p_transposed) {
        f = fopen(cnameT, "r");
        if (f) {
            printf("Reusing intermediary cache file %s\n", cnameT);
            unsigned int nc;
            unsigned int nr;
            size_t fsize;
            size_t rs;
            rs = fread(&nc, sizeof(unsigned int), 1, f);
            ASSERT_ALWAYS(rs == 1);
            rs = fread(&nr, sizeof(unsigned int), 1, f);
            ASSERT_ALWAYS(rs == 1);
            rs = fread(&fsize, sizeof(size_t), 1, f);
            ASSERT_ALWAYS(rs == 1);
            if (p_nr) { *p_nr = nr; p_nr = NULL; }
            if (p_nc) { *p_nc = nc; p_nc = NULL; }
            *p_transposed = malloc(fsize * sizeof(uint32_t));
            rs = fread(*p_transposed, sizeof(uint32_t), fsize, f);
            ASSERT_ALWAYS(rs == fsize);
            p_transposed = NULL;
            fclose(f);
        }
    }

    size_t sz = 0, szT = 0;
    unsigned int nr = 0;
    unsigned int nc = 0;

    if (p_direct || p_transposed || p_nr || p_nc) {
        read_easy_inner(filename, p_direct, p_transposed, &nr, &nc, &sz, &szT);
        if (p_nr) *p_nr = nr;
        if (p_nc) *p_nc = nc;
    }
#if 0
    /* Enable this code fragment ONLY when debugging internal mm building
     * routines. This speeds up considerably the upstart time, thus
     * facilitating debugging.
     */
    if (p_direct) {
        f = fopen(cname, "w");
        ASSERT_ALWAYS(f != NULL);
        printf("Saving intermediary cache file %s\n", cname);
        fwrite(&nr, sizeof(unsigned int), 1, f);
        fwrite(&nc, sizeof(unsigned int), 1, f);
        fwrite(&sz, sizeof(size_t), 1, f);
        fwrite(*p_direct, sizeof(uint32_t), sz, f);
        fclose(f);
    }
    if (p_transposed) {
        f = fopen(cnameT, "w");
        ASSERT_ALWAYS(f != NULL);
        printf("Saving intermediary cache file %s\n", cnameT);
        fwrite(&nc, sizeof(unsigned int), 1, f);
        fwrite(&nr, sizeof(unsigned int), 1, f);
        fwrite(&szT, sizeof(size_t), 1, f);
        fwrite(*p_transposed, sizeof(uint32_t), szT, f);
        fclose(f);
    }
#endif
    free(cname);
    free(cnameT);
}

/* vim: set sw=4: */
