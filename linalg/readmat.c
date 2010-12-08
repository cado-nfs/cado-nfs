#include "cado.h"
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "macros.h"
#include "readmat.h"

typedef int (*sortfunc_t) (const void *, const void *);

static int uint_cmp(unsigned int * a, unsigned int * b)
{
    if (*a < *b) return -1;
    else if (*b < *a) return 1;
    return 0;
}

unsigned int next_alloc_size(unsigned int s, unsigned int needed)
{
    /* most probably only one spin */
    for( ; s < needed ; ) {
        s = MAX(s + 1024, s + s / 4);
    }
    return s;
}

unsigned int *read_matrix_row(FILE * file, filter_matrix_t mat,
			      unsigned int *dst, int compact)
{
    unsigned int nc;
    unsigned int j;
    int ret;


    long a;
    unsigned long b;

    if (compact == 0) {
	ret = fscanf(file, "%ld %lu", &a, &b);
	ASSERT_ALWAYS(ret==2);

        /* FIXME -- having a and b stored in the sparse mat struct simply
         * won't scale ! */

        ASSERT_ALWAYS(a <= (long) INT_MAX);
        ASSERT_ALWAYS(a >= (long) INT_MIN);
        ASSERT_ALWAYS(b <= (unsigned long) INT_MAX);
    }

    ret = fscanf(file, "%d", &nc);
    ASSERT_ALWAYS(ret == 1);

    unsigned int offset = dst - mat->data;
    unsigned int needed = offset + nc + 1;

    if (compact == 0)
        needed += 2;

    if (needed > mat->alloc) {
	unsigned int next = next_alloc_size(mat->alloc, needed);
	mat->data = (unsigned int *) realloc(mat->data,
					     next * sizeof(unsigned int));
	mat->alloc = next;
	dst = mat->data + offset;
    }

    if (compact == 0) {
        *dst++ = (int) a;
        *dst++ = (int) b;
    }

    *dst++ = nc;
    unsigned int last = 0;
    int sorted = 1;
    unsigned int * head = dst;
    for (j = 0; j < nc; ++j) {
	unsigned int x;
	ret = fscanf(file, "%u", &x);
	ASSERT_ALWAYS(ret == 1);
        if (x < last)
            sorted = 0;
        last = x;
	*dst++ = x;
    }
    if (!sorted)
        qsort(head, nc, sizeof(unsigned int), (sortfunc_t) uint_cmp);

    mat->wt += nc;
    return dst;
}

void read_matrix_header(FILE * file, filter_matrix_t mat)
{
    int ret;

    ret = fscanf(file, "%d %d", &(mat->nrows), &(mat->ncols));
    ASSERT_ALWAYS(ret==2);
    mat->alloc = mat->nrows * 10;
    mat->data = (unsigned int *) malloc(mat->alloc * sizeof(unsigned int));
    mat->wt = 0;
}

void readmat(FILE * file, filter_matrix_t mat, int compact)
{
    unsigned int i;
    unsigned int *dst;
    read_matrix_header(file, mat);
    dst = mat->data;
    for (i = 0; i < mat->nrows; ++i) {
	dst = read_matrix_row(file, mat, dst, compact);
    }
}

void filter_matrix_init(filter_matrix_t dst)
{
    memset(dst, 0, sizeof(filter_matrix_t));
}

void filter_matrix_clear(filter_matrix_t dst)
{
    free(dst->data);
    memset(dst, 0, sizeof(filter_matrix_t));
}

/* vim: set sw=4 sta et: */
