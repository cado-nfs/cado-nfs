#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <gmp.h>

#include "blockmatrix.h"
#include "bit_matrices.h"
#include "macros.h"
#include "cado-endian.h"

blockmatrix blockmatrix_alloc(unsigned int nrows, unsigned int ncols)
{
    blockmatrix res = malloc(sizeof(struct blockmatrix_s));
    res->nrows = nrows;
    res->ncols = ncols;
    res->nrblocks = iceildiv(nrows, 64);
    res->ncblocks = iceildiv(ncols, 64);
    res->stride = res->nrblocks;
    if (nrows == 0 || ncols == 0) {
        res->mb = NULL;
    } else {
        res->mb = malloc(res->nrblocks * res->ncblocks * sizeof(mat64));
    }
    res->sub = 0;
    return res;
}

void blockmatrix_free(blockmatrix b)
{
    ASSERT_ALWAYS(b->sub == 0);
    free(b->mb); b->mb = NULL;
    b->nrblocks = 0;
    b->ncblocks = 0;
    b->nrows = 0;
    b->ncols = 0;
    b->stride = 0;
    free(b);
}

uint64_t * blockmatrix_subrow_ptr(blockmatrix res, int i, int j)
{
    ASSERT_ALWAYS(j % 64 == 0);
    return &(res->mb[(i/64) + (j/64) * res->stride][(i)%64]);
}
void blockmatrix_copy_colrange(blockmatrix B, blockmatrix A, int j0, int j1)
{
    ASSERT_ALWAYS(A->nrows == B->nrows);
    ASSERT_ALWAYS(A->ncols == B->ncols);
    ASSERT_ALWAYS(A->ncblocks == B->ncblocks);

    int block0 = j0 / 64;
    uint64_t * masks = malloc(A->ncblocks * sizeof(uint64_t));
    for(int b = block0 ; b*64 < j1 ; b++) {
        uint64_t mask = -((uint64_t)1);
        int z0 = j0 - b*64;
        ASSERT_ALWAYS(z0 < 64);
        if (z0>=0) mask &= ((uint64_t)-1) << z0;
        int z1 = 64 - (j1 - b*64);
        ASSERT_ALWAYS(z1 < 64);
        if (z1>=0) mask &= ((uint64_t)-1) >> z1;
        masks[b] = mask;
    }
    for(unsigned int i0 = 0 ; i0 < A->nrows ; i0 += 64) {
        for(unsigned int i = 0 ; i0 + i < A->nrows && i < 64 ; i++) {
            for(int b = block0 ; b*64 < j1 ; b++) {
                uint64_t m = masks[b];
                uint64_t v = A->mb[i0/64 + b*A->stride][i];
                v&=m;
                B->mb[i0/64 + b*B->stride][i] &= ~m;
                B->mb[i0/64 + b*B->stride][i] |= v;
            }
        }
    }
    free(masks);
}

void blockmatrix_print(blockmatrix b, const char *vname)
{
    FILE * f = fopen("/tmp/debug", "a");
    if (f == NULL) {
        fprintf(stderr, "Cannot append to /tmp/debug");
        return;
    }
    fprintf(f,"%s:=Matrix(GF(2),%u,%u,[\n", vname, b->nrows, b->ncols);
    for (unsigned int i = 0; i < b->nrows; i += 64) {
	for (unsigned int ii = 0; ii < 64 && i + ii < b->nrows; ii++) {
	    if (i + ii)
		fprintf(f,",\n");
	    for (unsigned int j = 0; j < b->ncols; j += 64) {
                mat64_ptr bl = b->mb[i / 64 + (j / 64) * b->stride];
		for (unsigned int jj = 0; jj < 64 && j + jj < b->ncols; jj++) {
		    if (j + jj)
			fprintf(f,",");
		    fprintf(f,"%"PRIu64, (bl[ii] >> jj) & ((uint64_t)1));
		}
	    }
	}
    }
    fprintf(f,"]);\n");
    fclose(f);
}

void blockmatrix_set_zero(blockmatrix b)
{
    if (!b->sub) {
        memset(b->mb, 0, b->nrblocks * b->ncblocks * sizeof(mat64));
    } else {
        for(unsigned int j = 0 ; j < b->ncblocks ; j++)  
            memset(b->mb + j * b->stride, 0, b->nrblocks * sizeof(mat64));
    }
}

void blockmatrix_set_identity(blockmatrix S)
{
    blockmatrix_set_zero(S);
    for(unsigned int i = 0 ; i < S->ncols ; i++)
        S->mb[i/64 + (i/64)*S->stride][i%64] ^= ((uint64_t)1) << (i%64);
}

/* Computes transpose(a) * b */
void blockmatrix_mul_Ta_b(blockmatrix c,
        const blockmatrix a,
        const blockmatrix b)
{
    ASSERT_ALWAYS(c->nrows == a->ncols);
    ASSERT_ALWAYS(c->ncols == b->ncols);
    ASSERT_ALWAYS(a->nrows == b->nrows);

    ASSERT_ALWAYS(a != b);
    ASSERT_ALWAYS(a != c);

    blockmatrix_set_zero(c);

    for(unsigned int i = 0 ; i < c->nrblocks ; i++) {
        for(unsigned int j = 0 ; j < c->ncblocks ; j++) {
            addmul_TN64_N64(c->mb[i + j * c->stride],
                    (uint64_t *) (a->mb + i * a->stride),
                    (uint64_t *) (b->mb + j * b->stride),
                    a->nrows);
        }
    }
}

void blockmatrix_mul_smallb(blockmatrix c,
        const blockmatrix a,
        const blockmatrix b)
{
    ASSERT_ALWAYS(c->nrows == a->nrows);
    ASSERT_ALWAYS(c->ncols == b->ncols);
    ASSERT_ALWAYS(a->ncols == b->nrows);

    ASSERT_ALWAYS(a != b);
    ASSERT_ALWAYS(a != c);

    blockmatrix_set_zero(c);

    for(unsigned int j = 0 ; j < b->ncblocks ; j++) {
        for(unsigned int i = 0 ; i < b->nrblocks ; i++) {
            addmul_N64_6464(
                    (uint64_t *) (c->mb + j * c->stride),
                    (uint64_t *) (a->mb + i * a->stride),
                    b->mb[i + j * b->stride], a->nrows);
        }
    }
}

/* Takes the block matrix m, and copy its data at position (i0, j0) inside
 * the tiny matrix.
 */
void blockmatrix_copy_to_flat(uint64_t * tiny, unsigned int stride,
        int i0, int j0, blockmatrix m)
{
    for(unsigned int i = 0 ; i < m->nrblocks ; i++) {
        for(unsigned int j = 0 ; j < m->ncblocks ; j++) {
            mat64_ptr tm = m->mb[i + j * m->stride];
            /* if needed, easy to transpose tm here and swap (i,j) */
            uint64_t * tp = tiny + (i0+i*64) * stride + j0/64 + j;
            for(unsigned int k = 0 ; k < 64 ; k++)
                tp[k*stride] = tm[k];
        }
    }
}

void blockmatrix_copy_transpose_to_flat(uint64_t * tiny, unsigned int stride,
        int i0, int j0, blockmatrix m)
{
    for(unsigned int i = 0 ; i < m->nrblocks ; i++) {
        for(unsigned int j = 0 ; j < m->ncblocks ; j++) {
            mat64 tm;
            transp_6464(tm, m->mb[i + j * m->stride]);
            uint64_t * tp = tiny + (i0+j*64) * stride + j0/64 + i;
            /* Note that the tiny matrix must have space allocated for rows and
             * cols multiples of 64, otherwise disaster will occur */
            for(unsigned int k = 0 ; k < 64 ; k++)
                tp[k*stride] = tm[k];
        }
    }
}

/* swap characters in a 64-bit word if necessary */
static uint64_t
u64_convert_to_little_endian (uint64_t v)
{
#if CADO_BYTE_ORDER == 1234 /* little endian: nothing to do */
  return v;
#elif CADO_BYTE_ORDER == 4321
  return ((v & 255) << 56) + (((v >> 8) & 255) << 48)
    + (((v >> 16) & 255) << 40) + (((v >> 24) & 255) << 32)
    + (((v >> 32) & 255) << 24) + (((v >> 40) & 255) << 16)
    + (((v >> 48) & 255) << 8) + ((v >> 56) & 255);
#else
#error "neither little nor big endian: implement me"
#endif
}
static uint64_t
u64_convert_from_little_endian (uint64_t v)
{
    return u64_convert_to_little_endian(v);
}

/* if mp_limb_t has 32 bits and we are on a big-endian machine, swap
   32-bit words */
void
swap_words_if_needed (uint64_t *v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
{
#if CADO_BYTE_ORDER == 1234
  /* do nothing */
#elif CADO_BYTE_ORDER == 4321
  if (GMP_LIMB_BITS == 32)
    {
      unsigned long i;
      for (i = 0; i < n; i++)
        v[i] = (v[i] >> 32) + ((v[i] & 4294967295UL) << 32);
    }
#else
#error "neither little nor big endian: implement me"
#endif
}

void blockmatrix_copy_transpose_from_flat(blockmatrix m, uint64_t * tiny, unsigned int stride, int i0, int j0)
{
    for(unsigned int i = 0 ; i < m->nrblocks ; i++) {
        for(unsigned int j = 0 ; j < m->ncblocks ; j++) {
            mat64 tm;
            uint64_t * tp = tiny + (i0+j*64) * stride + j0/64 + i;
            for(unsigned int k = 0 ; k < 64 ; k++)
                tm[k] = tp[k*stride];
            transp_6464(m->mb[i + j * m->stride], tm);
        }
    }
}

void blockmatrix_copy_from_flat(blockmatrix m, uint64_t * tiny, unsigned int stride, int i0, int j0)
{
    for(unsigned int i = 0 ; i < m->nrblocks ; i++) {
        for(unsigned int j = 0 ; j < m->ncblocks ; j++) {
            mat64_ptr tm = m->mb[i + j * m->stride];
            uint64_t * tp = tiny + (i0+i*64) * stride + j0/64 + j;
            for(unsigned int k = 0 ; k < 64 ; k++)
                tm[k] = tp[k*stride];
        }
    }
}

/* reads matrix from file 'name',  considering the input as little endian */
void
blockmatrix_read_from_flat_file (blockmatrix k, int i0, int j0,
                                 const char * name, unsigned int fnrows,
                                 unsigned int fncols)
{
    FILE * f = fopen(name, "r");
    ASSERT_ALWAYS(f);
    ASSERT_ALWAYS(i0 + fnrows <= k->nrows);
    ASSERT_ALWAYS(j0 + fncols <= k->ncols);
    ASSERT_ALWAYS(j0 % 64 == 0);
    ASSERT_ALWAYS(fncols % 64 == 0);
    blockmatrix_set_zero(k);
    for(unsigned int r = 0 ; r < fnrows ; r++) {
        for(unsigned int g = 0 ; g < fncols ; g+=64) {
            uint64_t v;
            int rc = fread(&v, sizeof(uint64_t), 1, f);
            ASSERT_ALWAYS(rc == 1);
            k->mb[((i0+r)/64) + ((j0+g)/64) * k->stride][(i0+r)%64] = 
              u64_convert_from_little_endian (v);
        }
    }
    fclose(f);
}

#if 1
void blockmatrix_read_transpose_from_flat_file(blockmatrix k, int i0, int j0, const char * name, unsigned int fnrows, unsigned int fncols)
{
    FILE * f = fopen(name, "r");
    ASSERT_ALWAYS(f);
    ASSERT_ALWAYS(i0 + fncols <= k->nrows);
    ASSERT_ALWAYS(j0 + fnrows <= k->ncols);
    ASSERT_ALWAYS(i0 % 64 == 0);
    ASSERT_ALWAYS(j0 % 64 == 0);
    ASSERT_ALWAYS(fncols % 64 == 0);
    blockmatrix_set_zero(k);
    mat64 * tmp = malloc((fncols/64) * sizeof(mat64));
    ASSERT_ALWAYS(tmp);
    for(unsigned int g = 0 ; g < fnrows ; g+=64) {
        /* Fill next block. */
        memset(tmp, 0, (fncols/64) * sizeof(mat64));
        for(unsigned int i = 0 ; g + i < fnrows && i < 64 ; i++) {
            for(unsigned int s = 0 ; s < fncols ; s+=64) {
                uint64_t v;
                int rc = fread(&v, sizeof(uint64_t), 1, f);
                ASSERT_ALWAYS(rc == 1);
                tmp[s/64][i]=v;
            }
        }
        for(unsigned int s = 0 ; s < fncols ; s+=64) {
            transp_6464(k->mb[s/64 + (g/64) * k->stride], tmp[s/64]);
        }
    }
    free(tmp);
    fclose(f);
}
#endif

void blockmatrix_write_to_flat_file(const char * name, blockmatrix k, int i0, int j0, unsigned int fnrows, unsigned int fncols)
{
    FILE * f = fopen(name, "w");
    ASSERT_ALWAYS(f);
    ASSERT_ALWAYS(i0 + fnrows <= k->nrows);
    ASSERT_ALWAYS(j0 + fncols <= k->ncols);
    ASSERT_ALWAYS(j0 % 64 == 0);
    // for our convenience, we allow larger files.
    // ASSERT_ALWAYS(fncols % 64 == 0);
    for(unsigned int r = 0 ; r < fnrows ; r++) {
        for(unsigned int g = 0 ; g < fncols ; g+=64) {
            uint64_t v;
            v = k->mb[((i0+r)/64) + ((j0+g)/64) * k->stride][(i0+r)%64];
            v = u64_convert_to_little_endian (v);
            int rc = fwrite(&v, sizeof(uint64_t), 1, f);
            ASSERT_ALWAYS(rc == 1);
        }
    }
    fclose(f);
}

void blockmatrix_transpose(blockmatrix b, blockmatrix a)
{
    ASSERT_ALWAYS(b->nrows == a->ncols);
    ASSERT_ALWAYS(a->nrows == b->ncols);
    for(unsigned int i = 0 ; i < a->nrblocks ; i++) {
        for(unsigned int j = 0 ; j < a->ncblocks ; j++) {
            transp_6464(b->mb[j + i * b->stride], a->mb[i + j * a->stride]);
        }
    }
}

blockmatrix blockmatrix_submatrix(blockmatrix k, int i0, int j0, unsigned int nrows, unsigned int ncols)
{
    ASSERT_ALWAYS(i0 % 64 == 0);
    ASSERT_ALWAYS(j0 % 64 == 0);
    ASSERT_ALWAYS(i0 + nrows <= k->nrows);
    ASSERT_ALWAYS(j0 + ncols <= k->ncols);
    blockmatrix res = malloc(sizeof(struct blockmatrix_s));
    res->nrows = nrows;
    res->ncols = ncols;
    res->stride = k->nrblocks;
    res->nrblocks = iceildiv(nrows, 64);
    res->ncblocks = iceildiv(ncols, 64);
    /* Do ***NOT*** call blockmatrix_free on a submatrix !!! */
    res->mb = k->mb + (i0/64) + (j0/64) * k->stride;
    if (nrows == 0 || ncols == 0) res->mb = NULL;
    res->sub = 1;
    return res;
}
void blockmatrix_swap(blockmatrix B, blockmatrix A)
{
    struct blockmatrix_s tmp[1];
    memcpy(tmp, A, sizeof(struct blockmatrix_s));
    memcpy(A, B, sizeof(struct blockmatrix_s));
    memcpy(B, tmp, sizeof(struct blockmatrix_s));
}

