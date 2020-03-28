#include "cado.h"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cinttypes>
#include <gmp.h>
#include "blockmatrix.hpp"
#include "bblas.hpp"
#include "macros.h"
#include "portability.h"
#include "cado-endian.h"

blockmatrix::blockmatrix(unsigned int nrows, unsigned int ncols)
{
    this->nrows = nrows;
    this->ncols = ncols;
    nrblocks = iceildiv(nrows, 64);
    ncblocks = iceildiv(ncols, 64);
    stride = nrblocks;
    if (nrblocks || ncblocks) {
        mb = (mat64 *) malloc(nrblocks * ncblocks * sizeof(mat64));
    } else {
        mb = NULL;
    }
    owner = true;
}

blockmatrix::~blockmatrix()
{
    if (owner && mb) free(mb);
    nrblocks = 0;
    ncblocks = 0;
    nrows = 0;
    ncols = 0;
    stride = 0;
}

blockmatrix::blockmatrix(blockmatrix & k, int i0, int j0, unsigned int nrows, unsigned int ncols)
{
    ASSERT_ALWAYS(i0 % 64 == 0);
    ASSERT_ALWAYS(j0 % 64 == 0);
    ASSERT_ALWAYS(i0 + nrows <= k.nrows);
    ASSERT_ALWAYS(j0 + ncols <= k.ncols);
    this->nrows = nrows;
    this->ncols = ncols;
    stride = k.nrblocks;
    nrblocks = iceildiv(nrows, 64);
    ncblocks = iceildiv(ncols, 64);
    mb = k.mb + (i0/64) + (j0/64) * k.stride;
    if (nrows == 0 || ncols == 0) mb = NULL;
    owner = false;
}
uint64_t * blockmatrix::subrow_ptr(int i, int j)
{
    ASSERT_ALWAYS(j % 64 == 0);
    return &(mb[(i/64) + (j/64) * stride][i%64]);
}
void blockmatrix::copy_colrange(blockmatrix const & A, int j0, int j1)
{
    ASSERT_ALWAYS(nrows == A.nrows);
    ASSERT_ALWAYS(ncols == A.ncols);
    ASSERT_ALWAYS(ncblocks == A.ncblocks);

    int block0 = j0 / 64;
    uint64_t * masks = (uint64_t *) malloc(A.ncblocks * sizeof(uint64_t));
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
    for(unsigned int i0 = 0 ; i0 < A.nrows ; i0 += 64) {
        for(unsigned int i = 0 ; i0 + i < A.nrows && i < 64 ; i++) {
            for(int b = block0 ; b*64 < j1 ; b++) {
                uint64_t m = masks[b];
                uint64_t v = A.mb[i0/64 + b*A.stride][i];
                v&=m;
                mb[i0/64 + b*stride][i] &= ~m;
                mb[i0/64 + b*stride][i] |= v;
            }
        }
    }
    free(masks);
}

void blockmatrix::print(const char *vname) const
{
    FILE * f = fopen("/tmp/debug", "a");
    if (f == NULL) {
        fprintf(stderr, "Cannot append to /tmp/debug");
        return;
    }
    fprintf(f,"%s:=Matrix(GF(2),%u,%u,[\n", vname, nrows, ncols);
    for (unsigned int i = 0; i < nrows; i += 64) {
	for (unsigned int ii = 0; ii < 64 && i + ii < nrows; ii++) {
	    if (i + ii)
		fprintf(f,",\n");
	    for (unsigned int j = 0; j < ncols; j += 64) {
                mat64 & bl = mb[i / 64 + (j / 64) * stride];
		for (unsigned int jj = 0; jj < 64 && j + jj < ncols; jj++) {
		    if (j + jj)
			fprintf(f,",");
		    fprintf(f,"%" PRIu64, (bl[ii] >> jj) & ((uint64_t)1));
		}
	    }
	}
    }
    fprintf(f,"]);\n");
    fclose(f);
}

void blockmatrix::set_zero()
{
    if (owner) {
        std::fill_n(mb, nrblocks * ncblocks, 0);
    } else {
        for(unsigned int j = 0 ; j < ncblocks ; j++)  
            std::fill_n(mb + j * stride, nrblocks, 0);
    }
}

void blockmatrix::set_identity()
{
    set_zero();
    for(unsigned int i = 0 ; i < ncols ; i++)
        mb[i/64 + (i/64)*stride][i%64] ^= ((uint64_t)1) << (i%64);
}

/* Computes transpose(a) * b */
void blockmatrix::mul_Ta_b(
        blockmatrix const & a,
        blockmatrix const & b)
{
    ASSERT_ALWAYS(  nrows == a.ncols);
    ASSERT_ALWAYS(  ncols == b.ncols);
    ASSERT_ALWAYS(a.nrows == b.nrows);

    if (this == &a || this == &b) {
        blockmatrix tc(a.ncols, b.ncols);
        tc.mul_Ta_b(a, b);
        tc.swap(*this);
        return;
    }

    ASSERT_ALWAYS(this != &b);
    ASSERT_ALWAYS(this != &a);

    set_zero();

    for(unsigned int i = 0 ; i < nrblocks ; i++) {
        for(unsigned int j = 0 ; j < ncblocks ; j++) {
            addmul_TN64_N64(mb[i + j * stride],
                    (uint64_t *) (a.mb + i * a.stride),
                    (uint64_t *) (b.mb + j * b.stride),
                    a.nrows);
        }
    }
}

void blockmatrix::mul_smallb(
        blockmatrix const & a,
        blockmatrix const & b)
{
    ASSERT_ALWAYS(  nrows == a.nrows);
    ASSERT_ALWAYS(  ncols == b.ncols);
    ASSERT_ALWAYS(a.ncols == b.nrows);

    if (this == &a || this == &b) {
        blockmatrix tc(a.nrows, b.ncols);
        tc.mul_smallb(a, b);
        tc.swap(*this);
        return;
    }

    ASSERT_ALWAYS(this != &b);
    ASSERT_ALWAYS(this != &a);

    set_zero();

    for(unsigned int j = 0 ; j < b.ncblocks ; j++) {
        for(unsigned int i = 0 ; i < b.nrblocks ; i++) {
            addmul_N64_6464(
                    (uint64_t *) (  mb + j *   stride),
                    (uint64_t *) (a.mb + i * a.stride),
                    b.mb[i + j * b.stride], a.nrows);
        }
    }
}

void blockmatrix::reverse_columns()
{
    if (nrows == ncols) {
        transpose(*this);
        reverse_rows();
        transpose(*this);
    } else {
        /* different stride. Better to allocate */
        blockmatrix ta(ncols, nrows);
        ta.transpose(*this);
        ta.reverse_rows();
        transpose(ta);
    }
}

void blockmatrix::reverse_rows()
{
    for(unsigned int ib = 0 ; ib < nrblocks ; ib++) {
        unsigned int xib = nrblocks - 1 - ib;
        if (xib < ib) break;
        for(unsigned int i = 0 ; i < 64 ; i++) {
            unsigned int ii = ib * 64 + i;
            unsigned int xi = 63 - i;
            unsigned int xii = xib * 64 + xi;
            if (xii < ii) break;
            for(unsigned int j = 0 ; j < ncblocks ; j++) {
                uint64_t t = mb[ib + j * stride][i];
                mb[ib + j * stride][i] = mb[xib + j * stride][xi];
                mb[xib + j * stride][xi] = t;
            }
        }
    }
}

/* Takes the block matrix m, and copy its data at position (i0, j0) inside
 * the tiny matrix.
 */
void blockmatrix::copy_to_flat(uint64_t * tiny, unsigned int flat_stride,
        int i0, int j0) const
{
    for(unsigned int i = 0 ; i < nrblocks ; i++) {
        for(unsigned int j = 0 ; j < ncblocks ; j++) {
            mat64 & tm = mb[i + j * stride];
            /* if needed, easy to transpose tm here and swap (i,j) */
            uint64_t * tp = tiny + (i0+i*64) * flat_stride + j0/64 + j;
            for(unsigned int k = 0 ; k < 64 ; k++)
                tp[k*flat_stride] = tm[k];
        }
    }
}

void blockmatrix::copy_transpose_to_flat(uint64_t * tiny, unsigned int flat_stride,
        int i0, int j0) const
{
    for(unsigned int i = 0 ; i < nrblocks ; i++) {
        for(unsigned int j = 0 ; j < ncblocks ; j++) {
            mat64 tm;
            mat64_transpose(tm, mb[i + j * stride]);
            uint64_t * tp = tiny + (i0+j*64) * flat_stride + j0/64 + i;
            /* Note that the tiny matrix must have space allocated for rows and
             * cols multiples of 64, otherwise disaster will occur */
            for(unsigned int k = 0 ; k < 64 ; k++)
                tp[k*flat_stride] = tm[k];
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
blockmatrix::swap_words_if_needed (uint64_t *v MAYBE_UNUSED, unsigned long n MAYBE_UNUSED)
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

void blockmatrix::copy_transpose_from_flat(uint64_t * tiny, unsigned int stride, int i0, int j0)
{
    for(unsigned int i = 0 ; i < nrblocks ; i++) {
        for(unsigned int j = 0 ; j < ncblocks ; j++) {
            mat64 tm;
            uint64_t * tp = tiny + (i0+j*64) * stride + j0/64 + i;
            for(unsigned int k = 0 ; k < 64 ; k++)
                tm[k] = tp[k*stride];
            mat64_transpose(mb[i + j * stride], tm);
        }
    }
}

void blockmatrix::copy_from_flat(uint64_t * tiny, unsigned int stride, int i0, int j0)
{
    for(unsigned int i = 0 ; i < nrblocks ; i++) {
        for(unsigned int j = 0 ; j < ncblocks ; j++) {
            mat64 & tm = mb[i + j * stride];
            uint64_t * tp = tiny + (i0+i*64) * stride + j0/64 + j;
            for(unsigned int k = 0 ; k < 64 ; k++)
                tm[k] = tp[k*stride];
        }
    }
}

/* reads matrix from file 'name',  considering the input as little endian */
void
blockmatrix::read_from_flat_file (int i0, int j0,
                                 const char * name, unsigned int fnrows,
                                 unsigned int fncols)
{
    FILE * f = fopen(name, "rb");
    ASSERT_ALWAYS(f);
    ASSERT_ALWAYS(i0 + fnrows <= nrows);
    ASSERT_ALWAYS(j0 + fncols <= ncols);
    ASSERT_ALWAYS(j0 % 64 == 0);
    ASSERT_ALWAYS(fncols % 64 == 0);
    set_zero();
    for(unsigned int r = 0 ; r < fnrows ; r++) {
        for(unsigned int g = 0 ; g < fncols ; g+=64) {
            uint64_t v;
            int rc = fread(&v, sizeof(uint64_t), 1, f);
            ASSERT_ALWAYS(rc == 1);
            mb[((i0+r)/64) + ((j0+g)/64) * stride][(i0+r)%64] = 
              u64_convert_from_little_endian (v);
        }
    }
    fclose(f);
}

#if 1
void blockmatrix::read_transpose_from_flat_file(int i0, int j0, const char * name, unsigned int fnrows, unsigned int fncols)
{
    FILE * f = fopen(name, "rb");
    ASSERT_ALWAYS(f);
    ASSERT_ALWAYS(i0 + fncols <= nrows);
    ASSERT_ALWAYS(j0 + fnrows <= ncols);
    ASSERT_ALWAYS(i0 % 64 == 0);
    ASSERT_ALWAYS(j0 % 64 == 0);
    ASSERT_ALWAYS(fncols % 64 == 0);
    set_zero();
    mat64 * tmp = (mat64 *) malloc((fncols/64) * sizeof(mat64));
    ASSERT_ALWAYS(tmp);
    for(unsigned int g = 0 ; g < fnrows ; g+=64) {
        /* Fill next block. */
        std::fill_n(tmp, fncols/64, 0);
        for(unsigned int i = 0 ; g + i < fnrows && i < 64 ; i++) {
            for(unsigned int s = 0 ; s < fncols ; s+=64) {
                uint64_t v;
                int rc = fread(&v, sizeof(uint64_t), 1, f);
                ASSERT_ALWAYS(rc == 1);
                tmp[s/64][i]=v;
            }
        }
        for(unsigned int s = 0 ; s < fncols ; s+=64) {
            mat64_transpose(mb[s/64 + (g/64) * stride], tmp[s/64]);
        }
    }
    free(tmp);
    fclose(f);
}
#endif

void blockmatrix::write_to_flat_file(const char * name, int i0, int j0, unsigned int fnrows, unsigned int fncols) const
{
    FILE * f = fopen(name, "wb");
    ASSERT_ALWAYS(f);
    ASSERT_ALWAYS(i0 + fnrows <= nrows);
    ASSERT_ALWAYS(j0 + fncols <= ncols);
    ASSERT_ALWAYS(j0 % 64 == 0);
    // for our convenience, we allow larger files.
    // ASSERT_ALWAYS(fncols % 64 == 0);
    for(unsigned int r = 0 ; r < fnrows ; r++) {
        for(unsigned int g = 0 ; g < fncols ; g+=64) {
            uint64_t v;
            v = mb[((i0+r)/64) + ((j0+g)/64) * stride][(i0+r)%64];
            v = u64_convert_to_little_endian (v);
            int rc = fwrite(&v, sizeof(uint64_t), 1, f);
            ASSERT_ALWAYS(rc == 1);
        }
    }
    fclose(f);
}

void blockmatrix::transpose(blockmatrix const & a)
{
    ASSERT_ALWAYS(nrows == a.ncols);
    ASSERT_ALWAYS(a.nrows == ncols);
    for(unsigned int i = 0 ; i < a.nrblocks ; i++) {
        mat64_transpose(mb[i + i * stride], a.mb[i + i * a.stride]);
        for(unsigned int j = i + 1 ; j < a.ncblocks ; j++) {
            mat64 tmp;
            mat64_transpose(tmp, a.mb[i + j * a.stride]);
            mat64_transpose(mb[i + j * stride], a.mb[j + i * a.stride]);
            mb[j + i * stride] = tmp;
        }
    }
}

void blockmatrix::swap(blockmatrix & A)
{
    std::swap(A, *this);
}

