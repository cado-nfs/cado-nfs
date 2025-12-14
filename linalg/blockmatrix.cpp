#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cinttypes>
#include <algorithm>          // for fill_n
#include "blockmatrix.hpp"
#include "bblas.hpp"
#include "macros.h"
#include "cado-endian.h"

blockmatrix::blockmatrix(unsigned int _nrows, unsigned int _ncols)
{
    this->_nrows = _nrows;
    this->_ncols = _ncols;
    nrblocks = iceildiv(_nrows, 64);
    ncblocks = iceildiv(_ncols, 64);
    if (nrblocks || ncblocks) {
        mb = mat64::alloc(nrblocks * ncblocks);
    } else {
        mb = NULL;
    }
}

blockmatrix::~blockmatrix()
{
    if (mb) mat64::free(mb, nrblocks * ncblocks);
    nrblocks = 0;
    ncblocks = 0;
    _nrows = 0;
    _ncols = 0;
}

void blockmatrix::copy_colrange(blockmatrix const & A, unsigned int j0, unsigned int j1)
{
    ASSERT_ALWAYS(_nrows == A._nrows);
    ASSERT_ALWAYS(_ncols == A._ncols);
    ASSERT_ALWAYS(ncblocks == A.ncblocks);

    constexpr const unsigned int B = matrix::width;
    using U = matrix::datatype;

    unsigned int const block0 = j0 / B;
    U * masks = new U[A.ncblocks];
    for(unsigned int b = block0 ; b*B < j1 ; b++) {
        U mask = -U(1);
        if (j0 >= b * B) {
            unsigned int const z0 = j0 - b*B;
            ASSERT_ALWAYS(z0 < B);
            mask &= (U(-1)) << z0;
        }
        if (j1 - b * B <= B) {
            unsigned int const z1 = B - (j1 - b*B);
            /* z1 is between 0 and B-1 */
            mask &= U(-1) >> z1;
        }
        masks[b] = mask;
    }
    for(unsigned int ii = 0 ; ii < A.nrows() ; ii++) {
        for(unsigned int jj = block0 * B ; jj < j1 ; jj += B) {
            U const m = masks[jj / B];
            (*this)[ii][jj] &= ~m;
            (*this)[ii][jj] |= A[ii][jj] & m;
        }
    }

    delete[] masks;
}

void blockmatrix::print(const char *vname) const
{
    constexpr const unsigned int B = matrix::width;
    FILE * f = fopen("/tmp/debug", "a");
    if (f == NULL) {
        fprintf(stderr, "Cannot append to /tmp/debug");
        return;
    }
    fprintf(f,"%s:=Matrix(GF(2),%u,%u,[\n", vname, _nrows, _ncols);
    for (unsigned int i = 0; i < _nrows; i += B) {
	for (unsigned int ii = 0; ii < B && i + ii < _nrows; ii++) {
	    if (i + ii)
		fprintf(f,",\n");
	    for (unsigned int j = 0; j < _ncols; j += B) {
                matrix const & bl = getblock(i/B, j/B);
		for (unsigned int jj = 0; jj < B && j + jj < _ncols; jj++) {
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

void blockmatrix::view_t::set_zero()
{
    for(unsigned int bi = 0 ; bi < nrblocks() ; bi++)
    for(unsigned int bj = 0 ; bj < ncblocks() ; bj++)
        getblock(bi, bj) = 0;
}

void blockmatrix::set_zero()
{
    std::fill_n(mb, nrblocks * ncblocks, 0);
}

void blockmatrix::set_identity()
{
    for(unsigned int bi = 0 ; bi < nrblocks ; bi++)
    for(unsigned int bj = 0 ; bj < ncblocks ; bj++)
        getblock(bi, bj) = (bi == bj);
}

/* Computes transpose(a) * b */
void blockmatrix::mul_Ta_b(
        blockmatrix::view_t c,
        blockmatrix::const_view_t a,
        blockmatrix::const_view_t b)
{
    ASSERT_ALWAYS(c.nrows() == a.ncols());
    ASSERT_ALWAYS(c.ncols() == b.ncols());
    ASSERT_ALWAYS(a.nrows() == b.nrows());

    ASSERT_ALWAYS(&c.M != &a.M);
    ASSERT_ALWAYS(&c.M != &b.M);

    c.set_zero();

    constexpr const unsigned int B = matrix::width;

    for(unsigned int i = 0 ; i < c.nrblocks() ; i++) {
        for(unsigned int j = 0 ; j < c.ncblocks() ; j++) {
            if (0) { // c.M.rowstride() == 1 && a.M.rowstride() == 1) {
                /* This code is fine, too, but the one below is more
                 * generic and probably just as fast. We want it to have
                 * the same level of testing, so let's use the more
                 * generic code in all cases.
                 */
                addmul_TN64_N64(c.getblock(i, j),
                        (uint64_t *) &a.getblock(0, i),
                        (uint64_t *) &b.getblock(0, j),
                        a.nrows());
            } else {
                for(unsigned int k = 0 ; k < a.nrblocks() ; k++) {
                    addmul_TN64_N64(c.getblock(i, j),
                            (uint64_t *) &a.getblock(k, i),
                            (uint64_t *) &b.getblock(k, j),
                            B);
                }
            }
        }
    }
}

void blockmatrix::mul_smallb(view_t a, const_view_t b)
{
    ASSERT_ALWAYS(a.ncols() == b.nrows());

    if (b.ncblocks() == 1 && b.nrblocks() == 1) {
        for(unsigned int i = 0 ; i < a.nrblocks() ; i++) {
            mul_6464_6464(
                    a.getblock(i, 0),
                    a.getblock(i, 0),
                    b.getblock(0, 0));
        }
    } else {
            blockmatrix c(a.nrows(), b.ncols());
            blockmatrix::mul_smallb(c, a, b);
            a.M.swap(c);
    }
}

void blockmatrix::mul_smallb(blockmatrix & c, const_view_t a, const_view_t b)
{
    ASSERT_ALWAYS(a.ncols() == b.nrows());
    ASSERT_ALWAYS(c.ncols() == b.ncols());
    ASSERT_ALWAYS(&c != &a.M && &c != &b.M);

    c.set_zero();
    for(unsigned int j = 0 ; j < b.ncblocks() ; j++) {
        ASSERT_ALWAYS(c.nrblocks == a.nrblocks());
        for(unsigned int k = 0 ; k < b.nrblocks() ; k++) {
            for(unsigned int i = 0 ; i < a.nrblocks() ; i++) {
                addmul_6464_6464(
                        c.getblock(i, j),
                        a.getblock(i, k),
                        b.getblock(k, j));
            }
        }
    }
}

void blockmatrix::reverse_columns()
{
    if (_nrows == _ncols) {
        transpose(*this);
        reverse_rows();
        transpose(*this);
    } else {
        /* different stride. Better to allocate */
        blockmatrix ta(_ncols, _nrows);
        ta.transpose(*this);
        ta.reverse_rows();
        transpose(ta);
    }
}

void blockmatrix::reverse_rows()
{
    for(unsigned int ib = 0 ; ib < nrblocks ; ib++) {
        unsigned int const xib = nrblocks - 1 - ib;
        if (xib < ib) break;
        for(unsigned int i = 0 ; i < 64 ; i++) {
            unsigned int const ii = ib * 64 + i;
            unsigned int const xi = 63 - i;
            unsigned int const xii = xib * 64 + xi;
            if (xii < ii) break;
            for(unsigned int j = 0 ; j < ncblocks ; j++) {
                uint64_t const t = getblock(ib, j)[i];
                getblock(ib, j)[i] = getblock(xib, j)[xi];
                getblock(xib, j)[xi] = t;
            }
        }
    }
}

/* Takes the block matrix m, and copy its data at position (i0, j0) inside
 * the tiny matrix.
 */
void blockmatrix::copy_to_flat(flat_area F, const_view_t V)
{
    using U = matrix::datatype;
    constexpr const unsigned int B = matrix::width;
    unsigned int const i0 = 0;
    unsigned int const j0 = 0;
    for(unsigned int i = 0 ; i < V.nrblocks() ; i++) {
        for(unsigned int j = 0 ; j < V.ncblocks() ; j++) {
            matrix const & tm = V.getblock(i, j);
            /* if needed, easy to transpose tm here and swap (i,j) */
            U * tp = F.p + (i0+i*B) * F.row_stride + j0/B + j;
            for(unsigned int k = 0 ; k < B ; k++)
                tp[k*F.row_stride] = tm[k];
        }
    }
}

void blockmatrix::copy_transpose_to_flat(flat_area F, const_view_t V)
{
    using U = matrix::datatype;
    constexpr const unsigned int B = matrix::width;
    unsigned int const i0 = 0;
    unsigned int const j0 = 0;
    for(unsigned int i = 0 ; i < V.nrblocks() ; i++) {
        for(unsigned int j = 0 ; j < V.ncblocks() ; j++) {
            matrix tm;
            matrix::transpose(tm, V.getblock(i, j));
            U * tp = F.p + (i0+j*B) * F.row_stride + j0/B + i;
            /* Note that the tiny matrix must have space allocated for rows and
             * cols multiples of B, otherwise disaster will occur */
            for(unsigned int k = 0 ; k < B ; k++)
                tp[k*F.row_stride] = tm[k];
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

void blockmatrix::copy_transpose_from_flat(view_t V, const_flat_area F)
{
    using U = matrix::datatype;
    constexpr const unsigned int B = matrix::width;
    unsigned int const i0 = 0;
    unsigned int const j0 = 0;
    for(unsigned int i = 0 ; i < V.nrblocks() ; i++) {
        for(unsigned int j = 0 ; j < V.ncblocks() ; j++) {
            matrix tm;
            const U * tp = F.p + (i0+j*B) * F.row_stride + j0/B + i;
            for(unsigned int k = 0 ; k < B ; k++)
                tm[k] = tp[k*F.row_stride];
            matrix::transpose(V.getblock(i, j), tm);
        }
    }
}

void blockmatrix::copy_from_flat(view_t V, const_flat_area F)
{
    using U = matrix::datatype;
    constexpr const unsigned int B = matrix::width;
    unsigned int const i0 = 0;
    unsigned int const j0 = 0;
    for(unsigned int i = 0 ; i < V.nrblocks() ; i++) {
        for(unsigned int j = 0 ; j < V.ncblocks() ; j++) {
            matrix & tm = V.getblock(i, j);
            const U * tp = F.p + (i0+i*B) * F.row_stride + j0/B + j;
            for(unsigned int k = 0 ; k < B ; k++)
                tm[k] = tp[k*F.row_stride];
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
    ASSERT_ALWAYS(i0 + fnrows <= _nrows);
    ASSERT_ALWAYS(j0 + fncols <= _ncols);
    ASSERT_ALWAYS(j0 % 64 == 0);
    ASSERT_ALWAYS(fncols % 64 == 0);
    set_zero();
    for(unsigned int r = 0 ; r < fnrows ; r++) {
        for(unsigned int g = 0 ; g < fncols ; g+=64) {
            uint64_t v;
            int const rc = fread(&v, sizeof(uint64_t), 1, f);
            ASSERT_ALWAYS(rc == 1);
            (*this)[i0+r][j0+g] = u64_convert_from_little_endian (v);
        }
    }
    fclose(f);
}

#if 1
void blockmatrix::read_transpose_from_flat_file(int i0, int j0, const char * name, unsigned int fnrows, unsigned int fncols)
{
    FILE * f = fopen(name, "rb");
    ASSERT_ALWAYS(f);
    ASSERT_ALWAYS(i0 + fncols <= _nrows);
    ASSERT_ALWAYS(j0 + fnrows <= _ncols);
    ASSERT_ALWAYS(i0 % 64 == 0);
    ASSERT_ALWAYS(j0 % 64 == 0);
    ASSERT_ALWAYS(fncols % 64 == 0);
    set_zero();
    mat64 * tmp = mat64::alloc(fncols/64);
    ASSERT_ALWAYS(tmp);
    for(unsigned int g = 0 ; g < fnrows ; g+=64) {
        /* Fill next block. */
        std::fill_n(tmp, fncols/64, 0);
        for(unsigned int i = 0 ; g + i < fnrows && i < 64 ; i++) {
            for(unsigned int s = 0 ; s < fncols ; s+=64) {
                uint64_t v;
                int const rc = fread(&v, sizeof(uint64_t), 1, f);
                ASSERT_ALWAYS(rc == 1);
                tmp[s/64][i]=v;
            }
        }
        for(unsigned int s = 0 ; s < fncols ; s+=64) {
            mat64_transpose(getblock((s/64), (g/64)), tmp[s/64]);
        }
    }
    mat64::free(tmp, fncols/64);
    fclose(f);
}
#endif

void blockmatrix::write_to_flat_file(const char * name, int i0, int j0, unsigned int fnrows, unsigned int fncols) const
{
    FILE * f = fopen(name, "wb");
    ASSERT_ALWAYS(f);
    ASSERT_ALWAYS(i0 + fnrows <= _nrows);
    ASSERT_ALWAYS(j0 + fncols <= _ncols);
    ASSERT_ALWAYS(j0 % 64 == 0);
    // for our convenience, we allow larger files.
    // ASSERT_ALWAYS(fncols % 64 == 0);
    for(unsigned int r = 0 ; r < fnrows ; r++) {
        for(unsigned int g = 0 ; g < fncols ; g+=64) {
            uint64_t v;
            v = getblock(((i0+r)/64), ((j0+g)/64))[(i0+r)%64];
            v = u64_convert_to_little_endian (v);
            int const rc = fwrite(&v, sizeof(uint64_t), 1, f);
            ASSERT_ALWAYS(rc == 1);
        }
    }
    fclose(f);
}

void blockmatrix::transpose(blockmatrix const & a)
{
    ASSERT_ALWAYS(_nrows == a._ncols);
    ASSERT_ALWAYS(a._nrows == _ncols);
    for(unsigned int i = 0 ; i < a.nrblocks ; i++) {
        mat64_transpose(getblock(i, i), a.getblock(i, i));
        for(unsigned int j = i + 1 ; j < a.ncblocks ; j++) {
            mat64 tmp;
            mat64_transpose(tmp, a.getblock(i, j));
            mat64_transpose(getblock(i, j), a.getblock(j, i));
            getblock(j, i) = tmp;
        }
    }
}

void blockmatrix::swap(blockmatrix & A)
{
    std::swap(A, *this);
}

