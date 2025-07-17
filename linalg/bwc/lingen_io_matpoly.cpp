#include "cado.h" // IWYU pragma: keep

#ifdef LINGEN_BINARY
#include <climits>
#include <cstdlib>
#endif
#include <cstdio>

#include <fstream>
#include <ostream>

#ifndef LINGEN_BINARY
#include <gmp.h>
#endif
#ifdef LINGEN_BINARY
#include <sys/types.h>
#include <unistd.h>     // pread
#endif

#include "lingen_io_matpoly.hpp"
#include "lingen_matpoly_select.hpp"
#include "macros.h"
#include "params.h"
#ifdef HAVE_OPENMP
#include "omp_proxy.h"
#endif



/* {{{ I/O helpers */

/* This is an indication of the number of bytes we read at a time for A
 * (input) and F (output) */
template<bool is_binary>
unsigned int lingen_io_matpoly<is_binary>::block_size = 1 << 20;

template<bool is_binary>
void lingen_io_matpoly<is_binary>::decl_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "io-block-size",
            "chunk size for reading the input or writing the output");
}

template<bool is_binary>
void lingen_io_matpoly<is_binary>::lookup_parameters(cxx_param_list & pl)
{
    param_list_lookup_string(pl, "io-block-size");
}

template<bool is_binary>
void lingen_io_matpoly<is_binary>::interpret_parameters(cxx_param_list & pl)
{
    param_list_parse_uint(pl, "io-block-size", &block_size);
}

/* {{{ matpoly_write
 * writes some of the matpoly data to f, either in ascii or binary
 * format. This can be used to write only part of the data (degrees
 * [k0..k1[). Returns the number of coefficients (i.e., matrices, so at
 * most k1-k0) successfully written, or
 * -1 on error (e.g. when some matrix was only partially written).
 */

#ifndef LINGEN_BINARY
template<>
int lingen_io_matpoly<false>::write(matpoly<false>::arith_hard * ab, std::ostream& os, matpoly<false> const & M, unsigned int k0, unsigned int k1, int ascii, int transpose)
{
    unsigned int const m = transpose ? M.n : M.m;
    unsigned int const n = transpose ? M.m : M.n;
    ASSERT_ALWAYS(k0 == k1 || (k0 < M.get_size() && k1 <= M.get_size()));
    for(unsigned int k = k0 ; k < k1 ; k++) {
        int err = 0;
        int matnb = 0;
        for(unsigned int i = 0 ; !err && i < m ; i++) {
            for(unsigned int j = 0 ; !err && j < n ; j++) {
                matpoly<false>::arith_hard::elt const & x = M.coeff(transpose?j:i,transpose?i:j,k);
                if (ascii) {
                    if (j) err = !(os << " ");
                    if (!err) err = !(ab->cxx_out(os, x));
                } else {
                    err = !ab->write(os, x);
                }
                if (!err) matnb++;
            }
            if (!err && ascii) err = !(os << "\n");
        }
        if (ascii) err = err || !(os << "\n");
        if (err) {
            return (matnb == 0) ? (int) (k - k0) : -1;
        }
    }
    return k1 - k0;
}
#else
template<>
int lingen_io_matpoly<true>::write(matpoly<true>::arith_hard *, std::ostream& os, matpoly<true> const & M, unsigned int k0, unsigned int k1, int ascii, int transpose)
{
    unsigned int const m = M.m;
    unsigned int const n = M.n;
    ASSERT_ALWAYS(k0 == k1 || (k0 < M.get_size() && k1 <= M.get_size()));
    unsigned int const mc = iceildiv(M.m, ULONG_BITS);
    unsigned int const nc = iceildiv(M.n, ULONG_BITS);
    size_t const ulongs_per_mat = transpose ? (M.n * mc) : (M.m * nc);
    std::vector<unsigned long> buf(ulongs_per_mat);
    for(unsigned int k = k0 ; k < k1 ; k++) {
        buf.assign(ulongs_per_mat, 0);
        bool err = false;
        size_t const kq = k / ULONG_BITS;
        size_t const kr = k % ULONG_BITS;
        unsigned long const km = 1UL << kr;
        if (!transpose) {
            for(unsigned int i = 0 ; i < m ; i++) {
                unsigned long * v = &(buf[i * nc]);
                for(unsigned int j = 0 ; j < n ; j++) {
                    unsigned int const jq = j / ULONG_BITS;
                    unsigned int const jr = j % ULONG_BITS;
                    unsigned long const bit = (M.part(i, j)[kq] & km) != 0;
                    v[jq] |= bit << jr;
                }
            }
        } else {
            for(unsigned int j = 0 ; j < n ; j++) {
                unsigned long * v = &(buf[j * mc]);
                for(unsigned int i = 0 ; i < m ; i++) {
                    unsigned int const iq = i / ULONG_BITS;
                    unsigned int const ir = i % ULONG_BITS;
                    unsigned long const bit = (M.part(i, j)[kq] & km) != 0;
                    v[iq] |= bit << ir;
                }
            }
        }
        if (ascii) {
            /* do we have an endian-robust wordsize-robust convention for
             * printing bitstrings in hex ?
             *
             * it's not even clear that we should care -- after all, as
             * long as mksol follows a consistent convention too, we
             * should be fine.
             */
            abort();
        } else {
            err = !(os.write((const char*) buf.data(), ulongs_per_mat * sizeof(unsigned long)));
        }
        if (err) return -1;
    }
    return k1 - k0;
}
#endif
/* }}} */

/* fw must be an array of ofstreams of exactly the same size as the
 * matrix to be written.
 */

#ifdef LINGEN_BINARY
template<>
int lingen_io_matpoly<true>::write_split(matpoly<true>::arith_hard * ab MAYBE_UNUSED, std::vector<std::ofstream> & fw, matpoly<true> const & M, unsigned int k0, unsigned int k1, int ascii)
{
    constexpr const unsigned int splitwidth = 64;
    ASSERT_ALWAYS(k0 == k1 || (k0 < M.get_size() && k1 <= M.get_size()));
    size_t const ulongs_per_mat = splitwidth * splitwidth / ULONG_BITS;
    std::vector<unsigned long> buf(ulongs_per_mat);
    for(unsigned int k = k0 ; k < k1 ; k++) {
        int err = 0;
        int matnb = 0;
        for(unsigned int i = 0 ; !err && i < M.m ; i += splitwidth) {
            for(unsigned int j = 0 ; !err && j < M.n ; j += splitwidth) {
                std::ostream& os = fw[i/splitwidth*M.n/splitwidth+j/splitwidth];
                /* In fact, we _do_ have code that we might want to use
                 * in the b64 layer.
                 */
                if (ascii)
                    abort();
                buf.assign(ulongs_per_mat, 0);
                size_t const kq = k / ULONG_BITS;
                size_t const kr = k % ULONG_BITS;
                for(unsigned int di = 0 ; di < splitwidth ; di++) {
                    unsigned int const ii = i + di;
                    unsigned int const ulongs_per_row = splitwidth / ULONG_BITS;
                    for(unsigned int dj0 = 0 ; dj0 < ulongs_per_row ; dj0++) {
                        for(unsigned int dj = 0 ; dj < ULONG_BITS ; dj++) {
                            unsigned int const jj = j + dj0 * ULONG_BITS + dj;
                            const unsigned long * mij = M.part(ii, jj);
                            unsigned long const bit = (mij[kq] >> kr) & 1;
                            buf[di * ulongs_per_row + dj0] ^= bit << dj;
                        }
                    }
                }
                err = !(os.write((const char *) buf.data(), ulongs_per_mat * sizeof(unsigned long)));
                if (!err) matnb++;
            }
        }
        if (err) {
            return (matnb == 0) ? (int) (k - k0) : -1;
        }
    }
    return k1 - k0;
}
#endif

#ifndef LINGEN_BINARY
template<>
int lingen_io_matpoly<false>::write_split(matpoly<false>::arith_hard * ab MAYBE_UNUSED, std::vector<std::ofstream> & fw, matpoly<false> const & M, unsigned int k0, unsigned int k1, int ascii)
{
    constexpr const unsigned int splitwidth = 1;
    ASSERT_ALWAYS(k0 == k1 || (k0 < M.get_size() && k1 <= M.get_size()));
    for(unsigned int k = k0 ; k < k1 ; k++) {
        int err = 0;
        int matnb = 0;
        for(unsigned int i = 0 ; !err && i < M.m ; i += splitwidth) {
            for(unsigned int j = 0 ; !err && j < M.n ; j += splitwidth) {
                std::ostream& os = fw[i/splitwidth*M.n/splitwidth+j/splitwidth];
                auto const & x = M.coeff(i, j, k);
                if (ascii) {
                    err = !(ab->cxx_out(os, x));
                    if (!err) err = !(os << "\n");
                } else {
                    err = !ab->write(os, x);
                }
                if (!err) matnb++;
            }
        }
        if (err) {
            return (matnb == 0) ? (int) (k - k0) : -1;
        }
    }
    return k1 - k0;
}
#endif

/* }}} */

/* {{{ matpoly_read
 * reads some of the matpoly data from f, either in ascii or binary
 * format. This can be used to parse only part of the data (degrees
 * [k0..k1[, k1 being an upper bound). Returns the number of coefficients
 * (i.e., matrices, so at most k1-k0) successfully read, or
 * -1 on error (e.g. when some matrix was only partially read).
 *
 * Note that the matrix must *not* be in pre-init state. It must have
 * been already allocated.
 */

#ifndef LINGEN_BINARY
template<>
int lingen_io_matpoly<false>::read(matpoly<false>::arith_hard * ab, FILE * f, matpoly<false> & M, unsigned int k0, unsigned int k1, int ascii, int transpose)
{
    ASSERT_ALWAYS(!M.check_pre_init());
    unsigned int const m = transpose ? M.n : M.m;
    unsigned int const n = transpose ? M.m : M.n;
    ASSERT_ALWAYS(k0 == k1 || (k0 < M.get_size() && k1 <= M.get_size()));
    for(unsigned int k = k0 ; k < k1 ; k++) {
        int err = 0;
        int matnb = 0;
        for(unsigned int i = 0 ; !err && i < m ; i++) {
            for(unsigned int j = 0 ; !err && j < n ; j++) {
                auto & x = M.coeff(transpose?j:i,transpose?i:j,k);
                if (ascii) {
                    err = ab->fscan(f, x) == 0;
                } else {
                    err = ab->fread(f, x) < 1;
                }
                if (!err) matnb++;
            }
        }
        if (err) {
            if (matnb == 0) 
                return (int) (k - k0);
            fprintf(stderr, "I/O error after %d,%d coefficients (matrix size %u*%u)\n",
                    (int) (k-k0), matnb, m, n);
            return -1;
        }
    }
    return k1 - k0;
}
#else
static int matpoly_read_inner(matpoly<true>::arith_hard *, FILE * f, matpoly<true> & M, unsigned int k0, unsigned int k1, int ascii, int transpose, off_t base, unsigned int batch = 1)
{
    /* Internally, the dimension of the matrix that are most packed
     * (i.e., rows when in row-major order) are padded to multiples of
     * ULONG_BITS
     */
    unsigned int const m = M.m;
    unsigned int const n = M.n;
    ASSERT_ALWAYS(k0 == k1 || (k0 < M.get_size() && k1 <= M.get_size()));
    unsigned int const mc = iceildiv(M.m, ULONG_BITS);
    unsigned int const nc = iceildiv(M.n, ULONG_BITS);
    size_t const ulongs_per_mat = transpose ? (M.n * mc) : (M.m * nc);
    batch = MIN(batch, k1 - k0);
    std::vector<unsigned long> buf(ulongs_per_mat * batch);
    for(unsigned int k = k0 ; k < k1 ; k+=batch) {
        if (k + batch > k1) batch = k1 - k;
        unsigned int good;
        if (ascii) {
            /* do we have an endian-robust wordsize-robust convention for
             * printing bitstrings in hex ?
             *
             * it's not even clear that we should care -- after all, as long as
             * mksol follows a consistent convention too, we should be fine.
             */
            abort();
        } else {
            int rc;
            if (base < 0) {
                rc = fread(buf.data(), sizeof(unsigned long), ulongs_per_mat * batch, f);
                if (rc != (int) (ulongs_per_mat * batch) && ferror(f))
                    return k - k0;
            } else {
                /* use pread -- good for multithreading */
                size_t const one = ulongs_per_mat * sizeof(unsigned long);
                off_t const off = base + (k - k0) * one;
                ssize_t const r = pread(fileno(f), buf.data(), one * batch, off);
                if (r < 0) rc = 0;
                else rc = r / sizeof(unsigned long);
            }
            good = rc / ulongs_per_mat;
        }
        for(unsigned int b = 0 ; b < good ; b++) {
            size_t const kq = (k + b) / ULONG_BITS;
            size_t const kr = (k + b) % ULONG_BITS;
            if (!transpose) {
                for(unsigned int i = 0 ; i < m ; i++) {
                    unsigned long * v = &(buf[b * ulongs_per_mat + i * nc]);
                    for(unsigned int j = 0 ; j < n ; j++) {
                        unsigned int const jq = j / ULONG_BITS;
                        unsigned int const jr = j % ULONG_BITS;
                        unsigned long const bit = (v[jq] >> jr) & 1;
                        M.part(i, j)[kq] &= ~(1UL << kr);
                        M.part(i, j)[kq] |= bit << kr;
                    }
                }
            } else {
                for(unsigned int j = 0 ; j < n ; j++) {
                    unsigned long * v = &(buf[b * ulongs_per_mat + j * mc]);
                    for(unsigned int i = 0 ; i < m ; i++) {
                        unsigned int const iq = i / ULONG_BITS;
                        unsigned int const ir = i % ULONG_BITS;
                        unsigned long const bit = (v[iq] >> ir) & 1;
                        M.part(i, j)[kq] &= ~(1UL << kr);
                        M.part(i, j)[kq] |= bit << kr;
                    }
                }
            }
        }
        if (good < batch) {
            return k - k0 + good;
        }
    }
    return k1 - k0;
}
template<>
int lingen_io_matpoly<true>::read(matpoly<true>::arith_hard * ab, FILE * f, matpoly<true> & M, unsigned int k0, unsigned int k1, int ascii, int transpose)
{
    int rc = 0;
    if (k0 % ULONG_BITS) {
        unsigned int const fk0 = MIN(k1, k0 + ULONG_BITS - (k0 % ULONG_BITS));
        rc += matpoly_read_inner(ab, f, M, k0, fk0, ascii, transpose, -1);
        if (rc < (int) (fk0 - k0) || fk0 == k1)
            return rc;
        return rc + read(ab, f, M, fk0, k1, ascii, transpose);
    }
    if (k1 % ULONG_BITS) {
        unsigned int const fk1 = MAX(k0, k1 - (k1 % ULONG_BITS));
        if (k0 < fk1) {
            /* recurse and to the bulk of the processing on aligned
             * values */
            rc += read(ab, f, M, k0, fk1, ascii, transpose);
            if (rc < (int) (fk1 - k0))
                return rc;
        }
        return rc + matpoly_read_inner(ab, f, M, fk1, k1, ascii, transpose, -1);
    }

    ASSERT_ALWAYS(!(k0 % ULONG_BITS));
    ASSERT_ALWAYS(!(k1 % ULONG_BITS));

    off_t const pos0 = ftell(f);
    unsigned int const mc = iceildiv(M.m, ULONG_BITS);
    unsigned int const nc = iceildiv(M.n, ULONG_BITS);
    size_t const ulongs_per_mat = transpose ? (M.n * mc) : (M.m * nc);
    size_t const one = ulongs_per_mat * sizeof(unsigned long);

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
    {
#ifdef HAVE_OPENMP
        unsigned int nth = omp_get_num_threads();
#else
        unsigned int const nth = 1;
#endif
        unsigned int const dk = ((k1 - k0)/ULONG_BITS) / nth;
        unsigned int const mk = ((k1 - k0)/ULONG_BITS) % nth;
        unsigned int idx = 0;
#ifdef HAVE_OPENMP
        idx = omp_get_thread_num();
#endif
        unsigned int lk0 = k0 / ULONG_BITS + idx * dk + MIN(idx, mk);
        unsigned int lk1 = lk0 + dk + (idx < mk);
        lk0 *= ULONG_BITS;
        lk1 *= ULONG_BITS;

        off_t const base = pos0 + (lk0 - k0) * one;
        int const my_rc = matpoly_read_inner(ab, f, M, lk0, lk1, ascii, transpose, base, UINT_MAX);

#ifdef HAVE_OPENMP
#pragma omp critical
#endif
        rc += my_rc;
    }

    fseek(f, pos0 + rc * one, SEEK_SET);
    return rc;
}

#endif
/* }}} */

#ifdef LINGEN_BINARY
template struct lingen_io_matpoly<true>;
#else
template struct lingen_io_matpoly<false>;
#endif
