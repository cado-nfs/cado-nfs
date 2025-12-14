#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdint>
#include <climits>

#include <algorithm>
#include <fstream>
#include <istream>
#include <ostream>
#include <ios>
#include <stdexcept>
#include <utility>

#include <sys/stat.h>

#include "cado-endian.h"
#include "fix-endianness.h"
#include "matrix_file.hpp"
#include "omp_proxy.h"
#include "parallelizing_info.hpp"
#include "subdivision.hpp"
#include "macros.h"

int matrix_file::lookup()/*{{{*/
{
    size_t rowcol_sizes[2] = { 0, 0 };
    size_t rowcol_offsets_sizes[2] = { 0, 0 };
    size_t rowcol_coeffs_sizes[2] = { 0, 0 };
    lookup_diagnostic.clear();

    for(int d = 0 ; d < 2 ; d++) {
        struct stat sbuf[1];
        std::string t;
        int rc;

        t = mfile + dotrowcols(d);
        rc = stat(t.c_str(), sbuf);
        has_rowcols[d] = (rc == 0);
        if (rc == 0)
            rowcol_sizes[d] = sbuf->st_size;

        t = mfile + dotrowcol_offsets(d);
        rc = stat(t.c_str(), sbuf);
        has_rowcol_offsets[d] = (rc == 0);
        if (rc == 0)
            rowcol_offsets_sizes[d] = sbuf->st_size;

        t = mfile + dotrowcol_coeffs(d);
        rc = stat(t.c_str(), sbuf);
        has_rowcol_coeffs[d] = (rc == 0);
        if (rc == 0)
            rowcol_coeffs_sizes[d] = sbuf->st_size;
    }

    /* Now do some checks */
    if (!has_rowcols[0] && !has_rowcols[1]) {
        lookup_diagnostic = "could not find any of .rows and .cols";
        return 0;
    }

    if (!has_rowcol_offsets[0] || !has_rowcol_offsets[1]) {
        lookup_diagnostic = "could not find both .row_offsets and .col_offsets";
        return 0;
    }

    if (has_rowcols[0] && has_rowcols[1]) {
        if (rowcol_sizes[0] != rowcol_sizes[1]) {
            lookup_diagnostic = "inconsistent sizes in .rows and .cols";
            return 0;
        }
    }

    /* Let's compute a default setting which we think makes sense for
     * this sepcific set of files. It doesn't mean that we'll read the
     * file exactly in this way. Reading in either direction remains
     * possible.
     */
    direction = !has_rowcols[0];

    withcoeffs = has_rowcol_coeffs[0] || has_rowcol_coeffs[1];

    for(int d = 0 ; d < 2 ; d++) {
        nrowcols[d] = rowcol_offsets_sizes[d] / sizeof(uint32_t);

        if (has_rowcols[d]) ncoeffs = rowcol_sizes[d] / sizeof(uint32_t);

        if (has_rowcol_coeffs[d] && !has_rowcols[d]) {
            lookup_diagnostic = dotrowcol_coeffs(d) + " without " + dotrowcols(d);
            return 0;
        }

        if (has_rowcol_coeffs[d] && rowcol_coeffs_sizes[d] != rowcol_sizes[d]) {
            lookup_diagnostic = dotrowcol_coeffs(d) + " has inconsistent size";
            return 0;
        }
    }

    return 1;
}/*}}}*/

/* {{{ helpers around pi_bcast and friends.
 * Have I coded this somewhere else already ? */
template<typename T> struct bwc_pi_type
{ static const pi_datatype_ptr value; };
template<> const pi_datatype_ptr bwc_pi_type<int>::value = BWC_PI_INT;
template<> const pi_datatype_ptr bwc_pi_type<unsigned int>::value = BWC_PI_UNSIGNED;
template<> const pi_datatype_ptr bwc_pi_type<unsigned long>::value = BWC_PI_UNSIGNED_LONG;
template<> const pi_datatype_ptr bwc_pi_type<unsigned long long>::value = BWC_PI_UNSIGNED_LONG_LONG;

template<typename T> static void bcast(parallelizing_info_ptr pi, T & s, unsigned int jrank, unsigned int trank)
{
    pi_bcast(&s, 1, bwc_pi_type<T>::value, jrank, trank, pi->m);
}
/* this is a full specialization */
template<> void bcast<bool>(parallelizing_info_ptr pi, bool & s, unsigned int jrank, unsigned int trank)
{
    int c = s;
    pi_bcast(&c, 1, BWC_PI_INT, jrank, trank, pi->m);
    s = c;
}
/* below are two overloads */
template<typename T, size_t N> 
static void
bcast(parallelizing_info_ptr pi, std::array<T, N> & s, unsigned int jrank, unsigned int trank)
{
    pi_bcast(s.data(), N, bwc_pi_type<T>::value, jrank, trank, pi->m);
}
static void bcast(parallelizing_info_ptr pi, std::string & s, unsigned int jrank, unsigned int trank)
{
    size_t sz = s.size();
    pi_bcast(&sz, 1, BWC_PI_SIZE_T, jrank, trank, pi->m);
    std::vector<char> foo(sz, ' ');
    if (pi->m->jrank == jrank && pi->m->trank == trank)
        std::ranges::copy(s, foo.begin());
    pi_bcast(foo.data(), sz, BWC_PI_BYTE, jrank, trank, pi->m);
    s = std::string(foo.begin(), foo.end());
}
template<typename T> static void allreduce(parallelizing_info_ptr pi, T & s, pi_op_ptr op)
{
    pi_allreduce(nullptr, &s, 1, bwc_pi_type<T>::value, op, pi->m);
}
/* }}} */

int matrix_file::lookup(parallelizing_info_ptr pi)/*{{{*/
{
    /* A priori, each thread / node has the same filename in memory. We
     * rely on that.
     * All threads on a given node have access to the same filesystem, so
     * they should be able to see the same file.
     * What is likely to differ, however, is how the different nodes see
     * the files. It could well be that only one node can see them.
     */
    
    /* begin with local lookups */
    int c = 0;
    if (pi->m->trank == 0)
        c = lookup();

    unsigned int j = pi->m->jrank;
    if (!c) j = UINT_MAX;
    allreduce(pi, j, BWC_PI_MIN);
    allreduce(pi, c, BWC_PI_MAX);

    if (c == 0) {
        /* the local lookup failed *everywhere*. Let's return the
         * diagnostic found by job 0, thread 0
         */
        bcast(pi, lookup_diagnostic, 0, 0);
        return 0;
    }

    /* We're only sharing this. The rest (has_rowcol etc) is best kept
     * local.
     */
    bcast(pi, nrowcols, j, 0);
    bcast(pi, ncoeffs, j, 0);
    bcast(pi, withcoeffs, j, 0);
    bcast(pi, direction, j, 0);

    return 1;
}/*}}}*/

static inline void vanished(std::string & filename)
{
    throw std::runtime_error(filename + " has vanished");
}

static inline void read32(std::istream& is, uint32_t * x, size_t n = 1)
{
    is.read((char *) x, n * sizeof(uint32_t));
    for(size_t i = 0 ; i < n ; i++) x[i] = bfix32(x[i]);
}

static inline uint32_t read32(std::istream& is)
{
    uint32_t x;
    read32(is, &x);
    return x;
}

static inline void write32(std::ostream& os, uint32_t const * x, size_t n = 1)
{
#ifdef CADO_LITTLE_ENDIAN
    os.write((char const *) x, n * sizeof(uint32_t));
#else
    for(auto x : *this) {
        x = bfix32(x);
        os.write((char const *) &x, sizeof(uint32_t));
    }
#endif
}

static inline void write32(std::ostream& os, uint32_t const & x)
{
    write32(os, &x);
}


template<bool withcoeffs>
static void inner_loop(uint32_t * p, uint32_t offset, unsigned int i0, unsigned int i1, std::ifstream & RW, std::ifstream& DATA, std::ifstream& COEFFS)
{
    /* TODO: display progress */
    p += i0 + offset * (1 + withcoeffs);
    for(unsigned int i = i0 ; i < i1 ; i++) {
        uint32_t const next_offset = read32(RW);
        unsigned int const w = *p++ = next_offset - offset;
        offset = next_offset;
        auto * q = p;
        if (!withcoeffs) {
            read32(DATA, p, w);
            if (!std::ranges::is_sorted(p, p+w))
                std::ranges::sort(p, p + w);
            p += w;
        } else {
            for(unsigned int k = 0 ; k < w ; k++) {
                *p++ = read32(DATA);
                *p++ = read32(COEFFS);
            }
            auto * h = reinterpret_cast<std::pair<uint32_t, int32_t> *>(q);
            if (!std::ranges::is_sorted(h, h+w)) {
                std::ranges::sort(h, h + w);
            }
        }
    }
}

void matrix_file::read(int direction, std::string const & sanity_check_vector MAYBE_UNUSED)
{
    /* if direction=0 we need .rows .row_offsets .row_coeffs */
    /* if direction=1 we need .cols .col_offsets .col_coeffs */

    int const d = direction;
    
    require_lookup();

    if (has_rowcols[d] && has_rowcol_offsets[d] && (!withcoeffs || has_rowcol_coeffs[d])) {
        /* We can read this data as is */
        std::string rwname = mfile + dotrowcol_offsets(d);
        std::string dataname = mfile + dotrowcols(d);
        std::string coeffsname = mfile + dotrowcol_coeffs(d);

        std::vector<uint32_t> & p(*this);
        p.assign(nrowcols[d] + ncoeffs * (1 + withcoeffs), 0);

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
        {
            subdivision const D(nrowcols[d], omp_get_num_threads());

            unsigned int const i0 = D.nth_block_start(omp_get_thread_num());
            unsigned int const i1 = D.nth_block_end(omp_get_thread_num());

            std::ifstream RW(rwname, std::ios::in | std::ios::binary);
            if (!RW) vanished(rwname);
            uint32_t z;
            if (i0) {
                RW.seekg((i0 - 1) * sizeof(uint32_t));
                z = read32(RW);
            } else {
                z = 0;
            }

            std::ifstream DATA(dataname, std::ios::in | std::ios::binary);
            if (!DATA) vanished(dataname);
            DATA.seekg(z * sizeof(uint32_t));

            std::ifstream COEFFS;
            if (withcoeffs) {
                COEFFS.open(coeffsname, std::ios::in | std::ios::binary);
                if (!COEFFS) vanished(coeffsname);
                COEFFS.seekg(z * sizeof(uint32_t));
            }

            if (withcoeffs) {
                inner_loop<true>(data(), z, i0, i1, RW, DATA, COEFFS);
            } else {
                inner_loop<false>(data(), z, i0, i1, RW, DATA, COEFFS);
            }
        }
    } else {
        int const other = d ^ 1;
        /* We will read the matrix in the wrong order, but use the known
         * offsets (in the correct direction) in order to place data
         * correctly.
         *
         * Let's suppose that d = direction = 1, and other = 0: we read
         * the .rows file, but we store in column-major form. It also
         * works the other way around if we so wish.
         */
        if (!(has_rowcols[other] && has_rowcol_offsets[d] && has_rowcol_offsets[other] && (!withcoeffs || has_rowcol_coeffs[other]))) {
            throw std::runtime_error(cmdline_filename + " cannot be read (some or all of the required files are missing for a read in " + rowcol[d] + "-major layout)");
        }
        std::string cwname = mfile + dotrowcol_offsets(d);
        std::string rwname = mfile + dotrowcol_offsets(other);
        std::string const dataname = mfile + dotrowcols(other);

        /* We need to load all offsets in the final direction beforehand
         * because it's a transposed read.
         */
        std::vector<uint32_t> offsets;
        offsets.assign(1 + nrowcols[d], 0);
        std::ifstream CW(cwname);
        if (!CW) vanished(cwname);
        read32(CW, offsets.data() + 1, nrowcols[d]);

        std::ifstream RW(rwname, std::ios::in | std::ios::binary);
        if (!RW) vanished(rwname);
        uint32_t z = 0;

        /* This is a relatively dumb, single-threaded version. How we
         * could expand it to multi-thread isn't totally clear
         */

        std::ifstream DATA(dataname, std::ios::in | std::ios::binary);
        if (!DATA) throw std::runtime_error(dataname + " has vanished");
        assign(nrowcols[d] + ncoeffs * (1 + withcoeffs), 0);
        uint32_t * p = data();

        std::ifstream COEFFS;
        if (withcoeffs) {
            std::string coeffsname = mfile + dotrowcol_coeffs(other);
            COEFFS.open(coeffsname, std::ios::in | std::ios::binary);
            if (!COEFFS) vanished(coeffsname);
        }

        /* First, adjust offsets to reflect more exactly the position in
         * he data array, with the line length and (when present)
         * coefficient space both taken into account. We write the
         * beginning-of-row markers to memory at the same time.
         */
        for(unsigned int j = 0 ; j < nrowcols[d] ; j++) {
            uint32_t const w = offsets[j + 1] - offsets[j];
            offsets[j] = j + offsets[j] * (1 + withcoeffs);
            p[offsets[j]++] = w;
        }

        /* This loop cannot be parallelied since it is synchronized with
         * the offsets[] vector
         */
        for(unsigned int i = 0 ; i < nrowcols[other] ; i++) {
            uint32_t const next_z = read32(RW);
            uint32_t const w = next_z - z;
            z = next_z;
            for(uint32_t k = 0 ; k < w ; k++) {
                uint32_t const j = read32(DATA);
                ASSERT_ALWAYS(j < nrowcols[d]);
                p[offsets[j]++] = i;
                if (withcoeffs)
                    p[offsets[j]++] = read32(COEFFS);
            }
        }
    }
}


void matrix_file::dump_data(std::ostream& os) const
{
    for(auto p = begin() ; p != end() ; ) {
        uint32_t const w = *p++;
        for(unsigned int i = 0 ; i < w ; i++, p += 1 + withcoeffs) {
            write32(os, bfix32(p[0]));
        }
    }
}

void matrix_file::dump_offsets(std::ostream& os) const
{
    uint32_t z = 0;
    for(auto p = begin() ; p != end() ; ) {
        uint32_t const w = *p++;
        z += w;
        write32(os, bfix32(z));
        p += w;
        if (withcoeffs) p += w;
    }
}

void matrix_file::dump_coeffs(std::ostream& os) const
{
    if (!withcoeffs) return;

    for(auto p = begin() ; p != end() ; ) {
        uint32_t const w = *p++;
        for(unsigned int i = 0 ; i < w ; i++, p += 1 + withcoeffs) {
            write32(os, bfix32(p[1]));
        }
    }
}

void matrix_file::dump_mixed(std::ostream& os) const
{
    write32(os, data(), size());
}

