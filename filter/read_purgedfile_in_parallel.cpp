#include "cado.h" // IWYU pragma: keep

#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "fmt/base.h"

#include "macros.h"
#include "merge_heap.h"
#include "merge_replay_matrix.h"
#include "omp_proxy.h"
#include "read_purgedfile_in_parallel.hpp"
#include "sparse.h"
#include "timing.h"
#include "typedefs.h"
#include "utils_cxx.hpp"

namespace {

/* This file sticks to stdio and to bare char pointers on purpose. The
 * C++ stream machinery was tried here and is not usable at this scale:
 * merely constructing an istringstream per row divides the throughput by
 * about ten, and imbuing it with a locale that treats ',' and ':' as
 * separators divides it by another four.
 */
using unique_file = std::unique_ptr<FILE, delete_FILE>;

/* Value of each character as a hexadecimal digit, or 255 if it is not
 * one. */
constexpr std::array<unsigned char, 256> hex_digits = []() {
    std::array<unsigned char, 256> T {};
    T.fill(255);
    for (unsigned char i = 0; i < 10; i++)
        T[std::size_t('0') + i] = i;
    for (unsigned char i = 0; i < 6; i++)
        T[std::size_t('a') + i] = T[std::size_t('A') + i] = 10 + i;
    return T;
}();

/* Parse a hexadecimal integer at p, and return a pointer to the first
 * character that is not a hexadecimal digit. Same thing as
 *      char * q; x = strtoul(p, &q, 16); return q;
 *
 * This much plumbing for a strtoul call deserves a justification, so
 * here is the measurement (34MB purged file, 28-bit indices, gcc 16.2,
 * glibc 2.43). On the parsing alone: 2000MB/s here, 1500MB/s with
 * std::from_chars, 260MB/s with strtoul. That is not lost in the noise
 * of the surrounding I/O: reading the same file end to end goes at
 * 370MB/s here and 160MB/s with strtoul on one thread, 1300MB/s versus
 * 560MB/s on four.
 */
inline char const * parse_hex(index_t & x, char const * p)
{
    index_t v = 0;
    for (unsigned char c; (c = hex_digits[(unsigned char) *p++]) != 255;)
        v = v * 16 + c;
    x = v;
    return p - 1;
}

/* Read the rows of a purged file with several threads.
 *
 * The file is cut in one contiguous byte range per thread. Thread i
 * reads whole lines: it begins with the first line that starts at or
 * after the beginning of its range, and stops when it reaches the
 * beginning of the range of thread i+1. The adjustment of the cut
 * positions to line boundaries is done by the threads themselves, so
 * that the cut positions have to be exchanged before any thread can tell
 * where it must stop.
 */
class parallel_reader {
    std::string filename;
    uint64_t skip;
    off_t endpos;

    /* Number of threads that we _really_ have, which is only known once
     * we are inside the parallel region. */
    unsigned int nthreads = 0;

    /* nthreads+1 cut positions. spos[i] is adjusted by thread i so that
     * it sits at the beginning of a line. */
    std::vector<off_t> spos;
    std::vector<std::size_t> rows_per_thread;

    /* Progress report. These are only ever touched from within the
     * critical section in read_range(), or outside the parallel region.
     */
    double tt0 = wct_seconds();
    std::size_t nrows_read = 0;
    std::size_t bytes_read = 0;
    std::size_t next_report = 1024;

    off_t file_size() const
    {
        unique_file const f = fopen_helper(filename, "r");
        int const rc = fseek(f.get(), 0, SEEK_END);
        DIE_ERRNO_DIAG(rc < 0, "fseek(%s)", filename.c_str());
        return ftell(f.get());
    }

    void print_report()
    {
        double const dt = wct_seconds() - tt0;
        fmt::print("# Read {} relations in {:.1f}s"
                   " -- {:.1f} MB/s -- {:.1f} rels/s\n",
                   nrows_read, dt,
                   double(bytes_read >> 20U) / dt, double(nrows_read) / dt);
        fflush(stdout);
        /* next time we report, we'll have read twice as many rows */
        next_report *= 2;
    }

    /* Cut the file in n pieces. Called from an omp single construct. */
    void setup(unsigned int n)
    {
        nthreads = n;
        fmt::print(stderr, "# {}: Doing I/O with {} threads\n",
                   filename, nthreads);
        rows_per_thread.assign(nthreads, 0);
        spos.resize(nthreads + 1);
        for (unsigned int i = 0; i < nthreads; i++)
            spos[i] = (endpos * i) / nthreads;
        spos[nthreads] = endpos;
    }

    std::vector<typerow_t *> read_range(FILE * f, off_t bytes_to_read);

  public:
    parallel_reader(std::string filename, uint64_t skip)
        : filename(std::move(filename))
        , skip(skip)
        , endpos(file_size())
    {
    }

    std::size_t read(filter_matrix_t * mat);
};

std::vector<typerow_t *> parallel_reader::read_range(FILE * f,
                                                        off_t bytes_to_read)
{
    std::vector<typerow_t *> rows;

    /* These two are reused from one row to the next. This avoids
     * frequent roundtrips to the malloc layer. */
    std::array<char, 4096> line;
    std::vector<typerow_t> primes;

    /* We add our contribution to the shared counters only every so
     * often. The step is chosen so that the threads collectively report
     * about twice per doubling of the global row count. */
    std::size_t local_next_report = 256;
    std::size_t reported_rows = 0;
    std::size_t reported_bytes = 0;

    long const start = ftell(f);

    for (long pos; (pos = ftell(f)) - start < bytes_to_read;) {
        if (fgets(line.data(), line.size(), f) == nullptr)
            break;

        if (line[0] == '#')
            continue;

        std::size_t const len = strlen(line.data());
        /* otherwise 4096 is not enough ! */
        ASSERT_ALWAYS(len < line.size() - 1);
        /* Where the interesting part of the line stops: at the
         * end-of-line delimiter, or at the end of the string if the last
         * line of the file happens not to be terminated. (A line that is
         * too long for our buffer would look the same, but the assertion
         * above has ruled that out.) */
        char const * eol = line.data() + len;
        if (len > 0 && eol[-1] == '\n')
            eol--;

        primes.clear();
        /* primes[0] is a placeholder: we will eventually store the
         * number of ideals of the row there. */
        typerow_t marker;
        setCell(&marker, 0, 0, 0);
        primes.push_back(marker);

        /* A row reads "a,b:i0,i1,...,ik". Only the hexadecimal ideal
         * indices after the colon are of interest to us. */
        char const * p = line.data();
        for (; *p && *p != ':'; p++)
            ;
        /* note the comparison: p may well jump _over_ eol, e.g. if there
         * is no colon at all in the line. We must not run past the end
         * of the buffer in that case. */
        for (; p < eol;) {
            /* skip the ':' or the ',' */
            p++;
            index_t x;
            p = parse_hex(x, p);
            if (x < skip)
                continue;
            typerow_t c;
            setCell(&c, 0, x, 1);
            primes.push_back(c);
        }

        /* Sort the ideals, and collapse the repeated ones. */
#ifdef FOR_DL
        std::sort(primes.begin() + 1, primes.end(),
                  [](typerow_t const & a, typerow_t const & b) {
                      return a.id < b.id;
                  });
#else
        std::sort(primes.begin() + 1, primes.end());
#endif
        std::size_t j = 1;
        for (std::size_t i = 1; i != primes.size();) {
            primes[j] = primes[i];
            std::size_t k = i + 1;
#ifdef FOR_DL
            /* accumulate the exponents of the ideals that are equal */
            for (; k != primes.size() && primes[k].id == primes[j].id; k++)
                primes[j].e += primes[k].e;
            j++;
#else
            /* we work modulo 2, so an ideal that appears an even number
             * of times simply goes away */
            for (; k != primes.size() && primes[k] == primes[j]; k++)
                ;
            j += (k - i) & 1;
#endif
            i = k;
        }
        primes.resize(j);

        /* Pay attention to the placeholder ! and fill it in, too. */
        index_t const n = primes.size() - 1;
        setCell(primes.data(), 0, n, 0);

        /* The 0 below must eventually become the row index, but we can't
         * write it right now: we do not know yet how many rows the
         * threads before us have read. See read(). */
        typerow_t * row = heap_alloc_row(0, n);
        compressRow(row, primes.data(), n);
        rows.push_back(row);

        if (rows.size() >= local_next_report) {
            std::size_t const local_bytes = pos - start;
#pragma omp critical
            {
                nrows_read += rows.size() - reported_rows;
                bytes_read += local_bytes - reported_bytes;
                reported_rows = rows.size();
                reported_bytes = local_bytes;
                if (nrows_read >= next_report)
                    print_report();
                local_next_report += next_report / 2 / nthreads;
            }
        }
    }

#pragma omp critical
    {
        std::size_t const local_bytes = ftell(f) - start;
        nrows_read += rows.size() - reported_rows;
        bytes_read += local_bytes - reported_bytes;
    }

    return rows;
}

std::size_t parallel_reader::read(filter_matrix_t * mat)
{
    unsigned int const max_threads =
        std::min<unsigned int>(omp_get_max_threads(), MAX_IO_THREADS);

    /* All threads get their private reading head. */
#pragma omp parallel num_threads(max_threads)
    {
        /* num_threads() is only an upper bound. When dynamic adjustment
         * of the number of threads is enabled (OMP_DYNAMIC=true, which
         * our test suite sets), the runtime may hand us a smaller team,
         * and libgomp routinely hands us a single thread. We must
         * therefore cut the file in as many pieces as we have threads
         * *for real*. Deciding on the cut before the parallel region
         * would leave the pieces of the threads that we did not get
         * entirely unread, and merge would then see only part of the
         * rows.
         */
#pragma omp single
        setup(omp_get_num_threads());
        /* the omp single construct above ends with an implicit barrier,
         * so that spos[] is set and visible to everyone here */

        int const i = omp_get_thread_num();

        /* the buffer must outlive the FILE that uses it, hence the
         * declaration order */
        std::array<char, 1U << 16U> buffer;
        unique_file const f = fopen_helper(filename, "r");
        setbuffer(f.get(), buffer.data(), buffer.size());

        int const rc = fseek(f.get(), spos[i], SEEK_SET);
        DIE_ERRNO_DIAG(rc < 0, "fseek(%s)", filename.c_str());

        /* Except when we're at the beginning of the stream, the line
         * that straddles the beginning of our range belongs to the
         * previous thread: read until we get a newline. Note that we
         * must stop at end of file too. A file whose last line is not
         * terminated would otherwise have us spin forever, since fgetc()
         * keeps returning EOF. */
        if (i > 0) {
            for (int c; (c = fgetc(f.get())) != '\n' && c != EOF;)
                ;
        }
        spos[i] = ftell(f.get());

#pragma omp barrier

        /* Cut positions are only ever moved forwards, and they start out
         * sorted, so this holds. read_range() would read the same rows
         * twice, or skip some, if it did not. */
        ASSERT_ALWAYS(spos[i] <= spos[i + 1]);

        std::vector<typerow_t *> const rows =
            read_range(f.get(), spos[i + 1] - spos[i]);
        rows_per_thread[i] = rows.size();

#pragma omp barrier

        /* Now that we know how many rows the threads before us have
         * read, we can renumber our rows and store them in order in
         * mat->rows. */
        std::size_t index = 0;
        for (int j = 0; j < i; j++)
            index += rows_per_thread[j];
        for (typerow_t * row: rows) {
            rowCell(row - 1, 0) = index;
            mat->rows[index] = row;
            index++;
        }
    }

    print_report();

    return nrows_read;
}

} // namespace

uint64_t read_purgedfile_in_parallel(filter_matrix_t * mat,
                                     std::string const & filename)
{
    return parallel_reader(filename, mat->skip).read(mat);
}
