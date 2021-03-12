#include "cado.h"
#include <stdio.h>
#include <istream>
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <locale>
#include "omp_proxy.h"
#include "timing.h"
#include "merge_heap.h"
#include "sparse.h"
#include "read_purgedfile_in_parallel.h"
#include "fmt/printf.h"
#include "fmt/format.h"


struct global_tracking {
    std::vector < uint64_t > rows_per_thread;
    std::vector < off_t > spos_tab;
    size_t bytes = 0;
    size_t next_report = 1024;
    double tt0;
    size_t nrows = 0;
    size_t nthreads() const { return rows_per_thread.size(); }
    global_tracking(int nthreads, off_t endpos) : rows_per_thread(nthreads,0) {
        tt0 = wct_seconds();
        for (int i = 0; i < nthreads; i++)
            spos_tab.push_back((endpos * i) / nthreads);
        spos_tab.push_back(endpos);
    }
    void print_report() {
        double dt = wct_seconds() - tt0;
        fmt::printf
            ("# Read %zu relations in %.1fs -- %.1f MB/s -- %.1f rels/s\n",
             nrows, dt, (bytes >> 20) / dt,
             nrows / dt);

        next_report *= 2;
    }
};

template<typename T> inline T hacked_strtoul16(char * & p)/*{{{*/
{
    /* functionally equivalent to:
         char * q;
         T x = strtoul(p, &q, 16);
         p = q;
         return x;
     */
    static const unsigned char ugly[256] = {
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        0,   1,   2,   3,   4,   5,   6,   7,  
        8,   9,   255, 255, 255, 255, 255, 255,
        255, 10,  11,  12,  13,  14,  15,  255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 10,  11,  12,  13,  14,  15,  255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255 
    };
    T x = 0;
    for(int c ; (c = ugly[(int) (unsigned char) *p++]) != 255 ; x = x*16+c) ;
    p--;
    return x;
}
/*}}}*/

#ifdef FOR_DL
template<>
struct std::less<ideal_merge_t>
{
    inline bool operator() (ideal_merge_t const &a, ideal_merge_t const &b) const {
        return a.id < b.id;
    }
};
#endif

struct local_row_reader {
    global_tracking & G;
    size_t local_next_report = 256;
    size_t local_nrows_at_last_report = 0;
    size_t local_bytes_at_last_report = 0;
    std::vector < typerow_t * > local_rows;
    filter_matrix_t * mat;
    int i;
    local_row_reader(global_tracking & G, filter_matrix_t * mat) : G(G), mat(mat) {
        i = omp_get_thread_num();
    }
    void read(std::istream & fi, off_t bytes_to_read, global_tracking & G) {

        std::string s;

        off_t start = fi.tellg();

        for (off_t pos ; ((pos = fi.tellg()) - start)  < bytes_to_read ; ) {
            std::vector < typerow_t > primes;
            /* Insert a temporary marker. We'll use it for storing the
             * size, eventually */
            typerow_t zz;
            setCell(&zz, 0, 0, 0);
            primes.push_back(zz);
            s.clear();
            {
                /* this is going to be a bit ugly, I know */
                std::getline(fi, s);
                if (s[0] == '#') continue;

                /* see "BAD IDEAS FOR PARSING LOOP" below for things that
                 * I tried and didn't play out well.  */
                char * p = &s[0];
                char * z = p + s.size();
                for( ; *p && *p != ':' ; p++);
                for( ; p++ != z ; ) {
                    index_t x = hacked_strtoul16<index_t>(p);
                    if (x < mat->skip)
                        continue;
                    typerow_t xx;
                    setCell(&xx, 0, x, 1);
                    primes.push_back(xx);
                }
            }

            std::sort(primes.begin() + 1, primes.end(), std::less<typerow_t>());

            auto jt = primes.begin() + 1;
            for (auto it = primes.begin() + 1; it != primes.end();) {
                *jt = *it;
                auto kt = it;
                ++kt;
#ifdef FOR_DL
                for (; kt != primes.end() && kt->id == jt->id; ++kt)
                    jt->e += kt->e;
                jt++;
#else
                for (; kt != primes.end() && *kt == *jt; ++kt);
                jt += ((kt - it) & 1);
#endif
                it = kt;
            }
            primes.erase(jt, primes.end());

            /* Pay attention to the special marker ! and update it, too. */
            unsigned int z = primes.size() - 1;
            setCell(primes, 0, z, 0);

            /* 0 here must eventually become the row index, but we can't
             * write it right now. We'll do so later on.  */
            typerow_t *newrow = heap_alloc_row(0, z);
            compressRow(newrow, &primes[0], z);
            local_rows.push_back(newrow);

            /* At this point we should consider reporting. */
            size_t local_bytes = pos - start;
            if (local_rows.size() >= local_next_report) {
#pragma omp critical
                {
                    G.nrows += local_rows.size() - local_nrows_at_last_report;
                    G.bytes += local_bytes - local_bytes_at_last_report;
                    local_nrows_at_last_report = local_rows.size();
                    local_bytes_at_last_report = local_bytes;
                    if (G.nrows >= G.next_report)
                        G.print_report();
                    local_next_report += G.next_report / 2 / G.nthreads();
                }
            }
        }
#pragma omp critical
        {
            G.nrows += local_rows.size() - local_nrows_at_last_report;
            size_t local_bytes = fi.tellg() - start;
            G.bytes += local_bytes - local_bytes_at_last_report;
            G.rows_per_thread[i] = local_rows.size();
        }
    }
};

uint64_t read_purgedfile_in_parallel(filter_matrix_t * mat,
				     const char *filename)
{
    off_t endpos;

    {
	std::ifstream f(filename);
	endpos = f.seekg(0, std::ios_base::end).tellg();
    }

    /* Find accurate starting positions for everyone */
    unsigned int nthreads = omp_get_max_threads();

    /* cap the number of I/O threads */
    if (nthreads > MAX_IO_THREADS)
	nthreads = MAX_IO_THREADS;

    fmt::fprintf(stderr, "# %s: Doing I/O with %u threads\n", filename,
		 nthreads);

    global_tracking G(nthreads, endpos);

    /* All threads get their private reading head. */
#pragma omp parallel num_threads(nthreads)
    {
        int i = omp_get_thread_num();
        std::ifstream fi;
        std::vector<char> buffer(1<<16);
        fi.rdbuf()->pubsetbuf(&buffer[0], buffer.size());

        fi.open(filename, std::ios_base::in);

        fi.seekg(G.spos_tab[i], std::ios_base::beg);
        /* Except when we're at the beginning of the stream, read until
         * we get a newline */
        for (; i > 0 && fi.get() != '\n';);
        G.spos_tab[i] = fi.tellg();

#pragma omp barrier

        off_t bytes_to_read = G.spos_tab[i + 1] - G.spos_tab[i];
        local_row_reader L(G, mat);
        L.read(fi, bytes_to_read, G);

#pragma omp barrier

        /* now we want to collect the number of rows for everyone,
         * renumber all relations, and store them in order in mat->rows.
         */
        uint64_t index = 0;
        for(int j = 0 ; j < i ; j++)
            index += G.rows_per_thread[j];
        for (auto & r : L.local_rows) {
            rowCell((r - 1), 0) = index;
            mat->rows[index] = r;
            index++;
        }
    }
    G.print_report();

    return G.nrows;
}


/* BAD IDEAS FOR PARSING LOOP */

/* In the critical parsing loop, there are several ways to parse the
 * string s and make a relation out of it.  The whole point of reading in
 * parallel is that we're going to have several threads that keep the I/O
 * layer happy, so that what happens CPU-bound shouldn't be *that*
 * important.
 *
 * However, we do care about lock contention. And apparently, there are
 * several bad ideas that one can come up with, which eventually destroy
 * performance.
 */

#if 0
    /* The #1 outrageously bad idea is this one.  Doing _just_ this
     * initialization of an istringstream, and nothing else, divides the
     * 16-thread reading performance by about 10 (from ~ 900MB/s to
     * 90MB/s from beegfs, with really nothing beyond progress printing
     * down the line.
     * */
    std::istringstream ss(s);
    /* This is another terrible idea, with catastrophic impact on
     * performance. On the same metric as above, we drop from 90MB/s to
     * 20MB/s.
     */
    ss.imbue(std::locale(std::locale(), new csv_reader()));
#if 0   /* corresponding companion structure */
        // https://stackoverflow.com/questions/1894886/parsing-a-comma-delimited-stdstring
        struct csv_reader: std::ctype<char> {
            csv_reader(): std::ctype<char>(get_table()) {}
            static std::ctype_base::mask const* get_table() {
                static std::vector<std::ctype_base::mask> rc(table_size, std::ctype_base::mask());

                rc[','] = std::ctype_base::space;
                rc['\n'] = std::ctype_base::space;
                rc[' '] = std::ctype_base::space;
                rc[':'] = std::ctype_base::space;
                return &rc[0];
            }
        }; 
#endif

    /* This first step of parsing over a,b doesn't degrade performance
     * too much.  */
    ss >> std::hex;
    {
        std::string a,b;
        ss >> a >> b;
    }
    /* this for loop, however, seems to tax the CPU somewhat (so yet
     * another small degradation, but we've already reached a pretty
     * degraded situation at this point).  I also see many futex calls, I
     * suspect that some libstc++ locking is at stake here. Maybe it's a
     * stupid thing with reference counting of a locale that goes with a
     * shared ptr, or something like this.
     */
    for( ; ss.good() ; ) {
        index_t x;
        ss >> x;
        typerow_t xx;
        setCell(&xx, 0, x, 1);
        primes.push_back(xx);
    }
#endif

#if 0
    /* it is also worth noting that hacked_strtoul16 is moderately, but
     * still measurably slower than its strtoul-based equivalent
     */
#endif

#if 0
/* We also see many lseek(fd, 0, SEEK_CUR) in strace that correspond to
 * our tellg() calls. We tried to kill them as follows, but it isn't
 * great, and raises double buffering questions that I don't really want
 * to go into.
 */
class countbuf : public std::streambuf {
    std::streambuf & sbuf;
    std::streamsize size = 0;
    char ch;
public:
    countbuf(std::streambuf & sbuf): sbuf(sbuf) {}
protected:
    virtual int underflow()
    {
        size++;
        ch = sbuf.sbumpc();
        setg(&ch, &ch, &ch+1);
        return ch;
    }
public:
    std::streamsize count() { return this->size; }
};
// this would be used as follows, the ifstream fi being replaced by fi2
// as follows.
        // countbuf cc(*fi.rdbuf());
        // std::istream fi2(&cc);

#endif

