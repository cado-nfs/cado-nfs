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


#if 0
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
#endif


struct local_row_reader {
    global_tracking & G;
    size_t local_next_report = 256;
    size_t local_nrows_at_last_report = 0;
    size_t local_bytes_at_last_report = 0;
    size_t local_coeffs = 0;
    std::vector < typerow_t * > local_rows;
    filter_matrix_t * mat;
    int i;
    local_row_reader(global_tracking & G, filter_matrix_t * mat) : G(G), mat(mat) {
        i = omp_get_thread_num();
    }
    void read(std::istream & fi, off_t bytes_to_read, global_tracking & G) {

        std::string s;

        // countbuf cc(*fi.rdbuf());
        // std::istream fi2(&cc);

        off_t start = fi.tellg();

        for (off_t pos ; ((pos = fi.tellg()) - start)  < bytes_to_read ; ) {
            std::vector < typerow_t > primes;
            /* Insert a temporary marker. We'll use it for storing the
             * size, eventually */
#ifdef FOR_DL
            ideal_merge_t zz { 0, 0 };
            primes.insert(primes.begin(), zz);
#else
            primes.emplace(primes.begin(), 0);
#endif
            s.clear();
            {
                /* this is going to be a bit ugly, I know */
                std::getline(fi, s);
                // pos += s.size();
                if (s[0] == '#') continue;

                /* At this point, there are several ways to parse the
                 * string s and make a relation out of it.
                 * The whole point of reading in parallel is that we're
                 * going to have several threads that keep the I/O layer
                 * happy, so that what happens CPU-bound shouldn't be
                 * *that* important.
                 *
                 * However, we do care about lock contention. And
                 * apparently, there are several bad ideas that one can
                 * come up with, which eventually destroy performance.
                 */
#if 0
                /* The #1 outrageously bad idea is this one. */
                std::istringstream ss(s);
#if 0
                /* This is just catastrophic in terms of performance, for
                 * reasons that I don't understand completely.
                 */
                ss.imbue(std::locale(std::locale(), new csv_reader()));
                ss >> std::hex;
                {
                    std::string a,b;
                    ss >> a >> b;
                }
                for( ; ss.good() ; ) {
                    index_t x;
                    ss >> x;
#if 0
#ifdef FOR_DL
                    ideal_merge_t xx { x, 1};
                    primes.push_back(xx);
#else
                    primes.emplace_back(x);
#endif
#endif
                }
#endif
#else
                char * p = &s[0];
                char * z = p + s.size();
                for( ; *p && *p != ':' ; p++);
                for( ; p++ != z ; ) {
#if 0
                    char * q;
                    index_t x MAYBE_UNUSED = strtoul(p, &q, 16);
#else
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
                    // index_t x MAYBE_UNUSED = strtoul(p, &q, 16);
                    index_t x = 0;
                    for(int c ; (c = ugly[(int) (unsigned char) *p++]) != 255 ; x = x*16+c) ;
                    p--;
#endif
                    if (x < mat->skip)
                        continue;
#if 1
#ifdef FOR_DL
                    ideal_merge_t xx { x, 1};
                    primes.push_back(xx);
#else
                    primes.emplace_back(x);
#endif
#endif
                }
#endif
            }
#if 1
            /* sqrt has a weirdo check in here. Not sure we want this 
               if ((rc != 2 && feof(fi)) || pos >= spos_tab[i+1])
               break;
               */
#ifdef FOR_DL
            struct cmp {
                bool operator() (typerow_t const &a, typerow_t const &b) const {
                    return a.id < b.id;
                }
            };
            std::sort(primes.begin() + 1, primes.end(), cmp());
#else
            std::sort(primes.begin() + 1, primes.end());
#endif
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

            /* Pay attention to the special marker ! */
            unsigned int z = primes.size() - 1;
            local_coeffs += z;

#ifdef FOR_DL
            primes[0].id = z;
#else
            primes[0] = z;
#endif

            /* 0 here must eventually become the row index, but we can't
             * write it right now.
             */
            typerow_t *newrow = heap_alloc_row(0, z);
            compressRow(newrow, &primes[0], z);
            local_rows.push_back(newrow);
#endif

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
        /* right, these are ugly pointers. well, this is how the
         * structures in merge.c work, and we have to live with it.
         * Anyway all the memory is owned by merge_heap.c
         */
        off_t bytes_to_read = G.spos_tab[i + 1] - G.spos_tab[i];
        local_row_reader L(G, mat);
        L.read(fi, bytes_to_read, G);
#pragma omp barrier
#pragma omp single
        {
            fmt::fprintf(stderr, "# All local buffers are ready\n");
        }

        /* now we want to collect the number of rows for
         * everyone, renumber all relations, and store them in order in
         * mat->rows.
         */
#pragma omp barrier
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
