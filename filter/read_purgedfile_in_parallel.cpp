#include "cado.h" // IWYU pragma: keep

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "filter_config.h"
#include "macros.h"
#include "merge_heap.h"
#include "merge_replay_matrix.h"
#include "omp_proxy.h"
#include "read_purgedfile_in_parallel.hpp"
#include "sparse.h"
#include "timing.h"
#include "typedefs.h"

#include "utils_cxx.hpp"

static uint64_t * rows_per_thread;
static off_t * spos_tab;
static size_t global_bytes;
static size_t global_next_report;
static double tt0;
static size_t global_nrows;
static size_t global_nthreads;

static void global_init(int nthreads, off_t endpos)
{
    global_nthreads = nthreads;
    rows_per_thread = (uint64_t *) malloc(nthreads * sizeof(uint64_t));
    tt0 = wct_seconds();
    /* these are accumulators, and this function may conceivably be
     * called more than once in the same process */
    global_nrows = 0;
    global_bytes = 0;
    global_next_report = 1024;
    spos_tab = (off_t *) malloc((nthreads + 1) * sizeof(off_t));

    for (int i = 0; i < nthreads; i++)
        spos_tab[i] = (endpos * i) / nthreads;
    spos_tab[nthreads] = endpos;
}

static void global_clear()
{
    free(rows_per_thread);
    free(spos_tab);
}

static void global_print_report() {
    double dt = wct_seconds() - tt0;
    printf
        ("# Read %zu relations in %.1fs -- %.1f MB/s -- %.1f rels/s\n",
         global_nrows, dt, (global_bytes >> 20) / dt,
         global_nrows / dt);
    fflush(stdout);

    global_next_report *= 2;
}

static inline char * hacked_strtoul16(index_t * px, char * p)/*{{{*/
{
    /* functionally equivalent to:
         char * q;
         * px = strtoul(p, &q, 16);
         return q;
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
    index_t x = 0;
    for(int c ; (c = ugly[(int) (unsigned char) *p++]) != 255 ; x = x*16+c) ;
    p--;
    *px = x;
    return p;
}
/*}}}*/

struct vector_of_typerow_pointer_s {
    typerow_t ** x;
    size_t size;
    size_t alloc;
};

typedef struct vector_of_typerow_pointer_s vector_of_typerow_pointer[1];
typedef struct vector_of_typerow_pointer_s * vector_of_typerow_pointer_ptr;
typedef const struct vector_of_typerow_pointer_s * vector_of_typerow_pointer_srcptr;

static void vector_of_typerow_pointer_init(vector_of_typerow_pointer_ptr V)
{
    V->x = NULL;
    V->size = V->alloc = 0;
}

static void vector_of_typerow_pointer_clear(vector_of_typerow_pointer_ptr V)
{
    free(V->x);
    V->x = NULL;
    V->size = V->alloc = 0;
}

static void vector_of_typerow_pointer_push_back(vector_of_typerow_pointer_ptr V, typerow_t * p)
{
    if (V->size >= V->alloc) {
        size_t newalloc = MAX(V->size * 2, 16);
        checked_realloc(V->x, newalloc);
        V->alloc = newalloc;
    }
    V->x[V->size++] = p;
}

struct vector_of_typerow_s {
    typerow_t * x;
    size_t size;
    size_t alloc;
};

typedef struct vector_of_typerow_s vector_of_typerow[1];
typedef struct vector_of_typerow_s * vector_of_typerow_ptr;
typedef const struct vector_of_typerow_s * vector_of_typerow_srcptr;

static void vector_of_typerow_init(vector_of_typerow_ptr V)
{
    V->x = NULL;
    V->size = V->alloc = 0;
}

static void vector_of_typerow_clear(vector_of_typerow_ptr V)
{
    free(V->x);
    V->x = NULL;
    V->size = V->alloc = 0;
}

static void vector_of_typerow_empty(vector_of_typerow_ptr V)
{
    V->size = 0;
}

static void vector_of_typerow_push_back(vector_of_typerow_ptr V, const typerow_t * p)
{
    if (V->size >= V->alloc) {
        size_t newalloc = MAX(V->size * 2, 16);
        checked_realloc(V->x, newalloc);
        V->alloc = newalloc;
    }
    memcpy(V->x + V->size++, p, sizeof(typerow_t));
}

static void read_local_rows(vector_of_typerow_pointer_ptr V, FILE * fi, off_t bytes_to_read, uint64_t skip)
{
    size_t local_next_report = 256;
    size_t local_nrows_at_last_report = 0;
    size_t local_bytes_at_last_report = 0;

    /* reuse local variables ; this avoids frequent roundtrips to the
     * malloc layer */

    char line[4096];
    vector_of_typerow primes;
    vector_of_typerow_init(primes);

    long start = ftell(fi);

    for (long pos ; ((pos = ftell(fi)) - start)  < bytes_to_read ; ) {

        vector_of_typerow_empty(primes);
        line[0]='\0';

        /* Insert a temporary marker. We'll use it for storing the
         * size, eventually */
        typerow_t zz;
        setCell(&zz, 0, 0, 0);
        vector_of_typerow_push_back(primes, & zz);
        {
            /* this is going to be a bit ugly, I know */
            if (fgets(line, sizeof(line), fi) == NULL) break;

            if (line[0] == '#') continue;

            /* see "BAD IDEAS FOR PARSING LOOP" below for things that
             * I tried and didn't play out well.  */
            char * p = line;
            size_t const len = strlen(line);
            /* otherwise 4096 is not enough ! */
            ASSERT_ALWAYS(len < (sizeof(line) - 1));
            /* Where the interesting part of the line stops: at the EOL
             * delimiter, or at the end of the string if the last line of
             * the file happens not to be terminated. (A line that is too
             * long for our buffer would look the same, but the assertion
             * above has ruled that out.) */
            char * z = line + len;
            if (len > 0 && z[-1] == '\n')
                z--;

            for( ; *p && *p != ':' ; p++);
            /* note the comparison: p may well jump _over_ z, e.g. if
             * there is no colon at all in the line. We must not run past
             * the end of the buffer in that case. */
            for( ; p < z ; ) {
                p++;
                index_t x;
                p = hacked_strtoul16(&x, p);
                if (x < skip)
                    continue;
                typerow_t xx;
                setCell(&xx, 0, x, 1);
                vector_of_typerow_push_back(primes, &xx);
            }
        }

        qsort(primes->x + 1, primes->size - 1, sizeof(typerow_t), cmp_typerow_t);

        unsigned int jt = 1;
        for (unsigned int it = 1; it != primes->size;) {
            primes->x[jt] = primes->x[it];
            unsigned int kt = it;
            ++kt;
#ifdef FOR_DL
            for (; kt != primes->size && primes->x[kt].id == primes->x[jt].id; ++kt)
                primes->x[jt].e += primes->x[kt].e;
            jt++;
#else
            for (; kt != primes->size && primes->x[kt] == primes->x[jt]; ++kt);
            jt += ((kt - it) & 1);
#endif
            it = kt;
        }
        primes->size = jt;

        /* Pay attention to the special marker ! and update it, too. */
        unsigned int z = primes->size - 1;
        setCell(primes->x, 0, z, 0);

        /* 0 here must eventually become the row index, but we can't
         * write it right now. We'll do so later on.  */
        typerow_t *newrow = heap_alloc_row(0, z);
        compressRow(newrow, &primes->x[0], z);
        vector_of_typerow_pointer_push_back(V, newrow);

        /* At this point we should consider reporting. */
        size_t local_bytes = pos - start;
        if (V->size >= local_next_report) {
#pragma omp critical
            {
                global_nrows += V->size - local_nrows_at_last_report;
                global_bytes += local_bytes - local_bytes_at_last_report;
                local_nrows_at_last_report = V->size;
                local_bytes_at_last_report = local_bytes;
                if (global_nrows >= global_next_report)
                    global_print_report();
                local_next_report += global_next_report / 2 / global_nthreads;
            }
        }
    }
#pragma omp critical
    {
        global_nrows += V->size - local_nrows_at_last_report;
        size_t local_bytes = ftell(fi) - start;
        global_bytes += local_bytes - local_bytes_at_last_report;
    }
    vector_of_typerow_clear(primes);
}

uint64_t read_purgedfile_in_parallel(filter_matrix_t * mat,
				     const char *filename)
{
    off_t endpos;

    {
	FILE * f = fopen(filename, "r");
        DIE_ERRNO_DIAG(f == NULL, "fopen(%s)", filename);
	int rc = fseek(f, 0, SEEK_END);
        DIE_ERRNO_DIAG(rc < 0, "fseek(%s)", filename);
        endpos = ftell(f);
        fclose(f);
    }

    unsigned int nthreads = omp_get_max_threads();

    /* cap the number of I/O threads */
    if (nthreads > MAX_IO_THREADS)
	nthreads = MAX_IO_THREADS;

    /* All threads get their private reading head. */
#pragma omp parallel num_threads(nthreads)
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
        {
            fprintf(stderr, "# %s: Doing I/O with %d threads\n", filename,
                    omp_get_num_threads());
            global_init(omp_get_num_threads(), endpos);
        }
        /* the omp single construct above ends with an implicit barrier,
         * so that spos_tab[] is set and visible to everyone here */

        int i = omp_get_thread_num();
        FILE * fi;
        char buffer[1 << 16];

        fi = fopen(filename, "r");
        DIE_ERRNO_DIAG(fi == NULL, "fopen(%s)", filename);

        setbuffer(fi, buffer, sizeof(buffer));

	int rc = fseek(fi, spos_tab[i], SEEK_SET);
        DIE_ERRNO_DIAG(rc < 0, "fseek(%s)", filename);

        /* Except when we're at the beginning of the stream, read until
         * we get a newline. Note that we must stop at end of file too. A
         * file whose last line is not terminated would otherwise have us
         * spin forever, since fgetc() keeps returning EOF. */
        if (i > 0) {
            for(int c ; (c = fgetc(fi)) != '\n' && c != EOF ; );
        }
        spos_tab[i] = ftell(fi);

#pragma omp barrier

        /* Cut positions are only ever moved forwards, and they start out
         * sorted, so this holds. read_local_rows() would read the same
         * rows twice, or skip some, if it did not. */
        ASSERT_ALWAYS(spos_tab[i] <= spos_tab[i + 1]);

        off_t bytes_to_read = spos_tab[i + 1] - spos_tab[i];
        vector_of_typerow_pointer V;
        vector_of_typerow_pointer_init(V);
        read_local_rows(V, fi, bytes_to_read, mat->skip);
        rows_per_thread[i] = V->size;

#pragma omp barrier

        /* now we want to collect the number of rows for everyone,
         * renumber all relations, and store them in order in mat->rows.
         */
        uint64_t index = 0;
        for(int j = 0 ; j < i ; j++)
            index += rows_per_thread[j];
        for (size_t s = 0 ; s < V->size ; s++) {
            typerow_t * r = V->x[s];
            rowCell((r - 1), 0) = index;
            mat->rows[index] = r;
            index++;
        }
        vector_of_typerow_pointer_clear(V);

        fclose(fi);
    }
    global_print_report();
    global_clear();

    return global_nrows;
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

