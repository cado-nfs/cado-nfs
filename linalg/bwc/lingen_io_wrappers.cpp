#include "cado.h"
#include <sys/types.h>
#include <sys/stat.h>
#include "lingen_io_wrappers.hpp"
#include "timing.h"
#include "lingen_average_matsize.hpp"
#include "lingen_io_matpoly.hpp"
#include "fmt/format.h"

#pragma GCC diagnostic ignored "-Wunused-parameter"

/* This one is currently in lingen_io_matpoly.cpp, but that file will go
 * away eventually */

extern unsigned int io_matpoly_block_size;

static int mpi_rank()
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
}


/* XXX It seems likely that lingen_average_matsize.?pp will go away
 */
double lingen_io_wrapper_base::average_matsize() const
{
    return ::average_matsize(ab, m, n, false);
}

unsigned int lingen_io_wrapper_base::preferred_window() const
{
    return iceildiv(io_matpoly_block_size, average_matsize());
}

double lingen_file_input::average_matsize() const
{
    return ::average_matsize(ab, m, n, ascii);
}

void lingen_file_input::open_file()
{
    if (mpi_rank()) return;
    fopen(filename.c_str(), ascii ? "r" : "rb");
}

void lingen_file_input::close_file()
{
    if (mpi_rank()) return;
    fclose(f);
}

size_t lingen_file_input::guessed_length() const
{
    unsigned long guess;

    if (mpi_rank() == 0) {
        struct stat sbuf[1];
        int rc = fstat(fileno(f), sbuf);
        if (rc < 0) {
            perror(filename.c_str());
            exit(EXIT_FAILURE);
        }

        size_t filesize = sbuf->st_size;

        size_t avg = average_matsize();

        if (!ascii) {
            if (filesize % avg) {
                fprintf(stderr, "File %s has %zu bytes, while its size"
                        " should be a multiple of %zu bytes "
                        "(assuming binary input; "
                        "perhaps --ascii is missing ?).\n",
                        filename.c_str(), filesize, avg);
                exit(EXIT_FAILURE);
            }
            guess = filesize / avg;
        } else {
            double expected_length = filesize / avg;
            printf("# Expect roughly %.2f items in the sequence.\n", expected_length);

            /* First coefficient is always lighter, so we add a +1. */
            guess = 1 + expected_length;
        }
    }
    MPI_Bcast(&guess, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
    return guess;
}

ssize_t lingen_file_input::read_to_matpoly(matpoly & dst, unsigned int k0, unsigned int k1)
{
    return matpoly_read(ab, f, dst, k0, k1, ascii, 0);
}

ssize_t lingen_random_input::read_to_matpoly(matpoly & dst, unsigned int k0, unsigned int k1)
{
    dst.fill_random(k0, k1, rstate);
    return k1 - k0;
}

void lingen_E_from_A::initial_read()
{
    // see bm_io::begin_read
    // and bm_io::compute_initial_F

    /* Our goal is to determine t0 and F0, and ultimately ensure that the
     * constant coefficient of E at t=t0 is a full-rank matrix.
     *
     * The initial constant coefficient of E is in fact the
     * coefficient of degree t0 of A'*F, where column j of A' is
     * column j of A, divided by X if j >= bm.d.nrhs.
     *
     * (recall that 0<=bm.d.nrhs<=n)
     *
     * We need m independent columns.
     */

}

ssize_t lingen_E_from_A::read_to_matpoly(matpoly & dst, unsigned int k0, unsigned int k1)
{
    // see bm_io::read1
    return 0;
}

template<>
ssize_t lingen_scatter<matpoly>::write_from_matpoly(matpoly const & src, unsigned int k0, unsigned int k1)
{
    /* copy coeffs k0 to k1 from src to E */
    return 0;
}

template<>
ssize_t lingen_scatter<bigmatpoly>::write_from_matpoly(matpoly const & src, unsigned int k0, unsigned int k1)
{
    // see bigmatpoly_producer_task
    return 0;
}

template<>
ssize_t lingen_gather<matpoly>::read_to_matpoly(matpoly & dst, unsigned int k0, unsigned int k1)
{
    /* copy coeffs */
    return 0;
}

template<>
ssize_t lingen_gather<bigmatpoly>::read_to_matpoly(matpoly & dst, unsigned int k0, unsigned int k1)
{
    /* copy coeffs */
    return 0;
}

ssize_t lingen_F_from_PI::read_to_matpoly(matpoly & dst, unsigned int k0, unsigned int k1)
{
    /* see compute_final_F -- this is where some rewrite is needed for
     * the binary case */
    return 0;
}

void lingen_output_to_singlefile::open_file()
{
    if (mpi_rank()) return;
    std::ios_base::openmode mode = std::ios_base::out;
    if (!ascii) mode |= std::ios_base::binary;
    os = std::unique_ptr<std::ofstream>(new std::ofstream(filename, mode));
    done_open = true;
}

ssize_t lingen_output_to_singlefile::write_from_matpoly(matpoly const & src, unsigned int k0, unsigned int k1)
{
    /* This should be fixed. We seem to be used to writing this
     * transposed. That's slightly weird.
     */
    return matpoly_write(ab, *os, src, k0, k1, ascii, 1);
}

ssize_t lingen_output_to_splitfile::write_from_matpoly(matpoly const & src, unsigned int k0, unsigned int k1)
{
    /* This should be fixed. We seem to be used to writing this
     * transposed. That's slightly weird.
     */
    if (!done_open)
        open_file();
    return matpoly_write_split(ab, fw, src, k0, k1, ascii);
}

void lingen_output_to_splitfile::open_file()
{
    if (mpi_rank()) return;
    ASSERT_ALWAYS(!done_open);
    std::ios_base::openmode mode = std::ios_base::out;
    if (!ascii) mode |= std::ios_base::binary;
    for(unsigned int i = 0 ; i < m ; i++) {
        for(unsigned int j = 0 ; j < n ; j++) {
            std::string s = fmt::format(pattern, i, i+1, j, j+1);
            fw.emplace_back(std::ofstream { s, mode });
            DIE_ERRNO_DIAG(!fw.back(), "open", s.c_str());
        }
    }
    done_open = true;
}

lingen_output_to_splitfile::lingen_output_to_splitfile(abdst_field ab, unsigned int m, unsigned int n, std::string const & pattern, bool ascii)
    : lingen_output_wrapper_base(ab, m, n)
    , pattern(pattern)
    , ascii(ascii)
{
}

lingen_output_to_sha1sum::~lingen_output_to_sha1sum()
{
    char out[41];
    f.checksum(out);
    printf("checksum(%s): %s\n", who.c_str(), out);
}

ssize_t
lingen_output_to_sha1sum::write_from_matpoly(matpoly const& src,
                               unsigned int k0,
                               unsigned int k1)
{
    return 0;
}


void pipe(lingen_input_wrapper_base & in, lingen_output_wrapper_base & out, bool print_progress)
{
    unsigned int window = std::min(in.preferred_window(), out.preferred_window());
    if (window == UINT_MAX) {
        window=4096;
    }
    matpoly F(in.ab, in.m, in.n, window);

    double tt0 = wct_seconds();
    double next_report_t = tt0;
    size_t next_report_k = 0;
    size_t expected = in.guessed_length();
    for(size_t done = 0 ; ; ) {
        ssize_t n = in.read_to_matpoly(F, 0, window);
        if (n < 0) break;
        ssize_t nn = out.write_from_matpoly(F, 0, n);
        if (nn < n) {
            fprintf(stderr, "short write\n");
            exit(EXIT_FAILURE);
        }
        done += n;
        if (!print_progress) continue;
        if (done >= next_report_k) {
            double tt = wct_seconds();
            if (tt > next_report_t) {
                printf(
                        "Processed %zu coefficients (%.1f%%)"
                        " in %.1f s (%.1f MB/s)\n",
                        done, 100.0 * done / expected,
                        tt-tt0, done * in.average_matsize() / (tt-tt0)/1.0e6);
                next_report_t = tt + 10;
                next_report_k = done + expected / 100;
            }
        }
    }
}

#if 0
/* {{{ reading A and writing F ... */
struct bm_io {/*{{{*/
    bmstatus & bm;
    unsigned int t0 = 0;
    FILE ** fr = NULL; /* array of n files when split_input_file is supported
                   (which is not the case as of now),
                   or otherwise just a 1-element array */
    char * iobuf = NULL;
    const char * input_file = NULL;
    const char * output_file = NULL;
#ifdef SELECT_MPFQ_LAYER_u64k1
    static constexpr const unsigned int simd = ULONG_BITS;
#else
    static constexpr const unsigned int simd = 1;
#endif
    int ascii = 0;
    /* This is only a rolling window ! */
    matpoly A, F;
    unsigned int (*fdesc)[2] = NULL;
    /* This k is the coefficient in A(X) div X of the next coefficient to
     * be read. This is thus the total number of coefficients of A(X) div
     * X which have been read so far.
     * In writing mode, k is the number of coefficients of F which have
     * been written so far.
     */
    unsigned int next_coeff_to_fetch_from_source = 0;   // ROOT ONLY!
    unsigned int next_coeff_to_consume = 0;     // ROOT ONLY!

    unsigned int guessed_length = 0;

    bool leader() const {
        int rank;
        MPI_Comm_rank(bm.com[0], &rank);
        return rank == 0;
    }
    unsigned int set_write_behind_size();
    void zero1(unsigned int deg);
    unsigned int fetch_more_from_source(unsigned int io_window, unsigned int batch);
    bm_io(bm_io const&)=delete;
    bm_io(bmstatus & bm, const char * input_file, const char * output_file, int ascii);
    ~bm_io();
    void begin_read();
    void end_read();
    void guess_length();
    void compute_initial_F() ;

    template<class Consumer, class Sink>
        void compute_final_F(Sink & S, Consumer& pi);
    template<class Producer>
        void compute_E(Producer& E, unsigned int expected, unsigned int allocated);
    template<typename T, typename Sink>
        void output_flow(T & pi);
};
/*}}}*/

/* The reading mode of bm_io is streaming, but with a look-back
 * functionality: we want to be able to access coefficients a few places
 * earlier, so we keep them in memory.
 *
 * The writing mode has a write-ahead feature. Coefficients of the result
 * are written at various times, some earlier than others. The time span
 * is the same as for the reading mode.
 */


/* We write the coefficients of the reversed polynomial \hat{F*\pi}, in
 * increasing degree order. Thus the first coefficients written
 * corresponds to high degree coefficients in \pi.
 * This is mostly for historical reasons, since in fact, we would prefer
 * having the coefficients in the reverse order.
 */

/* let mindelta and maxdelta be (as their name suggests) the minimum and
 * maximum over all deltas corresponding to solution columns.
 *
 * For 0<=i<n, coeff i,j,k of pi becomes coeff i,j,k' of the final f,
 * with k'=k-(maxdelta-delta[j]).
 *
 * For 0<=i<m, coeff n+i,j,k of pi becomes coeff c[i],j,k' of the final f,
 * with k'=k-(maxdelta-delta[j])-(t0-e[j]).
 *
 * Therefore the maximum write-behind distance is (maxdelta-mindelta)+t0.
 * We need one coeff more (because offset goes from 0 to maxoffset).
 */
unsigned int bm_io::set_write_behind_size()/*{{{*/
{
    bw_dimensions & d = bm.d;
    unsigned int mindelta, maxdelta;
    std::tie(mindelta, maxdelta) = bm.get_minmax_delta_on_solutions();
    unsigned int window = maxdelta - mindelta + t0 + 1;
    if (d.nrhs) {
        /* in sm-outside-matrix mode for DLP, we form the matrix F
         * slightly differently, as some *rows* are shifted out before
         * writing */
        window++;
    }
    window = simd * iceildiv(window, simd);
    if (leader()) {
        // F.realloc(window);
        F.zero_pad(window);
        // F.set_size(window);
    }
    return window;
}/*}}}*/

/* {{{ bm_output_* classes: the link from the in-memory structures to the
 * filesystem.
 * Note that coefficients are written transposed */
struct bm_output_singlefile {/*{{{*/
    bm_io & aa;
    matpoly const & P;
    std::ofstream f;
    char * iobuf;
    char * filename;
    bm_output_singlefile(bm_io &aa, matpoly const & P, const char * suffix = "")
        : aa(aa), P(P)
    {
        if (!aa.leader()) return;
        int rc = asprintf(&filename, "%s%s", aa.output_file, suffix);
        ASSERT_ALWAYS(rc >= 0);
        std::ios_base::openmode mode = std::ios_base::out;
        if (!aa.ascii) mode |= std::ios_base::binary;  
        f.open(filename, mode);
        DIE_ERRNO_DIAG(!f, "fopen", filename);
        iobuf = (char*) malloc(2 * io_matpoly_block_size);
        f.rdbuf()->pubsetbuf(iobuf, 2 * io_matpoly_block_size);
    }
    ~bm_output_singlefile()
    {
        if (!aa.leader()) return;
        printf("Saved %s\n", filename);
        free(filename);
        f.close();
        free(iobuf);
    }
    void write1(unsigned int deg)
    {
        if (!aa.leader()) return;
        deg = deg % P.capacity();
        matpoly_write(aa.bm.d.ab, f, P, deg, deg + 1, aa.ascii, 1);
    }
};/*}}}*/

struct bm_output_splitfile {/*{{{*/
    bm_io & aa;
    matpoly const & P;
    std::vector<std::ofstream> fw;
    bm_output_splitfile(bm_io & aa, matpoly const & P, const char * suffix = "")
        : aa(aa), P(P)
    {
        if (!aa.leader()) return;
        std::ios_base::openmode mode = std::ios_base::out;
        if (!aa.ascii) mode |= std::ios_base::binary;  
        for(unsigned int i = 0 ; i < P.m ; i++) {
            for(unsigned int j = 0 ; j < P.n ; j++) {
                char * str;
                int rc = asprintf(&str, "%s.sols%d-%d.%d-%d%s",
                        aa.output_file,
                        j, j + 1,
                        i, i + 1, suffix);
                ASSERT_ALWAYS(rc >= 0);
                fw.emplace_back(std::ofstream { str, mode } );
                DIE_ERRNO_DIAG(!fw.back(), "fopen", str);
                free(str);
            }
        }
        /* Do we want specific caching bufs here ? I doubt it */
        /*
           iobuf = (char*) malloc(2 * io_matpoly_block_size);
           setbuffer(f, iobuf, 2 * io_matpoly_block_size);
           */
    }
    ~bm_output_splitfile() {
    }
    void write1(unsigned int deg)
    {
        if (!aa.leader()) return;
        deg = deg % P.capacity();
        matpoly_write_split(aa.bm.d.ab, fw, P, deg, deg + 1, aa.ascii);
    }

};/*}}}*/
struct bm_output_checksum {/*{{{*/
    bm_io & aa;
    matpoly const & P;
    sha1_checksumming_stream f;
    const char * suffix;
    bm_output_checksum(bm_io & aa, matpoly const & P, const char * suffix = NULL)
        : aa(aa), P(P), suffix(suffix) { }
    ~bm_output_checksum() {
        char out[41];
        f.checksum(out);
        if (suffix)
            printf("checksum(%s): %s\n", suffix, out);
        else
            printf("checksum: %s\n", out);
    }
    void write1(unsigned int deg)
    {
        if (!aa.leader()) return;
        deg = deg % P.capacity();
        matpoly_write(aa.bm.d.ab, f, P, deg, deg + 1, aa.ascii, 1);
    }
};/*}}}*/
/* }}} */

void bm_io::zero1(unsigned int deg)/*{{{*/
{
    bw_dimensions & d = bm.d;
    unsigned int n = d.n;
    deg = deg % F.capacity();
    for(unsigned int j = 0 ; j < n ; j++) {
        for(unsigned int i = 0 ; i < n ; i++) {
#ifndef SELECT_MPFQ_LAYER_u64k1
            abset_zero(d.ab, F.coeff(i, j, deg));
#else
            F.coeff_accessor(i, j, deg) += F.coeff(i, j, deg);
#endif
        }
    }
}/*}}}*/

/* Transfers of matrix entries, from memory to memory ; this is possibly
 * almost a transparent layer, and the matpoly_* instances really _are_
 * transparent layers. However we want a unified interface for the case
 * where we read an MPI-scattered matrix by chunks.
 *
 * Given that we are talking memory-to-memory, the classes below are used
 * as follows. The *_write_task classes are at the beginning of the
 * computation, where the matrix A is read from disk, and the matrix E is
 * written to (possibly distributed) memory. The one which matters is E
 * here. Symmetrically, the *_read_task are for reading the computed
 * matrix pi from (possibly distributed) memory, and writing the final
 * data F to disk.
 *
 * Pay attention to the fact that the calls to these structures are
 * collective calls.
 */
#ifdef ENABLE_MPI_LINGEN
class bigmatpoly_consumer_task { /* {{{ */
    /* This reads a bigmatpoly, by chunks, so that the memory footprint
     * remains reasonable. */
    size_t simd;
    bigmatpoly const & xpi;
    matpoly pi; /* This is only temp storage ! */
    unsigned int B;
    unsigned int k0;
    int rank;

    public:
    bigmatpoly_consumer_task(bm_io & aa, bigmatpoly const & xpi) : simd(aa.simd), xpi(xpi) {
        bmstatus & bm = aa.bm;
        bw_dimensions & d = bm.d;
        unsigned int m = d.m;
        unsigned int n = d.n;
        unsigned int b = m + n;
        MPI_Comm_rank(bm.com[0], &rank);

        /* Decide on the temp storage size */
        double avg = average_matsize(d.ab, n, n, aa.ascii);
        B = iceildiv(io_matpoly_block_size, avg);
        B = simd * iceildiv(B, simd);
        if (!rank && !random_input_length) {
            printf("Writing F to %s\n", aa.output_file);
            printf("Writing F by blocks of %u coefficients"
                    " (%.1f bytes each)\n", B, avg);
        }
        pi = matpoly(d.ab, b, b, B);

        k0 = UINT_MAX;
    }

    inline unsigned int chunk_size() const { return B; }

    inline unsigned int size() { return xpi.get_size(); }

    /* locked means that we don't want to change the current access window */
    inline absrc_elt coeff_const_locked(unsigned int i, unsigned int j, unsigned int k) {
        ASSERT(!rank);
        ASSERT(k0 != UINT_MAX && k - k0 < B);
        /* Note that only the const variant is ok to use in the binary case */
        return pi.coeff(i, j, k - k0);
    }

    /* Here, we adjust the access window so as to access coefficient k */
    absrc_elt coeff_const(unsigned int i, unsigned int j, unsigned int k) {
        if (k0 == UINT_MAX || k - k0 >= B) {
            k0 = k - (k % B);
            pi.zero();
            pi.set_size(B);
            xpi.gather_mat_partial(pi, k0, MIN(xpi.get_size() - k0, B));
        }
        if (rank) return NULL;
        return coeff_const_locked(i, j, k);
    }
};      /* }}} */
class bigmatpoly_producer_task { /* {{{ */
    /* This writes a bigmatpoly, by chunks, so that the memory footprint
     * remains reasonable.
     * Note that in any case, the coefficient indices must be progressive
     * in the write.
     */
    size_t simd;
    bigmatpoly & xE;
    matpoly E; /* This is only temp storage ! */
    unsigned int B;
    unsigned int k0;
    int rank;

    /* forbid copies */
    bigmatpoly_producer_task(bigmatpoly_producer_task const &) = delete;

    public:
    bigmatpoly_producer_task(bm_io & aa, bigmatpoly & xE) : simd(aa.simd), xE(xE) {
        bmstatus & bm = aa.bm;
        bw_dimensions & d = bm.d;
        unsigned int m = d.m;
        unsigned int n = d.n;
        unsigned int b = m + n;
        MPI_Comm_rank(aa.bm.com[0], &rank);


        /* Decide on the MPI chunk size */
        double avg = average_matsize(d.ab, m, n, 0);
        B = iceildiv(io_matpoly_block_size, avg);
        B = simd * iceildiv(B, simd);

        if (!rank) {
            /* TODO: move out of here */
            if (aa.input_file)
                printf("Reading A from %s\n", aa.input_file);
            printf("Computing E by chunks of %u coefficients (%.1f bytes each)\n",
                    B, avg);
        }

        E = matpoly(d.ab, m, b, B);

        /* Setting E.size is rather artificial, since we use E essentially
         * as an area where we may write coefficients freely. The only aim is
         * to escape some safety checks involving ->size in matpoly_part */
        E.zero();

        k0 = UINT_MAX;
    }

    inline unsigned int chunk_size() const { return B; }

    // inline unsigned int size() { return xE.size; }
    inline void set_size(unsigned int s) {
        xE.set_size(s);
    }

    inline abdst_elt coeff_locked(unsigned int i, unsigned int j, unsigned int k) {
        ASSERT(k0 != UINT_MAX && k - k0 < B);
        ASSERT(!rank);

        E.set_size(E.get_size() + (E.get_size() == k - k0));
        ASSERT(E.get_size() == k - k0 + 1);

        return E.coeff(i, j, k - k0);
    }

    abdst_elt coeff(unsigned int i, unsigned int j, unsigned int k) {
        if (k0 == UINT_MAX) {
            E.zero();
            E.set_size(0);
            ASSERT(k == 0);
            k0 = 0;
        } else {
            if (k >= (k0 + B)) {
                /* We require progressive reads */
                ASSERT(k == k0 + B);
                xE.scatter_mat_partial(E, k0, B);
                E.zero();
                k0 += B;
            }
        }
        ASSERT(k0 != UINT_MAX && k - k0 < B);
        if (rank) return NULL;
        return coeff_locked(i, j, k);
    }

    void finalize(unsigned int length) {
        if (length > k0) {
            /* Probably the last chunk hasn't been dispatched yet */
            ASSERT(length < k0 + B);
            xE.scatter_mat_partial(E, k0, length - k0);
        }
        set_size(length);
    }
    friend void matpoly_extract_column(
        bigmatpoly_producer_task& dst, unsigned int jdst, unsigned int kdst,
        matpoly & src, unsigned int jsrc, unsigned int ksrc);
};

void matpoly_extract_column(
        bigmatpoly_producer_task& dst, unsigned int jdst, unsigned int kdst,
        matpoly & src, unsigned int jsrc, unsigned int ksrc)
{
    dst.coeff_locked(0, jdst, kdst);
    dst.E.extract_column(jdst, kdst - dst.k0, src, jsrc, ksrc);
}

/* }}} */
#endif  /* ENABLE_MPI_LINGEN */
class matpoly_consumer_task {/*{{{*/
    /* This does the same, but for a simple matpoly. Of course this is
     * much simpler ! */
    matpoly const & pi;

    public:
    matpoly_consumer_task(bm_io & aa, matpoly const & pi) :
            pi(pi)
    {
        if (!random_input_length) {
            printf("Writing F to %s\n", aa.output_file);
        }
    }

    inline unsigned int chunk_size() const { return 1; }

    inline unsigned int size() { return pi.get_size(); }

    absrc_elt coeff_const_locked(unsigned int i, unsigned int j, unsigned int k) {
        return pi.coeff(i, j, k);
    }
    inline absrc_elt coeff_const(unsigned int i, unsigned int j, unsigned int k) {
        return coeff_const_locked(i, j, k);
    }

    ~matpoly_consumer_task() { }
};/*}}}*/
class matpoly_producer_task { /* {{{ */
    size_t simd;
    matpoly & E;

    /* forbid copies */
    matpoly_producer_task(matpoly_producer_task const &) = delete;

    public:
    matpoly_producer_task(bm_io & aa, matpoly & E) : simd(aa.simd), E(E)
    {
        /* TODO: move out of here */
        if (aa.input_file)
            printf("Reading A from %s\n", aa.input_file);
    }

    inline unsigned int chunk_size() const { return simd; }
    inline void set_size(unsigned int s) {
        E.set_size(s);
    }
    // inline unsigned int size() { return E.size; }

#if 0
    abdst_elt coeff_locked(unsigned int i, unsigned int j, unsigned int k) {
        E.size += (E.size == k);
        ASSERT(E.size == k + 1);
        return E.coeff(i, j, k);
    }
    inline abdst_elt coeff(unsigned int i, unsigned int j, unsigned int k) {
        return coeff_locked(i, j, k);
    }
#endif
    void mark_coeff_as_read(unsigned int, unsigned int, unsigned int k) {
        E.set_size(E.get_size() + (E.get_size() == k));
        ASSERT(E.get_size() == k + 1);
    }
    inline absrc_elt coeff(unsigned int i, unsigned int j, unsigned int k) {
        mark_coeff_as_read(i, j, k);
        // note that for the binary case, this is a const accessor.
        return E.coeff(i, j, k);
    }
    void finalize(unsigned int length) { set_size(length); }
    ~matpoly_producer_task() { }
    friend void matpoly_extract_column(
        matpoly_producer_task& dst, unsigned int jdst, unsigned int kdst,
        matpoly & src, unsigned int jsrc, unsigned int ksrc);
};
void matpoly_extract_column(
        matpoly_producer_task& dst, unsigned int jdst, unsigned int kdst,
        matpoly & src, unsigned int jsrc, unsigned int ksrc)
{
    dst.E.extract_column(jdst, kdst, src, jsrc, ksrc);
}

/* }}} */

template<class Consumer, class Sink>
void bm_io::compute_final_F(Sink & S, Consumer& pi)/*{{{*/
{
    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    abdst_field ab = d.ab;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    int leader = rank == 0;


    /* We are not interested by pi.size, but really by the number of
     * coefficients for the columns which give solutions. */
    unsigned int maxdelta = bm.get_max_delta_on_solutions();

    if (leader) printf("Final f(X)=f0(X)pi(X) has degree %u\n", maxdelta);

    unsigned int window = F.capacity();

    /* Which columns of F*pi will make the final generator ? */
    std::vector<unsigned int> sols(n, 0);
    for(unsigned int j = 0, jj=0 ; j < m + n ; j++) {
        if (bm.lucky[j] <= 0)
            continue;
        sols[jj++]=j;
    }

    double tt0 = wct_seconds();
    double next_report_t = tt0 + 10;
    unsigned next_report_s = pi.size() / 100;

    /*
     * first compute the rhscontribs. Use that to decide on a renumbering
     * of the columns, because we'd like to get the "primary" solutions
     * first, in a sense. Those are the ones with fewer zeroes in the
     * rhscontrib part. So we would want that matrix to have its columns
     * sorted in decreasing weight order
     *
     * An alternative, possibly easier, is to have a function which
     * decides the solution ordering precisely based on the inspection of
     * this rhscoeffs matrix (?). But how should we spell that info when
     * we give it to mksol ??
     */

    /* This **modifies** the "sols" array */
    if (d.nrhs) {
        /* The data is only gathered at rank 0. So the stuff we compute
         * can only be meaningful there. On the other hand, it is
         * necessary that we call *collectively* the pi.coeff() routine,
         * because we need synchronisation of all ranks.
         */

        matpoly rhs(ab, d.nrhs, n, 1);
        rhs.zero_pad(1);

        {
            Sink Srhs(*this, rhs, ".rhs");

            /* adding the contributions to the rhs matrix. */
            for(unsigned int ipi = 0 ; ipi < m + n ; ipi++) {
                for(unsigned int jF = 0 ; jF < n ; jF++) {
                    unsigned int jpi = sols[jF];
                    unsigned int iF, offset;
                    if (ipi < n) {
                        iF = ipi;
                        offset = 0;
                    } else {
                        iF = fdesc[ipi-n][1];
                        offset = t0 - fdesc[ipi-n][0];
                    }
                    if (iF >= d.nrhs) continue;
                    ASSERT_ALWAYS(bm.delta[jpi] >= offset);
                    unsigned kpi = bm.delta[jpi] - offset;

                    ASSERT_ALWAYS(d.nrhs);
                    ASSERT_ALWAYS(iF < d.nrhs);
                    ASSERT_ALWAYS(jF < n);
                    absrc_elt src = pi.coeff_const(ipi, jpi, kpi);

                    if (leader) {
                        rhs.coeff_accessor(iF, jF, 0) += src;
                    }
                }
            }

            if (leader) {
                /* Now comes the time to prioritize the different solutions. Our
                 * goal is to get the unessential solutions last ! */
                std::vector<std::array<int, 2>> sol_score(n, {{0, 0}});
                /* score per solution is the number of non-zero coefficients,
                 * that's it. Since we have access to lexcmp2, we want to use it.
                 * Therefore, desiring the highest scoring solutions first, we
                 * negate the hamming weight.
                 */
                for(unsigned int jF = 0 ; jF < n ; jF++) {
                    sol_score[jF][1] = jF;
                    for(unsigned int iF = 0 ; iF < d.nrhs ; iF++) {
                        int z = !abis_zero(ab, rhs.coeff(iF, jF, 0));
                        sol_score[jF][0] -= z;
                    }
                }
                std::sort(sol_score.begin(), sol_score.end());

                if (leader) {
                    printf("Reordered solutions:\n");
                    for(unsigned int i = 0 ; i < n ; i++) {
                        printf(" %d (col %d in pi, weight %d on rhs vectors)\n", sol_score[i][1], sols[sol_score[i][1]], -sol_score[i][0]);
                    }
                }

                /* We'll now modify the sols[] array, so that we get a reordered
                 * F, too (and mksol/gather don't have to care about our little
                 * tricks */
                {
                    matpoly rhs2(ab, d.nrhs, n, 1);
                    rhs2.zero_pad(1);
                    for(unsigned int i = 0 ; i < n ; i++) {
                        rhs2.extract_column(i, 0, rhs, sol_score[i][1], 0);
                    }
                    rhs = std::move(rhs2);
                    if (sol_score[0][0] == 0) {
                        if (allow_zero_on_rhs) {
                            printf("Note: all solutions have zero contribution on the RHS vectors -- we will just output right kernel vectors (ok because of allow_zero_on_rhs=1)\n");
                        } else {
                            fprintf(stderr, "ERROR: all solutions have zero contribution on the RHS vectors -- we will just output right kernel vectors (maybe use allow_zero_on_rhs ?)\n");
                            rank0_exit_code = EXIT_FAILURE;
                        }
                    }
                    /* ugly: use sol_score[i][0] now to provide the future
                     * "sols" array. We'll get rid of sol_score right afterwards
                     * anyway.
                     */
                    for(unsigned int i = 0 ; i < n ; i++) {
                        sol_score[i][0] = sols[sol_score[i][1]];
                    }
                    for(unsigned int i = 0 ; i < n ; i++) {
                        sols[i] = sol_score[i][0];
                    }
                }
                Srhs.write1(0);
            }
        }
    }

    /* we need to read pi backwards. The number of coefficients in pi is
     * pilen = maxdelta + 1 - t0. Hence the first interesting index is
     * maxdelta - t0. However, for notational ease, we'll access
     * coefficients from index pi.size downwards. The latter is always
     * large enough.
     */

    ASSERT_ALWAYS(pi.size() >= maxdelta + 1 - t0);

    for(unsigned int s = 0 ; s < pi.size() ; s++) {
        unsigned int kpi = pi.size() - 1 - s;

        /* as above, we can't just have ranks > 0 do nothing. We need
         * synchronization of the calls to pi.coeff()
         */

        /* This call is here only to trigger the gather call. This one
         * must therefore be a collective call. Afterwards we'll use a
         * _locked call which is okay to call only at rank 0.
         */
        pi.coeff_const(0, 0, kpi);

        if (rank) continue;

        /* Coefficient kpi + window of F has been totally computed,
         * because of previous runs of this loop (which reads the
         * coefficients of pi).
         */
        if (kpi + window == maxdelta) {
            /* Avoid writing zero coefficients. This can occur !
             * Example:
             * tt=(2*4*1200) mod 1009, a = (tt cat tt cat * tt[0..10])
             */
            for(unsigned int j = 0 ; j < n ; j++) {
                int z = 1;
                for(unsigned int i = 0 ; z && i < n ; i++) {
                    absrc_elt src = F.coeff(i, j, 0);
                    z = abis_zero(ab, src);
                }

                if (z) {
                    /* This is a bit ugly. Given that we're going to
                     * shift one column of F, we'll have a
                     * potentially deeper write-back buffer. Columns
                     * which seemed to be ready still are, but they
                     * will now be said so only at the next step.
                     */
                    printf("Reduced solution column #%u from"
                            " delta=%u to delta=%u\n",
                            sols[j], bm.delta[sols[j]], bm.delta[sols[j]]-1);
                    window++;
                    F.realloc(window);
                    F.set_size(window);
                    bm.delta[sols[j]]--;
                    /* shift this column */
                    for(unsigned int k = 1 ; k < window ; k++) {
                        F.extract_column(j, k-1, F, j, k);
                    }
                    F.zero_column(j, window - 1);
                    break;
                }
            }
        }
        /* coefficient of degree maxdelta-kpi-window is now complete */
        if (kpi + window <= maxdelta) {
            S.write1((maxdelta-kpi) - window);
            zero1((maxdelta-kpi) - window);
        }
        /* Now use our coefficient. This might tinker with
         * coefficients up to degree kpi-(window-1) in the file F */

        if (kpi > maxdelta + t0 ) {
            /* this implies that we always have kF > delta[jpi],
             * whence we expect a zero contribution */
            for(unsigned int i = 0 ; i < m + n ; i++) {
                for(unsigned int jF = 0 ; jF < n ; jF++) {
                    unsigned int jpi = sols[jF];
                    absrc_elt src = pi.coeff_const_locked(i, jpi, kpi);
                    ASSERT_ALWAYS(abis_zero(ab, src));
                }
            }
            continue;
        }

        for(unsigned int ipi = 0 ; ipi < m + n ; ipi++) {
            for(unsigned int jF = 0 ; jF < n ; jF++) {
                unsigned int jpi = sols[jF];
                unsigned int iF, offset;
                if (ipi < n) {
                    /* Left part of the initial F is x^0 times
                     * identity. Therefore, the first n rows of pi
                     * get multiplied by this identity matrix, this
                     * is pretty simple.
                     */
                    iF = ipi;
                    offset = 0;
                } else {
                    /* next m rows of the initial F are of the form
                     * x^(some value) times some canonical basis
                     * vector. Therefore, the corresponding row in pi
                     * ends up contributing to some precise row in F,
                     * and with an offset which is dictated by the
                     * exponent of x.
                     */
                    iF = fdesc[ipi-n][1];
                    offset = t0 - fdesc[ipi-n][0];
                }
                unsigned int subtract = maxdelta - bm.delta[jpi] + offset;
                ASSERT(subtract < window);
                if (maxdelta < kpi + subtract) continue;
                unsigned int kF = (maxdelta - kpi) - subtract;
                unsigned int kF1 = kF - (iF < d.nrhs);
                if (kF1 == UINT_MAX) {
                    /* this has been addressed in the first pass,
                     * earlier.
                     */
                    continue;
                }
                absrc_elt src = pi.coeff_const_locked(ipi, jpi, kpi);
                ASSERT_ALWAYS(kF <= bm.delta[jpi] || abis_zero(ab, src));
                F.coeff_accessor(iF, jF, kF1 % window) += src;
            }
        }

        if (leader && s > next_report_s) {
            double tt = wct_seconds();
            if (tt > next_report_t) {
                printf( "Written %u coefficients (%.1f%%) in %.1f s\n",
                        s, 100.0 * s / pi.size(), tt-tt0);
                next_report_t = tt + 10;
                next_report_s = s + pi.size() / 100;
            }
        }
    }
    /* flush the pipe */
    if (leader && window <= maxdelta) {
        for(unsigned int s = window ; s-- > 0 ; )
            S.write1(maxdelta - s);
    }
}/*}}}*/

/* read 1 (or (batch)) coefficient(s) into the sliding window of input
 * coefficients of the input series A. The io_window parameter controls
 * the size of the sliding window. There are in fact two behaviours:
 *  - io_window == 0: there is no sliding window, really, and the new
 *    coefficient is appended as the last coefficient of A.
 *  - io_window > 0: there, we really have a sliding window. Coeffs
 *    occupy places in a circular fashion within the buffer.
 */
unsigned int bm_io::fetch_more_from_source(unsigned int io_window, unsigned int batch)/*{{{*/
{
    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    static unsigned int generated_random_coefficients = 0;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    ASSERT_ALWAYS(rank == 0);

    unsigned int pos = next_coeff_to_fetch_from_source;
    ASSERT_ALWAYS(A.get_size() % simd == 0);
    ASSERT_ALWAYS(next_coeff_to_fetch_from_source % simd == 0);
    ASSERT_ALWAYS(batch % simd == 0);
    ASSERT_ALWAYS(n % simd == 0);
    if (io_window) {
        pos = pos % io_window;
        ASSERT_ALWAYS(A.get_size() == io_window);
        ASSERT_ALWAYS(pos / io_window == (pos + batch - 1) / io_window);
    } else {
        ASSERT_ALWAYS(A.get_size() == next_coeff_to_fetch_from_source);
        if (next_coeff_to_fetch_from_source >= A.capacity())
            A.realloc(A.capacity() + batch);
        ASSERT_ALWAYS(next_coeff_to_fetch_from_source < A.capacity());
        A.set_size(A.get_size() + batch);
    }
    if (random_input_length) {
        if (generated_random_coefficients >= random_input_length)
            return 0;

        for (unsigned int i = 0; i < m ; i++) {
            for (unsigned int j = 0; j < n ; j++) {
                for (unsigned int b = 0; b < batch ; b += simd) {
#ifdef SELECT_MPFQ_LAYER_u64k1
                    unsigned int sq = (pos + b) / simd;
                    A.coeff(i, j)[sq] = gmp_urandomb_ui(rstate, simd);
#else
                    abrandom(d.ab, A.coeff(i, j, pos + b), rstate);
#endif
                }
            }
        }
        generated_random_coefficients += batch;
        next_coeff_to_fetch_from_source += batch;
        return batch;
    }

    for (unsigned int i = 0; i < m ; i++) {
        for (unsigned int j = 0; j < n ; j++) {
#ifdef SELECT_MPFQ_LAYER_u64k1
            memset(A.part(i, j) + (pos / simd), 0, abvec_elt_stride(d.ab, batch));
#else
            memset(A.coeff(i, j, pos), 0, abvec_elt_stride(d.ab, batch));
#endif
        }
    }

    for (unsigned int b = 0; b < batch ; b++ ) {
        for (unsigned int i = 0; i < m ; i++) {
            for (unsigned int j = 0; j < n ; j+= simd) {
                int rc;
#ifdef SELECT_MPFQ_LAYER_u64k1
                if (ascii)
                    abort();
                unsigned long data = 0;
                rc = fread(&data, sizeof(unsigned long), 1, fr[0]);
                rc = rc == 1;
                unsigned int sq = (pos + b) / simd;
                unsigned int sr = b % simd;
                for(unsigned int jr = 0 ; jr < simd ; jr++) {
                    unsigned long bit = (data & (1UL << jr)) != 0;
                    A.coeff(i, j + jr)[sq] ^= bit << sr;
                }
#else
                abdst_elt x = A.coeff(i, j, pos + b);
                if (ascii) {
                    rc = abfscan(d.ab, fr[0], x);
                    /* rc is the number of bytes read -- non-zero on success */
                } else {
                    size_t elemsize = abvec_elt_stride(d.ab, 1);
                    rc = fread(x, elemsize, 1, fr[0]);
                    rc = rc == 1;
                    abnormalize(d.ab, x);
                }
#endif
                if (!rc) {
                    if (i == 0 && j == 0) {
                        next_coeff_to_fetch_from_source += b;
                        return b;
                    }
                    fprintf(stderr,
                            "Parse error while reading coefficient (%d,%d,%d)%s\n",
                            i, j, 1 + next_coeff_to_fetch_from_source,
                            ascii ? "" : " (forgot --ascii?)");
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    next_coeff_to_fetch_from_source += batch;
    return batch;
}/*}}}*/

bm_io::bm_io(bmstatus & bm, const char * input_file, const char * output_file, int ascii)/*{{{*/
    : bm(bm)
    , input_file(input_file)
    , output_file(output_file)
    , ascii(ascii)
    , A(bm.d.ab, bm.d.m, bm.d.n, 1)
{
}/*}}}*/

bm_io::~bm_io()/*{{{*/
{
    if (fdesc) free(fdesc);
}/*}}}*/

void bm_io::begin_read()/*{{{*/
{
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    if (rank) return;

    if (random_input_length) {
        /* see below. I think it would be a bug to not do that */
        fetch_more_from_source(0, simd);
        next_coeff_to_consume++;
        return;
    }

    if (split_input_file) {
        fprintf(stderr, "--split-input-file not supported yet\n");
        exit(EXIT_FAILURE);
    }
    fr = (FILE**) malloc(sizeof(FILE*));
    fr[0] = fopen(input_file, ascii ? "r" : "rb");

    DIE_ERRNO_DIAG(fr[0] == NULL, "fopen", input_file);
    iobuf = (char*) malloc(2 * io_matpoly_block_size);
    setbuffer(fr[0], iobuf, 2 * io_matpoly_block_size);

    /* read the first coefficient ahead of time. This is because in most
     * cases, we'll discard it. Only in the DL case, we will consider the
     * first coefficient as being part of the series. Which means that
     * the coefficient reads in the I/O loop will sometimes correspond to
     * the coefficient needed at that point in time, while we will also
     * (in the DL case) need data from the previous read.
     */
    if (fetch_more_from_source(0, simd) < simd) {
        fprintf(stderr, "Read error from %s\n", input_file);
        exit(EXIT_FAILURE);
    }
    next_coeff_to_consume++;
}/*}}}*/

void bm_io::end_read()/*{{{*/
{
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);
    if (rank) return;
    if (random_input_length) return;
    fclose(fr[0]);
    free(fr);
    fr = NULL;
    free(iobuf);
    iobuf = 0;
}/*}}}*/

void bm_io::guess_length()/*{{{*/
{
    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    abdst_field ab = d.ab;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);

    if (random_input_length) {
        guessed_length = random_input_length;
        return;
    }

    if (!rank) {
        struct stat sbuf[1];
        int rc = fstat(fileno(fr[0]), sbuf);
        if (rc < 0) {
            perror(input_file);
            exit(EXIT_FAILURE);
        }

        size_t filesize = sbuf->st_size;

        if (!ascii) {
            size_t avg = average_matsize(ab, m, n, ascii);
            if (filesize % avg) {
                fprintf(stderr, "File %s has %zu bytes, while its size should be amultiple of %zu bytes (assuming binary input; perhaps --ascii is missing ?).\n", input_file, filesize, avg);
                exit(EXIT_FAILURE);
            }
            guessed_length = filesize / avg;
        } else {
            double avg = average_matsize(ab, m, n, ascii);
            double expected_length = filesize / avg;
            if (!rank)
                printf("# Expect roughly %.2f items in the sequence.\n", expected_length);

            /* First coefficient is always lighter, so we add a +1. */
            guessed_length = 1 + expected_length;
        }
    }
    MPI_Bcast(&(guessed_length), 1, MPI_UNSIGNED, 0, bm.com[0]);
}/*}}}*/

void bm_io::compute_initial_F() /*{{{ */
{
    bw_dimensions & d = bm.d;
    abdst_field ab = d.ab;
    unsigned int m = d.m;
    unsigned int n = d.n;
    int rank;
    fdesc = (unsigned int(*)[2])malloc(2 * m * sizeof(unsigned int));
    MPI_Comm_rank(bm.com[0], &rank);
    if (!rank) {
        /* read the first few coefficients. Expand A accordingly as we are
         * doing the read */

        ASSERT(A.m == m);
        ASSERT(A.n == n);

        abelt tmp MAYBE_UNUSED;
        abinit(ab, &tmp);

        /* First try to create the initial F matrix */
        printf("Computing t0\n");

        /* We want to create a full rank m*m matrix M, by extracting columns
         * from the first coefficients of A */

        matpoly M(ab, m, m, 1);
        M.zero_pad(1);

        /* For each integer i between 0 and m-1, we have a column, picked
         * from column cnum[i] of coeff exponent[i] of A which, once reduced modulo
         * the other ones, has coefficient at row pivots[i] unequal to zero.
         */
        std::vector<unsigned int> pivots(m, 0);
        std::vector<unsigned int> exponents(m, 0);
        std::vector<unsigned int> cnum(m, 0);
        unsigned int r = 0;

        for (unsigned int k = 1; r < m ; k++) {
            /*
             * k is the candidate for becoming the value t0.
             *
             * The initial constant coefficient of E is in fact the
             * coefficient of degree t0 of A'*F, where column j of A' is
             * column j of A, divided by X if j >= bm.d.nrhs.
             *
             * (recall that 0<=bm.d.nrhs<=n)
             *
             * Therefore, for j >= bm.d.nrhs, the contribution to
             * coefficient of degree k (tentative t0) of A'*F comes from
             * the following data from A:
             *  for cols of A with j < bm.d.rhs:  coeffs 0 <= deg <= k-1 
             *  for cols of A with j >= bm.d.rhs: coeffs 1 <= deg <= k
             *
             * This means that here, we're going to read data from the
             * following coefficient of A
             *  k   if bm.d.nrhs < n
             *  k-1 if bm.d.nrhs = n
             */
            unsigned int k_access = k - (bm.d.nrhs == n);

            if (rank == 0) {
                if (k_access >= next_coeff_to_consume)
                    next_coeff_to_consume++;
                ASSERT_ALWAYS(k_access < next_coeff_to_consume);
                if (k_access >= next_coeff_to_fetch_from_source) {
                    /* read a new coefficient into A */
                    fetch_more_from_source(0, simd);
                }
                ASSERT_ALWAYS(k_access <= next_coeff_to_fetch_from_source);
            }

            for (unsigned int j = 0; r < m && j < n; j++) {
                /* Extract a full column into M (column j, degree k in A) */
                /* adjust the coefficient degree to take into account the
                 * fact that for SM columns, we might in fact be
                 * interested by the _previous_ coefficient */
                M.extract_column(r, 0, A, j, k - (j < bm.d.nrhs));

                /* Now reduce it modulo all other columns */
                for (unsigned int v = 0; v < r; v++) {
                    unsigned int u = pivots[v];
                    /* the v-th column in the M is known to
                     * kill coefficient u (more exactly, to have a -1 as u-th
                     * coefficient, and zeroes for the other coefficients
                     * referenced in the pivots[0] to pivots[v-1] indices).
                     */
                    /* add M[u,r]*column v of M to column r of M */
                    for(unsigned int i = 0 ; i < m ; i++) {
#ifndef SELECT_MPFQ_LAYER_u64k1
                        abmul(ab, tmp, M.coeff(i, v, 0), M.coeff(u, r, 0));
                        abadd(ab, M.coeff(i, r, 0), M.coeff(i, r, 0), tmp);
#else
                        if (i == u) continue;
                        abelt x = { M.coeff(i, v, 0)[0] & M.coeff(u, r, 0)[0] };
                        M.coeff_accessor(i, r, 0) += x;
#endif
                    }
#ifndef SELECT_MPFQ_LAYER_u64k1
                    abset_zero(ab, M.coeff(u, r, 0));
#endif
                }
                unsigned int u = 0;
                for( ; u < m ; u++) {
                    if (!abis_zero(ab, M.coeff(u, r, 0)))
                        break;
                }
                if (u == m) {
                    printf("[X^%d] A, col %d does not increase rank (still %d)\n",
                           k - (j < bm.d.nrhs), j, r);

                    /* we need at least m columns to get as starting matrix
                     * with full rank. Given that we have n columns per
                     * coefficient, this means at least m/n matrices.
                     */

                    if (k * n > m + 40) {
                        printf("The choice of starting vectors was bad. "
                               "Cannot find %u independent cols within A\n", m);
                        exit(EXIT_FAILURE);
                    }
                    continue;
                }

                /* Bingo, it's a new independent col. */
                pivots[r] = u;
                cnum[r] = j;
                exponents[r] = k - 1;

                /* Multiply the column so that the pivot becomes -1 */
#ifndef SELECT_MPFQ_LAYER_u64k1
                /* this is all trivial in characteristic two, of course
                 */
                int rc = abinv(ab, tmp, M.coeff(u, r, 0));
                if (!rc) {
                    fprintf(stderr, "Error, found a factor of the modulus: ");
                    abfprint(ab, stderr, tmp);
                    fprintf(stderr, "\n");
                    exit(EXIT_FAILURE);
                }
                abneg(ab, tmp, tmp);
                for(unsigned int i = 0 ; i < m ; i++) {
                    abmul(ab, M.coeff(i, r, 0),
                              M.coeff(i, r, 0),
                              tmp);
                }
#endif

                r++;

                // if (r == m)
                    printf
                        ("[X^%d] A, col %d increases rank to %d (head row %d)%s\n",
                         k - (j < bm.d.nrhs), j, r, u,
                         (j < bm.d.nrhs) ? " (column not shifted because of the RHS)":"");
            }
        }

        if (r != m) {
            printf("This amount of data is insufficient. "
                   "Cannot find %u independent cols within A\n", m);
            exit(EXIT_FAILURE);
        }

        /* t0 is the k value for the loop index when we exited the loop.
         */
        t0 = exponents[r - 1] + 1;

        /* Coefficients of degree up to t0-1 of A' are read. Unless
         * bm.d.nrhs == n, for at least one of the columns of A, this
         * means up to degree t0.
         */
        if (rank == 0)
            ASSERT_ALWAYS(bm_io::next_coeff_to_consume == t0 + (bm.d.nrhs < n));

        printf("Found satisfactory init data for t0=%d\n", t0);

        /* We've also got some adjustments to make: room for one extra
         * coefficient is needed in A. Reading of further coefficients will
         * pickup where they stopped, and will always leave the last
         * t0+1+simd coefficients readable. */
        unsigned int window = simd + simd * iceildiv(t0 + 1, simd);
        A.realloc(window);
        A.zero_pad(window);

        for(unsigned int j = 0 ; j < m ; j++) {
            fdesc[j][0] = exponents[j];
            fdesc[j][1] = cnum[j];
            ASSERT_ALWAYS(exponents[j] < t0);
        }
        // free(pivots);
        // free(exponents);
        // free(cnum);
        // matpoly_clear(ab, M);
        abclear(ab, &tmp);
    }
    MPI_Bcast(fdesc, 2*m, MPI_UNSIGNED, 0, bm.com[0]);
    MPI_Bcast(&(t0), 1, MPI_UNSIGNED, 0, bm.com[0]);
    bm.set_t0(t0);
}				/*}}} */

template<class Writer>
void bm_io::compute_E(Writer& E, unsigned int expected, unsigned int allocated)/*{{{*/
{
    // F0 is exactly the n x n identity matrix, plus the
    // X^(s-exponent)e_{cnum} vectors. fdesc has the (exponent, cnum)
    // pairs
    bw_dimensions & d = bm.d;
    unsigned int m = d.m;
    unsigned int n = d.n;
    abdst_field ab = d.ab;
    int rank;
    MPI_Comm_rank(bm.com[0], &rank);

    unsigned int window = simd + simd * iceildiv(t0 + 1, simd);

    double tt0 = wct_seconds();
    double next_report_t = tt0 + 10;
    unsigned next_report_k = expected / 100;

    int over = 0;
    unsigned int final_length = allocated - t0;

    /* For the moment, despite the simd nature of the data, we'll compute
     * E one coeff at a time.
     */
    for(unsigned int kE = 0 ; kE + t0 < allocated ; kE++) {
        /* See the discussion in compute_initial_F ; to form coefficient
         * of degree kE of E, which is coefficient of degree t0+kE of
         * A'*F, we need to access the following coefficient of A:
         *
         *      t0+kE+1   if bm.d.nrhs < n
         *      t0+kE     if bm.d.nrhs = n
         *
         * more specifically, multiplying by row j of F will access
         * column j of A, and only coefficients of deg <= t0+kE+(j>=bm.d.nrhs)
         */
        
        unsigned int k_access = kE + t0 + (bm.d.nrhs < n);

        if (k_access % E.chunk_size() == 0 || kE + t0 >= expected) {
            /* check periodically */
            MPI_Bcast(&over, 1, MPI_INT, 0, bm.com[0]);
            if (over) break;
        }

        if (rank == 0) {
            if (!over && k_access >= next_coeff_to_fetch_from_source) {
                unsigned int b = fetch_more_from_source(window, simd);
                over = b < simd;
                if (over) {
                    printf("EOF met after reading %u coefficients\n", next_coeff_to_fetch_from_source);
                    final_length = kE + b;
                }
            }
            if (!over) {
                if (k_access >= next_coeff_to_consume) next_coeff_to_consume++;
                ASSERT_ALWAYS(k_access < next_coeff_to_consume);
                ASSERT_ALWAYS(k_access < next_coeff_to_fetch_from_source);
            }
        }

        /* This merely makes sure that coefficient E is writable: this
         * call may change the view window for E, and in
         * the case of an MPI run, this view window will be eventually
         * pushed to other nodes
         */
        E.coeff(0, 0, kE);

        if (rank || kE >= final_length)
            continue;

        if (kE + t0 > allocated) {
            fprintf(stderr, "Going way past guessed length%s ???\n", ascii ? " (more than 5%%)" : "");
        }

        for(unsigned int j = 0 ; j < n ; j++) {
            /* If the first columns of F are the identity matrix, then
             * in E we get data from coefficient kE+t0 in A', i.e.
             * coefficient of degree kE+t0+(j>=nrhs) in column j of A. More
             * generally, if it's x^q*identity, we read
             * coeficient of index kE + t0 + (j>=nrhs) - q.
             *
             * Note that we prefer to take q=0 anyway, since a
             * choice like q=t0 would create duplicate rows in E,
             * and that would be bad.
             */
            unsigned int kA = kE + t0 + (j >= bm.d.nrhs);
            ASSERT_ALWAYS(!rank || kA < next_coeff_to_consume);
            matpoly_extract_column(E, j, kE, A, j, kA % window);
        }

        for(unsigned int jE = n ; jE < m + n ; jE++) {
            unsigned int jA = fdesc[jE-n][1];
            unsigned int offset = fdesc[jE-n][0];
            unsigned int kA = kE + offset + (jA >= bm.d.nrhs);
            ASSERT_ALWAYS(!rank || kA < next_coeff_to_consume);
            matpoly_extract_column(E, jE, kE, A, jA, kA % window);
        }
        ASSERT_ALWAYS(!rank || k_access + 1 ==  next_coeff_to_consume);
        if (k_access > next_report_k) {
            double tt = wct_seconds();
            if (tt > next_report_t) {
                printf(
                        "Read %u coefficients (%.1f%%)"
                        " in %.1f s (%.1f MB/s)\n",
                        k_access, 100.0 * k_access / expected,
                        tt-tt0, k_access * average_matsize(ab, m, n, ascii) / (tt-tt0)/1.0e6);
                next_report_t = tt + 10;
                next_report_k = k_access + expected / 100;
            }
        }
    }
    MPI_Bcast(&final_length, 1, MPI_UNSIGNED, 0, bm.com[0]);
    E.finalize(final_length);
}/*}}}*/

template<typename T> struct matpoly_factory {};

#ifdef ENABLE_MPI_LINGEN
template<> struct matpoly_factory<bigmatpoly> {
    typedef bigmatpoly T;
    typedef bigmatpoly_producer_task producer_task;
    typedef bigmatpoly_consumer_task consumer_task;
    bigmatpoly_model model;
    matpoly_factory(MPI_Comm * comm, unsigned int m, unsigned int n) : model(comm, m, n) {}
    T init(abdst_field ab, unsigned int m, unsigned int n, int len) {
        return bigmatpoly(ab, model, m, n, len);
    }
    static T bw_lingen(bmstatus & bm, T & E) {
        return bw_biglingen_collective(bm, E);
    }
    static size_t capacity(T const & p) { return p.my_cell().capacity(); }
};
#endif

template<> struct matpoly_factory<matpoly> {
    typedef matpoly T;
    typedef matpoly_producer_task producer_task;
    typedef matpoly_consumer_task consumer_task;
    matpoly_factory() {}
    T init(abdst_field ab, unsigned int m, unsigned int n, int len) {
        return matpoly(ab, m, n, len);
    }
    static T bw_lingen(bmstatus & bm, T & E) {
        return bw_lingen_single(bm, E);
    }
    static size_t capacity(T const & p) { return p.capacity(); }
};

template<typename T, typename Sink>
void bm_io::output_flow(T & pi)
{
    unsigned int n = bm.d.n;

    matpoly::memory_guard dummy(SIZE_MAX);

    /* This object will store the rolling coefficients of the resulting
     * F. Since we compute coefficients on the fly, using F0 which has
     * degree t0, we need a memory of t0+1 coefficients in order to
     * always have one correct coefficient.
     *
     * Also, in the binary case where we want to store coefficients in
     * windows of 64, we need some extra room.
     */
    F = matpoly(bm.d.ab, n, n, t0 + 1);

    set_write_behind_size();

    Sink S(*this, F);

    typename matpoly_factory<T>::consumer_task pi_consumer(*this, pi);

    compute_final_F(S, pi_consumer);

    /* We need this because we want all our deallocation to happen before
     * the guard's dtor gets called */
    F = matpoly();
}


/*}}}*/
#endif
