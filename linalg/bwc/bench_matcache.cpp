#include "cado.h" // IWYU pragma: keep

/* This program is the simplest interface to the bare matrix
 * multiplication routine. It's meant to provide an easy way of benching,
 * and comparing, different matrix product implementations.
 *
 * It must be used with the INNER matrix, not the .info one. So either use
 * the vanilla cado format (before balance), or one of the .h<i>.v<j>
 * matrices as output by balance.
 */

#include <cstdio>
#include <cstdint>             // for uint64_t
#include <cstdlib>
#include <ctime>
#include <cerrno>
#include <climits>
#include <cstring>

#include <algorithm>
#include <string>
#include <memory>
#include <vector>
#include <utility>

#include <pthread.h>            // for pthread_mutex_lock, pthread_mutex_unlock

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "arith-cross.hpp"
#include "arith-generic.hpp"
#include "crc.h"        // cado_crc_lfsr
#include "cxx_mpz.hpp"    // for cxx_mpz
#include "gmp_aux.h"
#include "macros.h"
#include "matmul.hpp"
#include "matrix_u32.hpp"     // for matrix_u32
#include "params.h"
#include "portability.h" // asprintf // IWYU pragma: keep
#include "utils_cxx.hpp"        // for unique_ptr<FILE, delete_FILE>
#include "version_info.h" // cado_revision_string
#include "worker-threads.h"

static void usage()
{
    fmt::print(stderr,
            "Usage: ./bench [--impl <implementation>] [--tmax <time>] [--nmax <n_iter>] [--nchecks <number> | --nocheck] [-r|--rebuild] [-t|--transpose] [--nthreads <number>] [--cycles <frequency>] -- <file0> [<file1> ... ]\n");
    exit(1);
}


struct bench_args;

struct private_args {// {{{
    bench_args const * ba = nullptr;
    worker_threads_group * tg = nullptr;
    int tnum = -1;
    std::shared_ptr<matmul_interface> mm;
    arith_generic::owned_vector colvec;
    arith_generic::owned_vector rowvec;
    private_args() = default;
    private_args(worker_threads_group &, int, bench_args const &);
    void check();
    void fill_both_vectors_random(cxx_gmp_randstate & rstate) const;
    void fill_both_vectors_zero() const;
};// }}}

struct bench_args {// {{{
    std::unique_ptr<arith_generic> A;
    cxx_param_list & pl;
    std::string impl { "bucket" };
    int nthreads = 1;
    int nchecks = 4;
    bool transpose;
    int rebuild;
    cxx_mpz prime = 2;
    double freq = 1;
    const char * unit = nullptr;
    std::vector<std::string> mfiles;
    std::string source_vec;
    std::vector<private_args> p;
    worker_threads_group * tg;
    uint64_t ncoeffs_total = 0;
    double tmax = 100.0;
    unsigned int nmax = UINT_MAX;

    bench_args(bench_args const &) = delete;
    bench_args(bench_args &&) = delete;
    bench_args& operator=(bench_args const &) = delete;
    bench_args& operator=(bench_args &&) = delete;

    static void configure_aliases(cxx_param_list & pl) {
        param_list_configure_alias(pl, "--transpose", "-t");
        param_list_configure_alias(pl, "--rebuild", "-r");
    }

    static void configure_switches(cxx_param_list & pl) {
        param_list_configure_switch(pl, "--nocheck", nullptr);
        param_list_configure_switch(pl, "--transpose", nullptr);
        param_list_configure_switch(pl, "--rebuild", nullptr);
    }

    bench_args(cxx_param_list & pl, std::vector<std::string> const & mfiles)// {{{
        : pl(pl)
        , mfiles(mfiles)
    {
        param_list_parse(pl, "prime", prime);
        param_list_parse(pl, "nthreads", nthreads);
        if (param_list_parse(pl, "cycles", freq)) {
            unit = "cy/c";
        } else {
            unit = "ns/c";
        }
        param_list_parse(pl, "nchecks", nchecks);
        param_list_parse(pl, "tmax", tmax);
        param_list_parse(pl, "nmax", nmax);

        if (param_list_parse_switch(pl, "--nocheck"))
            nchecks = 0;

        param_list_parse(pl, "impl", impl);

        /* The api mandates that we set the desired value for nbys. Here,
         * that's not really our intent, since we really want to bench
         * the layer in its favorite working context. Most of the time,
         * setting nbys is pointless.
         */
        int nbys = 64;
        /* may leave nbys at the default value ! */
        param_list_parse(pl, "nbys", nbys);

        A.reset(arith_generic::instance(prime, nbys));

        fmt::print(stderr, "Using implementation \"{}\"\n", impl);
        transpose = param_list_parse_switch(pl, "--transpose");
        rebuild   = param_list_parse_switch(pl, "--rebuild");

        if (mfiles.size() != (size_t) nthreads) {
            fmt::print(stderr, "{} threads requested, but {} files given on the command line.\n", nthreads, mfiles.size());
            exit(EXIT_FAILURE);
        }

        p.resize(nthreads);
        tg = worker_threads_init(nthreads);
    }// }}}

    ~bench_args() { worker_threads_clear(tg); }

    void do_for_all_threads(void (*f)(worker_threads_group *, int, bench_args *)) {// {{{
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
        worker_threads_do(tg, (worker_func_t) f, this);
    }// }}}

    uint64_t get_ncoeffs_total();
    void display_per_thread_info() const;

    void fill_all_vectors_random(cxx_gmp_randstate & rstate)// {{{
    {
        for(auto & P : p)
            P.fill_both_vectors_random(rstate);
    }// }}}

    void do_simple_matmul_if_requested();

    void do_timing_run();
};// }}}

static void init_func(struct worker_threads_group * tg, int tnum, struct bench_args * ba)/*{{{*/
{
    ba->p[tnum] = private_args(*tg, tnum, *ba);
}/*}}}*/

static void check_func(struct worker_threads_group * tg MAYBE_UNUSED, int tnum, struct bench_args * ba)/*{{{*/
{
    ba->p[tnum].check();
}/*}}}*/

static void mul_func(struct worker_threads_group * tg MAYBE_UNUSED, int tnum, struct bench_args * ba)/*{{{*/
{
    private_args const & p(ba->p[tnum]);
    auto * dstvec = (ba->transpose ? p.colvec : p.rowvec).get();
    auto * srcvec = (ba->transpose ? p.rowvec : p.colvec).get();
    p.mm->mul(dstvec, srcvec, !ba->transpose);
}/*}}}*/

static void clear_func(struct worker_threads_group * tg MAYBE_UNUSED, int tnum, struct bench_args * ba)/*{{{*/
{
    private_args const & p(ba->p[tnum]);
    p.mm->report(ba->freq);

    pthread_mutex_lock(&tg->mu);
    fmt::print("{}\n", p.mm->report_string);
    pthread_mutex_unlock(&tg->mu);
}/*}}}*/

uint64_t bench_args::get_ncoeffs_total()// {{{
{
    ncoeffs_total = 0;
    for(auto const & P : p) {
        ncoeffs_total += P.mm->ncoeffs;
    }
    return ncoeffs_total;
}// }}}

void bench_args::display_per_thread_info() const// {{{
{
    uint64_t ncoeffs_total = 0;
    for(auto const & P : p) {
        fmt::print (stderr, "T{} {}: {} rows {} cols {} coeffs\n",
                P.tnum, mfiles[P.tnum],
                P.mm->dim[0],
                P.mm->dim[1],
                P.mm->ncoeffs);
    }
    fmt::print (stderr, "total {} coeffs\n", ncoeffs_total);
}// }}}

void bench_args::do_simple_matmul_if_requested()// {{{
{
    std::string srcvecname;
    if (!param_list_parse(pl, "srcvec", srcvecname))
        return;

    private_args const & P(p[0]);

    if (nthreads > 1) {
        fmt::print(stderr, "srcvec incompatible with multithread\n");
        exit(EXIT_FAILURE);
    }

    P.fill_both_vectors_zero();

    std::unique_ptr<FILE, delete_FILE> f(fopen(srcvecname.c_str(), "rb"));
    if (!f) {
        fmt::print(stderr, "fopen({}): {}\n", srcvecname, strerror(errno));
        exit(EXIT_FAILURE);
    }
    /* with transpose, we're doing vector times matrix (rowvec
     * times matrix -> colvec) */
    arith_generic::elt * srcvec = (transpose ? P.rowvec : P.colvec).get();
    arith_generic::elt * dstvec = (transpose ? P.colvec : P.rowvec).get();
    unsigned int const n = P.mm->dim[1 ^ transpose];
    fmt::print(stderr, "reading {} bytes from {}\n",
            n * sizeof(uint64_t), srcvecname);
    size_t const nread = fread(srcvec, sizeof(uint64_t), n, f.get());
    if (nread != size_t(n)) {
        fmt::print(stderr, "short read ({} < {})\n", nread, n);
        exit(1);
    }

    do_for_all_threads(mul_func);

    std::string const dstvecname = srcvecname + ".dst";
    f.reset(fopen(dstvecname.c_str(), "wb"));
    unsigned int const nw = P.mm->dim[0 ^ transpose];
    fmt::print(stderr, "writing {} bytes to {}\n",
            nw * sizeof(uint64_t), dstvecname);
    size_t const nwritten = fwrite(dstvec, sizeof(uint64_t), nw, f.get());
    if (nwritten != size_t(nw)) {
        fmt::print(stderr, "short write ({} < {})\n", nwritten, nw);
        exit(1);
    }
    fmt::print(stderr, "Saved [{}]{} * [{}] to [{}]\n",
            P.mm->cachefile_name,
            transpose ? "^T" : "",
            srcvecname,
            dstvecname);
}// }}}

void bench_args::do_timing_run()// {{{
{
#define NLAST   10
    double last[10]={0,};
    double sum_last = 0;

    clock_t const t0 = clock();
    clock_t next = 0.25 * CLOCKS_PER_SEC;
    clock_t t1 = 0;

    auto timecap = clock_t(tmax * CLOCKS_PER_SEC);

    next = std::min(next, timecap);

    fmt::print("Note: timings in seconds per iterations are given in cpu seconds ;\n");
    fmt::print("Note: with {} threads, this means {:.2f} walltime seconds.\n", nthreads, 1.0 / nthreads);
    for(unsigned int n = 0 ; n < nmax ; ) {
        do_for_all_threads(mul_func);

        clock_t const t = clock();
        clock_t dt = t - t0;
        sum_last -= last[n % NLAST];
        last[n % NLAST] = double(t - t1) / CLOCKS_PER_SEC;
        sum_last += double(t - t1) / CLOCKS_PER_SEC;
        t1 = t;
        unsigned int nlast = NLAST;
        n++;
        nlast = std::min(nlast, n);
        if (dt > next || n == nmax - 1) {
            do { next += 0.25 * CLOCKS_PER_SEC; } while (dt > next);
            next = std::min(next, timecap);
            dt /= CLOCKS_PER_SEC;

            auto fromstart = fmt::format("{} iters in {}s",
                    n, dt / CLOCKS_PER_SEC);

            auto average = fmt::format("{:.3f}/1, {:.2f} {}",
                    double(dt) / CLOCKS_PER_SEC / n,
                    freq * 1.0e9 * double(dt)/n/double(ncoeffs_total), unit);
            auto window = fmt::format("last {} : {:.3f}/1, {:.2f} {}",
                    nlast,
                    sum_last / nlast,
                    freq * 1.0e9 * sum_last / nlast/double(ncoeffs_total), unit);
            fmt::print("{}, {} ({})        \r", fromstart, average, window);
            if (dt > timecap)
                break;
        }
    }
    fmt::print("\n");
}// }}}

private_args::private_args(worker_threads_group & tg, int tnum, bench_args const & ba)// {{{
    : ba(&ba)
    , tg(&tg)
    , tnum(tnum)
    , mm(matmul_interface::create(
                ba.A.get(), 0, 0, ba.mfiles[tnum], ba.impl, ba.pl, !ba.transpose))
{
    time_t t0 = time(nullptr);
    fmt::print(stderr, "Expect to find or create cache file {}\n", mm->cachefile_name);
    if (!ba.rebuild && mm->reload_cache()) {
        /* reload_cache will set the matrix dimensions in
         * (matmul_public&) mm
         */
        pthread_mutex_lock(&tg.mu);
        fmt::print(stderr, "T{} Reusing cache file for {}\n",
                tnum, ba.mfiles[tnum]);
        fmt::print(stderr, "T{} Cache load time {} wct\n",
                tnum, time(nullptr) - t0);
        pthread_mutex_unlock(&tg.mu);
    } else {
        clock_t const ct0 = clock();
        pthread_mutex_lock(&tg.mu);
        fmt::print(stderr, "T{} Building cache file for {}\n",
                tnum, ba.mfiles[tnum]);
        pthread_mutex_unlock(&tg.mu);
        ASSERT_ALWAYS(mm->store_transposed == ba.transpose);

        /* TODO: do we need to pass ba.transpose all the way down to
         * matmul_build_cache? I think it should be okayish now. */
        auto matrix = matrix_u32::from_file(
                ba.mfiles[tnum],
                matrix_u32::transpose_option { ba.transpose },
                matrix_u32::withcoeffs_option { !ba.A->is_characteristic_two() });

        mm->dim = std::get<1>(matrix);
        mm->ncoeffs = std::get<2>(matrix);

        /* maybe the interface to build_cache would need to require the
         * dimensions, after all */
        mm->build_cache(std::move(std::get<0>(matrix)));

        pthread_mutex_lock(&tg.mu);
        fmt::print(stderr, "T{} Cache build time {:.2f}s cpu\n",
                tnum, (double) (clock()-ct0) / CLOCKS_PER_SEC);
        fmt::print(stderr, "T{} Saving cache file for {}\n",
                tnum, ba.mfiles[tnum]);
        pthread_mutex_unlock(&tg.mu);
        t0 = time(nullptr);
        mm->save_cache();
        pthread_mutex_lock(&tg.mu);
        fmt::print(stderr, "T{} Cache save time {} wct\n",
                tnum, time(nullptr) - t0);
        pthread_mutex_unlock(&tg.mu);
    }

    pthread_mutex_lock(&tg.mu);
    fmt::print(stderr, "T{} uses cache file {}\n",
            tnum,
            /* cache for mmt->locfile, */
            mm->cachefile_name);
    pthread_mutex_unlock(&tg.mu);

    unsigned int nr = mm->dim[0];
    unsigned int nc = mm->dim[1];
    mm->aux(MATMUL_AUX_GET_READAHEAD, &nr);
    mm->aux(MATMUL_AUX_GET_READAHEAD, &nc);

    colvec = ba.A->alloc_vector(nc, ALIGNMENT_ON_ALL_BWC_VECTORS);
    rowvec = ba.A->alloc_vector(nr, ALIGNMENT_ON_ALL_BWC_VECTORS);
    ba.A->vec_set_zero(colvec.get(), nc);
    ba.A->vec_set_zero(rowvec.get(), nr);
}// }}}

static uint32_t crc32(arith_generic * A, arith_generic::owned_vector const & v, unsigned int n)// {{{
{
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
    return crc32((unsigned long*) v.get(), A->vec_elt_stride(n));
}// }}}

void private_args::check()// {{{
{
    arith_generic * A = ba->A.get();

    std::unique_ptr<arith_cross_generic> AxA(arith_cross_generic::instance(A,A));

    unsigned int nr = mm->dim[0];
    unsigned int nc = mm->dim[1];
    unsigned int const nr0 = nr;
    unsigned int const nc0 = nc;
    mm->aux(MATMUL_AUX_GET_READAHEAD, &nr);
    mm->aux(MATMUL_AUX_GET_READAHEAD, &nc);

    auto colvec_bis = A->alloc_vector(nc, ALIGNMENT_ON_ALL_BWC_VECTORS);
    auto rowvec_bis = A->alloc_vector(nr, ALIGNMENT_ON_ALL_BWC_VECTORS);
    A->vec_set(colvec_bis.get(), colvec.get(), nc);
    A->vec_set(rowvec_bis.get(), rowvec.get(), nr);

    auto check0 = A->alloc_vector(A->simd_groupsize(), ALIGNMENT_ON_ALL_BWC_VECTORS);
    auto check1 = A->alloc_vector(A->simd_groupsize(), ALIGNMENT_ON_ALL_BWC_VECTORS);

    /* See the comment in matmul_mul about the direction argument and the
     * number of coordinates of source/destination vectors */

    fmt::print("T{} colvec({}): {}\n", tnum, nc, crc32(A, colvec, nc0));
    mm->mul(rowvec.get(), colvec.get(), 1);
    fmt::print("T{} rowvec({}): {}\n", tnum, nr, crc32(A, rowvec, nr0));

    A->vec_set_zero(check0.get(), A->simd_groupsize());
    AxA->add_dotprod(check0.get(), rowvec.get(), rowvec_bis.get(), nr0);

    fmt::print("T{} rowvec_bis({}): {}\n", tnum, nr, crc32(A, rowvec_bis, nr0));
    mm->mul(colvec_bis.get(), rowvec_bis.get(), 0);
    fmt::print("T{} colvec_bis({}): {}\n", tnum, nc, crc32(A, colvec_bis, nc0));

    A->vec_set_zero(check1.get(), A->simd_groupsize());
    AxA->add_dotprod(check1.get(), colvec.get(), colvec_bis.get(), nc0);

    if (A->vec_cmp(check0.get(), check1.get(), A->simd_groupsize()) != 0) {
        pthread_mutex_lock(&tg->mu);
        fmt::print(stderr, "T{} : Check failed\n", tnum);
        pthread_mutex_unlock(&tg->mu);
        abort();
    }
    mm->aux(MATMUL_AUX_ZERO_STATS);
}// }}}

void private_args::fill_both_vectors_random(cxx_gmp_randstate & rstate) const// {{{
{
    arith_generic * A = ba->A.get();
    A->vec_set_random(colvec.get(), mm->dim[1], rstate);
    A->vec_set_random(rowvec.get(), mm->dim[0], rstate);
    /* If we want shared vectors, we have a bit of a problem here. We
     * would need to steal thread 0's vector, more or less like so.
     * But it's really ugly, right?
     * if (ba->transpose)
     * p->rowvec = ba->p[0].rowvec;
     * else
     * p->colvec = ba->p[0].colvec;
     */
}// }}}

void private_args::fill_both_vectors_zero() const// {{{
{
    arith_generic * A = ba->A.get();
    A->vec_set_zero(colvec.get(), mm->dim[1]);
    A->vec_set_zero(rowvec.get(), mm->dim[0]);
}// }}}

static void banner(int argc, char const * argv[])// {{{
{
    /* print command line */
    fmt::print (stderr, "# ({}) {}", cado_revision_string, (argv)[0]);
    for (int i = 1; i < (argc); i++)
        fmt::print (stderr, " {}", (argv)[i]);
    fmt::print (stderr, "\n");

#ifdef  __GNUC__
    fmt::print(stderr, "# Compiled with gcc " __VERSION__ "\n");
#endif
    fmt::print(stderr, "# Compilation flags " CFLAGS "\n");
}// }}}

int main(int argc, char const * argv[])
{
    setvbuf(stdout, nullptr, _IONBF, 0);
    setvbuf(stderr, nullptr, _IONBF, 0);

    banner(argc, argv);

    cxx_param_list pl;

    bench_args::configure_aliases(pl);
    bench_args::configure_switches(pl);

    /* {{{ */
    argv++,argc--;

    std::vector<std::string> mfiles;

    mfiles.reserve(argc);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        if (strcmp(argv[0],"--") == 0) {
            argv++, argc--;
            for( ; argc ; argv++, argc--)
                mfiles.emplace_back(argv[0]);
            break;
        }
        if (argv[0][0] != '-') {
            mfiles.emplace_back(argv[0]);
            argv++,argc--;
            continue;
        }
        fmt::print (stderr, "Unknown option: {}\n", argv[0]);
        usage();
    }

    bench_args ba(pl, mfiles);

    param_list_lookup_string(pl, "matmul_bucket_methods");
    param_list_lookup_string(pl, "l1_cache_size");
    param_list_lookup_string(pl, "l2_cache_size");
    param_list_lookup_string(pl, "srcvec");



    if (param_list_warn_unused(pl)) {
        usage();
    }
    /*}}}*/

    ba.do_for_all_threads(&init_func);


    /* {{{ display some info */

    ba.display_per_thread_info();

    /* }}} */

    cxx_gmp_randstate rstate;
    gmp_randseed_ui(rstate, 1);

    /* {{{ Do checks if requested */
    for(int t = 0; t < ba.nchecks ; t++) {
        /* create deterministic test values */
        ba.fill_all_vectors_random(rstate);

        ba.do_for_all_threads(check_func);
        fmt::print(stderr, "Check {} ok\n", t);
    }
    if (ba.nchecks)
        fmt::print("All {} checks passed\n", ba.nchecks);
    /* }}} */


    ba.fill_all_vectors_random(rstate);
    ba.do_simple_matmul_if_requested();
    

    ba.do_timing_run();
    // fmt::print("Scanned {} coeffs in total\n", n * ncoeffs_total);

    ba.do_for_all_threads(clear_func);

    return 0;
}
