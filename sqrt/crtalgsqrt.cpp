/*
 * Program: crtalgsqrt
 * Authors: E. Thomé.
 * Purpose: computing the squareroots and finishing the factorization
 *
 * The algorithm implemented here is described in section 4 of:
 *
 *      Emmanuel Thomé, Square root algorithms for the number field sieve,
 *      4th International Workshop on Arithmetic in Finite Fields - WAIFI 2012,
 *      Jul 2012, Bochum, Germany. pp.208-224,
 *      DOI : 10.1007/978-3-642-31662-3_15
 *      https://hal.inria.fr/hal-00756838
 */

/*
  Usage within CADO-NFS:
  1) run: crtalgsqrt -v -depfile c75.dep.000 -polyfile c75.poly
     and let "alg" be the last integer value printed, for example 271...279:
# [83.72] c7 (+++---++) -2 [271185940941113750637336882505937475494764983427230684073069569288946725279]
     (note that the programm also accepts a ratdepfile parameter,
     but this one unused for the moment).
  2) let "rat" be the value of the rational square root (given by sqrt)
  3) compute gcd(alg-rat, n) and gcd(alg+rat, n)
  4) if this fails, try another dependency
 */

/* TODO list.
 *
 * This program can be used as a replacement for the _algebraic_ part of
 * the square root. At the moment it doesn't do the rational root, which
 * is a pity.
 *
 * There are several outstanding things which could/should be done with
 * this program.
 *
 * For correctness:
 *
 * - finish binding with the rest: check also the rational square root,
 *   and produce the factorization. Not necessarily independent from the
 *   above, since the elementary check is whether we get x^2=a mod N.
 * - update this TODO list.
 *
 * For speed / memory:
 *
 * - The program is probably overzealous with mpz_mod's sometimes. Some
 *   can be s(h)aved.
 * - The peak memory usage is apparently quite high. It is linked to the
 *   fact that we may have several threads doing full-length
 *   multiplications. Notwithstanding this dominant effect for scratch
 *   space usage, the amount of useful data should not exceed 3 times the
 *   ram_gb parameter.
 * - parallelize multiply_all_shares. It's not hard. Presently it
 *   accounts for 40% of the WCT, which is not acceptable.
 *
 * Important in terms of functionality, but not critical:
 *
 * - Use fpLLL for solving the knapsack. Should make it possible to go to
 *   12 primes or so in degree 6 (at least).
 * - Be a lot more flexible with regard to MPI setup. In particular, it
 *   should be possible to fix the s,t parameters not exactly equal to
 *   their computed minimal value.
 *
 * Not important:
 *
 * - understand and, if relevant, fix the memory leak diagnosed by
 *   valgrind. I don't understand.
 */

#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>
#include <cstdarg>
#include <cctype>
#include <cerrno>

#include <algorithm>
#include <complex>
#include <utility>
#include <vector>

#include <unistd.h>
#include <sys/stat.h>
#include <gmp.h>
#include "fmt/base.h"
#include "fmt/std.h"

#include "double_poly.h"
#include "polynomial.hpp"

#include "barrier.h"
#include "cado_poly.h"
#include "gmp-hacks.h"
#include "gmp_aux.h"
#include "knapsack.h"
#include "abfiles.hpp"
#include "macros.h"
#include "misc.h"
#include "arith/mod_ul.h"
#include "arith/modul_poly.h"
#include "mpz_poly.h"
#include "params.h"
#include "powers_of_p.hpp"
#include "rootfinder.h"
#include "select_mpi.h"
#include "timing.h"
#include "version_info.h"
#include "portability.h"
#include "utils_cxx.hpp"
#include "cxx_mpz.hpp"
#include "number_context.hpp"
#include "mpi_proxies.hpp"
#include "sqrt_wq.hpp"
#include "sqrt_cachefiles.hpp"
#include "cado_math_aux.hpp"
#include "cado_mp_conversions.hpp"

using cado_mpi::mpi_data_agrees;
using cado_mpi::allreduce;
using cado_mpi::broadcast;

/* {{{ time */
static double program_starttime;

#define WCT     (wct_seconds() - program_starttime)

#define STOPWATCH_DECL	        					\
    double t0 MAYBE_UNUSED, t1 MAYBE_UNUSED;				\
    double w0 MAYBE_UNUSED, w1 MAYBE_UNUSED;				\
    double rate MAYBE_UNUSED

#define STOPWATCH_GO()	        					\
    t0 = seconds();							\
    w0 = WCT;

#define STOPWATCH_GET()							\
    t1 = seconds();							\
    w1 = WCT;								\
    rate = (w1-w0) ? ((t1-t0)/(w1-w0) * 100.0) : 0;

#define log_step(step) logprint("%s%s\n", __func__, step)
#define log_begin() log_step(" starts")
#define log_step_time(step)                                             \
    logprint("%s%s. %.2lf, wct %.2lf, %.1f%%\n",                        \
            __func__, step, t1-t0, w1-w0, rate)
#define log_end() log_step_time(" ends")

/* }}} */

static int verbose = 0;

static double print_delay = 1;
static double ram_gb = 3.0;    // Number of gigabytes. Note that this is the
// maximum size of product that are obtained.
// The actual memory footprint can be quite
// considerably larger due to FFT allocation (by
// a constant factor, though).

static void usage()
{
    fprintf(stderr, "usage: crtalgsqrt algdepfile ratdepfile polyfile\n");
    exit(1);
}

/* {{{ logging */
static int max_loglevel=99;
static char prefix[20]={'\0'};

int
#ifndef HAVE_MINGW
/* Don't check format under MinGW as it still contains %zu here */
ATTR_PRINTF(1, 2)
#endif
logprint(const char * fmt, ...)
{
    va_list ap;
    int level=0;
    int s = strlen(fmt);
    const char * pfmt = fmt;
    if (s >= 3 && fmt[0] == '<' && isdigit(fmt[1]) && fmt[2]=='>') {
        level=fmt[1]-'0';
        pfmt += 3;
    }
    if (level > max_loglevel)
        return 0;

    static pthread_mutex_t obuf_lock = PTHREAD_MUTEX_INITIALIZER;
    static size_t st = 0;
    static char * t = NULL;
    pthread_mutex_lock(&obuf_lock);

    size_t wt = strlen(prefix) + strlen(fmt) + 80;
    if (wt > st) {
        checked_realloc(t, wt);
        st = wt;
    }
    snprintf(t, st, "# [%2.2lf] %s%s", WCT, prefix, pfmt);

    va_start(ap, fmt);
    int rc = vprintf(t, ap);
    va_end(ap);

    pthread_mutex_unlock(&obuf_lock);

    return rc;
}
/* }}} */

/* {{{ wrappers for some gmp operations, so as to report timings */
// above this threshold, we report each multiplication we do.
#define MUL_REPORT_THRESHOLD    8000000

#define REPORT_THIS(na, nb)     \
    (((na) > 10) && ((nb) > 10) && ((na) + (nb) > MUL_REPORT_THRESHOLD))

static void WRAP_mpz_mul(mpz_ptr c, mpz_srcptr a, mpz_srcptr b)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(b);
    mpz_mul(c,a,b);
    STOPWATCH_GET();
    if (REPORT_THIS(na, nb)) {
        logprint("<9> mpz_mul %d %d (%.1f) %.1f %.1f (%.1f%%)\n",
                (int) na, (int) nb, (double)na/nb, t1-t0, w1-w0, rate);
    }
}

#if 0
static void WRAP_mpz_invert(mpz_ptr c, mpz_srcptr a, mpz_srcptr b)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(b);
    mpz_invert(c,a,b);
    STOPWATCH_GET();
    if (REPORT_THIS(na, nb)) {
        logprint("<9> mpz_inv %d %d (%.1f) %.1f %.1f (%.1f%%)\n",
                (int) na, (int) nb, (double)na/nb, t1-t0, w1-w0, rate);
    }
}
#endif

#if 0 /* unused */
static void WRAP_mpz_addmul(mpz_ptr c, mpz_srcptr a, mpz_srcptr b)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(b);
    mpz_addmul(c,a,b);
    STOPWATCH_GET();
    if (REPORT_THIS(na, nb)) {
        logprint("<9> mpz_mul %d %d (%.1f) %.1f %.1f (%.1f%%)\n",
                (int) na, (int) nb, (double)na/nb, t1-t0, w1-w0, rate);
    }
}
#endif

static void WRAP_mpz_submul(mpz_ptr c, mpz_srcptr a, mpz_srcptr b)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(b);
    mpz_submul(c,a,b);
    STOPWATCH_GET();
    if (REPORT_THIS(na, nb)) {
        logprint("<9> mpz_mul %d %d (%.1f) %.1f %.1f (%.1f%%)\n",
                (int) na, (int) nb, (double)na/nb, t1-t0, w1-w0, rate);
    }
}

static void WRAP_mpz_mod(mpz_ptr c, mpz_srcptr a, mpz_srcptr p)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    mp_size_t na = mpz_size(a);
    mp_size_t nb = mpz_size(p);
    mpz_mod(c,a,p);
    STOPWATCH_GET();
    if (REPORT_THIS(na, nb) && na > nb + 10) {
        logprint("<9> mpz_mod %d %d (%.1f) %.1f %.1f (%.1f%%)\n",
                (int) na, (int) nb, (double)na/nb, t1-t0, w1-w0, rate);
    }
}
/* }}} */

/* {{{ mpi-gmp helpers */

/* {{{ share_1mpz -- two nodes communicate with eachother. integer z is
 * present on input for both nodes. Node s sends it to node r, which of
 * course receives it in integer t. The integer t is unused by node s. It
 * is an error to call this routine if the current node is not one of
 * {r,s}
 */
static void share_1mpz(mpz_ptr z, mpz_ptr t, int r, int s, MPI_Comm comm)
{
    // r <-- s
    int k;
    MPI_Comm_rank(comm, &k);
    mp_size_t nlimbs = mpz_size(z);

    if (k == r) {
        MPI_Recv(&nlimbs, 1, CADO_MPI_MP_SIZE_T, s, (r<<4), comm, MPI_STATUS_IGNORE);
        if (nlimbs == 0) {
            return;
        }

        _mpz_realloc(t, nlimbs);
        MPI_Recv(t->_mp_d, nlimbs, CADO_MPI_MP_LIMB_T, s, 1+(r<<4), comm, MPI_STATUS_IGNORE);
        MPI_Recv(&(t->_mp_size), 1, CADO_MPI_MPZ_INTERNAL_SIZE_T, s, 2+(r<<4), comm, MPI_STATUS_IGNORE);
    } else {
        MPI_Send(&nlimbs, 1, CADO_MPI_MP_SIZE_T, r, (r<<4), comm);
        if (nlimbs == 0)
            return;

        MPI_Send(z->_mp_d, nlimbs, CADO_MPI_MP_LIMB_T, r, 1+(r<<4), comm);
        MPI_Send(&(z->_mp_size), 1, CADO_MPI_MPZ_INTERNAL_SIZE_T, r, 2+(r<<4), comm);
    }
}/*}}}*/
#if 0
/* {{{ share_2mpz is the more complete version, where t is filled with
 * the peer integer on both nodes */
static void share_2mpz(mpz_ptr z, mpz_ptr t, int r, int s, MPI_Comm comm)
{
    int k;
    MPI_Comm_rank(comm, &k);
    int peer = r^s^k;

    mp_size_t nlimbs = mpz_size(z);
    mp_size_t nlimbs_peer = 0;

    MPI_Sendrecv(&nlimbs, 1, CADO_MPI_MP_SIZE_T, peer, k<<4,
            &nlimbs_peer, 1, CADO_MPI_MP_SIZE_T, peer, peer<<4,
            comm, MPI_STATUS_IGNORE);
    if (nlimbs == 0 || nlimbs_peer == 0) {
        return;
    }
    _mpz_realloc(t, nlimbs_peer);
    MPI_Sendrecv(z->_mp_d, nlimbs, CADO_MPI_MP_LIMB_T, peer, 1 + (k<<4),
            t->_mp_d, nlimbs_peer, CADO_MPI_MP_LIMB_T, peer, 1 + (peer<<4),
            comm, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&(z->_mp_size), 1, CADO_MPI_MPZ_INTERNAL_SIZE_T, peer, 2 + (k<<4),
            &(t->_mp_size), 1, CADO_MPI_MPZ_INTERNAL_SIZE_T, peer, 2 + (peer<<4),
            comm, MPI_STATUS_IGNORE);
}
/* }}} */
#endif
void reduce_mulmod_mpz(mpz_ptr z, int recv, MPI_Comm comm, mpz_srcptr px)/*{{{*/
{
    int me;
    int s;
    MPI_Comm_rank(comm, &me);
    MPI_Comm_size(comm, &s);
    ASSERT_ALWAYS(recv == 0);
    mpz_t t;
    mpz_init(t);
    for(int done_size = 1 ; done_size < s ; done_size <<= 1) {
        // XXX try to release the barrier, I think it's unnecessary.
        // MPI_Barrier(comm);
        if (me % done_size) continue;
        int receiver = me & ~done_size;
        int sender = me | done_size;
        if (sender >= s)
            continue;
        share_1mpz(z, t, receiver, sender, comm);
        // could be off-loaded to a working queue.
        if (me == receiver) {
            if (!mpz_size(t) || !mpz_size(z)) {
                /* It's not fundamentally wrong, but in our case, it
                 * indicates a bug, for sure. */
                fprintf(stderr, "Warning: received 0 from job %d\n", s);
                mpz_set_ui(z,0);
            }
            WRAP_mpz_mul(z,z,t);
            WRAP_mpz_mod(z,z, px);
        }
    }
    mpz_clear(t);
}/*}}}*/
void allreduce_mulmod_mpz(mpz_ptr z, MPI_Comm comm, mpz_srcptr px)/*{{{*/
{
#if 0
    /* This works only when the communicator size is a power of two.
     * If it isn't, then some parts don't get the result. Examples:
     *
     * 4 : (1 2 3 4) --> (12 12 34 34) --> (1234 1234 1234 1234).
     * 3 : (1 2 3) --> (12 12 3) --> (123 12 123)
     *
     * So we would need 3 rounds for a communicator of size 3, while the
     * circular algorithm needs two rounds only.
     */
    int me;
    int s;
    MPI_Comm_rank(comm, &me);
    MPI_Comm_size(comm, &s);
    mpz_t t;
    mpz_init(t);
    for(int done_size = 1 ; done_size < s ; done_size <<= 1) {
        // XXX try to release the barrier, I think it's unnecessary.
        // MPI_Barrier(comm);
        int receiver = me & ~done_size;
        int sender = me | done_size;
        if (sender >= s)
            continue;
        share_2mpz(z, t, receiver, sender, comm);
        if (!mpz_size(t) || !mpz_size(z)) {
            /* It's not fundamentally wrong, but in our case, it
             * indicates a bug, for sure. */
            fprintf(stderr, "Warning: received 0 from job %d\n", s);
            mpz_set_ui(z,0);
        }
        // could be off-loaded to a working queue.
        WRAP_mpz_mul(z,z,t);
        WRAP_mpz_mod(z,z, px);
    }
    mpz_clear(t);
#else
    /* completely stupid, but works */
    reduce_mulmod_mpz(z, 0, comm, px);
    broadcast(z, 0, comm);
#endif
}/*}}}*/

void reduce_mpz_poly_mul_mod_f(mpz_poly P, int recv, MPI_Comm comm, mpz_poly F)/*{{{*/
{
    int me;
    int s;
    MPI_Comm_rank(comm, &me);
    MPI_Comm_size(comm, &s);
    ASSERT_ALWAYS(recv == 0);
    mpz_poly Q;
    mpz_poly_init(Q, F->deg-1);
    ASSERT_ALWAYS(P->deg < F->deg);
    for(int done_size = 1 ; done_size < s ; done_size <<= 1) {
        // XXX try to release the barrier, I think it's unnecessary.
        // MPI_Barrier(comm);
        if (me % done_size) continue;
        int receiver = me & ~done_size;
        int sender = me | done_size;
        if (sender >= s)
            continue;
        for(int j = 0 ; j < F->deg ; j++) {
            if (j > P->deg)
                mpz_set_ui(mpz_poly_coeff(Q, j),0);
            share_1mpz(mpz_poly_coeff(P, j), mpz_poly_coeff(Q, j), receiver, sender, comm);
        }
        // could be off-loaded to a working queue.
        if (me == receiver) {
            mpz_poly_cleandeg(Q, F->deg - 1);
            mpz_poly_mul_mod_f(P, P, Q, F);
        }
    }
    mpz_poly_clear(Q);
}/*}}}*/
static void broadcast_poly(mpz_poly P, int maxdeg, int root, MPI_Comm comm) /*{{{*/
{
    /* maxdeg must be <= all allocation degrees. */
    ASSERT_ALWAYS(maxdeg + 1 >= 0);
    ASSERT_ALWAYS(((unsigned int) maxdeg + 1) <= P->alloc);
    for(int j = 0 ; j < maxdeg + 1 ; j++) {
        mpz_ptr z = mpz_poly_coeff(P, j);
        if (j > P->deg)
            mpz_set_ui(z,0);
        broadcast(z, root, comm);
    }
    mpz_poly_cleandeg(P, maxdeg);
}/*}}}*/
void allreduce_mpz_poly_mul_mod_f(mpz_poly P, MPI_Comm comm, mpz_poly F)/*{{{*/
{
    reduce_mpz_poly_mul_mod_f(P, 0, comm, F);
    broadcast_poly(P, F->deg - 1, 0, comm);
}
/* }}} */

/* }}} */

// some global variables (sheeh !)

struct sqrt_globals {
    int m;      // nprimes
    int n;      // degree
    int s;      // number of shares of A
    int t;      // number of prime groups
    int r;      // number of primes in each prime group
    unsigned long prec;
    cxx_mpz P;    // prime product (not to the power prec)
    size_t nbits_sqrt;
    // size_t nbits_a;
    cxx_mpz_poly f_hat;
    cxx_mpz_poly f_hat_diff;
    double f_hat_coeffs;
    cxx_cado_poly cpoly;
    mpz_t root_m;
    ab_source ab;
    int lll_maxdim = 50;
    int rank;
    int nprocs;
    int ncores = 2;
    struct work_queue wq[1];
    MPI_Comm acomm;     // same share of A
    MPI_Comm pcomm;     // same sub-product tree
    int arank, asize;
    int prank, psize;
    barrier_t barrier[1];
};

static struct sqrt_globals glob;

// {{{ TODO: Now that the v field is gone, replace the polymodF layer.
// Here's the only fragments which need to remain.
static void
mpz_poly_from_ab_monic(cxx_mpz_poly & tmp, long a, unsigned long b) {
    tmp->deg = b != 0;
    mpz_set_ui (mpz_poly_coeff(tmp, 1), b);
    mpz_neg (mpz_poly_coeff(tmp, 1), mpz_poly_coeff_const(tmp, 1));
    mpz_set_si (mpz_poly_coeff(tmp, 0), a);
    mpz_mul(mpz_poly_coeff(tmp, 0), mpz_poly_coeff_const(tmp, 0), mpz_poly_coeff_const(glob.cpoly->pols[1], glob.n));
}


// }}}

// {{{ floating point stuff
// {{{ getting the coefficients of the lagrange interpolation matrix.
polynomial<long double> lagrange_polynomial_abs(polynomial<cxx_mpz> const & f, std::complex<long double> r)
{
    auto q = f.div_q_xminusr(r) / f.derivative()(r);
    polynomial<long double> qa;
    for(int i = 0 ; i <= q.degree() ; i++)
        qa[i] = std::abs(decltype(r)(q[i]));
    return qa;
}
// }}}

void estimate_nbits_sqrt(size_t * sbits, ab_source_ptr ab) // , int guess)
{
    size_t abits[1];
    /*
       if (guess) {
     *abits = ab->digitbytes_estim * M_LN10 / M_LN2;
     *sbits = *abits / 2;
    // when doing gross estimates like this, we can hardly avoid
    // taking a safety margin.
     *abits += *abits / 10;
     *sbits += *sbits / 10;
     logprint("coefficients of A"
     " have at most %zu bits (ESTIMATED)\n", WCT, *abits);
     logprint("square root coefficients"
     " have at most %zu bits (ESTIMATED)\n", WCT, *sbits);
     return;
     }
     */

    int n = glob.cpoly->pols[1]->deg;

    /* gather roots of f and log|bigproduct(x)| at these roots */
    std::vector<std::pair<std::complex<long double>, double>>
        evals;

    evals.reserve(n);
    // take the roots of f, and multiply later on to obtain the roots of
    // f_hat. Otherwise we encounter precision issues.

    // {{{ compress the list of roots.
    int nreal = 0, ncomplex = 0;
    for(auto const & e : polynomial<cxx_mpz>(glob.cpoly->pols[1]).roots(cado::number_context<std::complex<long double>>())) {
        if (e.imag() > 0) { 
            evals.emplace_back(e, 0);
            ncomplex++;
        } else if (e.imag() < 0) {
            continue;
        } else {
            evals.emplace_back(e.real(), 0);
            nreal++;
        }
    }
    // }}}

    // post-scale to roots of f_hat, and store f_hat instead of f
    for(auto & [ x, logfx ] : evals)
        x *= mpz_get_d(mpz_poly_coeff_const(glob.cpoly->pols[1], n));
    
    // {{{ print the roots.
    if (nreal) {
        fmt::print("# [{:2.2f}]", WCT);
        fmt::print(" real");
        for(auto const & [ x, logfx ] : evals) {
            if (x.imag() == 0)
                fmt::print(" {}", x);
        }
    }
    if (ncomplex) {
        fmt::print(" complex");
        for(auto const & [ x, logfx ] : evals) {
            if (x.imag() > 0)
                fmt::print(" {}", x);
        }
    }
    fmt::print("\n");
    // }}}

    // {{{ now evaluate the product.

    // Consider the product A of all a-b\alpha's, which is a
    // polynomial in \alpha. This polynomial has _rational_ coefficients,
    // since each reduction involves dividing out by f_d.
    //
    // Easier to handle is f_d^nab*A (product of all f_d*a-b*f_d*\alpha),
    // which is an element of the order Z[f_d alpha] (f_d
    // alpha is an algebraic integer). It can be expressed with integer
    // coefficients in the powers of f_d\alpha. Its square root, however,
    // is not necessarily in this sub-order of the ring of integers.
    // Therefore we multiply by f_hat'(f_d\alpha), where f_hat is the
    // minimal polynomial of f_d\alpha.
    //
    // We ensure that we've rounded up nab to the next even multiple.

    int64_t a;
    uint64_t b;
    double w1,wt;
    w1 = WCT;
    ab_source_rewind(ab);
    for( ; ab_source_next(ab, &a, &b) ; ) {
        for(auto & [ x, logfx ] : evals) {
            std::complex<long double> y = a * mpz_get_d(mpz_poly_coeff_const(glob.cpoly->pols[1], n));
            std::complex<long double> w = x * b;
            y = y - w;
            logfx += cado_math_aux::log(cado_math_aux::abs(y));
        }
        wt = WCT;
        if (wt > w1 + print_delay || !(ab->nab % 10000000)) {
            w1 = wt;
            printf("# [%2.2lf] floating point evaluation: %zu (%.1f%%)\n",
                    WCT, ab->nab, 100.0*(double)ab->nab/ab->nab_estim);
        }
    }
    // }}}
    // note that now that we've read everything, we know the precise
    // number of (a,b)'s. Thus we can replace the estimation.
    ab->nab_estim = ab->nab;
// {{{ post-process evaluation: f'(alpha), and even nab. print.

    if (ab->nab & 1) {
        printf("# [%2.2lf] odd number of pairs !\n", WCT);
        for(auto & [ x, logfx] : evals)
            logfx += log(fabs(mpz_get_d(mpz_poly_coeff_const(glob.cpoly->pols[1], n))));
    }

    // multiply by the square of f_hat'(f_d\alpha).
    for(auto & [x, logfx] : evals) {
        auto s = polynomial<cxx_mpz>(glob.f_hat_diff)(x);
        logfx += 2 * std::log(s).real();
    }
    printf("# [%2.2lf] Log_2(A)", WCT);
    for(auto const & [x, logfx] : evals) {
        printf(" %.4g", logfx / M_LN2);
        if (x.imag() > 0)
            printf("*2");
    }
    printf("\n");
    // }}}

    // {{{ deduce the lognorm. print.
    double lognorm = 0;
    for(auto const & [x, logfx] : evals) {
        if (x.imag() > 0) {
            lognorm += 2*logfx;
        } else {
            lognorm += logfx;
        }
    }
    printf("# [%2.2lf] log_2(norm(A)) %.4g\n", WCT, lognorm / M_LN2);
    // }}}

    // {{{ now multiply this by the Lagrange matrix.
    std::vector<double> a_bounds(n, 0);
    std::vector<double> sqrt_bounds(n, 0);

    auto fz = polynomial<cxx_mpz>(glob.f_hat);
    for(auto const & [x, logfx] : evals) {
        auto q = lagrange_polynomial_abs(fz, x);
        for(int j = 0 ; j < n ; j++) {
            double za, zs;
            za = zs = log(q[j]);
            za += logfx;
            zs += logfx / 2;
            if (x.imag() > 0) {
                za += log(2);
                zs += log(2);
            }
            a_bounds[j] = std::max(a_bounds[j], za);
            sqrt_bounds[j] = std::max(sqrt_bounds[j], zs);
        }
    }
    // }}}

    // {{{ get global logbounds
    double logbound_sqrt = 0;
    double logbound_a = 0;
    for(int j = 0 ; j < n ; j++) {
        // note that we might have added up to n times the same thing.
        // (inequality a+b < 2max(a,b) )
        sqrt_bounds[j] += log(n);
        a_bounds[j] += log(n);

        // safety margin for inaccuracies ?
        sqrt_bounds[j] += 100 * M_LN2;
        a_bounds[j] += 100 * M_LN2;
        logbound_sqrt = std::max(logbound_sqrt, sqrt_bounds[j]);
        logbound_a = std::max(logbound_a, a_bounds[j]);
    }
    // }}}

#if 0
    // {{{ print logbounds
    printf("# [%2.2lf] logbounds per coeff of A", WCT);
    for(int j = 0 ; j < n ; j++) {
        printf(" %.4Lg", a_bounds[j] / M_LN2);
    }
    printf("\n");

    printf("# [%2.2lf] logbounds per coeff of sqrt(A)", WCT);
    for(int j = 0 ; j < n ; j++) {
        printf(" %.4Lg", sqrt_bounds[j] / M_LN2);
    }
    printf("\n");
    // }}}
#endif

    // Since we're being very gross, we use a stupid bound on the last
    // coefficient.

    *sbits = ceil(logbound_sqrt / M_LN2);
    *abits = ceil(logbound_a / M_LN2);
    printf("# [%2.2lf] coefficients of A"
            " have at most %zu bits\n", WCT, *abits);
    printf("# [%2.2lf] square root coefficients"
            " have at most %zu bits\n", WCT, *sbits);
}
/* }}} */

/*{{{ have to pre-declare prime_data */

#if 0
struct individual_contribution {
    uint64_t ratio;
    mpz_t modN;
};
#endif

struct prime_data {
    unsigned long p = 0;
    std::vector<unsigned long> r;
    // unsigned long rj;
    power_lookup_table powers;

    // computed somewhat late.
    cxx_mpz iHx;

    // struct alg_ptree_s * T = nullptr;     // always NULL.

    // after the rational ptree reduction, this contains the share of A
    // with coefficients reduced. Short-lived.
    cxx_mpz_poly A;

    std::vector<cxx_mpz> evals;       // set of evaluations.
    std::vector<cxx_mpz> lroots;      // (lifted) roots.
    std::vector<cxx_mpz> sqrts;       // square roots of A(x)     (only in the end !)

    prime_data(unsigned long p, std::vector<unsigned long> r)
        : p(p)
        , r(std::move(r))
        , powers(p)
        , evals(glob.n)
        , lroots(glob.n)
        , sqrts(glob.n)

    {
    }
};/* }}} */

/* {{{ product tree (algebraic) */
struct alg_ptree_s {
    struct alg_ptree_s * t0;
    struct alg_ptree_s * t1;
    cxx_mpz_poly s;
};
using alg_ptree_t = struct alg_ptree_s;

#if 0
alg_ptree_t * alg_ptree_build(struct prime_data * p, int i0, int i1)
{
    /* Everything being done the naive way, the algebraic ptree wins
     * nothing: the count of multiplications is exactly the same (n^2-n),
     * and the algebraic ptree incurs an extra storage for n
     * coefficients. It's thus disabled.
     */
    return NULL;

    mpz_srcptr px = p->powers(glob.prec);
    ASSERT_ALWAYS(i0 < i1);
    alg_ptree_t * res = (alg_ptree_t *) malloc(sizeof(alg_ptree_t));
    memset(res, 0, sizeof(alg_ptree_t));
    mpz_poly_init(res->s, i1-i0);
    if (i1-i0 == 1) {
        res->s->deg = 1;
        mpz_set_ui(res->s->coeff[1], 1);
        mpz_neg(res->s->coeff[0], p->rx);
        return res;
    }
    int d = (i1-i0)/2;
    res->t0 = alg_ptree_build(p, i0, i0+d);
    res->t1 = alg_ptree_build(p, i0+d, i1);
    mpz_poly_mul(res->s, res->t0->s, res->t1->s);
    mpz_poly_div_q_z(res->s, res->s, glob.F);
    mpz_poly_mod_mpz (res->s, res->s, px);
    return res;
}
#endif

void alg_ptree_clear(alg_ptree_t * T)
{
    if (T == NULL) return;
    alg_ptree_clear(T->t0);
    alg_ptree_clear(T->t1);
    free(T);
}

/* }}} */

/* {{{ finding the CRT primes */

std::vector<prime_data> suitable_crt_primes()
{
    unsigned int m = glob.m;

    std::vector<prime_data> res;
    res.reserve(m);

    // Note that p0 has to be well above the factor base bound. Thus on
    // 32-bit machines, it's possibly problematic.
    unsigned long p0 = ((~0UL)>>1)+1;    // 2^31 or 2^63
    unsigned long p = p0;

    if (glob.rank == 0)
        printf("# [%2.2lf] Searching for CRT primes\n", WCT);
    // printf("# [%2.2lf] p0=%lu\n", WCT, p0);

    cxx_gmp_randstate rstate;

    for( ; res.size() < m ; ) {
        p = ulong_nextprime(p);
        auto r = mpz_poly_roots(glob.cpoly->pols[1], p, rstate);
        if (r.size() != (size_t) glob.n) continue;
        for(auto & e : r) {
            cxx_mpz fd_e = mpz_poly_coeff_const(glob.cpoly->pols[1], glob.n);
            fd_e *= e;
            mpz_fdiv_r_ui(fd_e, fd_e, p);
            e = cado_math_aux::mpz_get<unsigned long>(fd_e);
        }
        std::ranges::sort(r);
        res.emplace_back(p, r);
        // printf("\n");
    }

    if (glob.rank == 0)
        printf("# [%2.2lf] Found all CRT primes\n", WCT);

    return res;
}
/* }}} */

/* {{{ everything that happens only modulo one prime */

/* {{{ tonelli-shanks */
void modul_find_ts_gen(residueul_t z, modulusul_t p)
{
    unsigned long pp = modul_getmod_ul(p)-1;
    int e = cado_ctzl(pp);
    pp >>= e;
    unsigned long s = 1UL << (e-1);
    residueul_t r;
    modul_init(r, p);
    do {
        modul_set_ul(z, rand (), p);
        modul_pow_ul(z, z, pp, p);
        modul_pow_ul(r, z, s, p);
        modul_add_ul(r, r, 1, p);
    } while (!modul_is0(r, p));
    modul_clear(r, p);
}

int modul_field_sqrt(residueul_t z, residueul_t a, residueul_t g, modulusul_t p)
{
    unsigned long pp = modul_getmod_ul(p)-1;
    int e = cado_ctzl(pp);
    pp >>= e+1;
    if (modul_is0(a, p)) {
        modul_set0(z, p);
        return 1;
    }
    if (modul_is0(g, p)) {
        modul_find_ts_gen(g, p);
    }
    residueul_t b, x, y, t;
    modul_init(b, p);
    modul_init(x, p);
    modul_init(y, p);
    modul_init(t, p);
    int r = e;
    unsigned long s = 1UL << (e-1);

    // modul_set(x, a, p);
    modul_set(y, g, p);

    modul_pow_ul(x, a, pp, p);
    modul_sqr(b, x, p);
    modul_mul(x, x, a, p);
    modul_mul(b, b, a, p);

    int m;
    for(;;) {
        modul_set(t, b, p);
        for(m=0; !modul_is1(t, p); m++)
            modul_sqr(t, t, p);
        ASSERT(m<=r);

        if (m==0 || m==r)
            break;

        s = 1UL << (r-m-1);
        r = m;

        modul_pow_ul(t, y, s, p);
        modul_sqr(y, t, p);
        modul_mul(x, x, t, p);
        modul_mul(b, b, y, p);
    }

    modul_set(z, x, p);
    modul_clear(t, p);
    modul_clear(x, p);
    modul_clear(y, p);
    modul_clear(b, p);
    return (m==0);
}
/* }}} */
/* {{{ sqrt / invsqrt */
void invsqrt_lift(struct prime_data * p, mpz_ptr A, mpz_ptr sx, int precision)
{
    double w0 = WCT;
    ASSERT(precision > 0);

    if (precision == 1) {
        gmp_printf("A is %Zd\n", A);
        logprint("computing a root of A mod %lu\n", p->p);
        residueul_t z, a;
        modulusul_t q;
        modul_initmod_ul(q, p->p);
        modul_init(z, q);
        modul_init(a, q);
        modul_set_ul_reduced(a, mpz_get_ui(A), q);
        int issquare = modul_field_sqrt(a, a, z, q);
        ASSERT_ALWAYS(issquare);
        modul_inv(a, a, q);
        if (modul_get_ul(a, q) & 1) modul_neg(a, a, q);
        mpz_set_ui(sx, modul_get_ul(a, q));
        modul_clear(a, q);
        modul_clear(z, q);
        modul_clearmod(q);
        return;
    }
    int lower = precision - precision / 2;

    mpz_srcptr pl = p->powers(lower);

    // we're going to recurse ; arrange so that we recurse on something
    // acceptably small.  The problem with this approach is that we store
    // many reductions in memory.
    mpz_t A_save;
    mpz_init(A_save);
    WRAP_mpz_mod(A_save, A, pl);
    // recurse.
    invsqrt_lift(p, A_save, sx, lower);
    mpz_clear(A_save);

    if (WCT > w0 + print_delay)
        logprint("precision %d\n", precision);

    mpz_srcptr pk = p->powers(precision);

    mpz_t ta;
    mpz_init(ta);

    WRAP_mpz_mul(ta, sx, sx);
    WRAP_mpz_mod(ta, ta, pk);
    WRAP_mpz_mul(ta, ta, A);
    WRAP_mpz_mod(ta, ta, pk);
    mpz_sub_ui(ta, ta, 1);
    if (mpz_odd_p(ta)) mpz_add(ta, ta, pk);
    mpz_div_2exp(ta, ta, 1);
    mpz_submul(sx, ta, sx);
    WRAP_mpz_mod(sx, sx, pk);

    mpz_clear(ta);
}

void sqrt_lift(struct prime_data * p, mpz_ptr A, mpz_ptr sx, int precision)
{
    double w0 = WCT;
    int lower = precision - precision / 2;
    mpz_srcptr pk = p->powers(precision);
    mpz_srcptr pl = p->powers(lower);
    mpz_t A_save;
    mpz_init(A_save);
    WRAP_mpz_mod(A_save, A, pl);
    invsqrt_lift(p, A_save, sx, lower);
    mpz_clear(A_save);

    if (WCT > w0 + print_delay)
        logprint("precision %d\n", precision);
    // inverse square root now in sx.

    mpz_t tmp;
    mpz_init(tmp);

    WRAP_mpz_mul(tmp, A, sx);
    WRAP_mpz_mod(tmp, tmp, pl);
    // XXX This destroys A !!!
    mpz_submul(A, tmp, tmp);
    if (mpz_odd_p(sx)) mpz_add(sx, sx, pk);
    mpz_div_2exp(sx, sx, 1);
    mpz_addmul(tmp, sx, A);
    WRAP_mpz_mod(sx, tmp, pk);

    mpz_clear(tmp);
}
/* }}} */

void root_lift(struct prime_data * p, mpz_ptr rx, mpz_ptr irx, int precision)/* {{{ */
{
    double w0 = WCT;
    ASSERT(precision > 0);

    if (precision == 1) {
        return;
    }
    int lower = precision - precision / 2;

    // recurse.
    root_lift(p, rx, irx, lower);

    if (WCT > w0 + print_delay)
        logprint("precision %d\n", precision);
    mpz_srcptr pk = p->powers(precision);
    mpz_srcptr pl = p->powers(lower);

    mpz_t ta, tb;
    mpz_init(ta);
    mpz_init(tb);

    mpz_t fprime;
    mpz_init(fprime);

    // we know r to half-precision, 1/f'(r) to quarter-precision. Compute
    // f(r) (to full precision), and obtain 1/f'(r) to half-precision. We
    // compute f'(r) to half-precision as a side-effect of the
    // computation of f(r).

    mpz_ptr rr[2] = { tb, fprime };
    mpz_poly_srcptr ff[2] = { glob.f_hat, glob.f_hat_diff };
    mpz_poly_eval_several_mod_mpz (rr, ff, 2, rx, pk);
    /* use irx. only one iteration of newton.  */
    WRAP_mpz_mod(fprime, fprime, pl);
    WRAP_mpz_mul(ta, irx, fprime);
    WRAP_mpz_mod(ta, ta, pl);
    mpz_sub_ui(ta,ta,1);
    WRAP_mpz_submul(irx, irx, ta);
    WRAP_mpz_mod(irx, irx, pl);

    mpz_clear(fprime);

    WRAP_mpz_mul(tb, irx, tb);
    mpz_sub(rx, rx, tb);
    WRAP_mpz_mod(rx, rx, pk);

    mpz_clear(ta);
    mpz_clear(tb);
}/* }}} */

/* }}} */

/* {{{ knapsack stuff */
struct crtalgsqrt_knapsack {
    knapsack_object ks;
    mpz_t Px;
    mpz_t lcx;
    mpz_t fhdiff_modN;
    mpz_t sqrt_modN;
    const mp_limb_t * tabN;
};

void crtalgsqrt_knapsack_init(struct crtalgsqrt_knapsack * cks)
{
    knapsack_object_init(cks->ks);
    mpz_init(cks->Px);
    mpz_init(cks->lcx);
    mpz_ptr z = cks->fhdiff_modN;
    mpz_init(z);
    mpz_init(cks->sqrt_modN);
}

int crtalgsqrt_knapsack_callback(struct crtalgsqrt_knapsack * cks,
        unsigned long v, int64_t x)
{
    knapsack_object_srcptr ks = cks->ks;
    const int64_t * c64 = cks->ks->tab;
    const mp_limb_t * cN = cks->tabN;
    unsigned int s;
    char * signs = (char *) malloc(ks->nelems+1);

    memset(signs, 0, ks->nelems+1);
    for(s = 0 ; s < ks->nelems ; s++) {
        signs[s]=(v & (1UL << s)) ? '+' : '-';
    }

    // now try to look at the solution more closely */
    mpz_t e;
    mpz_init_set_ui(e, 0);
    int spurious = 0;
    for(int k = glob.n - 1 ; k >= 0 ; k--) {
        mpz_t z, t, zN, qz, rz;
        mpz_init_set_ui(z, cks->ks->bound);
        mpz_init_set_ui(zN, 0);
        mpz_init_set_ui(t, 0);
        mpz_init(qz);
        mpz_init(rz);
        int64_t sk = ks->bound;
        for(int j = 0 ; j < glob.m * glob.n ; j++) {
            mpz_set_uint64(t, c64[j * glob.n + k]);
            mpz_t w;
            mp_size_t sN = mpz_size(glob.cpoly->n);
            MPZ_INIT_SET_MPN(w, cN + sN * (j * glob.n + k), sN);
            if (v & (1UL << j)) {
                sk+=c64[j * glob.n + k];
                mpz_add(z, z, t);
                mpz_add(zN, zN, w);
            } else {
                sk-=c64[j * glob.n + k];
                mpz_sub(z, z, t);
                mpz_sub(zN, zN, w);
            }
            mpz_clear(w);
        }
        if (sk >= 2 * ks->bound || sk < 0) {
            printf("[%s] recombination of coeff in X^%d yields noise (%" PRIu64 ")\n",signs, k, sk);
            spurious++;
            // break;
        }
        // so we have (sum of r's) mod p^k = something * p^k + small
        // the ``something'' is in the quotient of the division.
        // FIXME
        // The problem is that the ``small'' thing might be small
        // enough that we won't be able to tell apart s and p-s.
        // This implies that we will have to try 2^degree
        // combinations: those with the quotients as given, and
        // the other combinations with one or several quotients
        // lowered by one unit.
        mpz_set(qz,z);
        if (mpz_cmp_ui(qz,0) >= 0) {
            mpz_fdiv_q_2exp(qz, qz, 63);
            mpz_add_ui(qz,qz,mpz_odd_p(qz) != 0);
            mpz_fdiv_q_2exp(qz, qz, 1);
        } else {
            mpz_neg(qz,qz);
            mpz_fdiv_q_2exp(qz, qz, 63);
            mpz_add_ui(qz,qz,mpz_odd_p(qz) != 0);
            mpz_fdiv_q_2exp(qz, qz, 1);
            mpz_neg(qz,qz);
        }
        mpz_mul_2exp(rz, qz, 64);
        mpz_sub(rz, z, rz);
        // gmp_printf("[%d] %Zd = %Zd * 2^64 + %Zd\n", k, z, qz, rz);
        mpz_submul(zN, qz, cks->Px);
        mpz_mod(zN, zN, glob.cpoly->n);
        // gmp_printf("[X^%d] %Zd\n", k, zN);
        // good. we have the coefficient !
        mpz_mul(e, e, glob.root_m);
        mpz_mul(e, e, mpz_poly_coeff_const(glob.cpoly->pols[1], glob.n));
        mpz_add(e, e, zN);
        mpz_mod(e, e, glob.cpoly->n);
        mpz_clear(z);
        mpz_clear(t);
        mpz_clear(zN);
        mpz_clear(qz);
        mpz_clear(rz);

    }

    mpz_mul(e,e,cks->lcx);
    mpz_mod(e,e,glob.cpoly->n);

    mpz_mul(e,e,cks->fhdiff_modN);
    mpz_mod(e,e,glob.cpoly->n);

    mpz_set(cks->sqrt_modN, e);

    mpz_clear(e);

    if (spurious) {
        printf("# [%2.2lf] %lx (%s) %" PRId64 " [SPURIOUS]\n",
                WCT, v, signs, (int64_t) x);
    } else {
        gmp_printf("# [%2.2lf] %lx (%s) %" PRId64 " [%Zd]\n",
                WCT, v, signs, (int64_t) x, cks->sqrt_modN);
    }
    free(signs);

    return spurious == 0;
}

/* Compute all the stuff which will be needed for checking the solutions
*/
void crtalgsqrt_knapsack_prepare(struct crtalgsqrt_knapsack * cks, size_t lc_exp)
{
#if 0 /* {{{ display the knapsack contents */
    {
        for(int i = 0 ; i < glob.m ; i++) {
            for(int j =  0 ; j < glob.n ; j++) {
                for(int k = 0 ; k < glob.n ; k++) {
                    int64_t c = contribs64[ ((i * glob.n) + j) * glob.n + k];
                    printf(" %" PRId64, c);
                }
                for(int s = 0 ; s < glob.m * glob.n ; s++) {
                    printf(" %d", s == (i * glob.n + j));
                }
                printf("\n");
            }
        }
        mpz_t z;
        mpz_init(z);
        mpz_ui_pow_ui(z, 2, 64);
        for(int k = 0 ; k < glob.n ; k++) {
            for(int s = 0 ; s < glob.n + glob.m * glob.n ; s++) {
                if (s == k) {
                    gmp_printf(" %Zd", z);
                } else {
                    printf(" 0");
                }
            }
            printf("\n");
        }
        mpz_clear(z);
    }
#endif/*}}}*/

    mpz_powm_ui(cks->Px, glob.P, glob.prec, glob.cpoly->n);

    mpz_set(cks->lcx, mpz_poly_coeff_const(glob.cpoly->pols[1], glob.n));
    mpz_powm_ui(cks->lcx, mpz_poly_coeff_const(glob.cpoly->pols[1], glob.n), lc_exp, glob.cpoly->n);
    mpz_invert(cks->lcx, cks->lcx, glob.cpoly->n);

    // evaluate the derivative of f_hat in alpha_hat mod N, that is lc*m.
    {
        mpz_t alpha_hat;
        mpz_init(alpha_hat);
        mpz_mul(alpha_hat, glob.root_m, mpz_poly_coeff_const(glob.cpoly->pols[1], glob.n));
        mpz_poly_eval_mod_mpz(cks->fhdiff_modN,
                glob.f_hat_diff, alpha_hat, glob.cpoly->n);
        mpz_clear(alpha_hat);
    }
    mpz_invert(cks->fhdiff_modN, cks->fhdiff_modN, glob.cpoly->n);

    cks->ks->cb_arg = cks;
    cks->ks->cb = (knapsack_object_callback_t) crtalgsqrt_knapsack_callback;

    cks->ks->bound = 1+ceil(log(glob.m*glob.n)/log(2));

    cks->ks->stride = glob.n;

    unsigned int nelems = cks->ks->nelems = glob.m * glob.n;
    unsigned int k1 = nelems / 2;
    unsigned int k2 = nelems - k1;
    uint64_t n1 = UINT64_C(1) << k1;
    uint64_t n2 = UINT64_C(1) << k2;
    char buf[16];

    printf(
            "# [%2.2lf] Recombination: dimension %u = %u + %u, %s needed\n",
            WCT,
            nelems, k1, k2,
            size_disp((n1+n2) * sizeof(uint64_t), buf)
           );

}

void crtalgsqrt_knapsack_clear(struct crtalgsqrt_knapsack * cks)
{
    mpz_clear(cks->sqrt_modN);
    mpz_clear(cks->Px);
    mpz_clear(cks->lcx);
    mpz_clear(cks->fhdiff_modN);
    knapsack_object_clear(cks->ks);
}
/* }}} */

void mpi_custom_types_check()/*{{{*/
{
    int ret;
    /* If this ever fails, we have to do something to the defines above
     * (and maybe touch the cmake files) */
    MPI_Type_size(CADO_MPI_SIZE_T,&ret);   ASSERT_ALWAYS(ret==sizeof(size_t));
    MPI_Type_size(CADO_MPI_MP_SIZE_T,&ret);ASSERT_ALWAYS(ret==sizeof(mp_size_t));
    MPI_Type_size(CADO_MPI_MP_LIMB_T,&ret);ASSERT_ALWAYS(ret==sizeof(mp_limb_t));
    MPI_Type_size(CADO_MPI_UINT64_T,&ret);ASSERT_ALWAYS(ret==sizeof(uint64_t));
    MPI_Type_size(CADO_MPI_INT64_T,&ret);ASSERT_ALWAYS(ret==sizeof(int64_t));
    {
        mpz_t z;
        MPI_Type_size(CADO_MPI_MPZ_INTERNAL_SIZE_T,&ret);
        ASSERT_ALWAYS(ret==sizeof(z->_mp_size));
    }
}/*}}}*/

// {{{ parallelism parameters selection
// how many parts of A shall we consider ?
// now how do we choose the appropriate number of primes ? Let:
//
// R: the number of bits of the ram guide amount given.
// n: the degree
// T: the number of bits of one coefficient of the square root.
//
// Note thus that A occupies a memory space of 2Tn bits.
//
// and
//
// s: the number of parts of A considered
// B: the number of bits of lifted primes p^lambda
// r: the number of primes stored together in a product tree.
// t: the number of product trees considered.
//
// We have:
//
// [1] sR >= 2Tn (parts of A must form the totality of A)
// [2] nrB <= R  (sets of ``reduced'' A's must fit in RAM)
// [3] rtB >= T  (lifted primes must be large enough for the reconstruction)
//
// A further condition is given by the maximum knapsack lattice
// dimension which can be efficiently handled. This is somewhat
// arbitrary, let's say:
//
// [4] rtn <= 80
//
// Note that [1] readily gives s, and that [2+3] give t. Then per
// [2,3], the product rB is placed within an interval [T/t, R/n].
// Obeying [3], either by choosing the smallest/largest admissible
// power of two for r, or simply the smallest/largest admissible
// value, we deduce r and B. (large r will yield more
// parallelism-multithreading, but OTOH small r will yield less
// time).

static void get_parameters(int * pr, int * ps, int * pt, int asked_r)
{
    double R = ram_gb * 8UL * (double) (1UL << 30);
    int n = glob.n;
    size_t T = glob.nbits_sqrt + 80;    // some margin for knapsack.
    double ramfraction = 2*T*n / R;
    int s = ceil(ramfraction);
    int t = ceil((double)n*T / R);
    int dmax = glob.lll_maxdim - n;
    double rmax = (double) dmax / n / t;
    int r = rmax;

    if (r == 0) {
        r = 1;
        fprintf(stderr, "Warning: resorting to choosing r=1, which will give a reconstruction effort of dimension %d\n",
                n*t);
    }
    // for(int mask = ~0 ; r & mask ; mask <<= 1) r &= mask;
    // try r=1 always.
    if (asked_r > r) {
        printf("value %d requested for r is too large, using %d instead\n", asked_r, r);
    } else if (asked_r) {
        r = asked_r;
    } else {
        /* Do we provide a default ? */
        r = 1;
    }
    size_t B = ceil((double)T/t/r);

    char buf[16];
    char buf2[16];
    printf("# [%2.2lf] A is %s, %.1f%% of %s."
            " Splitting in %d parts\n",
            WCT,
            size_disp(2*T*n/8.0, buf),
            100.0 * ramfraction, size_disp(ram_gb*1024.0*1048576.0, buf2), s);

    printf("# [%2.2lf] %d groups of %d prime powers,"
            " each of %zu bits\n",
            WCT, t, r, B);

    /* If these asserts fail, it means that the constraints are
     * impossible to satisfy */
    ASSERT_ALWAYS((double)T/t < R/n);
    ASSERT_ALWAYS(B*n <= R);

    printf("# [%2.2lf] number of pairs is %zu\n", WCT,
            glob.ab->nab);
    *pr = r;
    *ps = s;
    *pt = t;
}
/* }}} */

static int rat_red_caches_ok(std::vector<prime_data> const & primes, int i0, int i1, size_t off0, size_t off1)/*{{{*/
{
    if (!rcache) return 0;
    int nok = 0;
    for(int i = i0 ; i < i1 ; i++) {
        ASSERT_ALWAYS(i < (int) primes.size());
        nok += cachefile_exists("a_%zu_%zu_mod_%lu_%lu", off0, off1, primes[i].p, glob.prec);
    }
    allreduce(nok, MPI_SUM, MPI_COMM_WORLD);
    return nok == glob.m * glob.s;
}/*}}}*/

static int alg_red_caches_ok(std::vector<prime_data> const & primes, int i0, int i1, size_t off0, size_t off1)/*{{{*/
{
    if (!rcache) return 0;
    int nok = 0;
    for(int i = i0 ; i < i1 ; i++) {
        ASSERT_ALWAYS(i < (int) primes.size());
        auto const & p = primes[i];
        for(int j = 0 ; j < glob.n ; j++) {
            nok += cachefile_exists("a_%zu_%zu_mod_%lu_%lu_%lu",
                    off0, off1, p.p, p.r[j], glob.prec);
        }
    }
    allreduce(nok, MPI_SUM, MPI_COMM_WORLD);
    return nok == glob.n * glob.m * glob.s;
}/*}}}*/

static int a_shared_caches_ok(std::vector<prime_data> const & primes, int i0, int i1)/*{{{*/
{
    if (!rcache) return 0;
    int nok = 0;
    for(int i = i0 ; i < i1 ; i++) {
        ASSERT_ALWAYS(i < (int) primes.size());
        auto const & p = primes[i];
        for(int j = 0 ; j < glob.n ; j++) {
            nok += cachefile_exists("a_mod_%lu_%lu_%lu", p.p, p.r[j], glob.prec);
        }
    }
    allreduce(nok, MPI_SUM, MPI_COMM_WORLD);
    return nok == glob.n * glob.m * glob.s;
}/*}}}*/

static int sqrt_caches_ok(std::vector<prime_data> const & primes, int i0, int i1)/*{{{*/
{
    if (!rcache) return 0;
    int nok = 0;
    for(int i = i0 ; i < i1 ; i++) {
        ASSERT_ALWAYS(i < (int) primes.size());
        auto const & p = primes[i];
        for(int j = 0 ; j < glob.n ; j++) {
            nok += cachefile_exists("sqrt_%lu_%lu_%lu", p.p, p.r[j], glob.prec);
        }
    }
    allreduce(nok, MPI_SUM, MPI_COMM_WORLD);
    return nok == glob.n * glob.m * glob.s;
}/*}}}*/

struct subtask_info_t {
    struct prime_data * p;
    int i0, i1;
    int j;
    size_t off0, off1;
    mpz_poly_ptr P;
    mpz_poly_ptr P0;
    mpz_poly_ptr P1;
    size_t nab_loc;
    struct wq_task * handle;
};

/*{{{ a_poly_read_share and companion */
#define ABPOLY_OFFSET_THRESHOLD        65536
// NOTE: This does not depend on p (nor r of course).
// the accumulation is done for all data between:
// the first data line starting at offset >= off0 (inclusive)
// the first data line starting at offset >= off1 (exclusive)
size_t accumulate_ab_poly(mpz_poly_ptr P, ab_source_ptr ab, size_t off0, size_t off1, cxx_mpz_poly & tmp)
{
    size_t res = 0;
    mpz_poly_set_ui(P, 1);
    if (off1 - off0 < ABPOLY_OFFSET_THRESHOLD) {
        ab_source_move_afterpos(ab, off0);
        logprint("<4> (a,b) rewind to %s, pos %zu\n",
                ab->nfiles ? ab->sname : ab->fname0, ab->cpos);
        for( ; ab->tpos < off1 ; res++) {
            int64_t a;
            uint64_t b;
            int r = ab_source_next(ab, &a, &b);
            FATAL_ERROR_CHECK(!r, "dep file ended prematurely\n");
            mpz_poly_from_ab_monic(tmp, a, b);
            mpz_poly_mul_mod_f(P, P, tmp, glob.f_hat);
        }
        return res;
    }
    size_t d = (off1 - off0) / 2;
    cxx_mpz_poly Pl, Pr;
    res += accumulate_ab_poly(Pl, ab, off0, off0 + d, tmp);
    res += accumulate_ab_poly(Pr, ab, off0 + d, off1, tmp);
    mpz_poly_mul_mod_f(P, Pl, Pr, glob.f_hat);
    return res;
}


void * a_poly_read_share_child(struct subtask_info_t * info)
{
    size_t off0 = info->off0;
    size_t off1 = info->off1;
        int rc;
    mpz_poly_ptr P = info->P;
    cachefile c;

    cachefile_init(c, "a_%zu_%zu", off0, off1);

    if (rcache && cachefile_open_r(c)) {
        logprint("reading cache %s\n", c->basename);
        rc = fscanf(c->f, "%zu", &info->nab_loc);
        ASSERT_ALWAYS(rc == 1);
        for(int i = 0 ; i < glob.n ; i++) {
            rc = gmp_fscanf(c->f, "%Zx", mpz_poly_coeff_const(info->P, i));
            ASSERT_ALWAYS(rc == 1);
        }
        mpz_poly_cleandeg(info->P, glob.n - 1);
        cachefile_close(c);
        return NULL;
    }

    ab_source ab;
    ab_source_init_set(ab, glob.ab);
    cxx_mpz_poly tmp;
    info->nab_loc = accumulate_ab_poly(P, ab, off0, off1, tmp);
    ab_source_clear(ab);

    if (wcache && cachefile_open_w(c)) {
        logprint("writing cache %s\n", c->basename);
        fprintf(c->f, "%zu\n", info->nab_loc);
        for(int i = 0 ; i < glob.n ; i++) {
            gmp_fprintf(c->f, "%Zx\n", mpz_poly_coeff_const(info->P, i));
        }
        cachefile_close(c);
    }

    return NULL;
}

void * a_poly_read_share_child2(struct subtask_info_t * info)
{
    mpz_poly_mul_mod_f(info->P0, info->P0, info->P1, glob.f_hat);
    return NULL;
}

size_t a_poly_read_share(mpz_poly_ptr P, size_t off0, size_t off1)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();

    size_t nab_loc = 0;
    int rc;
    log_begin();

    cachefile c;

    cachefile_init(c, "a_%zu_%zu", off0, off1);

    if (rcache && cachefile_open_r(c)) {
        logprint("reading cache %s\n", c->basename);
        rc = fscanf(c->f, "%zu", &nab_loc);
        ASSERT_ALWAYS(rc == 1);
        for(int i = 0 ; i < glob.n ; i++) {
            rc = gmp_fscanf(c->f, "%Zx", mpz_poly_coeff(P, i));
            ASSERT_ALWAYS(rc == 1);
        }
        mpz_poly_cleandeg(P, glob.n - 1);
        cachefile_close(c);
        return nab_loc;
    }

    size_t nparts = glob.ncores * glob.t;

    std::vector<subtask_info_t> a_tasks(glob.ncores);

    int j0 = glob.ncores * glob.arank;
    int j1 = j0 + glob.ncores;

    std::vector<cxx_mpz_poly> pols(j1-j0);
    for(int j = j0 ; j < j1 ; j++) {
        subtask_info_t & task = a_tasks[j-j0];
        task.P = pols[j-j0];
        task.off0 = off0 + (off1 - off0) * j / nparts;
        task.off1 = off0 + (off1 - off0) * (j+1) / nparts;
        task.nab_loc = 0;
        auto f = (wq_func_t) &a_poly_read_share_child;
        task.handle = wq_push(glob.wq, f, &task);
    }
    /* we're doing nothing, only waiting. */
    for(int j = j0 ; j < j1 ; j++) {
        wq_join(a_tasks[j-j0].handle);
        nab_loc += a_tasks[j-j0].nab_loc;
    }
    /* XXX freeing a_tasks is deferred because of course we still need A ! */

    STOPWATCH_GET();
    log_step_time(" done on leaf threads");
    log_step(": sharing among threads");

    std::vector<subtask_info_t> a_tasks2(glob.ncores);
    for(int done = 1 ; done < glob.ncores ; done<<=1) {
        for(int j = 0 ; j < glob.ncores ; j += done << 1) {
            if (j + done >= glob.ncores)
                break;
            subtask_info_t & task = a_tasks2[j];
            task.P0 = a_tasks[j].P;
            task.P1 = a_tasks[j+done].P;
            auto f = (wq_func_t) &a_poly_read_share_child2;
            task.handle = wq_push(glob.wq, f, &task);
        }
        /* we're doing nothing, only waiting. */
        for(int j = 0 ; j < glob.ncores ; j += done << 1) {
            if (j + done >= glob.ncores)
                break;
            wq_join(a_tasks2[j].handle);
        }
    }

    mpz_poly_swap(P, a_tasks[0].P);


    STOPWATCH_GET();
    log_step_time(" done locally");
    MPI_Barrier(MPI_COMM_WORLD);

    log_step(": sharing");

    allreduce_mpz_poly_mul_mod_f(P, glob.acomm, glob.f_hat);
    MPI_Allreduce(MPI_IN_PLACE, &nab_loc, 1, CADO_MPI_SIZE_T, MPI_SUM, MPI_COMM_WORLD);

    STOPWATCH_GET();
    log_end();

    if (wcache && cachefile_open_w(c)) {
        logprint("writing cache %s\n", c->basename);
        fprintf(c->f, "%zu\n", nab_loc);
        for(int i = 0 ; i < glob.n ; i++) {
            gmp_fprintf(c->f, "%Zx\n", mpz_poly_coeff_const(P, i));
        }
        cachefile_close(c);
    }

    return nab_loc;
}/*}}}*/

#if 0
struct tree_like_subtask_info_t {
    struct prime_data * p;
    int j;
    int i0, i1;
    struct wq_task * handle;

    /* In tree-like mode, tasks populate their argument with some info on
     * the child tasks they've spawned. It is of course mandatory to join
     * these tasks as well.
     */
    struct wq_task * t0;
    struct wq_task * t1;
};
#endif

/* {{{ precompute_powers */

void * precompute_powers_child(struct subtask_info_t * info)/* {{{ */
{
    struct prime_data * p = info->p;

    logprint("Precomputing p^%lu, p=%lu\n", glob.prec, p->p);
    // this triggers the whole precomputation.
    p->powers(glob.prec);

    return NULL;
}

/* }}} */

void precompute_powers(std::vector<prime_data> & primes, int i0, int i1)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();

    logprint("precompute_powers starts\n");

    {
        std::vector<subtask_info_t> tasks(i1-i0);
        for(int i = i0 ; i < i1 ; i++) {
            int k = i-i0;
            // if (k % glob.psize != glob.prank) continue;
            subtask_info_t & task = tasks[k];
            task.p = &primes[i];
            task.j = 0;
            auto f = (wq_func_t) &precompute_powers_child;
            task.handle = wq_push(glob.wq, f, &task);
        }
        /* we're doing nothing, only waiting. */
        for(int i = i0 ; i < i1 ; i++) {
            int k = i-i0;
            // if (k % glob.psize != glob.prank) continue;
            wq_join(tasks[k].handle);
        }
    }

    STOPWATCH_GET();

    logprint("precompute_powers ends. wct %.2lf.\n", w1-w0);
}
/*}}}*/

/* {{{ lifting roots */

void * lifting_roots_child(struct subtask_info_t * info)/* {{{ */
{
    struct prime_data * p = info->p;
    int j = info->j;

    mpz_srcptr p1 = p->powers(1);
    mpz_ptr rx = p->lroots[j];

    cachefile c;

    cachefile_init(c, "lroot_%lu_%lu_%lu", p->p, p->r[j], glob.prec);


    if (rcache && cachefile_open_r(c)) {
        logprint("reading cache %s\n", c->basename);
        gmp_fscanf(c->f, "%Zx", rx);
        cachefile_close(c);
        return NULL;
    }

    mpz_set_ui(rx, p->r[j]);

    cxx_mpz irx;

    mpz_poly_eval_mod_mpz(irx, glob.f_hat_diff, rx, p1);
    mpz_invert(irx, irx, p1);

    logprint("lifting p=%lu, r=%lu\n", p->p, p->r[j]);
    root_lift(p, rx, irx, glob.prec);

    if (wcache && cachefile_open_w(c)) {
        logprint("writing cache %s\n", c->basename);
        gmp_fprintf(c->f, "%Zx\n", rx);
        cachefile_close(c);
    }

    return NULL;
}

/* }}} */

void lifting_roots(std::vector<prime_data> & primes, int i0, int i1)
{
    int n = glob.n;

    STOPWATCH_DECL;
    STOPWATCH_GO();

    //  lifting all roots mod all primes (per pgroup)
    log_begin();

    {
        std::vector<subtask_info_t> lift_tasks((i1 - i0) * n);
        for(int j = 0 ; j < n ; j++) {
            for(int i = i0 ; i < i1 ; i++) {
                int k = (i-i0) * n + j;
                if (k % glob.psize != glob.prank)
                    continue;
                subtask_info_t & task = lift_tasks[k];
                task.p = &primes[i];
                task.j = j;
                auto f = (wq_func_t) &lifting_roots_child;
                task.handle = wq_push(glob.wq, f, &task);
            }
        }
        /* we're doing nothing, only waiting. */
        for(int j = 0 ; j < n ; j++) {
            for(int i = i0 ; i < i1 ; i++) {
                int k = (i-i0) * n + j;
                if (k % glob.psize != glob.prank)
                    continue;
                wq_join(lift_tasks[k].handle);
            }
        }
    }

    STOPWATCH_GET();

    log_step_time(" done locally");

    MPI_Barrier(MPI_COMM_WORLD);

    log_step(": sharing");

    {
        for(int j = 0 ; j < n ; j++) {
            for(int i = i0 ; i < i1 ; i++) {
                int k = (i-i0) * n + j;
                broadcast(primes[i].lroots[j], k % glob.psize, glob.pcomm);
            }
        }
    }
    STOPWATCH_GET();
    log_end();

    for(int i = i0 ; i < i1 ; i++)
        for(int j = 0 ; j < n ; j++)
            ASSERT_ALWAYS(mpz_size(primes[i].lroots[j]));
}/*}}}*/

/* {{{ rational reduction */

/* {{{ building rational product tree */
struct rat_ptree_s {
    struct rat_ptree_s * t0;
    struct rat_ptree_s * t1;
    mpz_t z;
    mpz_srcptr zx;
    struct prime_data * p;
};
using rat_ptree_t = struct rat_ptree_s;

//  rational product tree: all prime powers.
// it's quite easy to set up.
//
// As such the code does not use the worker threads at all.
// It is possible to do it in a distributed manner, however it is not
// clear that a lot is to be gained.

rat_ptree_t * rat_ptree_build_inner(struct prime_data * p, int i0, int i1)
{
    ASSERT_ALWAYS(i0 < i1);
    rat_ptree_t * res = (rat_ptree_t *) malloc(sizeof(rat_ptree_t));
    memset(res, 0, sizeof(rat_ptree_t));
    if (i1-i0 == 1) {
        res->p = p + i0;
        res->zx = res->p->powers(glob.prec);
        return res;
    }
    int d = (i1-i0)/2;
    res->t0 = rat_ptree_build_inner(p, i0, i0+d);
    res->t1 = rat_ptree_build_inner(p, i0+d, i1);
    mpz_init(res->z);
    res->zx = res->z;
    WRAP_mpz_mul(res->z, res->t0->zx, res->t1->zx);
    return res;
}

rat_ptree_t * rat_ptree_build(struct prime_data * p, int i0, int i1)
{
    if (i1 == i0)
        return NULL;

    STOPWATCH_DECL;
    STOPWATCH_GO();
    log_begin();
    rat_ptree_t * res = rat_ptree_build_inner(p, i0, i1);
    STOPWATCH_GET();
    log_end();
    return res;
}

void rat_ptree_clear(rat_ptree_t * t)
{
    if (t == NULL) return;
    rat_ptree_clear(t->t0);
    rat_ptree_clear(t->t1);
    if (t->p == NULL) mpz_clear(t->z);
    free(t);
}

#if 0
unsigned int rat_ptree_nleaves(rat_ptree_t * t)
{
    if (t->p)
        return 1;
    unsigned int n0 = rat_ptree_nleaves(t->t0);
    unsigned int n1 = rat_ptree_nleaves(t->t1);
    return n0 + n1;
}
#endif
/* }}} */

void reduce_poly_mod_rat_ptree(mpz_poly_ptr P, rat_ptree_t * T)/*{{{*/
{
    if (T == NULL)
        return;

    // of course this destroys P.
    // we assume that P is reduced mod the top-level T.
    if (T->p) {
        // on leaves, we have ->p filled. Not much to be done.
        mpz_poly_swap(T->p->A, P);
        return;
    }
    cxx_mpz_poly temp;
    temp->deg = glob.n - 1;
    for(int i = 0 ; i < glob.n ; i++) {
        WRAP_mpz_mod(mpz_poly_coeff(temp, i), mpz_poly_coeff_const(P, i), T->t0->zx);
        WRAP_mpz_mod(mpz_poly_coeff(P, i), mpz_poly_coeff_const(P, i), T->t1->zx);
    }
    mpz_poly_cleandeg(temp, glob.n - 1);
    mpz_poly_cleandeg(P, glob.n - 1);
    reduce_poly_mod_rat_ptree(temp, T->t0);
    reduce_poly_mod_rat_ptree(P, T->t1);
}/*}}}*/

void * rational_reduction_child(struct subtask_info_t * info)
{
    struct prime_data * primes = info->p;
    int i0 = info->i0;
    int i1 = info->i1;
    size_t off0 = info->off0;
    size_t off1 = info->off1;

    ASSERT_ALWAYS(info->P->deg >= 0);

    rat_ptree_t * ptree = rat_ptree_build(primes, i0, i1);

    if (ptree == NULL)
        ASSERT_ALWAYS(glob.r < glob.ncores);

    if (glob.r <= glob.ncores) {
        ASSERT_ALWAYS(ptree == NULL || ptree->p != NULL);
        if (ptree) {
            for(int i = 0 ; i < glob.n ; i++) {
                WRAP_mpz_mod(mpz_poly_coeff(ptree->p->A, i),
                        mpz_poly_coeff_const(info->P, i), ptree->zx);
            }
            mpz_poly_cleandeg(ptree->p->A, glob.n-1);
        }
        if (barrier_wait(glob.barrier, NULL, NULL, NULL) == BARRIER_SERIAL_THREAD) {
            cxx_mpz_poly foo;
            mpz_poly_swap(info->P, foo);
        }
    } else {
        ASSERT_ALWAYS(ptree != NULL);

        cxx_mpz_poly temp;
        for(int i = 0 ; i < glob.n ; i++) {
            WRAP_mpz_mod(mpz_poly_coeff(temp, i), mpz_poly_coeff_const(info->P, i), ptree->zx);
        }
        mpz_poly_cleandeg(temp, glob.n-1);
        if (barrier_wait(glob.barrier, NULL, NULL, NULL) == BARRIER_SERIAL_THREAD) {
            cxx_mpz_poly foo;
            mpz_poly_swap(info->P, foo);
        }
        // This computes P mod p_i for all p_i, and stores it into the
        // relevant field at the ptree leaves (->a)
        reduce_poly_mod_rat_ptree(temp, ptree);
    }

    rat_ptree_clear(ptree);

    for(int i = i0 ; i < i1 ; i++) {
        ASSERT_ALWAYS(primes[i].A->deg >= 0);
    }
    if (wcache) {
        for(int i = i0 ; i < i1 ; i++) {
            cachefile c;
            cachefile_init(c, "a_%zu_%zu_mod_%lu_%lu",
                    off0, off1, primes[i].p, glob.prec);
            if (cachefile_open_w(c)) {
                logprint("writing cache %s\n", c->basename);
                /* XXX nab_loc is a misnomer here. In reality we have nab_total
                 * there in this context */
                fprintf(c->f, "%zu\n", info->nab_loc);
                for(int j = 0 ; j < glob.n ; j++)
                    gmp_fprintf(c->f, "%Zx\n", mpz_poly_coeff_const(primes[i].A, j));
                cachefile_close(c);
            }
        }
    }

    return NULL;
}

void rational_reduction(std::vector<prime_data> & primes, int i0, int i1, mpz_poly_ptr P, size_t off0, size_t off1, size_t * p_nab_total)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    log_begin();
    int rc;

    ASSERT_ALWAYS(P->deg >= 0);

    if (rat_red_caches_ok(primes, i0, i1, off0, off1)) {
        ASSERT_ALWAYS(*p_nab_total == 0);
        for(int i = i0 ; i < i1 ; i++) {
            cachefile c;
            cachefile_init(c, "a_%zu_%zu_mod_%lu_%lu",
                    off0, off1, primes[i].p, glob.prec);
            if (cachefile_open_r(c)) {
                logprint("reading cache %s\n", c->basename);
                rc = fscanf(c->f, "%zu", p_nab_total);
                ASSERT_ALWAYS(rc == 1);
                for(int j = 0 ; j < glob.n ; j++) {
                    rc = gmp_fscanf(c->f, "%Zx", mpz_poly_coeff(primes[i].A, j));
                    ASSERT_ALWAYS(rc == 1);
                }
                cachefile_close(c);
            }
        }
        return;
    }

    {
        std::vector<subtask_info_t> tasks(glob.ncores);
        for(int k = 0 ; k < glob.ncores ; k++) {
            subtask_info_t & task = tasks[k];
            task.p = primes.data();
            task.i0 = i0 + (i1-i0) * k / glob.ncores;
            task.i1 = i0 + (i1-i0) * (k+1) / glob.ncores;
            task.off0 = off0;
            task.off1 = off1;
            task.P = P;
            task.nab_loc = *p_nab_total;
            auto f = (wq_func_t) &rational_reduction_child;
            task.handle = wq_push(glob.wq, f, &task);
        }
        /* we're doing nothing, only waiting. */
        for(int k = 0 ; k < glob.ncores ; k++) wq_join(tasks[k].handle);
    }

    STOPWATCH_GET();
    log_end();

    for(int i = i0 ; i < i1 ; i++) {
        ASSERT_ALWAYS(primes[i].A->deg >= 0);
    }
}
/* }}} */

/* {{{ algebraic reduction */

void * algebraic_reduction_child(struct subtask_info_t * info)
{
    struct prime_data * p = info->p;
    int j = info->j;
    int rc;
    logprint("alg_red (%lu, x-%lu) starts\n", p->p, p->r[j]);

    cachefile c;

    cachefile_init(c, "a_%zu_%zu_mod_%lu_%lu_%lu",
            info->off0, info->off1, info->p->p, info->p->r[j], glob.prec);

    if (rcache && cachefile_open_r(c)) {
        logprint("reading cache %s\n", c->basename);
        rc = fscanf(c->f, "%zu", &info->nab_loc);
        ASSERT_ALWAYS(rc == 1);
        rc = gmp_fscanf(c->f, "%Zx", (mpz_ptr) p->evals[j]);
        ASSERT_ALWAYS(rc == 1);
        cachefile_close(c);
        return NULL;
    }

    cxx_mpz ta;
    mpz_srcptr px = p->powers(glob.prec);
    mpz_set(ta, mpz_poly_coeff_const(p->A, glob.n - 1));
    for(int k = glob.n - 2 ; k >= 0 ; k--) {
        WRAP_mpz_mul(ta, ta, p->lroots[j]);
        mpz_add(ta, ta, mpz_poly_coeff_const(p->A, k));
        WRAP_mpz_mod(ta, ta, px);
    }
    WRAP_mpz_mul(p->evals[j], p->evals[j], ta);
    if (mpz_size(p->evals[j]) >= mpz_size(px) * 3/2)
        WRAP_mpz_mod(p->evals[j], p->evals[j], px);

    if (wcache && cachefile_open_w(c)) {
        logprint("writing cache %s\n", c->basename);
        fprintf(c->f, "%zu\n", info->nab_loc);
        gmp_fprintf(c->f, "%Zx\n", (mpz_srcptr) p->evals[j]);
        cachefile_close(c);
    }
    return NULL;
}

void algebraic_reduction(std::vector<prime_data> & primes, int i0, int i1, size_t off0, size_t off1, size_t * p_nab_total)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    log_begin();

    if (!alg_red_caches_ok(primes, i0, i1, off0, off1)) {
        for(int i = i0 ; i < i1 ; i++) {
            ASSERT_ALWAYS(primes[i].A->deg >= 0);
        }
    }

    std::vector<subtask_info_t> tasks((i1-i0) * glob.n);
    for(int j = 0 ; j < glob.n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * glob.n + j;
            // if (k % glob.psize != glob.prank) continue;
            subtask_info_t & task = tasks[k];
            task.off0 = off0;
            task.off1 = off1;
            task.p = &primes[i];
            task.j = j;
            task.nab_loc = *p_nab_total;
            auto f = (wq_func_t) &algebraic_reduction_child;
            task.handle = wq_push(glob.wq, f, &task);
        }
    }
    /* we're doing nothing */
    for(int j = 0 ; j < glob.n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * glob.n + j;
            // if (k % glob.psize != glob.prank) continue;
            wq_join(tasks[k].handle);
            * p_nab_total = tasks[k].nab_loc;
        }
    }
    for(int i = i0 ; i < i1 ; i++) {
        for(int k = glob.n - 1 ; k >= 0 ; k--) {
            mpz_realloc(mpz_poly_coeff(primes[i].A, k), 0);
        }
    }

    STOPWATCH_GET();
    log_end();
}/*}}}*/

void multiply_all_shares(std::vector<prime_data> & primes, int i0, int i1, size_t * p_nab_total)/*{{{*/
{
    STOPWATCH_DECL;
    STOPWATCH_GO();
    log_begin();
    int rc;

    // ok. for all eval points, collect the reduced stuff. we'll do the
    // collection in a tree-like manner.
    // this all happens within pcomm, which has size s.

    if (a_shared_caches_ok(primes, i0, i1)) {
        ASSERT_ALWAYS(*p_nab_total == 0);
        for(int i = i0 ; i < i1 ; i++) {
            struct prime_data * p = &primes[i];
            for(int j = 0 ; j < glob.n ; j++) {
                cachefile c;
                cachefile_init(c, "a_mod_%lu_%lu_%lu", p->p, p->r[j], glob.prec);
                int ok = cachefile_open_r(c);
                ASSERT_ALWAYS(ok);
                logprint("reading cache %s\n", c->basename);
                rc = fscanf(c->f, "%zu", p_nab_total);
                ASSERT_ALWAYS(rc == 1);
                rc = gmp_fscanf(c->f, "%Zx", (mpz_ptr) primes[i].evals[j]);
                ASSERT_ALWAYS(rc == 1);
                cachefile_close(c);
            }
        }
        return ;
    }

    // stupid -- this must be multithreaded
    // XXX
    for(int i = i0 ; i < i1 ; i++) {
        mpz_srcptr px = primes[i].powers(glob.prec);
        for(int j = 0 ; j < glob.n ; j++) {
            allreduce_mulmod_mpz(primes[i].evals[j], glob.pcomm, px);
        }
    }

    STOPWATCH_GET();
    log_end();

    if (wcache) {
        logprint("writing caches %d..%d\n", i0, i1);
        for(int i = i0 ; i < i1 ; i++) {
            struct prime_data * p = &primes[i];
            for(int j = 0 ; j < glob.n ; j++) {
                cachefile c;
                cachefile_init(c, "a_mod_%lu_%lu_%lu", p->p, p->r[j], glob.prec);
                if (cachefile_open_w(c)) {
                    /* it's ok if we can't write the cache file on every
                     * rank, since several ranks can very well share the same
                     * storage.
                     */
                    logprint("writing cache %s\n", c->basename);
                    fprintf(c->f, "%zu\n", * p_nab_total);
                    gmp_fprintf(c->f, "%Zx\n", (mpz_srcptr) primes[i].evals[j]);
                    cachefile_close(c);
                }
            }
        }
    }
}/*}}}*/

/* {{{ local square roots*/

void * local_square_roots_child(struct subtask_info_t * info)
{
    struct prime_data * p = info->p;
    int j = info->j;
    int rc;
    cachefile c;
    cachefile_init(c, "sqrt_%lu_%lu_%lu", p->p, p->r[j], glob.prec);
    if (rcache && cachefile_open_r(c)) {
        rc = fscanf(c->f, "%zu", &info->nab_loc);
        ASSERT_ALWAYS(rc == 1);
        logprint("reading cache %s\n", c->basename);
        rc = gmp_fscanf(c->f, "%Zx", (mpz_ptr) p->sqrts[j]);
        ASSERT_ALWAYS(rc == 1);
        cachefile_close(c);
        return NULL;
    }
    printf("# [%2.2lf] [P%dA%d] lifting sqrt (%lu, x-%lu) (p->evals[%d] has size %zu)\n", WCT, glob.arank, glob.prank, p->p, p->r[j], j, mpz_size(p->evals[j]));
    sqrt_lift(p, p->evals[j], p->sqrts[j], glob.prec);
    // printf("# [%2.2lf] done\n", WCT);
    mpz_realloc(p->evals[j], 0);

    if (wcache && cachefile_open_w(c)) {
        logprint("writing cache %s\n", c->basename);
        fprintf(c->f, "%zu\n", info->nab_loc);
        gmp_fprintf(c->f, "%Zx\n", (mpz_srcptr) p->sqrts[j]);
        cachefile_close(c);
    }

    return NULL;
}

void local_square_roots(std::vector<prime_data> & primes, int i0, int i1, size_t * p_nab_total)
{
    int n = glob.n;

    STOPWATCH_DECL;
    STOPWATCH_GO();

    log_begin();

    std::vector<subtask_info_t> tasks((i1-i0)*glob.n);
    for(int j = 0 ; j < glob.n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * glob.n + j;
            if (k % glob.psize != glob.prank) continue;
            auto & task = tasks[k];
            task.p = &primes[i];
            task.j = j;
            task.nab_loc = * p_nab_total;
            auto f = (wq_func_t) &local_square_roots_child;
            task.handle = wq_push(glob.wq, f, &task);
        }
    }
    /* we're doing nothing */
    for(int j = 0 ; j < glob.n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * glob.n + j;
            if (k % glob.psize != glob.prank) continue;
            wq_join(tasks[k].handle);
            * p_nab_total = tasks[k].nab_loc;
        }
    }

    STOPWATCH_GET();
    log_step_time(" done locally");

    MPI_Barrier(MPI_COMM_WORLD);
    log_step(": sharing");

    for(int j = 0 ; j < n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * n + j;
            mpz_ptr z = primes[i].sqrts[j];
            broadcast(z, k % glob.psize, glob.pcomm);
        }
    }

    STOPWATCH_GET();
    log_end();
}

/* }}} */

/* {{{ prime inversion lifts */

cxx_mpz inversion_lift(struct prime_data * p, cxx_mpz const & Hx, int precision)/* {{{ */
{
    double w0 = WCT;

    auto const & pk = p->powers(precision);
    ASSERT(precision > 0);

    if (precision == 1) {
        cxx_mpz iHx;
        mpz_invert(iHx, Hx, pk);
        return iHx;
    }
    int lower = precision - precision / 2;

    auto const & pl = p->powers(lower);

    // we're going to recurse. change the long Hx_mod value stored by a
    // temporary small one. The problem with this approach is that we
    // store many reductions in memory.
    cxx_mpz Hx_save;
    WRAP_mpz_mod(Hx_save, Hx, pl);
    // recurse.
    auto iHx = inversion_lift(p, Hx_save, lower);

    if (WCT > w0 + print_delay)
        logprint("precision %d\n", precision);

    cxx_mpz ta;

    WRAP_mpz_mul(ta, iHx, Hx);
    WRAP_mpz_mod(ta, ta, pk);

    mpz_sub_ui(ta, ta, 1);
    WRAP_mpz_mul(ta, ta, iHx);
    mpz_sub(iHx, iHx, ta);
    WRAP_mpz_mod(iHx, iHx, pk);

    return iHx;
    // gmp_printf("# [%2.2lf] %Zd\n", WCT, p->iHx_mod);
}/* }}} */

void * prime_inversion_lifts_child(struct subtask_info_t * info)
{
    struct prime_data * p = info->p;
    mpz_srcptr px = p->powers(glob.prec);

    cxx_mpz Hx;
    mpz_divexact_ui(Hx, glob.P, p->p);
    mpz_powm_ui(Hx, Hx, glob.prec, px);
    // compute the inverse of H^prec modulo p^prec.
    // need a recursive function for computing the inverse.
    printf("# [%2.2lf] [P%dA%d] lifting H^-l\n", WCT, glob.arank, glob.prank);
    p->iHx = inversion_lift(p, Hx, glob.prec);
    return NULL;
}

void prime_inversion_lifts(std::vector<prime_data> & primes, int i0, int i1)
{
    STOPWATCH_DECL;
    STOPWATCH_GO();

    log_begin();

    {
        std::vector<subtask_info_t> tasks(i1-i0);
        for(int i = i0 ; i < i1 ; i++) {
            int k = i-i0;
            if (k % glob.psize != glob.prank)
                continue;
            auto & task = tasks[k];
            task.p = &primes[i];
            task.j = INT_MAX;
            auto f = (wq_func_t) &prime_inversion_lifts_child;
            task.handle = wq_push(glob.wq, f, &task);
        }
        /* we're doing nothing, only waiting. */
        for(int i = i0 ; i < i1 ; i++) {
            int k = i-i0;
            if (k % glob.psize != glob.prank)
                continue;
            wq_join(tasks[k].handle);
        }
    }

    STOPWATCH_GET();
    log_step_time(" done locally");

    MPI_Barrier(MPI_COMM_WORLD);

    log_step(": sharing");

    for(int i = i0 ; i < i1 ; i++) {
        int k = i-i0;
        broadcast(primes[i].iHx, k % glob.psize, glob.pcomm);
    }

    STOPWATCH_GET();
    log_end();
}
/* }}} */

/* {{{ prime postcomputations (lagrange reconstruction) */
struct postcomp_subtask_info_t {
    struct prime_data * p;
    int j;
    int64_t * c64;
    mp_limb_t * cN;

    struct wq_task * handle;
};

void * prime_postcomputations_child(struct postcomp_subtask_info_t * info)
{
    struct prime_data * p = info->p;
    int j = info->j;
    int64_t * c64 = info->c64;
    mp_limb_t * cN = info->cN;

    mpz_srcptr px = p->powers(glob.prec);
    // printf("# [%2.2lf] done\n", WCT);

    // Lagrange reconstruction.
    //
    // Normally each coefficient has to be divided by the evaluation of
    // the derivative. However we skip this division, effectively
    // reconstructing the polynomial multiplied by the square of the
    // derivative -- which is exactly what we're looking for, in fact.

    // recall that we're working with the number field sieve in mind. So
    // we don't really care about the whole reconstruction in the number
    // field, and from here on we are going to take wild shortcuts.
    // Indeed, even though the ``magical sign combination'' is not known
    // at this point, we do know that the eventual reconstruction will be
    // linear. Thus instead of storing n^2 full length modular integers
    // (n for each root),  and do this for each prime, we store only the
    // pair (quotient mod p^x, residue mod N).

    mpf_t pxf, ratio;
    cxx_mpz z;

    mpf_init2(pxf, 256);
    mpf_init2(ratio, 256);

    mpf_set_z(pxf, px);

    cxx_mpz Hxm;
    mpz_set_ui(Hxm, p->p);
    mpz_invert(Hxm, Hxm, glob.cpoly->n);
    mpz_mul(Hxm, Hxm, glob.P);
    mpz_mod(Hxm, Hxm, glob.cpoly->n);
    mpz_powm_ui(Hxm, Hxm, glob.prec, glob.cpoly->n);

    cxx_mpz ta, tb;

    /* work mod one root */
    mpz_srcptr rx = p->lroots[j];
    mpz_ptr sx = p->sqrts[j];

    ASSERT_ALWAYS(mpz_size(rx));
    ASSERT_ALWAYS(mpz_size(sx));

    // so we have this nice square root. The first thing we do on our
    // list is to scramble it by multiplying it with the inverse of
    // H^x...
    WRAP_mpz_mul(sx, sx, p->iHx);
    WRAP_mpz_mod(sx, sx, px);

    // Now use the evaluation of f_hat mod rx to obtain the lagrange
    // coefficients.
    mpz_set_ui(ta, 1);
    for(int k = glob.n - 1 ; k >= 0 ; k--) {
        if (k < glob.n - 1) {
            WRAP_mpz_mul(ta, ta, rx);
            mpz_add(ta, ta, mpz_poly_coeff_const(glob.f_hat, k+1));
            WRAP_mpz_mod(ta, ta, px);
        }
        // multiply directly with H^-x * sqrt
        WRAP_mpz_mul(tb, ta, sx);
        if (k < glob.n - 1) {
            WRAP_mpz_mod(tb, tb, px);
        }
        ASSERT_ALWAYS(mpz_cmp_ui(tb, 0) >= 0);
        ASSERT_ALWAYS(mpz_cmp(tb, px) < 0);

        // now the shortcuts.
        mpf_set_z(ratio, tb);
        mpf_div(ratio, ratio, pxf);
        mpf_mul_2exp(ratio, ratio, 64);
        mpz_set_f(z, ratio);

        uint64_t u;
#if GMP_LIMB_BITS == 64
        u = mpz_get_ui(z);
#else
        u = (uint64_t) mpz_getlimbn(z,1);
        u <<= 32;
        u |= (uint64_t) mpz_getlimbn(z,0);
#endif
        c64[k] = (int64_t) u;

        mpz_mul(tb, tb, Hxm);
        mpz_mod(tb, tb, glob.cpoly->n);
        mp_size_t sN = mpz_size(glob.cpoly->n);
        ASSERT_ALWAYS(mpz_sgn(tb) > 0);
        MPN_SET_MPZ(cN + k * sN, sN, tb);
    }

    mpf_clear(pxf);
    mpf_clear(ratio);
    return NULL;
}

void prime_postcomputations(std::vector<prime_data> & primes, int i0, int i1, int64_t * contribs64, mp_limb_t * contribsN)
{
    int n = glob.n;

    STOPWATCH_DECL;
    STOPWATCH_GO();

    log_begin();

    std::vector<postcomp_subtask_info_t> tasks((i1-i0)*n);;
    mp_size_t sN = mpz_size(glob.cpoly->n);
    for(int j = 0 ; j < n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * n + j;
            if (k % glob.psize != glob.prank) continue;
            auto & task = tasks[k];
            task.p = &primes[i];
            task.j = j;
            task.c64 = contribs64 + (i * n + j) * n;
            task.cN = contribsN + ((i * n + j) * n) * sN;
            auto f = (wq_func_t) &prime_postcomputations_child;
            task.handle = wq_push(glob.wq, f, &task);
        }
    }
    /* we're doing nothing */
    for(int j = 0 ; j < n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * n + j;
            if (k % glob.psize != glob.prank) continue;
            wq_join(tasks[k].handle);
        }
    }

    STOPWATCH_GET();
    log_step_time(" done locally");

    MPI_Barrier(MPI_COMM_WORLD);

    log_step(": sharing");

    for(int j = 0 ; j < n ; j++) {
        for(int i = i0 ; i < i1 ; i++) {
            int k = (i-i0) * n + j;
            int64_t * c64 = contribs64 + (i * n + j) * n;
            mp_limb_t * cN = contribsN + ((i * n + j) * n) * sN;
            MPI_Bcast(c64, n, CADO_MPI_INT64_T,  k % glob.psize, glob.pcomm);
            MPI_Bcast(cN, n*sN, CADO_MPI_MP_LIMB_T, k % glob.psize, glob.pcomm);
        }
    }

    STOPWATCH_GET();
    log_end();
}
void old_prime_postcomputations(int64_t * c64, mp_limb_t * cN, struct prime_data * p)/* {{{ */
{
    mpz_srcptr px = p->powers(glob.prec);

    // Lagrange reconstruction.
    //
    // Normally each coefficient has to be divided by the evaluation of
    // the derivative. However we skip this division, effectively
    // reconstructing the polynomial multiplied by the square of the
    // derivative -- which is exactly what we're looking for, in fact.

    // recall that we're working with the number field sieve in mind. So
    // we don't really care about the whole reconstruction in the number
    // field, and from here on we are going to take wild shortcuts.
    // Indeed, even though the ``magical sign combination'' is not known
    // at this point, we do know that the eventual reconstruction will be
    // linear. Thus instead of storing n^2 full length modular integers
    // (n for each root),  and do this for each prime, we store only the
    // pair (quotient mod p^x, residue mod N).

    mpf_t pxf, ratio;
    cxx_mpz z;

    mpf_init2(pxf, 256);
    mpf_init2(ratio, 256);

    mpf_set_z(pxf, px);

    cxx_mpz Hxm;
    mpz_set_ui(Hxm, p->p);
    mpz_invert(Hxm, Hxm, glob.cpoly->n);
    mpz_mul(Hxm, Hxm, glob.P);
    mpz_mod(Hxm, Hxm, glob.cpoly->n);
    mpz_powm_ui(Hxm, Hxm, glob.prec, glob.cpoly->n);

    cxx_mpz ta, tb;

    // XXX Eh ! mpi-me !
    for(int j = 0 ; j < glob.n ; j++) {
        mpz_srcptr rx = p->lroots[j];
        mpz_ptr sx = p->sqrts[j];

        // so we have this nice square root. The first thing we do on our
        // list is to scramble it by multiplying it with the inverse of
        // H^x...
        mpz_mul(sx, sx, p->iHx);
        mpz_mod(sx, sx, px);

        // Now use the evaluation of f_hat mod rx to obtain the lagrange
        // coefficients.
        mpz_set_ui(ta, 1);
        for(int k = glob.n - 1 ; k >= 0 ; k--) {
            if (k < glob.n - 1) {
                WRAP_mpz_mul(ta, ta, rx);
                mpz_add(ta, ta, mpz_poly_coeff_const(glob.f_hat, k+1));
                WRAP_mpz_mod(ta, ta, px);
            }
            // multiply directly with H^-x * sqrt
            WRAP_mpz_mul(tb, ta, sx);
            if (k < glob.n - 1) {
                WRAP_mpz_mod(tb, tb, px);
            }
            ASSERT_ALWAYS(mpz_cmp_ui(tb, 0) >= 0);
            ASSERT_ALWAYS(mpz_cmp(tb, px) < 0);

            // now the shortcuts.
            mpf_set_z(ratio, tb);
            mpf_div(ratio, ratio, pxf);
            mpf_mul_2exp(ratio, ratio, 64);
            mpz_set_f(z, ratio);

            uint64_t u;
#if GMP_LIMB_BITS == 64
            u = mpz_get_ui(z);
#else
            u = (uint64_t) mpz_getlimbn(z,1);
            u <<= 32;
            u |= (uint64_t) mpz_getlimbn(z,0);
#endif
            c64[j*glob.n+k] = (int64_t) u;

            mpz_mul(tb, tb, Hxm);
            mpz_mod(tb, tb, glob.cpoly->n);
            mp_size_t sN = mpz_size(glob.cpoly->n);
            ASSERT_ALWAYS(mpz_sgn(tb) > 0);
            MPN_SET_MPZ(cN + (j*glob.n+k) * sN, sN, tb);
        }
    }
    mpf_clear(pxf);
    mpf_clear(ratio);
}
/* }}} */

/* }}} */

void mpi_set_communicators()/*{{{*/
{
    int s = glob.s;
    int t = glob.t;
    if (s * t > glob.nprocs) {
        if (glob.rank == 0) {
            fprintf(stderr, "Error: need at least %d*%d=%d jobs (or more RAM)\n",
                    s,t,s*t);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    } else if (s * t < glob.nprocs) {
        if (glob.rank == 0) {
            fprintf(stderr, "Warning: only %d*%d=%d jobs needed, have %d\n",
                    s,t,s*t,glob.nprocs);
        }
    }

    MPI_Comm_split(MPI_COMM_WORLD, glob.rank / t, glob.rank, &glob.acomm);
    MPI_Comm_rank(glob.acomm, &glob.arank);
    MPI_Comm_size(glob.acomm, &glob.asize);
    ASSERT_ALWAYS(glob.asize == t);

    MPI_Comm_split(MPI_COMM_WORLD, glob.rank % t, glob.rank, &glob.pcomm);
    MPI_Comm_rank(glob.pcomm, &glob.prank);
    MPI_Comm_size(glob.pcomm, &glob.psize);
    ASSERT_ALWAYS(glob.psize == s);
}
/*}}}*/

static void banner()
{
    MPI_Barrier(MPI_COMM_WORLD);
    if (glob.rank == 0) {
        printf("###########################################################################\n");
    }
}

int main(int argc, char const ** argv)
{
    int ret, i;
    int asked_r = 0;
    // int size_guess = 0;

    MPI_Init(&argc, (char ***) &argv);
    program_starttime = wct_seconds();

    MPI_Comm_rank(MPI_COMM_WORLD, &glob.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &glob.nprocs);

    /* {{{ parameter parsing */
    /* print the command line */
    printf("%s.r%s", argv[0], cado_revision_string);
    for (i = 1; i < argc; i++)
        printf(" %s", argv[i]);
    printf("\n");

    cxx_param_list pl;
    int cache=0;
    // param_list_configure_switch(pl, "-v", &verbose);
    param_list_configure_switch(pl, "--cache", &cache);
    param_list_configure_switch(pl, "--rcache", &rcache);
    param_list_configure_switch(pl, "--wcache", &wcache);
    // param_list_configure_switch(pl, "--size-guess", &size_guess);
    int wild = 0;
    argv++, argc--;
    for (; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }
        if (argv[0][0] != '-' && wild == 0) {
            param_list_add_key(pl, "depfile", argv[0],
                    PARAMETER_FROM_CMDLINE);
            wild++;
            argv++, argc--;
            continue;
        }
        if (argv[0][0] != '-' && wild == 1) {
            param_list_add_key(pl, "ratdepfile", argv[0],
                    PARAMETER_FROM_CMDLINE);
            wild++;
            argv++, argc--;
            continue;
        }
        if (argv[0][0] != '-' && wild == 2) {
            param_list_add_key(pl, "polyfile", argv[0],
                    PARAMETER_FROM_CMDLINE);
            wild++;
            argv++, argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage();
    }

    if (cache)  rcache = wcache = 1;

    param_list_parse_double(pl, "print_delay", &print_delay);
    param_list_parse_double(pl, "ram", &ram_gb);

    param_list_parse_int(pl, "verbose", &verbose);
    param_list_configure_switch(pl, "--cache", &cache);
    param_list_parse_int(pl, "ncores", &glob.ncores);
    param_list_parse_int(pl, "r", &asked_r);
    param_list_parse_int(pl, "lll_maxdim", &glob.lll_maxdim);

    if (param_list_lookup_string(pl, "depfile") == NULL)
        usage();
    if (param_list_lookup_string(pl, "ratdepfile") == NULL)
        usage();
    if (param_list_lookup_string(pl, "polyfile") == NULL)
        usage();
    /* }}} */

    ret = cado_poly_read(glob.cpoly, param_list_lookup_string(pl, "polyfile"));
    /* This assumes that we have a rational side 0 and an algebraic side 1*/
    ASSERT_ALWAYS(cado_poly_get_ratside(glob.cpoly) == 0);
    glob.n = glob.cpoly->pols[1]->deg;
    ASSERT_ALWAYS(ret == 1);
    mpz_init(glob.root_m);
    cado_poly_getm(glob.root_m, glob.cpoly, glob.cpoly->n);

    /* {{{ create f_hat, the minimal polynomial of alpha_hat = lc(f) *
     * alpha */
    {
        mpz_poly_to_monic(glob.f_hat, glob.cpoly->pols[1]);
        mpz_poly_derivative(glob.f_hat_diff, glob.f_hat);

        if (glob.rank == 0)
            printf("# [%2.2lf] Note: all computations done with polynomial f_hat\n", WCT);
    }
    /* }}} */

    // printf("# [%2.2lf] A is f_d^%zu*f_hat'(alpha_hat)*prod(f_d a - b alpha_hat)\n", WCT, nab + (nab &1));

    ab_source_init(glob.ab, param_list_lookup_string(pl, "depfile"),
            glob.rank, 0, MPI_COMM_WORLD);

    // note that for rsa768, this estimation takes only 10 minutes, so
    // it's not a big trouble.
    unsigned long toto;
    if (!param_list_parse_ulong(pl, "sqrt_coeffs_bits", &toto)) {
        if (glob.rank == 0) {
            estimate_nbits_sqrt(&glob.nbits_sqrt, glob.ab); //, size_guess);
        }
        MPI_Bcast(&glob.nbits_sqrt, 1, CADO_MPI_SIZE_T, 0, MPI_COMM_WORLD);
    } else {
        glob.nbits_sqrt = toto;
    }
    // MPI_Bcast(&glob.nbits_a, 1, CADO_MPI_SIZE_T, 0, MPI_COMM_WORLD);
    // we no longer need to know nab, so let's drop it as a proof !
    glob.ab->nab = 0;

    glob.r=0;
    glob.s=0;
    glob.t=0;
    if (glob.rank == 0) {
        get_parameters(&glob.r,&glob.s,&glob.t,asked_r);
    }
    broadcast(glob.s, 0, MPI_COMM_WORLD);
    broadcast(glob.t, 0, MPI_COMM_WORLD);
    broadcast(glob.r, 0, MPI_COMM_WORLD);
    ASSERT_ALWAYS(glob.r >= 1);
    ASSERT_ALWAYS(glob.s >= 1);
    ASSERT_ALWAYS(glob.t >= 1);

    int r = glob.r;
    int s = glob.s;
    // int n = glob.n;
    int t = glob.t;
    glob.m = r * t;


    // this is sufficiently trivial to do.
    auto primes = suitable_crt_primes();

    glob.P = 1;
    double log2_P = 0;
    for(int i = 0 ; i < glob.m ; i++) {
        log2_P += log(primes[i].p)/M_LN2;
        mpz_mul_ui(glob.P, glob.P, primes[i].p);
    }
    glob.prec = ceil((glob.nbits_sqrt + 128)/ log2_P);

    ASSERT_ALWAYS(mpi_data_agrees(glob.prec, MPI_COMM_WORLD));

    if (glob.rank == 0) {
        char sbuf[32];
        printf("# [%2.2lf] Lifting to precision l=%lu (p^l is approx %s)\n", WCT, glob.prec, size_disp(glob.prec * log(primes[0].p)/M_LN2 / 8, sbuf));
    }

    mpi_set_communicators();

    if (glob.rank == 0) {
        printf("# [%2.2lf] starting %d worker threads on each node\n", WCT, glob.ncores);
    }
    wq_init(glob.wq, glob.ncores);
    barrier_init(glob.barrier, NULL, glob.ncores);

    int pgnum = glob.rank % t;
    int apnum = glob.rank / t;

    snprintf(prefix, sizeof(prefix), "[P%dA%d] ", glob.arank, glob.prank);

    ASSERT_ALWAYS(glob.arank == pgnum);
    ASSERT_ALWAYS(glob.prank == apnum);

    int i0 = pgnum*r;
    int i1 = i0 + r;

    size_t off0 = apnum * glob.ab->totalsize / s;
    size_t off1 = (apnum+1) * glob.ab->totalsize / s;

    for(int i = i0 ; i < i1 ; i++)
        for(int j = 0 ; j < glob.n ; j++)
            primes[i].evals[j] = 1;

    // unless we have really nothing to do, we will need these two steps:

    banner(); /*********************************************/
    precompute_powers(primes, i0, i1);

    banner(); /*********************************************/
    lifting_roots(primes, i0, i1);

    // find some info on the caches.
    logprint("Checking sqrt caches\n");
    int sqrt_caches = sqrt_caches_ok(primes, i0, i1);
    logprint("Checking a caches\n");
    int shared_caches = a_shared_caches_ok(primes, i0, i1);
    logprint("Checking alg_red caches\n");
    int alg_caches = alg_red_caches_ok(primes, i0, i1, off0, off1);
    logprint("Checking rat_red caches\n");
    int rat_caches = rat_red_caches_ok(primes, i0, i1, off0, off1);

    if (rcache) {
        logprint("Caches: sqrt %s, shared %s, alg %s, rat %s\n",
                sqrt_caches ? "ok" : "no",
                shared_caches ? "ok" : "no",
                alg_caches ? "ok" : "no",
                rat_caches ? "ok" : "no");
    }

    size_t nab = 0;


    int next = 0;
    if (sqrt_caches) {
        next = 4;
    } else if (shared_caches) {
        next = 3;
    } else if (alg_caches) {
        next = 2;
    } else if (rat_caches) {
        next = 1;
    } else {
        next = 0;
    }

    logprint("Next step is: %d\n", next);

    {
        cxx_mpz_poly P;
        if (next < 1) {
            banner(); /*********************************************/
            nab = a_poly_read_share(P, off0, off1);
        }

        if (next < 2) {
            banner(); /*********************************************/
            rational_reduction(primes, i0, i1, P, off0, off1, &nab);
        }
        /* we no longer need P, we can free it now */
    }

    if (next < 3) {
        banner(); /*********************************************/
        algebraic_reduction(primes, i0, i1, off0, off1, &nab);
    }

    if (next < 4) {
        banner(); /*********************************************/
        multiply_all_shares(primes, i0, i1, &nab);
        // XXX This is a hack. we are evaluating the products
        // f_d^\epsilon*A*f_hat'(\hat\alpha)
        //
        // where A is the value denoted by the variable named A. It
        // is defined as:
        //
        // A(\hat\alpha) = (\prod_{(a,b)}(f_da-b\hat\alpha)
        //
        // Instead of fixing what's missing in A, we use shortcuts.
        // Extra f_d coefficients are added to the evaluation array,
        // and the derivative is eliminated later on by avoiding the
        // normalization in the Lagrange step.
        if (nab & 1) {
            logprint("Correcting by one leading term\n");
            for(int i = i0 ; i < i1 ; i++) {
                for(int j =  0 ; j < glob.n ; j++) {
                    mpz_ptr x = primes[i].evals[j];
                    mpz_mul(x, x, mpz_poly_coeff_const(glob.cpoly->pols[1], glob.n));
                }
            }
        }
    }

    banner(); /*********************************************/
    local_square_roots(primes, i0, i1, &nab);

    logprint("Number of pairs is %zu\n", nab);

    banner(); /*********************************************/
    prime_inversion_lifts(primes, i0, i1);

    banner(); /*********************************************/
    size_t nc = glob.m * glob.n * glob.n;
    std::vector<int64_t> contribs64(nc, 0);
    mp_size_t sN = mpz_size(glob.cpoly->n);
    std::vector<mp_limb_t> contribsN(nc * sN, 0);
#if 1
    prime_postcomputations(primes, i0, i1, contribs64.data(), contribsN.data());
    // now share the contribs !
    for(int i = 0 ; i < glob.m ; i++) {
        int root = i / r;
        int d64 = glob.n * glob.n;
        int dN = d64 * mpz_size(glob.cpoly->n);
        MPI_Bcast(contribs64.data() + i * d64, d64, CADO_MPI_INT64_T, root, glob.acomm);
        MPI_Bcast(contribsN.data() + i * dN, dN, CADO_MPI_MP_LIMB_T, root, glob.acomm);
    }
#else
    if (glob.prank == 0) {
        for(int i = i0 ; i < i1 ; i++) {
            int disp = i * glob.n * glob.n;
            old_prime_postcomputations(contribs64 + disp, contribsN + disp * sN, &primes[i]);
        }
    }

    // now share the contribs !
    if (glob.prank == 0) {
        for(int i = 0 ; i < glob.m ; i++) {
            int root = i / r;
            int d64 = glob.n * glob.n;
            int dN = d64 * mpz_size(glob.cpoly->n);
            MPI_Bcast(contribs64 + i * d64, d64, CADO_MPI_INT64_T, root, glob.acomm);
            MPI_Bcast(contribsN + i * dN, dN, CADO_MPI_MP_LIMB_T, root, glob.acomm);
        }
    }
#endif

    if (glob.rank == 0) {
        printf("# [%2.2lf] clearing work queues\n", WCT);
    }
    wq_clear(glob.wq, glob.ncores);
    barrier_destroy(glob.barrier, NULL);

    banner(); /*********************************************/
    // and only the leader does the knapsack.
    if (glob.rank == 0) {
        struct crtalgsqrt_knapsack cks[1];
        crtalgsqrt_knapsack_init(cks);
        crtalgsqrt_knapsack_prepare(cks, (nab+(nab&1)) / 2);
        cks->ks->tab = contribs64.data();
        cks->tabN = contribsN.data();
        knapsack_solve(cks->ks);
        crtalgsqrt_knapsack_clear(cks);
    }

    /****************************************************************/

    ab_source_clear(glob.ab);

    MPI_Finalize();

    return 0;
}
