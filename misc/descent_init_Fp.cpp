#include "cado.h" // IWYU pragma: keep

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <map>
#include <algorithm>
#include <utility>
#include <vector>

#include <unistd.h>
#include <gmp.h>
#include "fmt/base.h"


#include "cado_poly.h"
#include "ecm.h"
#include "lll.h"
#include "macros.h"
#include "mpz_poly.h"
#include "smooth_detect.hpp"
#include "mpz_mat.h"
#include "cxx_mpz.hpp"
#include "params.h"

static double default_B1done;

/*
ECM=...
gcc -Wall -g -std=c99 -I${ECM}/include descent_init_Fp.c smooth_detect.c -o
descent_init_Fp -L${ECM}/lib -lecm -lgmp -lm -lpthread
*/

static std::mutex mut_found;
static std::condition_variable cond_found;
static std::deque<std::pair<std::thread::id, descent_init_candidate>> winners;

class semaphore {
    std::mutex lk;
    bool b = false;
    public:
    operator bool() {
        std::lock_guard<std::mutex> const dummy(lk);
        return b;
    }
    void raise() {
        std::lock_guard<std::mutex> const dummy(lk);
        b = true;
    }
};

static semaphore please_die;

static void
HalfGcd(cxx_mpz & a, cxx_mpz & b, cxx_mpz & u)
{
    cxx_mpz v, w, x, q, r;
    u = 1;
    w = v = 0;
    x = q;
    /* invariant: a = u*a0 + v*b0 */
    while (mpz_cmpabs(a, u) > 0) {
        mpz_tdiv_qr(q, r, a, b);
        mpz_swap(a, b);
        mpz_swap(b, r);
        mpz_submul(u, q, w);
        mpz_swap(u, w);
        mpz_submul(v, q, x);
        mpz_swap(v, x);
    }
}

struct Fp_param
{
    cxx_mpz p;
    int n;         // extension degree
    cxx_mpz z;       // element to descend (n=1)
    cxx_mpz_poly zpol; // element to descend (n>1)
    cxx_mpz_poly f;    // used only in JL mode
    cxx_mpz_poly phi;  // used only in the n>1 case
    cxx_mpz m;       // root of f mod p (the one common with g) when n=1
};

static int
next_cand_Fp_hgcd(descent_init_candidate & cand, const void * params)
{
    if (please_die) {
        return 0;
    }
    Fp_param const & param = * (Fp_param const *) params;
    unsigned long const e = random();
    cxx_mpz u0, v0;
    mpz_powm_ui(u0, param.z, e, param.p);
    cxx_mpz tmp = param.p;
    HalfGcd(u0, tmp, v0);
    mpz_abs(u0, u0);
    mpz_abs(v0, v0);
    cand = descent_init_candidate(u0, v0, e);
    return 1;
}

static int
is_probably_sqrfree(cxx_mpz const & z)
{
    for (unsigned long const p :
            { 2,  3,  5,  7,  11, 13, 17, 19,
              23, 29, 31, 37, 41, 43, 47 })
    {
        unsigned long const p2 = p * p;
        if (mpz_gcd_ui(nullptr, z, p2) == p2)
            return 0;
    }
    return 1;
}

// JL version
// Returns a boolean meaning "failure, try again".
static int
get_JL_candidate_from_e(unsigned long e,
                       cxx_mpz_poly & U,
                       cxx_mpz_poly & V,
                       cxx_mpz & u,
                       cxx_mpz & v,
                       Fp_param const & param)
{
    cxx_mpz ze;
    mpz_powm_ui(ze, param.z, e, param.p);

    //**** Create matrix (see Barbulescu's PhD thesis, section 8.4.2)
    const int d = param.f->deg;

    cxx_mpz_mat M(2 * d, 2 * d);

    mpz_set(mpz_mat_entry(M, 0,0), param.p);

    // Topleft d*d block
    for (int i = 1; i < d; ++i) {
        mpz_set_ui(mpz_mat_entry(M, i, i), 1);
        mpz_neg(mpz_mat_entry(M, i, i-1), param.m);
    }
    // Bottomleft d*d block
    for (int i = 0; i < d; ++i)
        mpz_set(mpz_mat_entry(M, d+i, i), ze);

    // Bottomright d*d block
    for (int i = 0; i < d; ++i)
        mpz_set_ui(mpz_mat_entry(M, d+i, d+i), 1);

    //**** Apply LLL.
    {
        cxx_mpz det, a, b;
        a = b = 1;
        mpz_mat_LLL(det, M, nullptr, a, b);
    }

    //**** Recover rational reconstruction
    // z^e = U(alpha)/V(alpha) mod (p, x-m)
    for (int i = 0; i < d; ++i) {
        mpz_poly_setcoeff(U, i, mpz_mat_entry_const(M, 0, i));
        mpz_poly_setcoeff(V, i, mpz_mat_entry_const(M, 0, d+i));
    }
    mpz_poly_cleandeg(U, d - 1);
    mpz_poly_cleandeg(V, d - 1);

    //**** Compute norms
    mpz_poly_resultant(u, U, param.f);
    mpz_poly_resultant(v, V, param.f);
    mpz_abs(u, u);
    mpz_abs(v, v);
    // TODO: we don't want to deal with ideals of degree > 1.
    // The simplest way is to forbid squares at all. This is not optimal,
    // and furthermore, we can do it only for small factors...
    int fail = 0;
    if (!is_probably_sqrfree(u) || !is_probably_sqrfree(v))
        fail = 1;

    return fail;
}

// GF(p^n) version.
// Returns a boolean meaning "failure, try again".
static int
get_Fpn_candidate_from_e(unsigned long e,
                        cxx_mpz_poly & U,
                        cxx_mpz_poly & V,
                        cxx_mpz & u,
                        cxx_mpz & v,
                        Fp_param const & param)
{
    const int d = param.f->deg;
    const int n = param.phi->deg;

    cxx_mpz_poly ze;
    mpz_poly_pow_ui_mod_f_mod_mpz(ze, param.zpol, param.f, e, param.p);

    cxx_mpz_poly y; // y represents x^i
    cxx_mpz_poly c; // c represents ze*x^i mod f

    //**** Create matrix (see Barbulescu's PhD thesis, section 8.4.2)
    // Warning: the LLL code count indices starting with 1 and not 0.
    // (this is a bug, imho)
    mat_Z M;
    LLL_init(&M, 2 * d, 2 * d); // allocate and set to 0.

    // Topleft d*d block
    for (int i = 1; i < n + 1; ++i)
        mpz_set(M.coeff[i][i], param.p);
    // Let phi be the degree n polynomial defining the extension field
    // put the coef of phi=phi[0]+phi[1]*x+...+phi[n-1]*x^(n-1)+x^n
    for (int i = n + 1; i < d + 1; ++i) {
        for (int j = 0; j < n + 1; ++j)
            mpz_set(M.coeff[i][i - n + j], mpz_poly_coeff_const(param.phi, j));
    }

    // Bottomleft d*d block
    for (int i = 0; i < d; ++i) {
        mpz_poly_set_xi(y, i); // y is the polynomial alpha^i
        mpz_poly_mul_mod_f_mod_mpz(c, ze, y, param.f, param.p, nullptr, nullptr);
        for (int j = 0; j < d; ++j) {
            if (j > c->deg)
                mpz_set_ui(M.coeff[d + i + 1][j + 1], 0);
            else
                mpz_set(M.coeff[d + i + 1][j + 1], mpz_poly_coeff_const(c, j));
        }
    }

    // Bottomright d*d block
    for (int i = 0; i < d; ++i) {
        mpz_set_ui(M.coeff[d + i + 1][d + i + 1], 1);
    }

    //**** Apply LLL.
    {
        cxx_mpz det, a, b;
        a = b = 1;
        LLL(det, M, nullptr, a, b);
    }

    //**** Recover rational reconstruction
    // z^e = U(alpha)/V(alpha) mod (p, x-m)
    for (int i = 0; i < d; ++i) {
        mpz_poly_setcoeff(U, i, M.coeff[1][i + 1]);
        mpz_poly_setcoeff(V, i, M.coeff[1][d + i + 1]);
    }
    mpz_poly_cleandeg(U, d - 1);
    mpz_poly_cleandeg(V, d - 1);

    //**** Compute norms
    mpz_poly_resultant(u, U, param.f);
    mpz_poly_resultant(v, V, param.f);

    mpz_abs(u, u);
    mpz_abs(v, v);
    // TODO: we don't want to deal with ideals of degree > 1.
    // The simplest way is to forbid squares at all. This is not optimal,
    // and furthermore, we can do it only for small factors...
    int fail = 0;
    if (!is_probably_sqrfree(u) || !is_probably_sqrfree(v))
        fail = 1;

    LLL_clear(&M);
    return fail;
}

static int
next_cand_nonrat(descent_init_candidate & cand, const void * params)
{
    if (please_die) {
        return 0;
    }

    Fp_param const & param = * (Fp_param const *)params;
    cxx_mpz u, v;
    cxx_mpz_poly U, V;
    unsigned long e;
    int fail;
    do {
        e = random();
        if (param.n == 1) {
            fail = get_JL_candidate_from_e(e, U, V, u, v, param);
        } else {
            fail = get_Fpn_candidate_from_e(e, U, V, u, v, param);
        }
    } while (fail);

    cand = descent_init_candidate(u, v, e);
    return 1;
}

// Full factorization of z0; non-optimized.
// Assume fac_z has been allocated.
// Returns the number of factors.
static void
full_factor(std::vector<cxx_mpz> & fac_z, cxx_mpz const & z0)
{
    double B1 = 100.0;
    long sig;
    cxx_mpz z = z0;

    // Remove small primes, ECM can't separate them
    for (unsigned long const p : {
            2,  3,  5,  7,  11, 13, 17, 19,
            23, 29, 31, 37, 41, 43, 47 })
    {
        while (mpz_divisible_ui_p(z, p)) {
            fac_z.emplace_back(p);
            mpz_divexact_ui(z, z, p);
        }
    }

    while (!mpz_probab_prime_p(z, 10)) {
        cxx_mpz f;
        int success = 0;
        while (!success) {
            ecm_params ecm_par;
            ecm_init(ecm_par);
            ecm_par->B1done = default_B1done; /* issue with ECM 6.4.x */
            sig = random();
            mpz_set_ui(ecm_par->sigma, sig);
            success = ecm_factor(f, z, B1, ecm_par);
            B1 += sqrt(B1);
            ecm_clear(ecm_par);
        }
        if (mpz_perfect_power_p(f)) {
            cxx_mpz ff;
            for (int pow = 2;; pow++) {
                if (mpz_root(ff, f, pow)) {
                    mpz_set(f, ff);
                    break;
                }
            }
        }
        if (!mpz_probab_prime_p(f, 10)) {
            B1 -= 3 * sqrt(B1);
            B1 = std::max(B1, 20.0);
            continue;
        }
        do {
            fac_z.push_back(f);
            mpz_divexact(z, z, f);
        } while (mpz_divisible_p(z, f));
    }
    fac_z.push_back(z);
    std::ranges::sort(fac_z);
}

// Check if there are multiple factors.
// This assumes that the factors are sorted, so that multiple factors are
// consecutive.
static int
has_distinct_factors(std::vector<cxx_mpz> const & P)
{
    for (size_t i = 1; i < P.size() ; i++) {
        if (P[i-1] == P[i])
            return false;
    }
    return true;
}

static cxx_mpz
find_root(cxx_mpz const & p, cxx_mpz_poly const & f1, cxx_mpz_poly const & f2)
{
    // Check if projective root
    cxx_mpz r;
    mpz_mod(r, mpz_poly_lc(f1), p);
    if (mpz_cmp_ui(r, 0) == 0) {
        mpz_mod(r, mpz_poly_lc(f2), p);
        if (mpz_cmp_ui(r, 0) == 0) {
            return p;
        }
    }

    // Non projective
    cxx_mpz_poly G;
    mpz_poly_gcd_mpz(G, f1, f2, p);
    ASSERT_ALWAYS(G->deg == 1);
    mpz_invert(r, mpz_poly_coeff_const(G, 1), p);
    mpz_mul(r, r, mpz_poly_coeff_const(G, 0));
    mpz_neg(r, r);
    mpz_mod(r, r, p);
    return r;
}

// Possible modes are
//   MODE_RAT  (when there is a rational side)
//   MODE_JL   (for GF(p), with Joux-Lercier)
//   MODE_FPN  (for GF(p^n))
#define MODE_RAT 1
#define MODE_JL 2
#define MODE_FPN 3

struct descent_thread_param {
    Fp_param const & params;
    smooth_detect_params smooth_param;
    unsigned long target;
    int mode;
    unsigned int id;

    void operator()() const;
};

static void one_descent_thread(
        Fp_param const & params,
        smooth_detect_params const & smooth_param,
        unsigned long target,
        int mode)
{
    descent_init_candidate C;
    if (mode == MODE_RAT) {
        C = smooth_detect(
                      next_cand_Fp_hgcd,
                      (const void *) &params,
                      target,
                      smooth_param);
    } else if (mode == MODE_JL || mode == MODE_FPN) {
        C = smooth_detect(
                      next_cand_nonrat,
                      (const void *) &params,
                      target,
                      smooth_param);
    }

    std::lock_guard<std::mutex> const dummy(mut_found);

    winners.emplace_back(std::this_thread::get_id(), C);
    cond_found.notify_one();
}

static void descent_declare_usage(cxx_param_list & pl)
{
    param_list_usage_header(pl,
            "[-poly polfile] [-side xxx] [-extdeg n] [-jl] [-mt n] [-mineff "
            "e] [-maxeff E] [-seed s] [-lpb t] [-v] p z\n"
            "  If extdeg > 1, then z must be a white-separated sequence of "
            "coefs z0 z1 ... z_{k-1}\n");

    param_list_decl_usage(pl, "seed", "random seed");
    param_list_decl_usage(pl, "mt", "number of threads");
    param_list_decl_usage(pl, "minB1", "start ECM with this B1");
    param_list_decl_usage(pl, "mineff", "minimum ECM effort");
    param_list_decl_usage(pl, "maxeff", "maximum ECM effort");
    param_list_decl_usage(pl, "side", "side for target");
    param_list_decl_usage(pl, "v", "verbose mode");
    param_list_decl_usage(pl, "extdeg", "do descent for GF(p^extdeg)");
    param_list_decl_usage(pl, "jl", "use Joux-Lercier sieving");
    param_list_decl_usage(pl, "lpb", "bound on large primes");
    param_list_decl_usage(pl, "poly", "cado-nfs polynomial");
}

namespace {
    int verbose;
    int jl;
}

static void descent_configure_switches(cxx_param_list & pl)
{
    param_list_configure_switch(pl, "-v", &verbose);
    param_list_configure_switch(pl, "-jl", &jl);
}


int
main(int argc0, char const * argv0[])
{
    int argc = argc0;
    char const ** argv = argv0;

    unsigned long seed = 0;
    unsigned long target = 0; // the target smoothness
    unsigned long nthread = 1;
    double mineff = 2000.0;
    double maxeff = 1e20;
    double minB1 = 100.0;
    unsigned int ext = 1; // extension degree
    int side = 1;
    clock_t const tm = clock();

    cxx_param_list pl;

    descent_declare_usage(pl);
    descent_configure_switches(pl);

    cxx_cado_poly cpoly;

    std::vector<cxx_mpz> wild;

    argv++, argc--;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        if (argv[0][0] != '-') {
            cxx_mpz z;
            mpz_set_str(z, argv[0], 0);
            wild.push_back(z);
            argv++, argc--;
            continue;
        }
        fmt::print(stderr, "Unhandled parameter {}\n", argv[0]);
        param_list_print_usage(pl, argv0[0], stderr);
        return EXIT_FAILURE;
    }

    param_list_parse(pl, "seed", seed);
    param_list_parse(pl, "mt", nthread);
    param_list_parse(pl, "minB1", minB1);
    param_list_parse(pl, "mineff", mineff);
    param_list_parse(pl, "maxeff", maxeff);
    param_list_parse(pl, "side", side);
    param_list_parse(pl, "extdeg", ext);
    param_list_parse(pl, "lpb", target);
    const char * polyfilename = param_list_lookup_string(pl, "poly");
    if (polyfilename) {
        cado_poly_read(cpoly, polyfilename);
    }

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, argv0[0], stderr);
        return EXIT_FAILURE;
    }

    if ((jl || (ext > 1)) && !polyfilename) {
        fmt::print(
          stderr,
          "Error, must provide -poly when extdeg > 1 or using -jl option\n");
        param_list_print_usage(pl, argv0[0], stderr);
        return EXIT_FAILURE;
    }

    if (ext > 1 && jl) {
        fmt::print(stderr, "Warning: ignoring the -jl option with extdeg > 1\n");
        jl = 0;
    }

    if (wild.size() != ext + 1) {
        fmt::print(stderr, "Error: for extension degree {}, we need {} tail arguments\n",
                ext, ext + 1);
        return EXIT_FAILURE;
    }

    if (seed == 0)
        seed = getpid() + (time(nullptr) << 16U);
    srandom(seed);

    Fp_param params;
    params.p = wild[0];
    params.n = ext;
    if (jl || ext > 1) {
        params.f = cpoly->pols[side];
    }
    if (jl) {
        int const ret = cado_poly_getm(params.m, cpoly, params.p);
        ASSERT_ALWAYS(ret);
    }
    if (ext > 1) {
        mpz_poly_gcd_mpz(params.phi, cpoly->pols[0], cpoly->pols[1], params.p);
        for (unsigned int i = 0; i < ext; i++)
            mpz_poly_setcoeff(params.zpol, i, wild[i+1]);
    } else {
        params.z = wild[1];
    }

    const smooth_detect_params smooth_param = {
        mineff, maxeff, 10, verbose, minB1
    };


    int const mode = (ext > 1) ? MODE_FPN : (jl ? MODE_JL : MODE_RAT);

    std::map<std::thread::id, std::thread> threads;

    for( ; threads.size() < nthread ; ) {
        auto th = std::thread(one_descent_thread, params, smooth_param, target, mode);
        std::thread::id const id = th.get_id();
        threads[id] = std::move(th);
    }

    // Wait for one thread to signal success.
    for( ;; ) {
        std::pair<std::thread::id, descent_init_candidate> winner;
        {
            std::unique_lock<std::mutex> lk(mut_found);
            cond_found.wait(lk, [](){return !winners.empty();});
            std::swap(winner, winners.front());
            winners.pop_front();
        }

        threads[winner.first].join();
        threads.erase(winner.first);

        if (jl || ext > 1) {
            cxx_mpz u, v;
            cxx_mpz_poly U, V;
            // do again the LLL thing and print result.
            int fail;
            if (ext == 1)
                fail = get_JL_candidate_from_e(winner.second.e, U, V, u, v, params);
            else
                fail = get_Fpn_candidate_from_e(winner.second.e, U, V, u, v, params);
            ASSERT_ALWAYS(!fail);

            std::vector<cxx_mpz> facu, facv;
            full_factor(facu, u);
            full_factor(facv, v);
            if (!has_distinct_factors(facu) ||
                !has_distinct_factors(facv))
            {
                fmt::print("Fail: non-squarefree norm\n");
                // one of them is not squarefree. Restart the thread and wait
                // for another candidate.
                auto th = std::thread(one_descent_thread, params, smooth_param, target, mode);
                std::thread::id const id = th.get_id();
                threads[id] = std::move(th);
                continue;
            }

            fmt::print("U = {}\nV = {}\nu = {}\nv = {}\n",
                    mpz_poly_coeff_list(U, ","),
                    mpz_poly_coeff_list(V, ","),
                    u, v);

            fmt::print("fac_u =");
            for(cxx_mpz const & p : facu) {
                fmt::print(" {},{}", p, find_root(p, U, params.f));
            }
            fmt::print("\n");

            fmt::print("fac_v =");
            for(cxx_mpz const & p : facv) {
                fmt::print(" {},{}", p, find_root(p, V, params.f));
            }
            fmt::print("\n");

        }
        fmt::print("Youpi: e = {} is a winner\n", winner.second.e);
        break;
    }

    // Cancel other threads
    please_die.raise();

    for(auto & t : threads)
        t.second.join();
    threads.clear();

    fmt::print("Total CPU time: {:.1f} s\n",
           ((double)(clock() - tm)) / CLOCKS_PER_SEC);

    return EXIT_SUCCESS;
}
