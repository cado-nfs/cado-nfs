/* This file is _NOT_ a standalone compilation unit. It is included by
 * bench_arith_vs_arithxx_old_*.cpp, after one of the mod*_default.h
 * headers of the old arith layer, and written with the generic mod_*
 * names.
 */
#include "cado.h" // IWYU pragma: keep

#ifndef BENCH_ARITH_VS_ARITHXX_READY_TO_INCLUDE_IMPL_CODE
#error "This file must not be used as a standalone compilation unit"
#endif

#include <climits>
#include <cstddef>
#include <cstdint>

#include <gmp.h>

#include "bench_arith_vs_arithxx.hpp"
#include "cxx_mpz.hpp"

// scan-headers: skip

/* low 64 bits of a residue, in its plain (not Montgomery) form */
static uint64_t MOD_APPEND_TYPE(low64)(residue_t const x, modulus_t const m)
{
    modint_t z;
    mod_intinit(z);
    mod_get_int(z, x, m);
    unsigned long w[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    size_t const n = mod_intget_uls(w, z);
    uint64_t r = 0;
    for (size_t i = 0; i < n && i * ULONG_BITS < 64; i++)
        r |= uint64_t(w[i]) << (i * ULONG_BITS);
    mod_intclear(z);
    return r;
}

kbench_results MOD_APPEND_TYPE(bench_old)(cxx_mpz const & Mz, kbench_iters const & it)
{
    kbench_results res;

    unsigned long w[8];
    size_t n;
    mpz_export(w, &n, -1, sizeof(unsigned long), 0, 0, Mz);
    modint_t M, g;
    modulus_t m;
    mod_intinit(M);
    mod_intinit(g);
    mod_intset_uls(M, w, n);
    mod_initmod_int(m, M);

    residue_t x, y, z, px, pz, qx, qz, dx, dz, u, v, t;
    residue_t * all[] = { &x, &y, &z, &px, &pz, &qx, &qz, &dx, &dz, &u, &v, &t };
    for (auto * r : all)
        mod_init(*r, m);

    auto reset = [&]() {
        /* y and z are full-size: small operands would favour mod_mpz,
         * whose residues shrink with their value */
        mod_set_ul(x, 3, m); mod_set_ul(y, 5, m); mod_set_ul(z, 7, m);
        mod_inv(y, y, m); mod_inv(z, z, m);
        mod_set_ul(px, 11, m); mod_set_ul(pz, 13, m);
        mod_set_ul(qx, 17, m); mod_set_ul(qz, 19, m);
        mod_set_ul(dx, 23, m); mod_set_ul(dz, 29, m);
    };
    auto low64 = [&](residue_t const a) { return MOD_APPEND_TYPE(low64)(a, m); };
    double ns;

    if (kbench_wanted(it, "mul")) {
        reset();
        ns = kbench_time(it.mul, [&]() {
                for (size_t i = 0; i < it.mul; i++)
                    mod_mul(x, x, y, m);
                });
        res.push_back({ "mul", ns, low64(x) });
    }

    if (kbench_wanted(it, "sqr")) {
        reset();
        ns = kbench_time(it.mul, [&]() {
                for (size_t i = 0; i < it.mul; i++)
                    mod_sqr(x, x, m);
                });
        res.push_back({ "sqr", ns, low64(x) });
    }

    if (kbench_wanted(it, "add+sub")) {
        reset();
        ns = kbench_time(it.mul, [&]() {
                for (size_t i = 0; i < it.mul; i++) {
                    mod_add(x, x, y, m);
                    mod_sub(x, x, z, m);
                }
                });
        res.push_back({ "add+sub", ns, low64(x) });
    }

    if (kbench_wanted(it, "dadd")) {
        reset();
        ns = kbench_time(it.mul, [&]() {
                for (size_t i = 0; i < it.mul; i++) {
                    /* differential addition on a Montgomery curve, 4M+2S */
                    mod_sub(u, px, pz, m);
                    mod_add(v, qx, qz, m);
                    mod_mul(u, u, v, m);
                    mod_add(t, px, pz, m);
                    mod_sub(v, qx, qz, m);
                    mod_mul(v, t, v, m);
                    mod_add(t, u, v, m);
                    mod_sub(v, u, v, m);
                    mod_sqr(t, t, m);
                    mod_sqr(v, v, m);
                    mod_mul(px, t, dz, m);
                    mod_mul(pz, v, dx, m);
                }
                });
        res.push_back({ "dadd", ns, low64(px) ^ low64(pz) });
    }

    if (kbench_wanted(it, "div3")) {
        reset();
        ns = kbench_time(it.mul, [&]() {
                for (size_t i = 0; i < it.mul; i++)
                    mod_div3(x, x, m);
                });
        res.push_back({ "div3", ns, low64(x) });
    }

    if (kbench_wanted(it, "pow")) {
        reset();
        ns = kbench_time(it.pow, [&]() {
                for (size_t i = 0; i < it.pow; i++) {
                    mod_pow_ul(x, x, kbench_pow_exponent, m);
                    mod_add1(x, x, m);
                }
                });
        res.push_back({ "pow", ns, low64(x) });
    }

    if (kbench_wanted(it, "2pow")) {
        reset();
        unsigned long e = kbench_pow2_exponent;
        ns = kbench_time(it.pow, [&]() {
                for (size_t i = 0; i < it.pow; i++) {
                    mod_2pow_ul(x, e, m);
                    e = kbench_next_exponent(e);
                }
                });
        res.push_back({ "2pow", ns, low64(x) });
    }

    if (kbench_wanted(it, "inv")) {
        reset();
        ns = kbench_time(it.inv, [&]() {
                for (size_t i = 0; i < it.inv; i++) {
                    mod_inv(x, x, m);
                    mod_add1(x, x, m);
                }
                });
        res.push_back({ "inv", ns, low64(x) });
    }

    if (kbench_wanted(it, "gcd")) {
        reset();
        unsigned long gacc = 0;
        ns = kbench_time(it.inv, [&]() {
                for (size_t i = 0; i < it.inv; i++) {
                    mod_gcd(g, x, m);
                    gacc += mod_intget_ul(g);
                    mod_add(x, x, y, m);
                }
                });
        res.push_back({ "gcd", ns, gacc });
    }

    if (kbench_wanted(it, "isprime")) {
        /* called through a volatile pointer: mpz_probab_prime_p is declared
         * pure, and would otherwise be hoisted out of the loop */
        int (* volatile isprime)(modulus_t const) = [](modulus_t const mm) {
            return int(mod_isprime(mm));
        };
        unsigned long pacc = 0;
        ns = kbench_time(it.prime, [&]() {
                for (size_t i = 0; i < it.prime; i++)
                    pacc += isprime(m);
                });
        res.push_back({ "isprime", ns, pacc });
    }

    for (auto * r : all)
        mod_clear(*r, m);
    mod_clearmod(m);
    mod_intclear(M);
    mod_intclear(g);

    return res;
}
