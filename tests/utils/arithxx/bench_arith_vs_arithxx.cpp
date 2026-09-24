/* Kernel benchmark of the old arith layer against arithxx.
 *
 * For a given modulus size, the old layer is the one that facul picks
 * for that size (modredc_ul, modredc_15ul, modredc_2ul2 or mod_mpz), and
 * the new layer is the smallest of modredc64, modredc96, modredc126 and
 * mod_mpz_new that fits. The modulus is the largest prime below 2^bits.
 *
 * Each kernel is a chain of dependent operations (mul, sqr, add+sub, a
 * Montgomery curve differential addition, div3, pow, 2pow, inv, gcd,
 * isprime), run identically on both sides. The final values must agree;
 * the program fails if they do not. With -quick, the iteration counts
 * are small enough to run it as a test.
 *
 * Usage: bench_arith_vs_arithxx [-quick] [-only <op>,<op>...] [-new-first]
 *                               [<bits> ...]
 *
 * -only restricts the run to the given kernels. Both layers are timed
 * twice, in the order old, new, new, old, and the best time of each is
 * kept: on some machines, what runs second in a process is consistently
 * slower (by 2.5% on a Skylake, for 64-bit mul). -new-first times new,
 * then old, once, to observe this.
 *
 * Without <bits>, every size at which either layer changes is tried.
 * Timings are best of 5, in nanoseconds per iteration, on one thread:
 * pin the process to a core (e.g. with taskset), and do not run several
 * copies on the same core.
 */

#include "cado.h" // IWYU pragma: keep

#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <string>
#include <vector>

#include <gmp.h>

#include "arith/mod_mpz.h"
#include "arith/modredc_15ul.h"
#include "arith/modredc_2ul2.h"
#include "arith/modredc_ul.h"
#include "arithxx/mod_mpz_new.hpp"
#include "arithxx/modredc126.hpp"
#include "arithxx/modredc64.hpp"
#include "arithxx/modredc96.hpp"
#include "bench_arith_vs_arithxx.hpp"
#include "cxx_mpz.hpp"
#include "macros.h"

uint64_t kbench_low64(mpz_srcptr z)
{
    uint64_t r = 0;
    for (int i = 0; i * GMP_LIMB_BITS < 64; i++)
        r |= uint64_t(mpz_getlimbn(z, i)) << (i * GMP_LIMB_BITS);
    return r;
}

/* low bits of an arithxx Integer, without going through a temporary */
static uint64_t low64(cxx_mpz const & z) { return kbench_low64(z); }
template<typename I> static uint64_t low64(I const & z) { return z[0]; }

template<typename layer>
static kbench_results bench_new(cxx_mpz const & Mz, kbench_iters const & it)
{
    using Modulus = typename layer::Modulus;
    using Residue = typename layer::Residue;
    using Integer = typename layer::Integer;

    kbench_results res;

    Modulus const m { Integer(Mz) };
    Residue x(m), y(m), z(m), px(m), pz(m), qx(m), qz(m), dx(m), dz(m),
        u(m), v(m), t(m);

    auto reset = [&]() {
        /* y and z are full-size: small operands would favour mod_mpz,
         * whose residues shrink with their value */
        m.set(x, uint64_t(3)); m.set(y, uint64_t(5)); m.set(z, uint64_t(7));
        m.inv(y, y); m.inv(z, z);
        m.set(px, uint64_t(11)); m.set(pz, uint64_t(13));
        m.set(qx, uint64_t(17)); m.set(qz, uint64_t(19));
        m.set(dx, uint64_t(23)); m.set(dz, uint64_t(29));
    };
    auto get64 = [&](Residue const & a) { return low64(m.get(a)); };
    double ns;

    if (kbench_wanted(it, "mul")) {
        reset();
        ns = kbench_time(it.mul, [&]() {
                for (size_t i = 0; i < it.mul; i++)
                    m.mul(x, x, y);
                });
        res.push_back({ "mul", ns, get64(x) });
    }

    if (kbench_wanted(it, "sqr")) {
        reset();
        ns = kbench_time(it.mul, [&]() {
                for (size_t i = 0; i < it.mul; i++)
                    m.sqr(x, x);
                });
        res.push_back({ "sqr", ns, get64(x) });
    }

    if (kbench_wanted(it, "add+sub")) {
        reset();
        ns = kbench_time(it.mul, [&]() {
                for (size_t i = 0; i < it.mul; i++) {
                    m.add(x, x, y);
                    m.sub(x, x, z);
                }
                });
        res.push_back({ "add+sub", ns, get64(x) });
    }

    if (kbench_wanted(it, "dadd")) {
        reset();
        ns = kbench_time(it.mul, [&]() {
                for (size_t i = 0; i < it.mul; i++) {
                    /* differential addition on a Montgomery curve, 4M+2S */
                    m.sub(u, px, pz);
                    m.add(v, qx, qz);
                    m.mul(u, u, v);
                    m.add(t, px, pz);
                    m.sub(v, qx, qz);
                    m.mul(v, t, v);
                    m.add(t, u, v);
                    m.sub(v, u, v);
                    m.sqr(t, t);
                    m.sqr(v, v);
                    m.mul(px, t, dz);
                    m.mul(pz, v, dx);
                }
                });
        res.push_back({ "dadd", ns, get64(px) ^ get64(pz) });
    }

    if (kbench_wanted(it, "div3")) {
        reset();
        ns = kbench_time(it.mul, [&]() {
                for (size_t i = 0; i < it.mul; i++)
                    m.div3(x, x);
                });
        res.push_back({ "div3", ns, get64(x) });
    }

    if (kbench_wanted(it, "pow")) {
        reset();
        ns = kbench_time(it.pow, [&]() {
                for (size_t i = 0; i < it.pow; i++) {
                    m.pow(x, x, uint64_t(kbench_pow_exponent));
                    m.add1(x, x);
                }
                });
        res.push_back({ "pow", ns, get64(x) });
    }

    if (kbench_wanted(it, "2pow")) {
        reset();
        unsigned long e = kbench_pow2_exponent;
        ns = kbench_time(it.pow, [&]() {
                for (size_t i = 0; i < it.pow; i++) {
                    m.pow2(x, uint64_t(e));
                    e = kbench_next_exponent(e);
                }
                });
        res.push_back({ "2pow", ns, get64(x) });
    }

    if (kbench_wanted(it, "inv")) {
        reset();
        ns = kbench_time(it.inv, [&]() {
                for (size_t i = 0; i < it.inv; i++) {
                    m.inv(x, x);
                    m.add1(x, x);
                }
                });
        res.push_back({ "inv", ns, get64(x) });
    }

    if (kbench_wanted(it, "gcd")) {
        reset();
        unsigned long gacc = 0;
        Integer g;
        ns = kbench_time(it.inv, [&]() {
                for (size_t i = 0; i < it.inv; i++) {
                    m.gcd(g, x);
                    gacc += (unsigned long) low64(g);
                    m.add(x, x, y);
                }
                });
        res.push_back({ "gcd", ns, gacc });
    }

    if (kbench_wanted(it, "isprime")) {
        unsigned long pacc = 0;
        ns = kbench_time(it.prime, [&]() {
                for (size_t i = 0; i < it.prime; i++)
                    pacc += m.is_prime();
                });
        res.push_back({ "isprime", ns, pacc });
    }

    return res;
}

struct side {
    char const * name;
    kbench_results (*run)(cxx_mpz const &, kbench_iters const &);
};

/* the choice that facul makes, see FaculModulusBase::init_mpz */
static side old_layer(unsigned int bits)
{
    if (bits <= MODREDCUL_MAXBITS)
        return { "modredc_ul", bench_old_ul };
    if (bits <= MODREDC15UL_MAXBITS)
        return { "modredc_15ul", bench_old_15ul };
    if (bits <= MODREDC2UL2_MAXBITS)
        return { "modredc_2ul2", bench_old_2ul2 };
    return { "mod_mpz", bench_old_mpz };
}

static side new_layer(unsigned int bits)
{
    if (bits <= 64)
        return { "modredc64", bench_new<arithxx_modredc64> };
    if (bits <= 96)
        return { "modredc96", bench_new<arithxx_modredc96> };
    if (bits <= 126)
        return { "modredc126", bench_new<arithxx_modredc126> };
    return { "mod_mpz_new", bench_new<arithxx_mod_mpz_new> };
}

// coverity[root_function]
int main(int argc, char const * argv[])
{
    bool quick = false;
    bool new_first = false;
    std::string only;
    std::vector<unsigned int> sizes;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-quick") == 0)
            quick = true;
        else if (strcmp(argv[i], "-new-first") == 0)
            new_first = true;
        else if (strcmp(argv[i], "-only") == 0 && i + 1 < argc)
            only = argv[++i];
        else
            sizes.push_back(strtoul(argv[i], nullptr, 10));
    }
    if (sizes.empty()) {
        sizes = { MODREDCUL_MAXBITS, MODREDC15UL_MAXBITS,
                  MODREDC2UL2_MAXBITS, 64, 96, 126, 150 };
        std::sort(sizes.begin(), sizes.end());
        sizes.erase(std::unique(sizes.begin(), sizes.end()), sizes.end());
    }

    kbench_iters const fast { 2000000, 20000, 200000, 2000, only };
    kbench_iters const slow { 200000, 2000, 20000, 200, only };
    kbench_iters const test { 2000, 20, 200, 2, only };

    int mismatches = 0;
    printf("%5s %-12s %-12s %-8s %10s %10s %8s %s\n", "bits", "old", "new",
           "op", "old(ns)", "new(ns)", "new/old", "check");
    for (auto const bits: sizes) {
        ASSERT_ALWAYS(bits >= 2);
        cxx_mpz M;
        mpz_setbit(M, bits);
        do {
            mpz_sub_ui(M, M, 1);
        } while (!mpz_probab_prime_p(M, 25));

        side const o = old_layer(bits);
        side const n = new_layer(bits);
        kbench_iters const & it =
            quick ? test : (bits > MODREDC2UL2_MAXBITS || bits > 126) ? slow : fast;
        kbench_results ro, rn;
        if (new_first) {
            rn = n.run(M, it);
            ro = o.run(M, it);
        } else {
            ro = o.run(M, it);
            rn = n.run(M, it);
            auto const rn2 = n.run(M, it);
            auto const ro2 = o.run(M, it);
            for (size_t k = 0; k < ro.size(); k++) {
                ro[k].ns = std::min(ro[k].ns, ro2[k].ns);
                rn[k].ns = std::min(rn[k].ns, rn2[k].ns);
            }
        }
        ASSERT_ALWAYS(ro.size() == rn.size());
        for (size_t k = 0; k < ro.size(); k++) {
            ASSERT_ALWAYS(ro[k].op == rn[k].op);
            bool const ok = ro[k].chk == rn[k].chk;
            mismatches += !ok;
            printf("%5u %-12s %-12s %-8s %10.2f %10.2f %8.3f %s\n", bits,
                   o.name, n.name, ro[k].op.c_str(), ro[k].ns, rn[k].ns,
                   rn[k].ns / ro[k].ns, ok ? "match" : "MISMATCH");
        }
        fflush(stdout);
    }
    return mismatches ? EXIT_FAILURE : EXIT_SUCCESS;
}
