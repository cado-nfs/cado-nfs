#include "cado.h"
#include <sys/types.h>
#include <unistd.h>
#include <stdint.h>
#include <vector>
#include <tuple>
#include <string>
#include <map>
#include <stdio.h>
#include "gcd.h"
#include "macros.h"
#include "getprime.h"
#include "fb-types.h"
#include "las-plattice.hpp"

/* see plattice.sage */

/* This new c code is at very least not slower, and possibly even mildly
 * faster than the old assembly or C code. However the measurement is not
 * easy. Note also that this code covers many more cases.
 *
 * make -j8 test-reduce-plattice && for i in {1..10} ; do ./build/houblon/tests/sieve/test-reduce-plattice -I 65536 -ntests 10000000 -B $((2**25)) -T -seed 255775 ; done
 *
 * new code / two_legs # 5.8378
 * new code / simplistic # 6.0189
 * new code / swapping_loop # 5.9685
 * reference # 6.3625 
 * reference2 # 5.8297
 * reference2_asm # 5.9573
 */

struct plattice : public plattice_info {
    /* use this default ctor only for the comparison with the reference
     * routines.
     */
    plattice () = default;

    using plattice_info::mi0;
    using plattice_info::j0;
    using plattice_info::i1;
    using plattice_info::j1;
    using plattice_info::check_post_conditions;

    void simplistic(uint32_t I) {
        /* This is the main reduce_plattice loop */
        for( ;; ) {
            if (i1 < I) {
                if (i1 == 0) {
                    // Lo=matrix([ (mi0, j1-j0), (i1, j1)])
                    j0 = j1 - j0;
                    lattice_with_vertical_vector(I);
                    return;
                }
                int a = (mi0 + i1 - I) / i1;
                mi0 -= a * i1;
                j0  += a * j1;
                return;
            }
            {
                int k = mi0 / i1; mi0 -= k * i1; j0 += k * j1;
            }
            if (mi0 < I) {
                if (mi0 == 0) {
                    mi0 = i1;
                    i1 = j0 ; j0 = j1 ; j1 = i1;
                    i1 = 0;
                    lattice_with_vertical_vector(I);
                    return;
                }
                int a = (mi0 + i1 - I) / mi0;
                i1 -= a * mi0;
                j1 += a * j0;
                return;
            }
            {
                int k = i1 / mi0; i1 -= k * mi0; j1 += k * j0;
            }
        }
    }

    void swapping_loop(uint32_t I)
    {
        /* This is the main reduce_plattice loop */
        int flip;
        for(flip = 0 ; i1 >= I; flip ^= 1 ) {
            /* do partial unrolling for the frequent case where the
             * quotient is either 1 or 2.
             * this has a significant overall impact
             */
#if 1
            if (mi0 < i1 * 3) {
                { mi0 -= i1; j0 += j1; }
                if (mi0 >= i1) { mi0 -= i1; j0 += j1; }
            } else
#endif
            {
                int k = mi0 / i1; mi0 -= k * i1; j0 += k * j1;
            }
            std::swap(mi0, i1);
            std::swap(j0, j1);
        }
        /* an "UNLIKELY" macro here actually has an adverse
         * effect...  */
        if (i1 == 0) {
            if (!flip) {
                // Lo=matrix([ (mi0, j1-j0), (i1, j1)])
                j0 = j1 - j0;
            } else {
                std::swap(mi0, i1);
                std::swap(j0, j1);
            }
            lattice_with_vertical_vector(I);
        } else {
            int a = (mi0 + i1 - I) / i1;
            mi0 -= a * i1;
            j0  += a * j1;
        }
    }

    /* XXX
     * beware: this constructor takes I, but it shadows a constructor in
     * the production code which takes only logI !!!
     */
    plattice(const unsigned long q, const unsigned long r, bool proj, uint32_t I) : plattice_info()
    {
        initial_basis(q, r, proj);
        /* At this point, (mi0,j0) represents itself, i.e. a vector with
         * two positive coordinates.
         * Note that j0==0
         */
        // ASSERT_ALWAYS(check_pre_conditions(I));
        bool needs_special_treatment = (i1 == 0 || (j1 > 1 && mi0 < I));
        if (needs_special_treatment) {
            lattice_with_vertical_vector(I);
            return;
        }
        two_legs(I);
        // simplistic(I);
        // swapping_loop(I);
    }

};

int
reference (plattice *pli, const fbprime_t p, const fbroot_t r, uint32_t I)
{
    int32_t i0 = -((int32_t) p), i1 = (int32_t) r, j0 = 0, j1 = 1, k;
    const int32_t hI = (int32_t) I;
    const int32_t mhI = -hI;
    while ((i1 >= hI)) {
        k = i0 / i1; i0 %= i1; j0 -= k * j1;
        if ((i0 > mhI)) break;
        k = i1 / i0; i1 %= i0; j1 -= k * j0;
        /* We may conceivably unroll a bit more, or a bit less, here. Just
         * tuck in as many copies of the following block as you wish. */
        if (UNLIKELY(i1 < hI )) break;
        k = i0 / i1; i0 %= i1; j0 -= k * j1;
        if (UNLIKELY(i0 > mhI)) break;
        k = i1 / i0; i1 %= i0; j1 -= k * j0;
    }
    k = i1 - hI - i0;
    if (i1 > -i0) {
        if (UNLIKELY(!i0)) return 0;
        k /= i0; i1 -= k * i0; j1 -= k * j0;
    } else {
        if (UNLIKELY(!i1)) return 0;
        k /= i1; i0 += k * i1; j0 += k * j1;
    }
    pli->mi0 = - (int32_t) i0; pli->j0 = (uint32_t) j0; pli->i1 = (int32_t) i1; pli->j1 = (uint32_t) j1;
    return 1;
}

int
reference2 (plattice *pli, const fbprime_t p, const fbroot_t r, const uint32_t I)
{
    const int32_t hI = (int32_t) I;
    int32_t i0 = - (int32_t) p, i1 = (int32_t) r, j0, j1;

#define RPA do {							\
    i0 += i1; j0 += j1;							\
    if (LIKELY(i0 + i1 * 4 > 0)) {					\
        int64_t c0 = i0, c1 = j0;					\
        c0 += i1; c1 += j1; if (LIKELY(c0 <= 0)) { i0 = c0; j0 = c1; }	\
        c0 += i1; c1 += j1; if (LIKELY(c0 <= 0)) { i0 = c0; j0 = c1; }	\
        c0 += i1; c1 += j1; if (LIKELY(c0 <= 0)) { i0 = c0; j0 = c1; }	\
    } else								\
    RPC;								\
} while (0)
#define RPB do {							\
    i1 += i0; j1 += j0;							\
    if (LIKELY(i1 + i0 * 4 < 0)) {					\
        int64_t c0 = i1, c1 = j1;						\
        c0 += i0; c1 += j0; if (LIKELY(c0 >= 0)) { i1 = c0; j1 = c1; }	\
        c0 += i0; c1 += j0; if (LIKELY(c0 >= 0)) { i1 = c0; j1 = c1; }	\
        c0 += i0; c1 += j0; if (LIKELY(c0 >= 0)) { i1 = c0; j1 = c1; }	\
    } else								\
    RPD;								\
} while (0)
#define RPC do {					\
    int64_t k = i0 / i1; i0 %= i1; j0 -= k * j1;	\
} while (0)
#define RPD do {					\
    int64_t k = i1 / i0; i1 %= i0; j1 -= k * j0;	\
} while (0)

    /* This code seems odd (looks after the i0 <= mhI loop),
       but gcc generates the fastest asm with it... */
    j0 = 0; j1 = 1;
    if (LIKELY(i1 >= hI)) {
        const int32_t mhI = -hI;
        RPC;
        while (UNLIKELY(i0 < -0X7FFFFFFF / 5)) {
            RPD;
            if (UNLIKELY(i1 < 0X7FFFFFFF / 5)) goto p15;
            RPC;
        }
        if (LIKELY(i0 <= mhI)) {
            do {
                RPB;
p15:
                if (UNLIKELY(i1 < hI)) break;
                RPA;
            } while (LIKELY(i0 <= mhI));
        }
    }

#undef RPA
#undef RPB
#undef RPC
#undef RPD

    int64_t k = i1 - hI - i0;
    if (i1 > -i0) {
        if (UNLIKELY(!i0)) return 0;
        k /= i0; i1 -= k * i0; j1 -= k * j0;
    } else {
        if (UNLIKELY(!i1)) return 0;
        k /= i1; i0 += k * i1; j0 += k * j1;
    }
    pli->mi0 = -(int32_t) i0; pli->j0 = (uint32_t) j0; pli->i1 = (int32_t) i1; pli->j1 = (uint32_t) j1;
    return 1;
}

#if defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM) && !(defined(__APPLE_CC__) && defined(__llvm__) && __APPLE_CC__ == 5621)
#define TEST_ASSEMBLY_CODE_DELETED_BY_5f258ce8b
#endif

#ifdef TEST_ASSEMBLY_CODE_DELETED_BY_5f258ce8b
    int
reference2_asm (plattice *pli, const fbprime_t p, const fbroot_t r, const uint32_t I)
{
    const int32_t hI = (int32_t) I;
    int32_t i0 = - (int32_t) p, i1 = (int32_t) r, j0, j1;

    /* Mac OS X 10.8 embarks a version of llvm which crashes on the code
     * below (could be that the constraints are exerting too much of the
     * compiler's behaviour).
     *
     * See tracker #16540
     */

#define RPA(LABEL)      \
    "addl %2, %0\n"						\
    "leal (%0,%2,4), %%edx\n"					\
    "addl %3, %1\n"						\
    "testl %%edx, %%edx\n"					\
    "movl %0, %%eax\n"						\
    "jng " LABEL "\n"						\
    "addl %2, %%eax\n"						\
    "leal (%1,%3,1), %%edx\n"					\
    "cmovngl %%eax, %0\n"						\
    "cmovngl %%edx, %1\n"						\
    "addl %3, %%edx\n"						\
    "addl %2, %%eax\n"						\
    "cmovngl %%edx, %1\n"						\
    "cmovngl %%eax, %0\n"						\
    "addl %3, %%edx\n"						\
    "addl %2, %%eax\n"						\
    "cmovngl %%edx, %1\n"						\
    "cmovngl %%eax, %0\n"
#define RPB(LABEL)      \
    "addl %0, %2\n"						\
    "leal (%2,%0,4), %%edx\n"					\
    "addl %1, %3\n"						\
    "testl %%edx, %%edx\n"					\
    "movl %2, %%eax\n"						\
    "jns " LABEL "\n"						\
    "addl %0, %%eax\n"						\
    "leal (%1,%3,1), %%edx\n"					\
    "cmovnsl %%eax, %2\n"						\
    "cmovnsl %%edx, %3\n"						\
    "addl %1, %%edx\n"						\
    "addl %0, %%eax\n"						\
    "cmovnsl %%edx, %3\n"						\
    "cmovnsl %%eax, %2\n"						\
    "addl %1, %%edx\n"						\
    "addl %0, %%eax\n"						\
    "cmovnsl %%edx, %3\n"						\
    "cmovnsl %%eax, %2\n"
#define RPC     \
    "cltd\n"							\
    "idivl %2\n"							\
    "imull %3, %%eax\n"						\
    "movl %%edx, %0\n"						\
    "subl %%eax, %1\n"
#define RPD "cltd\n"    \
    "idivl %0\n"							\
    "imull %1, %%eax\n"						\
    "movl %%edx, %2\n"						\
    "subl %%eax, %3\n"

    int32_t mhI;
    __asm__ __volatile__ (
            "xorl %1, %1\n"
            "cmpl %2, %5\n"
            "movl $0x1, %3\n"
            "jg 9f\n"
            "movl %5, %%eax\n"
            "negl %%eax\n"
            "movl %%eax, %4\n"
            "movl %0, %%eax\n"
            "cltd\n"
            "idivl %2\n"
            "subl %%eax, %1\n"
            "cmpl $0xe6666667, %%edx\n"
            "movl %%edx, %0\n"
            "jl 0f\n"
            ".balign 8\n"
            "1:\n"
            "cmpl %0, %4\n"
            "jl 9f\n"
            RPB("3f")
            "2:\n"
            "cmpl %2, %5\n"
            "jg 9f\n"
            RPA("4f")
            "jmp 1b\n"
            ".balign 8\n"
            "3:\n"
            RPD
            "jmp 2b\n"
            ".balign 8\n"
            "4:\n"
            RPC
            "jmp 1b\n"
            ".balign 8\n"
            "0:\n"
            "movl %2, %%eax\n"
            "cltd\n"
            "idivl %0\n"
            "imull %1, %%eax\n"
            "subl %%eax, %3\n"
            "cmpl $0x19999999, %%edx\n"
            "movl %%edx, %2\n"
            "jle 2b\n"
            "movl %0, %%eax\n"
            "cltd\n"
            "idivl %2\n"
            "imull %3, %%eax\n"
            "subl %%eax, %1\n"
            "cmpl $0xe6666667, %%edx\n"
            "movl %%edx, %0\n"
            "jge 1b\n"
            "jmp 0b\n"
            "9:\n"
            : "+&r"(i0), "=&r"(j0), "+&r"(i1), "=&r"(j1),
        "=&rm"(mhI) : "rm"(hI) : "%rax", "%rdx", "cc");

#undef RPA
#undef RPB
#undef RPC
#undef RPD
    int64_t k = i1 - hI - i0;
    if (i1 > -i0) {
        if (UNLIKELY(!i0)) return 0;
        k /= i0; i1 -= k * i0; j1 -= k * j0;
    } else {
        if (UNLIKELY(!i1)) return 0;
        k /= i1; i0 += k * i1; j0 += k * j1;
    }
    pli->mi0 = -(int32_t) i0; pli->j0 = (uint32_t) j0; pli->i1 = (int32_t) i1; pli->j1 = (uint32_t) j1;
    return 1;
}
#endif

int main(int argc, char * argv[])
{
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    uint32_t I = 512;
    uint32_t B = 1e6;
    int ntests = 1000;
    unsigned long seed = 0;
    bool failed = false;
    bool quiet = false;
    bool timing = false;

    for( ; argc > 1 ; argv++,argc--) {
        if (strcmp(argv[1], "-B") == 0) {
            argv++,argc--;
            B = atoi(argv[1]);
        } else if (strcmp(argv[1], "-I") == 0) {
            argv++,argc--;
            I = atoi(argv[1]);
        } else if (strcmp(argv[1], "-ntests") == 0) {
            argv++,argc--;
            ntests = atoi(argv[1]);
        } else if (strcmp(argv[1], "-seed") == 0) {
            argv++,argc--;
            seed = atol(argv[1]);
        } else if (strcmp(argv[1], "-q") == 0) {
            quiet = true;
        } else if (strcmp(argv[1], "-T") == 0) {
            timing = true;
        }
    }

    std::vector<std::tuple<unsigned long, int>> prime_powers;
    /* count primes up to B, and add prime powers as well */
    prime_info pi;
    prime_info_init(pi);
    for (unsigned long p = 2 ; p < B ; p = getprime_mt (pi)) {
        int k = 1;
        unsigned long q = p;
        for( ; q < B ; q*=p, k++) {
            if (q < I) continue;
            prime_powers.emplace_back(p, k);
        }
    }
    prime_info_clear(pi);

    std::map<int, std::map<std::string, int> > stats;

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    if (!seed) {
        seed = getpid();
        if (!quiet) printf("seeding with pid = %lu\n", seed);
    }
    gmp_randseed_ui(rstate, seed);

    struct a_test {
        unsigned long p;
        unsigned long q;
        unsigned long r;
        int k;
    };
    std::vector<a_test> tests;

    /* interesting with I=0x20000 */
    tests.emplace_back(a_test { 5, 390625, 771795, 8 });

    for(int i = 0 ; i < ntests ; i++) {
        unsigned long j = gmp_urandomm_ui(rstate, prime_powers.size());
        unsigned long p;
        int k;
        std::tie(p, k) = prime_powers[j];
        unsigned long q = 1;
        for(int s = k ; s-- ; q*=p) ;
        unsigned long r = gmp_urandomm_ui(rstate, q + q / p);
        bool proj = r > q;
        if (proj)
            r = q + p * (r - q);
        tests.emplace_back(a_test { p, q, r, k });
    }

    clock_t clk0 = clock();

    for(auto const & aa : tests) {
        unsigned long p = aa.p;
        unsigned long q = aa.q;
        unsigned long r = aa.r;
        int k = aa.k;
        bool proj = r > q;
        const char * when = "pre";

        // if (proj || k > 1 || r == 0) continue;

        std::string desc;
        if (!quiet) {
            desc = proj ? "proj" : "affine";
            if (p == 2)
                desc += "+even";
        }
        try {
            plattice L(q, proj ? (r - q) : r, proj, I);
            // plattice L; reference(&L, q, r, I);
            // plattice L; reference2(&L, q, r, I);
#ifdef TEST_ASSEMBLY_CODE_DELETED_BY_5f258ce8b
            // plattice L; reference2_asm(&L, q, r, I);
#endif
            when = "post";
            ASSERT_ALWAYS(L.check_post_conditions(I));
            if (!quiet)
                stats[k][desc]++;
        } catch (std::runtime_error const & e) {
            fprintf(stderr, "Failed check (%s) for p^k=%lu^%d r=%lu\n", when, p, k, r);
            failed = true;
        }
    }
    clock_t clk1 = clock();
    if (timing) {
        printf("# %d tests in %.4fs\n",
                ntests, ((double)(clk1-clk0))/CLOCKS_PER_SEC);
    }
    if (!quiet) {
        for(auto & kv : stats) {
            int k = kv.first;
            printf("%d:", k);
            for(auto & dn : kv.second) {
                printf(" %s:%d", dn.first.c_str(), dn.second);
            }
            printf("\n");
        }
    }
    gmp_randclear(rstate);
    return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}

