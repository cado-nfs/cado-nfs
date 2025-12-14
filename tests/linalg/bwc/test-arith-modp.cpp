#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdint>        // for uint64_t

#include <type_traits>

#include <gmp.h>

#include "arith-modp.hpp"
#include "cxx_mpz.hpp"
#include "gmp-hacks.h"     // for MPN_SET_MPZ, MPZ_SET_MPN
#include "gmp_aux.h"
#include "macros.h"
#include "tests_common.h"
#include "timing.h"

using namespace arith_modp;

template<typename F>
static void do_tests(unsigned long iter, int summands, int cbound )
{

    size_t const maximum_limbs = std::conditional<F::is_constant_width, std::integral_constant<size_t, F::constant_width>, std::integral_constant<int, 32>>::type::value;

    mp_size_t const m = maximum_limbs * mp_bits_per_limb;
    // mp_size_t ell = F::extra*mp_bits_per_limb;

    uint64_t tt = 0;

    for(unsigned long test = 0 ; test < iter ; test++) {
        int coeffs[summands];

        for(int i = 0 ; i < summands ; i++) {
            /* cheat a little bit so that we effectively have 90% of
             * +1/-1
             */
            coeffs[i] = gmp_urandomm_ui(state, 2 * 10*cbound);
            coeffs[i] -= 10*cbound - (coeffs[i] >= 10*cbound);
            if (coeffs[i] > cbound) coeffs[i]=1;
            if (coeffs[i] < -cbound) coeffs[i]=-1;
        }

        cxx_mpz pz;
        /* Take a random odd modulus, and fits
         * within m bits. We insist on having the right number of words,
         * but we're happy if the number of bits varies, which is why we
         * won't content ourselves with just the output of mpz_rrandomb
         * */
        do {
            mpz_rrandomb(pz, state, m);
            /* Support p not having its high bit set! */
            mpz_fdiv_q_2exp(pz, pz, gmp_urandomm_ui(state, GMP_LIMB_BITS));
        } while (mpz_size(pz) < maximum_limbs || mpz_even_p(pz) || mpz_cmp_ui(pz, 1) == 0);
        F field(pz, 1U);

        size_t const N = field.template nlimbs<typename F::elt>();

        using a_type = typename F::elt_ur_for_add;
        using m_type = typename F::elt_ur_for_addmul;
        const int a_extra = field.template overhead_limbs<a_type>();
        const int m_extra = field.template overhead_limbs<m_type>();

        auto x = field.alloc(summands);
        auto y = field.alloc(summands);
        ARITH_MODP_TEMPORARY_ALLOC(&field, elt, r);

        cxx_mpz xz[summands];
        cxx_mpz yz[summands];


        for(int i = 0 ; i < summands ; i++) {
            mpz_urandomm(xz[i], state, pz);
            field.set(field.vec_item(x, i), xz[i]);
            mpz_urandomm(yz[i], state, pz);
            field.set(field.vec_item(y, i), yz[i]);
        }


        /** additions **/
        {
            cxx_mpz az, rz;
            mpz_set_ui(az, 0);

            /* Do the computation with gmp */
            for(int i = 0 ; i < summands ; i++)
                mpz_addmul_si(az, xz[i], coeffs[i]);

            mpz_mod(az, az, pz);

            {
                tt -= microseconds();
                ARITH_MODP_TEMPORARY_ALLOC(&field, elt_ur_for_add, a);
                field.set_zero(a);
                for(int i = 0 ; i < summands ; i++) {
                    if (coeffs[i] == 1) {
                        field.add(a, field.vec_item(x, i));
                    } else if (coeffs[i] == -1) {
                        field.sub(a, field.vec_item(x, i));
                    } else if (coeffs[i] > 0) {
                        field.addmul_ui(a, field.vec_item(x, i), coeffs[i]);
                    } else {
                        field.submul_ui(a, field.vec_item(x, i), -coeffs[i]);
                    }
                }
                field.reduce(r, a);
                tt += microseconds();
            }
            MPZ_SET_MPN(rz, r.pointer(), N);

            ASSERT_ALWAYS(mpz_cmp(az, rz) == 0);

            if (tests_common_get_verbose()) {
                gmp_printf("ok [%d+%d] %Zd\n", N, a_extra, (mpz_srcptr) pz);
            }
        }
        /** multiplications **/
        {
            cxx_mpz mz, rz;
            mpz_set_ui(mz, 0);
            /* Do the computation with gmp */
            for(int i = 0 ; i < summands ; i++) {
                mpz_addmul(mz, xz[i], yz[i]);
            }
            mpz_mod(mz, mz, pz);

            {
                tt -= microseconds();
                ARITH_MODP_TEMPORARY_ALLOC(&field, elt_ur_for_addmul, m);
                field.set_zero(m);
                for(int i = 0 ; i < summands ; i++) {
                    field.addmul_ur(m, field.vec_item(x, i), field.vec_item(y, i));
                }
                field.reduce(r, m);
                tt += microseconds();
            }
            MPZ_SET_MPN(rz, r.pointer(), N);

            ASSERT_ALWAYS(mpz_cmp(mz, rz) == 0);

            if (tests_common_get_verbose()) {
                gmp_printf("ok [%d+%d] %Zd\n", N, m_extra, (mpz_srcptr) pz);
            }
        }

        field.free(x);
        field.free(y);
    }
    printf("%lu tests ok [%zu] in %.4f s\n", iter, maximum_limbs, 1.0e-6 * tt);
}

// coverity[root_function]
int main(int argc, char const * argv[])
{
    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER | PARSE_VERBOSE);

    unsigned long iter = 100;
    tests_common_get_iter(&iter);

#define SUMMANDS 200
#define CBOUND 1000
    do_tests< gfp<1> >(iter, SUMMANDS, CBOUND);
    do_tests< gfp<2> >(iter, SUMMANDS, CBOUND);
    do_tests< gfp<3> >(iter, SUMMANDS, CBOUND);
    do_tests< gfp<4> >(iter, SUMMANDS, CBOUND);
    do_tests< gfp<5> >(iter, SUMMANDS, CBOUND);
    do_tests< gfp<6> >(iter, SUMMANDS, CBOUND);
    do_tests< gfp<7> >(iter, SUMMANDS, CBOUND);
    do_tests< gfp<8> >(iter, SUMMANDS, CBOUND);
    do_tests< gfp<9> >(iter, SUMMANDS, CBOUND);
    do_tests< gfp<0> >(iter, SUMMANDS, CBOUND);
    tests_common_clear();
}

