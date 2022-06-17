#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdint.h>        // for uint64_t
#include <gmp.h>
#include "gmp-hacks.h"     // for MPN_SET_MPZ, MPZ_SET_MPN
#include "tests_common.h"
#include "arith-modp.hpp"
#include "timing.h"
#include "macros.h"

using namespace arith_modp;

template<typename F>
void do_tests(unsigned long iter, int summands, int cbound )
{
    mp_size_t m = F::nlimbs * mp_bits_per_limb;
    // mp_size_t ell = F::extra*mp_bits_per_limb;

    typedef typename F::elt_ur_for_add a_type;
    constexpr const int a_extra = a_type::nlimbs - F::nlimbs;
    typedef typename F::elt_ur_for_addmul m_type;
    constexpr const int m_extra = m_type::nlimbs - F::nlimbs;

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
        } while (mpz_size(pz) < F::nlimbs || mpz_even_p(pz));
        F field(pz, 1);

        auto x = field.alloc(summands);
        auto y = field.alloc(summands);
        auto r = field.alloc();

        cxx_mpz xz[summands];
        cxx_mpz yz[summands];


        for(int i = 0 ; i < summands ; i++) {
            mpz_urandomm(xz[i], state, pz);
            x[i] = xz[i];
            mpz_urandomm(yz[i], state, pz);
            y[i] = yz[i];
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
                a_type a;
                a.zero();
                for(int i = 0 ; i < summands ; i++) {
                    if (coeffs[i] == 1) {
                        field.add(a, x[i]);
                    } else if (coeffs[i] == -1) {
                        field.sub(a, x[i]);
                    } else if (coeffs[i] > 0) {
                        field.addmul_ui(a, x[i], coeffs[i]);
                    } else {
                        field.submul_ui(a, x[i], -coeffs[i]);
                    }
                }
                field.reduce(*r, a);
                tt += microseconds();
            }
            MPZ_SET_MPN(rz, r->x, F::nlimbs);

            ASSERT_ALWAYS(mpz_cmp(az, rz) == 0);

            if (tests_common_get_verbose()) {
                gmp_printf("ok [%d+%d] %Zd\n", F::nlimbs, a_extra, (mpz_srcptr) pz);
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
                m_type m;
                m.zero();
                for(int i = 0 ; i < summands ; i++) {
                    field.addmul_ur(m, x[i], y[i]);
                }
                field.reduce(*r, m);
                tt += microseconds();
            }
            MPZ_SET_MPN(rz, r->x, F::nlimbs);

            ASSERT_ALWAYS(mpz_cmp(mz, rz) == 0);

            if (tests_common_get_verbose()) {
                gmp_printf("ok [%d+%d] %Zd\n", F::nlimbs, m_extra, (mpz_srcptr) pz);
            }
        }

        field.free(x);
        field.free(y);
        field.free(r);
    }
    printf("%lu tests ok [%d] in %.4f s\n", iter, F::nlimbs, 1.0e-6 * tt);
}

// coverity[root_function]
int main(int argc, const char * argv[])
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
    tests_common_clear();
}

