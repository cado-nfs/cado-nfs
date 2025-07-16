#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cstdlib>

#include <iostream>
#include <sstream>
#include <type_traits>
#include <vector>
#include <stdexcept>

#include <gmp.h>

#include "gmp_aux.h"
#include "cxx_mpz.hpp"
#include "arithxx/modint.hpp"
#include "tests_common.h"
#include "macros.h"

template<typename INTEGER>
struct is_fixed_width : public std::integral_constant<bool, true> {};
template<>
struct is_fixed_width<cxx_mpz> : public std::integral_constant<bool, false> {};

template <typename INTEGER>
class Tests {
    using Integer = INTEGER;

    template<typename T = INTEGER>
    static size_t test_bits()
    requires std::is_same_v<T, cxx_mpz>
    {
        return 16 + gmp_urandomm_ui(state, 1024);
    }

    template<typename T = INTEGER>
    static size_t test_bits()
    requires(!std::is_same_v<T, cxx_mpz>)
    {
        return T::max_bits;
    }

    bool test_get_set_cmp_fits() {
        cxx_mpz z;

        const size_t nbits = test_bits();

        mpz_rrandomb(z, state, nbits);

        const Integer i(z);

        /* cut in 64-bit chunks _manually_, using an algorithm that isn't
         * exactly the copy-paste of the one we have in the code.
         */
        std::vector<uint64_t> limbs;
        for(cxx_mpz zt = z ; zt != 0 ; zt >>= 64)
            limbs.push_back(mpz_get_uint64(zt));

        if (nbits <= 64) {
            /* test construction from just one uint64_t */
            const uint64_t w = mpz_get_uint64(z);
            const Integer i2(w);
            if (i2 != i)
                return false;
            if (!i2.fits_uint64_t())
                return false;
            if (w) {
                if (!(i > (w-1)))
                    return false;
            }
            if (w != UINT64_MAX) {
                if (!(i < (w+1)))
                    return false;
            }
        }

        {
            /* test construction from an array */
            if (Integer(limbs.data(), limbs.size()) != i)
                return false;
        }

        {
            /* test copy */
            if (Integer(i) != i)
                return false;
        }

        {
            /* test comparison. Pay attention to possible overflows. */
            if ((Integer(i) < (i + 1)) != ((i + 1) != 0))
                return false;
            if ((Integer(i) > (i - 1)) != (i != 0))
                return false;
        }

        {
            if (double(i) != mpz_get_d(z))
                return false;
        }

        return true;
    }

    bool test_add_sub() {
        cxx_mpz z0;
        cxx_mpz z1;
        const size_t nbits = test_bits();
        mpz_rrandomb(z0, state, nbits);
        mpz_rrandomb(z1, state, nbits);
        Integer i0(z0);
        Integer i1(z1);

        {
            if (i0.bits() != mpz_sizeinbase(z0, 2))
                return false;
        }

        {
            cxx_mpz z = z0 + z1;
            if (is_fixed_width<Integer>::value) {
                cxx_mpz m = 1;
                mpz_mul_2exp(m, m, nbits);
                if (z > m)
                    z -= m;
            }
            if (i0 + i1 != z)
                return false;
            
            i0 = z0;
            i0 += z1;
            if (i0 + i1 != z)
                return false;
        }

        {
            cxx_mpz z = z0 - z1;
            if (is_fixed_width<Integer>::value) {
                cxx_mpz m = 1;
                mpz_mul_2exp(m, m, nbits);
                if (z < 0)
                    z += m;
            }
            if (i0 + i1 != z)
                return false;
            i0 = z0;
            i0 -= z1;
            if (i0 - i1 != z)
                return false;
        }

        {
            const Integer i = z0 + z0;
            if (i >> 1 != z0)
                return false;
        }

        {
            const Integer i = z0;
            if (i << 1 != z0 + z0)
            return false;
        }

        {
            const Integer i(z0);
            if (i % z1 != (z0 % z1))
                return false;
        }

        {
            ASSERT_ALWAYS(nbits > 2);
            size_t dbits = 1 + gmp_urandomm_ui(state, (unsigned long) (nbits-2));
            mpz_rrandomb(z0, state, dbits);
            mpz_rrandomb(z1, state, nbits - dbits);
            const cxx_mpz z = z0 * z1;
            const Integer i(z);
            if (i.divexact(z0) != z1)
                return false;
        }

        {
            const Integer mz0 = -z0;

            if (mz0 + z0 != 0)
                return false;
        }

        return true;
    }


    bool test_operator_assign_cxx_mpz() {
        Integer i;
        cxx_mpz z;
        
        i = 1;
        i = z;
        if (i != 0) {
            std::cerr << "Setting to 0 failed" << "\n";
            return false;
        }
        z = 5;
        i = z;
        if (i != 5) {
            std::cerr << "Setting to 5 failed" << "\n";
            return false;
        }
        return true;
    }
    
    bool test_stream_operator() {
        cxx_mpz V;
        Integer v;
        bool ok = true;
        for (int i = 0; i < 256; i++) {
            V = 1;
            V <<= i;
            if (!V.fits<Integer>())
                break;
            v = V;
            std::stringstream V_str, v_str;
            V_str << V;
            v_str << v;
            if (V_str.str() != v_str.str()) {
                std::cerr << "GMP output: " << V << ", Integer output: " << v << "\n";
                ok = false;;
            }
        }
        return ok;
    }
public:
    bool runTests() {
        bool ok = true;
        ok &= test_operator_assign_cxx_mpz();
        ok &= test_stream_operator();
        ok &= test_get_set_cmp_fits();
        return ok;
    }
};

int main(int argc, char const * argv[])
{
    bool ok = true;
  
    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);

    try {
        Tests<Integer64> test1;
        ok &= test1.runTests();
        
        Tests<Integer128> test2;
        ok &= test2.runTests();
    } catch (std::runtime_error const & e) {
        std::cerr << e.what();
        ok = false;
    }
      
    tests_common_clear();
    exit(ok ? EXIT_SUCCESS : EXIT_FAILURE);
}
