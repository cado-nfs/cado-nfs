#include "cado.h" // IWYU pragma: keep

#include <cstdlib>

#include <iostream>
#include <string>
#include <memory>
#include <fstream>

#include <gmp.h>
#include "fmt/base.h"

#include "arithxx/mod64.hpp"
#include "cxx_mpz.hpp"
#include "macros.h"

template <class layer>
class Lucas {
public:
    using Modulus = typename layer::Modulus;
    using Integer = typename layer::Integer;
    using Residue = typename layer::Residue;
    using ResidueOp = typename layer::Modulus::ResidueOp;

    /* This does one Lucas strong primality test on N, with starting
     * point b
     */
    static bool test_one_lucas(Modulus const & m, Residue const & P)
    {
        Integer ell = m.getmod() + 1;
        int k = 0;
        for( ; (ell & 1) == 0 ; k++, ell>>=1) ;
        ASSERT_ALWAYS(k);
        Residue v(m);
        m.V(v, nullptr, P, ell);
        /* if v is already -2 or +2, things are not going to change */
        Residue const two = m(2);
        Residue mtwo(m);
        m.neg(mtwo, two);
        if (m.equal(v, two) || m.equal(v, mtwo))
            return true;
        for(int i = 0 ; i < k ; i++) {
            if (m.is0(v))
                return true;
            /* invariant: v = v_n = alpha^n + beta^n is neither 2 nor -2,
             */
            m.V_dbl(v, v, two); // v_{2n} = v_n^2-2
            /* Here, if N is prime then v_{2n} can't be +2 because that would
             * mean that v_n^2=4 and thus v_n == +2 or -2, which should both
             * have been checked earlier. So v_{2n} == 2 is a sign of
             * compositeness. Likewise, v_{2n}=-2 can only happen if
             * v_n^2=0. Since we checked v_n==0 before, that would mean
             * that n is not squarefree, and thus not prime.
             */
            if (m.equal(v, two) || m.equal(v, mtwo))
                return false;
        }
        return m.equal(v, two);
    }

    static bool my_test(cxx_mpz const & N)
    {
        if (N == 2) return true;
        if (!(N & 1)) return false;
        Modulus const m { Integer(N) };
        Residue P = m(3);
        Residue D = m(5);
        for(int i = 0 ; m.jacobi(D) != -1 ; i++) {
            m.add1(P, P);
            m.add(D, D, P);
            m.add(D, D, P);
            m.add(D, D, P);
            m.add(D, D, P);
            m.add1(P, P);
            if (i == 20) {
                /* N might be a square, in which case we're going to loop
                 * forever here. */
                if (mpz_perfect_square_p(N))
                    return false;
                // return mpz_probab_prime_p(N, 2);
            }
        }
        return test_one_lucas(m, P);
    }
};


int main(int argc, char const * argv[])
{
    std::istream * input = &std::cin;
    std::unique_ptr<std::istream> file_input;
    if (argc == 2) {
        auto * ptr = new std::ifstream(argv[1]);
        if (!ptr)
            throw std::runtime_error(fmt::format("File not found: {}", argv[1]));
        file_input.reset(ptr);
        input = ptr;
    }
    for(std::string s; std::getline(*input, s); ) {
        cxx_mpz z;
        mpz_set_str(z, s.c_str(), 0);
        if (Lucas<arithxx_mod64>::my_test(z))
            fmt::print("{}\n", z);
    }
}

