#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <ext/alloc_traits.h>
#include <cstdlib>

#include <iostream>
#include <sstream> // istringstream // IWYU pragma: keep
#include <vector>
#include <utility>
#include <stdexcept>
#include <string>              // for string, basic_string

#include <gmp.h>               // for mpz_cmp, mpz_set_str, mpz_t

#include "cxx_mpz.hpp"
#include "mpz_poly.h"
#include "mpz_poly_bivariate.hpp"
#include "macros.h"

// coverity[root_function]
static void tests_univariate()
{
    const std::vector<std::pair<std::pair<std::string, std::string>, std::vector<cxx_mpz>>> examples {
        { { "X", "X^128+(X+1)^2+(X^3+X^2+1)*(X+1)-X^2^7" }, {2UL, 3UL, 2UL, 2UL, 1UL}},
        { { "X", "X+1-(X+1)" }, {} },
        { { "z", "z+1" }, {1UL,1UL} },
        { { "a0", "a0+a0^3+1" }, {1UL,1UL,0,1} },
        /* see whether we correctly parse the c85 polynomial */
        { { "x", "960*x^4+85585660*x^3+405578084588*x^2+4213006218645262637*x-2287975327041639106629845" }, { "-2287975327041639106629845"_mpz, "4213006218645262637"_mpz, "405578084588"_mpz, 85585660UL, 960UL}},
    };

    const std::vector<std::pair<std::string, std::string>> expected_failures {
        { "x", "x+y" },
        { "x", "t" },
        { "x", "x+y-y" },
    };

    for(auto const & example : examples) {
        cxx_mpz_poly f;
        if (!(std::istringstream(example.first.second) >> f.named(example.first.first)))
            throw std::runtime_error("cannot parse polynomial\n");
        std::ostringstream os;
        os << f;
        std::cout << os.str() << "\n";
        ASSERT_ALWAYS((size_t)(f->deg+1) == example.second.size());
        for(size_t i = 0 ; i < example.second.size() ; ++i)
            ASSERT_ALWAYS(mpz_cmp(mpz_poly_coeff_const(f, i), example.second[i]) == 0);
        std::istringstream is(os.str());
        decltype(f) g;
        ASSERT_ALWAYS(is >> g && f == g);
    }
    for(auto const & example : expected_failures) {
        cxx_mpz_poly f;
        if ((std::istringstream(example.second) >> f.named(example.first)))
            throw std::runtime_error("unexpected success while parsing bad polynomial\n");
    }
}

static void tests_bivariate()
{
    struct test2 {
        std::string x,y;
        std::string s;
        test2(std::string x, std::string y, std::string s)
        : x(std::move(x))
        , y(std::move(y))
        , s(std::move(s))
        {}
    };

    const std::vector<test2> examples {
        { "x", "y", "x+y" },
        { "x", "t", "(x-t)^3+t+x" },
    };

    /*
        named_proxy(typename std::enable_if<
                    std::is_lvalue_reference<T>::value,
                    named_proxy<typename std::remove_const<T>::type>
                >::type const & c)
                */
    for(auto const & e : examples) {
        cxx_mpz_poly_bivariate f;
        // static_assert(std::is_lvalue_reference<>::value);
        // static_assert(std::is_lvalue_reference<decltype(f.named("x", "y"))>::value);
        // static_assert(std::is_same<decltype(f.named("x", "y")), cxx_mpz_poly_bivariate::named_proxy<cxx_mpz_poly_bivariate &>>::value);
        if (!(std::istringstream(e.s) >> f.named(e.x, e.y))) {
            std::cerr << "cannot parse polynomial\n";
            exit(EXIT_FAILURE);
        }
        std::ostringstream os;
        os << f.named("a0", "a1");
        std::cout << os.str() << "\n";
        std::istringstream is(os.str());
        decltype(f) g;
        ASSERT_ALWAYS(is >> g.named("a0", "a1") && f == g);
    }
}


// coverity[root_function]
int main()
{
    tests_univariate();
    tests_bivariate();
}
