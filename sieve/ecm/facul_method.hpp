#ifndef FACUL_METHOD_HPP_
#define FACUL_METHOD_HPP_

#include <algorithm>
#include <tuple>
#include <utility>
#include <vector>
#include "facul_ecm.h"  // for ecm_parameterization_t

#define PM1_METHOD 1
#define PP1_27_METHOD 2
#define PP1_65_METHOD 3
#define EC_METHOD 4
#define MPQS_METHOD 5

struct facul_method {
    long method = 0;    /* Which method to use (P-1, P+1 or ECM) */
    void *plan = NULL;  /* Parameters for that method */

    struct parameters {
        int method = 0;
        unsigned long B1 = 0;
        unsigned long B2 = 0;

        // only used for ECM
        // we have default values, only to silence code analyzer
        // warnings.
        ec_parameterization_t parameterization = BRENT12;
        unsigned long parameter = 0;
        int extra_primes = 0;

        /* we had traces of a notion of null methods in the source
         * files, which were ignored. I don't know what this is
         * supposed to be good for.  */
        bool is_null() const { return B1 == 0 && B2 == 0; }

        /* Note that parsing for this is (for the moment) in
         * facul_strategies.cpp */
        parameters() = default;
        ~parameters() = default;
        parameters(parameters const&) = default;
        parameters(parameters &&) = default;
        parameters& operator=(parameters const&) = default;
        parameters& operator=(parameters &&) = default;

        parameters(int method, unsigned long B1, unsigned long B2)
            : method(method), B1(B1), B2(B2)
        {
            ASSERT_ALWAYS(method != EC_METHOD);
        }

        parameters(int method, unsigned long B1, unsigned long B2, ec_parameterization_t parameterization, unsigned long parameter, int extra_primes = 1)
            : method(method), B1(B1), B2(B2),
            parameterization(parameterization),
            parameter(parameter),
            extra_primes(extra_primes)
        {
            ASSERT_ALWAYS(method == EC_METHOD);
        }

        bool operator<(parameters const & a) const {
            if (method != EC_METHOD || a.method != EC_METHOD)
                return
                    std::forward_as_tuple(method, B1, B2)
                    <
                    std::forward_as_tuple(a.method, a.B1, a.B2);
            else
                return
                    std::forward_as_tuple(method, B1, B2, parameterization, parameter, extra_primes)
                    <
                    std::forward_as_tuple(a.method, a.B1, a.B2, a.parameterization, a.parameter, a.extra_primes);
        }
    };

    struct parameters_with_side : public parameters {
        int side;
        parameters_with_side() = default;
        parameters_with_side(parameters_with_side const &) = default;
        parameters_with_side& operator=(parameters_with_side const &) = default;
        parameters_with_side(parameters_with_side&&) = default;
        parameters_with_side& operator=(parameters_with_side&&) = default;
        /*
        parameters_with_side(parameters const & p, int side)
            : parameters(p), side(side)
        {}
        */
        template<typename... Args>
        parameters_with_side(int side, Args&&... args)
        : parameters(std::forward<Args>(args)...), side(side)
        {}
        bool operator<(parameters_with_side const & a) const {
                return
                    std::forward_as_tuple(side, (parameters const&) *this)
                    <
                    std::forward_as_tuple(a.side, (parameters const&) a);
        }
    };

    facul_method() = default;
    facul_method(facul_method const &) = delete;
    facul_method& operator=(facul_method const &) = delete;
    facul_method(facul_method&& o) {
        plan = o.plan;
        method = o.method;
        o.plan = NULL;
        o.method = 0;
    }
    facul_method& operator=(facul_method&& o) {
        std::swap(plan, o.plan);
        std::swap(method, o.method);
        return *this;
    }
    facul_method(parameters const &, const int verbose = 0);
    ~facul_method();
};

struct facul_method_side {
    /* This points to just _one_ method, which is in the book of
     * precomputed method. It's _not_ an array. */
    facul_method const * method;

    int side;           /* To know on which side this method will be applied */
    int is_last;        /* To know if this method is the last on its side
                           (used when you chain methods).  */
    static void fix_is_last(std::vector<facul_method_side>& v);

    facul_method_side(facul_method const * method, int side)
        : method(method), side(side), is_last(0)
    { }
};

#endif	/* FACUL_METHOD_HPP_ */
