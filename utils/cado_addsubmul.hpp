#ifndef UTILS_CADO_ADDSUBMUL_HPP_
#define UTILS_CADO_ADDSUBMUL_HPP_

#include "cado_config.h"        // IWYU pragma: keep

#include <gmp.h>

#ifdef HAVE_MPFR
#include <mpfr.h>
#include "cxx_mpfr.hpp"
#include "mpfr_aux.h"
#endif

#ifdef HAVE_MPC
#include <mpc.h>
#include "cxx_mpc.hpp"
#include "mpc_aux.h"
#endif

#include "number_context.hpp"

namespace cado_math_aux {
    /* {{{ addmul/submul. It's similar to fma/fms, but the operand order
     * differs. Also, addmul/submul are compound operations, while fma is
     * not. When working with arbitrary precision types, the working
     * precision is the precision of x.
     *
     * TODO ET: do we really need the three-type template?
     * TODO ET: use eval_type_t ?
     */
    template<typename T, typename U, typename V>
    inline T& addmul(T & x, U const & y, V const & z)
    {
        /* return x += y*z */
        auto const & tr = cado::number_context<T>(x);
        return x = fma(tr(y), z, x);
    }
    template<typename T, typename U, typename V>
    inline T& submul(T & x, U const & y, V const & z)
    {
        /* return x -= y*z */
        auto const & tr = cado::number_context<T>(x);
        return x = -fms(tr(y), z, x);
    }

#ifdef HAVE_MPFR
    inline cxx_mpfr& addmul(cxx_mpfr & x, cxx_mpfr const & y, unsigned long& z)
    {
        mpfr_addmul_ui(x, y, z, MPFR_RNDN);
        return x;
    }
    inline cxx_mpfr& addmul(cxx_mpfr & x, cxx_mpfr const & y, long& z)
    {
        mpfr_addmul_si(x, y, z, MPFR_RNDN);
        return x;
    }
    inline cxx_mpfr& submul(cxx_mpfr & x, cxx_mpfr const & y, unsigned long& z)
    {
        mpfr_submul_ui(x, y, z, MPFR_RNDN);
        return x;
    }
    inline cxx_mpfr& submul(cxx_mpfr & x, cxx_mpfr const & y, long& z)
    {
        mpfr_submul_si(x, y, z, MPFR_RNDN);
        return x;
    }
#endif


#ifdef HAVE_MPC
    template<typename U, typename V>
        inline cxx_mpc& addmul(cxx_mpc & x, U const & y, V const & z)
        {
            auto const & tr = cado::number_context<cxx_mpc>(x);
            cxx_mpc yy = tr(y);
            yy *= z;
            x += yy;
            return x;
        }
    template<typename U, typename V>
        inline cxx_mpc& submul(cxx_mpc & x, U const & y, V const & z)
        {
            auto const & tr = cado::number_context<cxx_mpc>(x);
            cxx_mpc yy = tr(y);
            yy *= z;
            return x -= yy;
        }
    inline cxx_mpc& addmul(cxx_mpc & x, cxx_mpc const & y, cxx_mpc const & z)
    {
        mpc_addmul(x, y, z, MPC_RNDNN);
        return x;
    }
    inline cxx_mpc& addmul(cxx_mpc & x, cxx_mpc const & y, unsigned long& z)
    {
        mpc_addmul_ui(x, y, z, MPC_RNDNN);
        return x;
    }
    inline cxx_mpc& addmul(cxx_mpc & x, cxx_mpc const & y, long& z)
    {
        mpc_addmul_si(x, y, z, MPC_RNDNN);
        return x;
    }
    inline cxx_mpc& submul(cxx_mpc & x, cxx_mpc const & y, cxx_mpc const & z)
    {
        mpc_submul(x, y, z, MPC_RNDNN);
        return x;
    }
    inline cxx_mpc& submul(cxx_mpc & x, cxx_mpc const & y, unsigned long& z)
    {
        mpc_submul_ui(x, y, z, MPC_RNDNN);
        return x;
    }
    inline cxx_mpc& submul(cxx_mpc & x, cxx_mpc const & y, long& z)
    {
        mpc_submul_si(x, y, z, MPC_RNDNN);
        return x;
    }
#endif
    /* }}} */
} /* namespace cado_math_aux */

#endif	/* UTILS_CADO_ADDSUBMUL_HPP_ */
