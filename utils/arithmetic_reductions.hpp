#ifndef UTILS_ARITHMETIC_REDUCTIONS_HPP_
#define UTILS_ARITHMETIC_REDUCTIONS_HPP_

#include <string>
#include <tuple>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "mpz_poly.h"
#include "mpz_poly_bivariate.hpp"
#include "matrix.hpp"
#include "utils_cxx.hpp"

/* our goal here is to factor out a few reduction operators that can
 * makes sense both for univariate and bivariate polynomials.
 */
namespace cado::arithmetic_reductions {

    struct noop { /*{{{*/
        template<typename T>
        void operator()(T & a, T const & b) const {
            a = b;
        }
        std::string print() const;
    }; /*}}}*/

    struct mod_mpz { /*{{{*/
        /* This reduction operator is mod an odd integer, but not
         * necessarily prime. Attention: we hold a pointer to the
         * corresponding prime!
         *
         * Coefficient-wise reduction is done either on univariate or
         * bivariate polynomials.
         */
        mpz_srcptr p;

        public:
        explicit mod_mpz(mpz_srcptr p)
            : p(p)
        {
        }
        // mod_mpz(mod_mpz const &) = default;
        void operator()(cado::matrix<cxx_mpz> & b, cado::matrix<cxx_mpz> const & a) const;
        void operator()(cado::matrix<cxx_mpz_poly> & b, cado::matrix<cxx_mpz_poly> const & a) const;
        void operator()(cxx_mpz_poly & b, cxx_mpz_poly const & a) const;
        void operator()(cxx_mpz_poly_bivariate & b, cxx_mpz_poly_bivariate const & a) const;
        std::string print() const;
        // {{{ barrett pre-inverse
        private:
        cado::cached_property<cxx_mpz> cached_invp;
        cxx_mpz compute_invp() const {
            cxx_mpz a;
            ::barrett_precompute_inverse(a, p);
            return a;
        }
        protected:
        bool invp_is_cached() const { return cached_invp.is_cached(); }
        mpz_srcptr invp_lazy() const {
            if (cached_invp.is_cached())
                return invp();
            else
                return nullptr;
        }
        cxx_mpz const & invp() const {
            return cached_invp([this](){ return compute_invp();});
        }
        // }}}
    }; /*}}}*/

    struct mod_p : public mod_mpz {
        /* This is really the same as mod_mpz, but we prefer to let the
         * name reflect the fact that here, we asume that the modulus is
         * prime.
         */
        using mod_mpz::mod_mpz;
    };

    struct mod_fx { /*{{{*/
        /* Think arithmetic in a relative number field, for example.
         * There aren't many uses of this reducer at the moment.
         */
        cxx_mpz_poly const & fx;
        mod_fx(cxx_mpz_poly const & fx)
            : fx(fx)
        {
            ASSERT_ALWAYS(mpz_poly_is_monic(fx));
        }
        void operator()(cxx_mpz_poly & B, cxx_mpz_poly const & A) const;
        void operator()(cxx_mpz_poly_bivariate & B, cxx_mpz_poly_bivariate const & A) const;
        std::string print() const;
    }; /*}}}*/

    struct mod_fx_mod_p : public mod_p { /*{{{*/
        /* (p, fx) defines an ideal of Z[x], which may or may not be
         * prime. Note that p _must_ be prime, however.
         */
        cxx_mpz_poly const & fx;
        mod_fx_mod_p(mod_p const & R,
                               cxx_mpz_poly const & fx)
            : mod_p(R)
            , fx(fx)
        {
            /* trigger the computation of the barrett inverse, as it
             * probably never makes sense to skip it */
            mod_p::invp();
        }
        mod_fx_mod_p(cxx_mpz_poly const & fx, mpz_srcptr p)
            : mod_p(p)
            , fx(fx)
        {
            /* trigger the computation of the barrett inverse, as it
             * probably never makes sense to skip it */
            mod_p::invp();
        }
        void operator()(cado::matrix<cxx_mpz_poly> & b, cado::matrix<cxx_mpz_poly> const & a) const;
        void operator()(cxx_mpz_poly & B, cxx_mpz_poly const & A) const;
        void operator()(cxx_mpz_poly_bivariate & B, cxx_mpz_poly_bivariate const & A) const;
        std::string print() const;
        bool make_monic(cxx_mpz_poly_bivariate & f) const;
    };/*}}}*/

    struct mod_q : public mod_fx_mod_p {
        /* _this_ is really modulo a prime ideal
         * In effect, when (p, fx) define a prime ideal, this implements
         * arithmetic of polynomials with coefficients in a finite field.
         */
        mod_q(mod_p const & R, cxx_mpz_poly const & fx)
            : mod_fx_mod_p(R, fx)
        {}
        mod_q(cxx_mpz_poly const & fx, mpz_srcptr p)
            : mod_fx_mod_p(fx, p)
        {}

        private:
        cado::matrix<cxx_mpz> compute_frobenius_matrix() const;
        cado::cached_property<cado::matrix<cxx_mpz>> cached_fm;
        protected:
        cado::matrix<cxx_mpz> const & frobenius_matrix() const {
            return cached_fm([this](){ return compute_frobenius_matrix();});
        }

        private:
        cado::matrix<cxx_mpz> compute_inverse_frobenius_matrix() const;
        cado::cached_property<cado::matrix<cxx_mpz>> cached_ifm;
        protected:
        cado::matrix<cxx_mpz> const & inverse_frobenius_matrix() const {
            return cached_ifm([this](){ return compute_inverse_frobenius_matrix();});
        }

        protected:
        void matrix_apply(cxx_mpz_poly & B, cxx_mpz_poly const & A, cado::matrix<cxx_mpz> const & M) const;
        public:
        void frobenius(cxx_mpz_poly & B, cxx_mpz_poly const & A, int i) const;
        void frobenius(cxx_mpz_poly_bivariate & B, cxx_mpz_poly_bivariate const & A, int i) const;

        friend struct mod_fy_mod_q;
    };/*}}}*/

    struct mod_fy_mod_q : public mod_q {/*{{{*/
        public:
        private:
        cado::matrix<cxx_mpz_poly> compute_small_frobenius_matrix() const;
        cado::cached_property<cado::matrix<cxx_mpz_poly>> cached_sfm;
        protected:
        cado::matrix<cxx_mpz_poly> const & small_frobenius_matrix() const {
            return cached_sfm([this](){ return compute_small_frobenius_matrix();});
        }

        /*
        private:
        cado::matrix<cxx_mpz_poly> compute_inverse_small_frobenius_matrix() const;
        cado::cached_property<cado::matrix<cxx_mpz_poly>> cached_isfm;
        protected:
        cado::matrix<cxx_mpz_poly> const & inverse_small_frobenius_matrix() const {
            return cached_isfm([this](){ return compute_inverse_small_frobenius_matrix();});
        }
        public:
        */

        private:
        auto matrix_galois_norm(int k) const
            -> std::tuple<cado::matrix<cxx_mpz_poly>, cado::matrix<cxx_mpz>>;
        cado::matrix<cxx_mpz_poly> compute_frobenius_matrix() const;
        cado::cached_property<cado::matrix<cxx_mpz_poly>> cached_fm;
        protected:
        cado::matrix<cxx_mpz_poly> const & frobenius_matrix() const {
            return cached_fm([this](){ return compute_frobenius_matrix();});
        }

        /*
        private:
        cado::matrix<cxx_mpz_poly> compute_inverse_frobenius_matrix() const;
        cado::cached_property<cado::matrix<cxx_mpz_poly>> cached_ifm;
        protected:
        cado::matrix<cxx_mpz_poly> const & inverse_frobenius_matrix() const {
            return cached_ifm([this](){ return compute_inverse_frobenius_matrix();});
        }
        */
        public:

        void frobenius(cxx_mpz_poly_bivariate & B, cxx_mpz_poly_bivariate const & A, int i) const;

        /* In effect, when (p, fx) define a prime ideal, this implements
         * arithmetic of polynomials with coefficients in a finite field
         */
        cxx_mpz_poly_bivariate const & fy;
        mod_fy_mod_q(mod_q const & R,
                                      cxx_mpz_poly_bivariate const & fy)
            : mod_q(R)
            , fy(fy)
        {
        }
        mod_fy_mod_q(cxx_mpz_poly_bivariate const & fy, cxx_mpz_poly const & fx,
                                      mpz_srcptr p)
            : mod_q(fx, p)
            , fy(fy)
        {
        }
        void operator()(cxx_mpz_poly_bivariate & B, cxx_mpz_poly_bivariate const & A) const;
        void operator()(cado::matrix<cxx_mpz_poly> & b, cado::matrix<cxx_mpz_poly> const & a) const
        {
            ((mod_q const &)(*this))(b, a);
        }
        void operator()(cxx_mpz_poly & B, cxx_mpz_poly const & A) const
        {
            ((mod_q const &)(*this))(B, A);
        }
        std::string print() const;
    }; /*}}}*/

} /* namespace cado::arithmetic_reductions */

#endif	/* UTILS_ARITHMETIC_REDUCTIONS_HPP_ */
