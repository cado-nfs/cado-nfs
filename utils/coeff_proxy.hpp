#ifndef CADO_UTILS_COEFF_PROXY_HPP
#define CADO_UTILS_COEFF_PROXY_HPP

/* This is typically used in polynomial types, at least those that have:
 *      typedef P::coefficient_type
 *      P::coeffs
 *      P::cleandeg()
 *      P::degree()
 */

namespace cado_details {
    template<typename P>
    struct coeff_proxy {
        P & p;
        typedef typename P::coefficient_type T;
        int i;
        // NOLINTNEXTLINE(hicpp-explicit-conversions)
        operator T() { return (i <= p.degree()) ? p.coeffs[i] : 0; }
        coeff_proxy& operator=(T x) {
            if (i < p.degree()) {
                p.coeffs[i] = x;
            } else {
                p.coeffs.reserve(i + 1);
                for(int j = p.degree() + 1 ; j <= i ; j++)
                    p.coeffs.push_back(0);
                p.coeffs[i] = x;
                p.cleandeg(p.degree());
            }
            return *this;
        }
        coeff_proxy& operator+=(T x) { return (*this)=(T(*this)+x); }
        coeff_proxy& operator-=(T x) { return (*this)=(T(*this)-x); }
        coeff_proxy& operator*=(T x) { return (*this)=(T(*this)*x); }
        coeff_proxy& operator/=(T x) { return (*this)=(T(*this)/x); }
    };
    template<typename P>
    struct const_coeff_proxy {
        P const & p;
        typedef typename P::coefficient_type T;
        int i;
        // NOLINTNEXTLINE(hicpp-explicit-conversions)
        operator T() { return (i <= p.degree()) ? p.coeffs[i] : 0; }
    };
}


#endif	/* CADO_UTILS_COEFF_PROXY_HPP */
