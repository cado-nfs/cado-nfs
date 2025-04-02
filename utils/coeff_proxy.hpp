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
        unsigned int i;
        // NOLINTNEXTLINE(hicpp-explicit-conversions)
        operator T() { return (i < p.coeffs.size()) ? p.coeffs[i] : 0; }
        coeff_proxy& operator=(T x) {
            if (i + 1 < p.coeffs.size()) {
                p.coeffs[i] = x;
            } else {
                p.coeffs.insert(p.coeffs.end(), (i + 1 - p.size()), T(0));
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
        unsigned int i;
        // NOLINTNEXTLINE(hicpp-explicit-conversions)
        operator T() { return (i < p.size()) ? p.coeffs[i] : 0; }
    };
}


#endif	/* CADO_UTILS_COEFF_PROXY_HPP */
