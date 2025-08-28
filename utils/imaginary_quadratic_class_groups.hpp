#ifndef CADO_IMAGINARY_QUADRATIC_CLASS_GROUPS_HPP
#define CADO_IMAGINARY_QUADRATIC_CLASS_GROUPS_HPP

#include <exception>
#include <iostream>

#include "fmt/base.h"

#include "cado_poly.h"
#include "cxx_mpz.hpp"

bool cado_poly_is_imaginary_quadratic(cxx_cado_poly const & cpoly);

class imaginary_quadratic_form; /* forward declaration */

class imaginary_quadratic_class_group
{
    cxx_mpz disc;

    public:

    explicit imaginary_quadratic_class_group(cxx_mpz const & d);

    cxx_mpz const & discriminant() const
    {
        return disc;
    }

    bool operator==(imaginary_quadratic_class_group const & o) const
    {
        return mpz_cmp(disc, o.disc) == 0;
    }

    imaginary_quadratic_form operator()(cxx_mpz const & a, cxx_mpz const & b,
                                        cxx_mpz const & c) const;

    unsigned long bound_generating_set_grh() const;

    imaginary_quadratic_form one() const;

    friend std::ostream & operator<<(std::ostream & o,
                                     imaginary_quadratic_class_group const & G);
};

/*
 * A form has always a positive coefficient a and is always primitive
 * (gcd(a,b,c) = 1) and reduced (-a < b <= a and and a < c or a = c and b >= 0).
 * Its discriminant is equal to cl.discriminant() and thus always negative.
 */
class imaginary_quadratic_form
{
    imaginary_quadratic_class_group const & cl;
    cxx_mpz a, b, c;

    void set_c_from_disc();

    void throw_if_a_is_nonpositive() const;
    void throw_if_not_primitive() const;
    void throw_if_wrong_discriminant() const;

    /* q and r are scratch variables */
    void normalize(cxx_mpz & q, cxx_mpz & r);
    void rho(cxx_mpz & q, cxx_mpz & r);

    /*
     * Reduce the qfi f.
     * The form f is reduced if the form is
     *      normal (-a < b <= a)
     *    and
     *      a < c   or    a = c and b >= 0.
     *
     * t0 and t1 are scratch variables.
     */
    void reduction(cxx_mpz & t0, cxx_mpz & t1);
    void reduction();

    struct ScratchVars {
        cxx_mpz g, d, u0, v0, v1, v, w, s, t, a2d;
    };

    /* Assumes r, f1 and f2 belong to the same class group */
    static void compose_mul(imaginary_quadratic_form & r,
                            imaginary_quadratic_form const & f1,
                            imaginary_quadratic_form const & f2,
                            ScratchVars & tmp);

    /* Assumes r and f belong to the same class group */
    static void compose_sqr(imaginary_quadratic_form & r,
                            imaginary_quadratic_form const & f,
                            ScratchVars & tmp);

    public:

    struct not_primitive: public std::exception {
        const char * what() const noexcept override { return "not primitive"; }
    };

    imaginary_quadratic_form() = delete;
    imaginary_quadratic_form(imaginary_quadratic_class_group const & cl,
                             cxx_mpz const & a,
                             cxx_mpz const & b);
    imaginary_quadratic_form(imaginary_quadratic_class_group const & cl,
                             cxx_mpz const & a,
                             cxx_mpz const & b,
                             cxx_mpz const & c);

    imaginary_quadratic_form(imaginary_quadratic_form const &) = default;
    imaginary_quadratic_form & operator=(imaginary_quadratic_form const &);

    imaginary_quadratic_form(imaginary_quadratic_form &&);
    imaginary_quadratic_form & operator=(imaginary_quadratic_form &&);

    bool is_one() const
    {
        return a == 1U; /* forms are always reduced */
    }

    bool operator==(imaginary_quadratic_form const & o) const
    {
        return a == o.a && b == o.b && c == o.c;
    }

    bool operator<(imaginary_quadratic_form const & o) const
    {
        int r = mpz_cmp(a, o.a);
        if (r == 0) {
            r = mpz_cmp(b, o.b);
            if (r == 0) {
                r = mpz_cmp(c, o.c);
            }
        }
        return r < 0;
    }

    imaginary_quadratic_form operator*(imaginary_quadratic_form const & f) const;
    imaginary_quadratic_form operator^(cxx_mpz const & e) const;

    friend std::ostream & operator<<(std::ostream & o,
                                     imaginary_quadratic_form const & f);
    friend struct std::hash<imaginary_quadratic_form>;
};

namespace fmt {
    template <> struct formatter<imaginary_quadratic_class_group> : ostream_formatter {};
    template <> struct formatter<imaginary_quadratic_form> : ostream_formatter {};
}

#endif /* CADO_IMAGINARY_QUADRATIC_CLASS_GROUPS_HPP */
