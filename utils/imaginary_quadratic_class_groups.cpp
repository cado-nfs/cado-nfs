#include "cado.h" // IWYU pragma: keep

#include "gmp_aux.h" // mpz_log
#include "imaginary_quadratic_class_groups.hpp"

/* check that cpoly has only one side of the form x^2-d with d < 0 */
bool cado_poly_is_imaginary_quadratic(cxx_cado_poly const & cpoly)
{
    return cpoly->nb_polys == 1
            && cpoly->pols[0]->deg == 2
            && mpz_poly_is_monic(cpoly->pols[0])
            && mpz_sgn(mpz_poly_coeff_const(cpoly->pols[0], 1)) == 0
            && mpz_sgn(mpz_poly_coeff_const(cpoly->pols[0], 0)) > 0;
}


imaginary_quadratic_class_group::imaginary_quadratic_class_group(
        cxx_mpz const & d)
    : disc(d)
{
    ASSERT_ALWAYS(mpz_sgn(disc) < 0);
    ASSERT_ALWAYS(mpz_fdiv_ui(disc, 4) == 0 || mpz_fdiv_ui(disc, 4) == 1);
}

imaginary_quadratic_form
imaginary_quadratic_class_group::operator()(cxx_mpz const & a,
                                            cxx_mpz const & b,
                                            cxx_mpz const & c) const
{
    return { *this, a, b, c };
}

unsigned long imaginary_quadratic_class_group::bound_generating_set_grh() const
{
    double lg = mpz_log(disc);
    return ceil(6.*lg*lg);

}

imaginary_quadratic_form imaginary_quadratic_class_group::one() const
{
    return { *this, 1U, (mpz_odd_p(disc) ? 1U : 0U) };
}

std::ostream & operator<<(std::ostream & o,
                          imaginary_quadratic_class_group const & G)
{
    return o << "Cl(" << G.disc << ")";
}


imaginary_quadratic_form::imaginary_quadratic_form(
        imaginary_quadratic_class_group const &cl,
        cxx_mpz const & a,
        cxx_mpz const & b)
    : cl(cl)
    , a(a)
    , b(b)
{
    throw_if_a_is_nonpositive();
    set_c_from_disc();
    throw_if_not_primitive();
    reduction();
}

imaginary_quadratic_form::imaginary_quadratic_form(
        imaginary_quadratic_class_group const &cl,
        cxx_mpz const & a,
        cxx_mpz const & b,
        cxx_mpz const & c)
    : cl(cl)
    , a(a)
    , b(b)
    , c(c)
{
    throw_if_a_is_nonpositive();
    throw_if_not_primitive();
    throw_if_wrong_discriminant();
    reduction();
}

imaginary_quadratic_form & imaginary_quadratic_form::operator=(
        imaginary_quadratic_form const & o)
{
        ASSERT_ALWAYS(cl == o.cl);
        a = o.a;
        b = o.b;
        c = o.c;
        return *this;
}

imaginary_quadratic_form::imaginary_quadratic_form(
        imaginary_quadratic_form && o)
    : cl(o.cl), a(std::move(o.a)), b(std::move(o.b)), c(std::move(o.c))
{
}

imaginary_quadratic_form & imaginary_quadratic_form::operator=(
        imaginary_quadratic_form && o)
{
        ASSERT_ALWAYS(cl == o.cl);
        a = std::move(o.a);
        b = std::move(o.b);
        c = std::move(o.c);
        return *this;
}

std::ostream & operator<<(std::ostream & o, imaginary_quadratic_form const & f)
{
    return o << "(" << f.a << ", " << f.b << ", " << f.c << ")";
}

void
imaginary_quadratic_form::set_c_from_disc()
{
    mpz_mul(c, b, b);
    mpz_sub(c, c, cl.discriminant());
    ASSERT_EXPENSIVE(mpz_divisible_p(c, a));
    mpz_divexact(c, c, a);
    ASSERT_EXPENSIVE(mpz_divisible_ui_p(c, 4U));
    mpz_divexact_ui(c, c, 4U);
}

void
imaginary_quadratic_form::throw_if_a_is_nonpositive() const
{
  if (mpz_sgn(a) <= 0) {
    throw std::domain_error("a coefficient is <= 0");
  }
}

void
imaginary_quadratic_form::throw_if_not_primitive() const
{
  cxx_mpz g;
  mpz_gcd(g, a, b);
  mpz_gcd(g, g, c);
  if (g != 1U) {
    throw not_primitive();
  }
}

void
imaginary_quadratic_form::throw_if_wrong_discriminant() const
{
    cxx_mpz disc;
    mpz_mul(disc, a, c);
    mpz_mul_2exp(disc, disc, 2U);
    mpz_neg(disc, disc);
    mpz_addmul(disc, b, b);
    if (disc != cl.discriminant()) {
        throw std::domain_error("wrong discriminant");
    }
}

void
imaginary_quadratic_form::normalize(cxx_mpz & q, cxx_mpz & r)
{
    mpz_cdiv_qr(q, r, b, a); /* b = q*a + r    and    -a < r <= 0 */
    if (mpz_odd_p(q)) {
      mpz_add(r, r, a);
    }
    mpz_fdiv_q_2exp(q, q, 1U); /* div by 2 */
    /* Now we have b = (2*a)*q + r    and   -a < r <= a */
    mpz_swap(b, r);
    mpz_add(r, b, r);
    mpz_divexact_ui(r, r, 2U);
    mpz_submul(c, q, r);
}

void imaginary_quadratic_form::rho(cxx_mpz & q, cxx_mpz & r)
{
    mpz_swap(a, c);
    mpz_neg(b, b);
    normalize(q, r);
}

void imaginary_quadratic_form::reduction(cxx_mpz & t0, cxx_mpz & t1)
{
    int cmp;
    normalize(t0, t1);
    /* We know a and c > 0, do not consider signs for comparisons */
    while ((cmp = mpz_cmpabs(a, c)) > 0) { /* while a larger than c */
        rho(t0, t1);
    }

    if (cmp == 0 && mpz_sgn(b) < 0) { /* if a == c, we need b positive */
        mpz_neg(b, b);
    }
}

void
imaginary_quadratic_form::reduction()
{
    cxx_mpz t0, t1;
    reduction(t0, t1);
}

/* Naive implementation of composition r <- f1*f2:
 *  compute r.a and r.b using composition formulae, compute r.c using the
 *  discriminant, then reduce the form.
 * Assumes r, f1 and f2 belong to the same class group.
 */
void
imaginary_quadratic_form::compose_mul(imaginary_quadratic_form & r,
                                      imaginary_quadratic_form const & f1,
                                      imaginary_quadratic_form const & f2,
                                      ScratchVars & tmp)
{
    mpz_add(tmp.s, f1.b, f2.b);
    mpz_divexact_ui(tmp.s, tmp.s, 2);
    mpz_gcdext(tmp.g, tmp.u0, tmp.v0, f1.a, f2.a);
    mpz_gcdext(tmp.d, tmp.v1, tmp.w, tmp.g, tmp.s);
    mpz_mul(tmp.v, tmp.v0, tmp.v1);

    mpz_divexact(tmp.a2d, f2.a, tmp.d);

    mpz_mul(r.a, f1.a, tmp.a2d);
    mpz_divexact(r.a, r.a, tmp.d); /* r.a = (f1.a*f2.a)/d^2 */

    mpz_sub(tmp.t, tmp.s, f2.b);
    mpz_mul(tmp.t, tmp.t, tmp.v);
    mpz_submul(tmp.t, tmp.w, f2.c);
    mpz_mul(tmp.t, tmp.t, tmp.a2d);
    mpz_mul_2exp(tmp.t, tmp.t, 1U);
    mpz_add(r.b, f2.b, tmp.t);

    r.set_c_from_disc();
    r.reduction(tmp.u0, tmp.v0);
}

/* Naive implementation of r <- f^2. Same as above but simplified for squaring.
 * Assumes r, f1 and f2 belong to the same class group.
 */
void
imaginary_quadratic_form::compose_sqr(imaginary_quadratic_form & r,
                                      imaginary_quadratic_form const & f,
                                      ScratchVars & tmp)
{
    mpz_gcdext(tmp.d, tmp.u0, tmp.v, f.a, f.b);

    mpz_divexact(tmp.a2d, f.a, tmp.d);

    mpz_mul(r.a, tmp.a2d, tmp.a2d); /* r.a = (f.a/d)^2 */

    mpz_mul(tmp.t, tmp.a2d, f.c);
    mpz_mul(tmp.t, tmp.t, tmp.v);
    mpz_mul_2exp(tmp.t, tmp.t, 1U);
    mpz_sub(r.b, f.b, tmp.t); /* r.b = f.b - 2*v*c*a/d */

    r.set_c_from_disc();
    r.reduction(tmp.u0, tmp.v0);
}

imaginary_quadratic_form
imaginary_quadratic_form::operator*(imaginary_quadratic_form const & f) const
{
    ASSERT_ALWAYS(cl == f.cl);
    ScratchVars scratch;
    imaginary_quadratic_form r = cl.one();
    compose_mul(r, *this, f, scratch);
    return r;
}

imaginary_quadratic_form
imaginary_quadratic_form::operator^(cxx_mpz const & e) const
{
    imaginary_quadratic_form g = *this;
    ScratchVars tmp;
    for (unsigned int i = mpz_sizeinbase(e, 2) - 1; i--; ) {
        compose_sqr(g, g, tmp);
        if (mpz_tstbit(e, i)) {
            compose_mul(g, g, *this, tmp);
        }
    }

    return g;
}
