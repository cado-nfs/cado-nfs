#include "cado.h" // IWYU pragma: keep

#include "gmp_aux.h"
#include "all_valuations_above_p.hpp"
#include "numbertheory_internals.hpp"
#include "number_field_element.hpp"
#include "mpz_poly.h"

numbertheory_internals::all_valuations_above_p::all_valuations_above_p(cxx_mpz_poly const& f, cxx_mpz const& p, cxx_gmp_randstate & state)
    : f(f)
    , p(p)
    , K(number_field(f))
    , O(K.p_maximal_order(p))
    , F(O.factor_radical(p, state))
    , jjinv(O.fractional_ideal({{ K(1), K.gen() }}))
{
    K.bless("K", "alpha");
    O.bless(fmt::format("OK{}", p));

    for(auto const & Ie : F) inertia.push_back(Ie.inertia_degree());
    val_base.assign(f->deg, 0);
    if (f->deg > 1) {
        val_base = (*this)(jjinv);
    } else {
        int v = 0;
        cxx_mpz lc = mpz_poly_lc(f);
        for( ; mpz_divisible_p(lc, p) ; mpz_divexact(lc, lc, p), v++);
        /* to be honest, I'm not sure if we should put v or -v, here.
         * Well, it doesn't really seem to make sense anyway to deal with
         * number fields of degree one, does it? Chances are we'll never
         * reach this code branch with p dividing the leading coefficient
         * anyway. If we ever do, we'll have a test case that allows us
         * to find out more.
         */
        ASSERT_ALWAYS(v == 0);
        val_base.push_back(-v);
    }
}

void numbertheory_internals::all_valuations_above_p::bless_side(int side)
{
    K.bless(fmt::format("K{}", side), fmt::format("alpha{}", side));
    O.bless(fmt::format("OK{}", side));
}

std::vector<int> numbertheory_internals::all_valuations_above_p::operator()(number_field_fractional_ideal const & I) const
{
    std::vector<int> res;
    res.reserve(F.size());
    for(auto const & fkp : F)
        res.push_back(I.valuation(fkp));
    return res;
}

std::string numbertheory_internals::all_valuations_above_p::sagemath_string(int k, int)
{
    return fmt::format("{:S}", F[k]);
}

std::string numbertheory_internals::all_valuations_above_p::machine_description(int k) {
    /* This is about the same as above, in that we also return a
     * two-element form of the ideal, but here we return it in a
     * machine readable way. First the prime above which our ideal
     * sits, then a denominator, then the coefficients of the
     * polynomial in alpha that define the second element.
     */
    return fmt::format("{:M}", F[k]);
    /*
    number_field_prime_ideal::two_element two = F[k];
    number_field_element theta = K(two.second);

    std::vector<cxx_mpz> res = numbertheory_internals::write_element_as_list_of_integers(theta.coefficients);
    res.insert(res.begin(), two.first);
    return res;
    */
}

void numbertheory_internals::all_valuations_above_p::print_info(std::ostream& o, int k, cxx_mpz const& r MAYBE_UNUSED, int side) const
{
    o << fmt::format("# I{}={:S}; # e={} f={}\n",
            k,
            F[k],
            F[k].ramification_index(),
            F[k].inertia_degree());
    o << fmt::format("# I_{}_{}_{}_{} # {} {}\n",
            p, r, side, k,
            F[k].ramification_index(),
            F[k].inertia_degree());
}

/* create ideal I=<p^k,p^k*alpha,v*alpha-u> and decompose I*J */
std::vector<int> numbertheory_internals::all_valuations_above_p::operator()(int k, cxx_mpz const& r) const
{
    cxx_mpz pk;
    mpz_pow_ui(pk, p, k);
    cxx_mpz_poly A;
    mpz_poly_set_xi(A, 1);
    if (mpz_cmp(r, pk) < 0) {
        mpz_neg(mpz_poly_coeff(A, 0), r);
    } else {
        mpz_set_si(mpz_poly_coeff(A,0),-1);
        mpz_sub(mpz_poly_coeff(A,1), r, pk);
        mpz_poly_cleandeg(A, 1);
    }
    /* hell, do I _really_ need p*alpha here ??? */
    auto v = (*this)(O.fractional_ideal({K(pk), K(A), K(pk) * K.gen()}));
    ASSERT_ALWAYS(v.size() == val_base.size());
    for(size_t i = 0 ; i < v.size() ; i++)
        v[i] -= val_base[i];
    return v;
}

std::vector<int> numbertheory_internals::all_valuations_above_p::multiply_inertia(std::vector<int> const& v) const {
    ASSERT_ALWAYS(v.size() == inertia.size());
    std::vector<int> res(v.size(),0);
    for(unsigned int i = 0 ; i < v.size() ; i++) {
        res[i] = v[i] * inertia[i];
    }
    return res;
}


