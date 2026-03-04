#include "cado.h" // IWYU pragma: keep

#include <cmath>
#include <cstdint>

#include <ranges>
#include <ostream>

#include <gmp.h>

#include "gmp_aux.h"
#include "las-qlattice.hpp"
#include "rootfinder.h"
#include "special-q.hpp"
#include "macros.h"
#include "verbose.h"

#include "fmt/ranges.h"

std::ostream& operator<<(std::ostream& os, special_q_data_base const & Q)
{
    return Q.print(os);
}

/* check that the double x fits into an int32_t */
#define fits_int32_t(x) \
  ((double) INT32_MIN <= (x)) && ((x) <= (double) INT32_MAX)

/* We work with two vectors v0=(a0,b0) and v1=(a1,b1). The quadratic form
 * is proportional to a0^2+skewness^2*b0^2 */
static int 
generic_skew_gauss(mpz_t a[2], mpz_t b[2], double skewness)
{
    // skewness^2 can be written as a one-word (actually, 53 bits)
    // integer times a power of two, which is presumably larger than
    // 2^-53 in the most common case, since we expect skewness > 1.
    double mantissa;
    int64_t mantissa_z;
    int exponent;
    /* These are provided by c99 */
    mantissa = frexp(skewness * skewness, &exponent);
    mantissa_z = ldexp(mantissa, 53);
    exponent -= 53;

    mpz_t N[2], S, q, tmp;
    mpz_init(N[0]);
    mpz_init(N[1]);
    mpz_init(S);
    mpz_init(q);
    mpz_init(tmp);

    /* Compute the two norms, and the dot products */
#define QUADFORM(dst, i0, i1) do {					\
    mpz_mul(dst, a[i0], a[i1]);						\
    mpz_mul(tmp, b[i0], b[i1]);						\
    if (exponent < 0) mpz_mul_2exp(dst, dst, -exponent);		\
    if (exponent > 0) mpz_mul_2exp(tmp, tmp, exponent);			\
    mpz_addmul_int64(dst, tmp, mantissa_z);				\
} while (0)
    
    QUADFORM(N[0], 0, 0);
    QUADFORM(N[1], 1, 1);
    QUADFORM(S,    0, 1);

    /* After a reduction step (e.g. v0-=q*v1), N[0], N[1], and S are
     * updated using the following algorithm.
     * new_N0 = old_N0 + q^2old_N1 - 2q old_S
     * new_S = old_S - q old_N1
     *
     * i.e.
     *
     * new_N0 = old_N0 - q old_S - q*(old_S-q old_N1)
     * new_N0 = old_N0 - q * (old_S + new_S)
     *
     * cost: two products only
     */
    for(;;) {
        /* reduce v0 with respect to v1 */
        mpz_rdiv_q(q, S, N[1]);
        if (mpz_cmp_ui(q, 0) == 0) break;
        mpz_submul(a[0], q, a[1]);
        mpz_submul(b[0], q, b[1]);
        /* update */
        mpz_set(tmp, S);
        mpz_submul(S, q, N[1]);
        mpz_add(tmp, tmp, S);
        mpz_submul(N[0], q, tmp);

        /* reduce v1 with respect to v0 */
        mpz_rdiv_q(q, S, N[0]);
        if (mpz_cmp_ui(q, 0) == 0) break;
        mpz_submul(a[1], q, a[0]);
        mpz_submul(b[1], q, b[0]);
        /* update */
        mpz_set(tmp, S);
        mpz_submul(S, q, N[0]);
        mpz_add(tmp, tmp, S);
        mpz_submul(N[1], q, tmp);
    }

    /* We don't care about the sign of b. Down the road, there's an
     * IJToAB function which guarantees positive b */

    /* However we do care about vector 0 being the ``smallest'' in some
     * sense. The trick is that the comparison criterion used previously
     * by the code was the skewed L^\infty norm max(a,b*skewness). We
     * want the reduction step here to be oblivious to this L^\infty /
     * L^2 mess, so we provide something which is L^2 minimal. Wrapper
     * code around this function is still guaranteeing L^\infinity
     * minimality for compatibility, though (see sieve_info_adjust_IJ)
     */
    if (mpz_cmp(N[0], N[1]) > 0) {
        mpz_swap(a[0], a[1]);
        mpz_swap(b[0], b[1]);
    }
    mpz_clear(N[0]);
    mpz_clear(N[1]);
    mpz_clear(S);
    mpz_clear(q);
    mpz_clear(tmp);
    return 1;
}

static int
SkewGauss (qlattice_basis &basis,  mpz_srcptr p, mpz_srcptr r,
           const double skewness)
{
    mpz_t a[2], b[2];
    int fits;

    mpz_init_set (a[0], p);
    mpz_init_set (a[1], r);
    mpz_init_set_ui (b[0], 0);
    mpz_init_set_ui (b[1], 1);
    generic_skew_gauss (a, b, skewness);
    fits = mpz_fits_sint64_p (a[0]);
    fits = fits && mpz_fits_sint64_p (b[0]);
    fits = fits && mpz_fits_sint64_p (a[1]);
    fits = fits && mpz_fits_sint64_p (b[1]);
    if (fits)
      {
        basis.a0 = mpz_get_int64 (a[0]);
        basis.a1 = mpz_get_int64 (a[1]);
        basis.b0 = mpz_get_int64 (b[0]);
        basis.b1 = mpz_get_int64 (b[1]);
      }
    mpz_clear (a[0]);
    mpz_clear (a[1]);
    mpz_clear (b[0]);
    mpz_clear (b[1]);
    return fits ? 1 : 0;
}

qlattice_basis::qlattice_basis(special_q const & doing, double skew) :
    special_q_data_base(doing)
{
    /* Currently requires prime or composite square-free special-q
     * values, For powers, the base prime would have to be determined and
     * stored in a variable, so that powers of that prime in the factor
     * base can be skipped over.  */
    ASSERT_ALWAYS(!mpz_perfect_power_p(doing.p));

    if (!SkewGauss (*this, doing.p, doing.r, skew)) {
        throw too_skewed();
    }
}

/*
 *      (a)   (a0 a1)   (i)
 *      (b) = (b0 b1) * (j)
 */
void
qlattice_basis::convert_ij_to_ab(
        int64_t & a,
        uint64_t & b,
        int i,
        unsigned int j) const
{
    sublat.adjustIJ(i, j);

    int64_t const s = (int64_t)i * (int64_t) a0 + (int64_t)j * (int64_t) a1;
    int64_t const t = (int64_t)i * (int64_t) b0 + (int64_t)j * (int64_t) b1;
    if (t >= 0) {
        a = s;
        b = t;
    } else {
        a = -s;
        b = -t;
    }
}

void
qlattice_basis::convert_ij_to_ab(
        cxx_mpz & a,
        cxx_mpz & b,
        int i,
        unsigned int j) const
{
    sublat.adjustIJ(i, j);

    cxx_mpz zi = i;
    cxx_mpz zj = j;

    mpz_mul_int64(a, zi, a0);
    mpz_addmul_int64(a, zj, a1);

    mpz_mul_int64(b, zi, b0);
    mpz_addmul_int64(b, zj, b1);

    if (mpz_sgn(b) < 0) {
        mpz_neg(a, a);
        mpz_neg(b, b);
    }
}

/*
 *      (i)   (a0 a1)^-1   (a)   (b1 -a1)         (a)
 *      (j) = (b0 b1)    * (b) = (-b0 a0) / det * (b)
 * with
 *      det = doing.p
 */
int
qlattice_basis::convert_ab_to_ij(
        int & i,
        unsigned int & j,
        int64_t a,
        uint64_t b) const
{
    /* Both a,b and the coordinates of the lattice basis can be quite
     * large. However the result should be small.
     */
    cxx_mpz za, zb;
    mpz_set_int64(za, a);
    mpz_set_uint64(zb, b);
    return convert_ab_to_ij(i, j, za, zb);
}

int
qlattice_basis::convert_ab_to_ij(
        int & i,
        unsigned int & j,
        cxx_mpz const & a,
        cxx_mpz const & b) const
{
    /* Both a,b and the coordinates of the lattice basis can be quite
     * large. However the result should be small.
     */
    cxx_mpz ii, jj;
    int ok = 1;
    mpz_mul_int64(ii, a, b1); mpz_submul_int64(ii, b, a1);
    mpz_mul_int64(jj, b, a0); mpz_submul_int64(jj, a, b0);
    /*
    int64_t ii =   a * (int64_t) b1 - b * (int64_t)a1;
    int64_t jj = - a * (int64_t) b0 + b * (int64_t)a0;
    */
    if (!mpz_divisible_p(ii, doing.p)) ok = 0;
    if (!mpz_divisible_p(jj, doing.p)) ok = 0;
    mpz_divexact(ii, ii, doing.p);
    mpz_divexact(jj, jj, doing.p);

    if (mpz_sgn(jj) < 0 || (mpz_sgn(jj) == 0 && mpz_sgn(ii) < 0)) {
        mpz_neg(ii, ii);
        mpz_neg(jj, jj);
    }
    i = mpz_get_si(ii);
    j = mpz_get_ui(jj);
    if (sublat.m != 0) {
        int64_t imodm = i % int64_t(sublat.m);
        if (imodm < 0) {
            imodm += sublat.m;
        }
        int64_t jmodm = j % int64_t(sublat.m);
        if (jmodm < 0)
            jmodm += sublat.m;
        if (imodm != sublat.i0 || jmodm != sublat.j0) {
            verbose_fmt_print(1, 0, "# TraceAB: (i,j)=({},{}) does not belong "
                                    "to the right congruence class\n", i, j);
            ok = 0;
        } else {
            i = (i - imodm) / int64_t(sublat.m);
            j = (j - jmodm) / int64_t(sublat.m);
        }
    }

    return ok;
}

std::ostream & qlattice_basis::print(std::ostream & os) const
{
    return os << *this;
}

std::ostream& operator<<(std::ostream& os, qlattice_basis const & Q)
{
    if (mpz_cmp_ui(Q.doing.p, 0) != 0)
        os << Q.doing << ";";

    os  << " a0=" << Q.a0 << ";"
        << " b0=" << Q.b0 << ";"
        << " a1=" << Q.a1 << ";"
        << " b1=" << Q.b1;
    return os;
}

siqs_special_q_data::siqs_special_q_data(
        special_q const & sq,
        cxx_cado_poly const & cpoly)
    : special_q_data_base(sq)
{
    //TODO assert cpoly is x^2+cst
    cxx_gmp_randstate rstate;
    crt_data_modq.clear();
    crt_data_modq.reserve(doing.prime_factors.size());
    doing.r = 0u;
    for (auto const & qk: doing.prime_factors) {
        /* crt = q/qk * (1/(q/qk) mod qk), i.e.,
         *    crt = 0 mod q/qk and crt = 1 mod qk
         */
        cxx_mpz crt = doing.p / qk;
        crt *= crt.invmod(qk);

        std::vector<uint64_t> rk = mpz_poly_roots(cpoly[sq.side], qk, rstate);
        ASSERT_ALWAYS(rk.size() == 2);

        /* XXX assume poly is even, only one root is considered, -Rk will
         * correspond to the other root.
         */
        cxx_mpz t0 = rk[0] * crt;
        mpz_mod(t0, t0, doing.p);
        if (2*t0 >= doing.p)
            t0 = doing.p - t0;
        crt_data_modq.emplace_back(t0);

        doing.r += t0;
    }
}

void
siqs_special_q_data::convert_ij_to_ab(
        int64_t & a,
        uint64_t & b,
        int i,
        unsigned int j) const
{
    /* b = 1; a = rj + q*i */
    b = 1u;
    a = mpz_get_si(root_from_j(j));
    a += mpz_get_si(doing.p) * i;
}

void
siqs_special_q_data::convert_ij_to_ab(
        cxx_mpz & a,
        cxx_mpz & b,
        int i,
        unsigned int j) const
{
    /* b = 1; a = rj + q*i */
    b = 1u;
    a = root_from_j(j);
    mpz_addmul_si(a, doing.p, i);
}

int
siqs_special_q_data::convert_ab_to_ij(
        int & i,
        unsigned int & j,
        int64_t a,
        uint64_t b) const
{
    cxx_mpz za, zb;
    mpz_set_int64(za, a);
    mpz_set_uint64(zb, b);
    return convert_ab_to_ij(i, j, za, zb);
}

int
siqs_special_q_data::convert_ab_to_ij(
        int & i,
        unsigned int & j,
        cxx_mpz const & a,
        cxx_mpz const & b) const
{
    if (mpz_cmp_ui(b, 1u) != 0) {
        return 0; /* only b == 1 is sieved */
    }

    /* Here we will computed j directly without computing the its gray code g,
     * using the following facts:
     *   the k-th bit of j = xor(l-th bit of g for l >= k)
     *   the l-th bit of g = 0 if a % qk == Rk % qk else 1
     *
     * We will compute the correspond rj at the same time as it is needed to
     * compute i = (a-rj)/q.
     */
    j = 0u;
    cxx_mpz rj = 0u;
    uint64_t mask = 1u; /* invariant: mask = 2^(k+1)-1 = 0b11..1 with k 1's */
    for(unsigned int k = 0; k < nfactors(); ++k, mask = (mask << 1u) + 1u) {
        const uint64_t qk = doing.prime_factors[k];
        cxx_mpz const & Rk = crt_data_modq[k];

        uint64_t amodqk = mpz_tdiv_uint64(a, qk);
        if (mpz_sgn(a) < 0) {
            amodqk = qk - amodqk;
        }
        if (mpz_tdiv_uint64(Rk, qk) == amodqk) {
            rj += Rk;
        } else {
            ASSERT((qk - mpz_tdiv_uint64(Rk, qk)) == amodqk);
            rj -= Rk;
            j = j xor mask;
        }
    }

    cxx_mpz t;
    mpz_sub(t, rj, a);
    ASSERT(mpz_divisible_p(t, doing.p));
    mpz_divexact(t, t, doing.p);
    mpz_neg(t, t); /* t = (a-rj)/q  which is i */
    if (t.fits<int>()) {
      i = mpz_get_si(t);
      return 1;
    } else {
      return 0;
    }
}

std::ostream & siqs_special_q_data::print(std::ostream & os) const
{
    return os << *this;
}

std::ostream& operator<<(std::ostream& os, siqs_special_q_data const & Q)
{
    if (mpz_cmp_ui(Q.doing.p, 0) != 0)
        os << Q.doing << "; ";

    auto p = [&, Q](size_t i) { return fmt::format("({}, {})",
                            Q.doing.prime_factors[i], Q.crt_data_modq[i]); };
    auto v = std::views::iota(0u, Q.nfactors()) | std::views::transform(p);
    os << "CRT data = [" << fmt::format("{}", fmt::join(v, ", ")) << "; "
       << "r0 = " << Q.doing.r;

    return os;
}
