#ifndef CADO_ARITH_MODP_MAIN_HPP
#define CADO_ARITH_MODP_MAIN_HPP

#include <cstddef>
#include <cstdio>
#include <cstdint>

#include <algorithm>
#include <limits>
#include <type_traits>
#include <utility>
#include <array>
#include <string>
#include <vector>
#include <memory>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "fmt/format.h"
#include "gmp-hacks.h"
#include "gmp_aux.h"
#include "macros.h"
#include "memory.h"
#include "runtime_numeric_cast.hpp"

#include "arith-concrete-base.hpp"

#include "mpn_compile_time.hpp"

/* the elt objects are always reduced between 0 and p-1, so that it is
 * fine to have p as large as N*sizeof(mp_limb_t).
 *
 * However the unreduced elements are implicitly signed, and the leading
 * bit is here to reflect that. The reduction function takes this into
 * account. All subtraction, negation functions on unreduced elements are
 * therefore avoiding all normalizations.
 */

#define xxxDEBUG_INFINITE_LOOPS

namespace arith_modp::details {

    /* There is a good deal of complication that comes from the fact that
     * we want the same code to work for fixed-size types as well as for
     * variable-size types. In the latter case, the sizeof() of the
     * element type cannot (of course) be the true size, so we have to
     * resort to somewhat ugly gymnastics
     */

    struct mp_limb_array_base {
        /* This type is used only as a parent of child structures, and
         * only has the ability to convert to an mp_limb_t * pointer
         */
        mp_limb_t* pointer() { return reinterpret_cast<mp_limb_t *>(this); }
        const mp_limb_t* pointer() const { return reinterpret_cast<mp_limb_t const *>(this); }

        /* This makes the .pointer() syntax mostly superfluous */
        operator mp_limb_t*() { return pointer(); }
        operator const mp_limb_t*() const { return pointer(); }
    };

    template<size_t N>
    class mp_limb_array : public mp_limb_array_base {
        /* We only add this so that the concrete types have the correct
         * width in the fixed-width case. In the variable case, the
         * object size is meaningless anyway.
         * Note that we do _not_ want to inherit, because of the
         * array::pointer and array::operator[] which get in the way.
         */
        std::array<mp_limb_t, N> dummy;
    };

/* compile-time table that returns the largest power of two (but at most
 * 16) that divides k */
template<unsigned int k,unsigned int ell = 1, int divide_out = (ell < 16 && !(k & 1))>
struct alignment_divisor;
template<unsigned int k, unsigned int ell> struct alignment_divisor<k,ell,0> : public std::integral_constant<unsigned int, ell> {};
template<unsigned int k, unsigned int ell> struct alignment_divisor<k,ell,1> : public alignment_divisor<k/2,ell*2> {};

#ifdef HAVE_ALLOCA
#define ARITH_MODP_TEMPORARY_ALLOC(parent, type__, var) \
    mp_limb_t * CPP_PAD(var, _allocation) = alloca((parent)->stride<typename std::remove_pointer_t<decltype(parent)>::type__>()); \
    auto & var(*reinterpret_cast<typename std::remove_pointer_t<decltype(parent)>::type__ *>(CPP_PAD(var, _allocation)))
#else
#define ARITH_MODP_TEMPORARY_ALLOC(parent, type__, var) \
    auto CPP_PAD(var, _allocation) = typename std::unique_ptr<mp_limb_t[]>(new mp_limb_t[(parent)->template nlimbs<typename std::remove_pointer_t<decltype(parent)>::type__>()]);                                                               \
    auto & var(*reinterpret_cast<typename std::remove_pointer_t<decltype(parent)>::type__ *>(CPP_PAD(var, _allocation).get()))
#endif


template<size_t NN, typename T>
struct gfp_base : public arith_concrete_base
{
    static constexpr bool is_binary = false;
    static constexpr const size_t constant_width = NN;
    static constexpr const bool is_constant_width = NN > 0;
    static constexpr const bool is_variable_width = NN == 0;
    static std::string impl_name() { return is_constant_width ? fmt::format("p{}", NN) : "pz"; }
    static constexpr const bool is_characteristic_two = false;
    static constexpr int simd_groupsize() { return 1; }

    /* These three are empty structs. We only ever use them as
     * references.
     *
     * XXX These types are _not_ constructible!
     * XXX These types are _not_ assignable!
     * XXX These types are _not_ copyable!
     *
     * For pointer arithmetic, use ab->vec_subvec
     * For element access, use ab->vec_item
     * For assignment, use ab->set
     * For vector move/copy, use ab->vec_set
     * For element construction, use ab->alloc (or ARITH_MODP_TEMPORARY_ALLOC, which sets a reference to something allocated on the stack)
     *
     */
    struct elt : public arith_concrete_base::elt, public mp_limb_array<NN> {
        static constexpr int classifier = 0;
        static constexpr int alignment = (is_constant_width ? alignment_divisor<NN>::value : 1) * sizeof(mp_limb_t);
        elt() = delete;
        elt& operator=(elt const &) = delete;
        elt(elt const &) = delete;
        elt& operator=(elt &&) = delete;
        elt(elt &&) = delete;
        ~elt() = default;
    };
    struct elt_ur_for_add : public mp_limb_array<NN + (NN > 0)> {
        static constexpr int classifier = 1;
        static constexpr int alignment = (is_constant_width ? alignment_divisor<(NN + 1)>::value : 1) * sizeof(mp_limb_t);
        elt_ur_for_add() = delete;
        elt_ur_for_add& operator=(elt_ur_for_add const &) = delete;
        elt_ur_for_add(elt_ur_for_add const &) = delete;
        elt_ur_for_add& operator=(elt_ur_for_add &&) = delete;
        elt_ur_for_add(elt_ur_for_add &&) = delete;
        ~elt_ur_for_add() = default;
    };
    struct elt_ur_for_addmul : public mp_limb_array<2 * NN + (NN > 0)> {
        static constexpr int classifier = 2;
        static constexpr int alignment = sizeof(mp_limb_t);
        elt_ur_for_addmul() = delete;
        elt_ur_for_addmul& operator=(elt_ur_for_addmul const &) = delete;
        elt_ur_for_addmul(elt_ur_for_addmul const &) = delete;
        elt_ur_for_addmul& operator=(elt_ur_for_addmul &&) = delete;
        elt_ur_for_addmul(elt_ur_for_addmul &&) = delete;
        ~elt_ur_for_addmul() = default;
    };

    /*{{{ element layout information */

    /* This is just very verbosely saying that elt has N limbs,
     * elt_ur_for_add has N+1 limbs, and elt_ur_for_addmul has 2N+1
     * limbs. We want to have as many std::integral_constant types as we
     * can.
     *
     */
    template<typename X>
        std::integral_constant<size_t, 0> overhead_limbs() const
        requires std::is_same_v<X, elt>
        { return {}; }
    template<typename X>
        static std::integral_constant<size_t, 1> overhead_limbs()
        requires std::is_same_v<X, elt_ur_for_add>
        { return {}; }
    template<typename X>
        std::integral_constant<size_t, NN+1> overhead_limbs() const
        requires is_constant_width && std::is_same_v<X, elt_ur_for_addmul>
        { return std::integral_constant<size_t, NN+1>(); }
    template<typename X>
        size_t overhead_limbs() const
        requires is_variable_width && std::is_same_v<X, elt_ur_for_addmul>
        { return nlimbs<elt>() + 1; }

    template<typename X>
        std::integral_constant<size_t, NN> nlimbs() const
        requires is_constant_width && std::is_same_v<X, elt>
        { return std::integral_constant<size_t, NN>(); }
    template<typename X>
        size_t nlimbs() const 
        requires is_variable_width && std::is_same_v<X, elt>
        { return mpz_size(p); }
    template<typename X>
        std::integral_constant<size_t, NN+1> nlimbs() const
        requires is_constant_width && std::is_same_v<X, elt_ur_for_add>
        { return std::integral_constant<size_t, NN+1>(); }
    template<typename X>
        std::integral_constant<size_t, 2 * NN+1> nlimbs() const
        requires is_constant_width && std::is_same_v<X, elt_ur_for_addmul>
        { return std::integral_constant<size_t, 2 * NN+1>(); }
    template<typename X>
        size_t nlimbs() const
        requires is_variable_width && std::is_same_v<X, elt_ur_for_add>
        { return nlimbs<elt>() + 1; }
    template<typename X>
        size_t nlimbs() const
        requires is_variable_width && std::is_same_v<X, elt_ur_for_addmul>
        { return 2*nlimbs<elt>() + 1; }

    template<typename X>
    size_t stride() const { return sizeof(mp_limb_t) * nlimbs<X>(); }
    size_t elt_stride() const { return stride<elt>(); }
    size_t vec_elt_stride(size_t s) const { return s * elt_stride(); }

    /*}}}*/

    private:
    cxx_mpz p;
    /* How many bits to he have to shift p left so that its
     * highest bit is set?
     */
    int preinv_shift = 0;

    /* prime_preinv_for_add has just one limb. prime_preinv_for_addmul
     * has N+1 limbs (which is why we store it as vector: in the variable
     * size case, N is a runtime value).
     *
     * There are corner cases where prime_preinv_for_addmul might have
     * cancellations at high words (p=2^192-237 is one example). We
     * choose to keep a full-length limb list in that case, though.
     */
    mp_limb_t prime_preinv_for_add[decltype(overhead_limbs<elt_ur_for_add>())::value];
    std::vector<mp_limb_t> prime_preinv_for_addmul;
    public:
    gfp_base()
      : p(0)
    {}
    gfp_base(gfp_base const&) = default;
    gfp_base(gfp_base&&) = default;
    gfp_base& operator=(gfp_base const&) = default;
    gfp_base& operator=(gfp_base&&) = default;
    ~gfp_base() = default;
    gfp_base(cxx_mpz const& p, unsigned int simd_groupsize)
      : p(p)
    {
        auto N = nlimbs<elt>();
        ASSERT_ALWAYS(simd_groupsize == 1);
        const size_t m = mpz_sizeinbase(p, 2);
        ASSERT_ALWAYS(NN == 0 || m <= (size_t) N * mp_bits_per_limb);
        preinv_shift = (mp_bits_per_limb - m) % mp_bits_per_limb;
        cxx_mpz blah;
        compute_preinv(blah, overhead_limbs<elt_ur_for_add>());
        mpn_compile_time::MPN_SET_MPZ(prime_preinv_for_add, overhead_limbs<elt_ur_for_add>(), blah);
        prime_preinv_for_add[0] = mpz_get_ui(blah);
        compute_preinv(blah, overhead_limbs<elt_ur_for_addmul>());
        prime_preinv_for_addmul.clear();
        prime_preinv_for_addmul.reserve(overhead_limbs<elt_ur_for_addmul>());
        for(size_t i = 0 ; i < overhead_limbs<elt_ur_for_addmul>() ; i++)
            prime_preinv_for_addmul.push_back(mpz_getlimbn(blah, int(i)));
    }
    mp_limb_t const* prime_limbs() const { return mpz_limbs_read(p); }
    // elt const& prime() const { return *reinterpret_cast<mpn<N> const
    // *>(prime_limbs()); }


    /* The helpers from the mpn_compile_time namespace will resolve to
     * compile-time loops whenever possible, so we may assume that the
     * default set(), is_zero(), or cmp() functions, in particular, are fine.
     */

    mpz_srcptr characteristic() const { return p; }

    /*{{{ allocation / deallocation of (vectors of) elements */
    /* These allocation interfaces seem a bit stupid. At the low
     * hard level, we know that they're just the same as new[]
     * anyway.
     */

    template<typename X = elt>
        X* alloc(size_t k = 1, size_t al = X::alignment) const
        requires std::is_base_of_v<mp_limb_array_base, X>
        {
            return reinterpret_cast<X *>(::malloc_aligned(k * stride<X>(), al));
        }
    template<typename X = elt>
        void free(X* u) const
        requires std::is_base_of_v<mp_limb_array_base, X>
        {
            ::free_aligned(reinterpret_cast<void *>(u));
        }
    template<typename X = elt>
        X* realloc(X* u, size_t k0, size_t k, size_t al = X::alignment) const
        requires std::is_base_of_v<mp_limb_array_base, X>
        {
            return reinterpret_cast<X *>(::realloc_aligned(u->pointer(), k0 * stride<X>(), k * stride<X>(), al));
            /*
               X* v = alloc<X>(k);
               std::copy_n(u, k0, v);
               free<X>(u);
               return v;
               */
        }

    /*}}}*/
    /*{{{ predicates */
    template<typename X>
        bool is_zero(X const& x) const
        requires (X::classifier >= 0)
        {
            return mpn_compile_time::mpn_zero_p(x, nlimbs<X>());
        }

    template<typename X>
        int cmp(X const& x, unsigned long a) const
        requires (X::classifier >= 0)
        {
            return mpn_compile_time::mpn_cmp_ui(x, a, nlimbs<X>());
        }

    template<typename X>
        int cmp(X const& x, X const& y) const
        requires (X::classifier >= 0)
        {
            return mpn_compile_time::mpn_cmp(x, y, nlimbs<X>());
        }

    /*}}}*/
    /* {{{ assignments */
    template<typename X>
        X& set(X& x, unsigned long a) const
        requires (X::classifier >= 0)
        {
            mpn_compile_time::mpn_set_ui(x, a, nlimbs<X>());
            return x;
        }

    template<typename X>
        X& set(X& x, X const& a) const
        requires (X::classifier >= 0)
        {
            mpn_compile_time::mpn_copyi(x, a, nlimbs<X>());
            return x;
        }

    template<typename X>
        X& set(X& a, elt const& x) const
        requires (X::classifier >= 1)
        {
            auto N = nlimbs<elt>();
            mp_limb_t * ax = a.pointer() + N;
            auto nx = overhead_limbs<X>();
            mpn_compile_time::mpn_copyi(a, x, N);
            mpn_compile_time::mpn_zero(ax, nx);
            return a;
        }

    template<typename X>
        X& set(X& x, cxx_mpz const& a) const
        requires (X::classifier >= 0)
        {
            mpn_compile_time::MPN_SET_MPZ(x.pointer(), nlimbs<X>(), a);
            return x;
        }

    template<typename X>
        void set_zero(X& x) const
        requires (X::classifier >= 0)
        {
            mpn_compile_time::mpn_zero(x, nlimbs<X>());
        }

    template<typename X>
        void set_random(X& x, cxx_gmp_randstate & rstate) const
        requires (X::classifier >= 0)
        {
            T const* tx = static_cast<T const*>(this);
            cxx_mpz xz;
            mpz_urandomm(xz, rstate, p);
            tx->set(x, xz);
        }

    bool upperlimbs_are_zero(elt_ur_for_add const & a) const
    {
        return a[nlimbs<elt>()] == 0;
    }

    bool upperlimbs_are_zero(elt_ur_for_addmul const & a) const
    {
        auto N = nlimbs<elt>();
        mp_limb_t const * ax = a.pointer() + N;
        auto nx = overhead_limbs<elt_ur_for_addmul>();
        return mpn_compile_time::mpn_zero_p(ax, nx);
    }

    void stream_store(elt* dst, elt const& src) const {
        T const* tx = static_cast<T const*>(this);
        tx->set(*dst, src);
    }
    /* }}} */

    /*{{{ addition, at the element level (elt or wider) */
    void propagate_carry(elt_ur_for_addmul& a, mp_limb_t cy) const
    {
        auto N = nlimbs<elt>();
        mp_limb_t * ax = a.pointer() + N;
        auto nx = overhead_limbs<elt_ur_for_addmul>();
        mpn_add_1(ax, ax, nx, cy);
    }
    void propagate_carry(elt_ur_for_add& a, mp_limb_t cy) const
    {
        auto N = nlimbs<elt>();
        a[N] += cy;
    }

    template<typename X>
        void add(X& dst, elt const& a, elt const& b) const
        requires (X::classifier >= 1)
        {
            T const* tx = static_cast<T const*>(this);
            auto N = nlimbs<elt>();
            mp_limb_t cy = mpn_add_n(dst, a, b, N);
            mp_limb_t * dx = dst.pointer() + N;
            auto nx = overhead_limbs<X>();
            mpn_compile_time::mpn_zero(dx, nx);
            tx->propagate_carry(dst, cy);
        }
    template<typename X>
        void add(X& dst, elt const& src) const
        requires (X::classifier >= 1)
        {
            T const* tx = static_cast<T const*>(this);
            auto N = nlimbs<elt>();
            mp_limb_t cy = mpn_add_n(dst, dst, src, N);
            tx->propagate_carry(dst, cy);
        }
    /* This addition is only for unreduced types. These types are
     * always considered wide enough so that overflows work.
     */
    template<typename X>
        void add_ur(X& dst, X const& src) const
        requires (X::classifier >= 1)
        {
            T const* tx = static_cast<T const*>(this);
            tx->add(dst, src);
        }
    template<typename X>
        void add(X& dst, X const& src) const
        requires (X::classifier >= 1)
        {
            mpn_add_n(dst, dst, src, nlimbs<X>());
        }

    /* an addmul is ok for to go for an unreduced type which is
     * still somewhat narrow (only one extra limb).
     */
    template<typename X>
        void addmul_ui(X& dst, elt const& src, mp_limb_t x) const
        requires (X::classifier >= 1)
        {
            T const* tx = static_cast<T const*>(this);
            auto N = nlimbs<elt>();
            mp_limb_t cy = mpn_addmul_1(dst, src, N, x);
            tx->propagate_carry(dst, cy);
        }
    /*}}}*/

    /*{{{ subtraction, at the element level (elt or wider) */
    void propagate_borrow(elt_ur_for_addmul& a, mp_limb_t cy) const
    {
        auto N = nlimbs<elt>();
        mp_limb_t * ax = a.pointer() + N;
        auto nx = overhead_limbs<elt_ur_for_addmul>();
        mpn_sub_1(ax, ax, nx, cy);
    }
    void propagate_borrow(elt_ur_for_add& a, mp_limb_t cy) const
    {
        auto N = nlimbs<elt>();
        a[N] -= cy;
    }

    template<typename X>
        void sub(X& dst, elt const& a, elt const& b) const
        requires (X::classifier >= 1)
        {
            T const* tx = static_cast<T const*>(this);
            auto N = nlimbs<elt>();
            mp_limb_t cy = mpn_sub_n(dst, a, b, N);
            mp_limb_t * dx = dst.pointer() + N;
            auto nx = overhead_limbs<X>();
            mpn_compile_time::mpn_zero(dx, nx);
            tx->propagate_borrow(dst, cy);
        }
    template<typename X>
        void sub(X& dst, elt const& src) const
        requires (X::classifier >= 1)
        {
            T const* tx = static_cast<T const*>(this);
            auto N = nlimbs<elt>();
            mp_limb_t cy = mpn_sub_n(dst, dst, src, N);
            tx->propagate_borrow(dst, cy);
        }
    /* This subtraction is only for unreduced types. These types are
     * always considered wide enough so that overflows work.
     */
    template<typename X>
        void sub_ur(X& dst, X const& src) const
        requires (X::classifier >= 1)
        {
            T const* tx = static_cast<T const*>(this);
            tx->sub(dst, src);
        }
    template<typename X>
        void sub(X& dst, X const& src) const
        requires (X::classifier >= 1)
        {
            mpn_sub_n(dst, dst, src, nlimbs<X>());
        }

    /* a submul is ok for to go for an unreduced type which is
     * still somewhat narrow (only one extra limb).
     */
    template<typename X>
        void submul_ui(X& dst, elt const& src, mp_limb_t x) const
        requires (X::classifier >= 1)
        {
            T const* tx = static_cast<T const*>(this);
            auto N = nlimbs<elt>();
            mp_limb_t cy = mpn_submul_1(dst, src, N, x);
            tx->propagate_borrow(dst, cy);
        }
    /*}}}*/

    /*{{{ neg (for elt and anything larger) */

    /* XXX It's a rocky road here.
     *
     * The elt negation keeps data between 0 and p-1.
     * On the other hand, negation of an elt_ur does not care.
     * Right now, I have the impression that this is very
     * dangerous. It would only make sense if we are able to
     * statically have an idea of the "sign" of the unreduced
     * elements we're dealing with.
     */
    void neg(elt& dst, elt const& src) const
    {
        T const* tx = static_cast<T const*>(this);
        auto N = nlimbs<elt>();
        mpn_neg(dst, src, N);
        if (!tx->is_zero(dst))
            mpn_add_n(dst, dst, tx->prime_limbs(), N);
    }

    template<typename X>
        void neg(X& dst, X const& src) const
        requires (X::classifier >= 1)
        {
            mpn_neg(dst, src, nlimbs<X>());
        }
    /*}}}*/
    /*{{{ signed barrett reduction */

    /*{{{ compute_preinv */

    /* Preinverse for Barrett reduction. See also the code for
     * reduction, which is further below.
     *
     * We want to apply this to reduce a mod p, with the following
     * constraints.
     *
     *             2^(m-1) <  p  < 2^m
     *        -2^(ell-1)*p <= a  < 2^(ell-1)*p
     *
     * Let I=floor(2^(m+ell)/p). Because of the bound on p, we have 2^ell
     * < I < 2^(ell+1), so that 0<J=I-2^ell<2^ell (which actually fits
     * within ell bits). The preinverse we compute is this J.
     */

    /* e is the number of extra limbs */
    /* XXX we should make this (as well as n) an unsigned int */
    void compute_preinv(cxx_mpz & big, int e)
    {
        mpz_set_ui(big, 1);
        size_t m = mpz_sizeinbase(p, 2);
        const size_t ell = e * mp_bits_per_limb;
        mpz_mul_2exp(big, big, m + ell);
        mpz_fdiv_q(big, big, p);
        ASSERT_ALWAYS(mpz_sizeinbase(big, 2) == (ell + 1));
        mpz_fdiv_r_2exp(big, big, ell);
    }
    /*}}}*/
    /*{{{ reduce */
    /* Signed Barrett reduction (extended from Brent-Zimmermann 2010,
     * theorem 2.4)
     */

    /* input: a = a0 + a1 * 2^m, with         2^(m-1) <  p  < 2^m
     *                                   -2^(ell-1)*p <= a  < 2^(ell-1)*p
     *                                              0 <= a0 < 2^m
     * which imply in particular:
     *                                     -2^(ell-1) <= a1 < 2^(ell-1)
     *
     * Case a1 >= 0.
     *
     * Let q0 = floor(a1*I/2^ell) = floor(a1*J/2^ell) + a1.
     * We have 0 <= q0 < 2^ell.
     *
     * Moreover: q0 <= a1*I/2^ell <= a1*2^m/p <= a/p, so that r0=a-p*q0>=0.
     * use p*I >= 2^(m+ell)-p and 2^ell*q0 >= a1*I-2^ell
     *
     * compute 2^ell*p*q0 >= 2^(m+ell)*a1-a1*p-2^ell*p
     *                    >= 2^ell*(a-a0)-p*(a1+2^ell)
     *                    >  2^ell*a - 4*2^ell*p
     *             a-p*q0 <  4p
     * where the third line used a1 < 2^(ell-1) and a0 <= 2^m <= 2*p.
     *
     * Case a1 < 0.
     *
     * We let b1 = a1 + 2^ell, which is the unsigned limb used to
     * represent a1.
     *
     * Let q0 = floor(a1*I/2^ell) = floor(b1*J/2^ell) + b1 - 2^ell - J.
     *
     * Since a1 < 0, we have q0 < 0. With a1 >= -2^(ell-1) and
     * I<2^(ell+1), we obtain q0 > -2^ell. Therefore q0 is well
     * represented by the machine word
     *  q'0 = q0+2^ell = floor(b1*J/2^ell) + b1 - J
     *
     * We have p*q0 <= p*a1*I/2^ell < p*a1/2^ell*(2^(m+ell)/p-1)
     *              <  a1*2^m - a1/2^ell * p
     *              <  p*q+r-a0-a1/2^ell * p
     *         q-q0 >  (a0-r)/p + a1/2^ell
     *              >  -1.5   since a0>0, r<p, and a1>-2^(ell-1).
     *              >= -1     since q and q0 are integers.
     * So that q-(q0-1) >= 0.
     *
     * Note that because we have -2^ell < q0 < 0, then q0-1 is properly
     * represented by the unsigned machine word 2^ell-1+q0.
     *
     * Also, we have p*q0 >= p*a1*I/2^ell-p
     *                    >= a1*2^m-p
     *                    >= a-a0-p
     *         a-p*(q0-1) <= a0 + 2p < 4p
     *
     * To compute a-p*(q0-1), we actually compute
     * a-p*(q0+2^ell-1)+2^ell*p, which is a submul_ui followed by one
     * addition.
     */

    /* this reduces a in place, and copies the result to r */
    private:
    template<typename X>
        mp_limb_t const * preinv() const
        requires std::is_same_v<X, elt_ur_for_add>
        {
            return prime_preinv_for_add;
        }

    template<typename X>
        mp_limb_t const * preinv() const
        requires std::is_same_v<X, elt_ur_for_addmul>
        {
            return prime_preinv_for_addmul.data();
        }

    public:
    /* We're only going to emit a reduce function for our
     * two statically defined ur types, and only for these
     *
     * TODO: in fact, it does make sense to reduce from elt to
     * elt (e.g. when we generate a random bit string, or
     * possibly even when we read from a file)
     *
     */

    /* The code below is for anything that has more than N+1 words.
     * elt_ur_for_add could, conceivably, be N+2 words, and then we would
     * use this code. Our latest modifications have hardcoded it to be
     * N+1 though.
     */
    template<typename X>
        void reduce(elt& r, X& a) const
        requires (X::classifier >= 2)
        {
            T const* tx = static_cast<T const*>(this);
            auto j = preinv<X>();
            auto N = nlimbs<elt>();
            auto nx = overhead_limbs<X>();
            mp_limb_t tmp[nx + 1];
            if (preinv_shift) {
                mpn_lshift(tmp, a + N - 1, nx + 1, preinv_shift);
            } else {
                mpn_copyi(tmp + 1, a + N, nx);
            }
            mp_limb_t a1I[2 * nx];
            mpn_mul_n(a1I, tmp + 1, j, nx);
            mpn_add_n(a1I + nx, a1I + nx, tmp + 1, nx);
            mp_limb_t* q0 = a1I + nx;
            std::make_signed_t<mp_limb_t> sa1 = (tmp + 1)[nx - 1];
            if (sa1 < 0) {
                mpn_sub_n(q0, q0, j, nx);
                mpn_sub_1(q0, q0, nx, 1);
                mpn_add_n(a + nx, a + nx, prime_limbs(), N);
            }
            /* emulate a submul_n ; need to do mul first, then sub... */
            mp_limb_t scratch[N + nx];
            /* nx is N+1, so it's larger than N. However, if we use
             * X = elt_ur_for_add, then nx can be less than N (this is under
             * multiple ifs, since presently we've constrained elt_ur_for_add
             * to be N+1 and thus not use this code.
             */
            // mpn_mul_caller<N, nx>()(scratch, prime_limbs(), q0);
            mpn_mul(scratch, q0, nx, prime_limbs(), N);

            mpn_sub_n(a, a, scratch, N + nx);
#if !defined(NDEBUG) && !defined(DEBUG_INFINITE_LOOPS)
            int spin = 0;
#endif
            while (!upperlimbs_are_zero(a) || mpn_cmp(a, prime_limbs(), N) >= 0) {
                mp_limb_t cy = mpn_sub_n(a, a, prime_limbs(), N);
                tx->propagate_borrow(a, cy);
#if !defined(NDEBUG) && !defined(DEBUG_INFINITE_LOOPS)
                spin++;
                ASSERT_ALWAYS(spin < 4);
#endif
            }
            mpn_copyi(r, a, N);
        }
    template<typename X>
        void reduce(elt& r, X& a) const
        requires (X::classifier == 1)
        {
            T const* tx = static_cast<T const*>(this);
            auto N = nlimbs<elt>();
            auto j = preinv<X>();
            mp_limb_t a1 = a[N] << preinv_shift;
            if (preinv_shift) {
                a1 |= a[N - 1] >> (mp_bits_per_limb - preinv_shift);
            }
            const std::make_signed_t<mp_limb_t> sa1 = a1;
            mp_limb_t tmp[2];
#ifdef umul_ppmm
            umul_ppmm(tmp[1], tmp[0], a1, j[0]);
#elif defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
            __asm__("mulq %3" : "=a"(tmp[0]), "=d"(tmp[1]) : "0"(a1), "rm"(j[0]));
#else
            mpn_mul_n(tmp, &a1, j, 1);
#endif
            mp_limb_t q0 = tmp[1] + a1;
            if (sa1 < 0) {
                /* see above for the specificities of the negative case */
                q0 -= j[0] + 1;
                mpn_add_n(a + 1, a + 1, prime_limbs(), N);
            }
            mp_limb_t cy;
            cy = mpn_submul_1(a, prime_limbs(), N, q0);
            tx->propagate_borrow(a, cy);
#if !defined(NDEBUG) && !defined(DEBUG_INFINITE_LOOPS)
            int spin = 0;
#endif
            while (a[N] || mpn_cmp(a, prime_limbs(), N) >= 0) {
                cy = mpn_sub_n(a, a, prime_limbs(), N);
                tx->propagate_borrow(a, cy);
#if !defined(NDEBUG) && !defined(DEBUG_INFINITE_LOOPS)
                spin++;
                ASSERT_ALWAYS(spin < 4);
#endif
            }
            mpn_copyi(r, a, N);
        }
    /* We do need something that normalizes in place. This is
     * used every now and then (random generation is an example).
     *
     * It's not entirely clear to me if we can use the code above
     * directly or not. The solution below is perhaps a bit
     * unsatisfactory, but it's simple enough.  And it is almost
     * surely identical in terms of performance.
     */
    void reduce(elt& a) const
    {
        ARITH_MODP_TEMPORARY_ALLOC(this, elt_ur_for_add, r);
        T const* tx = static_cast<T const*>(this);
        tx->set(r, a);
        tx->reduce(a, r);
    }
    /*}}}*/
    /*}}}*/
    /*{{{ add_and_reduce */
    void add_and_reduce(elt& a, elt const& x) const
    {
        add_and_reduce(a, a, x);
    }
    void add_and_reduce(elt& a, elt const& x, elt const& y) const
    {
        T const* tx = static_cast<T const*>(this);
        auto N = nlimbs<elt>();
        ARITH_MODP_TEMPORARY_ALLOC(this, elt_ur_for_add, b);
        tx->add(b, x, y);
        if (b[N] || mpn_cmp(b, prime_limbs(), N) >= 0)
            mpn_sub_n(a, b, prime_limbs(), N);
        else
            mpn_copyi(a, b, N);
    }
    /*}}}*/
    /*{{{ sub_and_reduce */
    void sub_and_reduce(elt& a, elt const& x) const
    {
        sub_and_reduce(a, a, x);
    }
    void sub_and_reduce(elt& a, elt const& x, elt const& y) const
    {
        T const* tx = static_cast<T const*>(this);
        auto N = nlimbs<elt>();
        ARITH_MODP_TEMPORARY_ALLOC(this, elt_ur_for_add, b);
        tx->sub(b, x, y);
        if (b[N])
            mpn_add_n(a, b, prime_limbs(), N);
        else
            mpn_copyi(a, b, N);
    }
    /*}}}*/
    /*{{{ I/O */
    std::ostream& cxx_out(std::ostream& o, elt const& x) const
    {
        auto N = nlimbs<elt>();
        ARITH_MODP_TEMPORARY_ALLOC(this, elt, c);
        set(c, x);
        cxx_mpz cz;
        MPZ_SET_MPN(cz, (mp_limb_t const*)c, N);
        return o << cz;
    }
    std::ostream& write(std::ostream& o, elt const& x) const
    {
        return o.write((const char*)x.pointer(), elt_stride());
    }
    int fread(FILE* f, elt& x) const
    {
        auto N = nlimbs<elt>();
        int ret = ::fread((char*)x.pointer(), elt_stride(), 1, f);
        if (ret < 1)
            return ret;
        if (mpn_cmp(x, prime_limbs(), N) >= 0) {
            reduce(x);
        }
        return ret;
    }
    int fscan(FILE* f, elt& x) const
    {
        cxx_mpz xz;
        const size_t ret = mpz_inp_str(xz, f, 10);
        if (!ret)
            return 0;
        if (xz < 0 || xz >= p) {
            mpz_mod(xz, xz, p);
        }
        set(x, xz);
        return runtime_numeric_cast<int>(ret);
    }
    /*}}}*/

    void
        mul_ur(elt_ur_for_addmul& w, elt const& u, elt const& v) const
        {
            auto N = nlimbs<elt>();
            mpn_mul_n(w, u, v, N);
            /* elt_ur_for_addmul should be exactly 2N+1 limbs */
            w[2 * N] = 0;
        }

    void mul(elt& w, elt const& u, elt const& v) const
    {
        /* Should we use an mpn<2*N> instead? Would our
         * preinv<N+1> work ?
         */
        T const* tx = static_cast<T const*>(this);
        ARITH_MODP_TEMPORARY_ALLOC(this, elt_ur_for_addmul, t);
        tx->mul_ur(t, u, v);
        tx->reduce(w, t);
    }

    void
        addmul_ur(elt_ur_for_addmul & w, elt const& u, elt const& v) const
        {
            T const* tx = static_cast<T const*>(this);
            /* We don't have a convenient default at the mpn level */
            ARITH_MODP_TEMPORARY_ALLOC(this, elt_ur_for_addmul, t);
            tx->mul_ur(t, u, v);
            tx->add(w, t);
        }
    void addmul(elt& w, elt const& u, elt const& v) const
    {
        T const* tx = static_cast<T const*>(this);
        ARITH_MODP_TEMPORARY_ALLOC(this, elt_ur_for_addmul, t);
        tx->mul_ur(t, u, v);
        tx->add(t, w);
        tx->reduce(w, t);
    }
    void addmul_and_reduce(elt& w, elt const& u, elt const& v) const
    {
        T const* tx = static_cast<T const*>(this);
        ARITH_MODP_TEMPORARY_ALLOC(this, elt_ur_for_addmul, t);
        tx->set(t, w);
        tx->addmul_ur(t, u, v);
        tx->reduce(w, t);
    }

    /*{{{ accessors inside vectors */
    template<typename X>
        X* vec_subvec(X* p, ssize_t k) const
        requires (X::classifier >= 0)
        {
            ASSERT_ALWAYS(k <= std::numeric_limits<ssize_t>::max() / (ssize_t) nlimbs<X>());
            ASSERT_ALWAYS(k >= std::numeric_limits<ssize_t>::min() / (ssize_t) nlimbs<X>());
            return reinterpret_cast<X*>(p->pointer() + k * nlimbs<X>());
        }

    template<typename X>
        X& vec_item(X* p, size_t k) const
        requires (X::classifier >= 0)
        {
            return *vec_subvec(p, k);
        }

    template<typename X>
        X const* vec_subvec(X const* p, ssize_t k) const
        requires (X::classifier >= 0)
        {
            ASSERT_ALWAYS(k <= std::numeric_limits<ssize_t>::max() / (ssize_t) nlimbs<X>());
            ASSERT_ALWAYS(k >= std::numeric_limits<ssize_t>::min() / (ssize_t) nlimbs<X>());
            return reinterpret_cast<X const*>(p->pointer() + k * nlimbs<X>());
        }

    template<typename X>
        X const& vec_item(X const* p, size_t k) const
        requires (X::classifier >= 0)
        {
            return *vec_subvec(p, k);
        }
    /*}}}*/
    /* {{{ predicates on vectors */
    template<typename X>
        int vec_cmp(X const* a, X const* b, size_t k) const
        requires (X::classifier >= 0)
        {
            T const* tx = static_cast<T const*>(this);
            for (size_t i = 0; i < k; i++) {
                int r = tx->cmp(vec_item(a, i), vec_item(b, i));
                if (r)
                    return r;
            }
            return 0;
        }
    template<typename X>
        bool vec_is_zero(X const* p, size_t n) const
        requires (X::classifier >= 0)
        {
            T const* tx = static_cast<T const*>(this);
            for (size_t i = 0; i < n; i++)
                if (!tx->is_zero(vec_item(p, i)))
                    return false;
            return true;
        }
    /* }}} */
    /* {{{ assignments on vectors */
    template<typename X>
        void vec_set_zero(X* p, size_t n) const
        requires (X::classifier >= 0)
        {
            ASSERT_ALWAYS(!n || nlimbs<X>() <= SIZE_MAX / n);
            std::fill_n(p->pointer(), n * nlimbs<X>(), 0);
        }

    template<typename X>
        void vec_set(X* q, X const* p, size_t n) const
        requires (X::classifier >= 0)
        {
            ASSERT_ALWAYS(!n || nlimbs<X>() <= SIZE_MAX / n);
            size_t nn = n * nlimbs<X>();
            ASSERT_ALWAYS(nn <= SIZE_MAX / sizeof(X));
            if (q < p)
                std::copy_n(p->pointer(), nn, q->pointer());
            else
                std::copy_backward(p->pointer(), p->pointer() + nn, q->pointer() + nn);
        }

    /* extension */
    template<typename X>
        void vec_set(X* q, elt const* p, size_t n) const
        requires (X::classifier >= 1)
        {
            T const* tx = static_cast<T const*>(this);
            for (size_t i = 0; i < n; i++)
                tx->set(vec_item(q, i), vec_item(p, i));
        }

    template<typename X>
        void vec_set_random(X* p, size_t k, cxx_gmp_randstate & rstate) const
        requires (X::classifier >= 0)
        {
            T const* tx = static_cast<T const*>(this);
            for (size_t i = 0; i < k; ++i)
                tx->set_random(vec_item(p, i), rstate);
        }
    /* }}} */
    /*{{{ simd*/
    void simd_set_ui_at(elt& p, size_t, int v) const
    {
        T const* tx = static_cast<T const*>(this);
        tx->set(p, v);
    }
    void simd_add_ui_at(elt& p, size_t, int v) const
    {
        auto N = nlimbs<elt>();
        if (v > 0)
            mpn_add_1(p, p, N, v);
        else
            mpn_sub_1(p, p, N, -v);
    }
    int simd_test_ui_at(elt const& p, size_t) const
    {
        T const* tx = static_cast<T const*>(this);
        return !tx->is_zero(p);
    }
    int simd_hamming_weight(elt const& p) const { return !is_zero(p); }
    /*}}}*/

    int inverse(elt& res, elt const& x) const
    {
        T const* tx = static_cast<T const*>(this);
        if (tx->is_zero(x)) {
            tx->set_zero(res);
            return 0;
        }

        auto N = nlimbs<elt>();

        ASSERT_ALWAYS(mpn_cmp(x, prime_limbs(), N) < 0);
        ARITH_MODP_TEMPORARY_ALLOC(this, elt, g);
        ARITH_MODP_TEMPORARY_ALLOC(this, elt, u);
        ARITH_MODP_TEMPORARY_ALLOC(this, elt, v);
        mp_limb_t s[N + 1];
        tx->set(u, x);
        tx->set(v, p);

        mp_size_t sn = N + 1;
        mp_size_t gn = mpn_gcdext(g, s, &sn, u, N, v, N);
        if (gn != 1 || g[0] != 1) {
            set(res, g);
            return 0;
        }
        if (sn < 0) {
            sn = -sn;
            for (; (size_t) sn < N; s[sn++] = 0)
                ;
            mpn_sub_n(s, prime_limbs(), s, N);
        } else {
            for (; (size_t) sn < N; s[sn++] = 0)
                ;
        }
        mpn_copyi(res, s, N);
        return 1;
    }

    /*{{{ vec addition */
    template<typename X>
        void vec_add(X* q, elt const* p, size_t n) const
        requires (X::classifier >= 1)
        {
            T const* tx = static_cast<T const*>(this);
            for (size_t i = 0; i < n; i++)
                tx->add(vec_item(q, i), vec_item(p, i));
        }
    template<typename X>
        void vec_add(X* q, X* p, size_t n) const
        requires (X::classifier >= 1)
        {
            T const* tx = static_cast<T const*>(this);
            for (size_t i = 0; i < n; i++)
                tx->add(vec_item(q, i), vec_item(p, i));
        }
    void vec_add_and_reduce(elt* q,
            elt const* a,
            elt const* b,
            size_t n) const
    {
        T const* tx = static_cast<T const*>(this);
        for (size_t i = 0; i < n; i++)
            tx->add_and_reduce(
                    vec_item(q, i), vec_item(a, i), vec_item(b, i));
    }
    void vec_add_and_reduce(elt* q, elt const* a, size_t n) const
    {
        vec_add_and_reduce(q, q, a, n);
    }
    /*}}}*/
    /*{{{ vec subtraction */
    template<typename X>
        void vec_sub(X* q, elt const* p, size_t n) const
        requires (X::classifier >= 1)
        {
            T const* tx = static_cast<T const*>(this);
            for (size_t i = 0; i < n; i++)
                tx->sub(vec_item(q, i), vec_item(p, i));
        }
    template<typename X>
        void vec_sub(X* q, X* p, size_t n) const
        requires (X::classifier >= 1)
        {
            T const* tx = static_cast<T const*>(this);
            for (size_t i = 0; i < n; i++)
                tx->sub(vec_item(q, i), vec_item(p, i));
        }
    void vec_sub_and_reduce(elt* q,
            elt const* a,
            elt const* b,
            size_t n) const
    {
        T const* tx = static_cast<T const*>(this);
        for (size_t i = 0; i < n; i++)
            tx->sub_and_reduce(
                    vec_item(q, i),
                    vec_item(a, i),
                    vec_item(b, i));
    }
    /*}}}*/

    /*{{{ vec negation*/
    void vec_neg(elt* q, size_t n) const
    {
        T const* tx = static_cast<T const*>(this);
        for (size_t i = 0; i < n; i++)
            tx->neg(vec_item(q, i), vec_item(q, i));
    }
    void vec_neg(elt* q, elt const* p, size_t n) const
    {
        T const* tx = static_cast<T const*>(this);
        /* used in matpoly_sub */
        for (size_t i = 0; i < n; i++)
            tx->neg(vec_item(q, i), vec_item(p, i));
    }
    template<typename X>
        void vec_neg(X* q, size_t n) const
        requires (X::classifier >= 1)
        {
            T const* tx = static_cast<T const*>(this);
            for (size_t i = 0; i < n; i++)
                tx->neg(vec_item(q, i), vec_item(q, i));
        }
    template<typename X>
        void vec_neg(X* q, X const* p, size_t n) const
        requires (X::classifier >= 1)
        {
            T const* tx = static_cast<T const*>(this);
            /* used in matpoly_sub */
            for (size_t i = 0; i < n; i++)
                tx->neg(vec_item(q, i), vec_item(p, i));
        }
    /*}}}*/

    /*{{{ vec reduction */
    /* Note that we depend on reduce() being available */
    template<typename X>
        void vec_reduce(elt* q, X* p, size_t n) const
        requires (X::classifier >= 0)
        {
            T const* tx = static_cast<T const*>(this);
            for (size_t i = 0; i < n; i++)
                tx->reduce(vec_item(q, i), vec_item(p, i));
        }
    /*}}}*/

    /*{{{ vec simd operations*/
    void vec_simd_set_ui_at(elt* p, size_t k, int v) const
    {
        T const* tx = static_cast<T const*>(this);
        tx->simd_set_ui_at(vec_item(p, k), 0, v);
    }
    void vec_simd_add_ui_at(elt* p, size_t k, int v) const
    {
        T const* tx = static_cast<T const*>(this);
        tx->simd_add_ui_at(vec_item(p, k), 0, v);
    }
    size_t vec_simd_hamming_weight(elt const* p, size_t n) const
    {
        T const* tx = static_cast<T const*>(this);
        size_t r = 0;
        for (size_t i = 0; i < n; ++i)
            r += tx->simd_hamming_weight(vec_item(p, i));
        return r;
    }
    size_t vec_simd_find_first_set(elt& x, elt const* p, size_t n) const
    {
        T const* tx = static_cast<T const*>(this);
        static_assert(T::simd_groupsize() == 1, "this code wants trivial simd");
        for (size_t i = 0; i < n; ++i) {
            if (!tx->is_zero(tx->set(x, vec_item(p, i)))) {
                return i;
            }
        }
        return SIZE_MAX;
    }

    /*}}}*/

    private:
    void vec_conv_ur_ks(elt_ur_for_addmul* w,
            elt const* u,
            size_t n,
            elt const* v,
            size_t m) const /*{{{*/
    {
        T const* tx = static_cast<T const*>(this);
        auto N = nlimbs<elt>();
        /* The only thing that makes this a member function and
         * not a class function is the fact that we use the
         * instances's p as a max bound on the coefficient
         * values...
         *
         * (and of course, in the variable width case, we depend on the
         * width of p which exists only at the instance level)
         */
        // compute base as a power 2^GMP_NUMB_BITS
        // This is the least number of words that can accomodate
        //     log_2( (p-1)^2 * min(n,m) )
        cxx_mpz q;
        mpz_sub_ui(q, p, 1);
        mpz_mul(q, q, q);
        mpz_mul_ui(q, q, std::min(m, n));

        size_t nbits = mpz_sizeinbase(q, 2);
        size_t nwords = 1 + ((nbits - 1) / GMP_NUMB_BITS);
        nbits = GMP_NUMB_BITS * nwords;

        cxx_mpz W;
        {
            // Create big integers. We don't use the fine-tuning
            // from mpz_limbs_finish, since this would require
            // annoying checks when we read back from W.

            cxx_mpz U; /*{{{*/
            mp_limb_t* pU = mpz_limbs_write(U, n * nwords);
            mpn_zero(pU, n * nwords);
            for (size_t i = 0; i < n; ++i)
                mpn_copyi(pU + i * nwords, vec_item(u, i), N);
            // mpz_limbs_finish(U, n*nwords);
            SIZ(U) = n * nwords;
            /*}}}*/
            cxx_mpz V; /*{{{*/
            mp_limb_t* pV = mpz_limbs_write(V, m * nwords);
            mpn_zero(pV, m * nwords);
            for (size_t i = 0; i < m; ++i)
                mpn_copyi(pV + i * nwords, vec_item(v, i), N);
            SIZ(V) = m * nwords;
            /*}}}*/

            // Multiply
            mpz_mul(W, U, V);
        }

        // Put coefficients in w
        tx->vec_set_zero(w, m + n - 1);

        assert(mpz_size(W) >= (size_t)((m + n - 1) * nwords));
        const mp_limb_t* pW = mpz_limbs_read(W);
        ASSERT_ALWAYS(nwords <= nlimbs<elt_ur_for_addmul>());
        for (size_t i = 0; i < m + n - 1; ++i)
            mpn_copyi(vec_item(w, i), pW + i * nwords, nwords);
    }                                                   /*}}}*/
    void vec_conv_ur_n(elt_ur_for_addmul* w,
            elt const* u,
            elt const* v,
            size_t n) const /*{{{*/
    {
        T const* tx = static_cast<T const*>(this);
        if (n == 0)
            return;
        if (n == 1) {
            tx->mul_ur(vec_item(w, 0), vec_item(u, 0), vec_item(v, 0));
            return;
        }
        if (n == 2) { // Kara 2
            ARITH_MODP_TEMPORARY_ALLOC(this, elt, t1);
            ARITH_MODP_TEMPORARY_ALLOC(this, elt, t2);
            tx->mul_ur(vec_item(w, 0), vec_item(u, 0), vec_item(v, 0));
            tx->mul_ur(vec_item(w, 2), vec_item(u, 1), vec_item(v, 1));
            tx->add_and_reduce(t1, vec_item(u, 0), vec_item(u, 1));
            tx->add_and_reduce(t2, vec_item(v, 0), vec_item(v, 1));
            tx->mul_ur(vec_item(w, 1), t1, t2);
            tx->sub(vec_item(w, 1), vec_item(w, 0));
            tx->sub(vec_item(w, 1), vec_item(w, 2));
            return;
        }
        if (n == 3) { // do it in 6
            ARITH_MODP_TEMPORARY_ALLOC(this, elt, t1);
            ARITH_MODP_TEMPORARY_ALLOC(this, elt, t2);
            ARITH_MODP_TEMPORARY_ALLOC(this, elt_ur_for_addmul, s);
            // a0*b0*(1 - X)
            tx->mul_ur(vec_item(w, 0), vec_item(u, 0), vec_item(v, 0));
            tx->neg(vec_item(w, 1), vec_item(w, 0));
            // a1*b1*(-X + 2*X^2 - X^3)
            tx->mul_ur(vec_item(w, 2), vec_item(u, 1), vec_item(v, 1));
            tx->neg(vec_item(w, 3), vec_item(w, 2));
            tx->add(vec_item(w, 2), vec_item(w, 2));
            tx->add(vec_item(w, 1), vec_item(w, 3));
            // a2*b2*(-X^3+X^4)
            tx->mul_ur(vec_item(w, 4), vec_item(u, 2), vec_item(v, 2));
            tx->sub(vec_item(w, 3), vec_item(w, 4));
            // (a0+a1)*(b0+b1)*(X - X^2)
            tx->add_and_reduce(t1, vec_item(u, 0), vec_item(u, 1));
            tx->add_and_reduce(t2, vec_item(v, 0), vec_item(v, 1));
            tx->mul_ur(s, t1, t2);
            tx->add(vec_item(w, 1), s);
            tx->sub(vec_item(w, 2), s);
            // (a1+a2)*(b1+b2)*(X^3 - X^2)
            tx->add_and_reduce(t1, vec_item(u, 1), vec_item(u, 2));
            tx->add_and_reduce(t2, vec_item(v, 1), vec_item(v, 2));
            tx->mul_ur(s, t1, t2);
            tx->add(vec_item(w, 3), s);
            tx->sub(vec_item(w, 2), s);
            // (a0+a1+a2)*(b0+b1+b2)* X^2
            tx->add_and_reduce(t1, vec_item(u, 0), t1);
            tx->add_and_reduce(t2, vec_item(v, 0), t2);
            tx->mul_ur(s, t1, t2);
            tx->add(vec_item(w, 2), s);
            return;
        }

        // generic Kara
        size_t n0, n1;
        n0 = n / 2;
        n1 = n - n0;
        tx->vec_conv_ur_n(w, u, v, n0);
        tx->vec_conv_ur_n(vec_subvec(w, 2 * n0), vec_subvec(u, n0), vec_subvec(v, n0), n1);
        tx->set_zero(vec_item(w, 2 * n0 - 1));

        /* This always goes to the heap */
        elt* tmpu = tx->alloc(n1);
        elt* tmpv = tx->alloc(n1);
        auto* tmpw = tx->template alloc<elt_ur_for_addmul>(2 * n1 - 1);

        tx->vec_set(tmpu, u, n0);
        if (n1 != n0)
            tx->set_zero(vec_item(tmpu, n0));
        tx->vec_add_and_reduce(tmpu, vec_subvec(u, n0), n1);

        tx->vec_set(tmpv, v, n0);
        if (n1 != n0)
            tx->set_zero(vec_item(tmpv, n0));
        tx->vec_add_and_reduce(tmpv, vec_subvec(v, n0), n1);

        tx->vec_conv_ur_n(tmpw, tmpu, tmpv, n1);

        tx->vec_sub(tmpw, w, 2 * n0 - 1);
        tx->vec_sub(tmpw, vec_subvec(w, 2 * n0), 2 * n1 - 1);
        tx->vec_add(vec_subvec(w, n0), tmpw, 2 * n1 - 1);

        tx->free(tmpu);
        tx->free(tmpv);
        tx->template free<elt_ur_for_addmul>(tmpw);
    } /*}}}*/

    public:
    void vec_conv(elt* w,
            elt const* u,
            size_t n,
            elt const* v,
            size_t m) const /*{{{*/
    {
        // This should be good for pz
        T const* tx = static_cast<T const*>(this);
        auto* tmp = tx->template alloc<elt_ur_for_addmul>(m + n - 1);
        tx->vec_conv_ur(tmp, u, n, v, m);
        tx->vec_reduce(w, tmp, m + n - 1);
        tx->template free<elt_ur_for_addmul>(tmp);
    } /*}}}*/
    void vec_conv_ur(elt_ur_for_addmul* w,
            elt const* u,
            size_t n,
            elt const* v,
            size_t m) const
    { /*{{{*/
        T const* tx = static_cast<T const*>(this);
        if ((n > 1) && (m > 1) && (n + m > 15)) {
            tx->vec_conv_ur_ks(w, u, n, v, m);
            return;
        }
        if (n == m) {
            tx->vec_conv_ur_n(w, u, v, n);
            return;
        }
        ARITH_MODP_TEMPORARY_ALLOC(this, elt_ur_for_addmul, acc);
        if (n > m) {
            std::swap(u, v);
            std::swap(n, m);
        }
        for (size_t k = 0; k < n; ++k) {
            tx->set_zero(acc);
            for (size_t i = 0; i <= k; ++i)
                tx->addmul_ur(acc, vec_item(u, i), vec_item(v, k - i));
            tx->set(vec_item(w, k), acc);
        }
        for (size_t k = n; k < m; ++k) {
            tx->set_zero(acc);
            for (size_t i = 0; i < n; ++i)
                tx->addmul_ur(acc, vec_item(u, i), vec_item(v, k - i));
            tx->set(vec_item(w, k), acc);
        }
        for (size_t k = m; k < n + m - 1; ++k) {
            tx->set_zero(acc);
            for (size_t i = k - m + 1; i < n; ++i)
                tx->addmul_ur(acc, vec_item(u, i), vec_item(v, k - i));
            tx->set(vec_item(w, k), acc);
        }
    } /*}}}*/

    void vec_add_dotprod(elt& w, elt const* u, elt const* v, size_t n) const
    {
        T const* tx = static_cast<T const*>(this);
        ARITH_MODP_TEMPORARY_ALLOC(this, elt_ur_for_addmul, t);
        tx->set(t, w);
        for (size_t i = 0; i < n; ++i)
            tx->addmul_ur(t, vec_item(u, i), vec_item(v, i));
        tx->reduce(w, t);
    }

    void vec_addmul_and_reduce(elt* w,
            elt const* u,
            elt const& v,
            size_t n) const
    {
        T const* tx = static_cast<T const*>(this);
        for (size_t i = 0; i < n; ++i)
            tx->addmul_and_reduce(vec_item(w, i), vec_item(u, i), v);
    }
};

template<int n>
struct gfp : public gfp_base<n, gfp<n>>
{
    using super = gfp_base<n, gfp<n>>;
    template<typename... Args>
        gfp(Args&&... args)
        : super { std::forward<Args>(args)... }
    {}
};

/* Now for some sizes, we see a clear interest in using auxiliary
 * vector types. We call these "fast" types. The general compromise
 * is that we accept types which may be a little wider, but generally
 * allow for better performance. The specs go typically as follows.
 *
 * - conversion to the "fast" types must be done for both operands
 *   (say, source vector as well as destination vector). We don't
 *   intend to go with the same kind of arithmetic that what we do
 *   with elt and elt_ur above, where a "mixed" add/sub function
 *   exists.
 *
 * - "fast" types are amenable to vector instructions
 *
 * - up to some number of additions or subtractions may be performed
 *   on the fast type before reduction.
 *
 * - type may be ambiguous (so comparison entails conversion).
 *
 * We have two natural choices:
 *
 *  - RNS representation
 *  - carry-save (aka nails).
 *
 * The specializations below work with nails. The idea is as follows.
 * For a p spanning three 64-bit words, we spread data into four
 * 48-bit words in an avx register. Then we can accumulate up to 2^16
 * of these at little cost.
 */

/* the default version is just making no difference, so that we'll
 * use the simple elt / elt_ur mechanism */
template<typename T>
struct fast_type : public T
{};

}

#endif /* ARITH_MODP_MAIN_HPP_ */
