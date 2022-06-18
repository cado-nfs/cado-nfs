#ifndef ARITH_MODP_HPP_
#define ARITH_MODP_HPP_

#include <utility>
#include <algorithm>
#include <type_traits>

#include <gmp.h>

#include "gmp-hacks.h"
#include "gmp_aux.h"
#include "macros.h"
#include "memory.h"
#include "cxx_mpz.hpp"
#include "fmt/format.h"

#include "arith-concrete-base.hpp"

#if defined(HAVE_AVX2) || defined(HAVE_SSSE3)
#include <x86intrin.h>
#endif

/* the elt objects are always reduced between 0 and p-1, so that it is
 * fine to have p as large as N*sizeof(mp_limb_t).
 *
 * However the unreduced elements are implicitly signed, and the leading
 * bit is here to reflect that. The reduction function takes this into
 * account. All subtraction, negation functions on unreduced elements are
 * therefore avoiding all normalizations.
 */

#define  xxxDEBUG_INFINITE_LOOPS

namespace arith_modp {
namespace details {
    template<int N> struct mpn {
        static const int alignment = sizeof(mp_limb_t);
        static constexpr const int nlimbs = N;
        static constexpr const bool is_flat_storage = true;
        typedef mpn<N> self;
        mp_limb_t x[N];
        mpn() { mpn_zero(x, N); }
        mpn(mpn const& a) = default;
        mpn& operator=(mpn const& a) = default;
        mpn(mpz_srcptr a) { MPN_SET_MPZ(x, N, a); }
        self& operator=(mpz_srcptr a) { MPN_SET_MPZ(x, N, a); return *this; }
        self& operator=(unsigned long a) { zero(); x[0] = a; return *this; }
        mpn(unsigned long a) { zero(); x[0] = a; }
        void zero() { mpn_zero(x, N); }
        template<int na>
            mpn(mpn<na> const & a,
                    typename std::enable_if<(na < N)>::type * = 0)
            {
                mpn_copyi(x, a.x, na);
                mpn_zero(x + na, N - na);
            }
        template<int na>
            typename std::enable_if<(na < N), mpn>::type &
            operator=(mpn<na> const & a)
            {
                mpn_copyi(x, a.x, na);
                mpn_zero(x + na, N - na);
                return *this;
            }

        operator mp_limb_t * () { return x; }
        operator const mp_limb_t * () const { return x; }
        // mp_limb_t & operator[](long k) { return x[k]; }
        // mp_limb_t const & operator[](long k) const { return x[k]; }

        static void zero(self * x, int n) {
            // see bug 21663
            mpn_zero(x->x, n * N);
        }
        static void copy(self * y, const self * x, int n) {
            // see bug 21663
            mpn_copyi(y->x, x->x, n * N);
        }
        int cmp(self const& a) const {
            return mpn_cmp(x, a.x, N);
        }
        int cmp(unsigned long a) const {
            int r = (a < x[0]) - (x[0] < a);
            if (r) return r;
            for(int i = 1 ; i < N ; i++)
                if (x[i]) return 1;
            return 0;
        }
        bool operator==(self const& a) const { return cmp(a) == 0; }
        bool operator<(self const& a) const { return cmp(a) < 0; }
        bool is_zero() const {
            for(int i = 0 ; i < N ; i++) if (x[i]) return false;
            return true;
        }
        /*
        bool operator<(self const& a) {
            return memcmp(x, a.x, N * sizeof(mp_limb_t)) < 0;
        }
        bool operator>(self const& a) {
            return memcmp(x, a.x, N * sizeof(mp_limb_t)) > 0;
        }
        */
    };

    template<int N, typename T>
        struct gfp_base : public arith_concrete_base {
            static std::string impl_name() {
                return fmt::format(FMT_STRING("p{}"), N);
            }
            static constexpr const bool is_characteristic_two = false;
            static constexpr const bool simd_groupsize = 1;
            struct elt : public arith_concrete_base::elt, public mpn<N> {
                template<typename... Args> elt(Args&&... args) : mpn<N>(std::forward<Args>(args)...) {}
            };
            typedef mpn<N+1> elt_ur_for_add;
            typedef mpn<2*N+1> elt_ur_for_addmul;
            static constexpr const int nlimbs = N;
        private:
            cxx_mpz p;
            /* How many bits to he have to shift p left so that its
             * highest bit is set?
             */
            int preinv_shift=0;
            mpn<elt_ur_for_add::nlimbs - N> prime_preinv_for_add;
            mpn<elt_ur_for_addmul::nlimbs - N> prime_preinv_for_addmul;
        public:
            gfp_base() : p(0) {}
            gfp_base(gfp_base const &) = default;
            gfp_base(gfp_base &&) = default;
            gfp_base& operator=(gfp_base const &) = default;
            gfp_base& operator=(gfp_base &&) = default;
            gfp_base(cxx_mpz const & p, unsigned int simd_groupsize) : p(p) {
                ASSERT_ALWAYS(simd_groupsize == 1);
                size_t m = mpz_sizeinbase(p, 2);
                ASSERT_ALWAYS(m <= (size_t) N * mp_bits_per_limb);
                preinv_shift = (mp_bits_per_limb - m) % mp_bits_per_limb;
                compute_preinv(prime_preinv_for_add);
                compute_preinv(prime_preinv_for_addmul);
            }
            mp_limb_t const * prime_limbs() const { return mpz_limbs_read(p); }
            // elt const& prime() const { return *reinterpret_cast<mpn<N> const *>(prime_limbs()); }

            /*{{{ element layout information */
            static inline size_t elt_stride() { return sizeof(elt); }
            static inline size_t vec_elt_stride(size_t s) {
                return s * elt_stride();
            }
            /*}}}*/

            mpz_srcptr characteristic() const { return p; }

            /*{{{ allocation / deallocation of (vectors of) elements */
            /* These allocation interfaces seem a bit stupid. At the low
             * hard level, we know that they're just the same as new[]
             * anyway.
             */

            /* This macro leverages enable_if and exposes a
             * specialization that knows about a type X (which should
             * probably appear in the arguments), and which could be
             * anything that is _wider_ than mpn<N> = elt. Note that the
             * default type of X is elt
             */
#define gfp_polymorphic_with_default(storage, ret_type)		        \
            template<typename X = elt>					\
                storage							\
                typename std::enable_if<(X::nlimbs >= N), ret_type>::type

#define gfp_polymorphic(storage, above, ret_type)			\
            template<typename X>					\
                storage							\
                typename std::enable_if<(X::nlimbs >= N + above), ret_type>::type

            gfp_polymorphic_with_default(static inline, X *)
                alloc(size_t k = 1) { return new X[k]; }
            gfp_polymorphic_with_default(static inline, void)
                free(X * u) { delete[] u; }
            gfp_polymorphic_with_default(static inline, X *)
                realloc(X * u, size_t k0, size_t k) {
                    X * v = alloc<X>(k);
                    std::copy_n(u, k0, v);
                    free<X>(u);
                    return v;
                }

            /*}}}*/
            /*{{{ predicates */
            gfp_polymorphic(static inline, 0, bool)
            is_zero(X const & x) { return x.is_zero(); }

            gfp_polymorphic(static inline, 0, int)
            cmp(X const & x, unsigned long a) { return x.cmp(a); }

            gfp_polymorphic(static inline, 0, int)
            cmp(X const & x, X const & y) { return x.cmp(y); }

            /*}}}*/
            /* {{{ assignments */
            gfp_polymorphic(static inline, 0, X &)
            set(X& x, unsigned long a) { return x = a; }

            gfp_polymorphic(static inline, 0, X &)
            set(X& x, X const & a) { return x = a; }

            gfp_polymorphic(static inline, 0, X &)
            set(X& x, cxx_mpz const & a) { return x = a; }

            gfp_polymorphic(static inline, 1, X &)
            set(X& x, elt const & a) { return x = a; }

            gfp_polymorphic(static inline, 0, void)
            set_zero(X& x) { x.zero(); }

            gfp_polymorphic(inline, 0, void)
            set_random(X& x, gmp_randstate_ptr rstate) const {
                cxx_mpz xz;
                mpz_urandomm(xz, rstate, p);
                x = xz;
            }
            gfp_polymorphic(static inline, 0, bool)
            upperlimbs_are_zero(X & a)
            {
                for(int i = X::nlimbs ; i-- > N ; )
                    if (a.x[i]) return false;
                return true;
            }

            static inline void stream_store(elt * dst, elt const& src) { *dst = src; }
            /* }}} */

            /*{{{ addition, at the element level (elt or wider) */
            gfp_polymorphic(static inline, 2, void)
            propagate_carry(X & a, mp_limb_t cy)
            {
                mpn_add_1(a + N, a + N, X::nlimbs - N, cy);
            }
            static inline void propagate_carry(mpn<N+1> & a, mp_limb_t cy)
            {
                a[N] += cy;
            }

            gfp_polymorphic(static inline, 1, void)
            add(X & dst, elt const & a, elt const & b)
            {
                mp_limb_t cy = mpn_add_n(dst, a, b, N);
                mpn_zero(dst + N, X::nlimbs-N);
                T::propagate_carry(dst, cy);
            }
            gfp_polymorphic(static inline, 1, void)
            add(X & dst, elt const & src)
            {
                mp_limb_t cy = mpn_add_n(dst, dst, src, N);
                T::propagate_carry(dst, cy);
            }
            /* This addition is only for unreduced types. These types are
             * always considered wide enough so that overflows work.
             */
            gfp_polymorphic(static inline, 1, void)
            add_ur(X & dst, X const & src)
            {
                //mpn_add_n(dst, dst, src, X::nlimbs);
                T::add(dst, src);
            }
            gfp_polymorphic(static inline, 1, void)
            add(X & dst, X const & src)
            {
                mpn_add_n(dst, dst, src, X::nlimbs);
            }

            /* an addmul is ok for to go for an unreduced type which is
             * still somewhat narrow (only one extra limb).
             */
            gfp_polymorphic(static inline, 1, void)
            addmul_ui(X & dst, elt const & src, mp_limb_t x)
            {
                mp_limb_t cy = mpn_addmul_1(dst, src, N, x);
                T::propagate_carry(dst, cy);
            }
            /*}}}*/

            /*{{{ subtraction, at the element level (elt or wider) */
            gfp_polymorphic(static inline, 2, void)
            propagate_borrow(X & a, mp_limb_t cy)
            {
                mpn_sub_1(a + N, a + N, X::nlimbs - N, cy);
            }
            static inline void propagate_borrow(mpn<N+1> & a, mp_limb_t cy)
            {
                a[N] -= cy;
            }

            gfp_polymorphic(static inline, 1, void)
            sub(X & dst, elt const & a, elt const & b)
            {
                mp_limb_t cy = mpn_sub_n(dst, a, b, N);
                mpn_zero(dst + N, X::nlimbs-N);
                T::propagate_borrow(dst, cy);
            }
            gfp_polymorphic(static inline, 1, void)
            sub(X & dst, elt const & src)
            {
                mp_limb_t cy = mpn_sub_n(dst, dst, src, N);
                T::propagate_borrow(dst, cy);
            }
            /* This subition is only for unreduced types. These types are
             * always considered wide enough so that overflows work.
             */
            gfp_polymorphic(static inline, 1, void)
            sub_ur(X & dst, X const & src)
            {
                //mpn_sub_n(dst, dst, src, X::nlimbs);
                T::sub(dst, src);
            }
            gfp_polymorphic(static inline, 1, void)
            sub(X & dst, X const & src)
            {
                mpn_sub_n(dst, dst, src, X::nlimbs);
            }

            /* a submul is ok for to go for an unreduced type which is
             * still somewhat narrow (only one extra limb).
             */
            gfp_polymorphic(static inline, 1, void)
            submul_ui(X & dst, elt const & src, mp_limb_t x)
            {
                mp_limb_t cy = mpn_submul_1(dst, src, N, x);
                T::propagate_borrow(dst, cy);
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
            inline void
            neg(elt & dst, elt const & src) const
            {
                T const * tx = static_cast<T const *>(this);
                mpn_neg(dst, src, N);
                if (!dst.is_zero())
                    mpn_add_n(dst, dst, tx->prime_limbs(), N);
            }

            gfp_polymorphic(static inline, 1, void)
            neg(X & dst, X const & src)
            {
                mpn_neg(dst, src, X::nlimbs);
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
            template<int e>
            void compute_preinv(mpn<e> & j)
            {
                cxx_mpz big;
                mpz_set_ui(big,1);
                size_t m = mpz_sizeinbase(p, 2);
                size_t ell = e * mp_bits_per_limb;
                mpz_mul_2exp(big, big, m+ell);
                mpz_fdiv_q(big, big, p);
                ASSERT_ALWAYS(mpz_sizeinbase(big, 2) == (ell + 1));
                mpz_fdiv_r_2exp(big, big, ell);
                MPN_SET_MPZ(j, e, big);
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
             * I<2^(ell+1), we obtaib q0 > -2^ell. Therefore q0 is well
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
                inline typename std::enable_if<
                std::is_same<X, elt_ur_for_add>::value,
                decltype(prime_preinv_for_add)>::type const &
                    preinv() const { return prime_preinv_for_add; }

            template<typename X>
                inline typename std::enable_if<
                std::is_same<X, elt_ur_for_addmul>::value,
                decltype(prime_preinv_for_addmul)>::type const &
                    preinv() const { return prime_preinv_for_addmul; }
        public:
            /* We're only going to emit a reduce function for our
             * two statically defined ur types.
             *
             * TODO: in fact, it does make sense to reduce from elt to
             * elt (e.g. when we generate a random bit string, or
             * possibly even when we read from a file)
             *
             */
            template<typename X>
                typename std::enable_if<
                (std::is_same<X, elt_ur_for_add>::value ||
                std::is_same<X, elt_ur_for_addmul>::value) &&
                (X::nlimbs > 1 + N)
                >::type
            reduce(elt & r, X & a) const
            {
                auto j = preinv<X>();
                constexpr int extra = X::nlimbs - N;
                mp_limb_t tmp[extra + 1];
                if (preinv_shift) {
                    mpn_lshift(tmp, a + N - 1, extra + 1, preinv_shift);
                } else {
                    mpn_copyi(tmp + 1, a + N, extra);
                }
                mp_limb_t a1I[2*extra];
                mpn_mul_n(a1I, tmp + 1, j, extra);
                mpn_add_n(a1I + extra, a1I + extra, tmp + 1, extra);
                mp_limb_t * q0 = a1I + extra;
                typename std::make_signed<mp_limb_t>::type sa1 = (tmp+1)[extra-1];
                if (sa1 < 0) {
                    mpn_sub_n(q0, q0, j, extra);
                    mpn_sub_1(q0, q0, extra, 1);
                    mpn_add_n(a + extra, a + extra, prime_limbs(), N);
                }
                /* emulate a submul_n ; need to do mul first, then sub... */
                mp_limb_t scratch[N + extra];
                mpn_mul(scratch, prime_limbs(), N, q0, extra);
                mpn_sub_n(a, a, scratch, N + extra);
#if !defined(NDEBUG) && !defined(DEBUG_INFINITE_LOOPS)
                int spin=0;
#endif
                while (!upperlimbs_are_zero(a) || mpn_cmp(a, prime_limbs(), N) >= 0) {
                    mp_limb_t cy = mpn_sub_n(a, a, prime_limbs(), N);
                    T::propagate_borrow(a, cy);
#if !defined(NDEBUG) && !defined(DEBUG_INFINITE_LOOPS)
                    spin++;
                    ASSERT_ALWAYS(spin < 4);
#endif
                }
                mpn_copyi(r, a, N);
            }
            template<typename X>
                typename std::enable_if<
                (std::is_same<X, elt_ur_for_add>::value ||
                std::is_same<X, elt_ur_for_addmul>::value) &&
                (X::nlimbs == 1 + N)
                >::type
            reduce(elt & r, X & a) const
            {
                auto j = preinv<X>();
                mp_limb_t a1 = a[N] << preinv_shift;
                if (preinv_shift) {
                    a1 |= a[N-1] >> (mp_bits_per_limb - preinv_shift);
                }
                typedef std::make_signed<mp_limb_t>::type signed_mp_limb_t;
                signed_mp_limb_t sa1 = a1;
                mp_limb_t tmp[2];
#ifdef  umul_ppmm
                umul_ppmm(tmp[1], tmp[0], a1, j[0]);
#elif defined(HAVE_GCC_STYLE_AMD64_INLINE_ASM)
                __asm__ ("mulq %3" : "=a"(tmp[0]), "=d" (tmp[1]) : "0" (a1), "rm" (j[0]));
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
                T::propagate_borrow(a, cy);
#if !defined(NDEBUG) && !defined(DEBUG_INFINITE_LOOPS)
                int spin=0;
#endif
                while (a[N] || mpn_cmp(a, prime_limbs(), N) >= 0) {
                    cy = mpn_sub_n(a, a, prime_limbs(), N);
                    T::propagate_borrow(a, cy);
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
             * unsatisfactory, but it's simple ebough.  And it is almost
             * surely identical in terms of performance.
             */
            void reduce(elt& a) const {
                elt_ur_for_add r = a;
                reduce(a, r);
            }
            /*}}}*/
            /*}}}*/
            /*{{{ add_and_reduce */
            inline void add_and_reduce(elt & a, elt const & x) const {
                add_and_reduce(a, a, x);
            }
            inline void add_and_reduce(elt & a, elt const & x, elt const & y) const {
                elt_ur_for_add b;
                T::add(b, x, y);
                if (b[N] || mpn_cmp(b, prime_limbs(), N) >= 0)
                    mpn_sub_n(a, b, prime_limbs(), N);
                else
                    mpn_copyi(a, b, N);
            }
            /*}}}*/
            /*{{{ sub_and_reduce */
            inline void sub_and_reduce(elt & a, elt const & x) const {
                sub_and_reduce(a, a, x);
            }
            inline void sub_and_reduce(elt & a, elt const & x, elt const & y) const {
                elt_ur_for_add b;
                T::sub(b, x, y);
                if (b[N])
                    mpn_add_n(a, b, prime_limbs(), N);
                else
                    mpn_copyi(a, b, N);
            }
            /*}}}*/
            /*{{{ I/O */
            static inline std::ostream& cxx_out(std::ostream& o, elt const & x) {
                elt c = x;
                cxx_mpz cz;
                MPZ_SET_MPN(cz, (mp_limb_t const *) c, N);
                return o << cz;
            }
            static inline std::ostream& write(std::ostream& o, elt const & x) {
                return o.write((const char *) x.x, elt_stride());
            }
            inline int fread(FILE * f, elt & x) const {
                int ret = ::fread((char *) x.x, elt_stride(), 1, f);
                if (ret < 1) return ret;
                if (mpn_cmp(x, prime_limbs(), N) >= 0) {
                    reduce(x);
                }
                return ret;
            }
            inline int fscan(FILE * f, elt& x) const {
                cxx_mpz xz;
                size_t ret = mpz_inp_str(xz, f, 10);
                if (!ret) return 0;
                if (xz < 0 || xz >= p) {
                    mpz_mod(xz, xz, p);
                }
                x = xz;
                return ret;
            }
            /*}}}*/

            gfp_polymorphic(static inline, N, void)
                mul_ur(X & w, elt const & u, elt const & v)
            {
                mpn_mul_n(w, u, v, N);
                if (X::nlimbs > 2 * N) mpn_zero(w + N + N, X::nlimbs - N - N);
            }

            void mul(elt & w, elt const & u, elt const & v) const
            {
                /* Should we use an mpn<2*N> instead? Would our
                 * preinv<N+1> work ?
                 */
                mpn<2 * N + 1> t;
                T::mul_ur(t, u, v);
                static_cast<T const *>(this)->reduce(w, t);
            }

            gfp_polymorphic(static inline, N, void)
                addmul_ur(X & w, elt const & u, elt const & v)
            {
                /* We don't have a convenient default at the mpn level
                 */
                elt_ur_for_addmul t;
                T::mul_ur(t, u, v);
                T::add(w, t);
            }
            void addmul(elt & w, elt const & u, elt const & v) const
            {
                mpn<2 * N + 1> t;
                T::mul_ur(t, u, v);
                T::add(t, w);
                static_cast<T const *>(this)->reduce(w, t);
            }
            void addmul_and_reduce(elt & w, elt const & u, elt const & v) const
            {
                T const * tx = static_cast<T const *>(this);
                elt_ur_for_addmul t;
                tx->set(t, w);
                T::addmul_ur(t, u, v);
                tx->reduce(w, t);
            }

            /*{{{ accessors inside vectors */
            /* Likewise, this seems a bit stupid. THe only reason to have
             * them is if we want to have a pz layer. But it's very
             * doubtful that we can manage to have that with flat layout.
             */
            gfp_polymorphic(static inline, 0, X *)
            vec_subvec(X * p, size_t k) { return p + k; }

            gfp_polymorphic(static inline, 0, X &)
            vec_item(X * p, size_t k) { return p[k]; }

            gfp_polymorphic(static inline, 0, X const *)
            vec_subvec(X const * p, size_t k) { return p + k; }

            gfp_polymorphic(static inline, 0, X const &)
            vec_item(X const * p, size_t k) { return p[k]; }
            /*}}}*/
            /* {{{ predicates on vectors */
            gfp_polymorphic(static inline, 0, int)
            vec_cmp(X const * a, X const * b, size_t k)
            {
                for(size_t i = 0 ; i < k ; i++) {
                    int r = T::cmp(a[i], b[i]);
                    if (r) return r;
                }
                return 0;
            }
            gfp_polymorphic(static inline, 0, bool)
            vec_is_zero(X const * p, size_t n) {
                for(size_t i = 0 ; i < n ; i++)
                    if (!p[i].is_zero()) return false;
                return true;
            }
            /* }}} */
            /* {{{ assignments on vectors */
            gfp_polymorphic(static inline, 0, void)
            vec_set_zero(X * p, size_t n) {
                static_assert(X::is_flat_storage, "X must be flat");
                std::fill_n(&p[0][0], n * X::nlimbs, 0);
            }

            gfp_polymorphic(static inline, 0, void)
            vec_set(X * q, X const * p, size_t n)
            {
                static_assert(X::is_flat_storage, "X must be flat");
                if (q < p)
                    std::copy_n(p, n, q);
                else
                    std::copy_backward(p, p + n, q + n);
            }

            /* extension */
            gfp_polymorphic(static inline, 1, void)
            vec_set(X * q, elt const * p, size_t n)
            {
                for(size_t i = 0 ; i < n ; i++)
                    T::set(q[i], p[i]);
            }

            gfp_polymorphic(inline, 0, void)
            vec_set_random(X * p, size_t k, gmp_randstate_ptr rstate) const {
                for(size_t i = 0 ; i < k ; ++i) 
                    static_cast<const T *>(this)->set_random(p[i], rstate);
            }
            /* }}} */
                /*{{{ simd*/
                static inline void
                    simd_set_ui_at(elt & p, size_t, int v)
                    {
                        p = v;
                    }
                static inline void
                    simd_add_ui_at(elt & p, size_t, int v)
                    {
                        if (v > 0)
                            mpn_add_1(p, p, N, v);
                        else
                            mpn_sub_1(p, p, N, -v);
                    }
                static inline int
                    simd_test_ui_at(elt const & p, size_t)
                    {
                        return !p.is_zero();
                    }
                static inline int
                    simd_hamming_weight(elt const & p)
                    {
                        return !p.is_zero();
                    }
                /*}}}*/

            int inverse(elt & res, elt const & x) const
            {
                if (x.is_zero()) {
                    res.zero();
                    return 0;
                }

                ASSERT_ALWAYS (mpn_cmp(x, prime_limbs(), N) < 0);
                elt g,u,v;
                mpn<N+1> s;
                u = x;
                v = p;

                mp_size_t sn = N+1;
                mp_size_t gn = mpn_gcdext(g, s, &sn, u, N, v, N);
                if (gn != 1 || g[0] != 1) {
                    res = g;
                    return 0;
                }
                if (sn < 0) {
                    sn = -sn;
                    for( ; sn < N; s[sn++] = 0);
                    mpn_sub_n(s, prime_limbs(), s, N);
                } else {
                    for( ; sn < N; s[sn++] = 0);
                }
                mpn_copyi(res, s, N);
                return 1;
            }

            /*{{{ vec addition */
            gfp_polymorphic(static inline, 1, void)
            vec_add(X * q, elt const * p, size_t n) {
                for(size_t i = 0 ; i < n ; i++)
                    T::add(q[i], p[i]);
            }
            gfp_polymorphic(static inline, 1, void)
            vec_add(X * q, X * p, size_t n) {
                for(size_t i = 0 ; i < n ; i++)
                    T::add(q[i], p[i]);
            }
            inline void vec_add_and_reduce(elt * q, elt const * a, elt const * b, size_t n) const {
                T const * tx = static_cast<T const *>(this);
                for(size_t i = 0 ; i < n ; i++)
                    tx->add_and_reduce(tx->vec_item(q, i), tx->vec_item(a, i), tx->vec_item(b, i));
            }
            inline void vec_add_and_reduce(elt * q, elt const * a, size_t n) const {
                vec_add_and_reduce(q, q, a, n);
            }
            /*}}}*/
            /*{{{ vec subtraction */
            gfp_polymorphic(static inline, 1, void)
            vec_sub(X * q, elt const * p, size_t n) {
                for(size_t i = 0 ; i < n ; i++)
                    T::sub(q[i], p[i]);
            }
            gfp_polymorphic(static inline, 1, void)
            vec_sub(X * q, X * p, size_t n) {
                for(size_t i = 0 ; i < n ; i++)
                    T::sub(q[i], p[i]);
            }
            inline void vec_sub_and_reduce(elt * q, elt const * a, elt const * b, size_t n) const {
                for(size_t i = 0 ; i < n ; i++)
                    static_cast<const T *>(this)->sub_and_reduce(q[i], a[i], b[i]);
            }
            /*}}}*/
            /*{{{ vec negation*/
            inline void
            vec_neg(elt * q, size_t n) const {
                T const * tx = static_cast<T const *>(this);
                for(size_t i = 0 ; i < n ; i++)
                    tx->neg(q[i], q[i]);
            }
            inline void
            vec_neg(elt * q, elt const * p, size_t n) const {
                T const * tx = static_cast<T const *>(this);
                /* used in matpoly_sub */
                for(size_t i = 0 ; i < n ; i++)
                    tx->neg(q[i], p[i]);
            }
            gfp_polymorphic(static inline, 1, void)
            vec_neg(X * q, size_t n) {
                for(size_t i = 0 ; i < n ; i++)
                    T::neg(q[i], q[i]);
            }
            gfp_polymorphic(static inline, 1, void)
            vec_neg(X * q, X const * p, size_t n) {
                /* used in matpoly_sub */
                for(size_t i = 0 ; i < n ; i++)
                    T::neg(q[i], p[i]);
            }
            /*}}}*/

            /*{{{ vec reduction */
            /* Note that we depend on reduce() being available */
            gfp_polymorphic(inline, 0, void)
            vec_reduce(elt * q, X * p, size_t n) const {
                for(size_t i = 0 ; i < n ; i++)
                    static_cast<const T *>(this)->reduce(q[i], p[i]);
            }
            /*}}}*/
                /*{{{ vec simd operations*/
                static inline void
                    vec_simd_set_ui_at(elt * p, size_t k, int v)
                    {
                        T::simd_set_ui_at(T::vec_item(p, k), 0, v);
                    }
                static inline void
                    vec_simd_add_ui_at(elt * p, size_t k, int v)
                    {
                        T::simd_add_ui_at(T::vec_item(p, k), 0, v);
                    }
                static inline int
                    vec_simd_hamming_weight(elt const * p, size_t n)
                    {
                        int r = 0;
                        for(size_t i = 0 ; i < n ; ++i)
                            r += T::simd_hamming_weight(T::vec_item(p, i));
                        return r;
                    }
                static inline int
                    vec_simd_find_first_set(elt & x, elt const * p, size_t n)
                    {
                        static_assert(T::simd_groupsize == 1, "this code wants trivial simd");
                        for(size_t i = 0 ; i < n ; ++i) {
                            if (!T::is_zero(T::set(x, T::vec_item(p, i)))) {
                                return i;
                            }
                        }
                        return -1;
                    }

                /*}}}*/

        private:
            void vec_conv_ur_ks(elt_ur_for_addmul * w, elt const * u, size_t n, elt const * v, size_t m) const/*{{{*/
            {
                /* The only thing that makes this a member function and
                 * not a class function is the fact that we use the
                 * instances's p as a max bound on the coefficient
                 * values...
                 */
                // compute base as a power 2^GMP_NUMB_BITS
                // This is the least number of words that can accomodate
                //     log_2( (p-1)^2 * min(n,m) )
                cxx_mpz q;
                mpz_sub_ui(q, p, 1);
                mpz_mul(q, q, q);
                mpz_mul_ui(q, q, std::min(m, n));
    
                long nbits = mpz_sizeinbase(q, 2);
                long nwords = 1 + ((nbits-1) / GMP_NUMB_BITS);
                nbits = GMP_NUMB_BITS*nwords;

                cxx_mpz W;
                {
                    // Create big integers. We don't use the fine-tuning
                    // from mpz_limbs_finish, since this would require
                    // annoying checks when we read back from W.

                    cxx_mpz U;/*{{{*/
                    mp_limb_t * pU = mpz_limbs_write(U, n*nwords);
                    mpn_zero(pU, n*nwords);
                    for (size_t i = 0; i < n; ++i)
                        mpn_copyi(pU + i*nwords, u[i], N);
                    // mpz_limbs_finish(U, n*nwords);
                    SIZ(U) = n*nwords;
                    /*}}}*/
                    cxx_mpz V;/*{{{*/
                    mp_limb_t * pV = mpz_limbs_write(V, m*nwords);
                    mpn_zero(pV, m*nwords);
                    for (size_t i = 0; i < m; ++i)
                        mpn_copyi(pV + i*nwords, v[i], N);
                    SIZ(V) = m*nwords;
                    /*}}}*/

                    // Multiply
                    mpz_mul(W, U, V);
                }
    
                // Put coefficients in w
                T::vec_set_zero(w, m + n - 1);

                assert (mpz_size(W) >= (size_t) ((m+n-1)*nwords));
                const mp_limb_t * pW = mpz_limbs_read(W);
                ASSERT_ALWAYS(nwords <= elt_ur_for_addmul::nlimbs);
                for (size_t i = 0; i < m+n-1; ++i) 
                    mpn_copyi(w[i].x, pW + i * nwords, nwords);
            }/*}}}*/
            inline void vec_conv_ur_n(elt_ur_for_addmul * w, elt const * u, elt const * v, size_t n) const /*{{{*/
            {
                T const * tx = static_cast<T const *>(this);
                if (n == 0)
                    return;
                if (n == 1) {
                    tx->mul_ur(w[0], u[0], v[0]);
                    return;
                }
                if (n == 2) {  // Kara 2
                    elt t1, t2;
                    tx->mul_ur(w[0], u[0], v[0]);
                    tx->mul_ur(w[2], u[1], v[1]);
                    tx->add_and_reduce(t1, u[0], u[1]);
                    tx->add_and_reduce(t2, v[0], v[1]);
                    tx->mul_ur(w[1], t1, t2);
                    tx->sub(w[1], w[0]);
                    tx->sub(w[1], w[2]);
                    return;
                }
                if (n == 3) {  // do it in 6
                    elt t1, t2;
                    elt_ur_for_addmul s;
                    // a0*b0*(1 - X)
                    tx->mul_ur(w[0], u[0], v[0]);
                    tx->neg(w[1], w[0]);
                    // a1*b1*(-X + 2*X^2 - X^3)
                    tx->mul_ur(w[2], u[1], v[1]);
                    tx->neg(w[3], w[2]);
                    tx->add(w[2], w[2]);
                    tx->add(w[1], w[3]);
                    // a2*b2*(-X^3+X^4)
                    tx->mul_ur(w[4], u[2], v[2]);
                    tx->sub(w[3], w[4]);
                    // (a0+a1)*(b0+b1)*(X - X^2)
                    tx->add_and_reduce(t1, u[0], u[1]);
                    tx->add_and_reduce(t2, v[0], v[1]);
                    tx->mul_ur(s, t1, t2);
                    tx->add(w[1], s);
                    tx->sub(w[2], s);
                    // (a1+a2)*(b1+b2)*(X^3 - X^2)
                    tx->add_and_reduce(t1, u[1], u[2]);
                    tx->add_and_reduce(t2, v[1], v[2]);
                    tx->mul_ur(s, t1, t2);
                    tx->add(w[3], s);
                    tx->sub(w[2], s);
                    // (a0+a1+a2)*(b0+b1+b2)* X^2
                    tx->add_and_reduce(t1, u[0], t1);
                    tx->add_and_reduce(t2, v[0], t2);
                    tx->mul_ur(s, t1, t2);
                    tx->add(w[2], s);
                    return;
                }

                // generic Kara
                size_t n0, n1;
                n0 = n / 2;
                n1 = n - n0;
                tx->vec_conv_ur_n(w, u, v, n0);
                tx->vec_conv_ur_n(w + 2*n0, u + n0, v + n0, n1);
                tx->set_zero(w[2*n0-1]);

                elt * tmpu = tx->alloc(n1);
                elt * tmpv = tx->alloc(n1);
                auto * tmpw = tx->template alloc<elt_ur_for_addmul>(2*n1-1);

                tx->vec_set(tmpu, u, n0);
                if (n1 != n0) tx->set_zero(tmpu[n0]);
                tx->vec_add_and_reduce(tmpu, u+n0, n1);

                tx->vec_set(tmpv, v, n0);
                if (n1 != n0) tx->set_zero(tmpv[n0]);
                tx->vec_add_and_reduce(tmpv, v+n0, n1);

                tx->vec_conv_ur_n(tmpw, tmpu, tmpv, n1);

                tx->vec_sub(tmpw, w, 2*n0-1);
                tx->vec_sub(tmpw, w + 2*n0, 2*n1-1);
                tx->vec_add(w + n0, tmpw, 2*n1-1);

                return;
            }/*}}}*/

        public:
            inline void vec_conv(elt * w, elt const * u, size_t n, elt const * v, size_t m) const/*{{{*/
            {
                auto * tmp = T::template alloc<elt_ur_for_addmul>(m + n - 1);
                T const * tx = static_cast<T const *>(this);
                tx->vec_conv_ur(tmp, u, n, v, m);
                tx->vec_reduce(w, tmp, m+n-1);
                T::free(tmp);
            }/*}}}*/
            inline void vec_conv_ur(elt_ur_for_addmul * w, elt const * u, size_t n, elt const * v, size_t m) const {/*{{{*/
                T const * tx = static_cast<T const *>(this);
                if ((n > 1) && (m > 1) && (n+m > 15)) {
                    tx->vec_conv_ur_ks(w, u, n, v, m);
                    return;
                }
                if (n == m) {
                    tx->vec_conv_ur_n(w, u, v, n);
                    return;
                }
                elt_ur_for_addmul acc;
                if (n > m) {
                    std::swap(u, v);
                    std::swap(n, m);
                }
                for(size_t k = 0; k < n; ++k) {
                    acc.zero();
                    for(size_t i = 0; i <= k; ++i)
                        T::addmul_ur(acc, u[i], v[k-i]);
                    w[k] = acc;
                }
                for(size_t k = n; k < m; ++k) {
                    acc.zero();
                    for(size_t i = 0; i < n; ++i)
                        T::addmul_ur(acc, u[i], v[k-i]);
                    w[k] = acc;
                }
                for(size_t k = m; k < n+m-1; ++k) {
                    acc.zero();
                    for(size_t i = k-m+1; i < n; ++i)
                        T::addmul_ur(acc, u[i], v[k-i]);
                    w[k] = acc;
                }
            }/*}}}*/

            void vec_add_dotprod(elt & w, elt const * u, elt const * v, size_t n) const
            {
                T const * tx = static_cast<T const *>(this);
                elt_ur_for_addmul t;
                t = w;
                for(size_t i = 0; i < n; ++i)
                    T::addmul_ur(t, tx->vec_item(u, i), tx->vec_item(v, i));
                tx->reduce(w, t);
            }

            void vec_addmul_and_reduce(elt * w, elt const * u, elt const & v, size_t n) const
            {
                T const * tx = static_cast<T const *>(this);
                for(size_t i = 0; i < n; ++i)
                    tx->addmul_and_reduce(tx->vec_item(w, i), tx->vec_item(u, i), v);
            }
        };


    template<int n>
        struct gfp : public gfp_base<n, gfp<n> >
    {
        typedef gfp_base<n, gfp<n> > super;
        template<typename... Args> gfp(Args&&... args) : super(std::forward<Args>(args)...) {}
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
    template<typename T> struct fast_type : public T { };


#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    /* Now some specializations */

    /* {{{ gfp<1,1> */
    template<> struct gfp<1> : public gfp_base<1,gfp<1> > {
        typedef gfp_base<1, gfp<1> > super;
        using typename super::elt;
        using typename super::elt_ur_for_add;
        /* We need this so that the templates at the base level
         * participate in the resolution here.
         */
        using super::add;
        using super::sub;

        template<typename... Args> gfp(Args&&... args) : super(std::forward<Args>(args)...) {}

        static inline void add_ur(elt_ur_for_add & dst, elt_ur_for_add const & src)
        {
            asm("# gfp<1, 1>::add\n"
                "addq %q2, %q0\n"
                "adcq %q3, %q1\n"
                : "+r"(dst[0]), "+r"(dst[1])
                : "rm"(src[0]), "rm"(src[1])
               );
        }

        static inline void add(elt_ur_for_add & dst, elt_ur_for_add const & src)
        {
            add_ur(dst, src);
        }

        static inline void add(elt_ur_for_add & dst, elt const & src)
        {
            asm("# gfp<1, 1>::add\n"
                "addq %q2, %q0\n"
                "adcq $0x0, %q1\n"
                : "+r"(dst[0]), "+r"(dst[1])
                : "rm"(src[0])
               );
        }
        static inline void add(elt_ur_for_add & dst, elt const & a, elt const & b)
        {
            asm("# gfp<1, 1>::add\n"
                "movq %q2, %q0\n"
                "xorq %q1, %q1\n"
                "addq %q3, %q0\n"
                "adcq $0x0, %q1\n"
                : "=&r"(dst[0]), "=&r"(dst[1])
                : "rm"(a[0]), "rm"(b[0])
               );
        }

        static inline void sub(elt_ur_for_add & dst, elt const & src)
        {
            asm("# gfp<1, 1>::sub\n"
                "subq %q2, %q0\n"
                "sbbq $0x0, %q1\n"
                : "+r"(dst[0]), "+r"(dst[1])
                : "rm"(src[0])
               );
        }

        static inline void sub_ur(elt_ur_for_add & dst, elt_ur_for_add const & src)
        {
            asm("# gfp<1, 1>::sub\n"
                "subq %q2, %q0\n"
                "sbbq %q3, %q1\n"
                : "+r"(dst[0]), "+r"(dst[1])
                : "rm"(src[0]), "rm"(src[1])
               );
        }
        static inline void sub(elt_ur_for_add & dst, elt_ur_for_add const & src)
        { sub_ur(dst, src); }

        static inline void sub(elt_ur_for_add & dst, elt const & a, elt const & b)
        {
            asm("# gfp<1, 1>::sub\n"
                "movq %q2, %q0\n"
                "xorq %q1, %q1\n"
                "subq %q3, %q0\n"
                "sbbq $0x0, %q1\n"
                : "=&r"(dst[0]), "=&r"(dst[1])
                : "rm"(a[0]), "rm"(b[0])
               );
        }

        static inline void addmul_ui(elt_ur_for_add & dst, elt const & src, mp_limb_t x)
        {
            mp_limb_t foo, bar;
            asm("# gfp<1, 1>::addmul_ui\n"
                "mulq   %[mult]\n"
                "addq   %%rax, %[z0]\n"
                "adcq   $0, %%rdx\n"
                "addq   %%rdx, %[z1]\n"
            : "=a"(foo), "=&d"(bar), [z0]"+rm"(dst[0]), [z1]"+rm"(dst[1])
            : "0"(src[0]), [mult]"r1m"(x)
            );
        }
        static inline void submul_ui(elt_ur_for_add & dst, elt const & src, mp_limb_t x)
        {
            mp_limb_t foo, bar;
            asm("# gfp<1, 1>::submul_ui\n"
                "mulq   %[mult]\n"
                "subq   %%rax, %[z0]\n"
                "adcq   $0, %%rdx\n"
                "subq   %%rdx, %[z1]\n"
            : "=a"(foo), "=&d"(bar), [z0]"+rm"(dst[0]), [z1]"+rm"(dst[1])
            : "0"(src[0]), [mult]"r1m"(x)
            );
        }
    };
    /* }}} */

    /* {{{ gfp<2,1> */
    template<> struct gfp<2> : public gfp_base<2,gfp<2> > {
        typedef gfp_base<2, gfp<2> > super;
        using typename super::elt;
        using typename super::elt_ur_for_add;
        /* We need this so that the templates at the base level
         * participate in the resolution here.
         */
        using super::add;
        using super::sub;
        template<typename... Args> gfp(Args&&... args) : super(std::forward<Args>(args)...) {}
        static inline void add(elt_ur_for_add & dst, elt const & src) {
            asm("# gfp<2, 1>::add\n"
                "addq %q[s0], %q[d0]\n"
                "adcq %q[s1], %q[d1]\n"
                "adcq $0x0, %q[d2]\n"
                : [d0]"+rm"(dst[0]), [d1]"+rm"(dst[1]), [d2]"+rm"(dst[2])
                : [s0]"r"(src[0]), [s1]"r"(src[1])
               );
        }

        static inline void add(elt_ur_for_add & dst, elt_ur_for_add const & src) {
            asm("# gfp<2, 1>::add\n"
                "addq %q[s0], %q[d0]\n"
                "adcq %q[s1], %q[d1]\n"
                "adcq %q[s2], %q[d2]\n"
                : [d0]"+rm"(dst[0]), [d1]"+rm"(dst[1]), [d2]"+rm"(dst[2])
                : [s0]"r"(src[0]), [s1]"r"(src[1]), [s2]"r"(src[2])
               );
        }

        static inline void sub(elt_ur_for_add & dst, elt const & src) {
            asm("# gfp<2, 1>::sub\n"
                "subq %q[s0], %q[d0]\n"
                "sbbq %q[s1], %q[d1]\n"
                "sbbq $0x0, %q[d2]\n"
                : [d0]"+rm"(dst[0]), [d1]"+rm"(dst[1]), [d2]"+rm"(dst[2])
                : [s0]"r"(src[0]), [s1]"r"(src[1])
               );
        }
        static inline void sub(elt_ur_for_add & dst, elt_ur_for_add const & src) {
            asm("# gfp<2, 1>::sub\n"
                "subq %q[s0], %q[d0]\n"
                "sbbq %q[s1], %q[d1]\n"
                "sbbq %q[s2], %q[d2]\n"
                : [d0]"+rm"(dst[0]), [d1]"+rm"(dst[1]), [d2]"+rm"(dst[2])
                : [s0]"r"(src[0]), [s1]"r"(src[1]), [s2]"r"(src[2])
               );
        }
        static inline void add(elt_ur_for_add & dst, elt const & a, elt const & b) {
            dst = a; add(dst, b);
        }
        static inline void sub(elt_ur_for_add & dst, elt const & a, elt const & b) {
            dst = a; sub(dst, b);
        }

        static inline void addmul_ui(elt_ur_for_add & dst, elt const & src, mp_limb_t x)
        {
            mp_limb_t foo, bar;
            asm("# gfp<2, 1>::addmul_ui\n"
                "mulq   %[mult]\n"
                "addq   %%rax, %[z0]\n"
                "adcq   $0, %%rdx\n"
                "movq   %%rdx, %%rcx\n"
                "movq   %[s1], %%rax\n"
                "mulq   %[mult]\n"
                "addq   %%rcx, %%rax\n"
                "adcq   $0, %%rdx\n"
                "addq   %%rax, %[z1]\n"
                "adcq   $0, %%rdx\n"
                "addq   %%rdx, %[z2]\n"
            : "=&a"(foo), "=&d"(bar),
            [z0]"+rm"(dst[0]),
            [z1]"+rm"(dst[1]),
            [z2]"+rm"(dst[2])
            : [s0]"0"(src[0]), [s1]"rm"(src[1]), [mult]"rm"(x)
            : "rcx"
            );
        }

        static inline void submul_ui(elt_ur_for_add & dst, elt const & src, mp_limb_t x)
        {
            mp_limb_t foo, bar;
            asm("# gfp<2, 1>::submul_ui\n"
                "mulq   %[mult]\n"
                "subq   %%rax, %[z0]\n"
                "adcq   $0, %%rdx\n"
                "movq   %%rdx, %%rcx\n"
                "movq   %[s1], %%rax\n"
                "mulq   %[mult]\n"
                "addq   %%rcx, %%rax\n"
                "adcq   $0, %%rdx\n"
                "subq   %%rax, %[z1]\n"
                "adcq   $0, %%rdx\n"
                "subq   %%rdx, %[z2]\n"
            : "=&a"(foo), "=&d"(bar),
            [z0]"+rm"(dst[0]),
            [z1]"+rm"(dst[1]),
            [z2]"+rm"(dst[2])
            : [s0]"0"(src[0]), [s1]"rm"(src[1]), [mult]"rm"(x)
            : "rcx"
            );
        }
    };
    /* }}} */

    /* {{{ macros for assembly for further specializations */
#define FEED_IN_WITH_S0_IN_RAX(in1, r0, r1) \
        /* status: s0 in rax */                                \
        "mulq   %[mult]\n"             /* rdx:rax = s0 * v */  \
        "movq   %%rax, %%" #r0 "\n"    /* lo contrib to d1 */  \
        "movq   " in1 ", %%rax\n"        /* load s1          */  \
        "movq   %%rdx, %%" #r1 "\n"    /* hi contrib to d1 */
#define FEED_IN(in0, in1, r0, r1) \
        "movq   " in0 ", %%rax\n"       \
        FEED_IN_WITH_S0_IN_RAX(in1, r0, r1)
#define INNER_MUL(op, out, in, r0, r1, r2)   \
        /* status: r1:r0 to be added to d_{i+1}:d_i, rax = s_{i+1} */     \
        "xorq   %%" #r2 ", %%" #r2 "\n"                                   \
        "mulq   %[mult]\n"                   /* rdx:rax = s_{i+1} * v */  \
        "" #op "q %%" #r0 ", " out "\n" /* store d_i             */   \
        "adcq   %%rax, %%" #r1 "\n"         /* lo contrib to d_{i+1} */   \
        "adcq   %%rdx, %%" #r2 "\n"         /* hi contrib to d_{i+2} */   \
        "movq   " in ", %%rax\n"       /* load s_{i+2}          */
#define FINISH(op, opc, out0, out1, out2, r0, r1) \
        /* r1:r0 to be added to d_{i+1}:d_i ; rax = s_{i+2} */	\
        "mulq   %[mult]\n"                   			\
        "" #op "q   %%" #r0 ", " out0 "\n"  			\
        "adcq   %%rax, %%" #r1 "\n"				\
        "adcq   $0x0, %%rdx\n"					\
        "" #op "q   %%" #r1 ", " out1 "\n" 			\
        "" #opc "q   %%rdx, " out2 "\n" 
    /* }}} */
    /* {{{ this macro actually exposes the specialization in itself */
#define EXPOSE_SPECIALIZATION(n)					\
    template<> struct gfp<n> : public gfp_base<n,gfp<n> > {		\
        typedef gfp_base<n, gfp<n> > super;                     	\
        using typename super::elt;					\
        using typename super::elt_ur_for_add;				\
        template<typename... Args> gfp(Args&&... args)                  \
        : super(std::forward<Args>(args)...) {}                         \
        using super::add_ur;                                            \
        using super::sub_ur;                                            \
        /* We need this so that the templates at the base level		\
         * participate in the resolution here.				\
         */								\
        using super::add;						\
        using super::sub;						\
        static inline void add(elt_ur_for_add & dst, elt const & src) {		\
            asm("# gfp<" #n ", 1>::add\n"				\
                    ADDSUB_CODE ## n(add, adc)				\
               );							\
        }								\
        static inline void sub(elt_ur_for_add & dst, elt const & src) {		\
            asm("# gfp<" #n ", 1>::sub\n"					\
                    ADDSUB_CODE ## n(sub, sbb)				\
               );							\
        }								\
        static inline void add(elt_ur_for_add & dst, elt const & a, elt const & b) { \
            dst = a; add(dst, b);                                       \
        }                                                               \
        static inline void sub(elt_ur_for_add & dst, elt const & a, elt const & b) { \
            dst = a; sub(dst, b);                                       \
        }                                                               \
        static inline void add(elt_ur_for_add & dst, elt_ur_for_add const & src)	\
        { add_ur(dst, src); }						\
        static inline void sub(elt_ur_for_add & dst, elt_ur_for_add const & src)	\
        { sub_ur(dst, src); }						\
        static inline void addmul_ui(elt_ur_for_add & dst, elt const & src, mp_limb_t x) { \
            mp_limb_t foo MAYBE_UNUSED;					\
            asm ("# gfp<" #n ", 1>::addmul_ui\n"				\
                    ADDSUBMUL_CODE ## n(add, adc)			\
            );								\
        }								\
        static inline void submul_ui(elt_ur_for_add & dst, elt const & src, mp_limb_t x) { \
            mp_limb_t foo MAYBE_UNUSED;					\
            asm("# gfp<" #n ", 1>::submul_ui\n"				\
                    ADDSUBMUL_CODE ## n(sub, sbb)			\
            );								\
        }								\
    }
    /* }}} */

    /* {{{ code for gfp<3, 1> */
#define ADDSUBMUL_CODE3(op, opc)					\
                FEED_IN_WITH_S0_IN_RAX("%[s1]", r8, r9)			\
                INNER_MUL(op, "%[z0]", "%[s2]", r8, r9, r10)		\
                FINISH(op, opc, "%[z1]", "%[z2]", "%[z3]", r9, r10)	\
                : "=&a"(foo),                                           \
                    [z0]"+rm"(dst[0]),				\
                    [z1]"+rm"(dst[1]),				\
                    [z2]"+rm"(dst[2]),				\
                    [z3]"+rm"(dst[3])					\
                :							\
                    [s0]"0"(src[0]),					\
                    [s1]"rm"(src[1]),					\
                    [s2]"rm"(src[2]),					\
                    [mult]"rm"(x)					\
                : "r8", "r9", "r10", "rdx"

#define ADDSUB_CODE3(op, opc)   \
        "" #op  "q %q[s0], %q[d0]\n"					\
        "" #opc "q %q[s1], %q[d1]\n"					\
        "" #opc "q %q[s2], %q[d2]\n"					\
        "" #opc "q $0x0, %q[d3]\n"				\
                :							\
                [d0]"+rm"(dst[0]),					\
                [d1]"+rm"(dst[1]),					\
                [d2]"+rm"(dst[2]),					\
                [d3]"+rm"(dst[3])					\
                :							\
                [s0]"r"(src[0]),					\
                [s1]"r"(src[1]),					\
                [s2]"r"(src[2])

    EXPOSE_SPECIALIZATION(3);
    /* }}} */

    /* {{{ code for gfp<4, 1> */
    /*
#define ADDSUBMUL_CODE4(op, opc)					\
                FEED_IN_WITH_S0_IN_RAX("%[s1]", r8, r9)			\
                INNER_MUL(op, "%[z0]", "%[s2]", r8, r9, r10)		\
                INNER_MUL(op, "%[z1]", "%[s3]", r9, r10, r11)		\
                FINISH(op, opc, "%[z2]", "%[z3]", "%[z4]", r10, r11)	\
                : "=&a"(foo),                                           \
                    [z0]"+rm"(dst[0]),				\
                    [z1]"+rm"(dst[1]),				\
                    [z2]"+rm"(dst[2]),				\
                    [z3]"+rm"(dst[3]),				\
                    [z4]"+rm"(dst[4])					\
                :							\
                    [s0]"0"(src[0]),					\
                    [s1]"rm"(src[1]),					\
                    [s2]"rm"(src[2]),					\
                    [s3]"rm"(src[3]),					\
                    [mult]"rm"(x)					\
                : "r8", "r9", "r10", "r11", "rdx"

#define xADDSUB_CODE4(op, opc)   \
"" #op  "q %q[s0], (%[z])\n"					\
"" #opc "q %q[s1], 0x8(%[z])\n"					\
"" #opc "q %q[s2], 0x10(%[z])\n"				\
"" #opc "q %q[s3], 0x18(%[z])\n"				\
"" #opc "q $0x0, 0x20(%[z])\n"					\
        :							\
        :							\
            [z]"r"(&dst[0]),				        \
            [s0]"r"(src[0]),					\
            [s1]"r"(src[1]),					\
            [s2]"r"(src[2]),					\
            [s3]"r"(src[3])                                   \
        : "memory"


                */
#define ADDSUB_CODE4(op, opc)   \
        "" #op  "q %q[s0], %q[d0]\n"					\
        "" #opc "q %q[s1], %q[d1]\n"					\
        "" #opc "q %q[s2], %q[d2]\n"					\
        "" #opc "q %q[s3], %q[d3]\n"					\
        "" #opc "q $0x0, %q[d4]\n"					\
                :							\
                [d0]"+rm"(dst[0]),					\
                [d1]"+rm"(dst[1]),					\
                [d2]"+rm"(dst[2]),					\
                [d3]"+rm"(dst[3]),					\
                [d4]"+rm"(dst[4])					\
                :							\
                [s0]"r"(src[0]),					\
                [s1]"r"(src[1]),					\
                [s2]"r"(src[2]),					\
                [s3]"r"(src[3])


#define ADDSUBMUL_CODE4(op, opc)					\
                FEED_IN("0x0(%[s])", "0x8(%[s])", r8, r9)		\
                INNER_MUL(op, "0x0(%[z])", "0x10(%[s])", r8, r9, r10)	\
                INNER_MUL(op, "0x8(%[z])", "0x18(%[s])", r9, r10, r11)	\
                FINISH(op, opc,                                         \
                        "0x10(%[z])", "0x18(%[z])", "0x20(%[z])",       \
                        r10, r11)                               	\
                :							\
                :							\
                    [z]"D"(&dst[0]),				        \
                    [s]"S"(&src[0]),					\
                    [mult]"rm"(x)					\
                : "r8", "r9", "r10", "r11", "rax", "rdx", "memory"

    EXPOSE_SPECIALIZATION(4);
    /* }}} */

    /* {{{ code for gfp<5, 1> */

#define ADDSUBMUL_CODE5(op, opc)					\
                FEED_IN("0x0(%[s])", "0x8(%[s])", r8, r9)		\
                INNER_MUL(op, "0x0(%[z])", "0x10(%[s])", r8, r9, r10)	\
                INNER_MUL(op, "0x8(%[z])", "0x18(%[s])", r9, r10, r11)	\
                INNER_MUL(op, "0x10(%[z])", "0x20(%[s])", r10, r11, r8)	\
                FINISH(op, opc,                                         \
                        "0x18(%[z])", "0x20(%[z])", "0x28(%[z])",       \
                        r11, r8)                               	\
                :							\
                :							\
                    [z]"D"(&dst[0]),				        \
                    [s]"S"(&src[0]),					\
                    [mult]"rm"(x)					\
                : "r8", "r9", "r10", "r11", "rax", "rdx", "memory"

#define ADDSUB_CODE5(op, opc)   \
        "" #op  "q %q[s0], (%[z])\n"					\
        "" #opc "q %q[s1], 0x8(%[z])\n"					\
        "" #opc "q %q[s2], 0x10(%[z])\n"				\
        "" #opc "q %q[s3], 0x18(%[z])\n"				\
        "" #opc "q %q[s4], 0x20(%[z])\n"				\
        "" #opc "q $0x0, 0x28(%[z])\n"					\
                :							\
                :							\
                    [z]"r"(&dst[0]),				        \
                    [s0]"r"(src[0]),					\
                    [s1]"r"(src[1]),					\
                    [s2]"r"(src[2]),					\
                    [s3]"r"(src[3]),					\
                    [s4]"r"(src[4])                                   \
                : "memory"

    EXPOSE_SPECIALIZATION(5);
    /* }}} */

    /* {{{ code for gfp<6, 1> */

#define ADDSUBMUL_CODE6(op, opc)					\
                FEED_IN("0x0(%[s])", "0x8(%[s])", r8, r9)		\
                INNER_MUL(op, "0x0(%[z])", "0x10(%[s])", r8, r9, r10)	\
                INNER_MUL(op, "0x8(%[z])", "0x18(%[s])", r9, r10, r11)	\
                INNER_MUL(op, "0x10(%[z])", "0x20(%[s])", r10, r11, r8)	\
                INNER_MUL(op, "0x18(%[z])", "0x28(%[s])", r11, r8, r9)	\
                FINISH(op, opc,                                         \
                        "0x20(%[z])", "0x28(%[z])", "0x30(%[z])",       \
                        r8, r9)                               	\
                :							\
                :							\
                    [z]"D"(&dst[0]),				        \
                    [s]"S"(&src[0]),					\
                    [mult]"rm"(x)					\
                : "r8", "r9", "r10", "r11", "rax", "rdx", "memory"

#define ADDSUB_CODE6(op, opc)   \
        "" #op  "q %q[s0], (%[z])\n"					\
        "" #opc "q %q[s1], 0x8(%[z])\n"					\
        "" #opc "q %q[s2], 0x10(%[z])\n"				\
        "" #opc "q %q[s3], 0x18(%[z])\n"				\
        "" #opc "q %q[s4], 0x20(%[z])\n"				\
        "" #opc "q %q[s5], 0x28(%[z])\n"				\
        "" #opc "q $0x0, 0x30(%[z])\n"					\
                :							\
                :							\
                    [z]"r"(&dst[0]),				        \
                    [s0]"r"(src[0]),					\
                    [s1]"r"(src[1]),					\
                    [s2]"r"(src[2]),					\
                    [s3]"r"(src[3]),					\
                    [s4]"r"(src[4]),                                  \
                    [s5]"r"(src[5])                                   \
                : "memory"

    EXPOSE_SPECIALIZATION(6);
    /* }}} */

    /* {{{ code for gfp<7, 1> */
#define ADDSUBMUL_CODE7(op, opc)					\
                FEED_IN("0x0(%[s])", "0x8(%[s])", r8, r9)		\
                INNER_MUL(op, "0x0(%[z])", "0x10(%[s])", r8, r9, r10)	\
                INNER_MUL(op, "0x8(%[z])", "0x18(%[s])", r9, r10, r11)	\
                INNER_MUL(op, "0x10(%[z])", "0x20(%[s])", r10, r11, r8)	\
                INNER_MUL(op, "0x18(%[z])", "0x28(%[s])", r11, r8, r9)	\
                INNER_MUL(op, "0x20(%[z])", "0x30(%[s])", r8, r9, r10)	\
                FINISH(op, opc,                                         \
                        "0x28(%[z])", "0x30(%[z])", "0x38(%[z])",       \
                        r9, r10)                               	\
                :							\
                :							\
                    [z]"D"(&dst[0]),				        \
                    [s]"S"(&src[0]),					\
                    [mult]"rm"(x)					\
                : "r8", "r9", "r10", "r11", "rax", "rdx", "memory"

#define ADDSUB_CODE7(op, opc)   \
        "" #op  "q %q[s0], (%[z])\n"					\
        "" #opc "q %q[s1], 0x8(%[z])\n"					\
        "" #opc "q %q[s2], 0x10(%[z])\n"				\
        "" #opc "q %q[s3], 0x18(%[z])\n"				\
        "" #opc "q %q[s4], 0x20(%[z])\n"				\
        "" #opc "q %q[s5], 0x28(%[z])\n"				\
        "" #opc "q %q[s6], 0x30(%[z])\n"				\
        "" #opc "q $0x0, 0x38(%[z])\n"					\
                :							\
                :							\
                    [z]"r"(&dst[0]),				        \
                    [s0]"r"(src[0]),					\
                    [s1]"r"(src[1]),					\
                    [s2]"r"(src[2]),					\
                    [s3]"r"(src[3]),					\
                    [s4]"r"(src[4]),                                  \
                    [s5]"r"(src[5]),                                  \
                    [s6]"r"(src[6])                                   \
                : "memory"

    EXPOSE_SPECIALIZATION(7);
    /* }}} */

    /* {{{ code for gfp<8, 1> */
#define ADDSUBMUL_CODE8(op, opc)					\
                FEED_IN("0x0(%[s])", "0x8(%[s])", r8, r9)		\
                INNER_MUL(op, "0x0(%[z])", "0x10(%[s])", r8, r9, r10)	\
                INNER_MUL(op, "0x8(%[z])", "0x18(%[s])", r9, r10, r11)	\
                INNER_MUL(op, "0x10(%[z])", "0x20(%[s])", r10, r11, r8)	\
                INNER_MUL(op, "0x18(%[z])", "0x28(%[s])", r11, r8, r9)	\
                INNER_MUL(op, "0x20(%[z])", "0x30(%[s])", r8, r9, r10)	\
                INNER_MUL(op, "0x28(%[z])", "0x38(%[s])", r9, r10, r11)	\
                FINISH(op, opc,                                         \
                        "0x30(%[z])", "0x38(%[z])", "0x40(%[z])",       \
                        r10, r11)                               	\
                :							\
                :							\
                    [z]"D"(&dst[0]),				        \
                    [s]"S"(&src[0]),					\
                    [mult]"rm"(x)					\
                : "r8", "r9", "r10", "r11", "rax", "rdx", "memory"

#define ADDSUB_CODE8(op, opc)   \
        "" #op  "q %q[s0], (%[z])\n"					\
        "" #opc "q %q[s1], 0x8(%[z])\n"					\
        "" #opc "q %q[s2], 0x10(%[z])\n"				\
        "" #opc "q %q[s3], 0x18(%[z])\n"				\
        "" #opc "q %q[s4], 0x20(%[z])\n"				\
        "" #opc "q %q[s5], 0x28(%[z])\n"				\
        "" #opc "q %q[s6], 0x30(%[z])\n"				\
        "" #opc "q %q[s7], 0x38(%[z])\n"				\
        "" #opc "q $0x0, 0x40(%[z])\n"					\
                :							\
                :							\
                    [z]"r"(&dst[0]),				        \
                    [s0]"r"(src[0]),					\
                    [s1]"r"(src[1]),					\
                    [s2]"r"(src[2]),					\
                    [s3]"r"(src[3]),					\
                    [s4]"r"(src[4]),                                  \
                    [s5]"r"(src[5]),                                  \
                    [s6]"r"(src[6]),                                   \
                    [s7]"r"(src[7])                                   \
                : "memory"

    EXPOSE_SPECIALIZATION(8);
    /* }}} */
    
    /* further specialization only seem to bring very marginal
     * improvements. This should probably go away. */

#if 0
    /* This is old AVX/SSE carry-save code. It hasn't been tested in a
     * long while, and it has never been put in production as far as I
     * can tell.
     */
    /* AVX/SSE code is a choice which is mostly unrelated to the C++
     * dialect we use, except for a fine point below. We use a union to
     * protect an off-bounds write. And pre-C++11 does not really like
     * unions with contructors for members. Given that anyway we have
     * little use for this code at the moment, since it is suboptimal,
     * let's protect it with HAVE_CXX11
     */
#ifdef HAVE_CXX11
#if defined(HAVE_AVX2) || defined(HAVE_SSSE3)
    template<> struct fast_type<gfp<3, 1> > {
        typedef gfp<3, 1> super;
        struct elt;
        typedef elt elt_ur_for_add;
        struct elt {
#ifdef  HAVE_AVX2
        static const int alignment = 32;
#else
        static const int alignment = 16;
#endif
            typedef elt self;
#ifdef  HAVE_AVX2
            __m256i data[1];
#else
            __m128i data[2];
#endif
            elt() { zero(); }
            elt(elt const& a) = default;
            elt& operator=(elt const& a) = default;

            /* we do not construct (nor affect) from mpz, because we're not
             * positional */
            void zero() {
#ifdef  HAVE_AVX2
                data[0] = _mm256_setzero_si256();
#else
                data[0] = _mm_setzero_si128();
                data[1] = _mm_setzero_si128();
#endif
            }
            static void zero(elt * x, int N) {
                // see bug 21663
                memset(x->data, 0, N * sizeof(data));
            }
            static void copy(elt * y, const elt * x, int N) {
                // see bug 21663
                memcpy(y->data, x->data, N * sizeof(data));
            }
            bool operator==(elt const& a) {
                return memcmp(data, a.data, sizeof(data)) == 0;
            }
            elt(super::elt const& a) {
                convert(*this, a);
            }

            operator super::elt_ur_for_add() const {
                super::elt_ur_for_add carries(conv_backend_get_carries(*this));
                super::add(carries, conv_backend_get_main(*this));
                return carries;
            }

            /* same, but we assume carry is zero */
            operator super::elt() const {
                ASSERT(conv_backend_get_carries(*this).is_zero());
                return conv_backend_get_main(*this);
            }
        };

        static inline void stream_store(elt * dst, elt const& src) {
            /* Do we want to stream that or not ? In fact it's slower
             * when streaming... */
#if 0
#ifdef  HAVE_AVX2
            _mm256_stream_si256(dst->data+0, src.data[0]);
#else
            _mm_stream_si128(dst->data+0, src.data[0]);
            _mm_stream_si128(dst->data+1, src.data[1]);
#endif
#else
#ifdef  HAVE_AVX2
            _mm256_storeu_si256(dst->data+0, src.data[0]);
#else
            _mm_storeu_si128(dst->data+0, src.data[0]);
            _mm_storeu_si128(dst->data+1, src.data[1]);
#endif
#endif
        }
        static inline void add(elt & dst, elt const & src)
        {
#ifdef  HAVE_AVX2
            dst.data[0] = _mm256_add_epi64 (dst.data[0], src.data[0]);
#else
            dst.data[0] = _mm_add_epi64 (dst.data[0], src.data[0]);
            dst.data[1] = _mm_add_epi64 (dst.data[1], src.data[1]);
#endif
        }

        static inline void sub(elt & dst, elt const & src)
        {
#ifdef  HAVE_AVX2
            dst.data[0] = _mm256_sub_epi64 (dst.data[0], src.data[0]);
#else
            dst.data[0] = _mm_sub_epi64 (dst.data[0], src.data[0]);
            dst.data[1] = _mm_sub_epi64 (dst.data[1], src.data[1]);
#endif
        }

        static inline void sub_ur(elt_ur_for_add & dst, elt_ur_for_add const & src)
        {
#ifdef  HAVE_AVX2
            dst.data[0] = _mm256_sub_epi64 (dst.data[0], src.data[0]);
#else
            dst.data[0] = _mm_sub_epi64 (dst.data[0], src.data[0]);
            dst.data[1] = _mm_sub_epi64 (dst.data[1], src.data[1]);
#endif
        }

        /* conversions are done as a combination of blend & shuffle */

#ifdef  HAVE_AVX2
        /* We grok only values for w_i which are integer immediates
         * within {-1} \cup {0..15}
         */
#define shuffle_16bit_words(out, in,        		        	\
            w0, w1, w2, w3,						\
            w4, w5, w6, w7,						\
            w8, w9, wa, wb,						\
            wc, wd, we, wf)						\
    do {                                                                \
        *out = _mm256_xor_si256(				        \
            _mm256_shuffle_epi8(					\
            in,								\
            _mm256_setr_epi8( 						\
                (w0 < 0) ? -1 : ((w0 < 8)  ? 2*(w0&7) + 0 : -1),	\
                (w0 < 0) ? -1 : ((w0 < 8)  ? 2*(w0&7) + 1 : -1),	\
                (w1 < 0) ? -1 : ((w1 < 8)  ? 2*(w1&7) + 0 : -1),	\
                (w1 < 0) ? -1 : ((w1 < 8)  ? 2*(w1&7) + 1 : -1),	\
                (w2 < 0) ? -1 : ((w2 < 8)  ? 2*(w2&7) + 0 : -1),	\
                (w2 < 0) ? -1 : ((w2 < 8)  ? 2*(w2&7) + 1 : -1),	\
                (w3 < 0) ? -1 : ((w3 < 8)  ? 2*(w3&7) + 0 : -1),	\
                (w3 < 0) ? -1 : ((w3 < 8)  ? 2*(w3&7) + 1 : -1),	\
                (w4 < 0) ? -1 : ((w4 < 8)  ? 2*(w4&7) + 0 : -1),	\
                (w4 < 0) ? -1 : ((w4 < 8)  ? 2*(w4&7) + 1 : -1),	\
                (w5 < 0) ? -1 : ((w5 < 8)  ? 2*(w5&7) + 0 : -1),	\
                (w5 < 0) ? -1 : ((w5 < 8)  ? 2*(w5&7) + 1 : -1),	\
                (w6 < 0) ? -1 : ((w6 < 8)  ? 2*(w6&7) + 0 : -1),	\
                (w6 < 0) ? -1 : ((w6 < 8)  ? 2*(w6&7) + 1 : -1),	\
                (w7 < 0) ? -1 : ((w7 < 8)  ? 2*(w7&7) + 0 : -1),	\
                (w7 < 0) ? -1 : ((w7 < 8)  ? 2*(w7&7) + 1 : -1),	\
                (w8 < 0) ? -1 : ((w8 >= 8) ? 2*(w8&7) + 0 : -1),	\
                (w8 < 0) ? -1 : ((w8 >= 8) ? 2*(w8&7) + 1 : -1),	\
                (w9 < 0) ? -1 : ((w9 >= 8) ? 2*(w9&7) + 0 : -1),	\
                (w9 < 0) ? -1 : ((w9 >= 8) ? 2*(w9&7) + 1 : -1),	\
                (wa < 0) ? -1 : ((wa >= 8) ? 2*(wa&7) + 0 : -1),	\
                (wa < 0) ? -1 : ((wa >= 8) ? 2*(wa&7) + 1 : -1),	\
                (wb < 0) ? -1 : ((wb >= 8) ? 2*(wb&7) + 0 : -1),	\
                (wb < 0) ? -1 : ((wb >= 8) ? 2*(wb&7) + 1 : -1),	\
                (wc < 0) ? -1 : ((wc >= 8) ? 2*(wc&7) + 0 : -1),	\
                (wc < 0) ? -1 : ((wc >= 8) ? 2*(wc&7) + 1 : -1),	\
                (wd < 0) ? -1 : ((wd >= 8) ? 2*(wd&7) + 0 : -1),	\
                (wd < 0) ? -1 : ((wd >= 8) ? 2*(wd&7) + 1 : -1),	\
                (we < 0) ? -1 : ((we >= 8) ? 2*(we&7) + 0 : -1),	\
                (we < 0) ? -1 : ((we >= 8) ? 2*(we&7) + 1 : -1),	\
                (wf < 0) ? -1 : ((wf >= 8) ? 2*(wf&7) + 0 : -1),	\
                (wf < 0) ? -1 : ((wf >= 8) ? 2*(wf&7) + 1 : -1))),	\
            _mm256_shuffle_epi8(					\
                /* 0x4E is 0b01001110 aka (1,0,3,2) */			\
                _mm256_permute4x64_epi64 (in, _MM_SHUFFLE(1,0,3,2)), 	\
            _mm256_setr_epi8( 						\
                (w0 < 0) ? -1 : ((w0 >= 8) ? 2*(w0&7) + 0 : -1),	\
                (w0 < 0) ? -1 : ((w0 >= 8) ? 2*(w0&7) + 1 : -1),	\
                (w1 < 0) ? -1 : ((w1 >= 8) ? 2*(w1&7) + 0 : -1),	\
                (w1 < 0) ? -1 : ((w1 >= 8) ? 2*(w1&7) + 1 : -1),	\
                (w2 < 0) ? -1 : ((w2 >= 8) ? 2*(w2&7) + 0 : -1),	\
                (w2 < 0) ? -1 : ((w2 >= 8) ? 2*(w2&7) + 1 : -1),	\
                (w3 < 0) ? -1 : ((w3 >= 8) ? 2*(w3&7) + 0 : -1),	\
                (w3 < 0) ? -1 : ((w3 >= 8) ? 2*(w3&7) + 1 : -1),	\
                (w4 < 0) ? -1 : ((w4 >= 8) ? 2*(w4&7) + 0 : -1),	\
                (w4 < 0) ? -1 : ((w4 >= 8) ? 2*(w4&7) + 1 : -1),	\
                (w5 < 0) ? -1 : ((w5 >= 8) ? 2*(w5&7) + 0 : -1),	\
                (w5 < 0) ? -1 : ((w5 >= 8) ? 2*(w5&7) + 1 : -1),	\
                (w6 < 0) ? -1 : ((w6 >= 8) ? 2*(w6&7) + 0 : -1),	\
                (w6 < 0) ? -1 : ((w6 >= 8) ? 2*(w6&7) + 1 : -1),	\
                (w7 < 0) ? -1 : ((w7 >= 8) ? 2*(w7&7) + 0 : -1),	\
                (w7 < 0) ? -1 : ((w7 >= 8) ? 2*(w7&7) + 1 : -1),	\
                (w8 < 0) ? -1 : ((w8 < 8)  ? 2*(w8&7) + 0 : -1),	\
                (w8 < 0) ? -1 : ((w8 < 8)  ? 2*(w8&7) + 1 : -1),	\
                (w9 < 0) ? -1 : ((w9 < 8)  ? 2*(w9&7) + 0 : -1),	\
                (w9 < 0) ? -1 : ((w9 < 8)  ? 2*(w9&7) + 1 : -1),	\
                (wa < 0) ? -1 : ((wa < 8)  ? 2*(wa&7) + 0 : -1),	\
                (wa < 0) ? -1 : ((wa < 8)  ? 2*(wa&7) + 1 : -1),	\
                (wb < 0) ? -1 : ((wb < 8)  ? 2*(wb&7) + 0 : -1),	\
                (wb < 0) ? -1 : ((wb < 8)  ? 2*(wb&7) + 1 : -1),	\
                (wc < 0) ? -1 : ((wc < 8)  ? 2*(wc&7) + 0 : -1),	\
                (wc < 0) ? -1 : ((wc < 8)  ? 2*(wc&7) + 1 : -1),	\
                (wd < 0) ? -1 : ((wd < 8)  ? 2*(wd&7) + 0 : -1),	\
                (wd < 0) ? -1 : ((wd < 8)  ? 2*(wd&7) + 1 : -1),	\
                (we < 0) ? -1 : ((we < 8)  ? 2*(we&7) + 0 : -1),	\
                (we < 0) ? -1 : ((we < 8)  ? 2*(we&7) + 1 : -1),	\
                (wf < 0) ? -1 : ((wf < 8)  ? 2*(wf&7) + 0 : -1),	\
                (wf < 0) ? -1 : ((wf < 8)  ? 2*(wf&7) + 1 : -1)))	\
        );                                                              \
    } while (0)
#else
#define shuffle_16bit_words(out, lo, hi,       		        	\
            w0, w1, w2, w3,						\
            w4, w5, w6, w7,						\
            w8, w9, wa, wb,						\
            wc, wd, we, wf)						\
    do {                                                                \
        out[0] = _mm_xor_si128(			        	        \
            _mm_shuffle_epi8(lo, _mm_setr_epi8( 			\
                (w0 < 0) ? -1 : ((w0 < 8)  ? 2*(w0&7) + 0 : -1),	\
                (w0 < 0) ? -1 : ((w0 < 8)  ? 2*(w0&7) + 1 : -1),	\
                (w1 < 0) ? -1 : ((w1 < 8)  ? 2*(w1&7) + 0 : -1),	\
                (w1 < 0) ? -1 : ((w1 < 8)  ? 2*(w1&7) + 1 : -1),	\
                (w2 < 0) ? -1 : ((w2 < 8)  ? 2*(w2&7) + 0 : -1),	\
                (w2 < 0) ? -1 : ((w2 < 8)  ? 2*(w2&7) + 1 : -1),	\
                (w3 < 0) ? -1 : ((w3 < 8)  ? 2*(w3&7) + 0 : -1),	\
                (w3 < 0) ? -1 : ((w3 < 8)  ? 2*(w3&7) + 1 : -1),	\
                (w4 < 0) ? -1 : ((w4 < 8)  ? 2*(w4&7) + 0 : -1),	\
                (w4 < 0) ? -1 : ((w4 < 8)  ? 2*(w4&7) + 1 : -1),	\
                (w5 < 0) ? -1 : ((w5 < 8)  ? 2*(w5&7) + 0 : -1),	\
                (w5 < 0) ? -1 : ((w5 < 8)  ? 2*(w5&7) + 1 : -1),	\
                (w6 < 0) ? -1 : ((w6 < 8)  ? 2*(w6&7) + 0 : -1),	\
                (w6 < 0) ? -1 : ((w6 < 8)  ? 2*(w6&7) + 1 : -1),	\
                (w7 < 0) ? -1 : ((w7 < 8)  ? 2*(w7&7) + 0 : -1),	\
                (w7 < 0) ? -1 : ((w7 < 8)  ? 2*(w7&7) + 1 : -1))),	\
            _mm_shuffle_epi8(hi, _mm_setr_epi8( 			\
                (w0 < 0) ? -1 : ((w0 >= 8) ? 2*(w0&7) + 0 : -1),	\
                (w0 < 0) ? -1 : ((w0 >= 8) ? 2*(w0&7) + 1 : -1),	\
                (w1 < 0) ? -1 : ((w1 >= 8) ? 2*(w1&7) + 0 : -1),	\
                (w1 < 0) ? -1 : ((w1 >= 8) ? 2*(w1&7) + 1 : -1),	\
                (w2 < 0) ? -1 : ((w2 >= 8) ? 2*(w2&7) + 0 : -1),	\
                (w2 < 0) ? -1 : ((w2 >= 8) ? 2*(w2&7) + 1 : -1),	\
                (w3 < 0) ? -1 : ((w3 >= 8) ? 2*(w3&7) + 0 : -1),	\
                (w3 < 0) ? -1 : ((w3 >= 8) ? 2*(w3&7) + 1 : -1),	\
                (w4 < 0) ? -1 : ((w4 >= 8) ? 2*(w4&7) + 0 : -1),	\
                (w4 < 0) ? -1 : ((w4 >= 8) ? 2*(w4&7) + 1 : -1),	\
                (w5 < 0) ? -1 : ((w5 >= 8) ? 2*(w5&7) + 0 : -1),	\
                (w5 < 0) ? -1 : ((w5 >= 8) ? 2*(w5&7) + 1 : -1),	\
                (w6 < 0) ? -1 : ((w6 >= 8) ? 2*(w6&7) + 0 : -1),	\
                (w6 < 0) ? -1 : ((w6 >= 8) ? 2*(w6&7) + 1 : -1),	\
                (w7 < 0) ? -1 : ((w7 >= 8) ? 2*(w7&7) + 0 : -1),	\
                (w7 < 0) ? -1 : ((w7 >= 8) ? 2*(w7&7) + 1 : -1))));	\
        out[1] = _mm_xor_si128(			        	        \
            _mm_shuffle_epi8(lo, _mm_setr_epi8( 			\
                (w8 < 0) ? -1 : ((w8 < 8)  ? 2*(w8&7) + 0 : -1),	\
                (w8 < 0) ? -1 : ((w8 < 8)  ? 2*(w8&7) + 1 : -1),	\
                (w9 < 0) ? -1 : ((w9 < 8)  ? 2*(w9&7) + 0 : -1),	\
                (w9 < 0) ? -1 : ((w9 < 8)  ? 2*(w9&7) + 1 : -1),	\
                (wa < 0) ? -1 : ((wa < 8)  ? 2*(wa&7) + 0 : -1),	\
                (wa < 0) ? -1 : ((wa < 8)  ? 2*(wa&7) + 1 : -1),	\
                (wb < 0) ? -1 : ((wb < 8)  ? 2*(wb&7) + 0 : -1),	\
                (wb < 0) ? -1 : ((wb < 8)  ? 2*(wb&7) + 1 : -1),	\
                (wc < 0) ? -1 : ((wc < 8)  ? 2*(wc&7) + 0 : -1),	\
                (wc < 0) ? -1 : ((wc < 8)  ? 2*(wc&7) + 1 : -1),	\
                (wd < 0) ? -1 : ((wd < 8)  ? 2*(wd&7) + 0 : -1),	\
                (wd < 0) ? -1 : ((wd < 8)  ? 2*(wd&7) + 1 : -1),	\
                (we < 0) ? -1 : ((we < 8)  ? 2*(we&7) + 0 : -1),	\
                (we < 0) ? -1 : ((we < 8)  ? 2*(we&7) + 1 : -1),	\
                (wf < 0) ? -1 : ((wf < 8)  ? 2*(wf&7) + 0 : -1),	\
                (wf < 0) ? -1 : ((wf < 8)  ? 2*(wf&7) + 1 : -1))),	\
            _mm_shuffle_epi8(hi, _mm_setr_epi8(                         \
                (w8 < 0) ? -1 : ((w8 >= 8) ? 2*(w8&7) + 0 : -1),	\
                (w8 < 0) ? -1 : ((w8 >= 8) ? 2*(w8&7) + 1 : -1),	\
                (w9 < 0) ? -1 : ((w9 >= 8) ? 2*(w9&7) + 0 : -1),	\
                (w9 < 0) ? -1 : ((w9 >= 8) ? 2*(w9&7) + 1 : -1),	\
                (wa < 0) ? -1 : ((wa >= 8) ? 2*(wa&7) + 0 : -1),	\
                (wa < 0) ? -1 : ((wa >= 8) ? 2*(wa&7) + 1 : -1),	\
                (wb < 0) ? -1 : ((wb >= 8) ? 2*(wb&7) + 0 : -1),	\
                (wb < 0) ? -1 : ((wb >= 8) ? 2*(wb&7) + 1 : -1),	\
                (wc < 0) ? -1 : ((wc >= 8) ? 2*(wc&7) + 0 : -1),	\
                (wc < 0) ? -1 : ((wc >= 8) ? 2*(wc&7) + 1 : -1),	\
                (wd < 0) ? -1 : ((wd >= 8) ? 2*(wd&7) + 0 : -1),	\
                (wd < 0) ? -1 : ((wd >= 8) ? 2*(wd&7) + 1 : -1),	\
                (we < 0) ? -1 : ((we >= 8) ? 2*(we&7) + 0 : -1),	\
                (we < 0) ? -1 : ((we >= 8) ? 2*(we&7) + 1 : -1),	\
                (wf < 0) ? -1 : ((wf >= 8) ? 2*(wf&7) + 0 : -1),	\
                (wf < 0) ? -1 : ((wf >= 8) ? 2*(wf&7) + 1 : -1))));	\
    } while (0)
#endif

        /* case of 192 bits within 256 bits. Three 64-bit words
         * split into four 48-bit words.
         */
        static void convert(elt& dst, const super::elt& a) {
            /* index of 16-bit word in destination, fetched from
             * which index of 16-bit word in the gfp::elt. This is
             * given for the 256-bit registers
             *
             * 0    0
             * 1    1
             * 2    2
             * 3    <empty>
             * 4    3
             * 5    4
             * 6    5
             * 7    <empty>
             * 8    6
             * 9    7
             * 10   8
             * 11   <empty>
             * 12   9
             * 13   10
             * 14   11
             * 15   <empty>
             */
#ifdef  HAVE_AVX2
            /* I'm really upset here. _mm256_shuffle_epi8, aka VPSHUFB,
             * reads only 4-byte immediates (and discards the rest). As a
             * consequence, the following does not work: the indices
             * 12,13,14,15 read off bounds, while the 16,17, etc actually
             * do what they want, but based on the fact that they're
             * reduced mod 16 + implicitly considered wrt the high part
             * of the operand...
            dst.data[0] = _mm256_shuffle_epi8(
                    _mm256_loadu_si256((__m256i*) a.x),
                    _mm256_setr_epi8( 
                        0,1,2,3,4,5,-1,-1,
                        6,7,8,9,10,11,-1,-1,
                        12,13,14,15,16,17,-1,-1,
                        18,19,20,21,22,23,-1,-1));
            */

#if 0
            __m256i in = _mm256_loadu_si256((__m256i*) a.x);
            dst.data[0] =
                    _mm256_xor_si256(
                        _mm256_shuffle_epi8(
                        in,
                        _mm256_setr_epi8( 
                            0,1,2,3,4,5,-1,-1,
                            6,7,8,9,10,11,-1,-1,
                            -1,-1,-1,-1,0,1,-1,-1,
                            2,3,4,5,6,7,-1,-1)),
                        _mm256_shuffle_epi8(
                            /* 0x4E is 0b01001110 aka (1,0,3,2) */
                            _mm256_permute4x64_epi64 (in, _MM_SHUFFLE(1,0,3,2)), 
                        _mm256_setr_epi8( 
                            -1,-1,-1,-1,-1,-1,-1,-1,
                            -1,-1,-1,-1,-1,-1,-1,-1,
                            12,13,14,15,-1,-1,-1,-1,
                            -1,-1,-1,-1,-1,-1,-1,-1)));
#endif
            __m256i in = _mm256_loadu_si256((__m256i*) a.x);
            shuffle_16bit_words(dst.data, in,
                    0,1,2,-1, 3,4,5,-1, 6,7,8,-1, 9,10,11,-1);
#else   /* SSSE3 !! */
            __m128i lo = _mm_loadu_si128((__m128i*) a.x);
            __m128i hi = _mm_loadu_si128((__m128i*) (a.x + 2));
            shuffle_16bit_words(dst.data, lo, hi,
                    0,1,2,-1, 3,4,5,-1, 6,7,8,-1, 9,10,11,-1);
            /* note that 16bit-wide shuffles use an 8-bit immediate,
             * but do not offer the option to selectively insert
             * zeroes. So we're probably better off shuffling bytes.
             */
#endif
        }

        static super::elt conv_backend_get_main(elt const& src) {
            /* This is needed because we knowingly write off bounds */
            union station {
                super::elt e;
#ifdef  HAVE_AVX2
                __m256i v[1];
#else
                __m128i v[2];
#endif
                station() {};
            } main;
#ifdef  HAVE_AVX2
            shuffle_16bit_words(main.v, src.data[0],
                    0,1,2, 4,5,6, 8,9,10, 12,13,14, -1,-1,-1,-1);
#else
            shuffle_16bit_words(main.v, src.data[0], src.data[1],
                    0,1,2, 4,5,6, 8,9,10, 12,13,14, -1,-1,-1,-1);
#endif
            return main.e;
        }
        static super::elt_ur_for_add conv_backend_get_carries(elt const& src) {
            union station {
                super::elt_ur_for_add e;
#ifdef  HAVE_AVX2
                __m256i v[1];
#else
                __m128i v[2];
#endif
                station() {};
            } carries, ncarries;

            /* It's slightly more complicated than it seems. The carry
             * words may be negative. So we must sign-extend them to the
             * full unreduced element size.
             */
#ifdef  HAVE_AVX2
            shuffle_16bit_words(carries.v, src.data[0],
                    -1,-1,-1,3,
                    -1,-1,7,-1,
                    -1,11,-1,-1,
                    15,-1,-1,-1);
            __m256i zero = _mm256_setzero_si256();
            shuffle_16bit_words(ncarries.v,
                    _mm256_sub_epi16(zero,
                        _mm256_cmpgt_epi16(zero, carries.v[0])),
                    -1,-1,-1,-1,
                    3,-1,-1,6,
                    -1,-1,9,-1,
                    -1,12,-1,-1);
#else
            shuffle_16bit_words(carries.v, src.data[0], src.data[1],
                    -1,-1,-1,3,
                    -1,-1,7,-1,
                    -1,11,-1,-1,
                    15,-1,-1,-1);
            __m128i zero = _mm_setzero_si128();
            shuffle_16bit_words(ncarries.v,
                    _mm_sub_epi16(zero, _mm_cmpgt_epi16(zero, carries.v[0])),
                    _mm_sub_epi16(zero, _mm_cmpgt_epi16(zero, carries.v[1])),
                    -1,-1,-1,-1,
                    3,-1,-1,6,
                    -1,-1,9,-1,
                    -1,12,-1,-1);
#endif
            super::sub_ur(carries.e, ncarries.e);
            return carries.e;
        }



        /* (add|sub)mul_ui go through convert, do naively and convert
         * back. Yes, it's slightly painful. Here we assume that src
         * has undergone little to no accumulated additions, so that
         * it can basically be converted lossless to a gfp::elt
         *
         * This prototype was changed recently in the "main"
         * implementation. Ideally, the way to go would be to make the
         * version below a member function rather than a static function
         * that receives the prime and preinv.
         */
        static inline void addmul_ui(elt & dst, elt const & src, mp_limb_t x, super::elt const & p, super::preinv const & j)
        {
            super::elt zr;
            super::elt_ur_for_add z(dst);
            super::addmul_ui(z, (super::elt) src, x, p, j);
            super::reduce(zr, z, p, j);
            dst = zr;
        }
        static inline void submul_ui(elt_ur_for_add & dst, elt const & src, mp_limb_t x, super::elt const & p, super::preinv const & j)
        {
            super::elt zr;
            super::elt_ur_for_add z(dst);
            super::submul_ui(z, (super::elt) src, x, p, j);
            super::reduce(zr, z, p, j);
            dst = zr;
        }

        /* we have *TWO* reduction functions here. One which assigns to a
         * standard gfp::elt, and one which assigns to fast_type::elt */
        static void reduce(super::elt & r, elt const & a, super::elt const & p, super::preinv const & j)
        {
            super::elt_ur_for_add z(a);
            super::reduce(r, z, p, j);
        }
        static void reduce(elt & r, elt const & a, super::elt const & p, super::preinv const & j)
        {
            super::elt zr;
            reduce(zr, a, p, j);
            r = zr;
        }
    };
#endif  /* defined(HAVE_AVX2) || defined(HAVE_SSSE3) */
#endif
#endif

#endif
    }

/* expose only what we have in our public interface */
using details::gfp;
// using details::fast_type;
}



#endif
