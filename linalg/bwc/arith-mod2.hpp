#ifndef ARITH_MOD2_HPP_
#define ARITH_MOD2_HPP_

#include <array>
#include <ostream>
#include <cstdint>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <string>
#if defined(HAVE_SSE41) && defined(HAVE_POPCNT)
#include <x86intrin.h>
#endif


#include <gmp.h>
#include "fmt/format.h"
#include "gmp_aux.h"
#include "misc.h"       // u64_random
#include "cxx_mpz.hpp"

#include "arith-concrete-base.hpp"

namespace arith_mod2 {
    namespace details {

        /* The methods in flat_elt are added to the class only if the
         * underlying type is flat (all bits are data, no pointer)
         */

        template<typename E, bool check = E::is_flat_storage>
        struct flat_elt {};

        template<typename E>
        struct flat_elt<E, true> {
                /*{{{ element layout information */
                static constexpr inline size_t elt_stride() { return sizeof(E); }
                static inline size_t vec_elt_stride(size_t s) { return s * elt_stride(); }
                /*}}}*/

                /*{{{ allocation / deallocation of (vectors of) elements */
                /* These allocation interfaces seem a bit stupid. At the low
                 * hard level, we know that they're just the same as new[]
                 * anyway.
                 */

                static inline E * alloc(size_t k = 1) { return new E[k]; }
                static inline void free(E * u) { delete[] u; }
                static inline E * realloc(E * u, size_t k0, size_t k) {
                    E * v = new E[k];
                    std::copy_n(u, k0, v);
                    delete[] u;
                    return v;
                }
                /*}}}*/

        };

        /* The base type. This can be overridden later on. */
        template<unsigned int G>
            struct gf2_base_type : public std::array<uint64_t, G/64>, public arith_concrete_base::elt {
                    static constexpr const bool is_flat_storage = true;
#if 0
                    /* nlimbs is only for flat types which have a
                     * subscript operator.
                     */
                    static_assert(sizeof(decltype(((elt*)nullptr)->operator[](0))) * elt::nlimbs == sizeof(elt), "nlimbs is wrong");
                    static constexpr const bool nlimbs = K;
#endif
                    static constexpr const unsigned int K = G/64;
                    typedef std::array<uint64_t, K> super;
                    template<typename T> struct compatible_with {
                        static constexpr const bool value = is_flat_storage && (sizeof(super) == sizeof(T));
                    };
                    template<typename... Args> gf2_base_type(Args&&... args) : super(std::forward<Args>(args)...) {}
                    inline void zero() { std::fill_n(super::begin(), K, 0); }
                };

        template<unsigned int G, typename T>
            struct gf2_base : public arith_concrete_base, public flat_elt<gf2_base_type<G>> {
                static constexpr const bool is_characteristic_two = true;
                static constexpr const unsigned int simd_groupsize = G;
                static constexpr const unsigned int K = G / 64;

                typedef gf2_base_type<G> elt;
                static inline void
                    simd_set_ui_at(elt & p, size_t k, int v)
                    {
                        uint64_t * xp = p.begin();
                        uint64_t mask = ((uint64_t)1) << (k%64);
                        xp[k/64] = (xp[k/64] & ~mask) | ((((uint64_t)v) << (k%64))&mask);
                    }
                static inline int
                    simd_test_ui_at(elt const & p, size_t k)
                    {
                        uint64_t const * xp = p.begin();
                        uint64_t mask = ((uint64_t)(1)) << (k%64);
                        return (xp[k/64] & mask) != 0;
                    }
                static inline void
                    simd_add_ui_at(elt & p, size_t k, int v)
                    {
                        uint64_t * xp = p.begin();
                        uint64_t mask = ((uint64_t)(v&1)) << (k%64);
                        xp[k/64] ^= mask;
                    }
                static inline int
                    simd_hamming_weight(elt const & p)
                    {
                        int w = 0;
#if GNUC_VERSION_ATLEAST(3,4,0)
                        const unsigned long * xp = reinterpret_cast<const unsigned long *>(&p);
                        static_assert(T::elt::is_flat_storage && ((sizeof(typename T::elt) % sizeof(unsigned long)) == 0), "bad types");
                        // parentheses in the divisor are mandatory here, see
                        // https://gcc.gnu.org/pipermail/gcc-patches/2020-September/553888.html
                        for(size_t b = 0 ; b < T::elt_stride() / (sizeof(unsigned long)) ; b++) {
                            w += __builtin_popcountl(xp[b]);
                        }
#else
                        int tab[16] = { 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4 };
                        uint8_t const * xp = reinterpret_cast<uint8_t const *> p;
                        for(size_t b = 0 ; b < T::elt_stride() ; b++) {
                            w += tab[xp[b]&15] + tab[xp[b]>>4];
                        }
#endif
                        return w;
                    }

#if 0
                /*{{{ element layout information */
                static inline size_t elt_stride() { return sizeof(elt); }
                static inline size_t vec_elt_stride(size_t s) { return s * elt_stride(); }
                /*}}}*/

#endif

                /*{{{ predicates */
                static inline bool
                    is_zero(elt const & x) { 
                        for(unsigned int i = 0 ; i < K ; i++) {
                            if (x[i]) return false;
                        }
                        return true;
                    }

                static inline int
                    cmp(elt const & x, unsigned long a) {
                        int r = (x[0] > a) - (a > x[0]);
                        if (r) return r;
                        for(unsigned int i = 1 ; i < K ; i++)
                            if (x[i]) return 1;
                        return 0;
                    }

                static inline int
                    cmp(elt const & x, elt const & y) {
                        for(unsigned int i = 0 ; i < K ; i++) {
                            int r = (x[i] > y[i]) - (y[i] > x[i]);
                            if (r) return r;
                        }
                        return 0;
                    }
                /*}}}*/

                /* {{{ assignments */
                static inline elt &
                    set(elt& x, unsigned long a) { x[0] = a; return x; }

                static inline elt &
                    set(elt& x, cxx_mpz const & a) { x[0] = mpz_get_ui(a); return x; }


                static inline elt & set(elt& x, elt const & a) { return x = a; }
                static inline elt & neg(elt& x, elt const & a) { return set(x, a); }

                static inline void
                    set_zero(elt & x) { x.zero(); }

                static inline void 
                    set_random(elt& x, gmp_randstate_ptr rstate) {
                        for(unsigned int i = 0 ; i < K ; i++)
                            x[i] = u64_random(rstate);
                    }

                static inline void stream_store(elt * dst, elt const& src) { *dst = src; }

                static inline void reduce(elt&) {}
                /* }}} */

                /*{{{ addition, at the element level */
                static inline void
                    add(elt & dst, elt const & a, elt const & b)
                    {
                        for(unsigned int i = 0 ; i < K ; i++)
                            dst[i] = a[i] ^ b[i];
                    }
                /* It seems that we do not have to implement:
                 *  addmul
                 *  addmul_ur
                 *  addmul_and_reduce
                 *  addmul_ui
                 *  sub
                 *
                 */
                /*}}}*/

                /*{{{ I/O */
                static inline std::ostream& cxx_out(std::ostream& o, elt const & x) {
                    for(unsigned int i = 0 ; i < K ; i++) {
                        // We _know_ that we're dealing with a 64-bit type
                        // anyway, so what follows is unnecessary.
                        // o << fmt::format(FMT_STRING("{0:0{1}}"), x[i], sizeof(x[i]) * CHAR_BIT / 4);
                        o << fmt::format(FMT_STRING("{0:016x}"), x[i]);
                    }
                    return o;
                }
                static inline std::ostream& write(std::ostream& o, elt const & x) {
                    return o.write((const char *) x.x, sizeof(elt));
                }
                inline int fread(FILE * f, elt & x) const {
                    int ret = ::fread((char *) x.x, sizeof(elt), 1, f);
                    if (ret < 1) return ret;
                    return ret;
                }
                inline int fscan(FILE * f, elt& x) const
                {
                    char buf[K * 16 + 1];
                    int ret = fscanf(f, "%*s", sizeof(buf), buf);
                    if (!ret) return 0;
                    char * b = buf;
                    for(unsigned int i = 0 ; i < K ; i++, b+=16) {
                        std::string s(b, b+16);
                        std::istringstream is(s);
                        is >> std::hex >> x[i];
                    }
                    return ret;
                }
                /*}}}*/
            };

#if defined(HAVE_SSE41) && defined(HAVE_POPCNT)
        template<>
            struct gf2_base_type<128> : public arith_concrete_base::elt {
                static constexpr const bool is_flat_storage = true;
                __m128i x;
                void zero() { x = _mm_setzero_si128(); }
            };

        template<typename T>
            struct gf2_base<128, T>: public arith_concrete_base, public flat_elt<gf2_base_type<128>> {
                static constexpr const unsigned int G = 128;
                static constexpr const bool is_characteristic_two = true;
                static constexpr const unsigned int simd_groupsize = G;
                typedef gf2_base_type<G> elt;
                /* {{{ assignments */
                static inline elt & set(elt& x, elt const & a) { return x = a; }

                static inline void set_zero(elt & x) { x.zero(); }

                static inline void 
                    set_random(elt& x, gmp_randstate_ptr rstate) {
                        __m64 lo =  _mm_cvtsi64_m64(u64_random(rstate));
                        __m64 hi =  _mm_cvtsi64_m64(u64_random(rstate));
                        x = _mm_setr_epi64(lo, hi);
                    }

                static inline void stream_store(elt * dst, elt const& src) { *dst = src; }

                static inline void reduce(elt&) {}
                /* }}} */
                static inline void
                    add(elt & dst, elt const & a, elt const & b)
                    {
                        dst.x = _mm_xor_si128(a.x, b.x);
                    }
                /*{{{ simd*/
                static inline void
                    simd_set_ui_at(elt & p, size_t k, int v)
                    {
                        uint64_t * xp = reinterpret_cast<uint64_t *>(&p);
                        uint64_t mask = ((uint64_t)1) << (k%64);
                        xp[k/64] = (xp[k/64] & ~mask) | ((((uint64_t)v) << (k%64))&mask);
                    }
                static inline void
                    simd_add_ui_at(elt & p, size_t k, int v)
                    {
                        uint64_t * xp = reinterpret_cast<uint64_t *>(&p);
                        uint64_t mask = ((uint64_t)(v&1)) << (k%64);
                        xp[k/64] ^= mask;
                    }
                static inline int
                    simd_hamming_weight(elt const & p)
                    {
                        return _mm_popcnt_u64(_mm_extract_epi64(p.x, 0)) +
                            _mm_popcnt_u64(_mm_extract_epi64(p.x, 1));
                    }
                /*}}}*/
            };

#endif

        /*************************************************************/
        /* code below this point does not dive into the details of
         * the elt type.
         */

        template<unsigned int G, typename T>
            struct gf2_middle : public gf2_base<G, T> {
                typedef gf2_base<G, T> super;
                using typename super::elt;

                mpz_srcptr characteristic() const {
                    static mp_limb_t limbs[1] = {2};
                    static __mpz_struct a = { 1, 1, limbs };
                    return &a;
                }

                /*{{{ a few trivial proxies */
                using super::add;
                static inline void
                    add(elt & dst, elt const & b)
                    {
                        T::add(dst, dst, b);
                    }
                static inline void
                    add_and_reduce(elt & dst, elt const & b)
                    {
                        T::add(dst, dst, b);
                    }
                static inline void
                    sub_and_reduce(elt & dst, elt const & b)
                    {
                        T::add(dst, dst, b);
                    }
                /*}}}*/

            /* TODO: mul / addmul interface. What do we have to do?
             */

                /*{{{ accessors inside vectors */
                /* Likewise, this seems a bit stupid. THe only reason to have
                 * them is if we want to have a pz layer. But it's very
                 * doubtful that we can manage to have that with flat layout.
                 */
                static inline elt *
                    vec_subvec(elt * p, size_t k) { return p + k; }

                static inline elt &
                    vec_item(elt * p, size_t k) { return p[k]; }

                static inline elt const *
                    vec_subvec(elt const * p, size_t k) { return p + k; }

                static inline elt const &
                    vec_item(elt const * p, size_t k) { return p[k]; }
                /*}}}*/
                /* {{{ predicates on vectors */
                static inline int
                    vec_cmp(elt const * a, elt const * b, size_t k)
                    {
                        for(size_t i = 0 ; i < k ; i++) {
                            int r = T::cmp(a[i], b[i]);
                            if (r) return r;
                        }
                        return 0;
                    }
                static inline bool
                    vec_is_zero(elt const * p, size_t n) {
                        for(size_t i = 0 ; i < n ; i++)
                            if (!T::is_zero(p[i])) return false;
                        return true;
                    }
                /* }}} */
                /* {{{ assignments on vectors */

                static inline void
                    vec_set_zero(elt * p, size_t n) {
                        /* XXX This _should_ resolve to a memset */
                        for(size_t i = 0 ; i < n ; i++)
                            p[i].zero();
                    }

                static inline void
                    vec_set(elt * q, elt const * p, size_t n)
                    {
                        static_assert(elt::is_flat_storage, "X must be flat");
                        std::copy_n(p, n, q);
                    }

                static inline void
                    vec_neg(elt * q, elt const * p, size_t n)
                    {
                        vec_set(q, p, n);
                    }

                static inline void
                    vec_set_random(elt * p, size_t k, gmp_randstate_ptr rstate) {
                        for(size_t i = 0 ; i < k ; ++i) 
                            T::set_random(p[i], rstate);
                    }
                /* }}} */

                /*{{{ vec addition */
                static inline void
                    vec_add(elt * q, elt const * p, size_t n) {
                        for(size_t i = 0 ; i < n ; i++)
                            T::add(q[i], p[i]);
                    }
                static inline void
                    vec_add_and_reduce(elt * q, elt const * p, size_t n) {
                        T::vec_add(q, p, n);
                    }
                /*}}}*/
                /*{{{ vec simd operations*/
                static inline void
                    vec_simd_set_ui_at(elt * p, size_t k, int v)
                    {
                        T::simd_set_ui_at(p[k / G], k%G, v);
                    }
                static inline void
                    vec_simd_add_ui_at(elt * p, size_t k, int v)
                    {
                        T::simd_add_ui_at(p[k / G], k%G, v);
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
                    vec_simd_find_first_set(elt &, elt const * p, size_t n)
                    {
                        int j = 0;
                        for(size_t i = 0 ; i < n ; ++i, j+=T::simd_groupsize) {
                            elt const & x = T::vec_item(p, i);
                            if (!T::is_zero(x)) {
                                for(size_t k = 0 ; k < T::simd_groupsize ; k++, j++) {
                                    if (T::simd_test_ui_at(x, k))
                                        return j;
                                }
                            }
                        }
                        return -1;
                    }
                /*}}}*/
            void vec_dotprod(elt &, elt const *, elt const *, size_t) const
            {
                /* We do not really want to call this function, as it
                 * mostly exists for the groupsize == 1 case. Over the
                 * binary field, we never have groupsize == 1 anyway
                 */
                abort();
            }

            void vec_addmul_and_reduce(elt *, elt const *, elt const &, size_t) const
            {
                /* same as above */
                abort();
            }

            };

        template<unsigned int G>
            struct gf2 : public gf2_middle<G, gf2<G>> {
                typedef gf2_middle<G, gf2<G>> super;
                // we have nothing specific to do with the ctor.
                // template<typename... Args> gf2(Args&&... args) : super(std::forward<Args>(args)...) {}
                gf2(mpz_srcptr p, unsigned int simd_groupsize) {
                    /* Do the required sanity checks */
                    ASSERT_ALWAYS(mpz_cmp_ui(p, 2) == 0);
                    ASSERT_ALWAYS(simd_groupsize == G);
                }
                static std::string impl_name() {
                    return fmt::format(FMT_STRING("b{}"), G);
                }
            };
    }
    /* expose only what we have in our public interface */
    using details::gf2;
}


#endif	/* ARITH_MOD2_HPP_ */
