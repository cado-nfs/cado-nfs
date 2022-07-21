#ifndef ARITH_MOD2_HPP_
#define ARITH_MOD2_HPP_

#include <algorithm>
#include <array>
#include <cstdint>
#include <iomanip>
#include <ostream>
#include <sstream>
#include <string>
#if defined(HAVE_SSE41) && defined(HAVE_POPCNT)
#include <x86intrin.h>
#endif

#include "cxx_mpz.hpp"
#include "fmt/format.h"
#include "gmp_aux.h"
#include "macros.h"
#include "memory.h" // malloc_aligned
#include "misc.h" // u64_random
#include <gmp.h>

#include "arith-concrete-base.hpp"

namespace arith_mod2 {
namespace details {

/* compile-time table that returns the largest power of two (but at most
 * 128) that divides k */
template<unsigned int k,unsigned int ell = 1, int divide_out = (ell < 128 && !(k & 1))>
struct alignment_divisor;
template<unsigned int k, unsigned int ell> struct alignment_divisor<k,ell,0> : public std::integral_constant<unsigned int, ell> {};
template<unsigned int k, unsigned int ell> struct alignment_divisor<k,ell,1> : public alignment_divisor<k/2,ell*2> {};


/* {{{ layout_traits
 *
 * There are a few things that we do differently depending on whether
 * the base type is constant or variable size.
 */

template<typename elt, bool check = elt::is_constant_size>
struct layout_traits
{};

/* {{{ constant size */
template<typename elt>
struct layout_traits<elt, true>
{
    /* We assume that an element of type E can be operated on with
     * E::number_of_limbs() elements of type E::value_type ; this is
     * typically what happens with a non-zero array, but E might also be
     * a different type, and still allow this sort of mechanism behind
     * the scenes.
     */
    static constexpr unsigned int simd_groupsize()
    {
        return elt::simd_groupsize();
    }
    static constexpr unsigned int number_of_limbs()
    {
        return elt::number_of_limbs();
    }
    typedef typename elt::value_type value_type;

    /*{{{ element layout information */
    static constexpr inline size_t elt_stride() { return sizeof(elt); }
    static inline size_t vec_elt_stride(size_t s) { return s * elt_stride(); }
    /*}}}*/

    /*{{{ allocation / deallocation of (vectors of) elements */
    static inline elt* alloc(size_t k = 1, size_t al = alignment_divisor<sizeof(elt)>::value) {
        return reinterpret_cast<elt *>(malloc_aligned(k * sizeof(elt), al));
    }
    static inline void free(elt* u) {
        free_aligned(reinterpret_cast<void *>(u));
    }
    static inline elt* realloc(elt* u, size_t k0, size_t k, size_t al = alignment_divisor<sizeof(elt)>::value)
    {
        return reinterpret_cast<elt *>(::realloc_aligned(reinterpret_cast<void *>(u), k0 * sizeof(elt), k * sizeof(elt), al));
    }
    /*}}}*/

    /* {{{ setting and setting to zero */
    static inline elt& set(elt& x, elt const& a) { return x = a; }
    static inline void set_zero(elt& x) { x.zero(); }
    /* }}} */

    /* {{{ vector arithmetic and copies */
    static inline elt* vec_subvec(elt* p, size_t k) { return p + k; }

    static inline elt const* vec_subvec(elt const* p, size_t k)
    {
        return p + k;
    }

    static inline void vec_set(elt* q, elt const* p, size_t n)
    {
        if (q < p)
            std::copy_n(p, n, q);
        else
            std::copy_backward(p, p + n, q + n);
    }
    /* }}} */

    layout_traits(unsigned int G) { ASSERT_ALWAYS(G == simd_groupsize()); }
};
/* }}} */

/* {{{ variable size */
template<typename elt>
struct layout_traits<elt, false>
{
    typedef typename elt::value_type X;
    unsigned int K;
    inline unsigned int number_of_limbs() const { return K; }
    inline unsigned int simd_groupsize() const { return K * 64; }
    /* {{{ element layout information */
    /* These are really really weird */
    inline size_t elt_stride() const { return K * sizeof(X); }
    inline size_t vec_elt_stride(size_t s) const { return s * elt_stride(); }
    /* }}} */

    /* {{{ allocation / deallocation of (vectors of) elements */
    inline elt* alloc(size_t k = 1, size_t al = sizeof(X)) const {
        /* Do we want a runtime determination of the alignment ? */
        return reinterpret_cast<elt *>(malloc_aligned(k * K * sizeof(X), al));
    }
    inline void free(elt* u) const {
        free_aligned(reinterpret_cast<void *>(u));
    }
    inline elt* realloc(elt* u, size_t k0, size_t k, size_t al = sizeof(X)) const
    {
        return reinterpret_cast<elt *>(::realloc_aligned(reinterpret_cast<void *>(u), k0 * K * sizeof(X), k * K * sizeof(X), al));
    }
    /* }}} */

    /* {{{ setting and setting to zero */
    inline elt& set(elt& x, elt const& a) const
    {
        std::copy_n(a.data(), number_of_limbs(), x.data());
        return x;
    }
    inline void set_zero(elt& x) const
    {
        std::fill_n(x.data(), number_of_limbs(), 0);
    }
    /* }}} */

    /* {{{ vector arithmetic and copies */
    inline elt* vec_subvec(elt* p, size_t k) const
    {
        return reinterpret_cast<elt*>(p->data() + k * number_of_limbs());
    }

    inline elt const* vec_subvec(elt const* p, size_t k) const
    {
        return reinterpret_cast<elt const*>(p->data() + k * number_of_limbs());
    }

    inline void vec_set(elt* q, elt const* p, size_t n) const
    {
        X* qx = q->data();
        X const* px = p->data();

        if (q < p)
            std::copy_n(px, n * number_of_limbs(), qx);
        else
            std::copy_backward(
              px, px + n * number_of_limbs(), qx + n * number_of_limbs());
    }
    /* }}} */

    layout_traits(unsigned int G)
      : K(G / elt::value_type_bits)
    {}
};
/* }}} */
/* }}} */

/* This base will be defined later on. */
template<unsigned int G, typename T>
struct gf2_base;

/* This can host type-specific overrides on top of the default routines.
 */
template<unsigned int G, typename T>
struct gf2_override : public gf2_base<G, T>
{
    template<typename... Args>
    gf2_override(Args&&... args)
      : gf2_base<G, T>(std::forward<Args>(args)...)
    {}
};

/* The base type. This can be overridden later on. */
template<unsigned int G>
struct gf2_base_type
  : public std::array<uint64_t, G / 64>
  , public arith_concrete_base::elt
{
    static constexpr const bool is_constant_size = true;
    static constexpr const unsigned int K = G / 64;
    static constexpr unsigned int simd_groupsize() { return G; }
    static constexpr unsigned int number_of_limbs() { return K; }
    typedef std::array<uint64_t, K> elt;
    static constexpr const unsigned int value_type_bits = 64;
    /*
    template<typename T> struct compatible_with {
        static constexpr const bool value = is_flat_storage && (sizeof(elt) ==
    sizeof(T));
    };
    */
    template<typename... Args>
    gf2_base_type(Args&&... args)
      : elt(std::forward<Args>(args)...)
    {}
    inline void zero() { std::fill_n(elt::begin(), K, 0); }
};

/* This one is for the variable width case.
 *
 * It's really a bit special because the underlying element type has no
 * meaning at compilation time. It only depends on runtime data.
 *
 * We could be tempted to say that we're working with a 0-size array, but
 * that doesn't work because data() is allowed to return NULL in that
 * case.
 */
template<>
struct gf2_base_type<0> : public arith_concrete_base::elt
{
    static constexpr const bool is_constant_size = false;
    typedef uint64_t value_type;
    static constexpr const unsigned int value_type_bits = 64;

    value_type* data() { return reinterpret_cast<value_type*>(this); }
    value_type const* data() const
    {
        return reinterpret_cast<value_type const*>(this);
    }
};

#if defined(HAVE_SSE41) && defined(HAVE_POPCNT)
template<>
struct gf2_base_type<128> : public arith_concrete_base::elt
{
    static constexpr const bool is_constant_size = true;
    static constexpr unsigned int simd_groupsize() { return 128; }

    /* This is only so that the trivial functions can work on this type
     * easily */
    typedef uint64_t value_type;
    static constexpr const unsigned int value_type_bits = 64;
    static constexpr unsigned int number_of_limbs() { return 2; }

    /* It's slightly awkward that we use composition and not inheritance
     * here. */
    __m128i x;
    typedef uint64_t value_type;
    value_type* data() { return reinterpret_cast<value_type*>(&x); }
    value_type const* data() const
    {
        return reinterpret_cast<value_type const*>(&x);
    }
    void zero() { x = _mm_setzero_si128(); }
};

template<typename T>
struct gf2_override<128, T> : public gf2_base<128, T>
{
    template<typename... Args>
    gf2_override(Args&&... args)
      : gf2_base<G, T>(std::forward<Args>(args)...)
    {}
    /* {{{ assignments */
    static inline void set_random(elt& x, gmp_randstate_ptr rstate)
    {
        __m64 lo = _mm_cvtsi64_m64(u64_random(rstate));
        __m64 hi = _mm_cvtsi64_m64(u64_random(rstate));
        x = _mm_setr_epi64(lo, hi);
    }
    /* }}} */
    static inline void add(elt& dst, elt const& a, elt const& b)
    {
        dst.x = _mm_xor_si128(a.x, b.x);
    }
    /*{{{ simd*/
    static inline int simd_hamming_weight(elt const& p)
    {
        return _mm_popcnt_u64(_mm_extract_epi64(p.x, 0)) +
               _mm_popcnt_u64(_mm_extract_epi64(p.x, 1));
    }
    /*}}}*/
};

#endif

template<unsigned int G, typename T>
struct gf2_base
  : public arith_concrete_base
  , public layout_traits<gf2_base_type<G>>
{
    template<typename... Args>
    gf2_base(Args&&... args)
      : layout_traits<gf2_base_type<G>>(std::forward<Args>(args)...)
    {}
    static constexpr const bool is_characteristic_two = true;

    typedef gf2_base_type<G> elt;
    typedef layout_traits<elt> layout;
    static inline void simd_set_ui_at(elt& p, size_t k, int v)
    {
        uint64_t* xp = p.data();
        uint64_t mask = ((uint64_t)1) << (k % 64);
        xp[k / 64] =
          (xp[k / 64] & ~mask) | ((((uint64_t)v) << (k % 64)) & mask);
    }
    static inline int simd_test_ui_at(elt const& p, size_t k)
    {
        uint64_t const* xp = p.data();
        uint64_t mask = ((uint64_t)(1)) << (k % 64);
        return (xp[k / 64] & mask) != 0;
    }
    static inline void simd_add_ui_at(elt& p, size_t k, int v)
    {
        uint64_t* xp = p.data();
        uint64_t mask = ((uint64_t)(v & 1)) << (k % 64);
        xp[k / 64] ^= mask;
    }
    inline int simd_hamming_weight(elt const& p) const
    {
        int w = 0;
        T const* tx = static_cast<T const*>(this);
#if GNUC_VERSION_ATLEAST(3, 4, 0)
        const unsigned long* xp =
          reinterpret_cast<const unsigned long*>(p.data());
        static_assert(!T::elt::is_constant_size ||
                        sizeof(typename T::elt) % sizeof(unsigned long) == 0,
                      "bad types");
        // parentheses in the divisor are mandatory here, see
        // https://gcc.gnu.org/pipermail/gcc-patches/2020-September/553888.html
        for (size_t b = 0; b < tx->elt_stride() / (sizeof(unsigned long));
             b++) {
            w += __builtin_popcountl(xp[b]);
        }
#else
        int tab[16] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
        uint8_t const* xp = reinterpret_cast<uint8_t const*>(p.data());
        for (size_t b = 0; b < tx->elt_stride(); b++) {
            w += tab[xp[b] & 15] + tab[xp[b] >> 4];
        }
#endif
        return w;
    }

    /*{{{ predicates */
    inline bool is_zero(elt const& x) const
    {
        T const* tx = static_cast<T const*>(this);
        unsigned int K = tx->number_of_limbs();
        for (unsigned int i = 0; i < K; i++) {
            if (x.data()[i])
                return false;
        }
        return true;
    }

    inline int cmp(elt const& x, unsigned long a) const
    {
        T const* tx = static_cast<T const*>(this);
        unsigned int K = tx->number_of_limbs();
        int r = (x.data()[0] > a) - (a > x.data()[0]);
        if (r)
            return r;
        for (unsigned int i = 1; i < K; i++)
            if (x.data()[i])
                return 1;
        return 0;
    }

    inline int cmp(elt const& x, elt const& y) const
    {
        T const* tx = static_cast<T const*>(this);
        unsigned int K = tx->number_of_limbs();
        for (unsigned int i = 0; i < K; i++) {
            int r = (x.data()[i] > y.data()[i]) - (y.data()[i] > x.data()[i]);
            if (r)
                return r;
        }
        return 0;
    }
    /*}}}*/

    /* {{{ assignments */
    inline elt& set(elt& x, elt const& a) const { return layout::set(x, a); }
    inline void set_zero(elt& x) const { return layout::set_zero(x); }

    inline elt& set(elt& x, unsigned long a) const
    {
        T const* tx = static_cast<T const*>(this);
        tx->set_zero(x);
        x.data()[0] = a;
        return x;
    }

    inline elt& set(elt& x, cxx_mpz const& a) const
    {
        T const* tx = static_cast<T const*>(this);
        tx->set_zero(x);
        x.data()[0] = mpz_get_ui(a);
        return x;
    }

    inline elt& neg(elt& x, elt const& a) const
    {
        T const* tx = static_cast<T const*>(this);
        return tx->set(x, a);
    }

    inline void set_random(elt& x, gmp_randstate_ptr rstate) const
    {
        T const* tx = static_cast<T const*>(this);
        unsigned int K = tx->number_of_limbs();
        for (unsigned int i = 0; i < K; i++)
            x.data()[i] = u64_random(rstate);
    }

    inline void stream_store(elt* dst, elt const& src) const { set(*dst, src); }

    static inline void reduce(elt&) {}
    /* }}} */

    /*{{{ addition, at the element level */
    inline void add(elt& dst, elt const& a, elt const& b) const
    {
        T const* tx = static_cast<T const*>(this);
        unsigned int K = tx->number_of_limbs();
        for (unsigned int i = 0; i < K; i++)
            dst.data()[i] = a.data()[i] ^ b.data()[i];
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
    inline std::ostream& cxx_out(std::ostream& o, elt const& x) const
    {
        T const* tx = static_cast<T const*>(this);
        unsigned int K = tx->number_of_limbs();
        for (unsigned int i = 0; i < K; i++) {
            // We _know_ that we're dealing with a 64-bit type
            // anyway, so what follows is unnecessary.
            // o << fmt::format(FMT_STRING("{0:0{1}}"), x.data()[i],
            // sizeof(x.data()[i]) * CHAR_BIT / 4);
            o << fmt::format(FMT_STRING("{0:016x}"), x.data()[i]);
        }
        return o;
    }
    inline std::ostream& write(std::ostream& o, elt const& x)
    {
        T const* tx = static_cast<T const*>(this);
        return o.write(reinterpret_cast<const char*>(x.data()),
                       tx->elt_stride());
    }
    inline int fread(FILE* f, elt& x) const
    {
        T const* tx = static_cast<T const*>(this);
        int ret =
          ::fread(reinterpret_cast<char*>(x.data()), tx->elt_stride(), 1, f);
        if (ret < 1)
            return ret;
        return ret;
    }
    inline int fscan(FILE* f, elt& x) const
    {
        T const* tx = static_cast<T const*>(this);
        unsigned int K = tx->number_of_limbs();
        char buf[K * 16 + 1];
        int ret = fscanf(f, "%*s", sizeof(buf), buf);
        if (!ret)
            return 0;
        char* b = buf;
        for (unsigned int i = 0; i < K; i++, b += 16) {
            std::string s(b, b + 16);
            std::istringstream is(s);
            is >> std::hex >> x.data()[i];
        }
        return ret;
    }
    /*}}}*/
};

/*************************************************************/
/* code below this point does not dive into the details of
 * the elt type.
 */

template<unsigned int G, typename T>
struct gf2_middle : public gf2_override<G, T>
{
    typedef gf2_override<G, T> super;
    using typename super::elt;

    template<typename... Args>
    gf2_middle(Args&&... args)
      : super(std::forward<Args>(args)...)
    {}
    mpz_srcptr characteristic() const
    {
        static mp_limb_t limbs[1] = { 2 };
        static __mpz_struct a = { 1, 1, limbs };
        return &a;
    }

    /*{{{ a few trivial proxies */
    using super::add;
    inline void add(elt& dst, elt const& b) const
    {
        T const* tx = static_cast<T const*>(this);
        tx->add(dst, dst, b);
    }
    inline void add_and_reduce(elt& dst, elt const& b) const
    {
        T const* tx = static_cast<T const*>(this);
        tx->add(dst, dst, b);
    }
    inline void sub_and_reduce(elt& dst, elt const& b) const
    {
        T const* tx = static_cast<T const*>(this);
        tx->add(dst, dst, b);
    }
    /*}}}*/

    /* TODO: mul / addmul interface. What do we have to do?
     */

    /*{{{ accessors inside vectors */
    inline elt& vec_item(elt* p, size_t k) const
    {
        T const* tx = static_cast<T const*>(this);
        return *tx->vec_subvec(p, k);
    }

    inline elt const& vec_item(elt const* p, size_t k) const
    {
        T const* tx = static_cast<T const*>(this);
        return *tx->vec_subvec(p, k);
    }
    /*}}}*/

    /* {{{ predicates on vectors */
    inline int vec_cmp(elt const* a, elt const* b, size_t k) const
    {
        T const* tx = static_cast<T const*>(this);
        for (size_t i = 0; i < k; i++) {
            int r = tx->cmp(tx->vec_item(a, i), tx->vec_item(b, i));
            if (r)
                return r;
        }
        return 0;
    }
    inline bool vec_is_zero(elt const* p, size_t n) const
    {
        T const* tx = static_cast<T const*>(this);
        for (size_t i = 0; i < n; i++)
            if (!tx->is_zero(tx->vec_item(p, i)))
                return false;
        return true;
    }
    /* }}} */
    /* {{{ assignments on vectors */

    inline void vec_set_zero(elt* p, size_t n) const
    {
        T const* tx = static_cast<T const*>(this);
        /* XXX This _should_ resolve to a memset */
        for (size_t i = 0; i < n; i++)
            tx->set_zero(tx->vec_item(p, i));
    }

    inline void vec_neg(elt* q, elt const* p, size_t n) const
    {
        T const* tx = static_cast<T const*>(this);
        tx->vec_set(q, p, n);
    }

    inline void vec_set_random(elt* p, size_t k, gmp_randstate_ptr rstate) const
    {
        T const* tx = static_cast<T const*>(this);
        for (size_t i = 0; i < k; ++i)
            tx->set_random(tx->vec_item(p, i), rstate);
    }
    /* }}} */

    /*{{{ vec addition */
    inline void vec_add(elt* q, elt const* p, size_t n) const
    {
        T const* tx = static_cast<T const*>(this);
        for (size_t i = 0; i < n; i++)
            tx->add(tx->vec_item(q, i), tx->vec_item(p, i));
    }
    inline void vec_add_and_reduce(elt* q, elt const* p, size_t n) const
    {
        T const* tx = static_cast<T const*>(this);
        tx->vec_add(q, p, n);
    }
    /*}}}*/
    /*{{{ vec simd operations*/
    static inline void vec_simd_set_ui_at(elt* p, size_t k, int v)
    {
        // in truth, we know that everything is
        // contiguous, so it makes no sense to do a
        // division by the group size.
        // T::simd_set_ui_at(p[k / G], k%G, v);
        T::simd_set_ui_at(*p, k, v);
    }
    static inline void vec_simd_add_ui_at(elt* p, size_t k, int v)
    {
        // T::simd_add_ui_at(p[k / G], k%G, v);
        T::simd_add_ui_at(*p, k, v);
    }
    static inline int vec_simd_test_ui_at(elt const* p, size_t k, int v)
    {
        return T::simd_test_ui_at(*p, k, v);
    }

    inline int vec_simd_hamming_weight(elt const* p, size_t n) const
    {
        T const* tx = static_cast<T const*>(this);
        int r = 0;
        for (size_t i = 0; i < n; ++i)
            r += tx->simd_hamming_weight(tx->vec_item(p, i));
        return r;
    }
    inline int vec_simd_find_first_set(elt&, elt const* p, size_t n) const
    {
        T const* tx = static_cast<T const*>(this);
        int j = 0;
        for (size_t i = 0; i < n; ++i, j += tx->simd_groupsize()) {
            elt const& x = tx->vec_item(p, i);
            if (!tx->is_zero(x)) {
                for (size_t k = 0; k < tx->simd_groupsize(); k++, j++) {
                    if (T::simd_test_ui_at(x, k))
                        return j;
                }
            }
        }
        return -1;
    }
    /*}}}*/
    void vec_add_dotprod(elt&, elt const*, elt const*, size_t) const
    {
        /* We do not really want to call this function, as it
         * mostly exists for the groupsize == 1 case. Over the
         * binary field, we never have groupsize == 1 anyway
         */
        abort();
    }

    void vec_addmul_and_reduce(elt*, elt const*, elt const&, size_t) const
    {
        /* same as above */
        abort();
    }
};

template<unsigned int G>
struct gf2 : public gf2_middle<G, gf2<G>>
{
    typedef gf2_middle<G, gf2<G>> super;
    // we have nothing specific to do with the ctor.
    // template<typename... Args> gf2(Args&&... args) :
    // super(std::forward<Args>(args)...) {}
    gf2(mpz_srcptr p, unsigned int simd_groupsize)
      : super(simd_groupsize)
    {
        /* Do the required sanity checks. The check on the group size is
         * done manny layers up the chain, in layout_traits */
        ASSERT_ALWAYS(mpz_cmp_ui(p, 2) == 0);
    }
    static std::string impl_name()
    {
        if (G == 0)
            return "bz";
        return fmt::format(FMT_STRING("b{}"), G);
    }
};
}
/* expose only what we have in our public interface */
using details::gf2;
}

#endif /* ARITH_MOD2_HPP_ */
