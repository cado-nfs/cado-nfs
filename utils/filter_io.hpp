#ifndef CADO_FILTER_IO_HPP
#define CADO_FILTER_IO_HPP

#include "cado_config.h"

#include <cstdint>
#include <cstddef>
#include <climits>

#include <string>
#include <vector>
#include <type_traits>
#include <algorithm>
#include <array>
#include <utility>
#include <functional>
#include <istream>
#include <barrier>
#include <memory>
#include <condition_variable>
#include <thread>
#include <mutex>
#include <variant>

#include "fmt/base.h"
#include "fmt/format.h"

#ifdef HAVE_NANOSLEEP
#include <ctime>
#endif

#include "timing.h"
#include "typedefs.h"
#include "vector_with_cache.hpp"
#include "runtime_numeric_cast.hpp"
#include "macros.h"
#include "cxx_mpz.hpp"
#include "cado_compile_time_hacks.hpp"
#include "cado_type_traits.hpp"
#include "ringbuf.hpp"
#include "utils_cxx.hpp"
#include "stats.h"


#define RELATION_MAX_BYTES 4096

/* Initial size of primes_data array in earlyparsed_relation_s,
   If more than NB_PRIMES_OPT is needed (should be rare), *primes is
   allocated
*/
#define NB_PRIMES_OPT 31



/* Size of relation buffer between parsing & processing.  */
#define SIZE_BUF_REL (1U << 15U)

namespace cado::relation_building_blocks {
    /* {{{ hard-coded lookup table for hex and decimal input */
    static constexpr unsigned char Z = 255;
    static constexpr unsigned char const ugly[256] = {
        Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        /* digits, 0x30 - 0x39 */
        0, 1, 2, 3, 4, 5, 6, 7, 8, 9, Z, Z, Z, Z, Z, Z,
        /* capitals, 0x41 - 0x46 */
        Z,10,11,12,13,14,15, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        /* lowercase, 0x61 - 0x66 */
        Z,10,11,12,13,14,15, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
        Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z, Z,
    };
    /* }}} */

    /* The generic relation type below supersedes the old
     * earlyparsed_relation and friends
     *
     * The general grammar is as follows. {{{
     * Each block is documented later on in this file.
     *
     * <relation> = <ab> <sides>? <primelist> <sm>?
     *
     * <ab> = <integer:a> , <integer:b>
     * <sides> = @ <decimal_integer:side0> , <decimal_integer:side1>
     * <primelist> = : <prime>*
     * <prime> = <deprecated_optional_minus> <hex_integer:p> <exponent_info>?
     * <deprecated_optional_minus> = -?
     * <exponent_info> = / <abbreviated_polynomial_in_s>
     * <abbreviated_polynomial_in_s> = ( -? <term> ) ( [+-] <term> )*
     * <term> = [0-9]* (*? s+ (^ <decimal_integer>)?)*
     * <decimal_integer> = [0-9]+
     * <hex_integer> = [0-9a-f]+
     * <sm> = : <decimal_integer> (, <decimal_integer>)*
     *
     * where integer can be base 10 or base 16. Minus signs are only allowed
     * when explicitly marked. }}}
     */

    /* rel_number {{{
     * ==========
     * all relations have a number. No parsing is involved here, but
     * we need a terminating item in the class hierarchy. The number
     * gets put by the filtering layer itself.
     * 
     * This parser returns 0. It means that it is *NOT* legit for parsers
     * in the middle of the line, such as primelist_block or
     * primecount_block, to be built atop rel_number. Use ab_ignore in
     * this case.
     */
    struct rel_number {
        size_t num = 0;
        template<typename iterator>
        int parse(iterator &) { return 0; }
    };
// }}}

    /* ab block {{{
     * ========
     *
     * ab_type can be either
     *      uint64_t, in which case a is uint64_t and b is int64_t
     *      cxx_mpz
     *
     * the base in which a,b are expected to be written is decided at compile
     * time, although it is not actually reflected in the data layout, of
     * course.
     */
    template<typename ab_type, unsigned int base>
        struct ab_block : public rel_number
        { /* {{{ */
            using a_type = cado::make_signed_t<ab_type>;
            using b_type = cado::make_unsigned_t<ab_type>;

            a_type a;
            b_type b;
            std::array<int, 2> active_sides = { 0, 1 };

            /* {{{ parse */
            /* read the full block that contains a,b (in the requested
             * base), + potentially the sides. The first character that
             * is outside the block (a priori a colon) is returned. it0
             * is set to point to the character that *follows* that
             * colon.
             */
            template<typename iterator>
                int parse(iterator & it0)
                requires requires { *it0++ == '\0'; }
            {
                a = 0;
                b = 0;
                using cado_math_aux::log2_ct;
                static_assert(1U <= base && base <= 16U, "base should be in [1, 16]");

                unsigned long v;
                iterator it = it0;

                int c = *it++;
                bool negative = false;
                if (c == '-') {
                    negative = true;
                    c = *it++;
                }

                cado::make_unsigned_t<ab_type> w;

                for (w = 0 ; (v = ugly[c]) < base; w += v, c = *it++) {
                    if constexpr ((base & (base - 1)) == 0) {
                        w <<= log2_ct(base);
                    } else {
                        w *= base;
                    }
                }

                if (negative)
                    a = -runtime_numeric_cast<a_type>(w);
                else
                    a = runtime_numeric_cast<a_type>(w);

                ASSERT_ALWAYS(c == ',');
                c = *it++;

                for (w = 0 ; (v = ugly[c]) < base; w += v, c = *it++) {
                    if constexpr ((base & (base - 1)) == 0) {
                        w <<= log2_ct(base);
                    } else {
                        w *= base;
                    }
                }
                b = runtime_numeric_cast<b_type>(w);

                if (c == '@') {
                    int v, w = 0;
                    c = *it++;
                    for (w = 0 ; (v = ugly[c]) < 10; w += v, c = *it++)
                        w *= 10;
                    active_sides[0] = w;
                    ASSERT_ALWAYS(c == ',');
                    c = *it++;
                    for (w = 0 ; (v = ugly[c]) < 10; w += v, c = *it++)
                        w *= 10;
                    active_sides[1] = w;
                    if (active_sides != std::array<int, 2>({0,1})) {
                        fmt::print(stderr, "# warning, relation with non-standard sides\n");
                    }
                } else {
                    active_sides[0] = 0;
                    active_sides[1] = 1;
                }
                it0 = it;
                return c;
            } /* }}} */
            template<typename OutputIt>
                requires (base == 16)
            auto format_to(OutputIt out) const -> OutputIt
            {
                /* output something that makes it clear that there's data
                 * that was removed */
                return fmt::format_to(out, "{:x},{:x}", a, b);
            }
            template<typename OutputIt>
                requires (base == 10)
            auto format_to(OutputIt out) const -> OutputIt
            {
                /* output something that makes it clear that there's data
                 * that was removed */
                return fmt::format_to(out, "{},{}", a, b);
            }
        }; /* }}} */
    /* }}} */

    /* ab_ignore {{{
     * =========
     *
     * This recognizes the exact same thing as ab_block, but here
     * we do _not_ interpret the data, we just skip over it
     */
    template<unsigned int base>
        struct ab_ignore : public rel_number
        { /* {{{ */

            /* {{{ parse */
            template<typename iterator>
                int parse(iterator & it0)
                requires requires { *it0++ == '\0'; }
            {
                using cado_math_aux::log2_ct;
                static_assert(1U <= base && base <= 16U, "base should be in [1, 16]");

                iterator it = it0;

                int c = *it++;
                if (c == '-')
                    c = *it++;

                for ( ; ugly[c] < base; c = *it++) ;

                ASSERT_ALWAYS(c == ',');
                c = *it++;

                for ( ; ugly[c] < base; c = *it++) ;

                if (c == '@') {
                    c = *it++;
                    for ( ; ugly[c] < 10; c = *it++) ;
                    ASSERT_ALWAYS(c == ',');
                    c = *it++;
                    for (; ugly[c] < 10; c = *it++) ;
                }
                it0 = it;
                return c;
            } /* }}} */
            template<typename OutputIt>
            auto format_to(OutputIt out) const -> OutputIt
            {
                /* output something that makes it clear that there's data
                 * that was removed */
                fmt::format_to(out, "*");
                return out;
            }

        }; /* }}} */
    /* }}} */

    /* primelist block {{{
     * ===============
     *
     * Ideals in relations exist in multiple flavors. Both prime numbers and
     * prime indices are given in hex, in all formats. Sorting is not
     * required.
     *
     *  - When primes and not prime indices are given, multiple sides are
     *  possible, separated by colons (a line with n sides will always
     *  have n colons in total). The prime "zero" is (of course) not
     *  allowed. An empty side is entirely possible.
     *  FIXME: the way this interacts with active_sides is unclear. I
     *  surmise that if active_sides is defined, then we should determine
     *  the side by dereferencing that array. But then of course it
     *  implies that we must have exactly two colon-separated blocks.
     *
     * Primes and prime indices may have exponents. Those can be
     * specified with the following mechanism.
     *
     *  - unary representation of positive numbers is allowed. I.e., if
     *  we have 2,2,2, it means the prime number (or the prime ideal with
     *  index 2), three times.
     *
     *  - negative signs before the index negate the exponent (XXX: removed)
     * 
     *  - an optional exponent can be specified after the prime (or prime
     *  index) after the / suffix, with the syntax described below.
     *
     */

    template<typename prime_t>
    struct import_p_or_h_type {
        using failed_type = void;
    };

    template<typename prime_t>
    requires requires { typename prime_t::p_type; }
    struct import_p_or_h_type<prime_t>
    {
        using p_type = prime_t::p_type;
    };
    template<typename prime_t>
    requires requires { typename prime_t::index_type; }
    struct import_p_or_h_type<prime_t>
    {
        using index_type = prime_t::index_type;
    };

    template<typename prime_t, typename parent>
        requires
            requires { typename prime_t::p_or_h_type; }
         && requires { typename prime_t::e_type; }
    struct primes_block
        : public import_p_or_h_type<prime_t>
        , public parent
    { /* {{{ */
        using prime_type = prime_t;
        using p_or_h_type = prime_t::p_or_h_type;
        using e_type = prime_t::e_type;

        parent const & parent_block() const { return *this; }

        static_assert(!std::is_same_v<parent, rel_number>,
                "use ab_ignore instead of rel_number");

        vector_with_cache<prime_type, NB_PRIMES_OPT> primes;

        /* the number of sides is of course just one for indexed
         * relations. We could choose to ignore it, but it's probably not
         * worth the trouble.
         */
        int nsides = 0;

        /* This is the "simple" exponent parser. */
        template<typename iterator>
            int parse_exponent(int c,
                    iterator & it, 
                    typename prime_type::e_type & e)
            requires (std::is_integral_v<typename prime_type::e_type>
                    && requires { *it++ == '\0'; })
            { /* {{{ */
                e = 1;
                if (c == '/') {
                    c = *it++;
                    bool negative = false;
                    if (c == '-') {
                        negative = true;
                        c = *it++;
                    }
                    int v, w;

                    e = typename prime_type::e_type(0);
                    /* a leading 0 is always a bug */
                    ASSERT_ALWAYS(ugly[c]);
                    for (w = 0 ; (v = ugly[c]) < 10 ; w += v, c = *it++) w *= 10;
                    if (w == 0) w = 1; /* / means 1, and /- means -1.  */
                    e = negative ? -w : w;
                }
                return c;
            } /* }}} */

        template<typename iterator>
            int parse_exponent(int c,
                    iterator & it, 
                    typename prime_type::e_type & e)
            requires (prime_type::e_type::is_exponent_vector
                    && requires { *it++ == '\0'; })
            { /* {{{ */
                e = typename prime_type::e_type(1);
                if (c == '/') {
                    c = *it++;
                    c = e.parse(c, it);
                }
                return c;
            } /* }}} */

        void append_or_collate(p_or_h_type p_or_h, e_type e, int side)
            requires requires { decltype(prime_type::side)(); }
        {
            if (!primes.empty() && p_or_h == primes.last().p_or_h() && primes.last().side == side)
                primes.last().e += e;
            else {
                primes.emplace_back(p_or_h, e, side);
            }
        }

        void append_or_collate(p_or_h_type h, e_type e, int)
            requires (!requires { decltype(prime_type::side)(); })
        {
            if (!primes.empty() && h == primes.last().h)
                primes.last().e += e;
            else {
                primes.emplace_back(h, e);
            }
        }

        /* next_delim is an obscure hack. We want this same code to be called
         * even when an SM block is expected after the indexed relation. In
         * which case obviously we mustn't parse multiple sides.
         *
         * Note that this approach of using the same code path, both for
         * sieve relations and indexed relations, clearly has its limits.
         * On the other hand there is so little to adjust that
         * duplicating 50+ lines of code is not worth it.
         */
        template<typename iterator>
            int parse(iterator & it0, int delim_after_primes = '\n')
            requires requires { *it0++ == '\0'; }
        { /* {{{ */
            primes.clear();

            int c = parent::parse(it0);
            iterator it = it0;

            /* the side is not meaningful for prime indices, and we'll ignore
             * it in that case.
             */
            int side = -1;
            bool sorted = true;

            typename prime_type::p_or_h_type last_prime = 0;

            /* the very first delimiter is _before_ primes, so it's always a
             * colon. So we must make sure that we accept the column on the
             * first pass.
             *
             * We add a provision to break at null characters. It isn't going
             * to be used much, since implementations will probably have
             * newline-terminated lines. But we use it for testing.
             */
            for(bool first = true ; first || (c && c != delim_after_primes && c != '\n') ; first = false) {
                typename prime_type::p_or_h_type pr;
                if (c == ':') {
                    sorted = false;
                    last_prime = 0;
                    /* Empty sides are allowed. So we must check for
                     * separators another time (another colon, or a \n)
                     * otherwise we would push the prime 0 to the list...
                     *
                     */
                    c = *it++, side++, first=false;;
                    for( ; c != delim_after_primes && c == ':' ; c = *it++)
                        side++;
                    if (c == delim_after_primes || c == '\n')
                        break;
                } else {
                    ASSERT_ALWAYS(c == ',');
                    c = *it++;
                }
                /* XXX deprecated_optional_minus
                 * filter_galois used to produce this.
                 */
                bool negative = false;
                if (c == '-') {
                    negative = true;
                    c = *it++;
                }
                /* XXX end deprecated block */

                unsigned long v, w;
                for (w = 0 ; (v = ugly[c]) < 16 ; w += v, c = *it++) w <<= 4;
                pr = w;
                typename prime_type::e_type e;
                c = parse_exponent(c, it, e);

                /* XXX more stuff to remove eventually */
                if (negative)
                    e = -e;
                /* XXX end deprecated block */

                sorted = sorted && pr >= last_prime;

                append_or_collate(pr, e, side);
            }
            if (!sorted)
                sort_and_compress();
            it0 = it;
            nsides = side + 1;
            return c;
        } /* }}} */

        template<typename OutputIt>
            requires requires { prime_type().side; }
        auto format_to(OutputIt out) const -> OutputIt
        {
            out = parent::format_to(out);
            /* it's annoying but we have no guarantee that our primes are
             * sorted according to the side! */
            std::vector<std::vector<prime_type>> Ps(nsides);
            for(auto const & p : primes)
                Ps[p.side].push_back(p);
            for(auto const & s : Ps) {
                *out++ = ':';
                fmt::format_to(out, "{}", join(s, ","));
            }
            return out;
        }

        template<typename OutputIt>
            requires requires { prime_type().h; }
        auto format_to(OutputIt out) const -> OutputIt
        {
            out = parent::format_to(out);
            *out++ = ':';
            fmt::format_to(out, "{}", join(primes.begin(), primes.end(), ","));
            return out;
        }

        void sort_and_compress() {
            std::vector<prime_type> A;
            for(auto const & pse : primes)
                A.push_back(pse);
            std::ranges::sort(A);
            primes.clear();
            size_t j = 0;
            for(size_t i = 0; i < A.size(); i++) {
                if (j && A[i].p_or_h() == A[j-1].p_or_h()) {
                    A[j-1].e += A[i].e;
                    if (A[j-1].e == 0)
                        j--;
                } else {
                    A[j] = A[i];
                    j++;
                }
            }
            for(size_t i = 0 ; i < j ; i++)
                primes.push_back(A[i]);
        }
        /*
         * <primelist> = : <prime>*
         * <prime> = <deprecated_optional_minus> <hex_integer:p> <exponent_information>?
         * <deprecated_optional_minus> = -?
         * <exponent_information> = / <abbreviated_polynomial_in_s>
         * <abbreviated_polynomial_in_s> = ( -? <term> ) ( [+-] <term> )*
         * <term> = [0-9]* (*? s+ (^ <decimal_integer>)?)*
         * <decimal_integer> = [0-9]+
         * <hex_integer> = [0-9a-f]+
         */
    }; /* }}} */
    /* }}} */

    /* primecount block {{{
     * ====================
     *
     * The same syntax as in primelist is recognized, but here we only
     * _count_ the number of primes.
     *
     * There is a very, very significant caveat, though: the primes are
     * counted as they are, with multiplicities if they happen to be
     * repeated! For consistent results, it therefore makes sense to use
     * this block only on collated prime lists, which is the case in the
     * normal workflow (we only need it duting purge, as a stats
     * measure).
     *
     * Since we do not care about parsing for real, we _only_ count
     * commas. We don't parse the exponents.
     */

    template<typename parent>
    struct primecount_block
        : public parent
    { /* {{{ */
        size_t weight = 0;

        static_assert(!std::is_same_v<parent, rel_number>,
                "use ab_ignore instead of rel_number");

        /* we follow the exact same control flow as in primelist_block.  */
        template<typename iterator>
            int parse(iterator & it0, int delim_after_primes = '\n')
            requires requires { *it0++ == '\0'; }
        { /* {{{ */
            int c = parent::parse(it0);
            iterator it = it0;

            for( ; c && c != ':' ; c = *it++) ;

            for(bool first = true ; first || (c && c != delim_after_primes) ; first = false) {
                if (c == ':') {
                    /* skip over sides */
                    for( ; c != delim_after_primes && c == ':' ; c = *it++) ;
                    if (c == delim_after_primes)
                        break;
                } else {
                    ASSERT_ALWAYS(c == ',');
                    c = *it++;
                }
                /* XXX deprecated_optional_minus */
                if (c == '-')
                    c = *it++;
                /* XXX end deprecated block */

                for ( ; c && c != ',' && c != ':' && c != '\n' ; c = *it++) ;
                weight++;
            }
            it0 = it;
            return c;
        } /* }}} */
    }; /* }}} */
    /* }}} */

    /* exponent block {{{
     * ==============
     *
     * Exponents are given as a polynomial expression in the
     * indeterminate s, after the / prefix, and terminated at the next
     * comma, colon, end-of-line, or null character. For example the
     * relation below has a few exponents.
     *
     * 17,42:b,13/s^2,13/s^2,25:7,11,35/s^4-s
     *
     * in which the specified exponents represent: on side 0, and twice
     * the result of applying sigma^2 to the ideal of norm 19, and on
     * side 1, the result of applying s^4-s to the ideal of norm 53.
     * Specifying /s^0 or /1 is abbreviated by just not having anything
     * as exponent.
     *
     * The exponent expression can be abbreviated in multiple ways: all *
     * signs are optional, and a lone 1 can always be removed.
     *
     * E.g., the following expressions are all equivalent.
     *
     * 11,2a:3,4/-,5,8/-s+2ss,b,10/s^4
     * 11,2a:3,4/-1,5,8/-s+2*s*s,b,10/s^2*s^2
     * 11,2a:3,4/-1,5,8/-s+s*s+s*s,b,10/s^2*s^2
     * 11,2a:3/,4/-s^0,5,8/-1*s+s*s,8/s*s,b,10/s^2*s^2
     *
     * Note that depending on the exponent type, not all exponent
     * expressions are supported.
     */

    template<int n> struct exponent_vector /* {{{ */
        : public std::array<int, n>
    {
        static constexpr bool is_exponent_vector = true;
        explicit exponent_vector(int c) : std::array<int, n> { c, } {}
        exponent_vector() : std::array<int, n> { 0, } {}
        std::array<int, n> const & super() const { return *this; }
        auto operator==(exponent_vector const & o) const {
            return std::equal(
                    super().begin(), super().end(), o.begin(), o.end());
        }
        auto operator<=>(exponent_vector const & o) const {
            return std::lexicographical_compare_three_way(
                    super().begin(), super().end(), o.begin(), o.end());
        }
        auto operator==(int c) const { return operator==(exponent_vector<n>(c)); }
        auto operator<=>(int c) const { return operator<=>(exponent_vector<n>(c)); }
        exponent_vector& operator+=(exponent_vector const & f)
        {
            for(int j = 0 ; j < n ; j++)
                (*this)[j] += f[j];
            return *this;
        }
        exponent_vector operator-() const {
            exponent_vector res;
            for(int j = 0 ; j < n ; j++)
                res[j] = -(*this)[j];
            return res;
        }
        template<typename iterator> int parse(int c, iterator & it) /* {{{ */
            requires requires { *it++ == '\0'; }
        {
            /* we recognize this grammar
             * <abbreviated_poly_in_s> = ( -? <term> ) ( [+-] <term> )*
             * <term> = [0-9]* (*? s+ (^ <decimal_integer>)?)*
             */
            for(int j = 0 ; j < n ; (*this)[j++] = 0) ;
            for(bool first = true ; ; first = false) {
                bool negative = false;
                if (c == '-') {
                    negative = true;
                    c = *it++;
                } else if (!first) {
                    if (c == '+') {
                        c = *it++;
                    } else {
                        return c;
                    }
                }
                int j = 0;
                int v, w;
                /* a leading 0 is always a bug */
                ASSERT_ALWAYS(ugly[c]);
                for (w = 0 ; (v = ugly[c]) < 10 ; w += v, c = *it++)
                    w *= 10;
                if (w == 0) w = 1; /* / means 1, and /- means -1.  */
                const int scalar = w;
                for( ;; ) {
                    const int saved_c = c;
                    const auto saved_it = it;
                    if (c == '*') c = *it++;
                    if (c != 's') {
                        /* A * sign is only valid if followed by one or
                         * several s */
                        c = saved_c;
                        it = saved_it;
                        break;
                    }
                    for( ; c == 's' ; j++, c = *it++) ;
                    if (c == '^') {
                        c = *it++;
                        for (w = 0 ; (v = ugly[c]) < 10 ; w += v, c = *it++)
                            w *= 10;
                        /* we already counted one s as meaning
                         * something, so we must take it off. */
                        j += w - 1;
                    }
                }
                (*this)[j % n] += negative ? -scalar : scalar;
            }
        }
        /* }}} */
        template<typename OutputIt>
        auto format_to(OutputIt out) const -> OutputIt
        {
            ASSERT_ALWAYS(*this != 0);
            for(int i = 0, j = 0 ; i < n ; i++) {
                int c = (*this)[i];
                if (c == 0)
                    continue;
                if (c < -1)
                    out = fmt::format_to(out, "{}", c);
                else if (c == -1)
                    out = fmt::format_to(out, "-");
                else if (c == 1 && j)
                    out = fmt::format_to(out, "+");
                else if (c > 1 && j)
                    out = fmt::format_to(out, "+{}", c);
                else if (c > 1)
                    out = fmt::format_to(out, "{}", c);
                j++;
                if (i == 1)
                    out = fmt::format_to(out, "s");
                else if (i == 2)
                    out = fmt::format_to(out, "ss");
                else if (i > 3)
                    out = fmt::format_to(out, "s^{}", i);
            }
            return out;
        }
    }; /* }}} */
    /* }}} */

    /* {{{ sm block
     *
     * sm info is rather easy.
     * <sm> = : <decimal_integer> (, <decimal_integer>)*
     */

    template<typename parent>
    struct sm_block : public parent {
        std::vector<cxx_mpz> sm;
        template<typename iterator>
            int parse(iterator & it0)
            requires requires { *it0++ == '\0'; }
        { /* {{{ */
            sm.clear();
            /* use the "next_delim" hack from primes_block */
            int c = parent::parse(it0, ':');
            iterator it = it0;

            for(bool first = true ; c != '\n' ; first = false) {
                ASSERT_ALWAYS(c == (first ? ':' : ','));
                c = *it++;
                unsigned long v;
                cxx_mpz w;
                for (w = 0 ; (v = ugly[c]) < 16 ; w += v, c = *it++) w <<= 4;
                sm.push_back(w);
            }
            it0 = it;
            return c;
        } /* }}} */
        template<typename OutputIt>
        auto format_to(OutputIt out) const -> OutputIt
        {
            out = parent::format_to(out);
            /* output something that makes it clear that there's data
             * that was removed */
            return fmt::format_to(out, join(sm, ","));
        }
    };
    /* }}} */

    /* {{{ line block
     * This one is just a trivial addition that copies over the line info
     * for later full inspection
     */
    template<typename parent>
    struct line_block : public parent {
        std::string line;
        template<typename iterator>
            int parse(iterator & it0)
            requires requires { *it0++ == '\0'; }
        {
            line.clear();
            iterator it = it0;
            int c0 = parent::parse(it0);
            for(int c ; (c = *it++) != 0 && c != '\n'; ) {
                line.push_back(static_cast<char>(c));
            }
            return c0;
        }
        template<typename OutputIt>
        auto format_to(OutputIt out) const -> OutputIt
        {
            return fmt::format_to(out, line);
        }

    };
    /* }}} */

} /* namespace cado::relation_building_blocks */


struct prime_type_for_sieve_relations {
    using p_type = p_r_values_t;
    using p_or_h_type = p_type;
    using e_type = int;
    p_type p;
    e_type e;
    /* we're probably not in a very sane situation here. primes_block
     * doesn't (and perhaps cannot) dereference the active_sides array to
     * store really the _side_ and not the _side index_ here. We have a
     * variety of existing relation types, and it's not clear at all
     * that we can accomodate all situations easily.
     *
     * For the moment, we have sort of an implicit assumption that
     * active_sides is [0,1] so that for sieve relations, side_index is
     * the same as side.
     */
    int side;
    p_or_h_type& p_or_h() { return p; }
    p_or_h_type const & p_or_h() const { return p; }
    auto operator<=>(prime_type_for_sieve_relations const & o) const {
        if (auto r = side <=> o.side ; r != 0)
            return r;
        if (auto r = p <=> o.p ; r != 0)
            return r;
        return e <=> o.e;
    }
    auto operator==(prime_type_for_sieve_relations const & o) const {
        return (side == o.side) && (p == o.p) && (e == o.e);
    }
    template<typename OutputIt>
    auto format_to(OutputIt out) const -> OutputIt
    {
        static_assert(std::is_same_v<e_type, int>);
        ASSERT_ALWAYS(e >= 0);
        if (e == 0)
            return out;
        out = fmt::format_to(out, "{}", p);
        for(e_type f = 1 ; f < e ; f++)
            out = fmt::format_to(out, ",{}", p);
        return out;
    }
};

struct prime_type_for_indexed_relations {
    using index_type = index_t;
    using p_or_h_type = index_type;
    using e_type = int;
    index_type h;
    e_type e;
    // int side;
    p_or_h_type& p_or_h() { return h; }
    p_or_h_type const & p_or_h() const { return h; }
    auto operator<=>(prime_type_for_indexed_relations const & o) const {
        if (auto r = h <=> o.h ; r != 0)
            return r;
        return e <=> o.e;
    }
    auto operator==(prime_type_for_indexed_relations const & o) const {
        return (h == o.h) && (e == o.e);
    }
    template<typename OutputIt>
    auto format_to(OutputIt out) const -> OutputIt
    {
        static_assert(std::is_same_v<e_type, int>);
        /* for backwards compatibility, we write exponents in unary.
         * However, the way forward is to use the post-exponent notation.
         */
        ASSERT_ALWAYS(e != 0);
        if (e > 0) {
            out = fmt::format_to(out, "{:x}", h);
            for(e_type f = 1 ; f < e ; f++)
                out = fmt::format_to(out, ",{:x}", h);
        } else if (e < 0) {
            out = fmt::format_to(out, "-{:x}", h);
            for(e_type f = 1 ; f < -e ; f++)
                out = fmt::format_to(out, ",-{:x}", h);
        }
        return out;
    }
};

template<int n = 6>
struct prime_type_for_tnfs_sieved_relations {
    using p_type = p_r_values_t;
    using p_or_h_type = p_type;
    using e_type = cado::relation_building_blocks::exponent_vector<n>;
    p_type p;
    e_type e;
    int side;
    p_or_h_type& p_or_h() { return p; }
    p_or_h_type const & p_or_h() const { return p; }
    auto operator<=>(prime_type_for_tnfs_sieved_relations<n> const & o) const {
        if (auto r = side <=> o.side ; r != 0)
            return r;
        if (auto r = p <=> o.p ; r != 0)
            return r;
        return e <=> o.e;
    }
    auto operator==(prime_type_for_tnfs_sieved_relations<n> const & o) const {
        return (side == o.side) && (p == o.p) && (e == o.e);
    }
    template<typename OutputIt>
    auto format_to(OutputIt out) const -> OutputIt
    {
        ASSERT_ALWAYS(e != 0);
        out = fmt::format_to(out, "{:x}", p);
        if (e != 1)
            out = fmt::format_to(out, "/{}", e);
        return out;
    }
};

template<int n = 6>
struct prime_type_for_tnfs_indexed_relations {
    using index_type = index_t;
    using p_or_h_type = index_type;
    using e_type = cado::relation_building_blocks::exponent_vector<n>;
    index_type h;
    e_type e;
    // int side;
    p_or_h_type& p_or_h() { return h; }
    p_or_h_type const & p_or_h() const { return h; }
    auto operator<=>(prime_type_for_tnfs_indexed_relations<n> const & o) const {
        if (auto r = h <=> o.h ; r != 0)
            return r;
        return e <=> o.e;
    }
    auto operator==(prime_type_for_tnfs_indexed_relations<n> const & o) const {
        return (h == o.h) && (e == o.e);
    }
    template<typename OutputIt>
    auto format_to(OutputIt out) const -> OutputIt
    {
        ASSERT_ALWAYS(e != 0);
        out = fmt::format_to(out, "{:x}", h);
        if (e != 1)
            out = fmt::format_to(out, "/{}", e);
        return out;
    }
};

using tnfs_relation =
    cado::relation_building_blocks::primes_block<
        prime_type_for_tnfs_sieved_relations<12>,
        cado::relation_building_blocks::ab_block<uint64_t, 10>
        >;

using tnfs_indexed_relation =
    cado::relation_building_blocks::primes_block<
        prime_type_for_tnfs_indexed_relations<12>,
        cado::relation_building_blocks::ab_block<uint64_t, 16>
        >;



namespace cado::filter_io_details {

    /* {{{ concept is_locking_layer_v */
    template<typename T>
    concept is_locking_layer_v = std::is_empty_v<T>
        && requires { typename T::template critical_datatype<int>; }
        && requires { typename T::lock_t; }
        && requires { typename T::cond_t; }
        && requires { T::lock; }
        && requires { T::unlock; }
        && requires { T::wait; }
        && requires { T::signal; }
        && requires { T::signal_broadcast; }
    ; /* }}} */

    struct ifb_locking_posix {/*{{{*/
        static const int max_supported_concurrent = INT_MAX;
        template<typename T> struct critical_datatype {
            class t {
                T x;
                public:
                T load() const { return x; }
                void store(T a) { x = a; }
                explicit t(T const& a) : x(a) {}
                explicit t() : x(0) {}
                t& operator=(T const& a) { x = a; return *this; }
                T increment() { return x++; }
            };
        };
        using lock_t = std::mutex ;
        using cond_t = std::condition_variable  ;
        static void lock(lock_t * m) { m->lock(); }
        static void unlock(lock_t * m) { m->unlock(); }
        static void wait(cond_t * c, lock_t * m) {
            std::unique_lock<std::mutex> foo(*m, std::adopt_lock);
            c->wait(foo);
            foo.release();
        }
        static void signal(cond_t * c) { c->notify_one(); }
        static void signal_broadcast(cond_t * c) { c->notify_all(); }
        static int isposix() { return 1; }
    };
    static_assert(cado::filter_io_details::is_locking_layer_v<ifb_locking_posix>);
    /*}}}*/

    struct ifb_locking_lightweight {/*{{{*/
        /* we don't support several threads wanting to write to the same
         * location (we could, if we were relying on atomic compare and
         * swap, for instance) */
        static const int max_supported_concurrent = 1;
        template<typename T> struct critical_datatype {
            /* See bug #30068
             *
             * The total store ordering on x86 implies that we can play
             * very dirty games with the completed[] and scheduled[]
             * arrays. The underlying assumptions need not be true in
             * general, and definitely do not hold on arm64.
             *
             * Ideally, there would be a way to qualify our operations on
             * the atomic type that resolve to no emitted code at all if
             * the hardware memory model is x86. But I can't find a way
             * to do that.
             */
#if !defined(__x86_64) && !defined(__i386)
            class t : private std::atomic<T> {
                using super = std::atomic<T>;
                public:
                T load() const { return super::load(std::memory_order_acquire); }
                explicit t(T const& a) : super(a) {}
                explicit t() : super(0) {}
                void store(T a) { super::store(a, std::memory_order_release); }
                t& operator=(T const& a) { store(a); return *this; }
                T increment() { return super::fetch_add(1, std::memory_order_acq_rel); }
            };
#else
            class t {
                volatile T x;
                public:
                T load() const { return x; }
                explicit t(T const& a) : x(a) {}
                explicit t() : x(0) {}
                void store(T a) { x = a; }
                t& operator=(T const& a) { store(a); return *this; }
                /* c++20 frowns upon volatile. Well, it kinda forces us to
                 * use atomics, in fact. Let's silence the issue for the
                 * moment. It may well be that the correct way to go is the
                 * code branch above. But I still long for the optimal way to
                 * write things so that there's a 1-to-1 correspondence with
                 * the simple and easy code here.
                 */
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-volatile"
#elif defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvolatile"
#endif
                T increment() { return x++; }
#ifdef __clang__
#pragma clang diagnostic pop
#elif defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
            };
#endif
        };
        using lock_t = int;
        using cond_t = int;
        template<typename T> static T next(T a, int) { return a; }
        static void lock(lock_t *) {}
        static void unlock(lock_t *) {}
        static void wait(cond_t *, lock_t *) {
            /* {{{ define NANOSLEEP */
            /* The realistic minimal non-CPU waiting with nanosleep is about
             * 10 to 40 microseconds (1<<13 for nanosleep).  But all the I/O
             * between the threads have been buffered, and a thread does a
             * nanosleep only if its buffer is empty.  So I use here ~2ms
             * (1<<21) to optimize CPU scheduler.  Max pause is about 4 to
             * 8ms (1<<22, 1<<23); above that, the program is slowed down.
             */
#ifndef HAVE_NANOSLEEP
#ifdef HAVE_USLEEP
            usleep((unsigned long) (1<<21 / 1000UL));
#else
            sleep(0);
#endif
#else
            struct timespec wait_classical = { .tv_sec=0, .tv_nsec=1<<21 };
            nanosleep(&wait_classical, nullptr);
#endif
            /*}}}*/
        }
        static void signal(cond_t *) {}
        static void signal_broadcast(cond_t *) {}
        static int isposix() { return 0; }
    };
    static_assert(cado::filter_io_details::is_locking_layer_v<ifb_locking_lightweight>);
    /*}}}*/

    /* {{{ status table (utility for inflight_rels_buffer).
     *
     * In fact, when we use simple busy waits, we are restricted to
     * scheduled[k]==completed[k]+(0 or 1), and keeping track of the
     * processing level is useless. So we provide a trimmed-down
     * specialization for this case.
     *
     * the status table depends on the maximum number of threads per step. We
     * make it depend on the locking backend instead, for simplicity. */
    template<typename locking>
    struct status_table {
        using csize_t = typename locking::template critical_datatype<size_t>::t;
        typename locking::template critical_datatype<int8_t>::t x[SIZE_BUF_REL];
        /* {{{ ::catchup() (for ::schedule() termination) */
        void catchup(csize_t & last_completed, size_t last_scheduled, int level) {
            size_t c = last_completed.load();
            for( ; c < last_scheduled ; c++) {
                if (x[c & (SIZE_BUF_REL-1)].load() < level)
                    break;
            }
            last_completed.store(c);
        }
        /*}}}*/
        /*{{{ ::catchup_until_mine_completed() (for ::complete()) */
        /* (me) is the absolute relation index of the relation I'm currently
         * processing (the value of schedule[k] when it was called prior to
         * giving me this relation to process).
         */
        void catchup_until_mine_completed(csize_t & last_completed, size_t me, int level) {
            const size_t slot = me & (SIZE_BUF_REL-1);
            size_t c = last_completed.load();
            ASSERT(x[slot].load() == (int8_t) (level-1));
            /* The big question is how far we should go. By not exactly answering
             * this question, we avoid the reading of scheduled[k], which is good
             * because it is protected by m[k-1]. And even if we could consider
             * doing a rwlock for reading this, it's too much burden. So we leave
             * open the possibility that many relation slots ahead of us already
             * have x[slot] set to k, yet we do not increment
             * completed[k] that far. This will be caught later on by further
             * processing at this level.
             *
             * This logic is problematic regarding termination, though. See the
             * termination code in ::complete()
             */
            for( ; c < me ; c++) {
                if (x[c & (SIZE_BUF_REL-1)].load() < level)
                    break;
            }
            last_completed.store(c + (c == me));
            x[slot].increment();
            ASSERT(x[slot].load() == (int8_t) (level));
        }
        /*}}}*/
        void update_shouldbealreadyok(size_t slot, int level) {
            if (level < 0) {
                x[slot & (SIZE_BUF_REL-1)].store(level);
            } else {
                ASSERT(x[slot & (SIZE_BUF_REL-1)].load() == level);
            }
        }
    };

    template<>
    struct status_table<ifb_locking_lightweight> {
        using csize_t = ifb_locking_lightweight::critical_datatype<size_t>::t;
        static void catchup(csize_t & last_completed, size_t last_scheduled, int) {
            ASSERT_ALWAYS(last_completed.load() == last_scheduled);
        }
        static void catchup_until_mine_completed(csize_t & last_completed, size_t, int) {
            last_completed.increment();
        }
        void update_shouldbealreadyok(size_t, int) {}
    };
    /* }}} */

    /* {{{ multithreaded_call: pass this instead of a bare lambda in
     * order to have filter_rels spawn multiple threads at that
     * processing level while parsing the input files.
     */
    template<typename F>
        requires (!requires { typename F::first_type; })
    struct multithreaded_call : std::pair<size_t, F> {
        size_t nthreads() const { return this->first; }
        multithreaded_call(F const & f)
            : std::pair<size_t, F> { 1, f }
        {}
        multithreaded_call(size_t k, F const & f)
            : std::pair<size_t, F> { k, f }
        {}
        /*
        multithreaded_call(multithreaded_call<F> const & fk)
            : std::pair<size_t, F> { fk.first, fk.second }
        {}
        multithreaded_call & operator=(multithreaded_call<F> const & fk)
        {
            std::pair<size_t, F>::first = fk.first;
            std::pair<size_t, F>::second = fk.second;
            return *this;
        }
        */
        template<typename... Args>
        void operator()(Args&& ...args) const { std::pair<size_t, F>::second(std::forward<Args>(args)...); }
    };
    /* }}} */

    template<typename X>
    struct number_of_threads_helper {
        size_t operator()(X const &) const { return 1; }
    };
    template<typename F>
    struct number_of_threads_helper<multithreaded_call<F>> {
        size_t operator()(multithreaded_call<F> const & m) const
        { return m.nthreads(); }
    };
    template<typename X>
    size_t number_of_threads(X const & x)
    {
        return number_of_threads_helper<X>()(x);
    }

    /* {{{ inflight_rels_buffer: n-level buffer, with underyling locking
     * mechanism specified by the template class.  */
    template<
        typename locking_layer,
        typename relation_type,
        /* n is the "depth" of the processing pipe. We have n-1 functions
         * in total, which move relation from one processing level to the
         * next. Thus n distinct processing levels exist, and we have n-1
         * possible transitions.  */
        size_t n>
    struct inflight_rels_buffer {
        std::barrier<cado::nop_function> sync_point;
        using csize_t = locking_layer::template critical_datatype<size_t>::t;
        using lock_t = locking_layer::lock_t;
        using cond_t = locking_layer::cond_t;

        std::unique_ptr<relation_type[]> rels;
        /* always malloc()-ed to SIZE_BUF_REL, which is a power of two */
        /* invariant:
         * scheduled_0 >= ... >= completed_{n-1} >= scheduled_0 - SIZE_BUF_REL
         */
        csize_t completed[n];
        csize_t scheduled[n];
        status_table<locking_layer> status;
        lock_t m[n];
        cond_t bored[n];
        int active[n] = { 0, } ;     /* number of active threads */

        /*{{{ ctor */
        explicit inflight_rels_buffer(int nthreads_total)
            : sync_point(nthreads_total)
            , rels(std::make_unique<relation_type[]>(SIZE_BUF_REL))
        {
            // fill_n doesn't work with atomics, at least not with apple
            // clang 17.
            // std::fill_n(completed, n, csize_t(0));
            // std::fill_n(scheduled, n, csize_t(0));
            for(size_t i = 0 ; i < n ; i++) {
                completed[i] = 0;
                scheduled[i] = 0;
            }
        }/*}}}*/

        /* compared to the "original" inflight_rels_buffer type, this version
         * assumes that the relation_type has well-defined ctors and ctors,
         * so that everything can be default-destructed.
         */
        ~inflight_rels_buffer() = default;
        inflight_rels_buffer(inflight_rels_buffer const&) = delete;
        inflight_rels_buffer(inflight_rels_buffer &&) = delete;
        inflight_rels_buffer& operator=(inflight_rels_buffer const&) = delete;
        inflight_rels_buffer& operator=(inflight_rels_buffer &&) = delete;

        void drain() /*{{{ */
            /* This belongs to the buffer closing process.  The out condition of
             * this call is that all X(k) for k>0 terminate.  This call (as well
             * as init/clear) must be called on the producer side (step 0) (in a
             * multi-producer context, only one thread is entitled to call this)
             */
        {
            // size_t c = completed[0];
            active[0]--;

            for(size_t k = 0 ; k < n ; k++) {
                locking_layer::lock(m + k);
                while(active[k]) {
                    locking_layer::wait(bored + k, m + k);
                }
                completed[k].store(SIZE_MAX);
                locking_layer::signal_broadcast(bored + k);
                locking_layer::unlock(m + k);
            }
        }
        /*}}}*/

        relation_type * schedule(size_t k) /* {{{ */
            /* Schedule a new relation slot for processing at level k.
             *
             * This call may block until a relation is processed by level k-1
             * (or, if k==0, until a slot is made available in the relation
             * buffer).
             *
             * The relation is free for use by the current (consumer) thread
             * until it calls inflight_rels_buffer. When the producing
             * stream ends, this function returns nullptr. */
        {
            size_t const prev = k ? (k-1) : (n-1);
            // coverity[result_independent_of_operands]
            ASSERT(active[k] <= locking_layer::max_supported_concurrent);
            size_t s;
            size_t const a = k ? 0 : SIZE_BUF_REL;
            /* in 1-thread scenario, scheduled[k] == completed[k] */
            locking_layer::lock(m + prev);
            if constexpr (locking_layer::max_supported_concurrent == 1) {
                /* can't change */
                s = scheduled[k].load();
                while(s == a + completed[prev].load()) {
                    locking_layer::wait(bored + prev, m + prev);
                }
            } else {
                while((s=scheduled[k].load()) == a + completed[prev].load()) {
                    locking_layer::wait(bored + prev, m + prev);
                }
            }
            /* when completed[prev] == SIZE_MAX, the previous-level workers
             * are creating spuriouss relation created to trigger termination.
             * In this case, scheduled[prev] is safe to read now. we use it
             * as a marker to tell whether there's still work ahead of us, or
             * not.  */
            if (UNLIKELY(completed[prev].load() == SIZE_MAX) && scheduled[prev].load() == s) {
                /* prepare to return */
                /* note that scheduled[k] is *not* bumped here */
                locking_layer::unlock(m + prev);
                /* we emulate the equivalent of ::complete(), and terminate */
                locking_layer::lock(m + k);
                status.catchup(completed[k], s, k);
                active[k]--;
                locking_layer::signal_broadcast(bored + k);
                locking_layer::unlock(m + k);
                return nullptr;
            }
            // ASSERT(scheduled[k] < a + completed[prev]);
            scheduled[k].increment();
            const size_t slot = s & (SIZE_BUF_REL - 1);
            relation_type * rel = &rels[slot];
            status.update_shouldbealreadyok(s, k-1);
            locking_layer::unlock(m + prev);
            return rel;
        }
        /*}}}*/

        void complete(int k, relation_type const * rel) /* {{{ */
        {
            // coverity[result_independent_of_operands]
            ASSERT(active[k] <= locking_layer::max_supported_concurrent);
            const int slot = rel - rels.get();

            locking_layer::lock(m + k);

            size_t my_absolute_index;
            if constexpr (locking_layer::max_supported_concurrent == 1) {
                my_absolute_index = completed[k].load();
            } else {
                /* recover the integer relation number being currently
                 * processed from the one modulo SIZE_BUF_REL.
                 *
                 * We have (using ck = completed[k]):
                 *          ck <= zs < ck + N
                 *          ck <= s+xN < ck + N <= s+(x+1)N
                 *          xN < ck-s + N <= (x+1) N
                 *
                 */
                const size_t c = completed[k].load();
                my_absolute_index = slot;
                my_absolute_index += ((c - slot + SIZE_BUF_REL - 1) & -SIZE_BUF_REL);
            }

            /* morally, this is completed[k]++ */
            status.catchup_until_mine_completed(completed[k], my_absolute_index, k);
            locking_layer::signal_broadcast(bored + k);
            locking_layer::unlock(m + k);
        }
        /*}}}*/

        /* computation threads joining the computation are calling these */
        void enter(int k) {
            locking_layer::lock(m+k);
            active[k]++;
            locking_layer::unlock(m+k);
            sync_point.arrive_and_wait();
        }
        /* leave() is a no-op, since active-- is performed as part of the
         * normal drain() call */
        void leave(int) { }

        /* The calling scenario is as follows.
         *
         * For the owner thread.
         *  - constructor
         *  - start workers.
         *  - enter(0)
         *  - some schedule(0) / complete(0) for relations which get fed in.
         *  - drain() once all are produced
         *  - leave(0)
         *
         * For the workers (there may be more at each level if
         * ifb_locking_posix is used):
         *  - enter(k)
         *  - a loop on with schedule(k) / complete(k), exiting when
         *    schedule() returns NULL.
         *  - leave(k)
         *
         * Currently the owner thread is weakly assumed to be the only
         * level-0 thread, but that does not seem to be an absolute necessity
         * from the design. Additional level-0 threads would induce a loop
         * similar to other worker threads, but the fine points haven't been
         * considered yet.
         *
         * The current implementation has leave() a no-op, and uses drain()
         * at the owner thread to to a shutdown. This could change.
         */

        /*{{{ filter_rels consumer thread */
        template<size_t level>
            void spawn_threads(std::vector<std::thread> &,
                    timingstats_dict_ptr)
            {}

        template<size_t level, typename F, typename... Functions>
            void spawn_threads(std::vector<std::thread> & pool,
                    timingstats_dict_ptr stats,
                    F const & fk, Functions const & ...ff)
            {
                constexpr size_t k = level + 1;
                for(auto t = number_of_threads(fk); t-- ; ) {
                    pool.emplace_back([&, stats]() {
                            enter(k);
                            relation_type * slot;
                            for( ; (slot = schedule(k)) != nullptr ; ) {
                                fk(*slot);
                                complete(k, slot);
                            }
                            leave(k);
                            if (stats)
                                timingstats_dict_add_mythread(stats, "consumer");
                            });
                }
                spawn_threads<level + 1, Functions...>(pool, stats, ff...);
                /*
                   double thread_times[2];
                   thread_seconds_user_sys(thread_times);
                   fprintf(stderr, "Consumer thread (level %d) ends after having spent %.2fs+%.2fs on cpu\n", k, thread_times[0], thread_times[1]);
                   */
                // if (stats) timingstats_dict_add_mythread(stats, "consumer");
            }

/*}}}*/
    };
    /* }}} */

    /* {{{ all_functions_are_invocable: make sure that we can the
     * relation type and the functions in the template parameter pack are
     * all compatible.
     */
    template<typename relation_type, typename... Functions>
    struct all_functions_are_invocable
        : public std::false_type
    {};

    /* "all" functions in an empty pack always meet the constraint :-) */
    template<typename relation_type>
    struct all_functions_are_invocable<relation_type>
        : public std::true_type
    {};

    /* multithreaded function reduce to the single-thread condition */
    template<typename relation_type, typename F, typename... Functions>
    struct all_functions_are_invocable<
                relation_type,
                multithreaded_call<F>,
                Functions...>
        : all_functions_are_invocable<relation_type, F, Functions...>
    {};

    template<typename relation_type, typename F, typename... Functions>
                    requires std::is_invocable_v<F, relation_type &>
    struct all_functions_are_invocable<relation_type, F, Functions...>
        : all_functions_are_invocable<relation_type, Functions...>
    {};

    template<typename relation_type, typename... Functions>
    inline constexpr bool all_functions_are_invocable_v = all_functions_are_invocable<relation_type, Functions...>::value;
    /* }}} */

    /* This is always the same, it can live in the other compilation unit.  */
    extern void filter_rels_producer_thread(
            ringbuf & r,
            std::vector<std::string> const & input_files,
            timingstats_dict_ptr stats);
    extern void filter_rels_producer_thread(
            ringbuf & r,
            std::istream& is,
            timingstats_dict_ptr stats);

    template<typename iterator_type, typename relation_type>
    void parse_helper(relation_type & rel, uint64_t num, iterator_type it)
    requires requires { rel.parse(it); }
    {
        rel.parse(it);
        rel.num = num;
    }
    template<typename iterator_type, typename ... variants>
    void parse_helper(std::variant<variants...> & rel, uint64_t num, iterator_type it)
    {
        rel = {};
        auto & rel0 = std::get<0>(rel);
        rel0.parse(it);
        rel0.num = num;
    }

} /* namespace cado::filter_io_details */


template<
    typename locking_t,
    typename rel_t>
struct filter_rels_obj {
    using locking_layer = locking_t;
    using relation_type = rel_t;

    /* as annoying as it may be, we need this intermediary template for
     * the std::thread to compile correctly.
     */
    template<typename InputDescriptionType>
    static void producer(
            ringbuf & r,
            InputDescriptionType input,
            timingstats_dict_ptr stats)
    {
        cado::filter_io_details::filter_rels_producer_thread(r, input, stats);
    }

    template<typename... Args>
    size_t operator()(
            std::vector<std::string> const & input_files,
            Args&& ...args)
    {
        if (input_files.empty())
            return 0;

        return backend_call(input_files, std::forward<Args>(args)...);
    }

    template<typename... Args>
    size_t operator()(
            std::string const & f,
            Args&& ...args)
    {
        const std::vector<std::string> fs = {f};
        return backend_call(fs, std::forward<Args>(args)...);
    }

    template<typename... Args>
    size_t operator()(
            std::istream & is,
            Args&& ...args)
    {
        return backend_call(is, std::forward<Args>(args)...);
    }

    template<typename InputDescriptionType, typename... Functions>
        requires (cado::filter_io_details::all_functions_are_invocable_v<relation_type, Functions...> && cado::filter_io_details::is_locking_layer_v<locking_layer>)
    size_t backend_call(
            InputDescriptionType input_files,
            std::vector<bool> const * active,
            timingstats_dict_ptr tstats,
            Functions&& ...FF) const
    {
        constexpr size_t N = sizeof...(FF);
        using cado::filter_io_details::inflight_rels_buffer;
        using inflight_t = inflight_rels_buffer<locking_layer, relation_type, N + 1>;

        stats_data_t infostats;  /* for displaying progress */
        uint64_t nrels = 0, nactive = 0;
        size_t nB = 0;

        ringbuf rb(0);

        std::thread P(
                producer<InputDescriptionType>,
                std::ref(rb),
                std::forward<InputDescriptionType>(input_files),
                tstats);

        /* variants that we know about:
         *
         * _(AB_DECIMAL)|_(LINE): // dup1
         * _(AB_HEXA)|_(LINE): // dup1 + -abhexa
         * _(AB_HEXA): // dup2 (for renumbered files)
         * _(AB_DECIMAL)|_(PRIMES): // dup2/pass2
         * _(AB_HEXA)|_(PRIMES): // dup2/pass2
         * _(INDEX)|_(SORTED): // ???
         * _(INDEX): // all binaries after dup2 that do not need a,b
         * _(INDEX) | _(AB_HEXA): // reconstructlog
         * _(INDEX) | _(AB_HEXA) | _(SM): // reconstructlog
         * _(LINE): // purge/2 
         * _(LINE) | _(INDEX): // ???
         */

        // using cado::filter_io_details::multithreaded_call;

        // size_t number_of_consumers = (multithreaded_call(FF).nthreads() + ... + 0);

        using cado::filter_io_details::number_of_threads;
        size_t number_of_consumers = (number_of_threads(FF) + ... + 0);

        /* now prepare the inflight buffer, and the appropriate threads */
        inflight_t inflight(number_of_consumers + 1);

        /* {{{ setup and start all the consumer threads (at all levels) */
        /* these first few linse are also found in the non-templated function
         * which instantiates and calls us.
         */

        std::vector<std::thread> consumers;
        consumers.reserve(number_of_consumers);
        inflight.template spawn_threads<0, Functions...>(
                consumers, tstats, std::forward<Functions>(FF)...);

        /* }}} */

        /* {{{ main loop */

        /* will print report at 2^10, 2^11, ... 2^23 read rels and every 2^23 rels
         * after that */
        if (!active)
            stats_init (infostats, stdout, &nrels, 23, "Read", "relations", "",
                    "rels");
        else
            stats_init (infostats, stdout, &nrels, 23, "Read", "active rels",
                    "read rels", "rels");

        inflight.enter(0);
        for(size_t avail_seen = 0 ; ; ) {
            {
                std::unique_lock ux(rb.mx);
                while(rb.avail_to_read == avail_seen && !rb.done)
                    rb.bored.wait(ux);
                avail_seen = rb.avail_to_read; /* must be before mutex unlock ! */
                if (avail_seen == 0 && rb.done) {
                    /* end of producer1 is with rb->done = 1 -- which is
                     * compatible with bytes still being in the pipe ! */
                    break;
                }
            }

            /* We may have one or several lines which have just been
             * produced. As long as we succeed reading complete lines, we
             * consume them, and feed the second pipe.
             */
            int nl;
            for(size_t avail_offset = 0; avail_offset < avail_seen && (nl = rb.strchr('\n', 0)) >= 0 ; ) {
                /* skip comments and blank lines */
                if (*rb.begin() != '#' && *rb.begin() != '\n') {
                    uint64_t const relnum = nrels++;
                    if (!active || (*active)[relnum]) {
                        auto * slot = inflight.schedule(0);
                        ASSERT_ALWAYS(slot);
                        cado::filter_io_details::parse_helper(*slot, relnum, rb.begin());
                        inflight.complete(0, slot);
                        nactive++;
                    }
                }
                /* skip the newline byte as well */
                nl++;
                nB += nl;
                rb.skip_get(nl);
                avail_seen -= nl;
                avail_offset += nl;
                if (stats_test_progress(infostats))
                {
                    if (!active)
                        stats_print_progress (infostats, nrels, 0, nB, 0);
                    else
                        stats_print_progress (infostats, nactive, nrels, nB, 0);
                }
            }
        }
        inflight.drain();
        inflight.leave(0);
        /*}}}*/

        /* {{{ join all threads */
        for(auto & t : consumers)
            t.join();
        P.join();
        /*}}}*/

        /* NOTE: the inflight dtor is called automatically */

        if (!active)
            stats_print_progress (infostats, nrels, 0, nB, 1);
        else
            stats_print_progress (infostats, nactive, nrels, nB, 1);

        return nactive;
    }
};

template<
    typename locking_type,
    typename relation_type,
    typename InputDescriptionType,
    typename... Functions>
size_t filter_rels(
        InputDescriptionType input_description,
        std::vector<bool> const * active,
        timingstats_dict_ptr tstats,
        Functions&& ...FF)
{
    return filter_rels_obj<locking_type, relation_type>()(
            std::forward<InputDescriptionType>(input_description),
            active, tstats,
            std::forward<Functions>(FF)...);
}

#endif /* CADO_FILTER_IO_HPP */
