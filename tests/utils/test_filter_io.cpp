#include "cado.h"       // IWYU pragma: keep

#include <concepts>
#include <exception>
#include <string>
#include <unordered_map>

#include "fmt/base.h"

#include "filter_io.hpp"
#include "utils_cxx.hpp"

using rel_dup1 = cado::relation_building_blocks::primes_block<
    prime_type_for_sieve_relations,
    cado::relation_building_blocks::ab_block<uint64_t, 16>>;

using rel_dup2 = cado::relation_building_blocks::primes_block<
    prime_type_for_indexed_relations,
    cado::relation_building_blocks::ab_block<uint64_t, 16>>;

using rel_purge = cado::relation_building_blocks::primes_block<
    prime_type_for_indexed_relations,
    cado::relation_building_blocks::ab_ignore<16>>;

using rel_reconstructlog =
        cado::relation_building_blocks::sm_block<
        cado::relation_building_blocks::primes_block<
            prime_type_for_indexed_relations,
        cado::relation_building_blocks::ab_block<uint64_t, 16>>>;


template<typename relation_type>
struct TestSuite
{
    struct Valid {};

    using e_type = relation_type::e_type;
    using p_or_h_type = relation_type::p_or_h_type;

    void _unexpected_exception_msg(
            std::string const & s, std::exception const & e)
    {
        fmt::print(stderr, "Error, relation '{}' threw unexpected exception of"
                           "type '{}' with message '{}'\n", s, typeid(e).name(),
                           e.what());
    }

    void _wrong_ret_msg(
            std::string const & s, int expected, int got)
    {
        fmt::print(stderr, "Error, relation '{}' was not parsed correctly, "
                           "expected parse to return '0x{:x}', got '0x{:x}'\n",
                           s, expected, got);
    }

    static
    auto gen_pred_exponents(std::unordered_map<p_or_h_type, e_type> const & V)
    {
        auto pred = [&V](std::string const & s, relation_type const & rel) {
            for(auto const & pe : rel.primes) {
                if (!V.contains(pe.p_or_h())) {
                    fmt::print(stderr, "Error, relation '{}' was not parsed "
                                       "correctly: prime/index 0x{:x} was not "
                                       "expected but appears with exponent "
                                       "{}\n", s, pe.p_or_h(), pe.e);
                    return false;
                }
                auto e = V.at(pe.p_or_h());
                if (pe.e != e) {
                    fmt::print(stderr, "Error, relation '{}' was not parsed "
                                       "correctly: got {} for exponent of "
                                       "0x{:x}, {} was expected\n",
                                       s, pe.e, pe.p_or_h(), e);
                    return false;
                }
            }
            return true;
        };
        return pred;
    };

    template<std::predicate<std::string const &, relation_type> P>
    bool test(std::string const s, P pred, int expected_ret = '\0')
    {
        bool res = true;
        try {
            auto it = s.begin();
            relation_type rel;
            int c = rel.parse(it);

            if (c != expected_ret) {
                _wrong_ret_msg(s, expected_ret, c);
                res = false;
            } else if (not pred(s, rel)) {
                /* failing predicate should print its own error message */
                res = false;
            }
        } catch (std::exception const & e) {
            _unexpected_exception_msg(s, e);
            res = false;
        }
        ok = ok and res;
        return res;
    };

    template<std::same_as<Valid>>
    bool test(std::string const s, int expected_ret = '\0')
    {
        auto p = [](std::string const &, relation_type const &) {return true;};
        return test(s, p, expected_ret);
    };

    template<std::derived_from<std::exception> E>
    bool test(std::string const s, int expected_ret = '\0')
    {
        bool res = true;
        try {
            auto it = s.begin();
            relation_type rel;
            int c = rel.parse(it);

            res = false;
            if (c != expected_ret)
                _wrong_ret_msg(s, expected_ret, c);
            else
                fmt::print(stderr, "Error, relation '{}' should throw\n", s);
        } catch (E const & e) {
            /* parsing is supposed to throw exception of type E */
        } catch (std::exception const & e) {
            _unexpected_exception_msg(s, e);
            res = false;
        }
        ok = ok and res;
        return res;
    };

    bool ok = true;
};


bool
test_too_large_input()
{
    TestSuite<rel_dup2> T;
    using Valid = decltype(T)::Valid;
    using cado_out_of_range = cado::filter_io::out_of_range;
    using cado_cast_failure = runtime_numeric_cast_details::failure;

    /* Test too large a and b; two possible errors: parsing or casting */
    T.test<Valid>            ("e,2a:");
    T.test<cado_out_of_range>("10000000000000000,2a:");
    T.test<cado_out_of_range>("-fedcba9876543210f,2a:");
    T.test<cado_out_of_range>("e,123456789123456789abc:");
    T.test<cado_cast_failure>("ffffffffffffffff,2a:");
    T.test<cado_cast_failure>("8000000000000000,2a:");
    T.test<cado_cast_failure>("-8000000000000001,2a:");
    T.test<Valid>            ("-8000000000000000,2a:");
    T.test<Valid>            ("7fffffffffffffff,2a:");
    T.test<Valid>            ("-7fffffffffffffff,2a:");
    /* Test too large side info */
    T.test<Valid>            ("e,2a@101,51:");
    T.test<cado_out_of_range>("e,2a@18446744073709551616,51:");
    T.test<cado_out_of_range>("e,2a@101,18446744073709551616:");
    /* Test too large indexes */
    T.test<cado_out_of_range>("e,2a:10000000000000000");
    T.test<cado_out_of_range>("e,2a:-fedcba9876543210f"); // deprecated syntax
    T.test<cado_out_of_range>("e,2a:123456789123456789abc/42");
    /* Test too large exponent */
    T.test<Valid>            ("e,2a:12a/-1");
    T.test<cado_out_of_range>("e,2a:12a/18446744073709551616");
    T.test<cado_out_of_range>("e,2a:12a/-18446744073709551616");

    return T.ok;
}

template<typename relation_type>
bool
test_parse_error()
{
    TestSuite<relation_type> T;
    using Valid = typename decltype(T)::Valid;
    using parse_error = cado::filter_io::parse_error;

    T.template test<Valid>      ("e,2a:");
    T.template test<Valid>      ("e,2A:");
    T.template test<parse_error>("e,2a");
    T.template test<parse_error>("e,5g:");
    T.template test<parse_error>("e:");
    /* not valid, according to the spec, but is parsed as 'e,0:' */
    //T.template test<parse_error>("e,:");
    /* not valid, according to the spec, but is parsed as '0,2a:' */
    //T.template test<parse_error>(",2a:");

    T.template test<Valid>      ("e,2a@0,1:");
    T.template test<parse_error>("e,2a@42:");
    /* not valid, according to the spec, but is parsed as 'e,2a@42,0:' */
    //T.template test<parse_error>("e,2a@42,:");
    /* not valid, according to the spec, but is parsed as 'e,2a@0,42:' */
    //T.template test<parse_error>("e,2a@,42:");

    T.template test<Valid>      ("e,2a:12a");
    T.template test<Valid>      ("e,2a:12A");
    /* not valid, according to the spec, but is parsed as 'e,2a:0/' */
    //T.template test<parse_error>("e,2a:/");
    /* not valid, according to the spec, but is parsed as 'e,2a:0,0,12a' */
    //T.template test<parse_error>("e,2a:,,12a");
    T.template test<parse_error>("e,2a:s");
    T.template test<parse_error>("e,2a:z");
    T.template test<parse_error>("e,2a:11;12a");
    T.template test<parse_error>("e,2a:5g");

    T.template test<Valid>      ("e,2a:12a/");
    T.template test<Valid>      ("e,2a:12a/-");
    T.template test<parse_error>("e,2a:12a/052");
    T.template test<parse_error>("e,2a:12a/2a");
    T.template test<parse_error>("e,2a:12a/;11");

    T.template test<Valid>      ("e,2a:12a/");

    return T.ok;
}

bool
test_parse_sm()
{
    /* note: relation with sm must end with '\n' */
    TestSuite<rel_reconstructlog> T;
    using Valid = typename decltype(T)::Valid;
    using parse_error = cado::filter_io::parse_error;

    T.test<Valid>      ("e,2a::12a\n", '\n');
    T.test<Valid>      ("e,2a::12A,11\n", '\n');
    T.test<Valid>      ("e,2a::\n", '\n'); /* parsed as having no SM */
    T.test<Valid>      ("e,2a::0\n", '\n');
    T.test<Valid>      ("e,2a::123456789abcdef123456789ABCDEF\n", '\n');
    T.test<Valid>      ("e,2a::01234\n", '\n'); /* not parsed as octal */
    /* not valid, according to the spec, but is parsed as 'e,2a::0,0\n' */
    //T.test<parse_error>("e,2a::,\n", '\n');
    /* not valid, according to the spec, but is parsed as 'e,2a::0,0\n' */
    //T.test<parse_error>("e,2a::,0\n", '\n');
    T.test<parse_error>("e,2a::5g\n", '\n');
    T.test<parse_error>("e,2a::s\n", '\n');
    T.test<parse_error>("e,2a::z\n", '\n');
    T.test<parse_error>("e,2a::12a;11\n", '\n');
    T.test<parse_error>("e,2a::-12a\n", '\n');

    return T.ok;
}
bool
test_integer_exponents()
{
    TestSuite<rel_dup2> T;

    T.test("e,2a:11/-,12a/-1", T.gen_pred_exponents({{0x11, -1}, {0x12a, -1}}));
    T.test("e,2a:11/,12a/1,101010",  T.gen_pred_exponents({{0x11, 1},
                                                           {0x12a, 1},
                                                           {0x101010, 1}}));
    T.test("e,2a:11/5,12a/2,12a/12", T.gen_pred_exponents({{0x11, 5},
                                                           {0x12a, 14}}));
    T.test("e,2a:11/7,12a/-3,11/-11", T.gen_pred_exponents({{0x11, -4},
                                                            {0x12a, -3}}));
    T.test("e,2a:-11", T.gen_pred_exponents({{0x11, -1}})); // deprecated syntax

    return T.ok;
}

bool
test_tnfs_exponents()
{
    TestSuite<tnfs_indexed_relation> T;

    using PH = decltype(T)::p_or_h_type;
    using Exp = decltype(T)::e_type;

    Exp one {1};                    /* [1, 0, 0, .... ] ~ 1*/
    Exp mone {-1};                  /* [-1, 0, 0, .... ] ~ -1 */
    Exp e8; e8[1] = -1; e8[2] = 2;  /* [0, -1, 2, 0, ....] ~ -s + 2*s^2 */
    Exp e10; e10[4] = 1;            /* [0, 0, 0, 0, 1, 0, ... ] ~ s^4 */

    std::unordered_map<PH, Exp> E { {0x3, one}, {0x4, mone}, {0x5, one},
                                    {0x8, e8}, {0xb, one}, {0x10, e10}};
    auto check = T.gen_pred_exponents(E);

    T.test("11,2a:3,4/-,5,8/-s+2ss,b,10/s^4", check);
    T.test("11,2a:3,4/-1,5,8/-s+2*s*s,b,10/s^2*s^2\n", check, '\n');
    T.test("11,2a:3/,4/-s^0,5,8/-1*s+s*s,8/s*s,b,10/s^2*s^2", check);

    return T.ok;
}

int main()
{
    bool ok = true;

    ok &= test_parse_error<rel_dup1>();
    ok &= test_parse_error<rel_dup2>();
    ok &= test_parse_error<rel_purge>();
    ok &= test_parse_sm();
    ok &= test_integer_exponents();
    ok &= test_tnfs_exponents();
    ok &= test_too_large_input();

    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}

