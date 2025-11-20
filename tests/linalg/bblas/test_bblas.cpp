#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <string>
#include <set>
#include <vector>

#include "params.h"

#include "test_bblas_base.hpp"
#include "test_bblas_level2.hpp"
#include "test_bblas_level3.hpp"
#include "test_bblas_level4.hpp"
#include "test_bblas_level5.hpp"
#include "test_bpack.hpp"

void
static print_features() /*{{{*/
{
    printf("## compile-time features\n");
#ifdef HAVE_M4RI
    printf("## HAVE_M4RI\n");
#endif /* HAVE_M4RI */
#ifdef HAVE_M4RIE
    printf("## HAVE_M4RIE\n");
#endif /* HAVE_M4RIE */
#ifdef HAVE_PCLMUL
    printf("## HAVE_PCLMUL\n");
#endif /* HAVE_PCLMUL */
#ifdef HAVE_SSE2
    printf("## HAVE_SSE2\n");
#endif /* HAVE_SSE2 */
#ifdef HAVE_SSE41
    printf("## HAVE_SSE41\n");
#endif /* HAVE_SSE41 */
#ifdef HAVE_AVX2
    printf("## HAVE_AVX2\n");
#endif /* HAVE_AVX2 */
#ifdef VALGRIND
    printf("## VALGRIND\n");
#endif /* VALGRIND */
    printf("## ULONG_BITS=%d\n", ULONG_BITS);
} /*}}}*/

static void
declare_usage(cxx_param_list& pl) /*{{{*/
{
    param_list_decl_usage(pl, "seed", "seed for random data generation");
    param_list_decl_usage(pl, "n", "n value for some size-dependent tests");
    param_list_decl_usage(pl, "fast", "do quick tests only\n");
    param_list_decl_usage(pl, "tests", "list of tests to perform\n");
} /*}}}*/

// coverity[root_function]
int main(int argc, char const * argv[])
{
    cxx_param_list pl;
    unsigned int n = 2 * 1000 * 1000;
    int seed = 0;
    declare_usage(pl);
    param_list_configure_switch(pl, "-fast", &test_bblas_base::test_accel);
    const char* argv0 = argv[0];
    argc--, argv++;
    for (; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    param_list_parse_uint(pl, "n", &n);

    std::vector<std::string> tests;
    if (!param_list_parse(pl, "tests", tests))
        tests.emplace_back("all");

    param_list_parse(pl, "seed", seed);

    if (!seed)
        seed = time(nullptr);

    setvbuf(stdout, nullptr, _IONBF, 0);
    setvbuf(stderr, nullptr, _IONBF, 0);

    print_features();

    printf("# seeding with seed %d\n", seed);

    std::set<std::string> seen;

    test_bblas_level2 A2(n);
    A2.set_seed(seed);
    A2(tests, seen);

    test_bblas_level3 A3(n);
    A3.set_seed(seed);
    A3(tests, seen);

    test_bblas_level4 A4(n);
    A4.set_seed(seed);
    A4(tests, seen);

    test_bblas_level5 A5(n);
    A5.set_seed(seed);
    A5(tests, seen);

    test_bpack A6(n);
    A6.set_seed(seed);
    A6(tests, seen);

    for (auto const& s : tests) {
        if (seen.find(s) == seen.end())
            fprintf(stderr, "## no test for key %s\n", s.c_str());
    }

    return 0;
}
