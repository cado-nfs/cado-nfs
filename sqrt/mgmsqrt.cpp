#include "cado.h" // IWYU pragma: keep

/* see montgomery.md for some documentation. */

#include <cstdlib>
#include <cstdio>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

#include "fmt/base.h"

#include "cado_poly.h"
#include "params.h"
#include "renumber.hpp"
#include "relation.hpp"
#include "typedefs.h"
#include "utils_cxx.hpp"
#include "numbertheory/number_field.hpp"
#include "numbertheory/number_field_order.hpp"

struct montgomery_eth_root {
    number_field const & K;
    number_field_order const & OK;
    std::vector<std::pair<index_t, int>> const & ideals;
    std::vector<std::pair<relation_ab, int>> const & abs;
    int prec;

    montgomery_eth_root(number_field const & K,
            number_field_order const & OK,
            std::vector<std::pair<index_t, int>> const & ideals,
            std::vector<std::pair<relation_ab, int>> const & abs,
            int prec)
        : K(K)
        , OK(OK)
        , ideals(ideals)
        , abs(abs)
        , prec(prec)
    {
    }
};

struct mgmsqrt {
    cxx_cado_poly cpoly;
    unsigned int prec = 1000;
    renumber_t R;
    int e;
    int side;

    number_field K;
    number_field_order OK;

    std::vector<std::pair<index_t, int>> ideals;
    std::vector<std::pair<relation_ab, int>> abs;

    static void declare_usage(cxx_param_list & pl) {
        param_list_decl_usage(pl, "prec", "initial precision");
        param_list_decl_usage(pl, "ideals", "ideal indices (and valuations)");
        param_list_decl_usage(pl, "renumber", "renumber table");
        param_list_decl_usage(pl, "ab", "(a,b) pairs (and exponents)");
        param_list_decl_usage(pl, "e", "root order");
        param_list_decl_usage(pl, "side", "side");
    }

    explicit mgmsqrt(cxx_param_list & pl)
        : cpoly(pl)
        , R(cpoly, param_list_parse_mandatory<std::string>(pl, "renumber"), true)
        , e(param_list_parse_mandatory<int>(pl, "e"))
        , side(param_list_parse_mandatory<int>(pl, "side"))
        , K(cpoly->pols[side])
        , OK(K.maximal_order())
    {
        K.bless("alpha");
        {
            fmt::print("# Reading (a,b) pairs and exponents\n");
            std::ifstream ab(param_list_parse_mandatory<std::string>(pl, "ab"));
            for(std::string s; std::getline(ab, s); ) {
                strip(s);
                if (s.empty() || s[0] == '#')
                    continue;
                decltype(abs)::value_type x;
                if (!(std::istringstream(s) >> x.first >> x.second))
                    break;
                abs.emplace_back(x);
            }
        }
        {
            fmt::print("# Reading ideal indices and exponents\n");
            std::ifstream Iv(param_list_parse_mandatory<std::string>(pl, "ideals"));
            for(std::string s; std::getline(Iv, s); ) {
                strip(s);
                if (s.empty() || s[0] == '#')
                    continue;
                decltype(ideals)::value_type x;
                if (!(std::istringstream(s) >> x.first >> x.second))
                    break;
                ideals.emplace_back(x);
            }
        }
    }
};

int main(int argc, char const * argv[])
{
    setvbuf(stdout, nullptr, _IONBF, 0);
    setvbuf(stderr, nullptr, _IONBF, 0);
    cxx_param_list pl;

    cxx_cado_poly::declare_usage(pl);
    mgmsqrt::declare_usage(pl);

    /*
    configure_switches(pl);
    configure_aliases(pl);
    */

    argv++, argc--;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fmt::print(stderr, "Unhandled parameter {}\n", argv[0]);
        param_list_print_usage(pl, argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    param_list_print_command_line(stdout, pl);

    mgmsqrt MTY(pl);


    return 0;
}
