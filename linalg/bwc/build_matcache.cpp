#include "cado.h" // IWYU pragma: keep

#include <cstring>
#include <cstdio>
#include <cstdlib>

#include <memory>
#include <string>
#include <stdexcept>
#include <utility>

#include <gmp.h>
#include "fmt/format.h"

#include "cxx_mpz.hpp"
#include "matmul.hpp"
#include "macros.h"
#include "arith-generic.hpp"
#include "portability.h" // asprintf // IWYU pragma: keep
#include "params.h"
#include "matrix_u32.hpp"   // for matrix_u32


static void declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "matrix-file", "matrix file to work with");
    param_list_decl_usage(pl, "prime", "characteristic of the base field [default=2]");
    param_list_decl_usage(pl, "direction", "direction of the product, left for v*M, right for M*v [default=left for p=2, right otherwise]");
    param_list_decl_usage(pl, "withcoeffs", "whether we have coefficients in the matrix (i.e. not only 1's). Defaults to 1 (true) for p>2");
    param_list_decl_usage(pl, "impl", "name of the implementation backend. Defaults to bucket for p==2, basicp for p>2");
    param_list_decl_usage(pl, "groupsize", "number of vectors to consider together (defaults to 64 for p==2, 1 for p>2)");
    param_list_decl_usage(pl, "tmpdir", "directory where matrix cache file is saved (defaults to /tmp)\n");
}

struct direction_flag {
    int value;
    bool operator==(int x) const { return value == x; }
};

template<>
int param_list_parse<direction_flag>(param_list_ptr pl,
        std::string const & arg,
        direction_flag & D)
{
    std::string r;
    if (!param_list_parse(pl, arg, r))
        return 0;
    if (r == "left" || r == "LEFT") {
        D.value = 0;
    } else if (r == "right" || r == "RIGHT") {
        D.value = 1;
    } else {
        throw std::runtime_error(fmt::format("Wrong argument for direction flag ({}), must be left or right", r));
    }
    return 1;
}



int main(int argc, char const * argv[])
{
    const char *argv0 = argv[0];

    cxx_mpz prime;
    mpz_set_ui(prime, 2);
    std::string matrixfile;
    std::string tmpdir = "/tmp";

    setvbuf(stdout, nullptr, _IONBF, 0);
    setvbuf(stderr, nullptr, _IONBF, 0);

    cxx_param_list pl;

    declare_usage(pl);


    for(argv++, argc-- ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    param_list_parse(pl, "prime", prime);

    int withcoeffs = mpz_cmp_ui(prime, 2) != 0;     /* 0 == no coeffs */
    int const groupsize = mpz_cmp_ui(prime, 2) == 0 ? 64 : 1;
    direction_flag direction { mpz_cmp_ui(prime, 2) != 0 };     // 0 = left
    std::string impl = mpz_cmp_ui(prime, 2) == 0 ? "bucket" : "basicp";

    param_list_parse(pl, "withcoeffs", withcoeffs);
    param_list_parse(pl, "direction", direction);
    param_list_parse(pl, "tmpdir", tmpdir);
    param_list_parse(pl, "impl", impl);

    if (!param_list_parse(pl, "matrix-file", matrixfile)) {
        fprintf(stderr, "Error: argument matrix-file is mandatory\n");
        exit(EXIT_FAILURE);
    }
    param_list_warn_unused(pl);

    std::unique_ptr<arith_generic> const xx(arith_generic::instance(prime, groupsize));

    if (direction == 1) {
        fprintf(stderr, "Saving cache for matrix-times-vector\n");
    } else {
        fprintf(stderr, "Saving cache for vector-times-matrix\n");
    }

    /* build a file name for the cache file */
    /* TODO balancing_write has almost identical code that we could refactor
     */
    std::string locfile;
    {
        auto it = matrixfile.rfind('/');
        it = (it == std::string::npos) ? 0 : (it + 1);
        auto basename = matrixfile.substr(it);
        if ((it = basename.rfind(".bin")) != std::string::npos) {
            basename.erase(it, basename.size());
        }
        locfile = fmt::format("{}/{}", tmpdir, basename);
    }

    auto mm = matmul_interface::create(
                xx.get(), 0, 0, locfile, impl, pl, direction.value);

    /* uh ? */
    ASSERT_ALWAYS(mm->store_transposed == !direction.value);

    auto matrix = matrix_u32::from_file(
            matrixfile,
            matrix_u32::transpose_option { mm->store_transposed },
            matrix_u32::withcoeffs_option { !xx->is_characteristic_two() });

    mm->dim = std::get<1>(matrix);
    mm->ncoeffs = std::get<2>(matrix);
    mm->build_cache(std::move(std::get<0>(matrix)));
    mm->save_cache();
}
