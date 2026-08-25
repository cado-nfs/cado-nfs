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
#include "params.hpp"
#include "matrix_u32.hpp"   // for matrix_u32


static void declare_usage(cxx_param_list & pl)
{
    pl.declare_usage("matrix-file", "matrix file to work with");
    pl.declare_usage("prime", "characteristic of the base field [default=2]");
    pl.declare_usage("direction", "direction of the product, left for v*M, right for M*v [default=left for p=2, right otherwise]");
    pl.declare_usage("withcoeffs", "whether we have coefficients in the matrix (i.e. not only 1's). Defaults to 1 (true) for p>2");
    pl.declare_usage("impl", "name of the implementation backend. Defaults to bucket for p==2, basicp for p>2");
    pl.declare_usage("groupsize", "number of vectors to consider together (defaults to 64 for p==2, 1 for p>2)");
    pl.declare_usage("tmpdir", "directory where matrix cache file is saved (defaults to /tmp)\n");
}

struct direction_flag {
    int value;
    bool operator==(int x) const { return value == x; }
};

template<>
struct cado::params::parser<direction_flag> {
    bool operator()(std::string const & r, direction_flag & D) {
        if (r == "left" || r == "LEFT") {
            D.value = 0;
        } else if (r == "right" || r == "RIGHT") {
            D.value = 1;
        } else {
            return false;
            // throw cado::error("Wrong argument for direction flag ({}), must be left or right", r);
        }
        return true;
    }
};

int main(int argc, char const * argv[])
{
    cxx_mpz prime;
    mpz_set_ui(prime, 2);
    std::string matrixfile;
    std::string tmpdir = "/tmp";

    setvbuf(stdout, nullptr, _IONBF, 0);
    setvbuf(stderr, nullptr, _IONBF, 0);

    cxx_param_list pl;

    declare_usage(pl);

    pl.process_command_line(argc, argv, false);

    pl.parse("prime", prime);

    int withcoeffs = mpz_cmp_ui(prime, 2) != 0;     /* 0 == no coeffs */
    int const groupsize = mpz_cmp_ui(prime, 2) == 0 ? 64 : 1;
    direction_flag direction { mpz_cmp_ui(prime, 2) != 0 };     // 0 = left
    std::string impl = mpz_cmp_ui(prime, 2) == 0 ? "bucket" : "basicp";

    pl.parse("withcoeffs", withcoeffs);
    pl.parse("direction", direction);
    pl.parse("tmpdir", tmpdir);
    pl.parse("impl", impl);

    if (!pl.parse("matrix-file", matrixfile))
        pl.fail("Error: argument matrix-file is mandatory\n");

    pl.warn_unused();

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
