#include "cado.h" // IWYU pragma: keep

#include <cstdint>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>

#include <vector>

#include <gmp.h>
#include "fmt/base.h"

#include "gmp_aux.h"
#include "macros.h"
#include "parallelizing_info.hpp"
#include "utils_cxx.hpp"
#include "xvectors.hpp"

std::vector<uint32_t>
setup_x_random(unsigned int m, unsigned int nx, unsigned int nr,
               parallelizing_info_ptr pi, cxx_gmp_randstate & rstate,
               std::vector<unsigned int> const & forced)
{
    std::vector<uint32_t> xs(nx * m);

    if (forced.size() > m) {
        fmt::print(stderr,
                "Warning: we will place coefficients in X "
                "to account for the {} zero columns, but we will "
                "not be able to do so in a way that completely "
                "the rank defect, since m={} is less than {}\n",
                forced.size(), m, forced.size());
    }

    /* Here, everybody has to agree on an array of random values. The xs
     * pointer is on the stack of each calling thread, so threads must
     * converge to a commmon point of view on the data.
     */
    // job 0 thread 0 decides for everybody.
    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        /* Some parts here deserve a comment. First, we want to make sure
         * that none of the m vectors that we generate has duplicate
         * non-zero positions (see below).
         * Furthermore, there's a tricky logic in actually placing the
         * "forced" non-zero coordinates: we don't want them to appear in
         * batches over the same vector. It's much better to spread them
         * over several. Hence the funny transposition (i+j*m <-->
         * i*nx+j)
         */
        for (bool collision = true ; collision ;) {
            collision = false;
            unsigned int c = 0;
            for (unsigned int i = 0; i < m && !collision ; i++) {
                for (unsigned int j = 0; j < nx && !collision ; j++) {
                    if (i + j * m < forced.size()) {
                        ASSERT_ALWAYS(c < forced.size());
                        xs[i * nx + j] = forced[c++];
                    } else {
                        xs[i * nx + j] = gmp_urandomm_ui(rstate, nr);
                    }
                }
                if (nx == 1)
                    break;
                /* Make sure that there's no collision. Not that it
                 * matters so much, but at times the X vector is set with
                 * set_ui, and later on used in an additive manner. Plus,
                 * it does not make a lot of sense to have duplicates,
                 * since that amounts to having nothing anyway...
                 */
                std::ranges::sort(xs.data() + i * nx, xs.data() + (i + 1) * nx);
                for (unsigned int j = 1; j < nx; j++) {
                    if (xs[i * nx + j] == xs[i * nx + j - 1]) {
                        collision = true;
                        break;
                    }
                }
            }
        }
    }
    pi_bcast(xs.data(), xs.size() * sizeof(uint32_t), BWC_PI_BYTE, 0, 0, pi->m);

    return xs;
}

std::vector<uint32_t> load_x(unsigned int m, unsigned int & nx,
                                   parallelizing_info_ptr pi)
{
    std::ifstream is;

    /* pretty much the same deal as above */
    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        is.open("X");
        FATAL_ERROR_CHECK(!is, "Cannot open X for reading");
        is >> nx;
        FATAL_ERROR_CHECK(!is, "short read in file X");
    }
    pi_bcast(&nx, 1, BWC_PI_UNSIGNED, 0, 0, pi->m);
#ifdef __COVERITY__
    __coverity_mark_pointee_as_sanitized__(pnx, LOOP_BOUND);
#endif
    std::vector<uint32_t> xs(nx * m);
    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        for (unsigned int i = 0, k = 0; i < m; i++) {
            for (unsigned int j = 0; j < nx; j++, k++) {
                is >> xs[k];
                if (!is) {
                    fprintf(stderr,
                            "Short read in X, after reading data for %u rows "
                            "(compared to expected %u)\n",
                            i, m);
                    abort();
                }
            }
        }
    }
    pi_bcast(xs.data(), xs.size() * sizeof(uint32_t), BWC_PI_BYTE, 0, 0, pi->m);

    return xs;
}

std::vector<uint32_t> set_x_fake(unsigned int m, unsigned int & nx,
                                       parallelizing_info_ptr pi)
{
    /* Don't bother. */
    nx = 3;
    std::vector<uint32_t> xs(nx * m);
    for (unsigned int i = 0; i < nx * m; i++)
        xs[i] = i;
    serialize(pi->m);
    return xs;
}

void save_x(std::vector<uint32_t> const & xs, unsigned int m, unsigned int nx,
            parallelizing_info_ptr pi)
{
    /* Here, we expect that the data is already available to everybody, so
     * no synchronization is necessary.
     */
    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        ASSERT_ALWAYS(xs.size() == m * nx);
        // write the X vector
        std::ofstream fx("X");
        FATAL_ERROR_CHECK(!fx, "Cannot open X for writing");
        fx << nx << "\n";
        for (unsigned int i = 0; i < m; i++) {
            fx << join(xs.data() + i * nx, xs.data() + (i + 1) * nx, " ") << "\n";
        }
    }
}
