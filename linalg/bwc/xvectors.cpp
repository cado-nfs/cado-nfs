#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cinttypes>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include "parallelizing_info.hpp"
#include "xvectors.hpp"
#include "macros.h"              // for FATAL_ERROR_CHECK

void setup_x_random(uint32_t * xs,
        unsigned int m, unsigned int nx, unsigned int nr,
        parallelizing_info_ptr pi, gmp_randstate_t rstate)
{
    /* Here, everybody has to agree on an array of random values. The xs
     * pointer is on the stack of each calling thread, so threads must
     * converge to a commmon point of view on the data.
     */
    // job 0 thread 0 decides for everybody.
    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        for(unsigned int i = 0 ; i < m ; i++) {
            for(;;) {
                for(unsigned int j = 0 ; j < nx ; j++)
                    xs[i*nx+j] = gmp_urandomm_ui(rstate, nr);
                if (nx == 1)
                    break;
                /* Make sure that there's no collision. Not that it
                 * matters so much, but at times the X vector is set with
                 * set_ui, and later on used in an additive manner. Plus,
                 * it does not make a lot of sense to have duplicates,
                 * since that amounts to having nothing anyway...
                 */
                std::sort(xs + i * nx, xs + (i + 1) * nx);
                int collision=0;
                for(unsigned int j = 1 ; j < nx ; j++) {
                    if (xs[i*nx+j] == xs[i*nx+j-1]) {
                        collision=1;
                        break;
                    }
                }
                if (!collision)
                    break;
            }
        }
    }
    pi_bcast(xs, nx * m * sizeof(uint32_t), BWC_PI_BYTE, 0, 0, pi->m);
}

std::unique_ptr<uint32_t[]>
load_x(unsigned int m, unsigned int & nx,
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
    std::unique_ptr<uint32_t[]> xs(new unsigned int[nx * m]);
    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        for (unsigned int i = 0, k = 0 ; i < m; i++) {
            for (unsigned int j = 0 ; j < nx; j++, k++) {
                is >> xs[k];
                if (!is) {
                    fprintf(stderr, "Short read in X, after reading data for %u rows (compared to expected %u)\n",
                            i, m);
                    abort();
                }
            }
        }
    }
    pi_bcast(xs.get(), nx * m * sizeof(uint32_t), BWC_PI_BYTE, 0, 0, pi->m);

    return xs;
}

std::unique_ptr<uint32_t[]>
set_x_fake(unsigned int m, unsigned int & nx,
        parallelizing_info_ptr pi)
{
    /* Don't bother. */
    nx=3;
    std::unique_ptr<uint32_t[]> xs(new unsigned int[nx * m]);
    for(unsigned int i = 0 ; i < nx*m ; i++)
        xs[i] = i;
    serialize(pi->m);
    return xs;
}

void save_x(const uint32_t * xs, unsigned int m, unsigned int nx, parallelizing_info_ptr pi)
{
    /* Here, we expect that the data is already available to everybody, so
     * no synchronization is necessary.
     */
    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        // write the X vector
        std::ofstream fx("X");
        FATAL_ERROR_CHECK(!fx, "Cannot open X for writing");
        fx << nx << "\n";
        for(unsigned int i = 0 ; i < m ; i++) {
            for(unsigned int k = 0 ; k < nx ; k++) {
                if (k) fx << " ";
                fx << xs[i*nx+k];
            }
            fx << "\n";
        }
    }
}

