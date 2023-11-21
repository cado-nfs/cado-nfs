#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <inttypes.h>
#include <stdlib.h>
#include "parallelizing_info.hpp"
#include "xvectors.hpp"
#include "macros.h"              // for FATAL_ERROR_CHECK

typedef int (*sortfunc_t) (const void *, const void *);

int uint32_cmp(const uint32_t * xa, const uint32_t * xb)
{
    if (*xa < *xb) {
        return -1;
    } else if (*xb < *xa) {
        return 1;
    }
    return 0;
}

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
                for(unsigned int j = 0 ; j < nx ; j++) {
                    xs[i*nx+j] = gmp_urandomm_ui(rstate, nr);
                }
                /* Make sure that there's no collision. Not that it
                 * matters so much, but at times the X vector is set with
                 * set_ui, and later on used in an additive manner. Plus,
                 * it does not make a lot of sense to have duplicates,
                 * since that amounts to having nothing anyway...
                 */
                qsort(xs+i*nx,nx,sizeof(uint32_t),(sortfunc_t)uint32_cmp);
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

void load_x(uint32_t ** xs, unsigned int m, unsigned int *pnx,
        parallelizing_info_ptr pi)
{
    FILE * f = NULL;
    int rc = 0;

    /* pretty much the same deal as above */
    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        f = fopen("X", "r");
        FATAL_ERROR_CHECK(f == NULL, "Cannot open X for reading");
        rc = fscanf(f, "%u", pnx);
        FATAL_ERROR_CHECK(rc != 1, "short read in file X");
    }
    pi_bcast(pnx, 1, BWC_PI_UNSIGNED, 0, 0, pi->m);
#ifdef __COVERITY__
    __coverity_mark_pointee_as_sanitized__(pnx, LOOP_BOUND);
#endif
    *xs = (unsigned int *) malloc(*pnx * m * sizeof(unsigned int));
    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        for (unsigned int i = 0, k = 0 ; i < m; i++) {
            for (unsigned int j = 0 ; j < *pnx; j++, k++) {
                rc = fscanf(f, "%" SCNu32, &((*xs)[k]));
                if (rc != 1) {
                    fprintf(stderr, "Short read in X, after reading data for %u rows (compared to expected %u)\n",
                            i, m);
                    abort();
                }
            }
        }
        fclose(f);
    }
    pi_bcast(*xs, *pnx * m * sizeof(uint32_t), BWC_PI_BYTE, 0, 0, pi->m);
}

void set_x_fake(uint32_t ** xs, unsigned int m, unsigned int *pnx,
        parallelizing_info_ptr pi)
{
    /* Don't bother. */
    *pnx=3;
    *xs = (unsigned int *) malloc(*pnx * m * sizeof(unsigned int));
    for(unsigned int i = 0 ; i < *pnx*m ; i++) {
        (*xs)[i] = i;
    }
    serialize(pi->m);
}

void save_x(uint32_t * xs, unsigned int m, unsigned int nx, parallelizing_info_ptr pi)
{
    /* Here, we expect that the data is already available to everybody, so
     * no synchronization is necessary.
     */
    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        // write the X vector
        FILE * fx = fopen("X", "w");
        FATAL_ERROR_CHECK(fx == NULL, "Cannot open X for writing");
        fprintf(fx,"%u\n",nx);
        for(unsigned int i = 0 ; i < m ; i++) {
            for(unsigned int k = 0 ; k < nx ; k++) {
                fprintf(fx,"%s%" PRIu32,k?" ":"",xs[i*nx+k]);
            }
            fprintf(fx,"\n");
        }
        fclose(fx);
    }
}

