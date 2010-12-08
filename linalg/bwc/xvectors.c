#include "cado.h"
#include <stdio.h>
#include <inttypes.h>
#include "bwc_config.h"
#include "parallelizing_info.h"
#include "xvectors.h"
#include "utils.h"
#include "filenames.h"
#include "balancing.h"

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
        parallelizing_info_ptr pi)
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
                    xs[i*nx+j] = myrand() % nr;
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
    global_broadcast(pi->m, xs, nx * m * sizeof(unsigned int), 0, 0);
}

void load_x(uint32_t * xs, unsigned int m, unsigned int nx,
        parallelizing_info_ptr pi, balancing_ptr bal)
{
    /* pretty much the same deal as above */
    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        char * xn;
        int rc;
        rc = asprintf(&xn, X_FILE_BASE_PATTERN, bal->h->checksum);
        ASSERT_ALWAYS(rc >= 0);
        FILE * f = fopen(xn, "r");
        FATAL_ERROR_CHECK(f == NULL, "Cannot open "X_FILE_BASE_PATTERN" for reading");
        unsigned int nx_file;
        rc = fscanf(f, "%u", &nx_file);
        FATAL_ERROR_CHECK(rc != 1, "short read in file X");
        if (nx == 0) {
            nx = nx_file;       // nx means auto-detect.
        }
        FATAL_ERROR_CHECK(nx != nx_file, X_FILE_BASE_PATTERN " has bad nx value");
        for (unsigned int i = 0 ; i < nx * m; i++) {
            rc = fscanf(f, "%" SCNu32, &(xs[i]));
            FATAL_ERROR_CHECK(rc != 1, "short read in " X_FILE_BASE_PATTERN);
        }
        fclose(f);
        free(xn);
    }
    global_broadcast(pi->m, xs, nx * m * sizeof(unsigned int), 0, 0);
}

void save_x(uint32_t * xs, unsigned int m, unsigned int nx, parallelizing_info_ptr pi, balancing_ptr bal)
{
    /* Here, we expect that the data is already available to eveybody, so
     * no synchronization is necessary.
     */
    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        // write the X vector
        char * xn;
        int rc;
        rc = asprintf(&xn, X_FILE_BASE_PATTERN, bal->h->checksum);
        ASSERT_ALWAYS(rc >= 0);
        FILE * fx = fopen(xn,"w");
        FATAL_ERROR_CHECK(fx == NULL, "Cannot open "X_FILE_BASE_PATTERN" for writing");
        fprintf(fx,"%u\n",nx);
        for(unsigned int i = 0 ; i < m ; i++) {
            for(unsigned int k = 0 ; k < nx ; k++) {
                fprintf(fx,"%s%" PRIu32,k?" ":"",xs[i*nx+k]);
            }
            fprintf(fx,"\n");
        }
        fclose(fx);
        free(xn);
    }
}
