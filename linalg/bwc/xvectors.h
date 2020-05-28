#ifndef XVECTORS_H_
#define XVECTORS_H_

#include <gmp.h>
#include <stdint.h>
#include "parallelizing_info.h"

#ifdef __cplusplus
extern "C" {
#endif

void setup_x_random(uint32_t * xs,
        unsigned int m, unsigned int nx, unsigned int nr,
        parallelizing_info_ptr pi, gmp_randstate_t rstate);
void load_x(uint32_t ** xs, unsigned int m, unsigned int *pnx,
        parallelizing_info_ptr pi);
void save_x(uint32_t * xs, unsigned int m, unsigned int nx,
        parallelizing_info_ptr pi);
void set_x_fake(uint32_t ** xs, unsigned int m, unsigned int *pnx,
        parallelizing_info_ptr pi);

#ifdef __cplusplus
}
#endif

#endif	/* XVECTORS_H_ */
