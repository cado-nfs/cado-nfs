#ifndef ROPT_STAGE2_H
#define ROPT_STAGE2_H


#include "ropt_param.h"
#include "ropt_tree.h"
#include "ropt_arith.h"


/**
 * Sieving array 
 */
typedef struct {
  int16_t *array;
  unsigned int len_i;
  unsigned int len_j;
} _sievearray_t;
typedef _sievearray_t sievearray_t[1];


/* -- declarations -- */

#ifdef __cplusplus
extern "C" {
#endif

void
ropt_stage2 ( ropt_poly_t poly,
              ropt_s2param_t s2param,
              ropt_param_t param,
              ropt_info_t info,
              MurphyE_pq *global_E_pqueue,
              int w );

#if DEBUG
void print_sievearray ( double **A,
                        long A0,
                        long A1,
                        long B0,
                        long B1,
                        unsigned long K_ST,
                        unsigned long J_ST,
                        unsigned long MOD );
#endif

#ifdef __cplusplus
}
#endif


#endif /* ROPT_STAGE2_H */
