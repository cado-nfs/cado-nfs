#ifndef CADO_ROPT_STAGE2_H
#define CADO_ROPT_STAGE2_H

#include <stdint.h>      // int16_t
#include "ropt_str.hpp"    // ropt_param ...
#include "ropt_tree.h" // MurphyE_pq


/**
 * Sieving array 
 */
struct sievearray_s {
  int16_t *array;
  unsigned int len_i;
  unsigned int len_j;
};
typedef struct sievearray_s sievearray_t[1];


/* -- declarations -- */

void
ropt_stage2 ( ropt_poly const & poly,
              ropt_s2param_ptr s2param,
              ropt_param_srcptr param,
              ropt_info_ptr info,
              MurphyE_pq *global_E_pqueue,
              int w );


#endif /* CADO_ROPT_STAGE2_H */
