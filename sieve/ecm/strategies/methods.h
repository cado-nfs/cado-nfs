#ifndef __METHODS_H__
#define __METHODS_H__

#include "makestrat.h"

cofac_method_srcptr get_method_naive(int *cofac_range, prior_srcptr prior,
     float *acc_fail1, float *acc_fail5, float *acc_fail7, float *acc_fail11, 
     ppm1_history_ptr ppm1_history);

#endif   /* __METHODS_H__ */
