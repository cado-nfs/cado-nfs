#ifndef CADO_SINGLETON_REMOVAL_H
#define CADO_SINGLETON_REMOVAL_H

#include "purge_matrix.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void singleton_removal_oneiter_mono (purge_matrix_ptr mat);
void singleton_removal_oneiter_mt (purge_matrix_ptr, unsigned int);
int64_t singleton_removal (purge_matrix_ptr, unsigned int, int);

#ifdef __cplusplus
}
#endif


#endif /* CADO_SINGLETON_REMOVAL_H */
