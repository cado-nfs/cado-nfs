#ifndef CADO_RELATION_TOOLS_H
#define CADO_RELATION_TOOLS_H

#include <stdint.h>    // for int64_t, uint64_t
#include "typedefs.h"  // for p_r_values_t

#ifdef __cplusplus
extern "C" {
#endif

p_r_values_t relation_compute_r (int64_t a, uint64_t b, p_r_values_t p);
extern char * u64toa16 (char *p, uint64_t m);
extern char * u64toa10 (char *p, uint64_t m);
extern char * d64toa10 (char *p, int64_t m);
extern char * d64toa16 (char *p, int64_t m);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_RELATION_TOOLS_H */
