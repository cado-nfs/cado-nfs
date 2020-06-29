#ifndef GPF_H_
#define GPF_H_

/* A look-up table of the largest prime factor of a non-negative integer i.

Before calling gpf_get(i), gpf_init(m) must have been called with m >= i.

gpf_get(0) = 0, gpf_get(1) = 1.

*/

#ifdef __cplusplus
extern "C" {
#endif

extern unsigned int *gpf;

void gpf_init(unsigned int);

unsigned int gpf_safe_get(const unsigned long i);
static inline unsigned int gpf_get(const unsigned long i) {
#ifndef NDEBUG
    return gpf_safe_get(i);
#else
    return gpf[i];
#endif
}

/* FIXME: we never ever use it, it seems */
void gpf_clear();

#ifdef __cplusplus
}
#endif

#endif
