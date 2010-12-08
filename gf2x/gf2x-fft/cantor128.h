#ifndef CANTOR128_H_
#define CANTOR128_H_

#include <stdint.h>
#include <stdlib.h>
#include "gf2x.h"
#include "gf2x/gf2x-impl.h"  /* for GF2X_WORDSIZE */

#ifdef __cplusplus
extern "C" {
#endif

extern void mulCantor128(unsigned long *H, unsigned long *F, size_t Fl, unsigned long *G, size_t Gl);



// now for the gf2x-fft interface.
//
// note that we must implement some sort of fallback on the fake_fft
// thing when the polynomials are too small, or something like this. The
// gf2x-fft interface requires that any size can be used.

struct c128_info_struct {
    unsigned int k;
    size_t n;
};
typedef struct c128_info_struct c128_info_t[1];
typedef struct c128_info_struct * c128_info_ptr;
typedef const struct c128_info_struct * c128_info_srcptr;

#if GF2X_WORDSIZE == 32
typedef unsigned long c128_t[4];
#elif GF2X_WORDSIZE == 64
typedef unsigned long c128_t[2];
#else 
#error "Do not know how to define c128_t"
#endif

typedef c128_t * c128_ptr;
/* aliasing is not our friend here ;-( */
typedef c128_t * c128_srcptr;

extern void c128_init(c128_info_t p, size_t nF, size_t nG, ...);
static inline void c128_clear(c128_info_t p MAYBE_UNUSED) {}
extern c128_ptr c128_alloc(const c128_info_t p, size_t n);
extern void c128_free(
        const c128_info_t p MAYBE_UNUSED,
        c128_ptr x,
        size_t n MAYBE_UNUSED);
extern c128_ptr c128_get(const c128_info_t p, c128_ptr x, size_t k);
extern void c128_zero(const c128_info_t p, c128_ptr x, size_t n);

extern void c128_dft(const c128_info_t p, c128_ptr x, unsigned long * F, size_t nF);
extern void c128_compose(const c128_info_t p,
		c128_ptr y, c128_srcptr x1, c128_srcptr x2);
extern void c128_addcompose(const c128_info_t p,
		c128_ptr y, c128_srcptr x1, c128_srcptr x2);
extern void c128_add(const c128_info_t p,
		c128_ptr y, c128_srcptr x1, c128_srcptr x2);
extern void c128_cpy(const c128_info_t p, c128_ptr y, c128_srcptr x);
extern void c128_ift(const c128_info_t p,
		unsigned long * H, size_t Hl, c128_srcptr h);
extern size_t c128_size(c128_info_srcptr p);
extern void c128_init_similar(c128_info_ptr o, size_t bits_a, size_t bits_b, c128_info_srcptr other);
extern int c128_compatible(c128_info_srcptr o1, c128_info_srcptr o2);

#ifdef __cplusplus
}
#endif

#endif	/* CANTOR128_H_ */
