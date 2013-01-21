#define _GNU_SOURCE
#include "cado.h"
#include "fake_fft.h"
#include "gf2x.h"
#include "utils.h"


/* Number of words holding B bits ; better naming sought. */
#define BITS_TO_WORDS(B,W)      iceildiv((B),(W))

/* nF, nG : number of coefficients */
void fake_init(fake_info_ptr p, size_t nF, size_t nG, ...)
{
    p->n1 = nF;
    p->n2 = nG;
    p->n3 = nF + nG - 1;
    size_t nc = 0;
    if (p->n1 > nc) nc = p->n1;
    if (p->n2 > nc) nc = p->n2;
    p->size = 2 * BITS_TO_WORDS(nc, ULONG_BITS);
}

/* n is a number of coefficients ! */
void fake_dft(fake_info_srcptr p MAYBE_UNUSED, fake_ptr dst, unsigned long * src, size_t n) {
    ASSERT(n <= p->n1 || n <= p->n2);
    size_t s = BITS_TO_WORDS(n, ULONG_BITS);
    memcpy(dst, src, s * sizeof(unsigned long));
    if (n % ULONG_BITS) {
        /* Just as we are computing this assertion, we could easily mask out
         * the bits ourselves. However, our interface mandates that the high
         * bits be cleared in any case. So make sure we properly enforce this
         * constraint.
         */
        ASSERT((src[s-1] & ~((1UL << (n % ULONG_BITS)) - 1)) == 0);
    }
    memset(dst + s, 0, (p->size - s) * sizeof(unsigned long));
}

/* n is a number of coefficients ! -- requiring less than the compute size is
 * okay, as long as it's understood as a means for truncating the data. So we
 * don't do checking for zero high bits.
 */
void fake_ift(fake_info_srcptr p MAYBE_UNUSED, unsigned long * dst, size_t n, fake_srcptr src) {
    ASSERT(n <= p->n3);
    size_t t = BITS_TO_WORDS(n, ULONG_BITS);
    memcpy(dst, src, t * sizeof(unsigned long));
}
void fake_compose(fake_info_srcptr p MAYBE_UNUSED, fake_ptr dst, fake_srcptr s1, fake_srcptr s2) {
    size_t n1 = BITS_TO_WORDS(p->n1, ULONG_BITS);
    size_t n2 = BITS_TO_WORDS(p->n2, ULONG_BITS);
    gf2x_mul(dst, s1, n1, s2, n2);
}
void fake_addcompose(fake_info_srcptr p MAYBE_UNUSED, fake_ptr dst, fake_srcptr s1, fake_srcptr s2) {
    size_t n1 = BITS_TO_WORDS(p->n1, ULONG_BITS);
    size_t n2 = BITS_TO_WORDS(p->n2, ULONG_BITS);
    unsigned long * h = malloc(p->size * sizeof(unsigned long));
    /* lacking addmul in gf2x, we do some extra allocation */
    memset(h, 0, p->size * sizeof(unsigned long));
    gf2x_mul(h, s1, n1, s2, n2);
    for(unsigned int i = 0 ; i < p->size ; i++) {
        dst[i] ^= h[i];
    }
    free(h);
}
void fake_add(fake_info_srcptr p MAYBE_UNUSED, fake_ptr dst, fake_srcptr s1, fake_srcptr s2) {
    size_t i;
    for(i = 0 ; i < p->size ; i++) {
        dst[i] = s1[i] ^ s2[i];
    }
}

void fake_cpy(fake_info_srcptr p MAYBE_UNUSED, fake_ptr dst, fake_srcptr s) {
    memcpy(dst, s, (p->size)*sizeof(unsigned long));
}

size_t fake_size(fake_info_srcptr p) {
    return p->size;
}

void fake_init_similar(fake_info_ptr o, size_t bits_a, size_t bits_b, fake_info_srcptr other MAYBE_UNUSED)
{
    fake_init(o, bits_a, bits_b);
}

int fake_compatible(fake_info_srcptr o1 MAYBE_UNUSED, fake_info_srcptr o2 MAYBE_UNUSED)
{
    return 1;
}

/* vim: set sw=4 sta et: */
