#include "cado.h" // IWYU pragma: keep
#include "bblas_level2b.hpp"
#include "bblas_mat64.hpp"
#include "misc.h"      // cado_parity64

/**********************************************************************/
/* level 2b: vector times matrices.
 *      mul_o64_6464    (several variants)
 *      mul_o64_T6464   (several variants)
 */

/* implemented here:
 *    - mul_o64_6464 ; multiply a vector by a matrix
 *    - mul_o64_T6464 ; multiply a vector by the transpose of a matrix
 */

/* implements mul_o64_6464 */
void mul_o64_6464_C_lsb(uint64_t * r, uint64_t a, mat64 const & w)/*{{{*/
{
    uint64_t c = 0;
    for (unsigned int i = 0; i < 64; i++) {
	c ^= (w[i] & -(a & UINT64_C(1)));
	a >>= 1;
    }
    *r = c;
}/*}}}*/

/* implements mul_o64_6464 */
void mul_o64_6464_C_msb(uint64_t *r, uint64_t a, mat64 const & w)/*{{{*/
{
    uint64_t c = 0;
    for (int i = 64 - 1; i >= 0; i--) {
        c ^= (w[i] & (((int64_t) a) >> (64 - 1)));
        a <<= 1;
    }
    *r = c;
}/*}}}*/

/* implements mul_o64_T6464 */
void mul_o64_T6464_C_parity(uint64_t * w, uint64_t a, mat64 const & b)/*{{{*/
{
    // Uses unoptimized __builtin_parityl function -- maybe better with gcc 4.3
    // note that popcnt is faster in asm than the more restricted parity
    // functions. So if it's available, it should be tested.
    uint64_t c = 0;
    for (unsigned int i = 0; i < 64; i++) {
        uint64_t p = cado_parity64(a & b[i]);
	c ^= p << i;
    }
    *w = c;
}/*}}}*/

/* implements mul_o64_T6464 */
/* {{{ mul_o64_T6464_C_parity3 */

/* This is stolen from code by D. Harvey. (GPL, thus can't stay here) */
static inline uint64_t _parity64_helper2(const uint64_t* buf, uint64_t a)
{
#define XMIX32(b, a) (((((a) << 32) ^ (a)) >> 32) + \
                     ((((b) >> 32) ^ (b)) << 32))
#define XMIX16(b, a) (((((a) >> 16) ^ (a)) & 0x0000FFFF0000FFFFll) +     \
                     ((((b) << 16) ^ (b)) & 0xFFFF0000FFFF0000ll));
#define XMIX8(b, a) (((((a) >> 8) ^ (a)) & 0x00FF00FF00FF00FFll) + \
                    ((((b) << 8) ^ (b)) & 0xFF00FF00FF00FF00ll));
#define XMIX4(b, a) (((((a) >> 4) ^ (a)) & 0x0F0F0F0F0F0F0F0Fll) + \
                    ((((b) << 4) ^ (b)) & 0xF0F0F0F0F0F0F0F0ll));
#define XMIX2(b, a) (((((a) >> 2) ^ (a)) & 0x3333333333333333ll) + \
                    ((((b) << 2) ^ (b)) & 0xCCCCCCCCCCCCCCCCll));
#define XMIX1(b, a) (((((a) >> 1) ^ (a)) & 0x5555555555555555ll) + \
                    ((((b) << 1) ^ (b)) & 0xAAAAAAAAAAAAAAAAll));
   uint64_t a0, a1, b0, b1, c0, c1;
   a0 = XMIX32(buf[0x20] & a, buf[0x00] & a);
   a1 = XMIX32(buf[0x30] & a, buf[0x10] & a);
   b0 = XMIX16(a1, a0);
   a0 = XMIX32(buf[0x28] & a, buf[0x08] & a);
   a1 = XMIX32(buf[0x38] & a, buf[0x18] & a);
   b1 = XMIX16(a1, a0);
   c0 = XMIX8(b1, b0);
   a0 = XMIX32(buf[0x24] & a, buf[0x04] & a);
   a1 = XMIX32(buf[0x34] & a, buf[0x14] & a);
   b0 = XMIX16(a1, a0);
   a0 = XMIX32(buf[0x2C] & a, buf[0x0C] & a);
   a1 = XMIX32(buf[0x3C] & a, buf[0x1C] & a);
   b1 = XMIX16(a1, a0);
   c1 = XMIX8(b1, b0);
   return XMIX4(c1, c0);
}

void mul_o64_T6464_C_parity3(uint64_t * w, uint64_t a, mat64 const & b)
{
   uint64_t d0, d1, e0, e1;

   d0 = _parity64_helper2(b.data(), a);
   d1 = _parity64_helper2(b.data() + 2, a);
   e0 = XMIX2(d1, d0);

   d0 = _parity64_helper2(b.data() + 1, a);
   d1 = _parity64_helper2(b.data() + 3, a);
   e1 = XMIX2(d1, d0);

   *w = XMIX1(e1, e0);
} /* }}} */

/* final choices. These are static choices at this point, but it should
 * be the result of some tuning, ideally */
/* Given a 64-bit vector a, and a 64*64 matrix w, compute the product.
 *
 * With respect to endianness, we match a&1 with w[0]
 */
void mul_o64_6464(uint64_t * r, uint64_t a, mat64 const & w)
{
    mul_o64_6464_C_lsb(r,a,w);
}
void mul_o64_T6464(uint64_t * w, uint64_t a, mat64 const & b)
{
    mul_o64_T6464_C_parity(w,a,b);
}
