#include <stdint.h>
#include <stdlib.h>
#include "gf2x.h"

/*
 *
 * gf2x did not, prior to version 1.3, have a GF2X_VERSION macro defined.
 *
 * In order to check for the older gf2x versions, you may try to
 * *compile* this source file with one of the CHECK macros enabled.
 *
 * In particular, with -DCHECK7, you are checking for gf2x version 1.2 or
 * later.
 *
 * Such checks may be used in configure-time checks (autoconf, cmake).
 *
version             <---   -DCHECK#    --->
                     1   2   3   4   5   6   7
0.9.1 to 0.9.5      OK  NOK NOK NOK OK  NOK NOK
0.9.6               NOK OK  NOK NOK OK  NOK NOK
1.0   to 1.1        NOK OK  OK  NOK OK  NOK NOK
1.2                 NOK NOK NOK OK  OK  NOK OK
1.2-fft             NOK NOK NOK OK  OK  OK  OK
1.2-LGPL            NOK NOK NOK OK  NOK NOK OK
1.2-LGPL-fft        NOK NOK NOK OK  NOK OK  OK

later versions can check GF2X_VERSION
*/

#ifdef CHECK1
int version_0_9_1_to_0_9_5_ok() {
    extern int gf2x_mul_fft0;
    return gf2x_mul_fft0 == 0;
}
#endif
#ifdef CHECK2
int version_0_9_6_to_1_1_ok() {
    extern int gf2x_tfft_alloc;
    return gf2x_tfft_alloc == 0;
}
#endif
#ifdef CHECK3
int version_1_0_to_1_1_ok() {
    extern int gf2x_tfft_compatible;
    return gf2x_tfft_compatible == 0;
}
#endif
#ifdef CHECK4
int version_1_2_ok() {
    extern int gf2x_ternary_fft_compatible;
    return gf2x_ternary_fft_compatible == 0;
}
#endif
#ifdef CHECK5
int version_1_2_is_GPL_flavor() {
    extern int gf2x_mul_tc3u;
    return gf2x_mul_tc3u == 0;
}
#endif
#ifdef CHECK6
int version_1_2_has_fft_interface() {
    extern int gf2x_fake_fft_init;
    return gf2x_fake_fft_init == 0;
}
#endif
#ifdef CHECK7
int version_1_2_or_later() {
#ifdef GF2X_VERSION_MAJOR
    /* gf2x 1.3 or later define this */
    return 0;
#else
    /* This will fail to link for gf2x versions prior to 1.2 */
    extern int gf2x_ternary_fft_compatible;
    return gf2x_ternary_fft_compatible == 0;
#endif
}
#endif
#ifdef CHECK8
/* This is a slightly more verbose version of the above, if you fear the
 * linker may be averse to dirty casts.  */
#ifndef GF2X_VERSION_MAJOR
extern "C" {
struct gf2x_ternary_fft_info_s;
typedef const struct gf2x_ternary_fft_info_s * gf2x_ternary_fft_info_srcptr;
int gf2x_ternary_fft_compatible(gf2x_ternary_fft_info_srcptr o1, gf2x_ternary_fft_info_srcptr o2); 
}
#endif
/* This will fail to link for versions prior to v1.2. */
void fun() { gf2x_ternary_fft_compatible(0, 0); }
#endif

int main() {
    unsigned long a = 0xdeadbeef, b=0xcafecafe, c[2];
    /* result: 0x5fff29e51e50684a; */
    gf2x_mul(c,&a,1,&b,1);
    return ((uint32_t) c[0] == 0x1e50684a) ? 0 : EXIT_FAILURE;
}

