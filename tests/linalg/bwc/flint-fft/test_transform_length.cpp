#include "cado.h"
#include <iostream>
#include "fft.h"
#include "macros.h"

using namespace std;


extern "C" {
    extern void fft_transform_info_init_mulmod_inner(struct fft_transform_info * fti, mp_bitcnt_t bits1, mp_bitcnt_t bits2, unsigned int nacc, mp_bitcnt_t minwrap);
}

void test_transform_length()
{
    /* imagine various multiplication lengths, and run the transform
     * length selection algorithm to see whether the corner cases are
     * triggered. Implicitly, we rely on the ASSERTs there to make sure
     * that all the required inequalities hold, because we're doing no
     * check here.
     */
    struct fft_transform_info fti[1];
    fft_transform_info_init_mulmod_inner(fti, 1e6, 8e5, 12, 1e6 + 4); std::cout << fti->explain() << "\n";
    fft_transform_info_init_mulmod_inner(fti, 1e6, 8e5, 4, 0); std::cout << fti->explain() << "\n";
    fft_transform_info_init_mulmod_inner(fti, 1345570, 706750, 4, 0); std::cout << fti->explain() << "\n";
    /* This one is a corner case: we have bits above the firstwrap
     * position, yet those do not really wrap */
    fft_transform_info_init_mulmod_inner(fti, 1351600, 721000, 1, 0); std::cout << fti->explain() << "\n";

    fft_transform_info_init_mulmod_inner(fti, 1e7, 8e6, 4, 0); std::cout << fti->explain() << "\n";
    fft_transform_info_init_mulmod_inner(fti, 14e7, 7e7, 8, 0); std::cout << fti->explain() << "\n";
    fft_transform_info_init_mulmod_inner(fti, 6e8, 4e8, 8, 0); std::cout << fti->explain() << "\n";
#if ULONG_BITS == 64
    /* The following two are not totally out of question of 32-bits,
     * since the result would fit. However, our preferred size choices
     * lead us to overflows.
     */
    fft_transform_info_init_mulmod_inner(fti, 8e8, 7e8, 8, 0); std::cout << fti->explain() << "\n";
    fft_transform_info_init_mulmod_inner(fti, 14e8, 14e8, 8, 0); std::cout << fti->explain() << "\n";
    /* These are definitely 64-bit territory */
    fft_transform_info_init_mulmod_inner(fti, 3e9, 5e9, 8, 0); std::cout << fti->explain() << "\n";
    fft_transform_info_init_mulmod_inner(fti, 6e9, 4e9, 8, 0); std::cout << fti->explain() << "\n";
#endif

    fft_transform_info_init_mulmod_inner(fti, 10780, 4900, 16, 0); std::cout << fti->explain() << "\n";
}

int main()
{
    test_transform_length();
    return 0;
}
