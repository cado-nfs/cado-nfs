#include "cado.h"
#include <iostream>
#include "fft.h"
#include "macros.h"

using namespace std;


void fti_explain(const struct fft_transform_info * fti)
{
    std::cout <<
        "Tranform info for accumulating " << fti->nacc << " ";
    if (fti->minwrap)
        std::cout << "modular products (" << fti->bits1 << " by " << fti->bits2 << ") mod 2^K\\pm 1 with K>=" << fti->minwrap;
    else
        std::cout << "integer products (" << fti->bits1 << " by " << fti->bits2 << ")";
    std::cout << "; inputs splits in " << fti->bits << "-bit pieces"
        << ", hence 2 polynomials"
        << " of length " << iceildiv(fti->bits1, fti->bits)
        << " and " << iceildiv(fti->bits2, fti->bits)
        << ", multiplied modulo X^" << (1 << (fti->depth + 2)) << "-1"
        << ", in the ring R=Z/(2^" << (fti->w << fti->depth) << "+1)";

    /* 2^w a n-th root of -1 mod 2^(nw)+1
     * sqrt(2)^w a 4n-th root of 1
     */

    std::cout << ", in which ";
    if (fti->w == 1)
        std::cout << "sqrt(2)";
    else if (fti->w == 2)
        std::cout << "2";
    else
        std::cout << "2^" << (fti->w/2);
    std::cout << " is a " << (1 << (fti->depth+2)) << "-th root of 1";
    std::cout << std::endl;

    ASSERT_ALWAYS(fft_transform_info_check(fti));
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
    fft_get_transform_info_mulmod(fti, 1e6, 8e5, 12, 1e6 + 4); fti_explain(fti);
    fft_get_transform_info_mulmod(fti, 1e6, 8e5, 4, 0); fti_explain(fti);
    fft_get_transform_info_mulmod(fti, 1345570, 706750, 4, 0); fti_explain(fti);
    /* This one is a corner case: we have bits above the firstwrap
     * position, yet those do not really wrap */
    fft_get_transform_info_mulmod(fti, 1351600, 721000, 1, 0); fti_explain(fti);

    fft_get_transform_info_mulmod(fti, 1e7, 8e6, 4, 0); fti_explain(fti);
    fft_get_transform_info_mulmod(fti, 14e7, 7e7, 8, 0); fti_explain(fti);
    fft_get_transform_info_mulmod(fti, 6e8, 4e8, 8, 0); fti_explain(fti);
    fft_get_transform_info_mulmod(fti, 8e8, 7e8, 8, 0); fti_explain(fti);
    fft_get_transform_info_mulmod(fti, 14e8, 14e8, 8, 0); fti_explain(fti);

#if ULONG_BITS == 64
    fft_get_transform_info_mulmod(fti, 3e9, 5e9, 8, 0); fti_explain(fti);
    fft_get_transform_info_mulmod(fti, 6e9, 4e9, 8, 0); fti_explain(fti);
#endif

    fft_get_transform_info_mulmod(fti, 10780, 4900, 16, 0); fti_explain(fti);
}

int main()
{
    test_transform_length();
    return 0;
}
