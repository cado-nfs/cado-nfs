#include "cado.h" // IWYU pragma: keep
#include <cmath>
#include <cstdlib>
#include <climits>
#include <mutex>
#include <vector>
#include "getprime.h"
#include "ularith.h"
#include "gpf.h"
#include "macros.h"     // ASSERT ASSERT_EXPENSIVE

/* clumsy: we want a c-only accessor, but an automatic c++ dtor */
unsigned int * gpf;
std::vector<unsigned int> gpf_cxx;

unsigned int gpf_safe_get(const unsigned long i)
{
    ASSERT(i < gpf_cxx.size());
    return gpf_cxx[i];
}

void gpf_init(unsigned int m)
{
    static std::mutex mm;
    std::lock_guard<std::mutex> foo(mm);
    if (m < gpf_cxx.size())
        return;

    std::vector<unsigned int> res;

    res.assign(m+1, 0);

    /* When we increase gpf_max we could and probably should keep
       the old gpf data and sieve only over the newly added area.
       For now, Q&D: sieve everything again. */
    for (size_t i = 0; i <= m; i++) {
        res[i] = i;
    }
    prime_info pi;
    prime_info_init (pi);
    const unsigned int max_sieve = sqrt(m + 1);

    /* Do 2 separately so the compiler can use bit shift for division */
    if (2 <= max_sieve) {
        for (unsigned int i = 2; i <= m ; i += 2) {
            ASSERT(res[i] % 2 == 0);
            while (res[i] > 2) {
                res[i] /= 2;
                if (res[i] % 2 != 0)
                    break;
            }
        }
    }

    for (unsigned int p = getprime_mt (pi) ; p <= max_sieve; p = getprime_mt (pi)) {
        const unsigned int inv_p = ularith_invmod(p);
        const unsigned int lim_p = UINT_MAX / p;
        for (unsigned int i = 2*p; i <= m ; i += p) {
            ASSERT_EXPENSIVE(res[i] % p == 0);
            ASSERT_EXPENSIVE(res[i] / p == res[i] * inv_p);

            while (res[i] > p) {
                unsigned int candidate = res[i] * inv_p;
                if (candidate > lim_p)
                    break;
                res[i] = candidate;
            }
        }
    }
    prime_info_clear (pi);

    /* modify the globals only now. gpf_get only uses the global pointer
     * gpf, so it should be safe to do so (maybe not as helgrind sees it,
     * but in the most common memory model, we're not frightened, as a
     * pointer-size value can be read in one atomic op anyway). */
    std::swap(res, gpf_cxx);
    gpf = gpf_cxx.data();
}
