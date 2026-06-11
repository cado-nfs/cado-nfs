#include <cfenv>
#include <cstdlib>
#include <cstdint>
#include <cmath>

struct temporary_round_mode {
    int saved;
    explicit temporary_round_mode(int mode)
        : saved(fegetround())
    {
        fesetround(mode);
    }
    temporary_round_mode(temporary_round_mode const &) = delete;
    temporary_round_mode(temporary_round_mode &&) = delete;
    temporary_round_mode & operator=(temporary_round_mode const &) = delete;
    temporary_round_mode & operator=(temporary_round_mode &&) = delete;
    ~temporary_round_mode() { fesetround(saved); }
};

static inline bool rounding_towards_zero_works()
{
    /* This runtime check can detect dysfunctional math environments.
     * valgrind is one of them, unfortunately.
     */
    temporary_round_mode dummy(FE_TOWARDZERO);
    const double d0 = (uint64_t(1) << 52) + 1;
    const double d1 = std::ldexp(d0, 53);
    /* We can't initialize it as d = d0 + d1 because that wouldn't
     * necessarily obey the rounding mode
     */
    volatile double d = d0;
    d = d + d1;
    return d == d1;
}

int main()
{
    return rounding_towards_zero_works() ? EXIT_SUCCESS : EXIT_FAILURE;
}
