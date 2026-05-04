#ifndef CADO_UTILS_GETPRIME_H
#define CADO_UTILS_GETPRIME_H

// scan-headers: skip

#ifndef MAIN
#include "macros.h"
#else
#define ATTRIBUTE_DEPRECATED
#endif

#ifdef __cplusplus
#include <climits>

#include <limits>
#include <utility>
#endif

struct prime_info_s {
  unsigned long offset;  /* offset for current primes */
  unsigned long current;  /* index of previous prime */
  unsigned int *primes;  /* small primes up to sqrt(p) */
  unsigned long nprimes; /* length of primes[] */
  unsigned char *sieve;  /* sieving table */
  unsigned long len;     /* length of sieving table */
  unsigned int *moduli;  /* offset for small primes */
};
typedef struct prime_info_s prime_info[1];
typedef struct prime_info_s * prime_info_ptr;
typedef const struct prime_info_s * prime_info_srcptr;

#ifdef __cplusplus
extern "C" {
#endif

/* The getprime function returns successive primes, by default only odd
 * primes starting with 3. */
extern void prime_info_init (prime_info_ptr);

/* this second intialization function makes it possible to specify a
 * lower bound that is different from 3. It is even possible to pass 2,
 * in which case getprime_mt(pi) will return 2, 3, 5, 7, 11, ...
 */
extern void prime_info_init_seek(prime_info_ptr pi, unsigned long lower_bound);
extern void prime_info_clear (prime_info_ptr);
extern void prime_info_seek (prime_info_ptr, unsigned long lower_bound);
extern unsigned long getprime_mt (prime_info_ptr);

extern unsigned long getprime (unsigned long) ATTRIBUTE_DEPRECATED;


#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
struct prime_range
{
    unsigned long lower_bound = 2;
    unsigned long upper_bound = ULONG_MAX;
    unsigned long last = 2;
    prime_range(unsigned long pmin = 2, unsigned long pmax = ULONG_MAX)
        : lower_bound(pmin)
        , upper_bound(pmax)
    {
    }

    struct const_iterator {
        prime_info pi;
        unsigned long last;
        explicit const_iterator(unsigned long pmin)
        {
            prime_info_init(pi);
            prime_info_seek(pi, pmin);
            last = getprime_mt(pi);
        }
        const_iterator(const_iterator const &) = delete;
        const_iterator(const_iterator && o) noexcept
        {
            if (&o == this)
                return;
            memcpy(&pi, &o.pi, sizeof(prime_info_s));
            memset(&o.pi, 0, sizeof(prime_info_s));
        }
        const_iterator& operator=(const_iterator const &) = delete;
        const_iterator& operator=(const_iterator &&) = delete;
        ~const_iterator() {
            prime_info_clear(pi);
        }
        unsigned long operator*() const {
	    return last;
        }
        const_iterator& operator++() {
            last = getprime_mt(pi);
            return *this;
        }
        /*
        const_iterator& operator++(int) {
            auto c = *this; ++*this; return c;
        }
        */
        /* when the iteator reaches the range limit, it should stop there
         * and no longer increase, and a priori we won't see it wrap
         * around.
         */
        auto operator<=>(unsigned long p) const { return **this <=> p; }
        
        /* this one is a hack.
         */
        auto operator==(unsigned long p) const { return **this >= p; }
    };
    const_iterator begin() const {
        return const_iterator(lower_bound);
    }
    unsigned long end() const {
        return upper_bound;
    }
};
#endif

#endif	/* CADO_UTILS_GETPRIME_H */
