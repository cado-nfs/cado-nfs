#ifndef CADO_UTILS_MISC_H
#define CADO_UTILS_MISC_H

#include "cado_config.h"  // for HAVE_GCC_STYLE_AMD64_INLINE_ASM, ULONGLONG_...
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
#include <type_traits>
#include <string>
#include <vector>
#include <utility>
#include <ios>
#endif

#include <gmp.h>
#include "macros.h"
#ifdef __cplusplus
#include "cxx_mpz.hpp"
#endif

// scan-headers: stop here

/* we prefer GMP 5 or later, but the history of the why and how seems
 * lost. It seems that the late GMP-4.3 versions are fine, and the few
 * missing functions are provided as fallbacks */
#if !GMP_VERSION_ATLEAST(4,3,0)
#error "GNU MP 4.3.0 (at least) is required to compile CADO-NFS"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* These are is in misc2.cpp */
double nprimes_interval(double p0, double p1);
double prime_pi_2exp(unsigned int n);
double random_along_prime_distribution(unsigned int bits, gmp_randstate_t rstate);

uint64_t u64_random(gmp_randstate_t buf);
int64_t i64_random(gmp_randstate_t buf);

#define UMAX(A) (0xffffffffffffffffULL >>((8-sizeof(A))<<3U))
#define SMAX(A) (0x7fffffffffffffffLL  >>((8-sizeof(A))<<3U))
#define UMIN(A) (0)
#define SMIN(A) (1ULL<< ((sizeof(A)<<3U)-1))

/* uintmax_t is guaranteed to be larger or equal to uint64_t */
#define strtouint64(nptr,endptr,base) (uint64_t) strtoumax(nptr,endptr,base)

static inline void* pointer_arith(void * a, ptrdiff_t q) {
    return (void*)(((char*)a)+q);
}
static inline const void* pointer_arith_const(const void * a, ptrdiff_t q) {
    return (const void*)(((const char*)a)+q);
}

/* strtoul(), but with const char ** for second argument.
   Otherwise it's not possible to do, e.g., strtoul(p, &p, 10) when p is
   of type const char *
*/
extern unsigned long int strtoul_const(const char *nptr, const char **endptr, const int base);
extern unsigned long long int strtoull_const(const char *nptr, const char **endptr, const int base);

extern char * derived_filename(const char * prefix, const char * what, const char * ext);
extern int has_suffix(const char * path, const char * sfx);

extern char const ** filelist_from_file(const char * basepath, const char * filename,
                                  int typ);
extern void filelist_clear(char const ** filelist);

long get_arg_max(void);
extern int mkdir_with_parents(const char * dir, int fatal);

extern char * path_resolve(const char * progname, char * resolved);

extern void bit_reverse(unsigned long *, const unsigned long *, size_t);

/* k must be a power of 2. Returns the smallest multiple of k greater
 * than or equal to n
 */
static inline unsigned long
next_multiple_of_powerof2(unsigned long n, unsigned long k)
{
    ASSERT((k & (k-1)) == 0);
    return ((n-1)|(k-1)) + 1;
}
static inline unsigned long
next_multiple_of(unsigned long n, unsigned long k)
{
    return iceildiv(n, k) * k;
}
static inline unsigned long integer_sqrt(unsigned long a)
{
    /* returns 2 for a==3, otherwise returns floor(sqrt(a)) */
    for(unsigned long x = a, y = 1, z; ; x = y, y = z) {
        z = (y + a/y) / 2;
        if (z == x || z == y) return y;
    }
}

#ifndef __cplusplus
/* we have a template variant of this for C++ code */
static inline unsigned long next_power_of_2(unsigned long x)
{
    /* round x to the next power of two */
    for( ; x & (x - 1) ; ) {
        unsigned long low = x ^ (x - 1);
        x += (low >> 1U) + 1;
    }
    return x;
}
#endif

/* Best X86 medium memcpy with pointers & size length already cache
   lines aligned on modern X86, so dst/src/lg & 0x3F = 0 */
static inline void aligned_medium_memcpy(void *dst, void *src, size_t lg) {
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
  size_t lg_8bytes = lg >> 3U;
  __asm__ __volatile__ ("cld\nrep movsq\n":"+D"(dst),"+S"(src),"+c"(lg_8bytes));
#else
  memcpy(dst,src,lg);
#endif
}
/*{{{ cado_clz and variants*/
/* those builtins seem to have appeared in 3.4 (April 2004) */
#if GNUC_VERSION_ATLEAST(3,4,0)
#define cado_clzll(x)        __builtin_clzll(x)
#define cado_clzl(x)         __builtin_clzl(x)
#define cado_clz(x)          __builtin_clz(x)
#else
/* provide slow fallbacks */
static inline unsigned int cado_clzll(unsigned long long x)
{
#if ULONGLONG_BITS == 64
        static const int t[4] = { 2, 1, 0, 0 };
        int a = 0;
        int res;
        if (x >> 32U) { a += 32; x >>= 32U; }
        if (x >> 16U) { a += 16; x >>= 16U; }
        if (x >>  8U) { a +=  8; x >>=  8U; }
        if (x >>  4U) { a +=  4; x >>=  4U; }
        if (x >>  2U) { a +=  2; x >>=  2U; }
        res = 64 - 2 - a + t[x];
        return res;
#else
#error "cado_clzll might be wrong here"
#endif
}

static inline unsigned int cado_clzl(unsigned long x)
{
        static const int t[4] = { 2, 1, 0, 0 };
        int a = 0;
        int res;
#if ULONG_BITS == 64
        if (x >> 32U) { a += 32; x >>= 32U; }
#endif
        if (x >> 16U) { a += 16; x >>= 16U; }
        if (x >>  8U) { a +=  8; x >>=  8U; }
        if (x >>  4U) { a +=  4; x >>=  4U; }
        if (x >>  2U) { a +=  2; x >>=  2U; }
        res = GMP_LIMB_BITS - 2 - a + t[x];
        return res;
}

static inline int cado_clz(unsigned int x)
{
        static const unsigned int t[4] = { 2, 1, 0, 0 };
        int a = 0;
        int res;
        if (x >> 16U) { a += 16; x >>= 16U; }
        if (x >>  8U) { a +=  8; x >>=  8U; }
        if (x >>  4U) { a +=  4; x >>=  4U; }
        if (x >>  2U) { a +=  2; x >>=  2U; }
        res = 30 - a + t[x];
        return res;
}
#endif
/*}}}*/
/*{{{ cado_ctz and variants */
#if GNUC_VERSION_ATLEAST(3,4,0)
#define cado_ctzll(x)        __builtin_ctzll(x)
#define cado_ctzl(x)         __builtin_ctzl(x)
#define cado_ctz(x)          __builtin_ctz(x)
#else
/* the following code is correct because if x = 0...0abc10...0, then
   -x = ~x + 1, where ~x = 1...1(1-a)(1-b)(1-c)01...1, thus
   -x = 1...1(1-a)(1-b)(1-c)10...0, and x & (-x) = 0...000010...0 */
static inline unsigned int cado_ctzll(unsigned long long x)
{
  return (ULONGLONG_BITS - 1) - cado_clzll(x & - x);
}
static inline unsigned int cado_ctzl(unsigned long x)
{
  return (ULONG_BITS - 1) - cado_clzl(x & - x);
}
static inline unsigned int cado_ctz(unsigned int x)
{
  return (ULONG_BITS - 1) - cado_clzl(x & - x);
}
#endif
/*}}}*/
/*{{{ cado_parity and variants */
#if GNUC_VERSION_ATLEAST(3,4,0)
#define cado_parityll(x)        __builtin_parityll(x)
#define cado_parityl(x)         __builtin_parityl(x)
#define cado_parity(x)          __builtin_parity(x)
#else
/* slow equivalent */
static inline int cado_parityll(unsigned long long x)
{
#if ULONGLONG_BITS == 64
    x ^= x >> 32U;
    x ^= x >> 16U;
    x ^= x >> 8U;
    x ^= x >> 4U;
    int t[16] = { 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0 };
    return t[x&15];
#else
#error "fix cado_parityll here"
#endif
}
#endif

#if ULONGLONG_BITS == 64
static inline unsigned int cado_ctz64(uint64_t x) { return cado_ctzll(x); }
static inline unsigned int cado_clz64(uint64_t x) { return cado_clzll(x); }
static inline unsigned int cado_parity64(uint64_t x) { return cado_parityll(x); }
#else
#error "need proper equivalents for cado_ctz64 & friends"
#endif
/*}}}*/

/* l{0,2}abs() are specified to result in undefined behaviour if the argument
   is the most negative value. Our safe_l{0,2}abs() fix this.
   Casting to an unsigned type and taking the minus then is valid as unary
   minus of an unsigned type returns an unsigned type and arithmetic on
   unsigned types is defined as working modulo [TYPE]_MAX+1 by C99 6.2.5c9 */
static inline unsigned long long
safe_llabs(const long long n) {
    return (n < 0) ? -(unsigned long long)n : (unsigned long long)n;
}

static inline unsigned long
safe_labs(const long n) {
    return (n < 0) ? -(unsigned long)n : (unsigned long)n;
}

static inline unsigned int
safe_abs(const int n) {
    return (n < 0) ? -(unsigned int)n : (unsigned int)n;
}

static inline uint64_t
safe_abs64(const int64_t n) {
    return (n < 0) ? -(uint64_t)n : (uint64_t)n;
}

size_t file_bytes(FILE * f);

const char *size_disp_fine(size_t s, char buf[16], double cutoff);
const char *size_disp(size_t s, char buf[16]);

/* see below. The C++ code is the first-class citizen, but this proxy can
 * be used in C as well.
 */
extern void subdivide_primes_interval_proxy(unsigned long * r, unsigned long p0, unsigned long p1, size_t n);

extern int mpz_set_from_expression(mpz_ptr f, const char * value);


#ifdef __cplusplus
}
#endif

static inline const char * ok_NOK(int t)
{
    return t ? "ok" : "NOK";
}

static inline const char * ok_NOKNOK(int t)
{
    return t ? "ok" : "NOK NOK NOK NOK";
}

#ifdef __cplusplus
static inline std::string size_disp(size_t s) {
    char buf[16];
    size_disp(s, buf);
    return { buf };
}
static inline std::string size_disp_fine(size_t s, double cutoff) {
    char buf[16];
    size_disp_fine(s, buf, cutoff);
    return { buf };
}

static inline bool ends_with(std::string const & name, std::string const & suffix)
{
    // c++20: return name.ends_with(suffix)
    if (name.size() < suffix.size())
        return false;
    return name.substr(name.size()-suffix.size()) == suffix;
}
#endif

#ifdef __cplusplus
template<typename T>
static inline T next_power_of_2(T x)
{
    /* it's a bit crazy. why not just do:
     *
     * previous_power_of_two(x) {
     *   for(T c ; (c = x & (x-1)) != 0 ; x = c);
     *   return x;
     * }
     * next_power_of_2(x) { return previous_power_of_two(x-1)<<1U; }
     */
    static_assert(
            std::is_same_v<T, unsigned long> ||
            std::is_same_v<T, unsigned int> ||
            std::is_same_v<T, size_t>,
            "not supported for this type");
    /* round x to the next power of two */
    for( ; x & (x - 1) ; ) {
        T low = x ^ (x - 1);
        x += (low >> 1U) + 1;
    }
    return x;
}

template<typename T>
static inline T log2_of_next_power_of_2(T x)
{
    static_assert(
            std::is_same_v<T, unsigned long> ||
            std::is_same_v<T, unsigned int> ||
            std::is_same_v<T, size_t>,
            "not supported for this type");
    T m = 1;
    int i = 0;
    for( ; m && x > m ; i++, m<<=1);
    return i;
}


std::vector<unsigned long> subdivide_primes_interval(unsigned long p0, unsigned long p1, size_t n);
#endif

#ifdef __cplusplus
/* Use in any function that uses iomanip temporarily.
 *
 * Hmm, how much different is it from std::istream::sentry
 */

class IoStreamFlagsRestorer
{
public:
    IoStreamFlagsRestorer(std::ios_base & ioStream)
        : ioStream_(ioStream)
        , flags_(ioStream_.flags())
    {
    }

    ~IoStreamFlagsRestorer()
    {
        ioStream_.flags(flags_);
    }

private:
    std::ios_base & ioStream_;
    std::ios_base::fmtflags const flags_;
};
#endif


#if 0
#ifdef __cplusplus
// declare c++ containers as vector<T,pagealigned_allocator<T>>
template < typename _Tp > class pagealigned_allocator {
      public:
        typedef size_t size_type;
        typedef ptrdiff_t difference_type;
        typedef _Tp *pointer;
        typedef const _Tp *const_pointer;
        typedef _Tp & reference;
        typedef const _Tp & const_reference;
        typedef _Tp value_type;

        template < typename _Tp1 > struct rebind {
                typedef pagealigned_allocator < _Tp1 > other;
        };

        pagealigned_allocator()throw() { }
        pagealigned_allocator(const pagealigned_allocator &) throw() { }
        ~pagealigned_allocator()throw() { }

        template < typename _Tp1 >
        pagealigned_allocator(const pagealigned_allocator < _Tp1 > &) throw() { }

        pointer address(reference __x) const { return &__x; }
        const_pointer address(const_reference __x) const { return &__x; }
        // NB: __n is permitted to be 0.  The C++ standard says nothing
        // about what the return value is when __n == 0.
        pointer allocate(size_type __n, const void * = 0) {
                return static_cast <
                    _Tp * >(malloc_pagealigned(__n * sizeof(_Tp)));
        }

        // __p is not permitted to be a null pointer.
        void deallocate(pointer __p, size_type __n) {
                free_pagealigned(__p, __n * sizeof(_Tp));
        }

        size_type max_size()const throw() {
                return size_t(-1) / sizeof(_Tp);
        }

        void construct(pointer __p, const _Tp & __val) {
                ::new(__p) _Tp(__val);
        }

        void destroy(pointer __p) { __p->~_Tp(); }
};
#endif
#endif

#ifdef __cplusplus
std::vector<std::pair<cxx_mpz, int> > trial_division(cxx_mpz const& n0, unsigned long B, cxx_mpz & cofactor);
cxx_mpz mpz_from_expression(const char *);
#endif

#endif	/* CADO_UTILS_MISC_H */
