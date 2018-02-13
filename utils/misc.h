#ifndef CADO_UTILS_MISC_H_
#define CADO_UTILS_MISC_H_

#include <stddef.h>
#include <stdint.h>
#include <stdarg.h>
#include <limits.h>
#include <errno.h>
#include <string.h>
#include <gmp.h>
#include "macros.h"
#include "portability.h"

/* we prefer GMP 5 or later, but the history of the why and how seems
 * lost. It seems that the late GMP-4.3 versions are fine, and the few
 * missing functions are provided as fallbacks */
#if !GMP_VERSION_ATLEAST(4,3,0)
#error "GNU MP 4.3.0 (at least) is required to compile CADO-NFS"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define UMAX(A) (0xffffffffffffffffULL >>((8-sizeof(A))<<3))
#define SMAX(A) (0x7fffffffffffffffLL  >>((8-sizeof(A))<<3))
#define UMIN(A) (0)
#define SMIN(A) (1ULL<< ((sizeof(A)<<3)-1))

/* uintmax_t is guaranteed to be larger or equal to uint64_t */
#define strtouint64(nptr,endptr,base) (uint64_t) strtoumax(nptr,endptr,base)

static inline void* pointer_arith(void * a, ptrdiff_t q) {
    return (void*)(((char*)a)+q);
}
static inline const void* pointer_arith_const(const void * a, ptrdiff_t q) {
    return (const void*)(((const char*)a)+q);
}

/* MinGW's string.h does not declare a prototype for strdup if __STRICT_ANSI__
   is defined */
#if !defined(HAVE_STRDUP) || (defined(__MINGW32__) && defined(__STRICT_ANSI__))
char * strdup(const char *s);
#endif
#ifndef HAVE_STRNDUP
char * strndup(const char * a, size_t n);
#endif

#ifndef HAVE_STRLCPY
size_t strlcpy(char *dst, const char *src, size_t size) ATTRIBUTE((__warn_unused_result__));
#endif
#ifndef HAVE_STRLCAT
size_t strlcat(char *dst, const char *src, size_t size) ATTRIBUTE((__warn_unused_result__));
#endif

/* strtoul(), but with const char ** for second argument.
   Otherwise it's not possible to do, e.g., strtoul(p, &p, 10) when p is
   of type const char *
*/
extern unsigned long int strtoul_const(const char *nptr, const char **endptr, const int base);
extern unsigned long long int strtoull_const(const char *nptr, const char **endptr, const int base);

extern char * derived_filename(const char * prefix, const char * what, const char * ext);
extern int has_suffix(const char * path, const char * sfx);

extern char ** filelist_from_file(const char * basepath, const char * filename,
                                  int typ);
extern void filelist_clear(char ** filelist);

long get_arg_max(void);
extern int mkdir_with_parents(const char * dir, int fatal);

extern char * path_resolve(const char * progname, char * resolved);

/* k must be a power of 2. Returns the smallest multiple of k greater
 * than or equal to n
 */
static inline unsigned long
next_multiple_of_powerof2(unsigned long n, unsigned long k)
{
    ASSERT((k & (k-1)) == 0);
    return ((n-1)|(k-1)) + 1;
}
static inline unsigned long integer_sqrt(unsigned long a)
{
    /* returns 2 for a==3, otherwise returns floor(sqrt(a)) */
    for(unsigned long x = a, y = 1, z; ; x = y, y = z) {
        z = (y + a/y) / 2;
        if (z == x || z == y) return y;
    }
}

/* Best X86 medium memcpy with pointers & size length already cache
   lines aligned on modern X86, so dst/src/lg & 0x3F = 0 */
static inline void aligned_medium_memcpy(void *dst, void *src, size_t lg) {
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
  size_t lg_8bytes = lg >> 3;
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
static inline int cado_clzll(unsigned long long x)
{
#if ULONGLONG_BITS == 64
        static const int t[4] = { 2, 1, 0, 0 };
        int a = 0;
        int res;
        if (x >> 32) { a += 32; x >>= 32; }
        if (x >> 16) { a += 16; x >>= 16; }
        if (x >>  8) { a +=  8; x >>=  8; }
        if (x >>  4) { a +=  4; x >>=  4; }
        if (x >>  2) { a +=  2; x >>=  2; }
        res = 64 - 2 - a + t[x];
        return res;
#else
#error "cado_clzll might be wrong here"
#endif
}

static inline int cado_clzl(unsigned long x)
{
        static const int t[4] = { 2, 1, 0, 0 };
        int a = 0;
        int res;
#if ULONG_BITS == 64
        if (x >> 32) { a += 32; x >>= 32; }
#endif
        if (x >> 16) { a += 16; x >>= 16; }
        if (x >>  8) { a +=  8; x >>=  8; }
        if (x >>  4) { a +=  4; x >>=  4; }
        if (x >>  2) { a +=  2; x >>=  2; }
        res = GMP_LIMB_BITS - 2 - a + t[x];
        return res;
}

static inline int cado_clz(unsigned int x)
{
        static const int t[4] = { 2, 1, 0, 0 };
        int a = 0;
        int res;
        if (x >> 16) { a += 16; x >>= 16; }
        if (x >>  8) { a +=  8; x >>=  8; }
        if (x >>  4) { a +=  4; x >>=  4; }
        if (x >>  2) { a +=  2; x >>=  2; }
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
static inline int cado_ctzll(unsigned long long x)
{
  return (ULONGLONG_BITS - 1) - cado_clzll(x & - x);
}
static inline int cado_ctzl(unsigned long x)
{
  return (ULONG_BITS - 1) - cado_clzl(x & - x);
}
static inline int cado_ctz(unsigned int x)
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
    x ^= x >> 32;
    x ^= x >> 16;
    x ^= x >> 8;
    x ^= x >> 4;
    int t[16] = { 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0 };
    return t[x&15];
#else
#error "fix cado_parityll here"
#endif
}
#endif

#if ULONGLONG_BITS == 64
static inline int cado_ctz64(uint64_t x) { return cado_ctzll(x); }
static inline int cado_clz64(uint64_t x) { return cado_clzll(x); }
static inline int cado_parity64(uint64_t x) { return cado_parityll(x); }
#else
#error "need proper equivalents for cado_ctz64 & friends"
#endif
/*}}}*/
const char *size_disp_fine(size_t s, char buf[16], double cutoff);
const char *size_disp(size_t s, char buf[16]);

#ifdef __cplusplus
}
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

#endif	/* CADO_UTILS_MISC_H_ */
