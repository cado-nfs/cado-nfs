#include <limits.h>
#ifndef TYPECAST_H_
#define TYPECAST_H_

#include "macros.h"

/* Typecasting inline functions that assert that cast values are within 
   range of the cast-to type.
   Format is cast_<source type>_<target type>()
   where each of <source type> and <target type> are one of
   ulong	for unsigned long int
   uint		for unsigned int
   long		for long int
   int		for int
   ushort	for unsigned short int
   short	for short int
   uchar	for unsigned char
   char		for char
   size		for size_t
   
   For each function, there is also a
   cast_<source type>_<target type>_fast()
   variant which asserts range validity only if WANT_ASSERT_CONVERSION is
   defined. Thus the non-fast variant should be used in code that is not
   speed critical, and the fast variant in code that is.
*/

#ifdef WANT_ASSERT_CONVERSION
#define ASSERT_CONVERSION(x) assert(x)
#else
#define ASSERT_CONVERSION(x)
#endif


#ifdef __cplusplus
extern "C" {
#endif

static inline unsigned int
cast_ulong_uint(unsigned long int x)
{
#if ULONG_MAX > UINT_MAX
    ASSERT(x <= UINT_MAX);
#endif
  return (unsigned int) x;
}

static inline unsigned int
cast_ulong_uint_fast(unsigned long int x)
{
#if ULONG_MAX > UINT_MAX
    ASSERT_CONVERSION(x <= UINT_MAX);
#endif
  return (unsigned int) x;
}

static inline long int
cast_ulong_long(unsigned long int x)
{
#if 0 < LONG_MIN
    ASSERT(x >= LONG_MIN);
#endif
#if ULONG_MAX > LONG_MAX
    ASSERT(x <= LONG_MAX);
#endif
  return (long int) x;
}

static inline long int
cast_ulong_long_fast(unsigned long int x)
{
#if 0 < LONG_MIN
    ASSERT_CONVERSION(x >= LONG_MIN);
#endif
#if ULONG_MAX > LONG_MAX
    ASSERT_CONVERSION(x <= LONG_MAX);
#endif
  return (long int) x;
}

static inline int
cast_ulong_int(unsigned long int x)
{
#if 0 < INT_MIN
    ASSERT(x >= INT_MIN);
#endif
#if ULONG_MAX > INT_MAX
    ASSERT(x <= INT_MAX);
#endif
  return (int) x;
}

static inline int
cast_ulong_int_fast(unsigned long int x)
{
#if 0 < INT_MIN
    ASSERT_CONVERSION(x >= INT_MIN);
#endif
#if ULONG_MAX > INT_MAX
    ASSERT_CONVERSION(x <= INT_MAX);
#endif
  return (int) x;
}

static inline unsigned short int
cast_ulong_ushort(unsigned long int x)
{
#if ULONG_MAX > USHRT_MAX
    ASSERT(x <= USHRT_MAX);
#endif
  return (unsigned short int) x;
}

static inline unsigned short int
cast_ulong_ushort_fast(unsigned long int x)
{
#if ULONG_MAX > USHRT_MAX
    ASSERT_CONVERSION(x <= USHRT_MAX);
#endif
  return (unsigned short int) x;
}

static inline short int
cast_ulong_short(unsigned long int x)
{
#if 0 < SHRT_MIN
    ASSERT(x >= SHRT_MIN);
#endif
#if ULONG_MAX > SHRT_MAX
    ASSERT(x <= SHRT_MAX);
#endif
  return (short int) x;
}

static inline short int
cast_ulong_short_fast(unsigned long int x)
{
#if 0 < SHRT_MIN
    ASSERT_CONVERSION(x >= SHRT_MIN);
#endif
#if ULONG_MAX > SHRT_MAX
    ASSERT_CONVERSION(x <= SHRT_MAX);
#endif
  return (short int) x;
}

static inline unsigned char
cast_ulong_uchar(unsigned long int x)
{
#if ULONG_MAX > UCHAR_MAX
    ASSERT(x <= UCHAR_MAX);
#endif
  return (unsigned char) x;
}

static inline unsigned char
cast_ulong_uchar_fast(unsigned long int x)
{
#if ULONG_MAX > UCHAR_MAX
    ASSERT_CONVERSION(x <= UCHAR_MAX);
#endif
  return (unsigned char) x;
}

static inline char
cast_ulong_char(unsigned long int x)
{
#if 0 < CHAR_MIN
    ASSERT(x >= CHAR_MIN);
#endif
#if ULONG_MAX > CHAR_MAX
    ASSERT(x <= CHAR_MAX);
#endif
  return (char) x;
}

static inline char
cast_ulong_char_fast(unsigned long int x)
{
#if 0 < CHAR_MIN
    ASSERT_CONVERSION(x >= CHAR_MIN);
#endif
#if ULONG_MAX > CHAR_MAX
    ASSERT_CONVERSION(x <= CHAR_MAX);
#endif
  return (char) x;
}

static inline size_t
cast_ulong_size(unsigned long int x)
{
#if ULONG_MAX > SIZE_MAX
    ASSERT(x <= SIZE_MAX);
#endif
  return (size_t) x;
}

static inline size_t
cast_ulong_size_fast(unsigned long int x)
{
#if ULONG_MAX > SIZE_MAX
    ASSERT_CONVERSION(x <= SIZE_MAX);
#endif
  return (size_t) x;
}

static inline unsigned long int
cast_uint_ulong(unsigned int x)
{
#if UINT_MAX > ULONG_MAX
    ASSERT(x <= ULONG_MAX);
#endif
  return (unsigned long int) x;
}

static inline unsigned long int
cast_uint_ulong_fast(unsigned int x)
{
#if UINT_MAX > ULONG_MAX
    ASSERT_CONVERSION(x <= ULONG_MAX);
#endif
  return (unsigned long int) x;
}

static inline long int
cast_uint_long(unsigned int x)
{
#if 0 < LONG_MIN
    ASSERT(x >= LONG_MIN);
#endif
#if UINT_MAX > LONG_MAX
    ASSERT(x <= LONG_MAX);
#endif
  return (long int) x;
}

static inline long int
cast_uint_long_fast(unsigned int x)
{
#if 0 < LONG_MIN
    ASSERT_CONVERSION(x >= LONG_MIN);
#endif
#if UINT_MAX > LONG_MAX
    ASSERT_CONVERSION(x <= LONG_MAX);
#endif
  return (long int) x;
}

static inline int
cast_uint_int(unsigned int x)
{
#if 0 < INT_MIN
    ASSERT(x >= INT_MIN);
#endif
#if UINT_MAX > INT_MAX
    ASSERT(x <= INT_MAX);
#endif
  return (int) x;
}

static inline int
cast_uint_int_fast(unsigned int x)
{
#if 0 < INT_MIN
    ASSERT_CONVERSION(x >= INT_MIN);
#endif
#if UINT_MAX > INT_MAX
    ASSERT_CONVERSION(x <= INT_MAX);
#endif
  return (int) x;
}

static inline unsigned short int
cast_uint_ushort(unsigned int x)
{
#if UINT_MAX > USHRT_MAX
    ASSERT(x <= USHRT_MAX);
#endif
  return (unsigned short int) x;
}

static inline unsigned short int
cast_uint_ushort_fast(unsigned int x)
{
#if UINT_MAX > USHRT_MAX
    ASSERT_CONVERSION(x <= USHRT_MAX);
#endif
  return (unsigned short int) x;
}

static inline short int
cast_uint_short(unsigned int x)
{
#if 0 < SHRT_MIN
    ASSERT(x >= SHRT_MIN);
#endif
#if UINT_MAX > SHRT_MAX
    ASSERT(x <= SHRT_MAX);
#endif
  return (short int) x;
}

static inline short int
cast_uint_short_fast(unsigned int x)
{
#if 0 < SHRT_MIN
    ASSERT_CONVERSION(x >= SHRT_MIN);
#endif
#if UINT_MAX > SHRT_MAX
    ASSERT_CONVERSION(x <= SHRT_MAX);
#endif
  return (short int) x;
}

static inline unsigned char
cast_uint_uchar(unsigned int x)
{
#if UINT_MAX > UCHAR_MAX
    ASSERT(x <= UCHAR_MAX);
#endif
  return (unsigned char) x;
}

static inline unsigned char
cast_uint_uchar_fast(unsigned int x)
{
#if UINT_MAX > UCHAR_MAX
    ASSERT_CONVERSION(x <= UCHAR_MAX);
#endif
  return (unsigned char) x;
}

static inline char
cast_uint_char(unsigned int x)
{
#if 0 < CHAR_MIN
    ASSERT(x >= CHAR_MIN);
#endif
#if UINT_MAX > CHAR_MAX
    ASSERT(x <= CHAR_MAX);
#endif
  return (char) x;
}

static inline char
cast_uint_char_fast(unsigned int x)
{
#if 0 < CHAR_MIN
    ASSERT_CONVERSION(x >= CHAR_MIN);
#endif
#if UINT_MAX > CHAR_MAX
    ASSERT_CONVERSION(x <= CHAR_MAX);
#endif
  return (char) x;
}

static inline size_t
cast_uint_size(unsigned int x)
{
#if UINT_MAX > SIZE_MAX
    ASSERT(x <= SIZE_MAX);
#endif
  return (size_t) x;
}

static inline size_t
cast_uint_size_fast(unsigned int x)
{
#if UINT_MAX > SIZE_MAX
    ASSERT_CONVERSION(x <= SIZE_MAX);
#endif
  return (size_t) x;
}

static inline unsigned long int
cast_long_ulong(long int x)
{
#if LONG_MIN < 0
    ASSERT(x >= 0);
#endif
#if LONG_MAX > ULONG_MAX
    ASSERT(x <= ULONG_MAX);
#endif
  return (unsigned long int) x;
}

static inline unsigned long int
cast_long_ulong_fast(long int x)
{
#if LONG_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if LONG_MAX > ULONG_MAX
    ASSERT_CONVERSION(x <= ULONG_MAX);
#endif
  return (unsigned long int) x;
}

static inline unsigned int
cast_long_uint(long int x)
{
#if LONG_MIN < 0
    ASSERT(x >= 0);
#endif
#if LONG_MAX > UINT_MAX
    ASSERT(x <= UINT_MAX);
#endif
  return (unsigned int) x;
}

static inline unsigned int
cast_long_uint_fast(long int x)
{
#if LONG_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if LONG_MAX > UINT_MAX
    ASSERT_CONVERSION(x <= UINT_MAX);
#endif
  return (unsigned int) x;
}

static inline int
cast_long_int(long int x)
{
#if LONG_MIN < INT_MIN
    ASSERT(x >= INT_MIN);
#endif
#if LONG_MAX > INT_MAX
    ASSERT(x <= INT_MAX);
#endif
  return (int) x;
}

static inline int
cast_long_int_fast(long int x)
{
#if LONG_MIN < INT_MIN
    ASSERT_CONVERSION(x >= INT_MIN);
#endif
#if LONG_MAX > INT_MAX
    ASSERT_CONVERSION(x <= INT_MAX);
#endif
  return (int) x;
}

static inline unsigned short int
cast_long_ushort(long int x)
{
#if LONG_MIN < 0
    ASSERT(x >= 0);
#endif
#if LONG_MAX > USHRT_MAX
    ASSERT(x <= USHRT_MAX);
#endif
  return (unsigned short int) x;
}

static inline unsigned short int
cast_long_ushort_fast(long int x)
{
#if LONG_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if LONG_MAX > USHRT_MAX
    ASSERT_CONVERSION(x <= USHRT_MAX);
#endif
  return (unsigned short int) x;
}

static inline short int
cast_long_short(long int x)
{
#if LONG_MIN < SHRT_MIN
    ASSERT(x >= SHRT_MIN);
#endif
#if LONG_MAX > SHRT_MAX
    ASSERT(x <= SHRT_MAX);
#endif
  return (short int) x;
}

static inline short int
cast_long_short_fast(long int x)
{
#if LONG_MIN < SHRT_MIN
    ASSERT_CONVERSION(x >= SHRT_MIN);
#endif
#if LONG_MAX > SHRT_MAX
    ASSERT_CONVERSION(x <= SHRT_MAX);
#endif
  return (short int) x;
}

static inline unsigned char
cast_long_uchar(long int x)
{
#if LONG_MIN < 0
    ASSERT(x >= 0);
#endif
#if LONG_MAX > UCHAR_MAX
    ASSERT(x <= UCHAR_MAX);
#endif
  return (unsigned char) x;
}

static inline unsigned char
cast_long_uchar_fast(long int x)
{
#if LONG_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if LONG_MAX > UCHAR_MAX
    ASSERT_CONVERSION(x <= UCHAR_MAX);
#endif
  return (unsigned char) x;
}

static inline char
cast_long_char(long int x)
{
#if LONG_MIN < CHAR_MIN
    ASSERT(x >= CHAR_MIN);
#endif
#if LONG_MAX > CHAR_MAX
    ASSERT(x <= CHAR_MAX);
#endif
  return (char) x;
}

static inline char
cast_long_char_fast(long int x)
{
#if LONG_MIN < CHAR_MIN
    ASSERT_CONVERSION(x >= CHAR_MIN);
#endif
#if LONG_MAX > CHAR_MAX
    ASSERT_CONVERSION(x <= CHAR_MAX);
#endif
  return (char) x;
}

static inline size_t
cast_long_size(long int x)
{
#if LONG_MIN < 0
    ASSERT(x >= 0);
#endif
#if LONG_MAX > SIZE_MAX
    ASSERT(x <= SIZE_MAX);
#endif
  return (size_t) x;
}

static inline size_t
cast_long_size_fast(long int x)
{
#if LONG_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if LONG_MAX > SIZE_MAX
    ASSERT_CONVERSION(x <= SIZE_MAX);
#endif
  return (size_t) x;
}

static inline unsigned long int
cast_int_ulong(int x)
{
#if INT_MIN < 0
    ASSERT(x >= 0);
#endif
#if INT_MAX > ULONG_MAX
    ASSERT(x <= ULONG_MAX);
#endif
  return (unsigned long int) x;
}

static inline unsigned long int
cast_int_ulong_fast(int x)
{
#if INT_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if INT_MAX > ULONG_MAX
    ASSERT_CONVERSION(x <= ULONG_MAX);
#endif
  return (unsigned long int) x;
}

static inline unsigned int
cast_int_uint(int x)
{
#if INT_MIN < 0
    ASSERT(x >= 0);
#endif
#if INT_MAX > UINT_MAX
    ASSERT(x <= UINT_MAX);
#endif
  return (unsigned int) x;
}

static inline unsigned int
cast_int_uint_fast(int x)
{
#if INT_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if INT_MAX > UINT_MAX
    ASSERT_CONVERSION(x <= UINT_MAX);
#endif
  return (unsigned int) x;
}

static inline long int
cast_int_long(int x)
{
#if INT_MIN < LONG_MIN
    ASSERT(x >= LONG_MIN);
#endif
#if INT_MAX > LONG_MAX
    ASSERT(x <= LONG_MAX);
#endif
  return (long int) x;
}

static inline long int
cast_int_long_fast(int x)
{
#if INT_MIN < LONG_MIN
    ASSERT_CONVERSION(x >= LONG_MIN);
#endif
#if INT_MAX > LONG_MAX
    ASSERT_CONVERSION(x <= LONG_MAX);
#endif
  return (long int) x;
}

static inline unsigned short int
cast_int_ushort(int x)
{
#if INT_MIN < 0
    ASSERT(x >= 0);
#endif
#if INT_MAX > USHRT_MAX
    ASSERT(x <= USHRT_MAX);
#endif
  return (unsigned short int) x;
}

static inline unsigned short int
cast_int_ushort_fast(int x)
{
#if INT_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if INT_MAX > USHRT_MAX
    ASSERT_CONVERSION(x <= USHRT_MAX);
#endif
  return (unsigned short int) x;
}

static inline short int
cast_int_short(int x)
{
#if INT_MIN < SHRT_MIN
    ASSERT(x >= SHRT_MIN);
#endif
#if INT_MAX > SHRT_MAX
    ASSERT(x <= SHRT_MAX);
#endif
  return (short int) x;
}

static inline short int
cast_int_short_fast(int x)
{
#if INT_MIN < SHRT_MIN
    ASSERT_CONVERSION(x >= SHRT_MIN);
#endif
#if INT_MAX > SHRT_MAX
    ASSERT_CONVERSION(x <= SHRT_MAX);
#endif
  return (short int) x;
}

static inline unsigned char
cast_int_uchar(int x)
{
#if INT_MIN < 0
    ASSERT(x >= 0);
#endif
#if INT_MAX > UCHAR_MAX
    ASSERT(x <= UCHAR_MAX);
#endif
  return (unsigned char) x;
}

static inline unsigned char
cast_int_uchar_fast(int x)
{
#if INT_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if INT_MAX > UCHAR_MAX
    ASSERT_CONVERSION(x <= UCHAR_MAX);
#endif
  return (unsigned char) x;
}

static inline char
cast_int_char(int x)
{
#if INT_MIN < CHAR_MIN
    ASSERT(x >= CHAR_MIN);
#endif
#if INT_MAX > CHAR_MAX
    ASSERT(x <= CHAR_MAX);
#endif
  return (char) x;
}

static inline char
cast_int_char_fast(int x)
{
#if INT_MIN < CHAR_MIN
    ASSERT_CONVERSION(x >= CHAR_MIN);
#endif
#if INT_MAX > CHAR_MAX
    ASSERT_CONVERSION(x <= CHAR_MAX);
#endif
  return (char) x;
}

static inline size_t
cast_int_size(int x)
{
#if INT_MIN < 0
    ASSERT(x >= 0);
#endif
#if INT_MAX > SIZE_MAX
    ASSERT(x <= SIZE_MAX);
#endif
  return (size_t) x;
}

static inline size_t
cast_int_size_fast(int x)
{
#if INT_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if INT_MAX > SIZE_MAX
    ASSERT_CONVERSION(x <= SIZE_MAX);
#endif
  return (size_t) x;
}

static inline unsigned long int
cast_ushort_ulong(unsigned short int x)
{
#if USHRT_MAX > ULONG_MAX
    ASSERT(x <= ULONG_MAX);
#endif
  return (unsigned long int) x;
}

static inline unsigned long int
cast_ushort_ulong_fast(unsigned short int x)
{
#if USHRT_MAX > ULONG_MAX
    ASSERT_CONVERSION(x <= ULONG_MAX);
#endif
  return (unsigned long int) x;
}

static inline unsigned int
cast_ushort_uint(unsigned short int x)
{
#if USHRT_MAX > UINT_MAX
    ASSERT(x <= UINT_MAX);
#endif
  return (unsigned int) x;
}

static inline unsigned int
cast_ushort_uint_fast(unsigned short int x)
{
#if USHRT_MAX > UINT_MAX
    ASSERT_CONVERSION(x <= UINT_MAX);
#endif
  return (unsigned int) x;
}

static inline long int
cast_ushort_long(unsigned short int x)
{
#if 0 < LONG_MIN
    ASSERT(x >= LONG_MIN);
#endif
#if USHRT_MAX > LONG_MAX
    ASSERT(x <= LONG_MAX);
#endif
  return (long int) x;
}

static inline long int
cast_ushort_long_fast(unsigned short int x)
{
#if 0 < LONG_MIN
    ASSERT_CONVERSION(x >= LONG_MIN);
#endif
#if USHRT_MAX > LONG_MAX
    ASSERT_CONVERSION(x <= LONG_MAX);
#endif
  return (long int) x;
}

static inline int
cast_ushort_int(unsigned short int x)
{
#if 0 < INT_MIN
    ASSERT(x >= INT_MIN);
#endif
#if USHRT_MAX > INT_MAX
    ASSERT(x <= INT_MAX);
#endif
  return (int) x;
}

static inline int
cast_ushort_int_fast(unsigned short int x)
{
#if 0 < INT_MIN
    ASSERT_CONVERSION(x >= INT_MIN);
#endif
#if USHRT_MAX > INT_MAX
    ASSERT_CONVERSION(x <= INT_MAX);
#endif
  return (int) x;
}

static inline short int
cast_ushort_short(unsigned short int x)
{
#if 0 < SHRT_MIN
    ASSERT(x >= SHRT_MIN);
#endif
#if USHRT_MAX > SHRT_MAX
    ASSERT(x <= SHRT_MAX);
#endif
  return (short int) x;
}

static inline short int
cast_ushort_short_fast(unsigned short int x)
{
#if 0 < SHRT_MIN
    ASSERT_CONVERSION(x >= SHRT_MIN);
#endif
#if USHRT_MAX > SHRT_MAX
    ASSERT_CONVERSION(x <= SHRT_MAX);
#endif
  return (short int) x;
}

static inline unsigned char
cast_ushort_uchar(unsigned short int x)
{
#if USHRT_MAX > UCHAR_MAX
    ASSERT(x <= UCHAR_MAX);
#endif
  return (unsigned char) x;
}

static inline unsigned char
cast_ushort_uchar_fast(unsigned short int x)
{
#if USHRT_MAX > UCHAR_MAX
    ASSERT_CONVERSION(x <= UCHAR_MAX);
#endif
  return (unsigned char) x;
}

static inline char
cast_ushort_char(unsigned short int x)
{
#if 0 < CHAR_MIN
    ASSERT(x >= CHAR_MIN);
#endif
#if USHRT_MAX > CHAR_MAX
    ASSERT(x <= CHAR_MAX);
#endif
  return (char) x;
}

static inline char
cast_ushort_char_fast(unsigned short int x)
{
#if 0 < CHAR_MIN
    ASSERT_CONVERSION(x >= CHAR_MIN);
#endif
#if USHRT_MAX > CHAR_MAX
    ASSERT_CONVERSION(x <= CHAR_MAX);
#endif
  return (char) x;
}

static inline size_t
cast_ushort_size(unsigned short int x)
{
#if USHRT_MAX > SIZE_MAX
    ASSERT(x <= SIZE_MAX);
#endif
  return (size_t) x;
}

static inline size_t
cast_ushort_size_fast(unsigned short int x)
{
#if USHRT_MAX > SIZE_MAX
    ASSERT_CONVERSION(x <= SIZE_MAX);
#endif
  return (size_t) x;
}

static inline unsigned long int
cast_short_ulong(short int x)
{
#if SHRT_MIN < 0
    ASSERT(x >= 0);
#endif
#if SHRT_MAX > ULONG_MAX
    ASSERT(x <= ULONG_MAX);
#endif
  return (unsigned long int) x;
}

static inline unsigned long int
cast_short_ulong_fast(short int x)
{
#if SHRT_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if SHRT_MAX > ULONG_MAX
    ASSERT_CONVERSION(x <= ULONG_MAX);
#endif
  return (unsigned long int) x;
}

static inline unsigned int
cast_short_uint(short int x)
{
#if SHRT_MIN < 0
    ASSERT(x >= 0);
#endif
#if SHRT_MAX > UINT_MAX
    ASSERT(x <= UINT_MAX);
#endif
  return (unsigned int) x;
}

static inline unsigned int
cast_short_uint_fast(short int x)
{
#if SHRT_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if SHRT_MAX > UINT_MAX
    ASSERT_CONVERSION(x <= UINT_MAX);
#endif
  return (unsigned int) x;
}

static inline long int
cast_short_long(short int x)
{
#if SHRT_MIN < LONG_MIN
    ASSERT(x >= LONG_MIN);
#endif
#if SHRT_MAX > LONG_MAX
    ASSERT(x <= LONG_MAX);
#endif
  return (long int) x;
}

static inline long int
cast_short_long_fast(short int x)
{
#if SHRT_MIN < LONG_MIN
    ASSERT_CONVERSION(x >= LONG_MIN);
#endif
#if SHRT_MAX > LONG_MAX
    ASSERT_CONVERSION(x <= LONG_MAX);
#endif
  return (long int) x;
}

static inline int
cast_short_int(short int x)
{
#if SHRT_MIN < INT_MIN
    ASSERT(x >= INT_MIN);
#endif
#if SHRT_MAX > INT_MAX
    ASSERT(x <= INT_MAX);
#endif
  return (int) x;
}

static inline int
cast_short_int_fast(short int x)
{
#if SHRT_MIN < INT_MIN
    ASSERT_CONVERSION(x >= INT_MIN);
#endif
#if SHRT_MAX > INT_MAX
    ASSERT_CONVERSION(x <= INT_MAX);
#endif
  return (int) x;
}

static inline unsigned short int
cast_short_ushort(short int x)
{
#if SHRT_MIN < 0
    ASSERT(x >= 0);
#endif
#if SHRT_MAX > USHRT_MAX
    ASSERT(x <= USHRT_MAX);
#endif
  return (unsigned short int) x;
}

static inline unsigned short int
cast_short_ushort_fast(short int x)
{
#if SHRT_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if SHRT_MAX > USHRT_MAX
    ASSERT_CONVERSION(x <= USHRT_MAX);
#endif
  return (unsigned short int) x;
}

static inline unsigned char
cast_short_uchar(short int x)
{
#if SHRT_MIN < 0
    ASSERT(x >= 0);
#endif
#if SHRT_MAX > UCHAR_MAX
    ASSERT(x <= UCHAR_MAX);
#endif
  return (unsigned char) x;
}

static inline unsigned char
cast_short_uchar_fast(short int x)
{
#if SHRT_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if SHRT_MAX > UCHAR_MAX
    ASSERT_CONVERSION(x <= UCHAR_MAX);
#endif
  return (unsigned char) x;
}

static inline char
cast_short_char(short int x)
{
#if SHRT_MIN < CHAR_MIN
    ASSERT(x >= CHAR_MIN);
#endif
#if SHRT_MAX > CHAR_MAX
    ASSERT(x <= CHAR_MAX);
#endif
  return (char) x;
}

static inline char
cast_short_char_fast(short int x)
{
#if SHRT_MIN < CHAR_MIN
    ASSERT_CONVERSION(x >= CHAR_MIN);
#endif
#if SHRT_MAX > CHAR_MAX
    ASSERT_CONVERSION(x <= CHAR_MAX);
#endif
  return (char) x;
}

static inline size_t
cast_short_size(short int x)
{
#if SHRT_MIN < 0
    ASSERT(x >= 0);
#endif
#if SHRT_MAX > SIZE_MAX
    ASSERT(x <= SIZE_MAX);
#endif
  return (size_t) x;
}

static inline size_t
cast_short_size_fast(short int x)
{
#if SHRT_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if SHRT_MAX > SIZE_MAX
    ASSERT_CONVERSION(x <= SIZE_MAX);
#endif
  return (size_t) x;
}

static inline unsigned long int
cast_uchar_ulong(unsigned char x)
{
#if UCHAR_MAX > ULONG_MAX
    ASSERT(x <= ULONG_MAX);
#endif
  return (unsigned long int) x;
}

static inline unsigned long int
cast_uchar_ulong_fast(unsigned char x)
{
#if UCHAR_MAX > ULONG_MAX
    ASSERT_CONVERSION(x <= ULONG_MAX);
#endif
  return (unsigned long int) x;
}

static inline unsigned int
cast_uchar_uint(unsigned char x)
{
#if UCHAR_MAX > UINT_MAX
    ASSERT(x <= UINT_MAX);
#endif
  return (unsigned int) x;
}

static inline unsigned int
cast_uchar_uint_fast(unsigned char x)
{
#if UCHAR_MAX > UINT_MAX
    ASSERT_CONVERSION(x <= UINT_MAX);
#endif
  return (unsigned int) x;
}

static inline long int
cast_uchar_long(unsigned char x)
{
#if 0 < LONG_MIN
    ASSERT(x >= LONG_MIN);
#endif
#if UCHAR_MAX > LONG_MAX
    ASSERT(x <= LONG_MAX);
#endif
  return (long int) x;
}

static inline long int
cast_uchar_long_fast(unsigned char x)
{
#if 0 < LONG_MIN
    ASSERT_CONVERSION(x >= LONG_MIN);
#endif
#if UCHAR_MAX > LONG_MAX
    ASSERT_CONVERSION(x <= LONG_MAX);
#endif
  return (long int) x;
}

static inline int
cast_uchar_int(unsigned char x)
{
#if 0 < INT_MIN
    ASSERT(x >= INT_MIN);
#endif
#if UCHAR_MAX > INT_MAX
    ASSERT(x <= INT_MAX);
#endif
  return (int) x;
}

static inline int
cast_uchar_int_fast(unsigned char x)
{
#if 0 < INT_MIN
    ASSERT_CONVERSION(x >= INT_MIN);
#endif
#if UCHAR_MAX > INT_MAX
    ASSERT_CONVERSION(x <= INT_MAX);
#endif
  return (int) x;
}

static inline unsigned short int
cast_uchar_ushort(unsigned char x)
{
#if UCHAR_MAX > USHRT_MAX
    ASSERT(x <= USHRT_MAX);
#endif
  return (unsigned short int) x;
}

static inline unsigned short int
cast_uchar_ushort_fast(unsigned char x)
{
#if UCHAR_MAX > USHRT_MAX
    ASSERT_CONVERSION(x <= USHRT_MAX);
#endif
  return (unsigned short int) x;
}

static inline short int
cast_uchar_short(unsigned char x)
{
#if 0 < SHRT_MIN
    ASSERT(x >= SHRT_MIN);
#endif
#if UCHAR_MAX > SHRT_MAX
    ASSERT(x <= SHRT_MAX);
#endif
  return (short int) x;
}

static inline short int
cast_uchar_short_fast(unsigned char x)
{
#if 0 < SHRT_MIN
    ASSERT_CONVERSION(x >= SHRT_MIN);
#endif
#if UCHAR_MAX > SHRT_MAX
    ASSERT_CONVERSION(x <= SHRT_MAX);
#endif
  return (short int) x;
}

static inline char
cast_uchar_char(unsigned char x)
{
#if 0 < CHAR_MIN
    ASSERT(x >= CHAR_MIN);
#endif
#if UCHAR_MAX > CHAR_MAX
    ASSERT(x <= CHAR_MAX);
#endif
  return (char) x;
}

static inline char
cast_uchar_char_fast(unsigned char x)
{
#if 0 < CHAR_MIN
    ASSERT_CONVERSION(x >= CHAR_MIN);
#endif
#if UCHAR_MAX > CHAR_MAX
    ASSERT_CONVERSION(x <= CHAR_MAX);
#endif
  return (char) x;
}

static inline size_t
cast_uchar_size(unsigned char x)
{
#if UCHAR_MAX > SIZE_MAX
    ASSERT(x <= SIZE_MAX);
#endif
  return (size_t) x;
}

static inline size_t
cast_uchar_size_fast(unsigned char x)
{
#if UCHAR_MAX > SIZE_MAX
    ASSERT_CONVERSION(x <= SIZE_MAX);
#endif
  return (size_t) x;
}

static inline unsigned long int
cast_char_ulong(char x)
{
#if CHAR_MIN < 0
    ASSERT(x >= 0);
#endif
#if CHAR_MAX > ULONG_MAX
    ASSERT(x <= ULONG_MAX);
#endif
  return (unsigned long int) x;
}

static inline unsigned long int
cast_char_ulong_fast(char x)
{
#if CHAR_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if CHAR_MAX > ULONG_MAX
    ASSERT_CONVERSION(x <= ULONG_MAX);
#endif
  return (unsigned long int) x;
}

static inline unsigned int
cast_char_uint(char x)
{
#if CHAR_MIN < 0
    ASSERT(x >= 0);
#endif
#if CHAR_MAX > UINT_MAX
    ASSERT(x <= UINT_MAX);
#endif
  return (unsigned int) x;
}

static inline unsigned int
cast_char_uint_fast(char x)
{
#if CHAR_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if CHAR_MAX > UINT_MAX
    ASSERT_CONVERSION(x <= UINT_MAX);
#endif
  return (unsigned int) x;
}

static inline long int
cast_char_long(char x)
{
#if CHAR_MIN < LONG_MIN
    ASSERT(x >= LONG_MIN);
#endif
#if CHAR_MAX > LONG_MAX
    ASSERT(x <= LONG_MAX);
#endif
  return (long int) x;
}

static inline long int
cast_char_long_fast(char x)
{
#if CHAR_MIN < LONG_MIN
    ASSERT_CONVERSION(x >= LONG_MIN);
#endif
#if CHAR_MAX > LONG_MAX
    ASSERT_CONVERSION(x <= LONG_MAX);
#endif
  return (long int) x;
}

static inline int
cast_char_int(char x)
{
#if CHAR_MIN < INT_MIN
    ASSERT(x >= INT_MIN);
#endif
#if CHAR_MAX > INT_MAX
    ASSERT(x <= INT_MAX);
#endif
  return (int) x;
}

static inline int
cast_char_int_fast(char x)
{
#if CHAR_MIN < INT_MIN
    ASSERT_CONVERSION(x >= INT_MIN);
#endif
#if CHAR_MAX > INT_MAX
    ASSERT_CONVERSION(x <= INT_MAX);
#endif
  return (int) x;
}

static inline unsigned short int
cast_char_ushort(char x)
{
#if CHAR_MIN < 0
    ASSERT(x >= 0);
#endif
#if CHAR_MAX > USHRT_MAX
    ASSERT(x <= USHRT_MAX);
#endif
  return (unsigned short int) x;
}

static inline unsigned short int
cast_char_ushort_fast(char x)
{
#if CHAR_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if CHAR_MAX > USHRT_MAX
    ASSERT_CONVERSION(x <= USHRT_MAX);
#endif
  return (unsigned short int) x;
}

static inline short int
cast_char_short(char x)
{
#if CHAR_MIN < SHRT_MIN
    ASSERT(x >= SHRT_MIN);
#endif
#if CHAR_MAX > SHRT_MAX
    ASSERT(x <= SHRT_MAX);
#endif
  return (short int) x;
}

static inline short int
cast_char_short_fast(char x)
{
#if CHAR_MIN < SHRT_MIN
    ASSERT_CONVERSION(x >= SHRT_MIN);
#endif
#if CHAR_MAX > SHRT_MAX
    ASSERT_CONVERSION(x <= SHRT_MAX);
#endif
  return (short int) x;
}

static inline unsigned char
cast_char_uchar(char x)
{
#if CHAR_MIN < 0
    ASSERT(x >= 0);
#endif
#if CHAR_MAX > UCHAR_MAX
    ASSERT(x <= UCHAR_MAX);
#endif
  return (unsigned char) x;
}

static inline unsigned char
cast_char_uchar_fast(char x)
{
#if CHAR_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if CHAR_MAX > UCHAR_MAX
    ASSERT_CONVERSION(x <= UCHAR_MAX);
#endif
  return (unsigned char) x;
}

static inline size_t
cast_char_size(char x)
{
#if CHAR_MIN < 0
    ASSERT(x >= 0);
#endif
#if CHAR_MAX > SIZE_MAX
    ASSERT(x <= SIZE_MAX);
#endif
  return (size_t) x;
}

static inline size_t
cast_char_size_fast(char x)
{
#if CHAR_MIN < 0
    ASSERT_CONVERSION(x >= 0);
#endif
#if CHAR_MAX > SIZE_MAX
    ASSERT_CONVERSION(x <= SIZE_MAX);
#endif
  return (size_t) x;
}

static inline unsigned long int
cast_size_ulong(size_t x)
{
#if SIZE_MAX > ULONG_MAX
    ASSERT(x <= ULONG_MAX);
#endif
  return (unsigned long int) x;
}

static inline unsigned long int
cast_size_ulong_fast(size_t x)
{
#if SIZE_MAX > ULONG_MAX
    ASSERT_CONVERSION(x <= ULONG_MAX);
#endif
  return (unsigned long int) x;
}

static inline unsigned int
cast_size_uint(size_t x)
{
#if SIZE_MAX > UINT_MAX
    ASSERT(x <= UINT_MAX);
#endif
  return (unsigned int) x;
}

static inline unsigned int
cast_size_uint_fast(size_t x)
{
#if SIZE_MAX > UINT_MAX
    ASSERT_CONVERSION(x <= UINT_MAX);
#endif
  return (unsigned int) x;
}

static inline long int
cast_size_long(size_t x)
{
#if 0 < LONG_MIN
    ASSERT(x >= LONG_MIN);
#endif
#if SIZE_MAX > LONG_MAX
    ASSERT(x <= LONG_MAX);
#endif
  return (long int) x;
}

static inline long int
cast_size_long_fast(size_t x)
{
#if 0 < LONG_MIN
    ASSERT_CONVERSION(x >= LONG_MIN);
#endif
#if SIZE_MAX > LONG_MAX
    ASSERT_CONVERSION(x <= LONG_MAX);
#endif
  return (long int) x;
}

static inline int
cast_size_int(size_t x)
{
#if 0 < INT_MIN
    ASSERT(x >= INT_MIN);
#endif
#if SIZE_MAX > INT_MAX
    ASSERT(x <= INT_MAX);
#endif
  return (int) x;
}

static inline int
cast_size_int_fast(size_t x)
{
#if 0 < INT_MIN
    ASSERT_CONVERSION(x >= INT_MIN);
#endif
#if SIZE_MAX > INT_MAX
    ASSERT_CONVERSION(x <= INT_MAX);
#endif
  return (int) x;
}

static inline unsigned short int
cast_size_ushort(size_t x)
{
#if SIZE_MAX > USHRT_MAX
    ASSERT(x <= USHRT_MAX);
#endif
  return (unsigned short int) x;
}

static inline unsigned short int
cast_size_ushort_fast(size_t x)
{
#if SIZE_MAX > USHRT_MAX
    ASSERT_CONVERSION(x <= USHRT_MAX);
#endif
  return (unsigned short int) x;
}

static inline short int
cast_size_short(size_t x)
{
#if 0 < SHRT_MIN
    ASSERT(x >= SHRT_MIN);
#endif
#if SIZE_MAX > SHRT_MAX
    ASSERT(x <= SHRT_MAX);
#endif
  return (short int) x;
}

static inline short int
cast_size_short_fast(size_t x)
{
#if 0 < SHRT_MIN
    ASSERT_CONVERSION(x >= SHRT_MIN);
#endif
#if SIZE_MAX > SHRT_MAX
    ASSERT_CONVERSION(x <= SHRT_MAX);
#endif
  return (short int) x;
}

static inline unsigned char
cast_size_uchar(size_t x)
{
#if SIZE_MAX > UCHAR_MAX
    ASSERT(x <= UCHAR_MAX);
#endif
  return (unsigned char) x;
}

static inline unsigned char
cast_size_uchar_fast(size_t x)
{
#if SIZE_MAX > UCHAR_MAX
    ASSERT_CONVERSION(x <= UCHAR_MAX);
#endif
  return (unsigned char) x;
}

static inline char
cast_size_char(size_t x)
{
#if 0 < CHAR_MIN
    ASSERT(x >= CHAR_MIN);
#endif
#if SIZE_MAX > CHAR_MAX
    ASSERT(x <= CHAR_MAX);
#endif
  return (char) x;
}

static inline char
cast_size_char_fast(size_t x)
{
#if 0 < CHAR_MIN
    ASSERT_CONVERSION(x >= CHAR_MIN);
#endif
#if SIZE_MAX > CHAR_MAX
    ASSERT_CONVERSION(x <= CHAR_MAX);
#endif
  return (char) x;
}

#ifdef __cplusplus
}
#endif

#endif	/* TYPECAST_H_ */
