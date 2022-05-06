/* Data struct used for polyselect */
#include "cado.h" // IWYU pragma: keep
#include "polyselect_special_q.h"

/* The following are primes used as factors in the special-q.
   Warning: if you add larger primes, you should ensure that any product
   still fits in an "unsigned long" (cf routine collision_on_sq).
   Here on a 32-bit machine the largest product of 4 primes is
   241*251*257*263 < 2^32, and on a 64-bit machine the largest product of 8
   primes is 233*239*241*251*257*263*269*271 < 2^64.
   LEN_SPECIAL_Q is defined in the header.
*/
#if ULONG_MAX == 4294967295UL
const unsigned int SPECIAL_Q[LEN_SPECIAL_Q] = {
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
  31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
  73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
  127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
  179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
  233, 239, 241, 251, 257, 263, 0 };
#else /* 64-bit */
const unsigned int SPECIAL_Q[LEN_SPECIAL_Q] = {
  2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
  31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
  73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
  127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
  179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
  233, 239, 241, 251, 257, 263, 269, 271, 0 };
#endif
