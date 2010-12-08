#ifndef ABASE_GENERIC_H_
#define ABASE_GENERIC_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#ifdef  HAVE_ALLOCA_H
#include <alloca.h>
#endif
#include "random_generation.h"

/* This header contains the stuff that pertains to the abase interface,
 * but which is known to depend only on the data striding. Hence simple
 * bindings may be used */

typedef void * abase_generic_ptr;
#define abase_generic_ptr_add(p,k) (((void*)(p))+(k))

// I know arithmetic on void pointers is not C99. So if you insist, that
// sh*t should work.
// #define abase_generic_ptr_add(p,k) ((void*) (((char*)(p)) + (k)))


#define abase_generic_init(b,n) alignable_malloc((b)*next_multiple_of_powerof2((n),ABASE_UNIVERSAL_READAHEAD_ITEMS))
#define abase_generic_clear(b,p,n) alignable_free((p),(b)*next_multiple_of_powerof2((n),ABASE_UNIVERSAL_READAHEAD_ITEMS))
#ifdef  HAVE_ALLOCA_H
#define abase_generic_initf(b,n) alloca((b)*next_multiple_of_powerof2((n),ABASE_UNIVERSAL_READAHEAD_ITEMS))
#define abase_generic_clearf(b,p,n) /**/
#else
#define abase_generic_initf(b,n) malloc((b)*next_multiple_of_powerof2((n),ABASE_UNIVERSAL_READAHEAD_ITEMS))
#define abase_generic_clearf(b,p,n) free(p)
#endif
#define abase_generic_zero(b,p,n) memset((p),0,(b)*(n))
#define abase_generic_copy(b,q,p,n) memcpy((q),(p),(b)*(n))
#define abase_generic_write(b,p,n) fwrite((p), (b), (n), (f))
#define abase_generic_read(b,p,n) fread((p), (b), (n), (f))
#define abase_generic_random(b,p,n) myrand_area((p), (b)*(n))

#endif	/* ABASE_GENERIC_H_ */
