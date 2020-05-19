#ifndef MPFQ_LAYER_H_
#define MPFQ_LAYER_H_
// pragma no prototypes

#if defined(SELECT_MPFQ_LAYER_u64k1) || defined(SELECT_MPFQ_LAYER_u64)
#include "mpfq/mpfq_u64k1.h"
#elif defined(SELECT_MPFQ_LAYER_u64k2)
#include "mpfq/mpfq_u64k2.h"
#elif defined(SELECT_MPFQ_LAYER_u64k3)
#include "mpfq/mpfq_u64k3.h"
#elif defined(SELECT_MPFQ_LAYER_u64k4)
#include "mpfq/mpfq_u64k4.h"
#elif defined(SELECT_MPFQ_LAYER_m128)
#include "mpfq/mpfq_m128.h"
#elif defined(SELECT_MPFQ_LAYER_u64n)
#error "argh"
#include "mpfq/mpfq_u64n.h"
#elif defined(SELECT_MPFQ_LAYER_u128)
#error "argh"
#include "mpfq/mpfq_u128.h"
#elif defined(SELECT_MPFQ_LAYER_p16) /* This is really the first non-gf2 try */
#define NOT_OVER_GF2
#include "mpfq/mpfq_p16.h"

/* In reality we don't have everything configured for the moment. Update
 * CMakeLists.txt when the need for other moduli arises */

#elif defined(SELECT_MPFQ_LAYER_p_1)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_1.h"
#elif defined(SELECT_MPFQ_LAYER_p_2)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_2.h"
#elif defined(SELECT_MPFQ_LAYER_p_3)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_3.h"
#elif defined(SELECT_MPFQ_LAYER_p_4)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_4.h"
#elif defined(SELECT_MPFQ_LAYER_p_5)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_5.h"
#elif defined(SELECT_MPFQ_LAYER_p_6)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_6.h"
#elif defined(SELECT_MPFQ_LAYER_p_7)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_7.h"
#elif defined(SELECT_MPFQ_LAYER_p_8)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_8.h"
#elif defined(SELECT_MPFQ_LAYER_p_9)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_9.h"
#elif defined(SELECT_MPFQ_LAYER_p_10)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_10.h"
#elif defined(SELECT_MPFQ_LAYER_p_11)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_11.h"
#elif defined(SELECT_MPFQ_LAYER_p_12)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_12.h"
#elif defined(SELECT_MPFQ_LAYER_p_13)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_13.h"
#elif defined(SELECT_MPFQ_LAYER_p_14)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_14.h"
#elif defined(SELECT_MPFQ_LAYER_p_15)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_15.h"
#elif defined(SELECT_MPFQ_LAYER_pz)
#define NOT_OVER_GF2
#include "mpfq/mpfq_pz.h"
#else
// #warning "Using default selection for abase"
#error "argh. This code must be compiled with some SELECT_MPFQ_LAYER_ macro defined"
// #include "mpfq/mpfq_u64.h"
#endif

/* This is used as a shorthand throughout in order to ease the access to
 * the _primary_ abase. Other abases have to be accessed via the OO
 * interface.
 */
#include "mpfq/mpfq_name_ab.h"

#endif	/* MPFQ_LAYER_H_ */
