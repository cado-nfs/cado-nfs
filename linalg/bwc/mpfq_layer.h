#ifndef MPFQ_LAYER_H_
#define MPFQ_LAYER_H_
// pragma no prototypes

#if defined(MPFQ_FAKE_HPP_) && defined(MPFQ_LAYER_H_)
#error "mpfq_layer.h and mpfq_fake.hpp are incomaptible"
#endif

// IWYU pragma: begin_exports
#include "mpfq/mpfq.h"

#if defined(SELECT_MPFQ_LAYER_u64k1) || defined(SELECT_MPFQ_LAYER_u64)
#include "mpfq/mpfq_u64k1.h"
#ifdef __cplusplus
#include "mpfq/mpfq_u64k1.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_u64k2)
#include "mpfq/mpfq_u64k2.h"
#ifdef __cplusplus
#include "mpfq/mpfq_u64k2.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_u64k3)
#include "mpfq/mpfq_u64k3.h"
#ifdef __cplusplus
#include "mpfq/mpfq_u64k3.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_u64k4)
#include "mpfq/mpfq_u64k4.h"
#ifdef __cplusplus
#include "mpfq/mpfq_u64k4.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_m128)
#include "mpfq/mpfq_m128.h"
#ifdef __cplusplus
#include "mpfq/mpfq_m128.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_p16) /* This is really the first non-gf2 try */
#define NOT_OVER_GF2
#include "mpfq/mpfq_p16.h"
#ifdef __cplusplus
#include "mpfq/mpfq_p16.hpp"
#endif /* __cplusplus */

/* In reality we don't have everything configured for the moment. Update
 * CMakeLists.txt when the need for other moduli arises */

#elif defined(SELECT_MPFQ_LAYER_p_1)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_1.h"
#ifdef __cplusplus
#include "mpfq/mpfq_p_1.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_p_2)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_2.h"
#ifdef __cplusplus
#include "mpfq/mpfq_p_2.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_p_3)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_3.h"
#ifdef __cplusplus
#include "mpfq/mpfq_p_3.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_p_4)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_4.h"
#ifdef __cplusplus
#include "mpfq/mpfq_p_4.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_p_5)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_5.h"
#ifdef __cplusplus
#include "mpfq/mpfq_p_5.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_p_6)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_6.h"
#ifdef __cplusplus
#include "mpfq/mpfq_p_6.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_p_7)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_7.h"
#ifdef __cplusplus
#include "mpfq/mpfq_p_7.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_p_8)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_8.h"
#ifdef __cplusplus
#include "mpfq/mpfq_p_8.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_p_9)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_9.h"
#ifdef __cplusplus
#include "mpfq/mpfq_p_9.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_p_10)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_10.h"
#ifdef __cplusplus
#include "mpfq/mpfq_p_10.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_p_11)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_11.h"
#ifdef __cplusplus
#include "mpfq/mpfq_p_11.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_p_12)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_12.h"
#ifdef __cplusplus
#include "mpfq/mpfq_p_12.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_p_13)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_13.h"
#ifdef __cplusplus
#include "mpfq/mpfq_p_13.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_p_14)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_14.h"
#ifdef __cplusplus
#include "mpfq/mpfq_p_14.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_p_15)
#define NOT_OVER_GF2
#include "mpfq/mpfq_p_15.h"
#ifdef __cplusplus
#include "mpfq/mpfq_p_15.hpp"
#endif /* __cplusplus */
#elif defined(SELECT_MPFQ_LAYER_pz)
#define NOT_OVER_GF2
#include "mpfq/mpfq_pz.h"
#ifdef __cplusplus
#include "mpfq/mpfq_pz.hpp"
#endif /* __cplusplus */
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
// IWYU pragma: end_exports

#endif	/* MPFQ_LAYER_H_ */
