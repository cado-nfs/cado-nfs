#ifndef CADO_MATMUL_FACADE_HPP
#define CADO_MATMUL_FACADE_HPP

/* This file is included by all matmul implementations so as to enable
 * the macro MY_MATMUL_INTERFACE
 *
 * It does not make sense to have it included from the higher-level code
 * (in fact, it will error out).
 */

#if !defined(ARITH_LAYER)
#error "Please compile this file with the ARITH_LAYER macro defined"
#endif

#if !defined(MM_IMPL)
#error "Please compile this file with the MM_IMPL macro defined"
#endif

#include <cstdarg>  // for va_list
#include <cstddef>  // for size_t
#include <cstdint>  // for uint32_t

#include "macros.h"
#include "params.h"  // for param_list
#include "matmul.hpp"

/* Now some cpp glue */
#define MY_MATMUL_INTERFACE CADO_CONCATENATE4(matmul_, ARITH_LAYER, _, MM_IMPL)

#endif	/* MATMUL_FACADE_HPP_ */
