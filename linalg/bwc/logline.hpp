#ifndef CADO_LOGLINE_HPP
#define CADO_LOGLINE_HPP

#include <cstdarg>  // for va_list
#include <cstdio>   // for FILE, size_t
#include "macros.h"  // for ATTR_PRINTF
#include "params.h"  // for cxx_param_list

void logline_decl_usage(cxx_param_list & pl);
void logline_init_timer();
/* This queries (or sets) the time since timer init. These two functions
 * may be useful in checkpoints */
double logline_serialize();
void logline_unserialize(double);
int logline_interpret_parameters(cxx_param_list & pl);
int logline_begin(FILE * f, size_t size, const char * fmt, ...) ATTR_PRINTF(3,4);
int logline_end(double *, const char * fmt, ...);
int logline_vprintf(int level, const char * fmt, va_list ap);
int logline_printf(int level, const char * fmt, ...) ATTR_PRINTF(2,3);

#endif	/* LOGLINE_HPP_ */
