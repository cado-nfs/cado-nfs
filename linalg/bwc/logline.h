#ifndef LOGLINE_H_
#define LOGLINE_H_

#include <stdio.h>
#include <stdarg.h>
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

void logline_decl_usage(param_list_ptr pl);
void logline_init_timer();
int logline_interpret_parameters(param_list pl);
int logline_begin(FILE * f, size_t size, const char * fmt, ...) ATTR_PRINTF(3,4);
int logline_end(double *, const char * fmt, ...);
int logline_vprintf(int level, const char * fmt, va_list ap);
int logline_printf(int level, const char * fmt, ...) ATTR_PRINTF(2,3);

#ifdef __cplusplus
}
#endif

#endif	/* LOGLINE_H_ */
