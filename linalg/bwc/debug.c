#include "cado.h"

#include <stdio.h>
#include <stdarg.h>
#include "debug.h"

void debug_write(const void * v,
        size_t n, const char * fmt, ...)
{
    char * tmp;
    va_list ap;

    va_start(ap, fmt);
    int rc = vasprintf(&tmp, fmt, ap);
    FILE * f = fopen(tmp, "w");
    rc = fwrite(v, 1, n, f);
    fclose(f);
    va_end(ap);
    free(tmp);
}
