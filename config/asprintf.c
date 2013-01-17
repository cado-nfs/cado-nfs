
/* These feature test macros, under most well-behaving unices at least,
 * expose some available form of asprintf and vasprintf, which are used
 * in several places in cado-nfs. Should the system fail to provide these
 * functions, it is possible to provide workalikes */

#ifdef __OpenBSD__
#define _BSD_SOURCE
#else
#define _POSIX_C_SOURCE 200112L
#define _XOPEN_SOURCE   600
#define _BSD_SOURCE
#define _ISOC99_SOURCE
#ifndef __cplusplus
#define _GNU_SOURCE
#endif
#define _DARWIN_C_SOURCE
#define _NETBSD_SOURCE
#endif

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

char * foo(const char * fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    char * tmp;
    vasprintf(&tmp, fmt, ap);
    va_end(ap);
    return tmp;
}

int main()
{
    char * tmp;
    asprintf(&tmp, "Hello %d\n", 42);
    free(tmp);
    tmp = foo("Hello %d\n", 17);
    free(tmp);
    return 0;
}

