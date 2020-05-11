#ifndef CADO_POPEN_H_
#define CADO_POPEN_H_

// IWYU pragma: private, include "utils.h"

#include <stdio.h>
#include <unistd.h>
#ifdef HAVE_GETRUSAGE
#include <sys/resource.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef HAVE_MINGW
FILE * cado_popen(const char * command, const char * mode);
#ifdef HAVE_GETRUSAGE
int cado_pclose2(FILE * stream, struct rusage * r);
#else
int cado_pclose2(FILE * stream, void * r);
#endif
static inline int cado_pclose(FILE * stream) { return cado_pclose2(stream, NULL); }
#else

static inline FILE * cado_popen(const char * command, const char * mode) { return popen(command, mode); }
static inline int cado_pclose(FILE * stream) { return pclose(stream); }
/* we don't even provide cado_pclose2 for mingw */
#endif

#ifdef __cplusplus
}
#endif

#endif	/* CADO_POPEN_H_ */
