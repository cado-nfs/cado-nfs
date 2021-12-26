#ifndef CADO_POPEN_H_
#define CADO_POPEN_H_
// IWYU pragma: no_include <bits/types/struct_rusage.h>
// IWYU pragma: no_forward_declare rusage
#include "cado_config.h"  // just because we're a header.
#include <stdio.h>
#ifdef HAVE_GETRUSAGE
#include <sys/resource.h> // IWYU pragma: keep
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifndef HAVE_MINGW
FILE * cado_popen(const char * command, const char * mode);
int cado_fd_popen(const char * command, const char * mode);
#if defined(HAVE_GETRUSAGE) && defined(HAVE_WAIT4)
int cado_pclose2(FILE * stream, struct rusage * r);
int cado_fd_pclose2(int fd, struct rusage * r);
#else
int cado_pclose2(FILE * stream, void * r);
int cado_fd_pclose2(int fd, void * r);
#endif
static inline int cado_pclose(FILE * stream) { return cado_pclose2(stream, NULL); }
static inline int cado_fd_pclose(int fd) { return cado_fd_pclose2(fd, NULL); }
#else

static inline FILE * cado_popen(const char * command, const char * mode) { return popen(command, mode); }
static inline int cado_pclose(FILE * stream) { return pclose(stream); }
/* we don't even provide cado_pclose2 for mingw */
#endif

#ifdef __cplusplus
}
#endif

#endif	/* CADO_POPEN_H_ */
