#ifndef CADO_POPEN_HPP
#define CADO_POPEN_HPP

#include "cado_config.h"  // just because we're a header.

#include <cstdio>
#ifdef HAVE_GETRUSAGE
#include <sys/resource.h> // IWYU pragma: keep
#endif

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

#endif	/* CADO_POPEN_HPP */
