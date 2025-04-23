#ifndef CADO_PURGEDFILE_H
#define CADO_PURGEDFILE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void purgedfile_read_firstline (const char *, uint64_t *, uint64_t *);

#ifdef __cplusplus
}
#endif

#endif	/* CADO_PURGEDFILE_H */
