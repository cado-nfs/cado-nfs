#ifndef PURGEDFILE_H_
#define PURGEDFILE_H_
#include <stdio.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void purgedfile_read_firstline (const char *, uint64_t *, uint64_t *);

#ifdef __cplusplus
}
#endif

#endif	/* PURGEDFILE_H_ */
