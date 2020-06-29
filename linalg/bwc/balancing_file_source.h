#ifndef BALANCING_FILE_SOURCE_H_
#define BALANCING_FILE_SOURCE_H_

#include <stdint.h>                 // for uint32_t
#include <stdio.h>
#include <sys/stat.h>

#include "balancing_data_source.h"

struct file_source_s { 
    data_source b;
    const char * filename;
    struct stat sbuf[1];
    FILE * f;
    uint32_t * buf;
};
typedef struct file_source_s file_source[1];
typedef struct file_source_s *file_source_ptr;

#ifdef __cplusplus
extern "C" {
#endif

extern size_t file_source_get(file_source_ptr f, uint32_t ** p, size_t avail);
extern data_source_ptr file_source_alloc(const char * filename, size_t esz);
extern void file_source_free(data_source_ptr q);


#ifdef __cplusplus
}
#endif

#endif	/* BALANCING_FILE_SOURCE_H_ */
