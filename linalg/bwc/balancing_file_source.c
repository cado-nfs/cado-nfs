#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "balancing_file_source.h"
#include "utils.h"

// there's presently no provision for tweaking the following parameter,
// but adding it is clearly a completely trivial matter.
const size_t disk_io_window_size = 1 << 14;

size_t file_source_get(file_source_ptr f, uint32_t ** p, size_t avail)
{
    if (*p == NULL || avail == 0) {
        if (f->buf == NULL) {
            f->buf = malloc(disk_io_window_size * sizeof(uint32_t));
        }
        avail = disk_io_window_size;
        *p = f->buf;
    }
    size_t n = fread32_little(*p, avail, f->f);
    f->b->pos += n;
    if (n == 0) {
        if (ferror(f->f)) {
            perror(f->filename);
            exit(1);
        } else {
            if (f->b->pos * sizeof(uint32_t) != (size_t) f->sbuf->st_size) {
                fprintf(stderr, "Ended at wrong position (%zu != %zu) in %s\n",
                        f->b->pos, (size_t) f->sbuf->st_size, f->filename);
                exit(1);
            }
        }
    }
    return n;
}

data_source_ptr file_source_alloc(const char * filename, size_t esz)
{
    file_source_ptr p = malloc(sizeof(file_source));
    memset(p, 0, sizeof(file_source)); 
    p->filename = filename;
    p->f = fopen(filename, "r");
    if (p->f == NULL) {
	perror(filename);
	exit(1);
    }
    fstat(fileno(p->f), p->sbuf);
    if (esz) {
        if ((size_t) p->sbuf->st_size != esz) {
            fprintf(stderr, "%s: expected size %zu, not %zu\n",
                    filename, esz, (size_t) p->sbuf->st_size);
            exit(1);
        }
    }
    p->b->get = (size_t(*)(void*,uint32_t **,size_t)) file_source_get;
    return (data_source_ptr) p;
}

void file_source_free(data_source_ptr q)
{
    file_source_ptr p = (file_source_ptr) q;
    fclose(p->f);
    if (p->buf) free(p->buf);
    free(q);
}

