#include "cado.h"
#include "readbuffer.h"
#include <string.h>
#include <stdlib.h>
#include <errno.h>

#define DIE_ERRNO_DIAG(tst, func, arg) do {                             \
    if ((tst)) {                                                        \
        fprintf(stderr, func "(%s): %s\n", arg, strerror(errno));       \
        exit(1);                                                        \
    }                                                                   \
} while (0)

void rb_open(reading_buffer b, const char * filename)
{
    b->filename = strdup(filename);
    b->f = fopen(filename, "r");
    DIE_ERRNO_DIAG(b->f == NULL, "fopen", b->filename);
}
void rb_close(reading_buffer b)
{
    fclose(b->f);
    free(b->filename);
}

void rb_read_line(reading_buffer b)
{
    char * ptr;
    b->o = ftello(b->f);                                                  
    ptr = fgets(b->buf, sizeof(b->buf), b->f);
    if (ptr == NULL) {
        fprintf(stderr, "Unexpected %s at position %ld in %s\n",    
                feof(b->f) ? "EOF" : "error", (long) b->o, b->filename);            
        exit(1);
    }
    b->siz = strlen(ptr);
}

void rb_gobble_long_line(reading_buffer b)
{
    char * ptr;
    while (!(b->siz < (int) (sizeof(b->buf)-1) ||
                b->buf[sizeof(b->buf)-2] == '\n'))
    {
        ptr = fgets(b->buf, sizeof(b->buf), b->f);
        if (ptr == NULL) {
            fprintf(stderr, "Unexpected %s at position %ld in %s\n",    
                    feof(b->f) ? "EOF" : "error",
                    (long) ftello(b->f), b->filename);    
            exit(1);
        }
        b->siz = strlen(ptr);
    }                                                                   
    b->o = ftello(b->f);
}

void rb_feed_buffer_again_if_lowwater(reading_buffer b,
        int * ppos, int low)
{
    char * ptr;
    if (*ppos+low <= (int) sizeof(b->buf) ||
            (b->siz < (int) (sizeof(b->buf)-1) ||
                b->buf[sizeof(b->buf)-2] == '\n')) {
        return;
    }

    memcpy(b->buf, b->buf+*ppos, b->siz-*ppos);
    b->siz = b->siz-*ppos;
    *ppos=0;
    ptr = fgets(b->buf+b->siz,sizeof(b->buf)-b->siz,b->f);
    if (ptr == NULL) {
        fprintf(stderr, "Unexpected %s at position %ld in %s\n",    
                feof(b->f) ? "EOF" : "error",
                (long) ftello(b->f), b->filename);    
        exit(1);
    }
    b->siz += strlen(b->buf + b->siz);
}
/*  */
