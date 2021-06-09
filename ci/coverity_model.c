// #include "decls.h"

/* This is just so that coverity recognizes posix_memalign */
int posix_memalign(void **memptr, size_t alignment, size_t size)
{
    
   __coverity_negative_sink__(size);
   __coverity_negative_sink__(alignment);

    void * mem = __coverity_alloc__(size);
    if (mem) {
        *memptr = mem;
        return 0;
    } else {
        return 1;
    }
}
