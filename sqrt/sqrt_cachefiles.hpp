#ifndef SQRT_SQRT_CACHEFILES_HPP_
#define SQRT_SQRT_CACHEFILES_HPP_

#include <cstdio>

/*  cache files */
extern int rcache;
extern int wcache;
#define CACHEDIR        "/tmp"
#define CACHEPREFIX        "CRTALGSQRT."

struct cachefile_s {
    char basename[128];
    FILE * f;
    int writing;
};

typedef struct cachefile_s cachefile[1];
typedef struct cachefile_s * cachefile_ptr;

extern void cachefile_vinit(cachefile_ptr c, const char * fmt, va_list ap);
extern void cachefile_init(cachefile_ptr c, const char * fmt, ...);
extern int cachefile_open_w(cachefile_ptr c);
extern int cachefile_open_r(cachefile_ptr c);
extern int cachefile_exists(const char * fmt, ...);
extern void cachefile_close(cachefile_ptr c);

#endif	/* SQRT_SQRT_CACHEFILES_HPP_ */
