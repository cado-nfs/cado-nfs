#include "cado.h" // IWYU pragma: keep

#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cstdarg>
#include <cerrno>

#include <fcntl.h>
#include <unistd.h>

#include "sqrt_cachefiles.hpp"

int rcache = 1;
int wcache = 1;

void cachefile_vinit(cachefile_ptr c, const char * fmt, va_list ap)
{
    memset(c, 0, sizeof(*c));
    vsnprintf(c->basename, sizeof(c->basename), fmt, ap);
}

void cachefile_init(cachefile_ptr c, const char * fmt, ...)
{
    va_list ap;
    va_start(ap, fmt);
    cachefile_vinit(c, fmt, ap);
    va_end(ap);
}

int cachefile_open_w(cachefile_ptr c)
{
    char tname[256];
    char fname[256];
    snprintf(fname, sizeof(fname), CACHEDIR "/" CACHEPREFIX "%s", c->basename);
    snprintf(tname, sizeof(tname), CACHEDIR "/.pre." CACHEPREFIX "%s", c->basename);
    /* first create the destination file atomically. If it already
     * exists, we do nothing.
     */
    const int fd = open(fname, O_WRONLY | O_CREAT | O_EXCL, 0666);
    if (fd < 0) {
        if (errno != EEXIST)
            fprintf(stderr, "open(%s): %s\n", tname, strerror(errno));
        return 0;
    }
    /* close it: it will be replaced by the contents of the temp file
     * when we're done */
    close(fd);

    c->f = fopen(tname, "w");
    c->writing = 0;
    if (c->f == nullptr) {
        fprintf(stderr, "%s: %s\n", tname, strerror(errno));
    }
    c->writing = 1;
    return c->f != nullptr;
}

int cachefile_open_r(cachefile_ptr c)
{
    char tname[256];
    snprintf(tname, sizeof(tname), CACHEDIR "/" CACHEPREFIX "%s", c->basename);
    c->f = fopen(tname, "r");
    c->writing = 0;
    return c->f != nullptr;
}

int cachefile_exists(const char * fmt, ...)
{
    cachefile c;
    char tname[256];
    va_list ap;
    va_start(ap, fmt);
    cachefile_vinit(c, fmt, ap);
    va_end(ap);
    snprintf(tname, sizeof(tname), CACHEDIR "/" CACHEPREFIX "%s", c->basename);
    return access(tname, R_OK) == 0;
}

void cachefile_close(cachefile_ptr c)
{
    char tname[256];
    char fname[256];
    fclose(c->f);
    if (c->writing == 0)
        return;
    snprintf(tname, sizeof(tname), CACHEDIR "/.pre." CACHEPREFIX "%s", c->basename);
    snprintf(fname, sizeof(fname), CACHEDIR "/" CACHEPREFIX "%s", c->basename);
    if (rename(tname, fname) < 0)
        fprintf(stderr, "Could not rename temporary cache file %s to %s: %s\n", tname, fname, strerror(errno));
}
