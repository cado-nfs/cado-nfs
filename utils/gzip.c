#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "gzip.h"
#include "misc.h"

struct suffix_handler {
    const char * suffix;
    const char * pfmt_in;
    const char * pfmt_out;
};

#if 0
const char * suffix = NULL;

const char * copy_suffix_noalloc(const char * name)
{
    const char * p = strrchr(name, '.');
    if (p == NULL)
        p = name + strlen(name);
    return strdup(p);
}

const char * copy_suffix_alloc(const char * name)
{
    return strdup(copy_suffix_noalloc(name));
}
const char * path_remove_suffix(char * name)
{
    char * p = strrchr(name, '.');
    if (p) *p=0;
    return name;
}

#endif

struct suffix_handler supported_compression_formats[] = {
    { ".gz", "gzip -dc %s", "gzip -c --best > %s", },
    { ".bz2", "bzip2 -dc %s", "bzip2 -c --fast > %s", },
    /* These two have to be present */
    { "", NULL, NULL },
    { NULL, NULL, NULL },
};

const char * path_basename(const char * path)
{
    const char *p = strrchr(path, '/');
    if (p == NULL) {
        p = path;
    } else {
        p = p + 1;
    }
    return p;
}

int is_supported_compression_format(const char * s)
{
    struct suffix_handler * r = supported_compression_formats;
    for( ; r->suffix ; r++) {
        if (strcmp(r->suffix, s) == 0)
            return 1;
    }
    return 0;
}

FILE * fopen_compressed_r(const char * name, int* p_pipeflag, char const ** suf)
{
    const struct suffix_handler * r = supported_compression_formats;
    FILE * f;
    for( ; r->suffix ; r++) {
        if (!has_suffix(name, r->suffix))
            continue;
        if (suf) *suf = r->suffix;
        if (r->pfmt_in) {
            char * command;
            int ret = asprintf(&command, r->pfmt_in, name);
            ASSERT_ALWAYS(ret >= 0);
            f = popen(command, "r");
            free(command);
            if (p_pipeflag) *p_pipeflag = 1;
            return f;
        } else {
            f = fopen(name, "r");
            if (p_pipeflag) *p_pipeflag = 0;
            return f;
        }
    }
    return NULL;
}

FILE * fopen_compressed_w(const char * name, int* p_pipeflag, char const ** suf)
{
    const struct suffix_handler * r = supported_compression_formats;
    FILE * f;
    for( ; r->suffix ; r++) {
        if (!has_suffix(name, r->suffix))
            continue;
        if (suf) *suf = r->suffix;
        if (r->pfmt_out) {
            char * command;
            int ret = asprintf(&command, r->pfmt_out, name);
            ASSERT_ALWAYS(ret >= 0);
            f = popen(command, "w");
            free(command);
            if (p_pipeflag) *p_pipeflag = 1;
            return f;
        } else {
            f = fopen(name, "w");
            if (p_pipeflag) *p_pipeflag = 0;
            return f;
        }
    }
    return NULL;
}

FILE * gzip_open(const char * name, const char * mode)
{
    if (strcmp(mode, "r") == 0) {
        return fopen_compressed_r(name, NULL, NULL);
    } else if (strcmp(mode, "w") == 0) {
        return fopen_compressed_w(name, NULL, NULL);
    } else {
        abort();
    }
}

void gzip_close(FILE * f, const char * name)
{
    const char * e = name + strlen(name);
    const struct suffix_handler * r = supported_compression_formats;
    for( ; r->suffix ; r++) {
        const char * ne = e - MIN(strlen(name), strlen(r->suffix));
        if (strcmp(r->suffix, ne) != 0)
            continue;
        if (r->pfmt_out) {
            pclose(f);
            return;
        } else {
            fclose(f);
            return;
        }
    }
}

