#include "cado.h"       /* feature macros, no includes */
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <ctype.h>
/* For MinGW Build */
#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#endif

#include "cado_config.h"
#include "macros.h"
#include "misc.h"

#ifndef HAVE_STRDUP
char *strdup(const char * const s)
{
    const size_t size = strlen(s) + 1;
    char *r = (char *) malloc(size * sizeof(char));
    if (r != NULL)
        memcpy(r, s, size);
    return r;
}
#endif


#ifndef HAVE_STRNDUP
/* Not every libc has this, and providing a workalike is very easy */
char *strndup(const char * const a, const size_t n)
{
    const size_t size = MIN(strlen(a), n) + 1;
    char *r = (char *) malloc(size * sizeof(char));
    if (r != NULL) {
        memcpy(r, a, size);
        r[size] = '\0';
    }
    return r;
}
#endif

void
*malloc_check (const size_t x)
{
    void *p;
    p = malloc (x);
    ASSERT_ALWAYS(p != NULL);
    return p;
}

/* Not everybody has posix_memalign. In order to provide a viable
 * alternative, we need an ``aligned free'' matching the ``aligned
 * malloc''. We rely on posix_memalign if it is available, or else fall
 * back on ugly pointer arithmetic so as to guarantee alignment. Note
 * that not providing the requested alignment can have some troublesome
 * consequences. At best, a performance hit, at worst a segv (sse-2
 * movdqa on a pentium4 causes a GPE if improperly aligned).
 */

void *malloc_aligned(size_t size, size_t alignment)
{
#ifdef HAVE_POSIX_MEMALIGN
    void *res = NULL;
    int rc = posix_memalign(&res, alignment, size);
    ASSERT_ALWAYS(rc == 0);
    return res;
#else
    char * res;
    res = malloc(size + sizeof(size_t) + alignment);
    res += sizeof(size_t);
    size_t displ = alignment - ((unsigned long) res) % alignment;
    res += displ;
    memcpy(res - sizeof(size_t), &displ, sizeof(size_t));
    ASSERT_ALWAYS((((unsigned long) res) % alignment) == 0);
    return (void*) res;
#endif
}

void free_aligned(void * p, size_t size MAYBE_UNUSED, size_t alignment MAYBE_UNUSED)
{
#ifdef HAVE_POSIX_MEMALIGN
    free(p);
#else
    char * res = (char *) p;
    ASSERT_ALWAYS((((unsigned long) res) % alignment) == 0);
    size_t displ;
    memcpy(&displ, res - sizeof(size_t), sizeof(size_t));
    res -= displ;
    ASSERT_ALWAYS(displ == alignment - ((unsigned long) res) % alignment);
    res -= sizeof(size_t);
    free(res);
#endif
}

static long
pagesize (void)
{
#if defined(_WIN32) || defined(_WIN64)
  /* cf http://en.wikipedia.org/wiki/Page_%28computer_memory%29 */
  SYSTEM_INFO si;
  GetSystemInfo(&si);
  return si.dwPageSize;
#else
  return sysconf (_SC_PAGESIZE);
#endif
}

void *malloc_pagealigned(size_t sz)
{
    void *p = malloc_aligned (sz, pagesize ());
    ASSERT_ALWAYS(p != NULL);
    return p;
}

void free_pagealigned(void * p, size_t sz)
{
    free_aligned(p, sz, pagesize ());
}

int has_suffix(const char * path, const char * sfx)
{
    unsigned int lp = strlen(path);
    unsigned int ls = strlen(sfx);
    if (lp < ls) return 0;
    return strcmp(path + lp - ls, sfx) == 0;
}

// given a path to a file (prefix), and a suffix called (what), returns
// if the ext parameter is NULL, a malloced string equal to
// (prefix).(what) ; if ext is non-null AND (ext) is already a suffix of
// (prefix), say we have (prefix)=(prefix0)(ext), then we return
// (prefix0).(what)(ext)
// It is typical to use ".bin" or ".txt" as ext parameters.
char * derived_filename(const char * prefix, const char * what, const char * ext)
{
    char * dup_prefix;
    dup_prefix=strdup(prefix);

    if (ext && has_suffix(dup_prefix, ext)) {
        dup_prefix[strlen(dup_prefix)-strlen(ext)]='\0';
    }
    char * str;
    int rc = asprintf(&str, "%s.%s%s", dup_prefix, what, ext ? ext : "");
    if (rc<0) abort();
    free(dup_prefix);
    return str;
}


void chomp(char *s) {
    char *p;
    if (s && (p = strrchr(s, '\n')) != NULL)
        *p = '\0';
}

char ** filelist_from_file(const char * basepath, const char * filename,
                           int typ)
{
    char ** files = NULL;
    int nfiles_alloc = 0;
    int nfiles = 0;
    FILE *f;
    f = fopen(filename, "r");
    if (f == NULL) {
      if (typ == 0)
        perror ("Problem opening filelist");
      else
        perror ("Problem opening subdirlist");
      exit (1);
    }
    char relfile[FILENAME_MAX + 10];
    while (fgets(relfile, FILENAME_MAX + 10, f) != NULL) {

        // skip leading blanks
        char *rfile = relfile;
        while (isspace((int)(unsigned char)rfile[0]))
            rfile++;
        // if empty line or comment line, continue
        if ((rfile[0] == '#') || (rfile[0] == '\0') || (rfile[0] == '\n'))
            continue;
        chomp(rfile);

        if (nfiles == nfiles_alloc) {
            nfiles_alloc += nfiles_alloc / 2 + 16;
            files = realloc(files, nfiles_alloc * sizeof(char*));
        }
        if (basepath) {
            char * name;
            int ret = asprintf(&name, "%s/%s", basepath, rfile);
            ASSERT_ALWAYS(ret >= 0);
            files[nfiles] = name;
        } else {
            files[nfiles] = strdup(rfile);
        }
        nfiles++;
    }
    fclose(f);

    if (nfiles == nfiles_alloc) {
        nfiles_alloc += nfiles_alloc / 2 + 16;
        files = realloc(files, nfiles_alloc * sizeof(char*));
    }
    files[nfiles++] = NULL;
    return files;
}

void filelist_clear(char ** filelist)
{
    for(char ** p = filelist ; *p ; p++)
        free(*p);
    free(filelist);
}

int mkdir_with_parents(const char * dir, int fatal)
{
    char * tmp = strdup(dir);
    int n = strlen(dir);
    int pos = 0;
    if (dir[0] == '/')
        pos++;
    for( ; pos < n ; ) {
        for( ; dir[pos] == '/' ; pos++) ;
        if (pos == n) break;
        const char * slash = strchr(dir + pos, '/');
        strncpy(tmp, dir, n);
        if (slash) {
            pos = slash - dir;
            tmp[pos]='\0';
        } else {
            pos = n;
        }
        struct stat sbuf[1];
        int rc = stat(tmp, sbuf);
        if (rc < 0) {
            if (errno != ENOENT) {
                fprintf(stderr, "accessing %s: %s\n", tmp, strerror(errno));
                free(tmp);
                if (fatal) exit(1);
                return -errno;
            }
/* MinGW's mkdir has only one argument,
   cf http://lists.gnu.org/archive/html/bug-gnulib/2008-04/msg00259.html */
#if (defined _WIN32 || defined __WIN32__) && ! defined __CYGWIN__
            rc = mkdir (tmp);
#else
            rc = mkdir (tmp, 0777);
#endif
            if (rc < 0) {
                fprintf(stderr, "mkdir(%s): %s\n", tmp, strerror(errno));
                free(tmp);
                if (fatal) exit(1);
                return -errno;
            }
        }
    }
    free(tmp);
    return 0;
}


#ifdef MINGW

#ifndef PRISIZ
#error MINGW is defined, but PRISIZ is not
#endif

/* We call only v*printf(), so no infinite recursion even with rename in 
   effect here */

static const char *
subst_zu(const char *s)
{
    char *r = strdup(s);
    const char *prisiz = PRISIZ;
    const size_t l = strlen(r);
    size_t i;
    
    ASSERT_ALWAYS(strlen(prisiz) == 2);
    ASSERT_ALWAYS(r != NULL);
    for (i = 0; i + 2 < l; i++)
        if (r[i] == '%' && r[i+1] == 'z' && r[i+2] == 'u') {
            r[i+1] = prisiz[0];
            r[i+2] = prisiz[1];
        }
    return r;
}

int
printf_subst_zu (const char *format, ...)
{
  va_list ap;
  const char *subst_format;
  int r;
  
  va_start (ap, format);
  subst_format = subst_zu (format);  
  r = vprintf (subst_format, ap);
  free ((void *)subst_format);
  va_end (ap);
  return r;
}

int
fprintf_subst_zu (FILE *stream, const char *format, ...)
{
  va_list ap;
  const char *subst_format;
  int r;
  
  va_start (ap, format);
  subst_format = subst_zu (format);  
  r = vfprintf (stream, subst_format, ap);
  free ((void *)subst_format);
  va_end (ap);
  return r;
}

int
sprintf_subst_zu (char *str, const char *format, ...)
{
  va_list ap;
  const char *subst_format;
  int r;
  
  va_start (ap, format);
  subst_format = subst_zu (format);  
  r = vsprintf (str, subst_format, ap);
  free ((void *)subst_format);
  va_end (ap);
  return r;
}

int
snprintf_subst_zu (char *str, const size_t size, const char *format, ...)
{
  va_list ap;
  const char *subst_format;
  int r;
  
  va_start (ap, format);
  subst_format = subst_zu (format);  
  r = vsnprintf (str, size, subst_format, ap);
  free ((void *)subst_format);
  va_end (ap);
  return r;
}

#endif

#ifndef HAVE_ASPRINTF
/* Copied and improved from
 * http://mingw-users.1079350.n2.nabble.com/Query-regarding-offered-alternative-to-asprintf-td6329481.html
 */
int vasprintf( char **sptr, const char *fmt, va_list argv )
{
    int wanted = vsnprintf( *sptr = NULL, 0, fmt, argv );
    if (wanted<0)
        return -1;
    *sptr = malloc(1 + wanted);
    if (!*sptr)
        return -1;
    /* MinGW (the primary user of this code) can't grok %zu, so we have
     * to rewrite the format */
    const char *subst_format;
    subst_format = subst_zu (format);  
    int rc = vsnprintf(*sptr, 1+wanted, subst_zu, argv );
    free ((void *)subst_format);
    return rc;
}

int asprintf( char **sptr, const char *fmt, ... )
{
    int retval;
    va_list argv;
    va_start(argv, fmt);
    retval = vasprintf(sptr, fmt, argv);
    va_end(argv);
    return retval;
}
#endif  /* HAVE_ASPRINTF */


