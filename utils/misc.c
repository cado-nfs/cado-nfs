#include "cado.h"       /* feature macros, no includes */
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdint.h>
#include <errno.h>
#include <limits.h>
#include <unistd.h>
#include <ctype.h>
#ifdef HAVE_LINUX_BINFMTS_H
/* linux/binfmts.h defines MAX_ARG_STRLEN in terms of PAGE_SIZE, but does not
   include a header where PAGE_SIZE is defined, so we include sys/user.h
   as well. Alas, on some systems, PAGE_SIZE is not defined even there. */
#include <sys/user.h>
#include <linux/binfmts.h>
#endif
/* For MinGW Build */
#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#endif

#include "macros.h"
#include "portability.h"
#include "misc.h"
#include "typecast.h"

/* Wrapper around sysconf(ARG_MAX) that deals with availability of sysconf()
   and additional constraints on command line length */
long get_arg_max(void)
{
  long arg_max;
#ifdef HAVE_SYSCONF
  arg_max = sysconf (_SC_ARG_MAX);
#elif defined(ARG_MAX)
  /* Use value from limits.h */
  arg_max = ARG_MAX;
#else
  /* POSIX requires ARG_MAX >= 4096, and all but prehistoric systems allow
     at least as much */
  arg_max = 4096;
#endif
  /* Linux since 2.6.23 does not allow more than MAX_ARG_STRLEN characers in a
     single word on the command line. Since we need to be able to run
     "sh" "-c" "actual_command", this limit is effectively the limit on the
     command line length. */
#ifdef MAX_ARG_STRLEN
  /* MAX_ARG_STRLEN may be defined in terms of PAGE_SIZE, but PAGE_SIZE may
     not actually be defined in any header.  */
#ifndef PAGE_SIZE
  const unsigned int PAGE_SIZE = pagesize();
#endif
  if ((size_t) arg_max > (size_t) MAX_ARG_STRLEN)
    arg_max = MAX_ARG_STRLEN;
#endif
  return arg_max;
}

int has_suffix(const char * path, const char * sfx)
{
    unsigned int lp = strlen(path);
    unsigned int ls = strlen(sfx);
    if (lp < ls) return 0;
    return strcmp(path + lp - ls, sfx) == 0;
}

// given a path to a file (prefix), and a suffix called (what), returns:
// - if the ext parameter is NULL, return (prefix).(what) ;
// - if ext is non-null AND (ext) is already a suffix of (prefix), say
//   we have (prefix)=(prefix0)(ext), then we return (prefix0).(what)(ext)
// - if ext is non-null AND (ext) is NOT a suffix of (prefix), 
//   we return (prefix).(what)(ext)
// In all cases the returned string is malloced, and must be freed by the
// caller later.
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


static void chomp(char *s) {
    char *p;
    if (s && (p = strrchr(s, '\n')) != NULL)
        *p = '\0';
}


#ifndef HAVE_STRLCPY
size_t
strlcpy(char *dst, const char *src, const size_t size)
{
  strncpy (dst, src, size); /* Copy at most 'size' bytes from src to dst;
                               if strlen(src) < size, then dst is null-
                               terminated, otherwise it may not be */
  if (size > 0)
      dst[size - 1] = '\0'; /* Guarantee null-termination; thus 
                               strlen(dst) < size */
  return strlen(src);
}
#endif

#ifndef HAVE_STRLCAT
size_t
strlcat(char *dst, const char *src, const size_t size)
{
  const size_t dst_len = strnlen(dst, size); /* 0 <= dst_len <= size */
  const size_t src_len = strlen(src);

  /* From man page: Note however, that if strlcat() traverses size characters
     without finding a NUL, the length of the string is considered to be size
     and the destination string will not be NUL-terminated (since there was 
     no space for the NUL). */
  if (dst_len == size)
      return dst_len + src_len;

  /* Here, 0 <= dst_len < size, thus no underflow */
  strncpy(dst + dst_len, src, size - dst_len - 1);
  
  /* If dst_len + src_len < size, then string is 0-terminated. Otherwise
     we need to put '\0' in dst[size-1] to truncate the string to length
     size-1. */
  dst[size-1] = '\0';
  
  return dst_len + src_len;
}
#endif

/* Return a NULL-terminated list of file names read from filename.
   Empty lines and comment lines (starting with '#') are skipped.
   If basepath != NULL, it is used as path before each read filename
*/
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
    if (!filelist) return;
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
            /* Test if it's an MSDOS drive specifier */
            if (strlen(tmp) == 2 && (isupper(tmp[0]) || islower(tmp[0])) && tmp[1] == ':')
              continue;
            rc = mkdir (tmp);
#else
            rc = mkdir (tmp, 0777);
#endif
            /* We have an obvious race condition between the check above
             * and the mkdir here. So failing with EEXIST can be a
             * legitimate event */
            if (rc < 0 && errno != EEXIST) {
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

char * path_resolve(const char * progname, char * resolved)
{
  const char * path = getenv("PATH");
  if (!path) return 0;
  const char * next_path;
  for( ; *path ; path = next_path) {
      next_path = strchr(path, ':');
      char * segment;
      if (next_path) {
          segment = strndup(path, next_path - path);
          next_path++;
      } else {
          segment = strdup(path);
          next_path = path + strlen(path);
      }
      char dummy2[PATH_MAX];
#ifdef EXECUTABLE_SUFFIX
      snprintf(dummy2, PATH_MAX, "%s/%s" EXECUTABLE_SUFFIX, segment, progname);
#else
      snprintf(dummy2, PATH_MAX, "%s/%s", segment, progname);
#endif
      free(segment);
      if (realpath(dummy2, resolved))
          return resolved;
  }
  return NULL;
}

//  trivial utility
const char *size_disp_fine(size_t s, char buf[16], double cutoff)
{
    char *prefixes = "bkMGT";
    double ds = s;
    const char *px = prefixes;
    for (; px[1] && ds > cutoff;) {
	ds /= 1024.0;
	px++;
    }
    snprintf(buf, 10, "%.1f%c", ds, *px);
    return buf;
}
const char *size_disp(size_t s, char buf[16])
{
    return size_disp_fine(s, buf, 500.0);
}

/* strtoul(), but with const char ** for second argument.
   Otherwise it's not possible to do, e.g., strtoul(p, &p, 10) when p is
   of type const char *
*/
unsigned long int
strtoul_const(const char *nptr, const char **endptr, const int base)
{
  char *end;
  unsigned long r;
  r = strtoul(nptr, &end, base);
  *endptr = end;
  return r;
}

unsigned long long int
strtoull_const(const char *nptr, const char **endptr, const int base)
{
  char *end;
  unsigned long long r;
  r = strtoull(nptr, &end, base);
  *endptr = end;
  return r;
}

