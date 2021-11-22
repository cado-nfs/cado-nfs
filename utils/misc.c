#include "cado.h" // IWYU pragma: keep
#include <stdio.h>          // for asprintf, fprintf, perror, snprintf, fclose
#include <sys/stat.h>   // mkdir
#include <unistd.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
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
#include "misc.h"
#include "gmp_aux.h"    // mpz_get_uint64
#include "portability.h" // asprintf // IWYU pragma: keep

// re-entrant random with GMP.
uint64_t u64_random(gmp_randstate_t buf) {
#if ULONG_BITS == 64
    return gmp_urandomb_ui(buf, 64);
#elif ULONG_BITS == 32
    mpz_t z;
    mpz_init(z);
    mpz_urandomb(z, buf, 64);
    uint64_t r = mpz_get_uint64(z);
    mpz_clear(z);
    return r;
#endif
}

int64_t i64_random(gmp_randstate_t buf) {
    return (int64_t) u64_random(buf);
}

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
  /* Linux since 2.6.23 does not allow more than MAX_ARG_STRLEN characters in a
     single word on the command line. Since we need to be able to run
     "sh" "-c" "actual_command", this limit is effectively the limit on the
     command line length. */
#ifdef MAX_ARG_STRLEN
  /* MAX_ARG_STRLEN may be defined in terms of PAGE_SIZE, but PAGE_SIZE may
     not actually be defined in any header.  */
#if !defined(PAGE_SIZE) && defined(HAVE_SYSCONF)
  /* If we have sysconf(), we can resolve any reference to PAGE_SIZE in
     MAX_ARG_STRLEN to a variable */
  const size_t PAGE_SIZE MAYBE_UNUSED = sysconf (_SC_PAGESIZE);
  /* If we don't have sysconf(), either, we're pretty much out of options */
#endif
  if ((size_t) arg_max > (size_t) MAX_ARG_STRLEN)
    arg_max = MAX_ARG_STRLEN;
#endif
  /* as discussed on
   * https://sympa.inria.fr/sympa/arc/cado-nfs/2019-10/msg00000.html
   * one should subtract from arg_max the length of the environment, but
   * even this seems not enough, so we divide arg_max by 2 */
  return arg_max / 2;
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
            files = (char**) realloc(files, nfiles_alloc * sizeof(char*));
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
        files = (char**) realloc(files, nfiles_alloc * sizeof(char*));
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
        strncpy(tmp, dir, n+1);
        if (slash) {
            pos = slash - dir;
            tmp[pos]='\0';
        } else {
            pos = n;
        }
        struct stat sbuf[1];
        // coverity[fs_check_call]
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
    const char *prefixes = "bkMGT";
    double ds = s;
    const char *px = prefixes;
    for (; px[1] && ds > cutoff;) {
	ds /= 1024.0;
	px++;
    }
    snprintf(buf, 10, "%.2f %c%s", ds, *px, px==prefixes ? "" : "B");
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

static inline unsigned long MAYBE_UNUSED bitrev1(unsigned long a)/*{{{*/
{
    unsigned long m;
#if ULONG_BITS == 64
    /* these three should be doable with simple bswap, right ? */
    a = (a >> 32) ^ (a << 32);
    m = UINT64_C(0x0000ffff0000ffff); a = ((a >> 16) & m) ^ ((a << 16) & ~m);
    m = UINT64_C(0x00ff00ff00ff00ff); a = ((a >> 8) & m) ^ ((a << 8) & ~m);
    
    m = UINT64_C(0x0f0f0f0f0f0f0f0f); a = ((a >> 4) & m) ^ ((a << 4) & ~m);
    m = UINT64_C(0x3333333333333333); a = ((a >> 2) & m) ^ ((a << 2) & ~m);
    m = UINT64_C(0x5555555555555555); a = ((a >> 1) & m) ^ ((a << 1) & ~m);
#else
    a = (a >> 16) ^ (a << 16);
    m = 0x00ff00ffUL; a = ((a >> 8) & m) ^ ((a << 8) & ~m);

    m = 0x0f0f0f0fUL; a = ((a >> 4) & m) ^ ((a << 4) & ~m);
    m = 0x33333333UL; a = ((a >> 2) & m) ^ ((a << 2) & ~m);
    m = 0x55555555UL; a = ((a >> 1) & m) ^ ((a << 1) & ~m);
#endif
    return a;
}

extern void bit_reverse(unsigned long * dst, const unsigned long * ptr, size_t n)
{
    unsigned int i = 0;
    unsigned int j = n-1;
    for(i = 0 ; i < j ; i++, j--) {
        unsigned long lo = ptr[i];
        unsigned long hi = ptr[j];
        dst[j] = bitrev1(lo);
        dst[i] = bitrev1(hi);
    }
    if (i == j)
        dst[i] = bitrev1(ptr[i]);
}
