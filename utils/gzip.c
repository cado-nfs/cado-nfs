#include "cado.h"
#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef HAVE_WAIT_H
#include <sys/wait.h>
#endif
#include <errno.h>


#include "portability.h"
#include "macros.h"
#include "gzip.h"
#include "misc.h"
#include "cado_popen.h"

struct suffix_handler {
    const char * suffix;
    const char * pfmt_in;
    const char * pfmt_out;
};

static char antebuffer[PATH_MAX];	/* "directory/antebuffer" or "cat" */
static int antebuffer_buffer_size = 24; /* default value 2^24 = 16 Mo */

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
    { ".gz", "gzip -dc %s", "gzip -c --fast > %s", },
    { ".bz2", "bzip2 -dc %s", "bzip2 -c --fast > %s", },
    { ".lzma", "lzma -dc  %s", "lzma -c -0 > %s", },
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

void get_suffix_from_filename (char *s, char const **sfx)
{
  const struct suffix_handler * r = supported_compression_formats;
  for( ; r->suffix ; r++)
  {
    if (has_suffix(s, r->suffix))
    {
      *sfx = r->suffix;
      return;
    }
  }

  /* If we arrive here, it's because "" is not among the suffixes */
  abort();
  return;
}

static int try_antebuffer_path()
{
    int rc = access(antebuffer, X_OK);
    if (rc >= 0) {
        fprintf(stderr, "antebuffer set to %s\n", antebuffer);
        return 1;
    }
    fprintf(stderr, "access to %s: %s\n", antebuffer, strerror(errno));
    *antebuffer = 0;
    return 0;
}

int set_antebuffer_path (const char *executable_filename, const char *path_antebuffer)
{
  *antebuffer = 0;
  /* First, if we have path_antebuffer, we must have antebuffer or error */
  if (path_antebuffer) {
      struct stat sbuf[1];
      int rc = stat(path_antebuffer, sbuf);
      if (rc < 0) {
          fprintf(stderr, "%s: path_antebuffer=\"%s\" access error: %s\n",
                  __func__, path_antebuffer, strerror(errno));
      } else {
          /* Older versions had path_antebuffer be a directory. We still
           * support this, but only as a compatibility measure. */
          if (S_ISDIR(sbuf->st_mode)) {
#ifdef EXECUTABLE_SUFFIX
              snprintf(antebuffer, PATH_MAX, "%s/antebuffer" EXECUTABLE_SUFFIX, path_antebuffer);
#else
              snprintf(antebuffer, PATH_MAX, "%s/antebuffer", path_antebuffer);
#endif
          } else {
              strncpy(antebuffer, path_antebuffer, PATH_MAX);
          }
          if (try_antebuffer_path()) return 1;
      }
  }
  /* Second option: if we failed for any reason, and if $0 was given to
   * us, use that as a potential fallback */
  if (executable_filename) {
      char dummy[PATH_MAX];
      char dummy2[PATH_MAX];
      char * slash = strrchr(executable_filename, '/');
      if (slash) {
          int len = MIN(PATH_MAX - 1, slash - executable_filename);
          strncpy(dummy, executable_filename, len);
          dummy[len]='\0';
      } else {
          dummy[0]='.';
          dummy[1]='\0';
      }
#ifdef EXECUTABLE_SUFFIX
      snprintf(dummy2, PATH_MAX, "%s/../utils/antebuffer" EXECUTABLE_SUFFIX, dummy);
#else
      snprintf(dummy2, PATH_MAX, "%s/../utils/antebuffer", dummy);
#endif
      if (realpath(dummy2, antebuffer) && try_antebuffer_path())
          return 1;
  }
  /* Third option: walk $PATH */
  if ((path_resolve("antebuffer", antebuffer)) != NULL && try_antebuffer_path()) {
      return 1;
  }
  *antebuffer = 0;
  fprintf(stderr, "No antebuffer configured\n");
  *antebuffer='\0';
  return 0;
}

/* Return a list of unix commands to _read_ a set of files. Consecutive
 * files sharing the same decompression mechanism are grouped into a
 * single command line.
 *
 * antebuffer file1.gz file2.gz file3.gz | gzip -dc
 * antebuffer file4.bz2 | gzip -dc
 * antebuffer file5.gz | gzip -dc
 * antebuffer file6.gz | cat    // useless use of cat, should be fixed.
 *
 * Note that antebuffer may also not be defined. In that case, the
 * simpler command formats like "gzip -dc file1.gz file2.gz file3.gz" are
 * used.
 *
 * All strings returned are meant to be passed to popen(). The return
 * value is a malloc()-ed array of malloc()-ed strings, and the caller is
 * in charge of freeing it (with filelist_clear, for instance).
 */
char **prepare_grouped_command_lines(char **list_of_files)
{
    const struct suffix_handler *r = supported_compression_formats;
    char ** new_commands = NULL;
    size_t n_new_commands = 0;
    
    /* Allow a few bytes extra for popen's "/bin/sh" "-c" prefix */
    ASSERT_ALWAYS(get_arg_max() >= 20);
    size_t arg_max = get_arg_max() - 20;
    
    for(char ** grouphead = list_of_files ; *grouphead ; ) {
        char *cmd_prefix = NULL, *cmd_postfix = "";
        size_t prefix_len, postfix_len;
        const struct suffix_handler * this_suffix = r;
        for (; this_suffix && this_suffix->suffix; this_suffix++)
            if (has_suffix(*grouphead, this_suffix->suffix))
                break;
        ASSERT_ALWAYS(this_suffix);
        size_t filenames_total_size = 0;
        char ** grouptail;

        if (*antebuffer) {
            if (this_suffix->pfmt_in) {
                /* antebuffer 24 file1.gz file2.gz file3.gz | gzip -dc - */
                int rc = asprintf(&cmd_prefix, "%s %d ", antebuffer, antebuffer_buffer_size);
                ASSERT_ALWAYS(rc >= 0);
                char *tmp;
                rc = asprintf(&tmp, this_suffix->pfmt_in, "-");
                ASSERT_ALWAYS(rc >= 0);
                rc = asprintf(&cmd_postfix, " | %s", tmp);
                ASSERT_ALWAYS(rc >= 0);
                free(tmp);
            } else {
                /* antebuffer 24 file1.txt file2.txt file3.txt */
                /* avoid piping through cat */
                int rc = asprintf(&cmd_prefix, "%s %d ", antebuffer, antebuffer_buffer_size);
                ASSERT_ALWAYS(rc >= 0);
            }
        } else {
            if (this_suffix->pfmt_in) {
                /* gzip -dc file1.gz file2.gz file3.gz */
                int rc = asprintf(&cmd_prefix, this_suffix->pfmt_in, "");
                ASSERT_ALWAYS(rc >= 0);
            } else {
                /* cat file1.txt file2.txt file3.txt */
                /* There's potential for this to qualify as a useless use
                 * of cat, but anyway we don't expect to meet this case
                 * often.
                 */
                int rc = asprintf(&cmd_prefix, "cat ");
                ASSERT_ALWAYS(rc >= 0);
            }
        }
        prefix_len = strlen(cmd_prefix);
        postfix_len = strlen(cmd_postfix);
        
        for(grouptail = grouphead ; *grouptail ; grouptail++) {
            const struct suffix_handler * other_suffix = r;
            for (; other_suffix && other_suffix->suffix; other_suffix++)
                if (has_suffix(*grouptail, other_suffix->suffix))
                    break;
            if (other_suffix != this_suffix)
                break;
            /* Add 1 for a space */
            size_t ds = strlen(*grouptail) + 1;
            if (filenames_total_size + prefix_len + postfix_len + ds > arg_max)
                break;
            filenames_total_size += ds;
        }
        /* Now all file names referenced by pointers in the interval
         * [grouphead..grouptail[ have the same suffix. Create a new
         * command for unpacking them.
         */
        new_commands = realloc(new_commands, ++n_new_commands * sizeof(char*));

        /* intermediary string for the list of file names */
        char * tmp = malloc(filenames_total_size + 1);
        size_t k = 0;
        for(char ** g = grouphead ; g != grouptail ; g++) {
            k += snprintf(tmp + k, filenames_total_size + 1 - k, "%s ", *g);
        }
        tmp[k-1]='\0';  /* turn final space to a null byte */
        filenames_total_size--; /* and adjust filenames_total_size for deleted space */
            
        char * cmd;
        int rc;

        rc = asprintf(&cmd, "%s%s%s", cmd_prefix, tmp, cmd_postfix);
        ASSERT_ALWAYS(rc >= 0);
        ASSERT_ALWAYS(strlen(cmd) <= arg_max);
        ASSERT_ALWAYS(strlen(cmd) == filenames_total_size + prefix_len + postfix_len);
        new_commands[n_new_commands-1] = cmd;
        free(tmp);
        free(cmd_prefix);
        if (strcmp(cmd_postfix, ""))
          free(cmd_postfix);
        grouphead = grouptail;
    }
    new_commands = realloc(new_commands, ++n_new_commands * sizeof(char*));
    new_commands[n_new_commands-1] = NULL;
    return new_commands;
}

FILE*
fopen_maybe_compressed2 (const char * name, const char * mode, int* p_pipeflag, char const ** suf)
{
    const struct suffix_handler * r = supported_compression_formats;
    FILE * f;

    if (strchr(mode, 'r') && access(name, R_OK) != 0)
        return NULL;

    for( ; r->suffix ; r++) {
        if (!has_suffix(name, r->suffix)) continue;
        if (suf) *suf = r->suffix;
        char * command = NULL;
        if (strchr(mode, 'r') && r->pfmt_in) {
            int ret = asprintf(&command, r->pfmt_in, name);
            ASSERT_ALWAYS(ret >= 0);
        } else if (strchr(mode, 'w') && r->pfmt_out) {
            int ret = asprintf(&command, r->pfmt_out, name);
            ASSERT_ALWAYS(ret >= 0);
        }

        if (command) {
          /* apparently popen() under Linux does not accept the 'b' modifier */
            char pmode[2] = "x";
            pmode[0] = mode[0];
            f = cado_popen(command, pmode);
            if (p_pipeflag) *p_pipeflag = 1;
#ifdef F_SETPIPE_SZxxx
            /* The pipe capacity is 2^16 by default; we can increase it,
             * but it does not seem to make a difference, thus we don't
             * change it by default (patch from Alain Filbois). */
            fcntl (fileno (f), F_SETPIPE_SZ, 1UL << 20);
#endif
            free(command);
        } else {
            f = fopen(name, mode);
            if (p_pipeflag) *p_pipeflag = 0;
        }
        return f;
    }
    /* If we arrive here, it's because "" is not among the suffixes */
    abort();
    return NULL;
}


FILE*
fopen_maybe_compressed (const char * name, const char * mode)
{
    return fopen_maybe_compressed2(name, mode, NULL, NULL);
}

#ifdef  HAVE_GETRUSAGE
int
fclose_maybe_compressed2 (FILE * f, const char * name, struct rusage * rr)
#else
/* if we don't even have getrusage, then no fclose_maybe_compressed2 is
 * exposed. Yet, we use one as a code shortcut
 */
static int
fclose_maybe_compressed2 (FILE * f, const char * name, void * rr MAYBE_UNUSED)
#endif
{
    const struct suffix_handler * r = supported_compression_formats;

    for( ; r->suffix ; r++) {
        if (!has_suffix(name, r->suffix)) continue;
        /* It doesn't really make sense to imagine that one of these two
         * may exist and not the other */
        ASSERT_ALWAYS((r->pfmt_out == NULL) == (r->pfmt_in == NULL));
        if (r->pfmt_in || r->pfmt_out) {
            int status;
#ifdef  HAVE_GETRUSAGE
            if (rr)
                status = cado_pclose2(f, rr);
            else
#endif
                status = cado_pclose(f);
#if defined(WIFEXITED) && defined(WEXITSTATUS)
            /* Unless child process finished normally and with exit status 0,
               we return an error */
            if (status == -1 || !WIFEXITED(status) || WEXITSTATUS(status) != 0)
                return EOF;
#else
            /* What do under MinGW? -1 definitely means an error, but how do
               we parse the other possible status codes? */
            return (status == -1) ? EOF : 0;
#endif
            return 0;
        } else {
#ifdef  HAVE_GETRUSAGE
            if (rr) memset(rr, 0, sizeof(*rr));
#endif
            return fclose(f);
        }
    }
    /* If we arrive here, it's because "" is not among the suffixes */
    abort();
    return EOF;
}

int
fclose_maybe_compressed (FILE * f, const char * name)
{
    return fclose_maybe_compressed2(f, name, NULL);
}
