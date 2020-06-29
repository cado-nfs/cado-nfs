/* Alternate implementation of singleton_simple.cpp, in pure C, stand-alone.

   Usage: singleton_simple -col-max-index nnn -out xxx x y z.

   (there is no -col-min-index option, all ideals are taken into account).
*/

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <unistd.h>
#include <fcntl.h>
#include "macros.h"
#include "portability.h"

#define MAX_SIZE 1024

#undef ASSERT_ALWAYS
#define ASSERT_ALWAYS assert

#define AB_BASE 16 /* base of a,b */

unsigned char *T;
unsigned long nrels = 0;

int has_suffix(const char * path, const char * sfx)
{
    unsigned int lp = strlen(path);
    unsigned int ls = strlen(sfx);
    if (lp < ls) return 0;
    return strcmp(path + lp - ls, sfx) == 0;
}

struct suffix_handler {
    const char * suffix;
    const char * pfmt_in;
    const char * pfmt_out;
};

struct suffix_handler supported_compression_formats[] = {
    { ".gz", "gzip -dc %s", "gzip -c --fast > %s", },
    { ".bz2", "bzip2 -dc %s", "bzip2 -c --fast > %s", },
    { ".lzma", "lzma -dc  %s", "lzma -c -0 > %s", },
    /* These two have to be present */
    { "", NULL, NULL },
    { NULL, NULL, NULL },
};

static inline FILE * cado_popen(const char * command, const char * mode) { return popen(command, mode); }

static inline int cado_pclose(FILE * stream) { return pclose(stream); }

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

static int
fclose_maybe_compressed2 (FILE * f, const char * name)
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
                status = cado_pclose(f);
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
    return fclose_maybe_compressed2(f, name);
}

/* Pass 1: read relations and count ideal weight in T[p] (capped to 255) */
void
pass1 (int argc, char *argv[], unsigned long col_max_index)
{
  char s[MAX_SIZE], *end;
  long a;
  unsigned long b, p, lastp, nideals = 0;
  int ret, e;
  FILE *fp;

  while (argc > 1)
    {
      printf ("Pass 1: read file %s\n", argv[1]);
      fflush (stdout);
      fp = fopen_maybe_compressed (argv[1], "r");
      argc --;
      argv ++;
      while (fgets (s, MAX_SIZE, fp) != NULL)
	{
	  ASSERT (s[0] != '\n');
          /* ensure line ends with \n */
          if (s[strlen(s)-1] != '\n')
            sprintf (s + strlen(s), "\n");
	  if (s[0] != '#') /* not a comment */
	    {
	      a = strtol (s, &end, AB_BASE);
	      b = strtoul (end + 1, &end, AB_BASE);
	      lastp = ULONG_MAX;
	      while (*end != '\n')
		{
		  p = strtoul (end + 1, &end, 16);
		  if (p != lastp)
		    {
		      if (lastp != ULONG_MAX && (e & 1))
			{
			  ASSERT (lastp < col_max_index);
			  T[lastp] += (T[lastp] < 255);
			}
		      e = 1;
		      lastp = p;
		    }
		}
	      if (e & 1)
		{
		  ASSERT (lastp < col_max_index);
		  T[lastp] += (T[lastp] < 255);
		}
	      nrels ++;
	      if ((nrels & 0x7fffff) == 0)
		{
		  printf ("Pass 1: read %lu relations\n", nrels);
		  fflush (stdout);
		}
	    }
	}
      fclose (fp);
    }
  for (p = 0; p < col_max_index; p++)
    nideals += (T[p] != 0);
  printf ("Pass 1: read %lu relations on %lu ideals\n", nrels, nideals);
  fflush (stdout);
}

/* Pass 2: output relations with no singleton */
void
pass2 (int argc, char *argv[], char *out_file)
{
  char s[MAX_SIZE], *end;
  long a;
  unsigned long b, p;
  int ret, e;
  FILE *fp, *out;
  unsigned long nrels_out = 0;

  out = fopen (out_file, "w");
  nrels = 0;
  while (argc > 1)
    {
      printf ("Pass 2: read file %s\n", argv[1]);
      fflush (stdout);
      fp = fopen_maybe_compressed (argv[1], "r");
      argc --;
      argv ++;
      while (fgets (s, MAX_SIZE, fp) != NULL)
	{
	  ASSERT (s[0] != '\n');
          /* ensure line ends with \n */
          if (s[strlen(s)-1] != '\n')
            sprintf (s + strlen(s), "\n");
	  if (s[0] != '#') /* not a comment */
	    {
	      int singleton = 0;
	      a = strtol (s, &end, AB_BASE);
	      b = strtoul (end + 1, &end, AB_BASE);
	      while (*end != '\n' && singleton == 0)
		{
		  p = strtoul (end + 1, &end, 16);
		  ASSERT (T[p] != 0);
		  if (T[p] == 1)
		    singleton = 1;
		}
	      nrels ++;
	      if (singleton == 0)
		{
		  fprintf (out, s);
		  nrels_out ++;
		}
	      if ((nrels & 0x7fffff) == 0)
		{
		  printf ("Pass 2: read %lu relations\n", nrels);
		  fflush (stdout);
		}
	    }
	}
      fclose (fp);
    }
  printf ("Pass 2: read %lu relations, output %lu\n", nrels, nrels_out);
  fflush (stdout);
  fclose (out);
}

int
main (int argc, char *argv[])
{
  unsigned long col_max_index = 0;
  char *out_file = NULL;

  while (argc > 2 && argv[1][0] == '-')
    {
      if (strcmp (argv[1], "-col-max-index") == 0)
        {
          col_max_index = strtoul (argv[2], NULL, 10);
          argv += 2;
          argc -= 2;
        }
      else if (strcmp (argv[1], "-out") == 0)
        {
          out_file = argv[2];
          argv += 2;
          argc -= 2;
        }
      else
        {
          fprintf (stderr, "Error, unknown option %s\n", argv[2]);
          exit (1);
        }
    }

  if (col_max_index == 0)
    {
      fprintf (stderr, "Error, missing -col-max-index option\n");
      exit (1);
    }

  T = malloc (col_max_index);
  memset (T, 0, col_max_index);
  pass1 (argc, argv, col_max_index);
  pass2 (argc, argv, out_file);
  free (T);
  return 0;
}
