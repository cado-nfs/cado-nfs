/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2013, 2015
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

   This program is free software; you can redistribute it and/or modify it
   under the terms of either:
    - If the archive contains a file named toom-gpl.c (not a trivial
    placeholder), the GNU General Public License as published by the Free
    Software Foundation; either version 3 of the License, or (at your
    option) any later version.
    - If the archive contains a file named toom-gpl.c which is a trivial
    placeholder, the GNU Lesser General Public License as published by
    the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the license text for more details.
   
   You should have received a copy of the GNU General Public License as
   well as the GNU Lesser General Public License along with this program;
   see the files COPYING and COPYING.LIB.  If not, write to the Free
   Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
   02110-1301, USA.
*/

#define _DEFAULT_SOURCE /* _BSD_SOURCE is deprecated */
#define _XOPEN_SOURCE 500 /* for random() */
#define _POSIX_C_SOURCE 200112L
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include "tuning-common.h"

double mulstep = 1.0;
FILE * rp;
const char * outfile = NULL;
double MINTIME = 0.5; /* timer resolution */

/* sets the clock() resolution in seconds */
void
set_clock_resolution ()
{
  clock_t c0, c1, c2;
  int i, iter = 1000; /* with 1000 iterations we expect an accuracy of 0.1% */
  double resolution;

  c0 = c1 = clock ();
  for (i = 0; i < iter; i++)
    {
      do { c2 = clock (); } while (c2 == c1);
      c1 = c2;
    }
  resolution = (double) (c2 - c0) / (double) CLOCKS_PER_SEC;
  MINTIME = (resolution > 0.5) ? 0.5 : resolution;
  printf ("Using MINTIME = %.2es with clock() resolution of %.2es\n",
          MINTIME, resolution / (double) iter);
}

void random_wordstring(unsigned long *a, long n)
{
    long i;
    for (i = 0; i < n; i++)
      {
        /* random () returns a value between 0 and RAND_MAX = 2^31-1 */
	a[i] = random () | (random () << 31);
        if (sizeof (long) > sizeof (int))
          a[i] = (a[1] << 31) | random ();
      }
}

void dump(const unsigned long *a, long m, const unsigned long *b, long n,
	  const unsigned long *c, const unsigned long *d)
{
    printf("failed:=[");
    printf("[");
    int j;
    for (j = 0; j < m; j++) {
	if (j)
	    printf(", ");
	printf("%lu", a[j]);
    }
    printf("],");
    printf("[");
    for (j = 0; j < n; j++) {
	if (j)
	    printf(", ");
	printf("%lu", b[j]);
    }
    printf("],");
    printf("[");
    for (j = 0; j < m + n; j++) {
	if (j)
	    printf(", ");
	printf("%lu", c[j]);
    }
    printf("],");
    printf("[");
    for (j = 0; j < m + n; j++) {
	if (j)
	    printf(", ");
	printf("%lu", d[j]);
    }
    printf("]");
    printf("];\n");
}

void check(const unsigned long *a, long m,
	   const unsigned long *b, long n,
	   const char * cname, const unsigned long *c,
           const char * dname, const unsigned long *d)
{
    long i = 0;
    int error = 0;
    for(i = 0 ; i < m + n ; i++) {
        if (c[i] != d[i]) {
            if (error++ < 10)
              {
                fprintf(stderr,
                        "Error: %s and %s differ for %ldx%ld at word %ld\n",
                        cname, dname, m, n, i);
                fprintf (stderr, "expected %lx, got %lx\n", c[i], d[i]);
              }
        }
    }
    if (error)
      {
        if (m + n < 1000) {
          dump(a, m, b, n, c, d);
        }
        abort();
      }
}


void set_tuning_output()
{
    if (outfile) {
        if ((rp = fopen(outfile, "w")) == NULL) {
            fprintf(stderr, "fopen(%s): %s\n", outfile, strerror(errno));
            exit(1);
        }
    } else {
#if 0
        /* Not really useful since there is -o anyway */
        /* If file descriptor # 3 is open for writing, use it. */
        if ((rp = fdopen(3, "w")) == NULL) {
            rp = stdout;
        }
#else
        rp = stdout;
#endif
    }

    setbuf(rp, NULL);
    setbuf(stdout, NULL);
}

void close_tuning_output()
{
  if (outfile)
    fclose (rp);
}

int handle_tuning_mulstep(int * p_argc, char *** p_argv)
{
    if (strcmp((*p_argv)[0], "--step") == 0 || strcmp((*p_argv)[0], "-s") == 0) {
        (*p_argc)--,(*p_argv)++;
        if (! (*p_argc)) {
            return -1;
        }
        mulstep = atof((*p_argv)[0]);
        return 1;
    }
    return 0;
}

int handle_tuning_outfile(int * p_argc, char *** p_argv)
{
    if (strcmp((*p_argv)[0], "--output") == 0 || strcmp((*p_argv)[0], "-o") == 0) {
        (*p_argc)--,(*p_argv)++;
        if (! (*p_argc)) {
            return -1;
        }
        outfile = (*p_argv)[0];
        return 1;
    }
    return 0;
}
