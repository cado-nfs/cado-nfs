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

#ifndef TUNING_COMMON_H_
#define TUNING_COMMON_H_

extern double mulstep;
extern FILE * rp;

extern double MINTIME;

#define TIME(x, i)				\
  { long j, k = 1;				\
    double s0 = seconds ();			\
    do {					\
      for (j = 0; j < k; j++) (i);		\
      k = 2 * k;				\
      x = seconds () - s0;			\
    } while (x < MINTIME);			\
    (x) = (x) / (k - 1);			\
  }

#ifdef __cplusplus
extern "C" {
#endif

extern void set_clock_resolution ();
extern void random_wordstring(unsigned long * a, long n);

extern void check(const unsigned long *a, long m,
	   const unsigned long *b, long n,
	   const char * cname, const unsigned long *c,
           const char * dname, const unsigned long *d);

extern void set_tuning_output();
extern void close_tuning_output();
extern int handle_tuning_mulstep(int * p_argc, char *** p_argv);
extern int handle_tuning_outfile(int * p_argc, char *** p_argv);
#ifdef __cplusplus
}
#endif

#endif	/* TUNING_COMMON_H_ */
