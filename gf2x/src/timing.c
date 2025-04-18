/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009, 2010, 2013, 2015, 2016
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

#include "gf2x/gf2x-config.h"

#include <sys/types.h>
#include <time.h>
#include <stdint.h>
#include "timing.h"

uint64_t microseconds()
{
  if (CLOCKS_PER_SEC < 1000000)
    return clock () * (1000000 / CLOCKS_PER_SEC);
  else
    return clock () / (CLOCKS_PER_SEC / 1000000);
}

/* returns the clock() resolution in seconds */
double
clock_resolution ()
{
  clock_t c0, c1, c2;
  int i, iter = 100; /* with 100 iterations we expect an accuracy of 1% */

  c0 = c1 = clock ();
  for (i = 0; i < iter; i++)
    {
      do { c2 = clock (); } while (c2 == c1);
      c1 = c2;
    }
  return (double) (c2 - c0) / iter;
}
