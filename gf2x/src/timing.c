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

#include "gf2x/gf2x-config.h"

#include <sys/types.h>
#ifdef HAVE_SYS_RESOURCE_H
#include <sys/resource.h>
#else
#include <time.h>
#endif  /* HAVE_SYS_RESOURCE_H */
#include <stdint.h>
#include "timing.h"

uint64_t microseconds()
{
#ifdef HAVE_SYS_RESOURCE_H
        struct rusage res[1];
        getrusage(RUSAGE_SELF,res);
        uint64_t r;
        r  = (uint64_t) res->ru_utime.tv_sec;
        r *= (uint64_t) 1000000UL;
        r += (uint64_t) res->ru_utime.tv_usec;
        return r;
#else
        return clock () / (CLOCKS_PER_SEC / 1000000);
#endif
}
