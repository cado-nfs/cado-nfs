/* This file is part of the gf2x library.

   Copyright 2007, 2008, 2009
   Richard Brent, Pierrick Gaudry, Emmanuel Thome', Paul Zimmermann

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation; either version 2.1 of the License, or (at
   your option) any later version.
   
   This program is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
   License for more details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with CADO-NFS; see the file COPYING.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
   Boston, MA 02110-1301, USA.
*/
#include <sys/types.h>    /* for cputime */
#ifdef HAVE_RESOURCE_H
#include <sys/resource.h> /* for cputime */
#endif
#include <stdint.h>
#include "timing.h"

uint64_t microseconds()
{
        struct rusage res[1];
        getrusage(RUSAGE_SELF,res);
        uint64_t r;
        r  = (uint64_t) res->ru_utime.tv_sec;
        r *= (uint64_t) 1000000UL;
        r += (uint64_t) res->ru_utime.tv_usec;
        return r;
}
