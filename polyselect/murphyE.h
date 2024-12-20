#ifndef MURPHYE_H_
#define MURPHYE_H_

/* Header file for murphyE.c.

Copyright 2010 Paul Zimmermann

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "cado_poly.h"
#include "gmp_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MURPHY_K 1000

double MurphyE (cado_poly_srcptr cpoly, double Bf, double Bg, double area, int K,
                unsigned long B);

#ifdef __cplusplus
}
#endif

#endif	/* MURPHYE_H_ */
