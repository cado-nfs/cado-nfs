/* Compute skewness of an imported polynomial

Copyright 2024 Paul Zimmermann

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

#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include "auxiliary.hpp"
#include "polyselect_norms.h"
#include "cado_poly.hpp"
#include "macros.h"

static void
compute_skewness (const char *input_file, const char *output_file)
{
  cxx_cado_poly p;
  if (!p.read(input_file))
    {
      fprintf (stderr, "Error reading polynomial file %s\n", input_file);
      exit (EXIT_FAILURE);
    }
  p.skew = L2_combined_skewness2 (p[RAT_SIDE], p[ALG_SIDE]);
  if (output_file == NULL) {
    printf("%g\n", p.skew);
  } else {
    FILE *of;
    of = fopen (output_file, "w");
    if (of == NULL)
      {
        fprintf (stderr, "Error writing polynomial file %s\n", output_file);
        exit (EXIT_FAILURE);
      }
    p.fprintf (of);
  }
} 

// usage: skewness input_file output_file
int main(int argc, char const * argv[])
{
  ASSERT_ALWAYS (argc == 2 || argc == 3);
  const char *input_file = argv[1];
  const char *output_file = (argc == 3) ? argv[2] : NULL;
  compute_skewness (input_file, output_file);
}
