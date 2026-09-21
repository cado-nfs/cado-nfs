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
#include "polyselect_norms.hpp"
#include "cado_poly.hpp"
#include "macros.h"
#include "cado_main.hpp"
#include "utils_cxx.hpp"

static void
compute_skewness (const char *input_file, const char *output_file)
{
  cxx_cado_poly p;
  if (!p.read(input_file))
    {
      throw cado::error("Error reading polynomial file {}", input_file);
    }
  p.skew = L2_combined_skewness2 (p[RAT_SIDE], p[ALG_SIDE]);
  if (output_file == NULL) {
    printf("%g\n", p.skew);
  } else {
    FILE *of;
    of = fopen (output_file, "w");
    if (of == NULL)
      {
        throw cado::error("Error writing polynomial file {}", output_file);
      }
    p.fprintf (of);
  }
} 

// usage: skewness input_file output_file
static int main_(int argc, char const * argv[]);

int main(int argc, char const * argv[])
{
    return cado::main_wrapper(main_, argc, argv);
}

static int main_(int argc, char const * argv[])
{
  ASSERT_ALWAYS (argc == 2 || argc == 3);
  const char *input_file = argv[1];
  const char *output_file = (argc == 3) ? argv[2] : NULL;
  compute_skewness (input_file, output_file);

  return EXIT_SUCCESS;
}
