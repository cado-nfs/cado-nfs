/* Compute Murphy E-value of an imported polynomial

Copyright 2024 Paul Zimmermann and Tyler Busby

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

#include <stdio.h>
#include <stdlib.h>

#include "auxiliary.h"

#include "cado_poly.h"
#include "params.h"
#include "murphyE.h"

static int
compute_murphyE (const char *input_file, const char *output_file)
{
    cado_poly p;
    double E;
    cado_poly_init(p);
    int rc = EXIT_SUCCESS;
    if (!cado_poly_read(p, input_file)) {
	fprintf(stderr, "Error reading polynomial file %s\n", input_file);
        rc = EXIT_FAILURE;
    } else if (output_file == NULL) {
	E = MurphyE(p, bound_f, bound_g, area, MURPHY_K, get_alpha_bound());
	printf("%.5g\n", E);
    } else {
	FILE *of;
	of = fopen(output_file, "w");
	if (of == NULL) {
	    fprintf(stderr, "Error writing polynomial file %s\n",
		    output_file);
            rc = EXIT_FAILURE;
	} else {
            cado_poly_fprintf(of, "", p);
        }
    }
    cado_poly_clear(p);
    return rc;
}

// usage: score input_file <output_file>
int main(int argc, char const * argv[])
{
    param_list pl;
    param_list_init(pl);

    const char * progname = argv[0];
    const char *input_file = NULL;
    const char *output_file = NULL;

    param_list_decl_usage(pl, "Bf", "smoothness bound on f side");
    param_list_decl_usage(pl, "Bg", "smoothness bound on g side");
    param_list_decl_usage(pl, "area", "estimate of sieve are");
    argv++,argc--;

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        if (!input_file) {
            input_file = argv[0];
            argv++, argc--;
            continue;
        }
        if (!output_file) {
            output_file = argv[0];
            argv++, argc--;
            continue;
        }
        fprintf(stderr, "Unknown option: %s\n", argv[0]);
        param_list_print_usage(pl, progname, stderr);
        return EXIT_FAILURE;
    }

    param_list_parse_double(pl, "Bf", &bound_f);
    param_list_parse_double(pl, "Bg", &bound_g);
    param_list_parse_double(pl, "area", &area);

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, progname, stderr);
        return EXIT_FAILURE;
     }

    int rc = compute_murphyE(input_file, output_file);
    param_list_clear(pl);
    return rc;
}
