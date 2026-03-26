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

#include <cstdio>
#include <cstdlib>

#include "auxiliary.hpp"

#include "cado_poly.hpp"
#include "params.h"
#include "murphyE.hpp"
#include "polyselect_norms.h"
#include "polyselect_alpha.h"

int fullscore = 0;

static int
compute_murphyE (const char *input_file, const char *output_file)
{
    cxx_cado_poly p;
    double E;
    int rc = EXIT_SUCCESS;

    FILE * of = NULL;

    if (!p.read(input_file)) {
	fprintf(stderr, "Error reading polynomial file %s\n", input_file);
        rc = EXIT_FAILURE;
    } else if (output_file == NULL) {
	E = MurphyE(p, bound_f, bound_g, area, MURPHY_K, get_alpha_bound());
        of = stdout;
        if (!fullscore)
            printf("%.5g\n", E);
    } else {
	FILE *of;
	of = fopen(output_file, "w");
	if (of == NULL) {
	    fprintf(stderr, "Error writing polynomial file %s\n",
		    output_file);
            rc = EXIT_FAILURE;
	} else {
            p.fprintf(of);
        }
    }

    if (fullscore) {
        double skew = p.skew;
        // note that there's some inconsistency regarding whether we take
        // the skewness of just f, or the skewness of f combined with the
        // skewness of g.
        // double skew = L2_skewness(p[1]);
        double logmu = L2_lognorm(p[1], skew);
        unsigned int nroots = mpz_poly_number_of_real_roots(p[1]);
        // double comb_skew = L2_combined_skewness2 (p[0], p[1]);

        // p.skew = L2_skewness(p[1]);
        p.skew = L2_combined_skewness2 (p[0], p[1]);

        double exp_E = logmu + expected_rotation_gain(p[1], p[0]);

        double E = MurphyE (p, bound_f, bound_g, area, MURPHY_K, get_alpha_bound());

        const char * prefix = "# ";

        if (output_file == NULL) {
            fprintf(of, "poly0="); mpz_poly_print_raw(p[0]);
            fprintf(of, "poly1="); mpz_poly_print_raw(p[1]);
        }
        fprintf(of, "%sskew_f=%1.2f, skew_fg=%1.2f\n",
                prefix,
                L2_skewness(p[1]),
                L2_combined_skewness2 (p[0], p[1]));
        fprintf(of, "%sexp_E %1.2f, lognorm %1.2f, skew %1.2f, %d rroots\n",
                prefix, exp_E, logmu, skew, nroots);
        fprintf(of, "%sMurphyE(Bf=%.3e, Bf=%.3e, area=%.3e)=%1.2e\n",
                prefix, bound_f, bound_g, area, E);

    }

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
    param_list_decl_usage(pl, "full", "give full score");
    param_list_configure_switch(pl, "--full", &fullscore);
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

    if (param_list_warn_unused(pl) || !input_file) {
        param_list_print_usage(pl, progname, stderr);
        return EXIT_FAILURE;
     }

    int rc = compute_murphyE(input_file, output_file);

    param_list_clear(pl);
    return rc;
}
