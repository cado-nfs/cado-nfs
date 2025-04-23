#ifndef CADO_ROPT_IO_H
#define CADO_ROPT_IO_H

#include <stdio.h> // FILE
#include <gmp.h>        // mpz_t
#include "ropt_str.h"    // ropt_param_t

/* -- declarations -- */

#ifdef __cplusplus
extern "C" {
#endif

int cachesize_guess ( int ); // from utils

void ropt_L1_cachesize ();


/* ropt on polys in formats cado or msieve or from stdin */
void ropt_on_stdin ( ropt_param_ptr param );

void ropt_on_cadopoly ( FILE *file,
                        ropt_param_ptr param );

void ropt_on_msievepoly ( FILE *file,
                          ropt_param_ptr param );

double print_poly_info_short ( mpz_poly_srcptr f,
                               mpz_poly_srcptr g,
                               mpz_srcptr N );

#if 0
/* parse stage 2 parameters from argv */
void ropt_parse_param ( int argc,
                        char **argv,
                        ropt_param_t param );
#endif


#ifdef __cplusplus
}
#endif

#endif
