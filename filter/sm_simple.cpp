/* Shirokauer maps 

   Given a list of a,b pairs, compute the corresponding SMs.

   */

#include "cado.h"

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>

#include "macros.h"
#include "utils.h"

static void my_sm(const char *outfile, const char *infile, 
		  sm_side_info *sm_info, int nb_polys)
{
  FILE *in;
  in = fopen(infile, "r");
  if (in == NULL) {
    fprintf(stderr, "Error: could not open %s for reading\n", infile);
    exit(EXIT_FAILURE);
  }
  FILE *out = stdout;
  if (outfile != NULL) {
    out = fopen(outfile, "w");
    if (out == NULL) {
      fprintf(stderr, "Error: could not open %s for writing\n", outfile);
      exit(EXIT_FAILURE);
    }
  }

  char buf[1024];
  int maxdeg = sm_info[0]->f->deg;
  for(int side = 1; side < nb_polys; side++)
      maxdeg = MAX(maxdeg, sm_info[side]->f->deg);

  cxx_mpz_poly pol(maxdeg);
  cxx_mpz_poly smpol(maxdeg);

  while (fgets(buf, 1024, in)) {
    if (buf[0] == '#')
      continue;


    if (buf[0] != 'X') {
        // we read a relation in legacy format.
        cxx_mpz a, b;
        int ret = gmp_sscanf(buf, "%Zd,%Zd:", (mpz_ptr) a, (mpz_ptr) b);
        ASSERT_ALWAYS(ret == 2);
        mpz_poly_set_zero(pol);
        mpz_poly_setcoeff(pol, 0, a);
        mpz_neg(b,b);
        mpz_poly_setcoeff(pol, 1, b);
    } else {
        // here we parse several coefficients.
        char * p = buf+1;
        for(unsigned int i = 0 ; ; i++) {
            for( ; isspace(*p) || *p == ',' ; p++);
            cxx_mpz x;
            int c;
            if (gmp_sscanf(p, "%Zx%n", (mpz_ptr) x, &c) < 1) {
                ASSERT_ALWAYS(*p++ == ':');
                break;
            }
            mpz_poly_setcoeff(pol, i, x);
            p += c;
        }
    }

    for (int side = 0; side < nb_polys; ++side) {
      compute_sm_piecewise(smpol, pol, sm_info[side]);
      print_sm(out, smpol, sm_info[side]->nsm, sm_info[side]->f->deg);
      if (side == 0 && sm_info[0]->nsm > 0 && sm_info[1]->nsm > 0)
          fprintf(out, " ");
    }
    fprintf(out, "\n");
  }

  if (out != NULL)
    fclose(out);
  fclose(in);
}


static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "poly", "(required) poly file");
  param_list_decl_usage(pl, "inp", "(required) input file containing relations");
  param_list_decl_usage(pl, "out", "output file");
  param_list_decl_usage(pl, "ell", "(required) group order");
  verbose_decl_usage(pl);
}

static void usage (const char *argv, const char * missing, param_list pl)
{
  if (missing) {
    fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
        missing);
  }
  param_list_print_usage(pl, argv, stderr);
  exit (EXIT_FAILURE);
}

/* -------------------------------------------------------------------------- */

int main (int argc, char **argv)
{
  char *argv0 = argv[0];

  const char *polyfile = NULL;
  const char *infile = NULL;
  const char *outfile = NULL;

  param_list pl;
  cado_poly pol;
  mpz_poly_ptr F[NB_POLYS_MAX];

  mpz_t ell, ell2;
  double t0;

  /* read params */
  param_list_init(pl);
  declare_usage(pl);

  if (argc == 1)
    usage (argv[0], NULL, pl);

  argc--,argv++;
  for ( ; argc ; ) {
    if (param_list_update_cmdline (pl, &argc, &argv)) { continue; }
    fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
    usage (argv0, NULL, pl);
  }

  /* Read poly filename from command line */
  if ((polyfile = param_list_lookup_string(pl, "poly")) == NULL) {
    fprintf(stderr, "Error: parameter -poly is mandatory\n");
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }

  /* Read purged filename from command line */
  if ((infile = param_list_lookup_string(pl, "inp")) == NULL) {
    fprintf(stderr, "Error: parameter -inp is mandatory\n");
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }

  /* Read outfile filename from command line ; defaults to stdout. */
  outfile = param_list_lookup_string(pl, "out");

  /* Read ell from command line (assuming radix 10) */
  mpz_init (ell);
  if (!param_list_parse_mpz(pl, "ell", ell)) {
    fprintf(stderr, "Error: parameter -ell is mandatory\n");
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }

  /* Init polynomial */
  cado_poly_init (pol);
  cado_poly_read(pol, polyfile);
  for(int side = 0; side < pol->nb_polys; side++)
      F[side] = pol->pols[side];

  if (param_list_warn_unused(pl))
    usage (argv0, NULL, pl);
  verbose_interpret_parameters(pl);
  param_list_print_command_line (stdout, pl);

  mpz_init(ell2);
  mpz_mul(ell2, ell, ell);

  sm_side_info sm_info[NB_POLYS_MAX];

  for(int side = 0 ; side < pol->nb_polys; side++) {
    sm_side_info_init(sm_info[side], F[side], ell);
  }

  for (int side = 0; side < pol->nb_polys; side++) {
    fprintf(stdout, "\n# Polynomial on side %d:\nF[%d] = ", side, side);
    mpz_poly_fprintf(stdout, F[side]);

    printf("# SM info on side %d:\n", side);
    sm_side_info_print(stdout, sm_info[side]);

    fflush(stdout);
  }

  t0 = seconds();

  my_sm(outfile, infile, sm_info, pol->nb_polys);

  fprintf(stdout, "\n# sm completed in %2.2lf seconds\n", seconds() - t0);
  fflush(stdout);

  for(int side = 0 ; side < pol->nb_polys ; side++) {
    sm_side_info_clear(sm_info[side]);
  }

  mpz_clear(ell);
  mpz_clear(ell2);
  cado_poly_clear(pol);
  param_list_clear(pl);

  return 0;
}
