/*
  Size-optimize all polynomial in a file.
  File must contains polynomials in CADO format and polynomials must be
  separated by a newline.
*/

#include "cado.h" // IWYU pragma: keep
#include <float.h> // for DBL_MAX
#include <stdio.h>
#include <stdlib.h>
#include "auxiliary.h"
#include "cado_poly.h"
#include "mpz_poly.h"
#include "params.h"
#include "size_optimization.h"
#include "usp.h"        // numberOfRealRoots
#include "verbose.h"             // verbose_decl_usage
#include "polyselect_norms.h"
#include "polyselect_alpha.h"


static void
declare_usage(param_list pl)
{
  char str[200];
  param_list_decl_usage(pl, "inputpolys", "size-optimized the polynomials "
                                          "given in this file");
  snprintf (str, 200, "size-optimization effort (default %d)", SOPT_DEFAULT_EFFORT);
  param_list_decl_usage(pl, "sopteffort", str);
  param_list_decl_usage(pl, "v", "verbose mode");
  param_list_decl_usage(pl, "translation-only", "(switch) do not use rotations");
  verbose_decl_usage(pl);
}

static void
usage (const char *argv, param_list pl)
{
  param_list_print_usage (pl, argv, stderr);
  param_list_clear (pl);
  exit (EXIT_FAILURE);
}

/*
   main function 
*/
int main (int argc, char **argv)
{
  char *argv0 = argv[0];
  int use_only_translation = 0; /* only use translation in optimization */
  unsigned int sopt_effort = SOPT_DEFAULT_EFFORT; /* size optimization effort */
  int verbose = 0;
  FILE *polys_file = NULL;
  const char *polys_filename = NULL;
  unsigned int nb_input_polys = 0; /* number of input polynomials */
  cado_poly cpoly;
  /* For statistics */
  double ave_raw_lognorm = 0.0, ave_raw_alpha = 0.0;
  double min_raw_lognorm = DBL_MAX, max_raw_lognorm = 0.0;
  double ave_sopt_lognorm = 0.0, ave_sopt_alpha = 0.0;
  double min_sopt_lognorm = DBL_MAX, max_sopt_lognorm = 0.0;


  /* read params */
  param_list pl;
  param_list_init (pl);

  declare_usage(pl);

  param_list_configure_switch (pl, "-v", &verbose);
  param_list_configure_switch (pl, "-translation-only", &use_only_translation);

  if (argc == 1)
    usage (argv0, pl);
  
  
  argv++, argc--;
  for ( ; argc; ) {
    if (param_list_update_cmdline (pl, &argc, &argv)) continue;
    fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
    usage (argv0, pl);
  }

  /* parse poly filename */
  polys_filename = param_list_lookup_string (pl, "inputpolys");

  /* parse size optimization effort that passed to size_optimization */
  param_list_parse_uint (pl, "sopteffort", &sopt_effort);

  if (param_list_warn_unused(pl))
    usage (argv0, pl);

  /* print out commands */
  verbose_interpret_parameters(pl);
  param_list_print_command_line (stdout, pl);

  /* Open file containing polynomials. */
  printf ("# Reading polynomials from %s\n", polys_filename);
  polys_file = fopen(polys_filename, "r");
  if (polys_file == NULL)
  {
    perror("Could not open file");
    exit(EXIT_FAILURE);
  }

  cado_poly_init (cpoly);

  /* Main loop: read all polynomials from file and do size-optimization. */
  while (cado_poly_read_next_poly_from_stream (cpoly, polys_file))
  {
    unsigned int nrroots;
    double lognorm, alpha, alpha_proj, exp_E;

    printf ("\n### Input raw polynomial (%u) ###\n", nb_input_polys);
    cpoly->skew = L2_skewness (cpoly->pols[ALG_SIDE], SKEWNESS_DEFAULT_PREC);
    nrroots = numberOfRealRoots ((const mpz_t *) cpoly->pols[ALG_SIDE]->coeff, cpoly->pols[ALG_SIDE]->deg, 0, 0, NULL);
    lognorm = L2_lognorm (cpoly->pols[ALG_SIDE], cpoly->skew);
    alpha = get_alpha (cpoly->pols[ALG_SIDE], get_alpha_bound ());
    alpha_proj = get_alpha_projective (cpoly->pols[ALG_SIDE], get_alpha_bound ());
    exp_E = lognorm
      + expected_rotation_gain (cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE]);
    cado_poly_fprintf (stdout, cpoly, "# ");
    cado_poly_fprintf_info (stdout, lognorm, exp_E, alpha, alpha_proj, nrroots,
                            "# ");

    ave_raw_lognorm += lognorm;
    min_raw_lognorm = (lognorm < min_raw_lognorm) ? lognorm : min_raw_lognorm;
    max_raw_lognorm = (lognorm > max_raw_lognorm) ? lognorm : max_raw_lognorm;
    ave_raw_alpha += alpha;

    /* Size-optimize */
    if (use_only_translation)
      sopt_local_descent (cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE], cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE], 
                                      1, -1, SOPT_DEFAULT_MAX_STEPS, verbose);
    else
      size_optimization (cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE], cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE],
                                                        sopt_effort, verbose);

    printf ("### Size-optimized polynomial (%u) ###\n", nb_input_polys);
    cpoly->skew = L2_skewness (cpoly->pols[ALG_SIDE], SKEWNESS_DEFAULT_PREC);
    nrroots = numberOfRealRoots ((const mpz_t *) cpoly->pols[ALG_SIDE]->coeff, cpoly->pols[ALG_SIDE]->deg, 0, 0, NULL);
    lognorm = L2_lognorm (cpoly->pols[ALG_SIDE], cpoly->skew);
    alpha = get_alpha (cpoly->pols[ALG_SIDE], get_alpha_bound ());
    alpha_proj = get_alpha_projective (cpoly->pols[ALG_SIDE], get_alpha_bound ());
    exp_E = lognorm
      + expected_rotation_gain (cpoly->pols[ALG_SIDE], cpoly->pols[RAT_SIDE]);
    cado_poly_fprintf (stdout, cpoly, NULL);
    cado_poly_fprintf_info (stdout, lognorm, exp_E, alpha, alpha_proj, nrroots,
                            NULL);


    ave_sopt_lognorm += lognorm;
    min_sopt_lognorm = (lognorm < min_sopt_lognorm) ? lognorm : min_sopt_lognorm;
    max_sopt_lognorm = (lognorm > max_sopt_lognorm) ? lognorm : max_sopt_lognorm;
    ave_sopt_alpha += alpha;

    nb_input_polys++;
  }
  /* Did we stop at the end of the file or was there an error. */
  if (!feof (polys_file))
  {
    fprintf (stderr, "Error while reading file %s.\n", polys_filename);
    abort ();
  }

  printf("\n# Number of polynomials read: %u\n", nb_input_polys);
  printf("# Average input lognorm: %3.3f\n", ave_raw_lognorm / nb_input_polys);
  printf("# Minimum input lognorm: %3.3f\n", min_raw_lognorm);
  printf("# Maximum input lognorm: %3.3f\n", max_raw_lognorm);
  printf("# Average input alpha value: %3.3f\n", ave_raw_alpha / nb_input_polys);
  printf("# Average sopt lognorm: %3.3f\n", ave_sopt_lognorm / nb_input_polys);
  printf("# Minimum sopt lognorm: %3.3f\n", min_sopt_lognorm);
  printf("# Maximum sopt lognorm: %3.3f\n", max_sopt_lognorm);
  printf("# Average sopt alpha value: %3.3f\n", ave_sopt_alpha / nb_input_polys);

  cado_poly_clear (cpoly);
  fclose (polys_file);
  param_list_clear (pl);
  return 0;
}
