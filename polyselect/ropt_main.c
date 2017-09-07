/**
   Root optimization for polynomials in the number field sieve.

   ./polyselect_ropt [options]

    -inputpolys    "input polynomial file"
    -ropteffort    "effort scales"
    -t             "number of threads"
    -area          "sieving area"
    -I             "I-value (alternative to -area)"
    -Bf            "algebraic smoothness bound"
    -Bg            "rational smoothness bound"
    -v             "toggle verbose"
    -boundmaxlognorm  "maximum lognorm to bound the rotation"

   Only the -inputpolys option (follwed by the FILENAME) is compulsory.
   Input polynomials are in "CADO". They are read from FILENAME.
   Output polynomials are in "CADO" format.

   Please report bugs to the public mailing list
       cado-nfs-discuss@lists.gforge.inria.fr
*/


#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif
#include "utils.h"
#include "rho.h"
#include "auxiliary.h"
#include "murphyE.h"
#include <pthread.h>
#include "portability.h"
#include "area.h"
#include "ropt.h"

pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER; /* used as mutual exclusion
                                                     lock for output */
unsigned int nthreads = 1;
int tot_found = 0; /* total number of polynomials */
cado_poly best_poly;
double best_MurphyE = 0.0; /* Murphy's E (the larger the better) */
ropt_param_t ropt_param; /* params for ropt algorithms */
double total_exp_E = 0.0; /* cumulated expected E for input polynomials */
double total_E = 0.0; /* cumulated E-value for input polynomials */
double nb_read = 0.0; /* number of read polynomials so far */
double nb_optimized = 0.0; /* number of optimized polynomials so far */

/**
 * Usage
 */
void
usage_adv (char **argv)
{
  fprintf (stderr, "Error: Unhandled parameter: %s\n", argv[1]);
  fprintf (stderr, "\n");
  fprintf (stderr, "Options followed by %s\n", argv[1]);
  fprintf (stderr, "Usage: %s %s -f fname [options]\n", argv[0], argv[1]);
  fprintf (stderr, "       %s %s -f fname --s2 [options]\n", argv[0], argv[1]);
  fprintf (stderr, "       %s %s -fm fname -n N -d D [options]\n", argv[0], argv[1]);
  fprintf (stderr, "       %s %s -fm fname --s2 -n N -d D [options]\n", argv[0], argv[1]);
  fprintf (stderr, "       %s %s [options]\n", argv[0], argv[1]);
  fprintf (stderr, "       %s %s --s2 [options]\n", argv[0], argv[1]);
  fprintf (stderr, "\n");
  fprintf (stderr, "Options:\n");
  fprintf (stderr, " -f fname      read polynomials in CADO format.\n");
  fprintf (stderr, " -fm fname     read polynomials in MSIEVE format (use together\n");
  fprintf (stderr, "                with options -n and -d).\n");
  fprintf (stderr, " -n N          (only in -fm mode) the integer to be factored.\n");
  fprintf (stderr, " -d D          (only in -fm mode) the degree of the polynomial.\n");
  fprintf (stderr, " -amin L       L is the lower bound for quadratic rotation.\n");
  fprintf (stderr, " -amax R       R is the upper bound for quadratic rotation.\n");
  fprintf (stderr, " -bmax U       U and -U are the sieving length for linear rotation.\n");
  fprintf (stderr, " -cmax V       V and -V are the sieving length for constant rotation.\n");
  fprintf (stderr, " -e N E_1 E_2 ... E_N\n");
  fprintf (stderr, "               N is the number of small primes in finding sublattices.\n");
  fprintf (stderr, "               E_1 to E_N are the prime powers in finding sublattices.\n");
  fprintf (stderr, " -norm M       M is the (estimated) lognorm upper bound in ropt.\n");
  fprintf (stderr, " -v            number of occurrences defines verbose level {0, 1, 2, 3} (default 0).\n");
  fprintf (stderr, " --s2          (switch) sieve-only mode (use together with the following\n");
  fprintf (stderr, "               four options: -a, -b, -c and -mod.\n");
  fprintf (stderr, " -a A          Fix the quadratic rotation of the sublattice by A.\n");
  fprintf (stderr, " -b B          Fix the linear rotation of the sublattice by B.\n");
  fprintf (stderr, " -c C          Fix the constant rotation of the sublattice by C.\n");
  fprintf (stderr, " -mod M        M is the sublattice modulus.\n");
  fprintf (stderr, " --skip_ropt   (switch) skip root sieve (use with option -fm).\n");
  fprintf (stderr, " --sopt        (switch) do size optimization for polynomial.\n");
  fprintf (stderr, " --gen_raw     (switch) regenerate raw polynomial and skip root sieve.\n");
  fprintf (stderr, " -Bf F         algebraic smoothness bound (default %.2e).\n", BOUND_F);
  fprintf (stderr, " -Bg G         rational smoothness bound (default %.2e).\n", BOUND_G);
  fprintf (stderr, " -area A       sieving area (default %.2e).\n", AREA);
  fprintf (stderr, " -I I          I-value (alternative to -area).\n");
  fprintf (stderr, " -ropteffort M   sieving effort ranging from 1 to 10 (default %.0f).\n", DEFAULT_ROPTEFFORT);

  fprintf (stderr, "\nExample 1: %s %s -f fname\n", argv[0], argv[1]);
  fprintf (stderr, "Root optimization for all CADO-formatted polynomials in 'fname'.\n");

  fprintf (stderr, "\nExample 2: %s %s -f fname -amin -512 -amax 512  -norm 70\n", argv[0], argv[1]);
  fprintf (stderr, "As above, but restricts the quadratic rotation between -512\n"
           "and 512 and the sieving region by norm 70.\n");

  fprintf (stderr, "\nExample 3: %s %s -f fname -amin -512 -amax 512 -e 5 7 4 3 2 2 -bmax 16 -cmax 10000000\n", argv[0], argv[1]);
  fprintf (stderr, "As above, but uses five prime factors (2, 3, 5, 7, 11) in the\n"
           "sublattice with powers (7, 4, 3, 2, 2). It also tells ropt to\n"
           "root sieve a region of 16 by 10000000.\n");

  fprintf (stderr, "\nExample 4: %s %s -fm rsa768.poly -amin -512 -amax 512 -e 5 7 4 3 2 2 -bmax 16 -cmax 10000000 -n $N -d 6\n", argv[0], argv[1]);
  fprintf (stderr, "As above, but reads msieve format where each line contains\n"
           "'c_d Y1 Y0'. The parameters -n and -d are compulsory.\n");

  fprintf (stderr, "\nExample 5: %s %s -f fname --s2 -a 12 -b 345 -c 6789 -mod 1814400 -bmax 16 -cmax 10000000\n", argv[0], argv[1]);
  fprintf (stderr, "Sieve-only mode. Assume that we know the polynomial has good\n"
           "root property at rotation (12*x^2 + 345*x + 6789), we want to\n"
           "search 12*x^2 + (345 + 1814400*i)*x + (6789 + 1814400*j) where\n"
           "i, j are bounded by -bmax and -cmax.\n");
  exit(1);
}


/**
 * parse manually input parameters to param.
 */
static void
ropt_parse_param ( int argc,
                   char **argv,
                   ropt_param_t param )
{
  int I = 0;
  int I_area = 0; /* 0 when neither -area nor I is given,
                     1 when -area is given
                     2 when -I is given */

  /* filename only */
  if (argc > 1) {
    /* stage 2 (root sieve) parameters only */
    if (strcmp (argv[1], "--s2") == 0) {

      param->stage_flag = 2;
      argv += 1;
      argc -= 1;
      while (argc >= 2 && argv[1][0] == '-') {
        if (strcmp (argv[1], "-v") == 0)
        {
          param->verbose ++;
          argv += 1;
          argc -= 1;
        }
        else if (argc >= 3 && strcmp (argv[1], "-a") == 0)
        {
          param->s2_w = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-b") == 0)
        {
          mpz_set_str (param->s2_u, argv[2], 10);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-c") == 0)
        {
          mpz_set_str (param->s2_v, argv[2], 10);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-mod") == 0)
        {
          mpz_set_str (param->s2_mod, argv[2], 10);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-bmax") == 0)
        {
          param->s2_Amax = atol (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-cmax") == 0)
        {
          param->s2_Bmax = atol (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-n") == 0)
        {
          mpz_set_str (param->n, argv[2], 10);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-norm") == 0)
        {
          param->bound_lognorm = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-d") == 0)
        {
          param->d = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-Bf") == 0)
        {
          bound_f = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-Bg") == 0)
        {
          bound_g = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-area") == 0)
        {
          if (I_area == 2)
            {
              fprintf (stderr, "Error, both -I and -area are given\n");
              exit (1);
            }
          area = atof (argv[2]);
          I_area = 1;
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-I") == 0)
        {
          if (I_area == 1)
            {
              fprintf (stderr, "Error, both -area and -I are given\n");
              exit (1);
            }
          I_area = 2;
          I = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-ropteffort") == 0)
        {
          param->effort = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "--skip_ropt") == 0)
        {
          fprintf (stderr, "Error: cannot use --skip_ropt with --s2.\n");
          exit(1);
        }
        else if (argc >= 3 && strcmp (argv[1], "--gen_raw") == 0)
        {
          fprintf (stderr, "Error: cannot use --gen_raw with --s2.\n");
          exit(1);
        }
        else if (argc >= 3 && strcmp (argv[1], "--sopt") == 0)
        {
          fprintf (stderr, "Error: cannot use --sopt with --s2.\n");
          exit(1);
        }
        else {
          usage_adv (argv);
        }
        if (I_area == 2)
          area = bound_f * pow (2.0, (double) (2*I-1));
      }
    }
    /* stage 1 and stage 2 parameters */
    else {
      while (argc >= 2 && argv[1][0] == '-') {
        if (strcmp (argv[1], "-v") == 0)
        {
          param->verbose ++;
          argv += 1;
          argc -= 1;
        }
        else if (argc >= 3 && strcmp (argv[1], "-amin") == 0)
        {
          param->w_left_bound = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-amax") == 0)
        {
          param->w_length = atoi (argv[2]) - param->w_left_bound + 1;
          if (param->w_length < 0) {
            fprintf (stderr, "Error in options -amin and/or -amax.\n");
            exit(1);
          }
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-bmax") == 0)
        {
          param->s2_Amax = atol (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-cmax") == 0)
        {
          param->s2_Bmax = atol (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-n") == 0)
        {
          mpz_set_str (param->n, argv[2], 10);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-norm") == 0)
        {
          param->bound_lognorm = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-d") == 0)
        {
          param->d = atoi (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-Bf") == 0)
        {
          bound_f = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-Bg") == 0)
        {
          bound_g = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-area") == 0)
        {
          area = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-ropteffort") == 0)
        {
          param->effort = atof (argv[2]);
          argv += 2;
          argc -= 2;
        }
        else if (argc >= 2 && strcmp (argv[1], "--skip_ropt") == 0)
        {
          param->skip_ropt = 1;
          argv += 1;
          argc -= 1;
        }
        else if (argc >= 2 && strcmp (argv[1], "--gen_raw") == 0)
        {
          param->gen_raw = 1;
          argv += 1;
          argc -= 1;
        }
        else if (argc >= 2 && strcmp (argv[1], "--sopt") == 0)
        {
          param->sopt = 1;
          argv += 1;
          argc -= 1;
        }
        else if (argc >= 3 && strcmp (argv[1], "-e") == 0)
        {
          param->stage_flag = 1;
          param->s1_num_e_sl = atoi (argv[2]);
          argv += 2;
          argc -= 2;
          int j = 0;
          for (j = 0; j < param->s1_num_e_sl; j ++) {
            param->s1_e_sl[j] = (unsigned short) atoi (argv[1]);
            argv += 1;
            argc -= 1;
          }
        }
        else {
          usage_adv (argv);
        }
      }
    }
  }
}


static void
ropt_wrapper (cado_poly_ptr input_poly, unsigned int poly_id,
              ropt_time_t tott)
{
  double curr_MurphyE, st, st1;
  mpz_t t;
  cado_poly ropt_poly;
  ropt_time_t eacht;

  cado_poly_init (ropt_poly);
  eacht->ropt_time = 0.0;
  eacht->ropt_time_stage1 = 0.0;
  eacht->ropt_time_tuning = 0.0;
  eacht->ropt_time_stage2 = 0.0;

  if (nthreads > 1)
    pthread_mutex_lock (&lock);
  printf ("\n### input polynomial %u ###\n", poly_id);
  total_exp_E += cado_poly_fprintf_with_info (stdout, input_poly, "# ", 0);
  nb_read += 1.0;
  fflush (stdout);
  if (nthreads > 1)
    pthread_mutex_unlock (&lock);

  /* If the content of the algebraic polynomial has content <> 1, then print a
     warning (this should not be frequent) and divide all coefficients of the
     polynomial by the content. */
  mpz_init (t);
  mpz_poly_content (t, input_poly->pols[ALG_SIDE]);
  if (mpz_cmp_ui (t, 1) != 0 && !mpz_divisible_p (input_poly->n, t))
  {
    gmp_printf ("# WARNING: the content of the algebraic side of polynomial %u "
                "is not 1 (%Zd). The input polynomial will be divided by its "
                "content.\n", poly_id, t);
    mpz_poly_divexact_mpz (input_poly->pols[ALG_SIDE], input_poly->pols[ALG_SIDE], t);
  }
  mpz_clear (t);

  st = seconds_thread ();
  ropt_polyselect (ropt_poly, input_poly, ropt_param, eacht);
  st1 = seconds_thread ();

  /* MurphyE */
  /* use the skewness minimizing the sum of the lognorms */
  ropt_poly->skew = L2_combined_skewness2 (ropt_poly->pols[0],
                                           ropt_poly->pols[1],
                                           SKEWNESS_DEFAULT_PREC);
  curr_MurphyE = MurphyE (ropt_poly, bound_f, bound_g, area, MURPHY_K);

  if (nthreads > 1)
    pthread_mutex_lock (&lock);

  /* update time */
  tott->ropt_time += st1 - st;
  tott->ropt_time_stage1 += eacht->ropt_time_stage1;
  tott->ropt_time_tuning += eacht->ropt_time_tuning;
  tott->ropt_time_stage2 += eacht->ropt_time_stage2;

  if (curr_MurphyE > best_MurphyE)
    {
      best_MurphyE = curr_MurphyE;
      cado_poly_set (best_poly, ropt_poly);
    }
  printf ("\n### root-optimized polynomial %u ###\n", poly_id);
  total_E += cado_poly_fprintf_with_info_and_MurphyE (stdout, ropt_poly,
                                                      curr_MurphyE,
                                                      bound_f, bound_g, area,
                                                      "# ");
  nb_optimized += 1.0;
  printf ("### Best MurphyE so far is %.2e, av. exp_E %.2f, av. E %.2f\n",
          best_MurphyE, total_exp_E / nb_read, total_E / nb_optimized);
  fflush (stdout);

  if (nthreads > 1)
    pthread_mutex_unlock (&lock);

  cado_poly_clear (ropt_poly);
}

/**
 * Interface main_adv(). This will call ropt_on_cadopoly().
 */
static int
main_adv (int argc, char **argv)
{

  char **argv0 = argv;
  argv += 1;
  argc -= 1;
  int i;

  ropt_param_t param;
  ropt_param_init (param);

  /* print command-line arguments */
  fprintf (stderr, "# %s.r%s", *argv, cado_revision_string);
  for (i = 1; i < argc; i++)
    fprintf (stderr, " %s", *(argv+i));
  fprintf (stderr, "\n");

  /* read polynomials in cado format */
  if (argc >= 3 && strcmp(argv[1], "-f") == 0) {

    FILE *file = NULL;
    char *filename = NULL;
    filename = argv[2];
    argv += 2;
    argc -= 2;

    file = fopen(filename, "r");
    if (file == NULL) {
      fprintf(stderr, "# Error in reading file\n");
      exit (1);
    }

    /* parse parameters */
    ropt_parse_param (argc, argv, param);

    /* call ropt_on_cadopoly() */
    ropt_on_cadopoly (file, param);

    fclose (file);
    ropt_param_free (param);

    return EXIT_SUCCESS;
  }
  /* read polynomials in msieve format */
  else if (argc >= 3 && strcmp(argv[1], "-fm") == 0) {

    FILE *file = NULL;
    char *filename = NULL;
    filename = argv[2];
    argv += 2;
    argc -= 2;

    file = fopen(filename, "r");
    if (file == NULL) {
      fprintf(stderr, "# Error in reading file\n");
      exit (1);
    }

    /* parse parameters */
    ropt_parse_param (argc, argv, param);

    if (mpz_cmp_ui(param->n, 0) == 0) {
      fprintf(stderr, "# Error: please input parameter \"-n number\"\n");
      exit (1);
    }
    if ( param->d == 0 ) {
      fprintf(stderr, "# Error: please input parameter \"-n degree\"\n");
      exit (1);
    }

    /* call ropt_on_msievepoly() */
    ropt_on_msievepoly (file, param);

    fclose (file);
    ropt_param_free (param);

    return EXIT_SUCCESS;
  }
  /* read polynomials from stdin */
  else if (argc >= 3)
  {

    /* parse parameters */
    ropt_parse_param (argc, argv, param);

    /* call ropt_on_stdin() */
    ropt_on_stdin (param);

    ropt_param_free (param);

    return EXIT_SUCCESS;
  }
  else
    usage_adv (argv0);

  return EXIT_SUCCESS;
}


static void
declare_usage_basic (param_list pl)
{
  char str[200];
  param_list_decl_usage(pl, "inputpolys", "root-sieve the size-optimized "
                                          "polynomials given in this file");
  snprintf (str, 200, "root-sieve effort (default %.0f)", DEFAULT_ROPTEFFORT);
  param_list_decl_usage(pl, "ropteffort", str);
  param_list_decl_usage(pl, "t", "number of threads to use (default 1)");
  snprintf (str, 200, "sieving area (default %.2e)", AREA);
  param_list_decl_usage(pl, "area", str);
  snprintf (str, 200, "I-value (alternative to -area)");
  param_list_decl_usage(pl, "I", str);
  snprintf (str, 200, "algebraic smoothness bound (default %.2e)", BOUND_F);
  param_list_decl_usage(pl, "Bf", str);
  snprintf (str, 200, "rational smoothness bound (default %.2e)", BOUND_G);
  param_list_decl_usage(pl, "Bg", str);
  param_list_decl_usage(pl, "v", "(switch) verbose mode");
  param_list_decl_usage(pl, "boundmaxlognorm", "Maximum lognorm. Used to compute"
                                               " bounds for rotations for the "
                                               "sieve stage of ropt");
  verbose_decl_usage(pl);
}


static void
usage_basic (const char *argv, const char * missing, param_list pl)
{
  if (missing) {
    fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
            missing);
  }
  param_list_print_usage (pl, argv, stderr);
  param_list_clear (pl);
  exit (EXIT_FAILURE);
}


/**
 * Interface main_basic()
 */
static int
main_basic (int argc, char **argv)
{
  char **argv0 = argv;
  const char *polys_filename = NULL;
  double st0 = seconds ();
  FILE *polys_file = NULL;
  cado_poly *input_polys = NULL;
  unsigned int nb_input_polys = 0; /* number of input polynomials */
  unsigned int size_input_polys = 16; /* Size of input_polys tab. */
  ropt_time_t tott;
  int I;

  tott->ropt_time = 0.0;
  tott->ropt_time_stage1 = 0.0;
  tott->ropt_time_tuning = 0.0;
  tott->ropt_time_stage2 = 0.0;

  cado_poly_init (best_poly);
  input_polys = (cado_poly *) malloc (size_input_polys * sizeof (cado_poly));
  ASSERT_ALWAYS (input_polys != NULL);
  for (unsigned int i = 0; i < size_input_polys; i++)
    cado_poly_init (input_polys[i]);
  ropt_param_init (ropt_param);
  ropt_param->effort = DEFAULT_ROPTEFFORT; /* ropt effort */

  /* read params */
  param_list pl;
  param_list_init (pl);
  declare_usage_basic(pl);
  param_list_configure_switch (pl, "-v", &(ropt_param->verbose));

  if (argc == 1)
    usage_basic (argv0[0], NULL, pl);

  argv++, argc--;
  for ( ; argc; ) {
    if (param_list_update_cmdline (pl, &argc, &argv)) continue;
    fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
    usage_basic (argv0[0], NULL, pl);
  }

  param_list_parse_uint (pl, "t", &nthreads);

  if (param_list_parse_double (pl, "Bf", &bound_f) == 0) /* no -Bf */
    bound_f = BOUND_F;
  if (param_list_parse_double (pl, "Bg", &bound_g) == 0) /* no -Bg */
    bound_g = BOUND_G;
  int has_area = param_list_parse_double (pl, "area", &area);
  int has_I = param_list_parse_int (pl, "I", &I);
  if (has_area == 0 && has_I == 0) /* no -area nor -I */
    area = AREA;
  else if (has_area && has_I)
    {
      fprintf (stderr, "Error, both -area and -I given\n");
      exit (1);
    }
  else if (has_I)
    area = bound_f * pow (2.0, (double) (2*I-1));

  /* sieving effort that passed to ropt */
  param_list_parse_double (pl, "ropteffort", &(ropt_param->effort));
  /* param for ropt */
  param_list_parse_double (pl, "boundmaxlognorm", &(ropt_param->bound_lognorm));
  /* filename of the file with the list of polynomials to root-sieve */
  polys_filename = param_list_lookup_string (pl, "inputpolys");

  if (param_list_warn_unused(pl))
    usage_basic (argv0[0], NULL, pl);

  /* print command line */
  verbose_interpret_parameters(pl);
  param_list_print_command_line (stdout, pl);

  /* Check that nthreads is >= 1 */
  if (nthreads == 0)
  {
    fprintf (stderr, "Error, -t must be non-zero.\n");
    usage_basic (argv0[0], NULL, pl);
  }

  if (polys_filename == NULL)
    usage_basic (argv0[0], "inputpolys", pl);

  printf ("# Info: Will use %u thread%s\n# Info: ropteffort = %.0f\n", nthreads,
          (nthreads > 1) ? "s": "", ropt_param->effort);

  printf ("# Info: L1_cachesize = %u, size_tune_sievearray = %u\n",
          L1_cachesize, size_tune_sievearray);

  /* Open file containing polynomials. */
  printf ("# Reading polynomials from %s\n", polys_filename);
  polys_file = fopen(polys_filename, "r");
  if (polys_file == NULL)
  {
    perror("Could not open file");
    exit(EXIT_FAILURE);
  }

  /* Remove initial empty lines */
  int c;
  while ((c = fgetc (polys_file)) == '\n');
  if (c != EOF)
    ungetc (c, polys_file);

  /* Read all polynomials from file. Store then in input_polys. */
  while (cado_poly_read_next_poly_from_stream (input_polys[nb_input_polys],
                                               polys_file))
  {
    nb_input_polys++;
    if (nb_input_polys == size_input_polys) /* Realloc if needed */
    {
      if (ropt_param->verbose > 0)
        fprintf (stderr, "# Reallocating input_polys\n");
      unsigned int new_size = 2 * size_input_polys;
      input_polys = (cado_poly *) realloc (input_polys, new_size * sizeof (cado_poly));
      ASSERT_ALWAYS (input_polys != NULL);
      for (unsigned int i = size_input_polys; i < new_size; i++)
        cado_poly_init (input_polys[i]);
      size_input_polys = new_size;
    }
  }
  /* Did we stop at the end of the file or was there an error. */
  if (!feof (polys_file))
  {
    fprintf (stderr, "Error while reading file %s.\n", polys_filename);
    abort ();
  }
  fclose (polys_file);
  printf ("# %u polynomial(s) read.\n", nb_input_polys);

  /* Main loop: do root-optimization on input_polys. */
#ifdef HAVE_OPENMP
  omp_set_num_threads (nthreads);
#pragma omp parallel
#pragma omp master
  printf ("# Info: Using OpenMP with %u thread(s)\n", omp_get_num_threads ());
#pragma omp parallel for schedule(dynamic)
#endif
  for (unsigned int i = 0; i < nb_input_polys; i++)
    ropt_wrapper (input_polys[i], i, tott);

  /* print total time and rootsieve time.
     These two lines gets parsed by the script. */
  printf ("\n\n# Stat: total phase took %.2fs\n", seconds () - st0);

#ifndef HAVE_RUSAGE_THREAD /* rootsieve_time is correct only if RUSAGE_THREAD
                              works or in mono-thread mode */
  if (nthreads == 1)
#endif
  {
    printf ("# Stat: rootsieve took %.2fs\n", tott->ropt_time);
    printf ("# Stat:  (stage 1 took %.2fs)\n", tott->ropt_time_stage1);
    printf ("# Stat:  (tuning took %.2fs)\n", tott->ropt_time_tuning);
    printf ("# Stat:  (stage 2 (sieving) took %.2fs)\n", tott->ropt_time_stage2);
  }

  if (best_MurphyE == 0.0)
  {
    if (nb_input_polys > 0)
    {
      /* At least one poly was parsed and the best MurphyE is still 0 => Error */
      fprintf (stderr, "Error, the best MurphyE value is still 0 at the end.\n");
      abort ();
    }
    else
    {
      /* No polynomials were parsed in the input files but the whole file was
         parsed without errors (see !feof above) => Warning */
      printf ("# WARNING: No polynomials were found in the input file %s\n",
              polys_filename);
    }
  }
  else
  {
    printf ("# Best polynomial found:\n");
    cado_poly_fprintf_with_info_and_MurphyE (stdout, best_poly, best_MurphyE,
                                             bound_f, bound_g, area, NULL);
    printf ("# Average exp_E: %.2f, average E: %.2f\n",
            total_exp_E / (double) nb_input_polys,
            total_E / (double) nb_input_polys);
  }

  ropt_param_free (ropt_param);
  cado_poly_clear (best_poly);
  for (unsigned int i = 0; i < size_input_polys; i++)
    cado_poly_clear (input_polys[i]);
  free (input_polys);
  param_list_clear (pl);

  return EXIT_SUCCESS;
}


/**
 * There are two interfaces: either call "main_adv" or "main_basic"
 * depending whether we have the option --adv.
 * Note the cadoprograms.py uses the main_basic() interface.
 */
int
main (int argc, char **argv)
{
  /* usage */

  /* detect L1 cache size */
  ropt_L1_cachesize ();

  if (argc > 1 && strcmp(argv[1], "--adv") == 0)
    return main_adv (argc, argv);
  else
    return main_basic (argc, argv);
}
