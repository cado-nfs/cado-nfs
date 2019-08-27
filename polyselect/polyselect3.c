/* Recomputes the Murphy-E value with a larger value of ALPHA_BOUND
   on the top polynomials found after the rootsieve:

   polyselect3 -poly cxxx.poly
*/

#include "cado.h"
#include "utils.h"
#include "murphyE.h"
#include "auxiliary.h"

static void
declare_usage (param_list pl)
{
  param_list_decl_usage(pl, "poly", "polynomial file");
  param_list_decl_usage(pl, "Bf", "factor base bound on the algebraic side");
  param_list_decl_usage(pl, "Bg", "factor base bound on the linear side");
  param_list_decl_usage(pl, "area", "sieving area");
  param_list_decl_usage(pl, "I", "sieving parameter I");
  verbose_decl_usage(pl);
}

int
main (int argc, char *argv[])
{
  param_list pl;
  cado_poly cpoly;
  FILE *f;
  char *argv0 = argv[0];
  double Bf, Bg, area;
  int I;

  param_list_init(pl);
  declare_usage(pl);
  cado_poly_init(cpoly);

  argv++, argc--;
  for( ; argc ; ) {
      if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }

      /* Could also be a file */
      if ((f = fopen(argv[0], "r")) != NULL) {
          param_list_read_stream(pl, f, 0);
          fclose(f);
          argv++,argc--;
          continue;
      }

      fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
      param_list_print_usage(pl, argv0, stderr);
      exit (EXIT_FAILURE);
  }
  verbose_interpret_parameters (pl);
  param_list_print_command_line (stdout, pl);

  const char * filename;
  if ((filename = param_list_lookup_string (pl, "poly")) == NULL)
    {
      fprintf (stderr, "Error: parameter -poly is mandatory\n");
      param_list_print_usage (pl, argv0, stderr);
      exit (EXIT_FAILURE);
    }

  if (param_list_parse_double (pl, "Bf", &Bf) == 0)
    {
      fprintf (stderr, "Error: parameter -Bf is mandatory\n");
      param_list_print_usage (pl, argv0, stderr);
      exit (EXIT_FAILURE);
    }

  if (param_list_parse_double (pl, "Bg", &Bg) == 0)
    {
      fprintf (stderr, "Error: parameter -Bg is mandatory\n");
      param_list_print_usage (pl, argv0, stderr);
      exit (EXIT_FAILURE);
    }

  int has_area = param_list_parse_double (pl, "area", &area);
  int has_I = param_list_parse_int (pl, "I", &I);
  if (has_area == 0 && has_I == 0)
    {
      fprintf (stderr, "Error: parameter -area or -I is mandatory\n");
      param_list_print_usage (pl, argv0, stderr);
      exit (EXIT_FAILURE);
    }

  if (!cado_poly_read (cpoly, filename))
    {
      fprintf (stderr, "Error reading polynomial file %s\n", filename);
      exit (EXIT_FAILURE);
    }

  if (has_I)
    area = Bf * pow (2.0, (double) (2 * I - 1));

  double e = MurphyE (cpoly, Bf, Bg, area, MURPHY_K, 10 * ALPHA_BOUND);
  f = fopen (filename, "a");
  fprintf (f, "# MurphyE (Bf=%.3e,Bg=%.3e,area=%.3e) = %.2e [B=%u]\n",
           Bf, Bg, area, e, 10 * ALPHA_BOUND);
  fclose (f);

  cado_poly_clear (cpoly);
  param_list_clear(pl);

  return 0;
}
