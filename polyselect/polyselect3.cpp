/* Recomputes the Murphy-E value with a larger value of ALPHA_BOUND
   on the top polynomials found after the rootsieve:

   polyselect3 -poly cxxx.poly -num 10 -Bf ... -Bg ... -area ...

   will process cxxx.poly.0, cxxx.poly.1, ..., cxxx.poly.9
   and add the new Murphy-E value at the end of each file.
*/

#include "cado.h" // IWYU pragma: keep
/* The following avoids to put #ifdef HAVE_OPENMP ... #endif around each
   OpenMP pragma. It should come after cado.h, which sets -Werror=all.
#ifdef  __GNUC__
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif
 *
 * unfortunately, while it looks like a reasonable thing to do in theory,
 * it's gcc specific. We can't expect such a thing to work with other
 * compilers.
 */
#include <cstdio>      // FILE
#include <cstdlib>     // exit ...
#include "murphyE.hpp"
#include "cado_poly.hpp"
#include "auxiliary.hpp"
#include "params.hpp"
#include "omp_proxy.h" // IWYU pragma: keep
#include "verbose.hpp"             // verbose_output_print
#include "polyselect_alpha.h"

static void
declare_usage (cxx_param_list & pl)
{
  pl.declare_usage("poly", "polynomial prefix");
  pl.declare_usage("t", "number of threads");
  pl.declare_usage("num", "number of files to process");
  pl.declare_usage("Bf", "factor base bound on the algebraic side");
  pl.declare_usage("Bg", "factor base bound on the linear side");
  pl.declare_usage("area", "sieving area");
  pl.declare_usage("v", "verbose level");
  verbose_decl_usage(pl);
}

int
main (int argc, char const * argv[])
{
  cxx_param_list pl;
  double Bf, Bg, area;
  int nthreads = 1;
  int num = 1; /* number of files to process */
  int verbose = 0;

  declare_usage(pl);
  pl.configure_switch_old("-v", &verbose);

  pl.process_command_line_and_extra_parameter_files(argc, argv);

  verbose_interpret_parameters (pl);
  pl.print_command_line(stdout);

  pl.parse("t", nthreads);
#ifdef HAVE_OPENMP
  omp_set_num_threads (nthreads);
#endif
  pl.parse("num", num);

  const char * filename = pl.lookup_old("poly");
  if (!filename)
      pl.fail("Error: parameter -poly is mandatory\n");

  pl.parse_mandatory("Bf", Bf);
  pl.parse_mandatory("Bg", Bg);
  pl.parse_mandatory("area", area);

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
  for (int i = 0; i < num; i++)
    {
      cxx_cado_poly cpoly;
      char s[1024];
      FILE *fp;

      sprintf (s, "%s.%d", filename, i);
      if (!cpoly.read(s))
        {
          fprintf (stderr, "Error reading polynomial file %s\n", s);
          exit (EXIT_FAILURE);
        }

      double e = MurphyE (cpoly, Bf, Bg, area, MURPHY_K, 10 * get_alpha_bound ());
      fp = fopen (s, "a");
      fprintf (fp, "# MurphyF (Bf=%.3e,Bg=%.3e,area=%.3e) = %.3e\n",
               Bf, Bg, area, e);
      fclose (fp);
    }

  return 0;
}
