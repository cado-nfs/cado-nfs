#include "cado.h" // IWYU pragma: keep
// TODO: did not work.

#include <stdlib.h>
#include "parameters.h"

int main()
{
  /*mpz_t p;*/
  /*mpz_init(p);*/
  /*mpz_set_str(p, "10822639589", 10);*/
  /*unsigned int n = 6;*/

  /*uint64_t number_a = 100;*/

  /*unsigned int lpb_min = 25;*/
  /*unsigned int lpb_max = 26;*/

  /*unsigned int t_min = 3;*/
  /*unsigned int t_max = 4;*/

  /*uint64_t size_start = 34359738368;*/

  /*mpz_poly h;*/
  /*mpz_poly_init(h, n);*/
  /*for (int k = 0; k <= (int)n; k++) {*/
    /*mpz_poly_setcoeff_si(h, k, 1);*/
  /*}*/

  /*int coeff0 = -2;*/
  /*int coeff1 = 2;*/

  /*unsigned int nb_times = 2;*/
  
  /*double weight_0 = 0.5;*/
  /*double weight_1 = 0.5;*/

  /*mpz_poly f0;*/
  /*mpz_poly_init(f0, n);*/
  /*mpz_poly f1;*/
  /*mpz_poly_init(f1, n);*/

  /*find_parameters_adapted(p, n, number_a, lpb_min, lpb_max, t_min, t_max,*/
      /*size_start, h, coeff0, coeff1, nb_times, weight_0, weight_1);*/

  /*mpz_poly_clear(h);*/
  /*mpz_clear(p);*/

  sieving_bound_t H;
  unsigned int t = 3;
  sieving_bound_init(H, t);
  unsigned int nb_elem = 60;
  mpz_t p;
  mpz_init(p);
  mpz_set_str(p, "10822639589", 10);
  unsigned int n = 6;
  sieving_region_square(H, nb_elem);
  sieving_bound_fprintf(stdout, H);
  printf("Size: %f\n", log2(size_sieving_region(H)));
  printf("%d\n", sieving_region_classical(H, p, n, nb_elem));
  sieving_bound_fprintf(stdout, H);
  printf("Size: %f\n", log2(size_sieving_region(H)));
  sieving_bound_clear(H);

  t = 4;
  sieving_bound_init(H, t);
  printf("%d\n", sieving_region_classical(H, p, n, nb_elem));
  sieving_bound_fprintf(stdout, H);
  printf("Size: %f\n", log2(size_sieving_region(H)));
  sieving_bound_clear(H);
  
  t = 5;
  sieving_bound_init(H, t);
  printf("%d\n", sieving_region_classical(H, p, n, nb_elem));
  sieving_bound_fprintf(stdout, H);
  printf("Size: %f\n", log2(size_sieving_region(H)));
  sieving_bound_clear(H);

  mpz_clear(p);
}
