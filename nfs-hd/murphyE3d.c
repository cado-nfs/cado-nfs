#include "cado.h"
#include "utils.h"
#include <math.h>
#include "alpha3d.h"
#include "murphyE3d.h"

typedef struct {
  //Array of spq with MqLLL.
  double x;
  double y;
  double z;
} s_point3d_t;

typedef s_point3d_t point3d_t[1];
typedef s_point3d_t * point3d_ptr;
typedef const s_point3d_t * point3d_srcptr;

//Based on https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
void fibonacci_sphere(point3d_t * point, unsigned int N)
{
  double offset = 2.0 / (double) N;
  double increment = M_PI * (3.0 - sqrt(5.0));
  double r = 0.0;
  double phi = 0.0;
  for (unsigned int i = 0; i < N; i++) {
    point[i]->y = (((double)i * offset) - 1) + (offset / 2.0);
    r = sqrt(1 - point[i]->y * point[i]->y);
    phi = (double)((i + 1) % N) * increment;
    point[i]->x = cos(phi) * r;
    point[i]->z = sin(phi) * r;
  }
}

double murphyE3d(cado_poly_srcptr f, double * lpb, double volume,
    unsigned int N, int q_side, double Q, double s,
    unsigned long p, gmp_randstate_t rstate, unsigned int N_alpha)
{
  double * alpha = (double *) malloc(sizeof(double) * f->nb_polys);
  for (int i = 0; i < f->nb_polys; i++) {
    alpha[i] = alpha3d(f->pols[i], p, rstate, N_alpha);
  }

  double E = 0.0;
  double sx = pow(volume, 1.0 / 3.0) / sqrt(s);
  double sy = pow(volume, 1.0 / 3.0);
  double sz = pow(volume, 1.0 / 3.0) * sqrt(s);

  point3d_t * points = (point3d_t *) malloc(sizeof(point3d_t) * N);
  fibonacci_sphere(points, N);

  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  mpz_poly a;
  mpz_poly_init(a, 2);
  mpz_t resultant;
  mpz_init(resultant);
  double * resultants = (double *) malloc(sizeof(double) * f->nb_polys);
  double * u = (double *) malloc(sizeof(double) * f->nb_polys);
  double prob = 0.0;

  for (unsigned int i = 0; i < N; i++) {
    x = sx * points[i]->x;
    y = sy * points[i]->y;
    z = sz * points[i]->z;

    //a($0) = x * $0^2 + y * $0 + z
    mpz_poly_setcoeff_double(a, 2, x);
    mpz_poly_setcoeff_double(a, 1, y);
    mpz_poly_setcoeff_double(a, 0, z);

    for (int j = 0; j < f->nb_polys; j++) {
      //TODO: use double_poly_resultant
      mpz_poly_resultant(resultant, f->pols[j], a);
      resultants[j] = mpz_get_d(resultant);
    }
    if (q_side > -1) {
      resultants[q_side] = resultants[q_side] / Q;
    }
    for (int j = 0; j < f->nb_polys; j++) {
      u[j] = (log(fabs(resultants[j])) + alpha[j]) / log(lpb[j]);
    }
    prob = 1.0;
    for (int j = 0; j < f->nb_polys; j++) {
      prob *= dickman_rho(u[j]);
    }
    E += prob;
  }

  free(u);
  free(resultants);
  mpz_clear(resultant);
  mpz_poly_clear(a);
  free(points);
  free(alpha);

  return E / (double)N;
}

#ifdef MAIN_MURPHYE3D
void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "poly", "path to the polynomial file");
  param_list_decl_usage(pl, "lpb", "large prime bounds");
  param_list_decl_usage(pl, "V", "size of the volume of the sphere");
  param_list_decl_usage(pl, "N", "number of points on the sphere");
  param_list_decl_usage(pl, "q", "size of special-q");
  param_list_decl_usage(pl, "q_side", "side of the special-q");
  param_list_decl_usage(pl, "s", "skewness");
  param_list_decl_usage(pl, "p", "bound on p for alpha");
  param_list_decl_usage(pl, "seed", "seed for Monte Carlo for alpha");
  param_list_decl_usage(pl, "N_alpha", "bound for Monte Carlo for alpha");
}

void initialise_parameters(int argc, char * argv[], cado_poly_ptr f,
    double ** lpb, double * volume, unsigned int * N, double * Q,
    int * q_side, double * s, unsigned long * p, gmp_randstate_t rstate, unsigned * N_alpha)
{
  param_list pl;
  param_list_init(pl);
  declare_usage(pl);
  FILE * fpl;
  char * argv0 = argv[0];

  argv++, argc--;
  for( ; argc ; ) {
    if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }

    /* Could also be a file */
    if ((fpl = fopen(argv[0], "r")) != NULL) {
      param_list_read_stream(pl, fpl, 0);
      fclose(fpl);
      argv++,argc--;
      continue;
    }

    fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
    param_list_print_usage(pl, argv0, stderr);
    exit (EXIT_FAILURE);
  }

  unsigned int size_path = 1024;
  char path [size_path];
  param_list_parse_string(pl, "poly", path, size_path);
  cado_poly_read(f, path);
  ASSERT(f->nb_polys >= 2);
  unsigned int V = (unsigned int) f->nb_polys;

  * lpb = realloc(*lpb, sizeof(double) * V);
  unsigned int * lpb_size = (unsigned int *) malloc(sizeof(unsigned int) * V);
  param_list_parse_uint_list(pl, "lpb", lpb_size, (size_t) V, ",");
  for (unsigned int i = 0; i < V; i++) {
    (*lpb)[i] = pow(2.0, (double)lpb_size[i]);
  }
  free(lpb_size);

  double volume_size = 40.0;
  param_list_parse_double(pl, "V", &volume_size);
  * volume = pow(2.0, (double) volume_size);

  param_list_parse_uint(pl, "N", N);

  double qsize = 0;
  param_list_parse_double(pl, "q", &qsize);
  * Q = pow(2.0, (double)qsize);

  * q_side = -1;
  param_list_parse_int(pl, "q_side", q_side);
  ASSERT(* q_side < (int)V);

  * s = 1.0;
  param_list_parse_double(pl, "s", s);

  * p = 10000;
  param_list_parse_ulong(pl, "p", p);

  unsigned long seed;
  if (param_list_parse_ulong(pl, "seed", &seed))
      gmp_randseed_ui(rstate, seed);

  * N_alpha = 10000;
  param_list_parse_uint(pl, "N_alpha", N_alpha);

  param_list_clear(pl);
}

int main(int argc, char ** argv)
{
  cado_poly f;
  double * lpb;
  double volume;
  unsigned int N;
  double Q;
  int q_side;
  double s;
  unsigned long p;
  gmp_randstate_t rstate;
  unsigned int N_alpha;
  cado_poly_init(f);

  gmp_randinit_default(rstate);
  lpb = NULL;

  initialise_parameters(argc, argv, f, &lpb, &volume, &N, &Q, &q_side, &s, &p,
      rstate, &N_alpha);

  printf("%le\n", murphyE3d(f, lpb, volume, N, q_side, Q, s, p, rstate, N_alpha));

  free(lpb);
  cado_poly_clear(f);
  return 0;
}
#endif // MAIN_MURPHYE3D
