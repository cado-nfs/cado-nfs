#include "cado.h" // IWYU pragma: keep

#define _GNU_SOURCE // for pthread_tryjoin_np()
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h> 
#include <unistd.h>
#include <pthread.h>
#include <sys/time.h>
#include <gmp.h>
#include "cado_poly.h"  // cado_poly
#include "ecm.h"
#include "lll.h"        // mat_Z LLL
#include "mpz_poly.h"   // mpz_poly
#include "smooth_detect.h"
#include "macros.h"

static double default_B1done;

/*
ECM=...
gcc -Wall -g -std=c99 -I${ECM}/include descent_init_Fp.c smooth_detect.c -o descent_init_Fp -L${ECM}/lib -lecm -lgmp -lm -lpthread
*/

pthread_cond_t cond_found = PTHREAD_COND_INITIALIZER;
pthread_mutex_t mut_found = PTHREAD_MUTEX_INITIALIZER;
static unsigned int thid_found = 0;

pthread_mutex_t mut_die = PTHREAD_MUTEX_INITIALIZER;
static unsigned int please_die = 0;

void HalfGcd(mpz_t a, mpz_t b, mpz_t u) {
  mpz_t v, w, x, q, r;
  mpz_set_ui(u, 1);
  mpz_init_set_ui(w, 0);
  mpz_init_set_ui(v, 0);
  mpz_init_set_ui(x, 1);
  mpz_init(q);
  mpz_init(r);
  /* invariant: a = u*a0 + v*b0 */
  while (mpz_cmpabs(a, u) > 0) {
    mpz_tdiv_qr(q, r, a, b);
    mpz_swap(a, b);
    mpz_swap(b, r);
    mpz_submul(u, q, w);
    mpz_swap(u, w);
    mpz_submul(v, q, x);
    mpz_swap(v, x);
  }
  mpz_clear(v);
  mpz_clear(w);
  mpz_clear(x);
  mpz_clear(q);
  mpz_clear(r);
}

typedef struct {
  mpz_t p;
  int n;          // extension degree
  mpz_t z;        // element to descend (n=1)
  mpz_poly zpol;  // element to descend (n>1)
  mpz_poly f;     // used only in JL mode
  mpz_poly phi;   // used only in the n>1 case
  mpz_t m;        // root of f mod p (the one common with g) when n=1
} Fp_param;

int next_cand_Fp_hgcd(cand_t cand, void *params) {
  pthread_mutex_lock(&mut_die);
  int die = please_die;
  pthread_mutex_unlock(&mut_die);
  if (die) {
    return 0;
  }
  Fp_param * param = (Fp_param *)params;
  unsigned long e = random();
  mpz_t u0, v0;
  mpz_init(u0);
  mpz_init(v0);
  mpz_powm_ui(u0, param->z, e, param->p);
  mpz_t tmp;
  mpz_init(tmp);
  mpz_set(tmp, param->p);
  HalfGcd(u0, tmp, v0);
  mpz_abs(u0, u0);
  mpz_abs(v0, v0);
  cand_set_original_values(cand, u0, v0, e);
  mpz_clear(tmp);
  mpz_clear(u0);
  mpz_clear(v0);
  return 1;
}

int is_probably_sqrfree(mpz_t z)
{
  unsigned long x;
  const int nb_p = 15;
  unsigned long tab_p[15] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47};

  for (int i = 0; i < nb_p; ++i) {
    unsigned long p2 = tab_p[i]*tab_p[i];
    x = mpz_gcd_ui(NULL, z, p2);
    if (x == p2)
      return 0;
  }
  return 1;
}

// JL version
// Returns a boolean meaning "failure, try again".
int get_JL_candiate_from_e(unsigned long e, mpz_poly_ptr UU, mpz_poly_ptr VV,
    mpz_t u, mpz_t v, Fp_param * param)
{
  mpz_t ze;
  mpz_init(ze);
  mpz_powm_ui(ze, param->z, e, param->p);

  //**** Create matrix (see Barbulescu's PhD thesis, section 8.4.2)
  // Warning: the LLL code count indices starting with 1 and not 0.
  // (this is a bug, imho)
  int d = param->f->deg;
  mat_Z M;
  LLL_init(&M, 2*d, 2*d);  // allocate and set to 0.
  // Topleft d*d block
  mpz_set(M.coeff[1][1], param->p);
  for (int i = 1; i < d; ++i) {
    mpz_set_ui(M.coeff[i+1][i+1], 1);
    mpz_neg(M.coeff[i+1][i], param->m);
  }
  // Bottomleft d*d block
  for (int i = 0; i < d; ++i) {
    mpz_set(M.coeff[d+i+1][i+1], ze);
  }
  // Bottomright d*d block
  for (int i = 0; i < d; ++i) {
    mpz_set_ui(M.coeff[d+i+1][d+i+1], 1);
  }

  //**** Apply LLL.
  { 
    mpz_t det, a, b;
    mpz_init(det);
    mpz_init_set_ui (a, 1);
    mpz_init_set_ui (b, 1);
    LLL(det, M, NULL, a, b);
    mpz_clear(det);
    mpz_clear(a);
    mpz_clear(b);
  }

  //**** Recover rational reconstruction
  // z^e = U(alpha)/V(alpha) mod (p, x-m)
  mpz_poly U, V;
  mpz_poly_init(U, d-1);
  mpz_poly_init(V, d-1);
  for (int i = 0; i < d; ++i) {
    mpz_poly_setcoeff(U, i, M.coeff[1][i+1]);
    mpz_poly_setcoeff(V, i, M.coeff[1][d+i+1]);
  }
  mpz_poly_cleandeg(U, d-1);
  mpz_poly_cleandeg(V, d-1);

  //**** Compute norms
  mpz_poly_resultant(u, U, param->f);
  mpz_poly_resultant(v, V, param->f);
  if (UU != NULL)
    mpz_poly_set(UU, U);
  if (VV != NULL)
    mpz_poly_set(VV, V);
  mpz_abs(u, u);
  mpz_abs(v, v);
  // TODO: we don't want to deal with ideals of degree > 1.
  // The simplest way is to forbid squares at all. This is not optimal,
  // and furthermore, we can do it only for small factors...
  int fail = 0;
  if (!is_probably_sqrfree(u) || !is_probably_sqrfree(v))
    fail = 1;

  LLL_clear(&M);
  mpz_poly_clear(U);
  mpz_poly_clear(V);
  mpz_clear(ze);
  return fail;
}


// GF(p^n) version.
// Returns a boolean meaning "failure, try again".
int get_Fpn_candiate_from_e(unsigned long e, mpz_poly_ptr UU, mpz_poly_ptr VV,
    mpz_t u, mpz_t v, Fp_param * param)
{
  int d = param->f->deg;
  int n = param->phi->deg;

  mpz_poly ze;
  mpz_poly_init(ze, d-1);
  mpz_poly_pow_ui_mod_f_mod_mpz(ze, param->zpol, param->f, e, param->p);

  mpz_poly y;           // y represents x^i
  mpz_poly_init(y,0);
  mpz_poly c;           // c represents ze*x^i mod f
  mpz_poly_init(c, d-1);

  //**** Create matrix (see Barbulescu's PhD thesis, section 8.4.2)
  // Warning: the LLL code count indices starting with 1 and not 0.
  // (this is a bug, imho)
  mat_Z M;
  LLL_init(&M, 2*d, 2*d);  // allocate and set to 0.

  // Topleft d*d block
  for (int i = 1; i < n+1; ++i)
    mpz_set(M.coeff[i][i], param->p);
  // Let phi be the degree n polynomial defining the extension field 
  // put the coef of phi=phi[0]+phi[1]*x+...+phi[n-1]*x^(n-1)+x^n
  for (int i = n+1; i < d+1; ++i) {
    for (int j = 0; j < n+1; ++j)
      mpz_set(M.coeff[i][i-n+j],param->phi->coeff[j]);
  }

  // Bottomleft d*d block
  for (int i = 0; i < d; ++i) {
    mpz_poly_set_xi(y,i);     // y is the polynomial alpha^i
    mpz_poly_mul_mod_f_mod_mpz(c, ze, y, param->f, param->p, NULL, NULL);
    for (int j = 0; j < d; ++j) {
      if (j > c->deg)
        mpz_set_ui(M.coeff[d+i+1][j+1], 0);
      else
        mpz_set(M.coeff[d+i+1][j+1], c->coeff[j]);
    }
  }

  // Bottomright d*d block
  for (int i = 0; i < d; ++i) {
    mpz_set_ui(M.coeff[d+i+1][d+i+1], 1);
  }

  //**** Apply LLL.
  {
    mpz_t det, a, b;
    mpz_init(det);
    mpz_init_set_ui (a, 1);
    mpz_init_set_ui (b, 1);
    LLL(det, M, NULL, a, b);
    mpz_clear(det);
    mpz_clear(a);
    mpz_clear(b);
  }

  //**** Recover rational reconstruction
  // z^e = U(alpha)/V(alpha) mod (p, x-m)
  mpz_poly U, V;
  mpz_poly_init(U, d-1);
  mpz_poly_init(V, d-1);
  for (int i = 0; i < d; ++i) {
    mpz_poly_setcoeff(U, i, M.coeff[1][i+1]);
    mpz_poly_setcoeff(V, i, M.coeff[1][d+i+1]);
  }
  mpz_poly_cleandeg(U, d-1);
  mpz_poly_cleandeg(V, d-1);

  //**** Compute norms
  mpz_poly_resultant(u, U, param->f);
  mpz_poly_resultant(v, V, param->f);
  if (UU != NULL)
    mpz_poly_set(UU, U);
  if (VV != NULL)
    mpz_poly_set(VV, V);
  mpz_abs(u, u);
  mpz_abs(v, v);
  // TODO: we don't want to deal with ideals of degree > 1.
  // The simplest way is to forbid squares at all. This is not optimal,
  // and furthermore, we can do it only for small factors...
  int fail = 0;
  if (!is_probably_sqrfree(u) || !is_probably_sqrfree(v))
    fail = 1;

  LLL_clear(&M);
  mpz_poly_clear(U);
  mpz_poly_clear(V);
  mpz_poly_clear(y);
  mpz_poly_clear(c);
  mpz_poly_clear(ze);
  return fail;
}


int next_cand_nonrat(cand_t cand, void *params) {
  pthread_mutex_lock(&mut_die);
  int die = please_die;
  pthread_mutex_unlock(&mut_die);
  if (die) {
    return 0;
  }

  Fp_param * param = (Fp_param *)params;
  mpz_t u, v;
  mpz_init(u);
  mpz_init(v);
  unsigned long e;
  int fail;
  do {
    e = random();
    if (param->n == 1) {
      fail = get_JL_candiate_from_e(e, NULL, NULL, u, v, param);
    } else {
      fail = get_Fpn_candiate_from_e (e, NULL, NULL, u, v, param);
    }
  } while (fail); 

  cand_set_original_values(cand, u, v, e);
  mpz_clear(u);
  mpz_clear(v);
  return 1;
}

int my_mpz_cmp(const void *a, const void *b) {
  mpz_t *pa, *pb;
  pa = (mpz_t *)a;
  pb = (mpz_t *)b;
  return mpz_cmp(pa[0], pb[0]);
}

// Full factorization of z0; non-optimized.
// Assume fac_z has been allocated.
// Returns the number of factors.
int full_factor(mpz_t *fac_z, mpz_t z0) {
  double B1 = 100.0;
  long sig;
  int nf = 0;
  mpz_t z, f;
  mpz_init_set(z, z0);
  mpz_init(f);

  // Remove small primes, ECM can't separate them
  const int nb_p = 15;
  unsigned long tab_p[15] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47};
  for (int i = 0; i < nb_p; ++i) {
    while (mpz_divisible_ui_p(z, tab_p[i])) {
      mpz_set_ui(fac_z[nf], tab_p[i]);
      nf++;
      mpz_divexact_ui(z, z, tab_p[i]);
    }
  }

  while (!mpz_probab_prime_p(z, 10)) {
    int success = 0;
    while (!success) {
      ecm_params ecm_par;
      ecm_init(ecm_par);
      ecm_par->B1done = default_B1done; /* issue with ECM 6.4.x */
      sig = random();
      mpz_set_ui(ecm_par->sigma, sig);
      success = ecm_factor(f, z, B1, ecm_par);
      B1 += sqrt(B1);
      ecm_clear(ecm_par);
    }
    if (mpz_perfect_power_p(f)) {
      mpz_t ff;
      mpz_init(ff);
      for (int pow = 2; ; pow++) {
        if (mpz_root(ff, f, pow)) {
          mpz_set(f, ff);
          break;
        }
      }
      mpz_clear(ff);
    }
    if (!mpz_probab_prime_p(f, 10)) {
      B1 -= 3*sqrt(B1);
      if (B1 <= 20.0)
        B1 = 20.0;
      continue;
    }
    do {
      mpz_set(fac_z[nf], f);
      nf++;
      mpz_divexact(z, z, f);
    } while (mpz_divisible_p(z, f));
  }
  mpz_set(fac_z[nf], z);
  nf++;
  qsort(fac_z, nf, sizeof(mpz_t), my_mpz_cmp);
  mpz_clear(z);
  mpz_clear(f);
  return nf;
}

// Check if there are multiple factors.
// This assumes that the factors are sorted, so that multiple factors are
// consecutive.
int has_distinct_factors(mpz_t *fac_z, int nf) {
  for (int i = 1; i <= nf; ++i) {
    if (mpz_cmp(fac_z[i], fac_z[i-1]) == 0)
      return 0;
  }
  return 1;
}

void find_root(mpz_t r, mpz_t p, mpz_poly f1, mpz_poly f2)
{
  // Check if projective root
  mpz_mod(r, f1->coeff[f1->deg], p);
  if (mpz_cmp_ui(r, 0) ==  0) {
    mpz_mod(r, f2->coeff[f2->deg], p);
    if (mpz_cmp_ui(r, 0) ==  0) {
      mpz_set(r, p);
      return;
    }
  }

  // Non projective
  mpz_poly G;
  mpz_poly_init(G, 1);
  mpz_poly_gcd_mpz(G, f1, f2, p);
  ASSERT_ALWAYS(G->deg == 1);
  mpz_invert(r, G->coeff[1], p);
  mpz_mul(r, r, G->coeff[0]);
  mpz_neg(r, r);
  mpz_mod(r, r, p);
  mpz_poly_clear(G);
}


void print_fac(mpz_t *fac, mpz_t *roots, int nf) {
  for (int i = 0; i < nf - 1; ++i)
    gmp_printf("%Zd,%Zd ", fac[i], roots[i]);
  gmp_printf("%Zd,%Zd", fac[nf-1], roots[nf-1]);
}


void usage(char *argv0) {
  fprintf(stderr, "%s [-poly polfile] [-side xxx] [-extdeg n] [-jl] [-mt n] [-mineff e] [-maxeff E] [-seed s] [-lpb t] [-v] p z\n", argv0);
  fprintf(stderr, "  If extdeg > 1, then z must be a white-separated sequence of coefs z0 z1 ... z_{k-1}\n");
  abort();
}

// Possible modes are
//   MODE_RAT  (when there is a rational side)
//   MODE_JL   (for GF(p), with Joux-Lercier)
//   MODE_FPN  (for GF(p^n))
#define MODE_RAT 1
#define MODE_JL  2
#define MODE_FPN 3
 
typedef struct {
  Fp_param * params;
  const smooth_detect_param_s * smooth_param;
  unsigned long target;
  int mode;
  unsigned int id;
} thparam_s;
typedef thparam_s* thparam_t;

void * process_one_thread(void *thparams) {
  thparam_t thparam = (thparam_s*) thparams;
  cand_t C;
  cand_init(C);
  if (thparam->mode == MODE_RAT) {
    smooth_detect(C, next_cand_Fp_hgcd, thparam->params, thparam->target,
        thparam->smooth_param);
  } else if (thparam->mode == MODE_JL || thparam->mode == MODE_FPN) {
    smooth_detect(C, next_cand_nonrat, thparam->params, thparam->target,
        thparam->smooth_param);
  }

  cand_s *winner = (cand_s*) malloc(sizeof(cand_s));
  cand_init(winner);
  cand_set(winner, C);
  cand_clear(C);
  pthread_mutex_lock(&mut_found);
  thid_found =  thparam->id;
  pthread_cond_signal(&cond_found);
  pthread_mutex_unlock(&mut_found);
  return (void *)winner;
}

int main(int argc, char **argv) {
  char *argv0 = argv[0];
  unsigned long seed = 0;
  unsigned long target = 0; // the target smoothness
  unsigned long nthread = 1;
  double mineff = 2000.0;
  double maxeff = 1e20;
  double minB1 = 100.0;
  int verbose = 0;
  int jl = 0;
  int ext = 1;  //extension degree
  int side = 1;
  cado_poly cpoly;
  int gotpoly = 0;
  clock_t tm = clock();
  
  cado_poly_init(cpoly);

  while (argc > 2 && argv[1][0] == '-') {
    if (strcmp(argv[1], "-seed") == 0) {
      seed = strtoul(argv[2], NULL, 10);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-mt") == 0) {
      nthread = strtoul(argv[2], NULL, 10);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-minB1") == 0) {
      minB1 = strtof(argv[2], NULL);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-mineff") == 0) {
      mineff = strtof(argv[2], NULL);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-side") == 0) {
      side = strtoul(argv[2], NULL, 10);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-maxeff") == 0) {
      maxeff = strtof(argv[2], NULL);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-v") == 0) {
      verbose = 1;
      argc -= 1;
      argv += 1;
    } else if (strcmp(argv[1], "-extdeg") == 0) {
      ext = strtoul(argv[2], NULL, 10);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-jl") == 0) {
      jl = 1;
      argc -= 1;
      argv += 1;
    } else if (strcmp(argv[1], "-lpb") == 0) {
      target = strtoul(argv[2], NULL, 10);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-poly") == 0) {
      if (!cado_poly_read(cpoly, argv[2])) {
        fprintf(stderr, "Error reading polynomial file %s\n", argv[2]);
        usage(argv0);
      }
      gotpoly = 1;
      argc -= 2;
      argv += 2;
    } else {
      fprintf(stderr, "Unknown option %s\n", argv[1]);
      usage(argv0);
    }
  }
  if (argc != 2 + ext)
    usage(argv0);

  if ((jl || (ext > 1)) && !gotpoly) {
    fprintf(stderr, "Error, must provide -poly when extdeg > 1 or using -jl option\n");
    usage(argv0);
  }

  if (ext > 1 && jl) {
    fprintf(stderr, "Warning: ignoring the -jl option with extdeg > 1\n");
    jl = 0;
  }

  if (seed == 0)
    seed = getpid() + (time(NULL)<<16);
  srandom(seed);

  Fp_param params;
  mpz_init_set_str(params.p, argv[1], 10);
  params.n = ext;
  if (jl || ext > 1) {
    mpz_poly_init(params.f, cpoly->pols[side]->deg);
    mpz_poly_set(params.f, cpoly->pols[side]);
  }
  if (jl) {
    mpz_init(params.m);
    int ret = cado_poly_getm(params.m, cpoly, params.p);
    ASSERT_ALWAYS(ret);
  }
  if (ext > 1) {
    mpz_poly_init(params.phi, 0);
    mpz_poly_gcd_mpz(params.phi, cpoly->pols[0], cpoly->pols[1], params.p);
    mpz_poly_init(params.zpol, ext-1);
    for(int i = 0; i < ext; i ++) {
      mpz_t tmp;
      mpz_init_set_str(tmp, argv[i+2],10);
      mpz_poly_setcoeff(params.zpol, i, tmp);
      mpz_clear(tmp);
    }
  }

  if (ext == 1) {
    mpz_init_set_str(params.z, argv[2], 10);
  }
  
  const smooth_detect_param_s smooth_param = {mineff, maxeff, 10, verbose, minB1};
  
  thparam_s * thparam;
  thparam = (thparam_s *)malloc(nthread*sizeof(thparam_s));

  pthread_t *thid = malloc(nthread*sizeof(pthread_t));
  for (unsigned int i = 0; i < nthread; ++i) {
    thparam[i].params = &params;
    thparam[i].smooth_param = &smooth_param;
    thparam[i].target = target;
    if (ext > 1) {
      thparam[i].mode = MODE_FPN;
    } else {
      if (jl)
        thparam[i].mode = MODE_JL;
      else
        thparam[i].mode = MODE_RAT;
    }
    thparam[i].id = i;
    int ret = pthread_create(&thid[i], NULL, process_one_thread, &thparam[i]);
    ASSERT_ALWAYS(ret == 0);
  }

  // Wait for one thread to signal success.
  pthread_mutex_lock(&mut_found);
  pthread_cond_wait(&cond_found, &mut_found);
  unsigned int i = thid_found;
  pthread_mutex_unlock(&mut_found);
  
  int found = 0;
  do {
    cand_s *winner;
    int ret = pthread_join(thid[i], (void *)(&winner));
    ASSERT_ALWAYS(ret == 0);
    if (jl || ext > 1) {
      mpz_t u, v;
      mpz_init(u); mpz_init(v);
      mpz_poly U, V;
      mpz_poly_init(U, params.f->deg-1);
      mpz_poly_init(V, params.f->deg-1);
      // do again the LLL thing and print result.
      int fail;
      if (ext == 1)
        fail = get_JL_candiate_from_e(winner->id, U, V, u, v, &params);
      else
        fail = get_Fpn_candiate_from_e(winner->id, U, V, u, v, &params);
      ASSERT_ALWAYS(!fail);

      mpz_t facu[50];
      mpz_t facv[50];
      mpz_t roots[50];
      int nfu, nfv;
      for (int j = 0; j < 50; ++j) {
        mpz_init(facu[j]);
        mpz_init(facv[j]);
        mpz_init(roots[j]);
      }
      nfu = full_factor(facu, u);
      nfv = full_factor(facv, v);
      if (!has_distinct_factors(facu, nfu) || !has_distinct_factors(facv, nfv)) {
        // one of them is not squarefree. Restart the thread and wait for
        // another candidate.
        ret = pthread_create(&thid[i], NULL, process_one_thread, &thparam[i]);
        ASSERT(ret == 0);
        pthread_mutex_lock(&mut_found);
        pthread_cond_wait(&cond_found, &mut_found);
        i = thid_found;
        pthread_mutex_unlock(&mut_found);
        printf("Fail: non-squarefree norm\n");
        mpz_clear(u); mpz_clear(v);
        mpz_poly_clear(U); mpz_poly_clear(V);
        for (int j = 0; j < 50; ++j) {
          mpz_clear(facu[j]);
          mpz_clear(facv[j]);
          mpz_clear(roots[j]);
        }
        cand_clear(winner);
        free(winner);
        continue;
      }

      printf("U = ");
      mpz_poly_fprintf_coeffs(stdout, U, ',');
      printf("V = ");
      mpz_poly_fprintf_coeffs(stdout, V, ',');
      gmp_printf("u = %Zd\nv = %Zd\n", u, v);

      for (int j = 0; j < nfu; ++j)
        find_root(roots[j], facu[j], U, params.f);
      printf("fac_u = "); print_fac(facu, roots, nfu); printf("\n");

      for (int j = 0; j < nfv; ++j)
        find_root(roots[j], facv[j], V, params.f);
      printf("fac_v = "); print_fac(facv, roots, nfv); printf("\n");

      for (int j = 0; j < 50; ++j) {
        mpz_clear(facu[j]);
        mpz_clear(facv[j]);
        mpz_clear(roots[j]);
      }
      mpz_clear(u); mpz_clear(v);
      mpz_poly_clear(U); mpz_poly_clear(V);
    }
    printf("Youpi: e = %lu is a winner\n", winner->id);
    cand_clear(winner);
    free(winner);

    found = 1;
  } while (!found);

  // Cancel other threads
  pthread_mutex_lock(&mut_die);
  please_die = 1;
  pthread_mutex_unlock(&mut_die);
  for (unsigned int j = 0; j < nthread; ++j) {
    if (i == j)
      continue;
    cand_s *winner;
    pthread_join(thid[j], (void *)(&winner));
    cand_clear(winner);
    free(winner);
  }

  printf("Total CPU time: %.1f s\n", ((double)(clock() - tm)) / CLOCKS_PER_SEC);
  free(thid);
  free(thparam);
  if (gotpoly)
    cado_poly_clear(cpoly);
  if (jl || ext > 1) 
    mpz_poly_clear(params.f);
  if (jl) {
    mpz_clear(params.m);
  }
  if (ext > 1) {
    mpz_poly_clear(params.phi);
    mpz_poly_clear(params.zpol);
  } else {
    mpz_clear(params.z);
  }

  mpz_clear(params.p);
  return EXIT_SUCCESS;
} 
