#define _DEFAULT_SOURCE  // for random()
#define _GNU_SOURCE // for pthread_tryjoin_np()
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <sys/types.h> 
#include <sys/resource.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <gmp.h>

#include "smooth_detect.h"

/*
ECM=...
gcc -Wall -g -std=c99 -I${ECM}/include descent_init_Fp.c smooth_detect.c -o descent_init_Fp -L${ECM}/lib -lecm -lgmp -lm -lpthread
*/


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

// Function and parameters to be passed to smooth_detect.
typedef struct {
  mpz_t p;
  mpz_t z;
} Fp_hgcd_param;

void next_cand_Fp_hgcd(cand_t cand, void *params) {
  Fp_hgcd_param * param = (Fp_hgcd_param *)params;
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
}

void usage(char *argv0) {
  fprintf(stderr, "./%s [-mt n] [-mineff e] [-maxeff E] [-seed s] [-target t] p z\n", argv0);
  abort();
}

typedef struct {
  Fp_hgcd_param * params;
  const smooth_detect_param_s * smooth_param;
  unsigned long target;
} thparam_s;
typedef thparam_s* thparam_t;

void * process_one_thread(void *thparams) {
  thparam_t thparam = (thparam_s*) thparams;
  cand_t C;
  cand_init(C);
  smooth_detect(C, next_cand_Fp_hgcd, thparam->params, thparam->target,
      thparam->smooth_param);
  cand_s *winner = (cand_s*) malloc(sizeof(cand_s));
  cand_init(winner);
  cand_set(winner, C);
  return (void *)winner;
}

int main(int argc, char **argv) {
  char *argv0 = argv[0];
  unsigned long seed = 0;
  unsigned long target = 0;
  unsigned long nthread = 1;
  double mineff = 2000.0;
  double maxeff = 1e20;
  
  while (argc > 2 && argv[1][0] == '-') {
    if (strcmp(argv[1], "-seed") == 0) {
      seed = strtoul(argv[2], NULL, 10);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-mt") == 0) {
      nthread = strtoul(argv[2], NULL, 10);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-mineff") == 0) {
      mineff = strtof(argv[2], NULL);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-maxeff") == 0) {
      maxeff = strtof(argv[2], NULL);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-target") == 0) {
      target = strtoul(argv[2], NULL, 10);
      argc -= 2;
      argv += 2;
    } else {
      fprintf(stderr, "Unknown option %s\n", argv[1]);
      usage(argv0);
    }
  }
  if (argc != 3)
    usage(argv0);

  if (seed == 0)
    seed = getpid() + (time(NULL)<<16);
  srandom(seed);

  Fp_hgcd_param params;
  mpz_init_set_str(params.p, argv[1], 10);
  mpz_init_set_str(params.z, argv[2], 10);
  
  const smooth_detect_param_s smooth_param = {mineff, maxeff, 10};

  thparam_s thparam[1];
  thparam->params = &params;
  thparam->smooth_param = &smooth_param;
  thparam->target = target;

  pthread_t *thid = malloc(nthread*sizeof(pthread_t));
  for (unsigned int i = 0; i < nthread; ++i) {
    int ret = pthread_create(&thid[i], NULL, process_one_thread, thparam);
    assert(ret == 0);
  }
  
  while (1) {
    for (unsigned int i = 0; i < nthread; ++i) {
      cand_s *winner;
      // Non-blocking join: non-portable, gnu-specific.
      int ret = pthread_tryjoin_np(thid[i], (void *)(&winner));
      if (ret == 0) {
        printf("Youpi: e = %lu is a winner\n", winner->id);
        goto end;
      }
    }
    sleep(1);
  }

end:
  free(thid);
  mpz_clear(params.p);
  mpz_clear(params.z);
  return EXIT_SUCCESS;
}
