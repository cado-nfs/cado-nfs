#define _DEFAULT_SOURCE  // for random()
#define _GNU_SOURCE // for pthread_tryjoin_np()
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/types.h> 
#include <sys/resource.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <sys/time.h>
#include <gmp.h>

#include "smooth_detect.h"
#include "macros.h"
#include "portability.h"

/*
 * Look for a bi-smooth pair among a set of pre-sieved pairs.
 *
 * Input must be given in a file, one per line, each line of the form:
 *   a b cofac_a cofac_b q
 * (only b, cofac_a, cofac_b are read)
 *
 * b is used as an identifier of the winner line.
 *
 */


/*
ECM=...
gcc -Wall -g -std=c99 -I${ECM}/include filter_sieved_pairs.c smooth_detect.c -o filter_sieved_pairs -L${ECM}/lib -lecm -lgmp -lm -lpthread
*/

typedef struct {
  pthread_mutex_t mutex;
  FILE *file;
} sieved_params;

int next_cand_sieved_pair(cand_t cand, void *params) {
  sieved_params *param = (sieved_params *)params;
  int ret = pthread_mutex_lock(&param->mutex);
  ASSERT(ret == 0);
  unsigned long id;
  mpz_t u, v, u0;
  mpz_init(u);
  mpz_init(v);
  // skip commented lines
  int c = fgetc(param->file);
  ASSERT (c != EOF);
  while (c == '#') {
    char str[2048];
    fgets(str, 2048, param->file);
    c = fgetc(param->file);
    ASSERT (c != EOF);
  }
  ungetc(c, param->file);
  // read the line
  ret = gmp_fscanf(param->file, "%*Zd %lu %Zd %Zd %*Zd\n", &id, u, v);
  ASSERT (ret == 3);
  mpz_init_set_ui(u0, 1);
  mpz_mul_2exp(u0, u0, 512);
  cand_set_presieved_values(cand, u0, u0, u, v, 0, 0, id);
  mpz_clear(u);
  mpz_clear(v);
  mpz_clear(u0);
  pthread_mutex_unlock(&param->mutex);
  return 1;
}

void usage(char *argv0) {
  fprintf(stderr, "./%s [-mt n] [-minB1 b1] [-mineff e] [-maxeff E] [-seed s] -target t <filename>\n", argv0);
  abort();
}

typedef struct {
  sieved_params * params;
  const smooth_detect_param_s * smooth_param;
  unsigned long target;
} thparam_s;
typedef thparam_s* thparam_t;

void * process_one_thread(void *thparams) {
  thparam_t thparam = (thparam_s*) thparams;
  cand_t C;
  cand_init(C);
  smooth_detect(C, next_cand_sieved_pair, thparam->params, thparam->target,
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
  double minB1 = 100.0;
  int verbose = 0;
  clock_t tm = clock();
  
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
    } else if (strcmp(argv[1], "-minB1") == 0) {
      minB1 = strtof(argv[2], NULL);
      argc -= 2;
      argv += 2;
    } else if (strcmp(argv[1], "-v") == 0) {
      verbose = 1;
      argc -= 1;
      argv += 1;
    } else if (strcmp(argv[1], "-target") == 0) {
      target = strtoul(argv[2], NULL, 10);
      argc -= 2;
      argv += 2;
    } else {
      fprintf(stderr, "Unknown option %s\n", argv[1]);
      usage(argv0);
    }
  }
  if (argc != 2)
    usage(argv0);

  if (seed == 0)
    seed = getpid() + (time(NULL)<<16);
  srandom(seed);

  sieved_params params[1];
  params->file = fopen(argv[1], "r");
  if (!params->file) {
    fprintf(stderr, "Could not open %s for reading: ", argv[1]);
    perror("\n");
    exit(EXIT_FAILURE);
  }
  pthread_mutex_init(&params->mutex, NULL);
  
  const smooth_detect_param_s smooth_param = {mineff, maxeff, 10, verbose, minB1};

  thparam_s thparam[1];
  thparam->params = &params[0];
  thparam->smooth_param = &smooth_param;
  thparam->target = target;

  pthread_t *thid = malloc(nthread*sizeof(pthread_t));
  for (unsigned int i = 0; i < nthread; ++i) {
    int ret = pthread_create(&thid[i], NULL, process_one_thread, thparam);
    ASSERT(ret == 0);
  }
  
  while (1) {
    for (unsigned int i = 0; i < nthread; ++i) {
      cand_s *winner;
      // Non-blocking join: non-portable, gnu-specific.
      int ret = pthread_tryjoin_np(thid[i], (void *)(&winner));
      if (ret == 0) {
        printf("Youpi: b = %lu is a winner\n", winner->id);
        goto end;
      }
    }
    sleep(1);
  }

end:
  printf("Total CPU time: %.0f s\n", ((double)(clock() - tm)) / CLOCKS_PER_SEC);
  free(thid);
  fclose(params->file);
  return EXIT_SUCCESS;
}
