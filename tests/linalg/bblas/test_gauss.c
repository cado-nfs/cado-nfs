#include "cado.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "macros.h"
#include "gauss.h"

#ifndef NB_TEST
#define NB_TEST 100
#endif

#ifdef VERBOSE
static void printVector(mp_limb_t *V, int n) {
  int i;
  printf("[");
  for (i = 0; i < n; ++i) {
    if (V[i / ULONG_BITS] & (1UL << (i % ULONG_BITS)))
      printf("1");
    else
      printf("0");
    if (i < n - 1) printf(",");
  }
  printf("]");
}

static void printMatrix(mp_limb_t *M, int nrows, int ncols,
			int limbs_per_row) {
  int i;
  printf("[\n");
  for (i = 0; i < nrows; ++i) {
    printVector(M + i*limbs_per_row, ncols);
    if (i < nrows - 1)
      printf(",\n");
  }
  printf("\n]\n");
}
#endif

int main(int argc, char **argv) {
  mp_limb_t* mat;
  mp_limb_t** ker;
  int nrows, ncols, limbs_per_row, limbs_per_col;
  int i, justrank = 0;
#ifdef VERBOSE
  int j;
#endif
  
  nrows = ncols = 0;

  argc--,argv++;
  for( ; argc ; argc--,argv++) {
      int u = atoi(*argv);
      if (u) {
          if (nrows == 0) {
              nrows = ncols = u;
          } else {
              ncols = u;
          }
      } else if (strcmp(argv[2], "justrank") == 0) {
          justrank = 1;
      } else {
	printf("usage: %s n [justrank]\n", argv[0]);
	printf("   where n is the size of the matrix\n");
	exit(1);
      }
  }
  if (nrows == 0) nrows = ncols = 163;

  limbs_per_row = iceildiv(ncols, ULONG_BITS);
  limbs_per_col = iceildiv(nrows, ULONG_BITS);

  mat = (mp_limb_t *)malloc(limbs_per_row*nrows*sizeof(mp_limb_t));
  ker = (mp_limb_t **)malloc(nrows*sizeof(mp_limb_t *));
  for (i = 0; i < nrows; ++i)
    ker[i] = (mp_limb_t *)malloc(limbs_per_col*sizeof(mp_limb_t));

  for (i = 0; i < NB_TEST; ++i) {
    int dim MAYBE_UNUSED;

    mpn_random(mat, limbs_per_row*nrows);

#ifdef VERBOSE
    printf("M:=\n");
    printMatrix(mat, nrows, ncols, limbs_per_row);
    printf(";\n");
#endif


#ifdef VERBOSE
    mp_limb_t * s = malloc(nrows * limbs_per_col * sizeof(mp_limb_t));
    dim = spanned_basis(s, mat, nrows, ncols, limbs_per_row, limbs_per_col);
    printf("rank:=%d\n", dim);
    printf("S:=\n");
    printMatrix(s,nrows,nrows,limbs_per_col);
    printf(";\n");
    free(s);
#endif


    if (!justrank)
      dim = kernel(mat, ker, nrows, ncols, limbs_per_row, limbs_per_col);
    else
      dim = kernel(mat, NULL, nrows, ncols, limbs_per_row, limbs_per_col);

#ifdef VERBOSE
    printf("dimker:=%d;\n", dim);

    if ((!justrank) && (dim > 0)) 
      printf("Bker:=[\n");
      for (j = 0; j < dim; ++j) {
	printVector(ker[j], nrows);
        if (j < dim-1) printf(",");
        printf("\n");
      }
      printf("];\n");
#endif

  }

  free(mat);
  free(ker);
}

