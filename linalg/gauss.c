/* gauss.c ; solve small linear systems over GF(2)
 * 
 * Algorithm: Gaussian elimination as described in ...
   
   Copyright 2007, 2008, 2009, 2010 Pierrick Gaudry, Francois Morain, Emmanuel Thom\'e, Paul Zimmermann
   
   This file is part of CADO-NFS.
   
   CADO-NFS is free software; you can redistribute it and/or modify it
   under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation; either version 2.1 of the License, or
   (at your option) any later version.
   
   CADO-NFS is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more
   details.
   
   You should have received a copy of the GNU Lesser General Public
   License along with CADO-NFS; see the file COPYING.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
   Boston, MA 02110-1301, USA.
*/


/*
Compile with one of these:

gcc -I../ -I../utils -I../build/x86_64 -L../build/x86_64/utils gauss.c -lgmp -lutils -lm

gcc -O4 -DNDEBUG -funroll-loops -I../ -I../utils -I../build/x86_64 -L../build/x86_64/utils -DMULTI_ROW=3 -DNB_TEST=10 gauss.c -lgmp -lutils -lm


Note: for small sizes (100-200), MULTI_ROW=2 seems to be fine,
      for larger sizes (also for small sizes on alpha), try MULTI_ROW=3

      With NO_MAIN option, a .o file is generated, exporting only the
      kernel function.  In that case, the code does not use GMP at all.
      
The exported function is the following:

 int kernel(mp_limb_t* mat, mp_limb_t** ker, int nrows, int ncols,
            int limbs_per_row, int limbs_per_col);

which returns the dimension of the kernel and fill-in the pointer ker.
If just the dimension of kernel is wanted, set ker=NULL.

*/

/*===========================================================================*/
/*    includes, defines, and declarations...                                 */
/* Warning: there are some (static) global variables                         */
/*===========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gmp.h"   /* only used for setting a random matrix */
#include "utils.h"  /* for seconds() */
#include "gauss.h"

#ifndef VERBOSE
#define VERBOSE 0
#endif
#ifndef NB_TEST
#define NB_TEST 100
#endif
#ifndef MULTI_ROW
#define MULTI_ROW 2
#endif


/* Look on which computer we are */
#define MACHINE_WORD_SIZE GMP_LIMB_BITS

/* If the flag below is set, the interface to kernel() changes slightly, as the
 * mp_limb_t ** ker argument is expected to hold space for pointers BUT
 * these pointers need not be initialized to anything. The number of pointers
 * eventually returned, and pointing to data, is exactly the dimension of
 * the kernel.
 *
 * For tiny matrices, this feature does not make sense, as a flat storage
 * for the kernel vectors is enough.
 */
#define SAVE_KERNEL_MEMORY 0

#ifndef __GMP_H__
typedef unsigned long int mp_limb_t;
#endif

/* some global variables storing the data */
static int NROWS = 0;
static int NCOLS = 0;
static int LIMBS_PER_ROW = 0;
static mp_limb_t* matrix = NULL;
static mp_limb_t** ptr_rows = NULL;

/* functions declaration */

static void check_soundness();

#if 0 /* ununsed currently */
static void addPartialRows(int row, int pivot, mp_limb_t mask,
			   int j_cur, mp_limb_t **ptr);
#endif
static inline void addRows(int row, int pivot, mp_limb_t mask,
			   /*int j_cur,*/ mp_limb_t **ptr);
static inline void add2Rows(int row, int row2, int pivot, mp_limb_t mask,
			    /*int j_cur,*/ mp_limb_t **ptr_current);
static inline void add3Rows(int row, int row2, int row3, int pivot,
			    mp_limb_t mask, /* int j_current, */
			    mp_limb_t **ptr_current);
static inline int getPivot(mp_limb_t **ptr, mp_limb_t mask);
#if VERBOSE
#ifndef NO_MAIN
static void printVector(mp_limb_t *V, int n);
static void printMatrix(mp_limb_t *M, int nrows, int ncols, int limbs_per_row);
#endif
#endif


/*===========================================================================*/
/*  functions definition                                                     */
/*===========================================================================*/

#ifndef NO_MAIN
int main(int argc, char **argv) {
  mp_limb_t* mat;
  mp_limb_t** ker;
  int nrows, ncols, limbs_per_row, limbs_per_col;
  int i, justrank = 0;
#if VERBOSE
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

  limbs_per_row = iceildiv(ncols, MACHINE_WORD_SIZE);
  limbs_per_col = iceildiv(nrows, MACHINE_WORD_SIZE);

  mat = (mp_limb_t *)malloc(limbs_per_row*nrows*sizeof(mp_limb_t));
  ker = (mp_limb_t **)malloc(nrows*sizeof(mp_limb_t *));
  for (i = 0; i < nrows; ++i)
    ker[i] = (mp_limb_t *)malloc(limbs_per_col*sizeof(mp_limb_t));

  for (i = 0; i < NB_TEST; ++i) {
    int dim;
    mpn_random(mat, limbs_per_row*nrows);

#if 0
    // special case for a 1024x128 matrix.
    {
        int i;
        FILE * f;
        f = fopen("/tmp/cado.ZxyM8F58ky//c59.bwc/K.1", "r");
        for(i = 0 ; i < nrows ; i++) {
            fread(&mat[i * limbs_per_row], sizeof(uint64_t), 1, f);
        }
        fclose(f);
        f = fopen("/tmp/cado.ZxyM8F58ky//c59.bwc/K.0", "r");
        for(i = 0 ; i < nrows ; i++) {
            fread(&mat[1 + i * limbs_per_row], sizeof(uint64_t), 1, f);
        }
        fclose(f);
    }
#endif


#if VERBOSE
    printf("M:=\n");
    printMatrix(mat, nrows, ncols, limbs_per_row);
    printf(";\n");
#endif


#if VERBOSE
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

    printf("dimker:=%d;\n", dim);

#if VERBOSE
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
#endif

static void check_soundness(){
#ifndef NDEBUG
  mp_limb_t x = 2342326;
  int i;
  ASSERT (LIMBS_PER_ROW*MACHINE_WORD_SIZE >= NCOLS);
  /* touch everywhere in the matrix to look for SEGV's */
  for (i = 0; i < NROWS; ++i) 
    x += matrix[LIMBS_PER_ROW*i];
  if (x == matrix[LIMBS_PER_ROW*NROWS-1])
    ASSERT(0);  /* very unlikely to happen, but forces no optimizations */
#endif
}

void simple_addrows(mp_limb_t * p, int d, int s, int stride)
{
    if (!p) return;
    mp_limb_t * src = p + s * stride;
    mp_limb_t * dst = p + d * stride;
    int j;
    for(j = 0 ; j < stride ; j++) {
        dst[j] ^= src[j];
    }
}

/* 
 * Kernel computation.
 *
 * Put the vectors of the kernel in the (preallocated) table of cols ker.
 * Return the dimension of the kernel.
 * If just the dimension is wanted, then pass ker = NULL, and in that case,
 * limbs_per_cols is meaningless.
 * 
 * Variables description:
 *   ptr_current : table containing for each row the address of
 *                 the limb which we are currently dealing with.
 *                 These are incremented after MACHINE_WORD_SIZE 
 *                 pivot elimination.
 *   col_current : index of the current column
 *   j_current   : the index of the limb containing the current column
 *   mask1       : mask with only one 1 set at the current column
 *   mask2       : mask with 1's set at column higher than the current
 *   set_pivot   : table. at index j, contains the index of the row that
 *                 used for pivoting column j. Initialized with -1.
 *   set_used    : table. at index i, contains 1 if the i-th row has
 *                 been used for pivoting, at some point, 0 otherwise.
 *   rank        : the rank of the matrix
 *   i, j        : used in for(;;) loops
 */

struct gaussian_elimination_data {
    int rank;
    int nrows;
    int ncols;
    mp_limb_t* mat;
    mp_limb_t* lmat;
    int limbs_per_row;
    int limbs_per_col;
    int *set_pivot;
    int *set_used;
};

#define ADDBIT_SLOW(r, l, i, j)    \
    r[i * l + j / MACHINE_WORD_SIZE] ^= 1UL << (j % MACHINE_WORD_SIZE)

void gaussian_elimination(struct gaussian_elimination_data * G)
{
  int i;

  /* store the data in the global variables */
  matrix = G->mat;
  NROWS = G->nrows;
  NCOLS = G->ncols;
  LIMBS_PER_ROW = G->limbs_per_row;
  int limbs_per_col = G->limbs_per_col;
  ptr_rows = (mp_limb_t **)malloc((NROWS+1)*sizeof(mp_limb_t *));
  for (i = 0; i < NROWS+1; ++i)
    ptr_rows[i] = matrix + LIMBS_PER_ROW*i;


  /* if G->lmat == NULL, then don't bug the user if he just gave 0 for the
   * otherwise unused limbs_per_col value */
  assert (G->lmat == NULL || limbs_per_col*MACHINE_WORD_SIZE >= NROWS);

  G->set_pivot = (int *)malloc(NCOLS*sizeof(int));
  G->set_used  = (int *)malloc(NROWS*sizeof(int));

  for (i = 0; i < NCOLS; ++i)
    G->set_pivot[i] = -1;
  for (i = 0; i < NROWS; ++i)
    G->set_used[i] = 0;

  mp_limb_t **ptr_current;
  int j_current, col_current;
  mp_limb_t mask1, mask2;
  double st0 = seconds(), st;

  // shut up.
#if VERBOSE
  fprintf (stderr, "Using MULTI_ROW=%d\n", MULTI_ROW);
#endif

  check_soundness();
  
  /* initialize ptr_current, set_used and set_pivot */
  ptr_current = (mp_limb_t **)malloc((NROWS)*sizeof(mp_limb_t *));
  for (i = 0; i < NROWS; ++i)
    ptr_current[i] = ptr_rows[i];

  /* Main Loop: for each column, find the pivot and eliminate */
  j_current = 0;
  mask1 = 1UL;
  mask2 = (-(1UL)) ^ (1UL);
  col_current = 0;
  while (col_current < NCOLS) {
    int pivot;

    /* Get the pivot */
    pivot = getPivot(ptr_current, mask1);
    if (pivot != -1) {
#ifndef NDEBUG
      for (i = 0; i < pivot; ++i)
	assert (!((*ptr_current[i]) & mask1 ));
      assert (((*ptr_current[pivot]) & mask1 ));
      assert (G->set_used[pivot] == 0);
#endif
      G->set_pivot[col_current] = pivot;
      G->set_used[pivot] = 1;

      /* Eliminate with the pivot */
#if MULTI_ROW == 1
      for (i = pivot + 1; i < NROWS; ++i) {
	/* is there a 1 in position (i, col_current) ? */
	if ( (*ptr_current[i]) & mask1 ) {
	  addRows(i, pivot, mask1, /*j_current,*/ ptr_current);
          simple_addrows(G->lmat, i, pivot, limbs_per_col);
        }
      }      
#elif MULTI_ROW == 2
      /* try with two rows at a time */
      { 
	int i1;
	i = pivot + 1;
	i1 = -1;
	while (i < NROWS) {
	  if ( (*ptr_current[i]) & mask1 ) {
	    if (i1 >= 0) {
	      add2Rows(i1, i, pivot, mask1, /*j_current,*/ ptr_current);
              simple_addrows(G->lmat, i, pivot, limbs_per_col);
              simple_addrows(G->lmat, i1, pivot, limbs_per_col);
	      i1 = -1;
	    } else {
	      i1 = i;
	    }
	  } /* end if */
	  ++i;
	} /* end while */
	if (i1 >= 0)  {
	  addRows(i1, pivot, mask1, /*j_current,*/ ptr_current);
          simple_addrows(G->lmat, i1, pivot, limbs_per_col);
        }
      } /* end block */
#elif MULTI_ROW == 3
      /* try with three rows at a time */
      { 
	int i1, i2;
	i = pivot + 1;
	i1 = -1;  i2 = -1;
	while (i < NROWS) {
	  if ( (*ptr_current[i]) & mask1 ) {
	    if (i1 >= 0) {
	      if (i2 >= 0) {
		add3Rows(i1, i2, i, pivot, mask1, /*j_current,*/ ptr_current);
                simple_addrows(G->lmat, i, pivot, limbs_per_col);
                simple_addrows(G->lmat, i1, pivot, limbs_per_col);
                simple_addrows(G->lmat, i2, pivot, limbs_per_col);
		i1 = -1; i2 = -1;
	      }
	      else {
		i2 = i;
	      }
	    } else {
	      i1 = i;
	    }
	  } /* end if */
	  ++i;
	} /* end while */
	if (i1 >= 0)  {
	  if (i2 >= 0) {
	    add2Rows(i1, i2, pivot, mask1, /*j_current,*/ ptr_current);
            simple_addrows(G->lmat, i1, pivot, limbs_per_col);
            simple_addrows(G->lmat, i2, pivot, limbs_per_col);
          } else {
	    addRows(i1, pivot, mask1, /*j_current,*/ ptr_current);
            simple_addrows(G->lmat, i1, pivot, limbs_per_col);
          }
	}
      } /* end block */
#else
#error MULTI_ROW must be 1, 2, or 3.
#endif 

      /* Purge the pivot line */
      addRows(pivot, pivot, mask1, /*j_current,*/ ptr_current);
#if 0
#if VERBOSE
      if (G->lmat) {
          printf("L:=\n");
          printMatrix(G->lmat, NROWS, NROWS, limbs_per_col);
          printf(";\n");
      }
#endif
#endif
    }

    /* increment */
    ++col_current;
    mask1 <<= 1;
    mask2 <<= 1;
    if (!mask1) { /* we have done MACHINE_WORD_SIZE operations */
      mask1 = 1UL;
      mask2 = (-(1UL)) ^ (1UL);
      j_current++;
      assert (j_current <= LIMBS_PER_ROW);
      for (i = 0; i < NROWS; ++i)
	ptr_current[i]++;
    }      

    /* some verbosity... */
    if ((col_current % 128) == 0)
      {
        st = (seconds () - st0); /* time in seconds */
#if VERBOSE
        fprintf (stderr, "done %d pivots in %1.0fs (est. %1.0fs)\n",
			col_current,
                 st, (double) NCOLS * st / (double) col_current);
#endif
      }

  } /* end while */
  free(ptr_current);
}


int kernel(mp_limb_t* mat, mp_limb_t** ker, int nrows, int ncols,
	   int limbs_per_row, int limbs_per_col)
{
    int i,j;

    struct gaussian_elimination_data G[1] = {{
        .mat = mat,
        .nrows = nrows, .ncols = ncols,
        .limbs_per_row = limbs_per_row,
        .limbs_per_col = limbs_per_col,
        .rank = 0, }};
    gaussian_elimination(G);

#if VERBOSE
  printf("Pivots:=[");
  for (i = 0; i < NCOLS-1; ++i)
    printf("%d,", G->set_pivot[i]);
  printf("%d];\n", G->set_pivot[NCOLS-1]);

  printf("usedrows:= [");
  for (i = 0; i < NROWS-1; ++i)
    printf("%d,", G->set_used[i]);
  printf("%d];\n", G->set_used[NROWS-1]);
#endif     

 
  /* if ker == NULL, then don't bug the user if he just gave 0 for the
   * otherwise unused limbs_per_col value */
  assert (ker == NULL || limbs_per_col*MACHINE_WORD_SIZE >= NROWS);
  /* Explore the set of unused rows, get the G->rank */
  G->rank = NROWS;
  for (i = 0; i < NROWS; ++i)
    if (!G->set_used[i]) {
      G->rank--;
      if (ker != NULL) {
	/* each unused row gives one element of the kernel */
#if SAVE_KERNEL_MEMORY
	ker[0] = (mp_limb_t *) malloc(limbs_per_col * sizeof(mp_limb_t));
#endif

	/* ker[0] contains the space for the new element of the kernel */
	mp_limb_t *ptr = ptr_rows[i];
	
	for (j = 0; j < limbs_per_col; ++j)
	  ker[0][j] = 0;
	ker[0][i / MACHINE_WORD_SIZE] = 1UL << (i % MACHINE_WORD_SIZE);

	mp_limb_t mask1 = 1UL;
	j = 0;
	while (j < NCOLS) {
	  if ((*ptr) & mask1) { /* have to recover the pivoting */
	    int pivot;
	    pivot = G->set_pivot[j];
	    (ker[0])[pivot / MACHINE_WORD_SIZE] |= 
	      (1UL << (pivot % MACHINE_WORD_SIZE));
	  } /* end if */
	  /* increment */
	  ++j;
	  mask1 <<= 1;
	  if (!mask1) { /* we have done MACHINE_WORD_SIZE operations */
	    mask1 = 1UL;
	    ptr++;
	  }
	} /* end while */
	ker++;
      } /* end if ker != NULL */
    } /* end if / for, loop over all unused */

  free(ptr_rows);
  free(G->set_pivot);
  free(G->set_used);

  return NROWS - G->rank;
} /* end function kernel */

/* puts in lmat a full rank matrix of size nrows*nrows, and return a rank
 * value r, such that the r first (least significant bits) rows of the
 * matrix lmat * mat form a basis of the row span of the matrix mat,
 * while the last (nrows-r) (most significant bits) rows are zero.
 *
 * lmat is expected to be already allocated, and to have a striding equal
 * to limbs_per_col.
 *
 * elim_table is an extra field, which should only be used when ncols is
 * small. If not null, one expects there a malloc()'ed area (see size
 * below), where each block of 16*limbs_per_col limbs contains
 * precomputed row combinations. Row combination of index i, in the k-th
 * such block, is such that the multiplication of this row by the
 * original matrix mat yields a row whose bits of indices k to k+4 are
 * exactly set to i. This can be used to do a projection.
 */
int spanned_basis(mp_limb_t * lmat, mp_limb_t * mat, int nrows, int ncols,
        int limbs_per_row, int limbs_per_col, mp_limb_t * elim_table)
{
    int i;

#if 0
#if VERBOSE
    printf("M:=\n");
    printMatrix(mat, nrows, ncols, limbs_per_row);
    printf(";\n");
#endif
#endif
    mp_limb_t * lmat2 = malloc(nrows * limbs_per_col * sizeof(mp_limb_t));
    memset(lmat2, 0, nrows * limbs_per_col * sizeof(mp_limb_t));
    for (i = 0; i < nrows; ++i)
        ADDBIT_SLOW(lmat2, limbs_per_col, i, i);

    struct gaussian_elimination_data G[1] = {{
            .mat = mat, .lmat = lmat2,
            .nrows = nrows, .ncols = ncols,
            .limbs_per_row = limbs_per_row,
            .limbs_per_col = limbs_per_col,
            .rank = 0, }};
    gaussian_elimination(G);

    int r = 0;
    int h = NROWS;

    for (i = 0; i < NROWS; ++i) {
        assert(r < h);
        if (G->set_used[i]) {
            memcpy(lmat + r * limbs_per_col, lmat2 + i * limbs_per_col, limbs_per_col * sizeof(mp_limb_t));
            r++;
        } else {
            h--;
            memcpy(lmat + h * limbs_per_col, lmat2 + i * limbs_per_col, limbs_per_col * sizeof(mp_limb_t));
        }
    }

    if (elim_table) {
        // clearing is done progressively below
        // memset(elim_table, 0, iceildiv(NCOLS, 4) * 16 * limbs_per_col);
        mp_limb_t * zcol = malloc(limbs_per_col * sizeof(mp_limb_t));
        memset(zcol, 0, limbs_per_col * sizeof(mp_limb_t));
        mp_limb_t * e = elim_table;
        for (i = 0; i < NCOLS ; i+=4, e+=16 * limbs_per_col) {
            memset(e, 0, 16 * limbs_per_col);
            mp_limb_t * killer[4];
            for(int di = 0;di<4 ; di++) {
                if (i+di >= NCOLS || G->set_pivot[i+di]<0)
                    killer[di] = zcol;
                else
                    killer[di] = lmat2 + G->set_pivot[i+di] * limbs_per_col;
            }
            /* Do a gray code walk to fill in the table. */
            mp_limb_t * ex, * ey;
            int d = 0;
            memset(e, 0, limbs_per_col); ex=e;
#define DO_BIT(z) do {							\
            d^=1<<z; ey=e+d*limbs_per_col;				\
            for(int j = 0 ; j < limbs_per_col ; j++) {			\
                ey[j] = ex[j] ^ killer[z][j];				\
            }								\
            ex=ey;							\
        } while (0)
            /* zero */ DO_BIT(0); DO_BIT(1); DO_BIT(0);
            DO_BIT(2); DO_BIT(0); DO_BIT(1); DO_BIT(0);
            DO_BIT(3); DO_BIT(0); DO_BIT(1); DO_BIT(0);
            DO_BIT(2); DO_BIT(0); DO_BIT(1); DO_BIT(0);
        }
#undef  DO_BIT
        free(zcol);
    }

    free(lmat2);
    free(ptr_rows);
    free(G->set_pivot);
    free(G->set_used);

    return r;
}


/* add the pivot row to the given row, expect for the given column 
 *   mask gives the bit we have to keep
 *   j_current is the number of the limb which contains the column
 */

static inline void addRows(int row, int pivot, mp_limb_t mask,
                           /*int j_current,*/  mp_limb_t **ptr_current) {
  int i;
  mp_limb_t *ptr1, *ptr2;

  ptr1 = ptr_rows[row];
  ptr2 = ptr_rows[pivot];

  for (i = 0; i < LIMBS_PER_ROW; ++i) 
    ptr1[i] ^= ptr2[i];
  *ptr_current[row] |= mask; 
} /* end function addRows */


/* add the pivot row to the two given rows, expect for the given column 
 *   mask gives the bit we have to keep
 *   j_current is the number of the limb which contains the column
 */

static inline void add2Rows(int row, int row2, int pivot, mp_limb_t mask,
			    /*int j_current,*/ mp_limb_t **ptr_current)
{
  // int i;
  mp_limb_t *ptr1, *ptr2, *ptr3, *ptr_lim;

  ptr1 = ptr_rows[row];
  ptr2 = ptr_rows[row2];
  ptr3 = ptr_rows[pivot];
  ptr_lim = ptr_rows[pivot + 1];
  while (ptr3 != ptr_lim) {
    *ptr1++ ^= *ptr3;
    *ptr2++ ^= *ptr3++;
  }

  *ptr_current[row] |= mask; 
  *ptr_current[row2] |= mask; 
} /* end function add2Rows */


static inline void add3Rows(int row, int row2, int row3, int pivot,
			    mp_limb_t mask, /* int j_current, */
			    mp_limb_t **ptr_current) {
  // int i;
  mp_limb_t *ptr1, *ptr2, *ptr3, *ptr_piv, *ptr_lim;

  ptr1 = ptr_rows[row];
  ptr2 = ptr_rows[row2];
  ptr3 = ptr_rows[row3];
  ptr_piv = ptr_rows[pivot];
  ptr_lim = ptr_rows[pivot + 1];
  while (ptr_piv != ptr_lim) {
    mp_limb_t aux = *ptr_piv++;
    *ptr1++ ^= aux;
    *ptr2++ ^= aux;
    *ptr3++ ^= aux;
  }

  *ptr_current[row] |= mask; 
  *ptr_current[row2] |= mask; 
  *ptr_current[row3] |= mask; 
} /* end function add3Rows */

#if 0 /* unused currently */
/* add the pivot row to the given row, starting with the given column+1
   Note: the mask and ptr_current allow to do this quickly. */
/*
  TODO: do several elimination simultaneously to reduce the cost of the loop
*/

static void addPartialRows(int row, int pivot, mp_limb_t mask, int j_current,
			   mp_limb_t **ptr_current) {
  // int i;
  mp_limb_t *ptr1, *ptr2, *ptrlim;

  *ptr_current[row] ^= (*ptr_current[pivot] & mask);
  
  ptr1 = ptr_current[row] + 1;
  ptr2 = ptr_current[pivot] + 1;
  ptrlim = ptr_rows[row+1];
  while (ptr1 < ptrlim)
    *ptr1++ ^= *ptr2++;
} /* end function addPartialRows */
#endif

/* Return the index of the first row with a one in the desired column. */
/* If the column is empty, return -1 */

static inline int getPivot(mp_limb_t **ptr, mp_limb_t mask) {
  int i;
  for (i = 0; i < NROWS; ++i)
    if (*ptr[i]&mask)
      return i;
 
  return -1;
}


#if	VERBOSE
#ifndef NO_MAIN
static void printVector(mp_limb_t *V, int n) {
  int i;
  printf("[");
  for (i = 0; i < n; ++i) {
    if (V[i / MACHINE_WORD_SIZE] & (1UL << (i % MACHINE_WORD_SIZE)))
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
#endif
