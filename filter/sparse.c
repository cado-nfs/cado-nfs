/*
 * Program: history
 * Author : F. Morain
 * Purpose: managing history of merges
 *
 * Algorithm:
 *
 */


#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "portability.h"
#include "macros.h"
#include "filter_config.h"
#include "utils.h"
#include "merge_replay_matrix.h"
#include "sparse.h"

#define DEBUG 0

void
fprintRow(FILE *file, typerow_t *row)
{
  index_t i;
#ifdef FOR_DL
    fprintf(file, "[%" PRid "]", row[0].id);
    for(i = 1; i <= row[0].id; i++)
        fprintf(file, " %" PRid "(%" PRId32 ")", row[i].id, row[i].e);
#else
    fprintf(file, "[%" PRid "]", rowCell (row, 0));
    for(i = 1; i <= rowCell (row, 0); i++)
        fprintf(file, " %" PRid "", rowCell (row, i));
#endif
}

// i1 += i2
// A row is row[0..max] where row[0] = max and the real components are
// row[1..max].
// If len != -1, then it is the real length of row[i1]+row[i2].
//
// If j is given, it is the index of the column that is used for
// pivoting in the case of DL. Then, the operation is
//   i1 = e2*i1 + e1*i2
// where e1 and e2 are adjusted so that the j-th column is zero in i1.
//
// Also update the data for the index file, if needed (i.e. if the given
// pointer is not NULL).
void
addRowsUpdateIndex(typerow_t **rows, index_data_t index_data,
                   index_t i1, index_t i2, MAYBE_UNUSED index_t j)
{
    uint32_t k1, k2, k, len;
    typerow_t *tmp;

    ASSERT(i1 != i2);

    ASSERT(rows[i1] != NULL);
    ASSERT(rows[i2] != NULL);
#if DEBUG >= 1
    fprintf(stderr, "R[%d] =", i1); fprintRow(stderr, rows[i1]);
    fprintf(stderr, "\n");
    fprintf(stderr, "R[%d] =", i2); fprintRow(stderr, rows[i2]);
    fprintf(stderr, "\n");
#endif
    len = 1 + rowLength(rows, i1) + rowLength(rows, i2);
    tmp = (typerow_t *) malloc (len * sizeof(typerow_t));

#ifdef FOR_DL /* look for the exponents of j in i1 i2*/
    int e1 = 0, e2 = 0;
    int d;
    unsigned int l;
    for (l = 1 ; l <= rowLength(rows, i1) ; l++)
        if (rowCell(rows[i1], l) == j)
            e1 = rows[i1][l].e;
    for (l = 1 ; l <= rowLength(rows, i2) ; l++)
        if (rowCell(rows[i2], l) == j)
            e2 = rows[i2][l].e;

    ASSERT (e1 != 0 && e2 != 0);

    d  = (int) gcd_int64 ((int64_t) e1, (int64_t) e2);
    e1 /= -d;
    e2 /= d;

    for (l = 1 ; l <= rowLength(rows, i2) ; l++)
        rows[i2][l].e *= e1;
    for (l = 1 ; l <= rowLength(rows, i1) ; l++)
        rows[i1][l].e *= e2;

#if DEBUG >= 1
    fprintf (stdout, "Computing %d*rows[%d] + %d*rows[%d] for j=%d\n",
                      e2, i1, e1, i2, j);
#endif

#endif


    // loop while everybody is here
    k = k1 = k2 = 1;
    while((k1 <= rowLength(rows, i1)) && (k2 <= rowLength(rows, i2))){
        if (rowCell(rows[i1], k1) < rowCell(rows[i2], k2))
          {
            tmp[k++] = rowFullCell(rows[i1], k1);
            k1 ++;
          }
        else if (rowCell(rows[i1], k1) > rowCell(rows[i2], k2))
          {
            tmp[k++] = rowFullCell(rows[i2], k2);
            k2 ++;
          }
        else{
#ifdef FOR_DL
            if (rows[i1][k1].e + rows[i2][k2].e != 0)
            {
                tmp[k].id = rows[i1][k1].id;
                tmp[k++].e = rows[i1][k1].e + rows[i2][k2].e;
            }
#else
#if DEBUG >= 1
            fprintf(stderr, "WARNING: j1=j2=%d in addRows\n", k1);
#endif
#endif
            k1++;
            k2++;
        }
    }
    // finish with k1
    for( ; k1 <= rowLength(rows, i1); k1++)
        tmp[k++] = rowFullCell(rows[i1], k1);
    // finish with k2
    for( ; k2 <= rowLength(rows, i2); k2++)
        tmp[k++] = rowFullCell(rows[i2], k2);
    ASSERT(k <= len);

    setRawCell(tmp, 0, k-1, 1);

    // copy back
#if defined(FOR_DL) || __SIZEOF_INDEX__ == 4 || __SIZEOF_INDEX__ == 8
    free(rows[i1]);
    if(k == len)
      rows[i1] = tmp;
    else
      rows[i1] = reallocRow (tmp, k);
#else
    if (k-1 != rowLength(rows, i1))
      rows[i1] = reallocRow (rows[i1], k);
    compressRow (rows[i1], tmp, k-1);
    free (tmp);
#endif
#ifdef FOR_DL
    /* restore old coeff for row i2 */
    for (l = 1 ; l <= rowLength(rows, i2) ; l++)
        rows[i2][l].e /= e1;
#endif


    // Now, deal with the index_data.
    if (index_data != NULL) {
        k = k1 = k2 = 0;   // in index_data_t, we count from 0...

        relset_t r1 = index_data[i1];
        relset_t r2 = index_data[i2];
        relset_t tmp;
        tmp.n = 0;
        tmp.rels = (multirel_t *) malloc((r1.n+r2.n)*sizeof(multirel_t));
        while ((k1 < r1.n) && (k2 < r2.n)) {
            if (r1.rels[k1].ind_row < r2.rels[k2].ind_row) {
#ifdef FOR_DL
                tmp.rels[k].e = e2*r1.rels[k1].e;
#endif
                tmp.rels[k++].ind_row = r1.rels[k1++].ind_row;
            } else if (r1.rels[k1].ind_row > r2.rels[k2].ind_row) {
#ifdef FOR_DL
                tmp.rels[k].e = e1*r2.rels[k2].e;
#endif
                tmp.rels[k++].ind_row = r2.rels[k2++].ind_row;
            } else {
#ifdef FOR_DL
                int32_t e = e2*r1.rels[k1].e + e1*r2.rels[k2].e;
                if (e != 0) {
                    tmp.rels[k].ind_row = r1.rels[k1].ind_row;
                    tmp.rels[k++].e = e;
                }
#endif
            k1++;
            k2++;
            }
        }
        // finish with k1 and k2
        for( ; k1 < r1.n; ) {
#ifdef FOR_DL
            tmp.rels[k].e = e2*r1.rels[k1].e;
#endif
            tmp.rels[k++].ind_row = r1.rels[k1++].ind_row;
        }
        for( ; k2 < r2.n; ) {
#ifdef FOR_DL
            tmp.rels[k].e = e1*r2.rels[k2].e;
#endif
            tmp.rels[k++].ind_row = r2.rels[k2++].ind_row;
        }
        ASSERT (k <= r1.n + r2.n);
        tmp.n = k;
	if (k < r1.n + r2.n) /* reallocate to save memory */
	  tmp.rels = realloc (tmp.rels, k * sizeof(multirel_t));

        // copy back to i1
        free (index_data[i1].rels);
        index_data[i1] = tmp;
    }

#if DEBUG >= 1
    fprintf(stderr, "row[%d]+row[%d] =", i1, i2);
    fprintRow(stderr, rows[i1]); fprintf(stderr, "\n");
#endif
}

// A line is "[-]i i1 ... ik [#j]"
int parse_hisfile_line (index_signed_t *ind, const char *t, index_t *j)
{
  int ni = 0;
  index_signed_t sg = 1;

  if(*t == '-')
    {
      sg = -1;
	    t++;
    }

  ind[0] = 0;
  while(1)
    {
      if((*t == '\n') || (*t == '#'))
          break;
      else if (*t == ' ')
        {
          ni++;
          ind[ni] = 0;
        }
	    else
          ind[ni] = 10 * ind[ni] + (*t - '0');

	    t++;
    }

  ind[0] = sg * ind[0];
  *j=0;

  if (*t == '#')
      for (t++ ; (*t != '\n') ; t++)
          *j = 10 * (*j) + (*t - '0');
  else
    {
      ni++;
      *j = -1;
    }

  return ni;
}

/* allocates memory for a row of n elements (including the length which is
   element 0) */
typerow_t*
mallocRow (uint32_t n)
{
  typerow_t *row;

#if defined(FOR_DL) || __SIZEOF_INDEX__ == 4 || __SIZEOF_INDEX__ == 8
  row = (typerow_t*) malloc (n * sizeof (typerow_t));
#else
  /* for 5 <= __SIZEOF_INDEX__ <= 7, we need n*__SIZEOF_INDEX__ bytes */
  row = (typerow_t*) malloc (n * __SIZEOF_INDEX__);
#endif
  FATAL_ERROR_CHECK(row == NULL, "Cannot allocate memory");
  return row;
}

/* reallocates memory for a row of n elements */
typerow_t*
reallocRow (typerow_t *row, uint32_t n)
{
#if defined(FOR_DL) || __SIZEOF_INDEX__ == 4 || __SIZEOF_INDEX__ == 8
  row = (typerow_t*) realloc (row, n * sizeof (typerow_t));
#else
  /* for 5 <= __SIZEOF_INDEX__ <= 7, we need n*__SIZEOF_INDEX__ bytes */
  row = (typerow_t*) realloc (row, n * __SIZEOF_INDEX__);
#endif
  FATAL_ERROR_CHECK(row == NULL, "Cannot allocate memory");
  if (n * __SIZEOF_INDEX__ == 50) printf ("row=%lx\n", (unsigned long) row);
  return row;
}


/* experimental code for factorization with 5 <= __SIZEOF_INDEX__ <= 7 */
#if !defined(FOR_DL) && 5 <= __SIZEOF_INDEX__ && __SIZEOF_INDEX__ <= 7

/* We use the following data structure for a matrix row:

 * each entry has __SIZEOF_INDEX__ bytes which are stored consecutively
   (little endian, i.e., smallest byte first)
 * the first entry represents the length (say n) of the row
 * thus we have in total (n+1)*__SIZEOF_INDEX__ bytes */

index_t
rowCell (index_t *row, int k)
{
  index_t res = 0;
  uint8_t *ptr = ((uint8_t*) row) + k * __SIZEOF_INDEX__;

  for (int t = 0; t < __SIZEOF_INDEX__; t++)
    res += ((index_t) ptr[t]) << (t * 8);
  return res;
}

/* this function is only used for factorization, where the exponent e is not
   used */
void setCell(index_t *row, int k, index_t j, exponent_t e MAYBE_UNUSED)
{
  uint8_t *ptr = ((uint8_t*) row) + k * __SIZEOF_INDEX__;
  index_t j0 = j;

  for (int t = 0; t < __SIZEOF_INDEX__; t++)
    {
      ptr[t] = (uint8_t) j;
      j >>= 8;
    }
  /* we should have used all bits of j */
  if (j != 0) printf ("j0=%lu\n", j0);
  ASSERT_ALWAYS(j == 0);
}

/* set (compressed) row[0..n] from buf[0..n] */
void compressRow (index_t *row, index_t *buf, int n)
{
  for (int k = 0; k <= n; k++)
    setCell (row, k, buf[k], 1);
}

#endif

