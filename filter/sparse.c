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
#include "macros.h"
#include "sparse.h"

#define DEBUG 0

void
fprintRow(FILE *file, int32_t *row)
{
    int i;

    fprintf(file, "[%d]", row[0]);
    for(i = 1; i <= row[0]; i++)
	fprintf(file, " %d", row[i]);
}

// row[0..row[0]] is of lenfth row[0]+1.
int32_t*
copyRow(int32_t *row)
{
    int32_t *tmp = (int32_t *)malloc((1+row[0]) * sizeof(int32_t));
    
    memcpy(tmp, row, (1+row[0]) * sizeof(int32_t));
    return tmp;
}

// i1 += i2
// A row is row[0..max] where row[0] = max and the real components are
// row[1..max].
// If len != -1, then it is the real length of row[i1]+row[i2].
void
addRows(int32_t **rows, int i1, int i2, int len0)
{
    int32_t k1, k2, k, len, *tmp;

    ASSERT(rows[i1] != NULL);
    ASSERT(rows[i2] != NULL);
#if DEBUG >= 1
    fprintf(stderr, "R[%d] =", i1); fprintRow(stderr, rows[i1]); 
    fprintf(stderr, "\n");
    fprintf(stderr, "R[%d] =", i2); fprintRow(stderr, rows[i2]);
    fprintf(stderr, "\n");
#endif
    len = 1 + (len0 != -1 ? len0 : rows[i1][0] + rows[i2][0]);
    tmp = (int32_t *)malloc(len * sizeof(int32_t));
    k = k1 = k2 = 1;

    // loop while everybody is here
    while((k1 <= rows[i1][0]) && (k2 <= rows[i2][0])){
	if(rows[i1][k1] < rows[i2][k2])
	    tmp[k++] = rows[i1][k1++];
	else if(rows[i1][k1] > rows[i2][k2])
            tmp[k++] = rows[i2][k2++];
	else{
#if DEBUG >= 1
	    fprintf(stderr, "WARNING: j1=j2=%d in addRows\n", k1);
#endif
	    k1++;
	    k2++;
	}
    }
    // finish with k1
    for( ; k1 <= rows[i1][0]; k1++)
	tmp[k++] = rows[i1][k1];
    // finish with k2
    for( ; k2 <= rows[i2][0]; k2++)
	tmp[k++] = rows[i2][k2];
    ASSERT(k <= len);
    // copy back
    free(rows[i1]);
    if(len0 != -1){
	tmp[0] = k-1;
	rows[i1] = tmp;
	ASSERT(tmp[0] == len0);
    }
    else{
        /* FIXME: why not use realloc here instead? Since k <= len,
           it suffices to shrink the tmp[] array to k entries.
           Also, we might detect the special case k = len. */
#if 0
	int *tmp2 = (int32_t *)malloc(k * sizeof(int32_t));
	memcpy(tmp2, tmp, k * sizeof(int32_t));
	tmp2[0] = k-1;
	rows[i1] = tmp2;
	free(tmp);
#else
	tmp[0] = k-1;
	if(k == len)
	    rows[i1] = tmp;
	else
	    rows[i1] = realloc(tmp, k * sizeof(int32_t));
#endif
    }
#if DEBUG >= 1
    fprintf(stderr, "row[%d]+row[%d] =", i1, i2);
    fprintRow(stderr, rows[i1]); fprintf(stderr, "\n");
#endif
}

void
removeWeight(int32_t **rows, int *wt, int i)
{
    int32_t k;

    for(k = 1; k <= rows[i][0]; k++)
	wt[rows[i][k]]--;
}

int
hasCol(int32_t **rows, int i, int32_t j)
{
    int32_t k;

    for(k = 1; k <= rows[i][0]; k++)
	if(rows[i][k] == j)
	    return 1;
    return 0;
}

int 
cmp (const void *p, const void *q)
{
  int x = *((int *)p);
  int y = *((int *)q);
  return (x <= y ? -1 : 1);
}
