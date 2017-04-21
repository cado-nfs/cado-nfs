#ifndef MST_H_
#define MST_H_

#ifdef __cplusplus
extern "C" {
#endif

// #define TIMINGS

#ifdef TIMINGS
extern double trecomputeR;
extern double tfill[MERGE_LEVEL_MAX], tmst[MERGE_LEVEL_MAX];
extern double nfill[MERGE_LEVEL_MAX];
#endif

void fillRowAddMatrix(int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], filter_matrix_t *mat, int m, index_t *ind, index_t ideal);
int minimalSpanningTree(int *start, int *end, int m, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX]);
int minCostUsingMST(filter_matrix_t *mat, int m, index_t *ind, index_t j);
void printMST(int *father, int *sons, int m, index_t *ind);


#ifdef __cplusplus
}
#endif

#endif	/* MST_H_ */
