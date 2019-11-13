void fillRowAddMatrix(int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], filter_matrix_t *mat, int m, index_t *ind, index_t ideal);
int minimalSpanningTree(int *start, int *end, int m, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX]);
int minCostUsingMST(filter_matrix_t *mat, int m, index_t *ind, index_t j);
