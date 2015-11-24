#ifndef MERGE_MONO_H_
#define MERGE_MONO_H_

#ifdef __cplusplus
extern "C" {
#endif

extern void addRowsAndUpdate(filter_matrix_t *mat, int i1, int i2, int32_t j);
extern void removeRowDefinitely(report_t *rep, filter_matrix_t *mat, int32_t i);
extern void MSTWithA(report_t *rep, filter_matrix_t *mat, int m, int32_t *ind, int32_t j, double *tMST, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX]);
extern int deleteHeavyColumns(report_t *rep, filter_matrix_t *mat);
extern int addFatherToSons(int history[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1], filter_matrix_t *mat, int m, int *ind, int32_t j,	int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX], int *father, int *height MAYBE_UNUSED, int hmax MAYBE_UNUSED, int sons[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX+1]);

extern int number_of_superfluous_rows(filter_matrix_t *mat);
extern void mergeOneByOne(report_t *rep, filter_matrix_t *mat, int maxlevel, double target_density);

extern void resume(report_t *rep, filter_matrix_t *mat, const char *resumename);

#ifdef __cplusplus
}
#endif

#endif	/* MERGE_MONO_H_ */
