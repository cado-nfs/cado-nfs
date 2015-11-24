#ifndef MERGE_MPI_H_
#define MERGE_MPI_H_

#ifdef __cplusplus
extern "C" {
#endif

extern void mpi_start_proc(char *outname, filter_matrix_t *mat, FILE *purgedfile, char *purgedname, double target_density, char *resumename);
extern void mpi_send_inactive_rows(int i);
extern void mpi_add_rows(filter_matrix_t *mat, int m, int32_t j, int32_t *ind);
extern void mpi_load_rows_for_j(filter_matrix_t *mat, int m, int32_t j);

#ifdef __cplusplus
}
#endif

#endif	/* MERGE_MPI_H_ */
