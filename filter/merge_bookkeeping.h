#ifndef MERGE_BOOKKEEPING_H_
#define MERGE_BOOKKEEPING_H_



/* a lot of verbosity */
// #define BIG_BROTHER

/* some more verbosity which requires additional operations */
// #define BIG_BROTHER_EXPENSIVE

#ifdef BIG_BROTHER
extern unsigned char *touched_columns;
#endif

/* define CANCEL to count column cancellations */
// #define CANCEL
#ifdef CANCEL
#define CANCEL_MAX 2048
extern unsigned long cancel_rows;
extern unsigned long cancel_cols[CANCEL_MAX];
#endif

/* define DEBUG if printRow or copy_matrix is needed */
// #define DEBUG

/* cbound_incr is the increment on the maximal cost of merges at each step.
   Setting it to 1 is optimal in terms of matrix size, but will take a very
   long time (typically 10 times more than with cbound_incr=10). */
#define CBOUND_INCR_DEFAULT 8


/* Note about variables used in the code:
 * cwmax is the (current) maximal weight of columns that will be considered
   for a merge. It starts at cwmax=2. Once we have performed *all* 2-merges,
   we increase cwmax to 3, and at each step of the algorithm, we increase it
   by 1 (not waiting for all 3-merges to be completed).
 * cbound is the maximum (current) fill-in that is allowed for a merge
   (in fact, it is a biased value to avoid negative values, one should subtract
    BIAS from cbound to get the actual value). It starts at 0, and once all
    the 2-merges have been performed (which all give a negative fill-in, thus
    they will all be allowed), we increase cbound by CBOUND_INCR at each step
    of the algorithm (where CBOUND_INCR differs for integer factorization and
    discrete logarithm).
 * j0 means that we assume that columns of index < j0 cannot have
   weight <= cwmax. It depends on cwmax (decreases when cwmax increases).
   At the first call to compute_weights(), the values j0(cwmax=2) up to
   j0(MERGE_LEVEL_MAX) are computed once for all (since the weight of a
   column usually does not decrease, the values of j0 should remain correct
   during the algorithm, but not optimal).
   In several places we use the fact that the rows are sorted by increasing
   columns: if we start from the end, we can stop at soon as j < j0.
*/

#define COMPUTE_W  0 /* compute_weights */
#define COMPUTE_R  1 /* compute_R */
#define COMPUTE_M  2 /* compute_merges */
#define APPLY_M    3 /* apply_merges */
#define PASS       4 /* pass */
#define RECOMPRESS 5 /* recompress */
#define FLUSH      6 /* buffer_flush */
#define GC         7 /* garbage collection */
extern double cpu_t[8];
extern double wct_t[8];

extern int merge_verbose; /* verbosity level */

#define MARGIN 5 /* reallocate dynamic lists with increment 1/MARGIN */

/* define TRACE_J to trace all occurrences of ideal TRACE_J in the matrix */
// #define TRACE_J 1438672



#ifdef __cplusplus
extern "C" {
#endif

void
print_timings (char *s, double cpu, double wct);

#ifdef __cplusplus
}
#endif

#endif	/* MERGE_BOOKKEEPING_H_ */
