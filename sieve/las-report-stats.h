#ifndef LAS_REPORT_STATS_H_
#define LAS_REPORT_STATS_H_

/* las_report: Structure for gathering reports and stats on sieving */

struct las_report_s {
    unsigned long reports;
  /* unsigned long survivors0; */ /* Obsolete */
    unsigned long survivors1;
    unsigned long survivors2;
    double tn[2];            /* norms */
    double ttbuckets_fill;
    double ttbuckets_apply;
    double ttf;                 /* factor_survivors */
    double ttcof;                /* cofactorisation */
    unsigned long (*survivor_sizes)[256]; /* First index: rational side */
    unsigned long (*report_sizes)[256];
    unsigned long both_even;
};
typedef struct las_report_s las_report[1];
typedef struct las_report_s * las_report_ptr;
typedef const struct las_report_s * las_report_srcptr;


void las_report_init(las_report_ptr p);
void las_report_clear(las_report_ptr p);
void las_report_accumulate(las_report_ptr p, las_report_ptr q);
void las_report_copy(las_report_ptr p, las_report_ptr q);
void las_report_stat_fprint(FILE *file, las_report_ptr p);

#endif	/* LAS_REPORT_STATS_H_ */
