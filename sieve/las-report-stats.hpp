#ifndef LAS_REPORT_STATS_HPP_
#define LAS_REPORT_STATS_HPP_

#include <string>
#include "tdict.hpp"

/* las_report: Structure for gathering reports and stats on sieving */

struct las_report_s {
    struct {
        unsigned long before_sieve;
        unsigned long after_sieve;
        unsigned long not_both_even;
        unsigned long not_both_multiples_of_p;
        unsigned long trial_divided_on_side[2];
        unsigned long check_leftover_norm_on_side[2];
        unsigned long enter_cofactoring;
        unsigned long cofactored;
        unsigned long smooth;
    } survivors;
    unsigned long reports;
    double tn[2];            /* norms */
    double ttbuckets_fill;
    double ttbuckets_apply;
    double ttf;                 /* factor_survivors */
    double ttcof;                /* cofactorisation */
    unsigned long (*survivor_sizes)[256]; /* First index: rational side */
    unsigned long (*report_sizes)[256];
};
typedef struct las_report_s las_report[1];
typedef struct las_report_s * las_report_ptr;
typedef const struct las_report_s * las_report_srcptr;


void las_report_init(las_report_ptr p);
void las_report_clear(las_report_ptr p);
void las_report_accumulate_and_clear(las_report_ptr p, las_report_ptr q);
void las_report_copy(las_report_ptr p, las_report_ptr q);

struct coarse_las_timers {
    static int bookkeeping() { return 0; }
    static int search_survivors() { return 1; }
    static int sieving(int side) { return 2 + side; }
    static int sieving_mixed() { return 4; }
    static int cofactoring(int side) { return 5 + side; }
    static int cofactoring_mixed() { return 7; }
    static int batch(int side) { return 8 + side; }
    static int batch_mixed() { return 10; }
    static int thread_wait() { return 11; }
    static std::string explain(int x) {
        switch(x) {
            case -1: return "uncounted";
            case 0: return "bookkeeping";
            case 1: return "search_survivors";
            case 2: return "sieving on side 0";
            case 3: return "sieving on side 1";
            case 4: return "sieving (not differentiated)";
            case 5: return "cofactoring on side 0";
            case 6: return "cofactoring on side 1";
            case 7: return "cofactoring (not differentiated)";
            case 8: return "product trees on side 0";
            case 9: return "product trees on side 1";
            case 10: return "product trees (not differentiated)";
            case 11: return "worker thread wait";
            default: ASSERT_ALWAYS(0);
        }
    }
};

#define TIMER_CATEGORY(timer, cat) \
    timer.current->coarse_flag = coarse_las_timers::cat

#endif	/* LAS_REPORT_STATS_HPP_ */
