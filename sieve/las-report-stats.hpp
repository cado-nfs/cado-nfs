#ifndef LAS_REPORT_STATS_HPP_
#define LAS_REPORT_STATS_HPP_

#include <string>
#include <array>
#include <memory>
#include <string.h>
#include "tdict.hpp"
#include "las-base.hpp"

/* las_report: Structure for gathering reports and stats on sieving */

struct las_report {
    struct survivors_t : public _padded_pod<survivors_t> {
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
    unsigned long reports=0;
    unsigned long duplicates=0;   /* used with -dup option */
    double tn[2]={0,0};           /* norms */
    double ttbuckets_fill=0;
    double ttbuckets_apply=0;
    double ttf=0;                 /* factor_survivors */
    double ttcof=0;               /* cofactorisation */
    class count_matrix {
        /* First index: rational side */
        unsigned long data[256][256];
        public:
        count_matrix() { memset(data, 0, sizeof(data)); }
        unsigned long * operator[](int x) { return data[x]; }
        unsigned long const * operator[](int x) const { return data[x]; }
    };
    std::shared_ptr<count_matrix> survivor_counts;
    std::shared_ptr<count_matrix> report_counts;
    void accumulate_and_clear(las_report && q)
    {
        {
            unsigned long * ps = (unsigned long*) & survivors;
            unsigned long * qs = (unsigned long*) & q.survivors;
            for(size_t i = 0 ; i < sizeof(survivors) / sizeof(unsigned long) ; i++) {
                ps[i] += qs[i];
            }
        }
        reports += q.reports;
        duplicates += q.duplicates;
        for(int side = 0 ; side < 2 ; side++) tn[side]  += q.tn[side];
        ttbuckets_fill  += q.ttbuckets_fill;
        ttbuckets_apply += q.ttbuckets_apply;
        ttf     += q.ttf;
        ttcof     += q.ttcof;
        if (survivor_counts) {
            count_matrix * ps = survivor_counts.get();
            count_matrix * qs = q.survivor_counts.get();
            count_matrix * pr = report_counts.get();
            count_matrix * qr = q.report_counts.get();
            for(size_t S0 = 0 ; S0 < 256 ; ++S0) {
                for(size_t S1 = 0 ; S1 < 256 ; ++S1) {
                    (*ps)[S0][S1] += (*qs)[S0][S1];
                    (*pr)[S0][S1] += (*qr)[S0][S1];
                }
            }
        }
        q = las_report();
    }
    void mark_survivor(uint8_t S0, uint8_t S1) {
        if (survivor_counts) {
            (*survivor_counts.get())[S0][S1]++;
        }
    }
    void mark_report(uint8_t S0, uint8_t S1) {
        if (survivor_counts) {
            (*report_counts.get())[S0][S1]++;
        }
    }
    /* Well, at this point it's never used... */
    void allocate_count_matrices() {
        survivor_counts = std::make_shared<count_matrix>();
        report_counts = std::make_shared<count_matrix>();
    }
    void display_survivor_counters() const;
};

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
