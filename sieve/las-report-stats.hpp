#ifndef LAS_REPORT_STATS_HPP_
#define LAS_REPORT_STATS_HPP_

#include <cstring> // for memset, size_t
#include <exception>
#include <memory>   // for shared_ptr, allocator, make_shared, __shared...
#include <sstream>  // for basic_ostream::operator<<, operator<<, ostri...
#include <stdint.h> // for uint8_t
#include <string>   // for string, operator<<

#include "macros.h" // for ASSERT_ALWAYS

/* las_report: Structure for gathering reports and stats on sieving */

struct las_report
{
    struct survivors_t {
        unsigned long before_sieve = 0;
        unsigned long after_sieve = 0;
        unsigned long not_both_even = 0;
        unsigned long not_both_multiples_of_p = 0;
        unsigned long trial_divided_on_side[2] { 0, 0 };
        unsigned long check_leftover_norm_on_side[2] { 0, 0 };
        unsigned long enter_cofactoring = 0;
        unsigned long cofactored = 0;
        unsigned long smooth = 0;
    } survivors;
    unsigned long nr_sq_processed = 0;
    unsigned long nr_sq_discarded = 0;
    unsigned long total_logI = 0;
    unsigned long total_J = 0;
    unsigned long reports = 0;
    unsigned long duplicates = 0;   /* used with -dup option */
    unsigned long multi_print = 0;  /* the ones that were printed several
                                       times because the special-q was
                                       restarted. For convenience, we
                                       _also_ count them in [reports].  */
    int nwaste = 0;                 /* number of restarted special-q's */
    double waste = 0;               /* restarted special-q's */
    double cumulated_wait_time = 0; /* wait time in threadpool */
    las_report() = default;
    las_report(las_report const&) = delete;
    las_report(las_report&&) = default;
    las_report operator=(las_report const&) = delete;
    las_report& operator=(las_report&&) = default;
    class count_matrix
    {
        /* First index: rational side */
        unsigned long data[256][256];

        public:
        count_matrix() { memset(data, 0, sizeof(data)); }
        unsigned long* operator[](int x) { return data[x]; }
        unsigned long const* operator[](int x) const { return data[x]; }
    };
    std::shared_ptr<count_matrix> survivor_counts;
    std::shared_ptr<count_matrix> report_counts;

    void accumulate_and_clear(las_report&& q)
    {
        {
            unsigned long* ps = (unsigned long*)&survivors;
            unsigned long* qs = (unsigned long*)&q.survivors;
            for (size_t i = 0; i < sizeof(survivors) / sizeof(unsigned long);
                 i++) {
                ps[i] += qs[i];
            }
        }
        nr_sq_processed += q.nr_sq_processed;
        nr_sq_discarded += q.nr_sq_discarded;
        total_logI += q.total_logI;
        total_J += q.total_J;
        reports += q.reports;
        duplicates += q.duplicates;
        multi_print += q.multi_print;
        nwaste += q.nwaste;
        waste += q.waste;
        if (survivor_counts) {
            count_matrix* ps = survivor_counts.get();
            count_matrix* qs = q.survivor_counts.get();
            count_matrix* pr = report_counts.get();
            count_matrix* qr = q.report_counts.get();
            for (size_t S0 = 0; S0 < 256; ++S0) {
                for (size_t S1 = 0; S1 < 256; ++S1) {
                    (*ps)[S0][S1] += (*qs)[S0][S1];
                    (*pr)[S0][S1] += (*qr)[S0][S1];
                }
            }
        }
        q = las_report();
    }
    void mark_survivor(uint8_t S0, uint8_t S1)
    {
        /* This *must* be called on a *private* las_report structure, so
         * that taking the mutex is not necessary */
        if (survivor_counts) {
            (*survivor_counts.get())[S0][S1]++;
        }
    }
    void mark_report(uint8_t S0, uint8_t S1)
    {
        /* This *must* be called on a *private* las_report structure, so
         * that taking the mutex is not necessary */
        if (survivor_counts) {
            (*report_counts.get())[S0][S1]++;
        }
    }
    /* Well, at this point it's never used... */
    void allocate_count_matrices()
    {
        survivor_counts = std::make_shared<count_matrix>();
        report_counts = std::make_shared<count_matrix>();
    }
    void display_survivor_counters() const;
};

struct coarse_las_timers
{
    static int bookkeeping() { return 0; }
    static int search_survivors() { return 1; }
    static int sieving(int side) { return 2 + (side << 8); }
    static int sieving_mixed() { return 3; }
    static int norms(int side) { return 4 + (side << 8); }
    static int cofactoring(int side) { return 5 + (side << 8); }
    static int cofactoring_mixed() { return 6; }
    static int batch(int side) { return 7 + (side << 8); }
    static int batch_mixed() { return 8; }
    static int thread_wait() { return 9; }
    static std::string explain_base(int x)
    {
        int side = x >> 8;
        x &= 255;
        switch (x) {
            case 255:
                return "uncategorized (top-level bookkeeping)";
            case 0:
                return "bookkeeping (lower levels)";
            case 1:
                return "search_survivors";
            case 2:
                return "sieving on side " + std::to_string(side);
            case 3:
                return "sieving (not differentiated)";
            case 4:
                return "norms on side " + std::to_string(side);
            case 5:
                return "cofactoring on side " + std::to_string(side);
            case 6:
                return "cofactoring (not differentiated)";
            case 7:
                return "product trees on side " + std::to_string(side);
            case 8:
                return "product trees (not differentiated)";
            case 9:
                return "worker thread wait";
            default:
                ASSERT_ALWAYS(0);
        }
    }
    static std::string explain(int x)
    {
        std::ostringstream os;
        os << x << " " << explain_base(x);
        return os.str();
    }
};

#ifndef DISABLE_TIMINGS
#define TIMER_CATEGORY(timer, cat)                                             \
    timer.set_current_category(coarse_las_timers::cat)
#else
#define TIMER_CATEGORY(timer, cat) /**/
#endif

#endif /* LAS_REPORT_STATS_HPP_ */
