#include "cado.h" // IWYU pragma: keep

#ifdef LINGEN_BINARY
#error "lingen_tune_cutoffs does not work with binary, at least for the moment"
#endif

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <cstdint>
#ifdef  HAVE_SIGHUP
#include <csignal>
#endif

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <gmp.h>

#include "arith-hard.hpp"
#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "lingen_bw_dimensions.hpp"
#include "lingen_fft_select.hpp"
#include "lingen_matpoly_ft.hpp"
#include "lingen_matpoly_select.hpp"
#ifndef LINGEN_BINARY
#include "lingen_polymat.hpp"
#endif
#include "lingen_tune_cutoffs.hpp"
#include "macros.h"
#include "params.h"
#include "runtime_numeric_cast.hpp"
#include "select_mpi.h"
#include "timing.h"
#include "tree_stats.hpp"
#include "utils_cxx.hpp"

using namespace std;

void lingen_tune_cutoffs_decl_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "B",        "minimum end of bench span window");
    param_list_decl_usage(pl, "catchsig", "enable intercepting ^C");
}

void lingen_tune_cutoffs_lookup_parameters(cxx_param_list & pl)
{
    param_list_lookup_string(pl, "B");
    param_list_lookup_string(pl, "catchsig");
}

/* interface to C programs for list of cutoffs we compute *//*{{{*/
/* This is the simple C-type cutoff table which is used down the line by
 * C programs */

struct cutoff_list_s {
    size_t k;     /* this is the number of coefficients (length) of
                   * pi.  Admittedly,  it is slightly bogus to have
                   * this as a main variable.  Length [k] for pi
                   * corresponds to an "input length" of [k*(m+n)/m].
                   * In the MP case, this means that E has length
                   * [input_length+k-1], and the output E_right has
                   * length [input_length]. */
    int choice;   /* semantics vary. For FFT-based stuff,
                   * this is the depth reduction from the
                   * flint default.
                   */
};

using cutoff_list = cutoff_list_s *;

/* cutoff list *ALWAYS* finish with UINT_MAX, UINT_MAX */
static int cutoff_list_get(cutoff_list cl, size_t k)
{
    if (!cl) return 0;
    if (k < cl->k) return 0;
    for( ; k >= cl->k ; cl++);
    return (--cl)->choice;
}
/*}}}*/

/* The -B argument may be used to request printing of timing results at
 * least up to this input length */
static size_t bench_atleast_uptothis = 0;

/* timer backends *//*{{{*/
#ifdef  HAVE_GCC_STYLE_AMD64_INLINE_ASM
struct timer_rdtsc {
    double operator()() const {
        uint64_t const c = cputicks();
        return c / 1.0e9;
    }
    static const char * timer_unit() { return "Gccyc"; }
};
#endif


struct timer_rusage {
    double operator()() const { return seconds(); }
    static const char * timer_unit() { return "seconds (rusage)"; }
};
struct timer_wct {
    double operator()() const { return wct_seconds(); }
    static const char * timer_unit() { return "seconds (wct)"; }
};
/* It is important, when doing MPI benches, that we use this algorithm */
struct timer_wct_synchronized {
    double operator()() const { 
        double d = wct_seconds();
        MPI_Allreduce(MPI_IN_PLACE, &d, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
        return d;
    }
};
/*}}}*/

/* how do we do one single measurements ? We'll do it many times, but how
 * long ? Here is the place for making the corresponding choices.
 */
struct measurement_choice {
    /* provide some defaults */
    int enough_repeats = 1000;
    double minimum_time = 0.5;
    double enough_time = 1.0;
};

template<typename T>
struct small_bench {
    measurement_choice mc;
    double& tfinal;
    double last_starting_tt;
    /* accumulated_time so far. Can be increased either with
     * - add_since_last() (which resets our timer)
     * - set_since_last() (which does not)
     * - inject() (which does not either, since we do not care about our
     *   own timer in this case).
     */
    double tt = 0;
    int n = 0;

    explicit small_bench(double& t, measurement_choice mc = measurement_choice())
        : mc(mc)
        , tfinal(t)
    {
        reset();
    }
    /* are we done ? when we are, we set tfinal to the average value. */
    int done() {
        if (tt < mc.enough_time) {
            if (tt < mc.minimum_time && n < mc.enough_repeats)
                return 0;
        }
        // fprintf(stderr, "%.2f %d\n", tt, n);
        tfinal = tt / n;
        return 1;
    }
    small_bench& operator++() { n++; return *this; }
    void reset() { last_starting_tt = T()(); }
    /* adds to tt the time since the last timer reset, and resets the
     * timer.  */
    void add_since_last(unsigned int weight = 1) { 
        double now = T()();
        tt += weight * (now - last_starting_tt);
        last_starting_tt = now;
    }
    /* *sets* tt to the time since the last timer reset. Do not reset the
     * timer.  */
    void set_since_last(unsigned int weight = 1) { 
        double const now = T()();
        tt = weight * (now - last_starting_tt);
    }
    void inject(double t, unsigned int weight = 1) { tt += weight*t; }
};

/* {{{ cutoff_finder and its child small_bench */
template<typename Timer_backend = timer_rusage>
struct cutoff_finder {
    measurement_choice mc;
    double scale = 1.1;
    int stable_cutoff_break = 4;
    unsigned int ntests;
    int stable = 0;
    vector<double> benches;
    vector<bool> meaningful;
    vector<pair<size_t, pair<vector<double>, int> > > all_results;
    map<int, string> method_names;
    double slowness_ratio = 2;
    unsigned int age_slow_discard = 5;

    explicit cutoff_finder(unsigned int ntests, measurement_choice mc = measurement_choice())
        : mc(mc)
        , ntests(ntests)
        , benches(ntests)
        , meaningful(ntests , true)
    {
    }

    void set_method_name(int k, string const& s) {
        method_names.insert(make_pair(k, s));
    }
    string method_name(int i) const {
        ostringstream v;
        const auto z = method_names.find(i);
        if (z == method_names.end()) {
            v << "method " << i;
            return v.str();
        }
        return z->second;
    }

    size_t next_length(size_t k) const {
        return std::max(k+1, size_t(scale*double(k)));
    }

    int done() const {
        int nactive=0;
        for(unsigned int i = 0 ; i < ntests ; i++) {
            nactive += meaningful[i];
        }
        if (nactive > 1) return 0;
        return !(stable < stable_cutoff_break);
    }

    bool still_meaningful_to_test(int i) { return meaningful[i]; }

    string summarize_for_this_length(size_t k)
    {
        std::ostringstream comments;
        int best = -1;
        for(unsigned int i = 0 ; i < ntests ; i++) {
            if (meaningful[i] && (best < 0 || benches[i] < benches[best])) best = i;
        }
        ASSERT_ALWAYS(best >= 0);
        all_results.emplace_back(k,
                    make_pair(benches, best));

        if (all_results.empty() || best != all_results.back().second.second) {
            stable = 0;
        } else {
            stable++;
        }

        vector<bool> was_meaningful = meaningful;

        /* try to discard those which are slow. We live on the assumption
         * that the higher-numbered strategies win eventually, therefore
         * we don't discard them if [best] is smaller currently */
        for(int i = 0 ; i < best ; i++) {
            /* already dead ? */
            if (!meaningful[i]) continue;

            /* for how long has it been slow ? */
            unsigned int age_slow = 0;
            for( ; age_slow < all_results.size() ; age_slow++) {
                vector<double> const& these(all_results[all_results.size()-1-age_slow].second.first);
                int const best_there = all_results[all_results.size()-1-age_slow].second.second;
                /* do not discard this method if a less advanced one is
                 * still alive.
                 */
                if (best_there < i) break;
                /* not so slow */
                if (these[i] <= these[best_there] * slowness_ratio) break;
            }
            if (age_slow >= age_slow_discard) {
                /* ok, discard */
                comments << "  ; discarding " << method_name(i);
                meaningful[i] = false;
            }
        }
        ostringstream s;
        s << k;
        for(unsigned int i = 0 ; i < ntests ; i++) {
            if (was_meaningful[i]) {
                s << " " << benches[i];
                if ((int) i==best) s << '#';
            } else {
                s << " *";
            }
            benches[i] = 999999;
        }
        s << " " << method_name(best) << comments.str();
        return s.str();
    }

    small_bench<Timer_backend> micro_bench(int i) {
        return small_bench<Timer_backend>(benches[i], mc);
    }
    small_bench<Timer_backend> micro_bench(int i, measurement_choice mc1) {
        return small_bench<Timer_backend>(benches[i], mc1);
    }

    /* This is really limited to karatsuba-like cutoffs. It's ugly */
    vector<pair<size_t, int> > export_best_table()
    {
        vector<pair<size_t, int> > steps;
        steps.emplace_back(1,0);

        for(auto const & it : all_results) {
            size_t const size = it.first;
            int const best = it.second.second;
            if (steps.empty() || best != steps.back().second) {
                steps.emplace_back(size, best);
            }
        }
        return steps;
    }
    vector<pair<size_t, int> > export_kara_cutoff_data(struct polymat_cutoff_info * dst)
    {
        /* This size will eventually feed the ->subdivide field for the
         * cutoff info */
        size_t first_kara_size = SIZE_MAX;

        /* This size will eventually feed the ->cut field for the
         * cutoff info. For sizes above this, we will _always_ use
         * karatsuba. */
        size_t first_alwayskara_size = SIZE_MAX;

        dst->cut = first_alwayskara_size;
        dst->subdivide = first_kara_size;
        dst->table.clear();
        dst->table.emplace_back(1,0);

        for(auto const & it : all_results) {
            size_t const size = it.first;
            vector<double> const& benches(it.second.first);
            int best = it.second.second;
            /* In case fft wins, we invent something which will use
             * karatsuba still */
            if (best > 1) best = benches[1] < benches[0];
            if (dst->table.empty() || best != dst->table.back().second) {
                dst->table.emplace_back(size, best);
                if (best == 1 && size < first_kara_size)
                    first_kara_size = size;
                /* assign it multiple times */
                first_alwayskara_size = size;
            }
        }
        return dst->table;
    }
    vector<pair<size_t, int> > export_kara_cutoff_data_force_kara_now(struct polymat_cutoff_info * dst, size_t size)
    {
        vector<double> const allz(ntests);
        all_results.emplace_back(size, make_pair(allz, 1));
        auto x = export_kara_cutoff_data(dst);
        all_results.pop_back();
        return x;
    }
    static string print_result(vector<pair<size_t, int> > const& tab) {
        ostringstream s;
        s << "{";
        for(auto const & y : tab)
            s << " { " << y.first << ", " << y.second << " },";
        s << " }";
        return s.str();
    }

};
/* }}} */

#ifdef  HAVE_SIGHUP
static double last_hup = 0;
static int hup_caught = 0;

static void sighandler(int sig MAYBE_UNUSED)
{
    double const t = wct_seconds();
    if (t < last_hup + 0.5) {
        fprintf(stderr, "Interrupt twice in half a second; exiting\n");
        exit(1);
    }
    last_hup = wct_seconds();
    hup_caught++;
}


static void catch_control_signals()
{
    struct sigaction sa[1];
    memset(sa, 0, sizeof(sa));
    sa->sa_handler = sighandler;
    sigaction(SIGHUP, sa, nullptr);
    sigaction(SIGINT, sa, nullptr);

    /* 
     * should play sigprocmask and friends
    sigset_t sset[1];
    sigemptyset(sset);
    sigaddset(sset, SIGHUP);
    */
}
#endif

/* This code will first try to bench the basic operations (first
 * local; global are later) of middle product E*pi and product pi*pi.
 * The goal is to see where is the good value for the thresholds:
 * polymat_mul_kara_threshold
 * polymat_mp_kara_threshold
 */

// It's buggy and this must be investigated. We want to keep it compiled.
template<bool is_binary>
static void lingen_tune_mul_fti_depth [[maybe_unused]] (typename matpoly<is_binary>::arith_hard * ab, unsigned int m, unsigned int n, cutoff_list *cl_out)/*{{{*/
{
    cxx_gmp_randstate rstate;

    gmp_randseed_ui(rstate, 1);

#ifdef  HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using timer_t = timer_rdtsc;
#else
    using timer_t = timer_rusage;
#endif

    int const nadjs=7;
    measurement_choice mc;
    mc.enough_time = 2.0;
    mc.minimum_time = 0.1;
    mc.enough_repeats = 100;

    measurement_choice mcout = mc;
    mcout.enough_repeats = 2;
    mc.enough_time = 1.0;

    cutoff_finder<timer_t> finder(nadjs, mcout);
    finder.slowness_ratio = 1.4;
    // finder.age_slow_discard = 5;
    finder.scale = 1.01;

    for(int i = 0 ; i < nadjs ; i++) {
        ostringstream o;
        o << "depth_adj=" << nadjs-1-i;
        finder.set_method_name(i, o.str());
    }

    cout << "# Tuning FFT depth adjustment (mul) for"
                << " m=" << m
                << ", n=" << n
                << "\n";

    cout << "# Timings reported in " << timer_t::timer_unit() << "\n";
    cout << "# inputlength, ncoeffs, various depth adjustments";
    cout << "\n";

    /* for multiplication */
    cout << "# Note: for input length k within lingen,"
        << " we use ncoeffs=k*m/(m+n) = "<<(double)m/(m+n)<<"*k\n";

    /* This is for forcing the bench to run until a large length. This is
     * unnecessary here */
    // unsigned int const min_bench = 0;

    /* Beware, k is the length of piL, not a degree. Hence length 1
     * clearly makes no sense */
    for(size_t k = 2 ; !hup_caught && (/* k < min_bench || */ !finder.done()) ; k=finder.next_length(k)) {
        size_t const input_length = (m+n) * k / m;

        mpz_t p;
        mpz_init_set(p, ab->characteristic());

        auto * A = ab->alloc(k);
        auto * B = ab->alloc(k);
        auto * C = ab->alloc(2*k-1);
        
        ab->vec_set_random(A, k, rstate);
        ab->vec_set_random(B, k, rstate);
        ab->vec_set_zero(C, 2*k-1);

        ostringstream extra_info;

        for(int index = 0 ; index < nadjs ; index++) {
            struct fft_transform_info fti[1];
            fft_transform_info_init_fppol(fti, p,
                    runtime_numeric_cast<mp_size_t>(k),
                    runtime_numeric_cast<mp_size_t>(k),
                    m+n);
            fft_transform_info_set_first_guess(fti);
            if (index == 0) {
                extra_info << "(from " << fti->depth << ")";
            }
            if (nadjs-1-index >= fti->depth) {
                finder.benches[index] = 999999999;
                continue;
            }
            if (fti->depth >= 11 && index > 0) {
                finder.benches[index] = 999999999;
                continue;
            }
            fft_transform_info_adjust_depth(fti, nadjs-1-index);

            size_t fft_alloc_sizes[3];
            fft_transform_info_get_alloc_sizes(fti, fft_alloc_sizes);

            using temp_t = std::unique_ptr<void, free_delete<void>>;
            const std::unique_ptr<void, free_delete<void>> test_diagnostic(malloc(fft_alloc_sizes[2]));
            const temp_t tt(malloc(fft_alloc_sizes[2]));
            const temp_t qt(malloc(fft_alloc_sizes[1]));
            const temp_t tA(malloc(fft_alloc_sizes[0]));
            const temp_t tB(malloc(fft_alloc_sizes[0]));
            const temp_t tC(malloc(fft_alloc_sizes[0]));

            using bt = small_bench<timer_t>;
            for(bt x = finder.micro_bench(index); !x.done(); ++x) {

                double t_dftA=0;
                for(bt y = bt(t_dftA, mc); !y.done(); ++y) {
                    fft_prepare(fti, tA.get());
                    fft_dft(fti, tA.get(),
                            (mp_limb_t*)A,
                            runtime_numeric_cast<mp_size_t>(k),
                            qt.get());
                    y.set_since_last();
                }
                x.inject(t_dftA, (m+n)*(m+n));

                double t_dftB=0;
                for(bt y = bt(t_dftB, mc); !y.done(); ++y) {
                    fft_prepare(fti, tB.get());
                    fft_dft(fti, tB.get(),
                            (mp_limb_t*)B,
                            runtime_numeric_cast<mp_size_t>(k),
                            qt.get());
                    y.set_since_last();
                }
                x.inject(t_dftB, (m+n)*(m+n));

                double t_conv=0;
                for(bt y = bt(t_conv, mc); !y.done(); ++y) {
                    fft_prepare(fti, tC.get());
                    fft_zero(fti, tC.get());
                    fft_addcompose(fti, tC.get(), tA.get(), tB.get(),
                            tt.get(), qt.get());
                    y.set_since_last();
                }
                x.inject(t_conv, (m+n)*(m+n)*(m+n));

                double const t_iftC=0;
                for(bt y = bt(t_conv, mc); !y.done(); ++y) {
                    fft_prepare(fti, tC.get());
                    fft_ift(fti, (mp_limb_t*)C, 2*k-1, tC.get(), qt.get());
                    y.set_since_last();
                }
                x.inject(t_iftC, (m+n)*(m+n));
            }
        }

        ab->free(A);
        ab->free(B);
        ab->free(C);
        
        mpz_clear(p);

        cout << input_length
            << " " << finder.summarize_for_this_length(k)
            << extra_info.str()
            << "\n";
    }
    hup_caught = 0;

    auto table = finder.export_best_table();

    for(auto & x : table)
        x.second = nadjs-1-x.second;

    cout << "/* FFT depth adjustments for "
                << (m)<<"*"<<(m+n)
                <<" times "
                << (m+n)<<"*"<<(m+n)<<" products */\n";
    cout << "#define MUL_FTI_DEPTH_ADJ_" <<m+n<<"_"<<(m+n)<<"_"<<(m+n)
        << " " << decltype(finder)::print_result(table) << "\n";

    if (cl_out) {
        *cl_out = (cutoff_list) malloc((table.size()+1)*sizeof(**cl_out));
        for(unsigned int i = 0 ; i < table.size(); i++) {
            (*cl_out)[i].k = table[i].first;
            (*cl_out)[i].choice = table[i].second;
        }
        (*cl_out)[table.size()].k      = SIZE_MAX;
        (*cl_out)[table.size()].choice = INT_MAX;
    }
}/*}}}*/

// It's buggy and this must be investigated. We want to keep it compiled.
template<bool is_binary>
static void lingen_tune_mp_fti_depth [[maybe_unused]] (typename matpoly<is_binary>::arith_hard * ab, unsigned int m, unsigned int n, cutoff_list * cl_out) /*{{{*/
{
    cxx_gmp_randstate rstate;
    gmp_randseed_ui(rstate, 1);

#ifdef  HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using timer_t = timer_rdtsc;
#else
    using timer_t = timer_rusage;
#endif

    int const nadjs=7;
    measurement_choice mc;
    mc.enough_time = 2.0;
    mc.minimum_time = 0.1;
    mc.enough_repeats = 100;

    measurement_choice mcout = mc;
    mcout.enough_repeats = 2;
    mc.enough_time = 1.0;

    cutoff_finder<timer_t> finder(nadjs, mcout);
    finder.slowness_ratio = 1.4;
    // finder.age_slow_discard = 5;
    finder.scale = 1.01;

    for(int i = 0 ; i < nadjs ; i++) {
        ostringstream o;
        o << "depth_adj=" << nadjs-1-i;
        finder.set_method_name(i, o.str());
    }

    cout << "# Tuning FFT depth adjustment (mp) for"
                << " m=" << m
                << ", n=" << n
                << "\n";

    cout << "# Timings reported in " << timer_t::timer_unit() << "\n";
    cout << "# inputlength, ncoeffs, various depth adjustments";
    cout << "\n";

    /* for multiplication */
    cout << "# Note: for input length k within lingen,"
        << " we use ncoeffs=k*m/(m+n) = "<<(double)m/(m+n)<<"*k\n";

    /* This is for forcing the bench to run until a large length. This is
     * unnecessary here */
    // unsigned int const min_bench = 0;

    /* Beware, k is the length of piL, not a degree. Hence length 1
     * clearly makes no sense */
    for(size_t k = 2 ; !hup_caught && (/* k < min_bench || */ !finder.done()) ; k=finder.next_length(k)) {
        size_t const input_length = (m+n) * k / m;

        size_t const E_length = k + input_length - 1;

        cxx_mpz p(ab->characteristic());

        auto * A = ab->alloc(E_length);
        auto * B = ab->alloc(k);
        auto * C = ab->alloc(input_length);
        
        ab->vec_set_random(A, E_length, rstate);
        ab->vec_set_random(B, k, rstate);
        ab->vec_set_zero(C, input_length);

        ostringstream extra_info;

        for(int index = 0 ; index < nadjs ; index++) {
            struct fft_transform_info fti[1];
            fft_transform_info_init_fppol_mp(fti, p,
                    runtime_numeric_cast<mp_size_t>(k),
                    runtime_numeric_cast<mp_size_t>(E_length),
                    m+n);
            fft_transform_info_set_first_guess(fti);
            if (index == 0) {
                extra_info << "(from " << fti->depth << ")";
            }
            if (nadjs-1-index >= fti->depth) {
                finder.benches[index] = 999999999;
                continue;
            }
            if (fti->depth >= 11 && index > 0) {
                finder.benches[index] = 999999999;
                continue;
            }
            fft_transform_info_adjust_depth(fti, nadjs-1-index);

            size_t fft_alloc_sizes[3];
            fft_transform_info_get_alloc_sizes(fti, fft_alloc_sizes);
            using temp_t = std::unique_ptr<void, free_delete<void>>;
            // NOLINTBEGIN(cppcoreguidelines-no-malloc,hicpp-no-malloc)
            const temp_t tt(malloc(fft_alloc_sizes[2]));
            const temp_t qt(malloc(fft_alloc_sizes[1]));
            const temp_t tA(malloc(fft_alloc_sizes[0]));
            const temp_t tB(malloc(fft_alloc_sizes[0]));
            const temp_t tC(malloc(fft_alloc_sizes[0]));
            // NOLINTEND(cppcoreguidelines-no-malloc,hicpp-no-malloc)

            using bt = small_bench<timer_t>;
            for(bt x = finder.micro_bench(index); !x.done(); ++x) {

                double t_dftA=0;
                for(bt y = bt(t_dftA, mc); !y.done(); ++y) {
                    fft_prepare(fti, tA.get());
                    fft_dft(fti, tA.get(), (mp_limb_t*)A,
                            runtime_numeric_cast<mp_size_t>(E_length),
                            qt.get());
                    y.set_since_last();
                }
                x.inject(t_dftA, m*(m+n));

                double t_dftB=0;
                for(bt y = bt(t_dftB, mc); !y.done(); ++y) {
                    fft_prepare(fti, tB.get());
                    fft_dft(fti, tB.get(), (mp_limb_t*)B,
                            runtime_numeric_cast<mp_size_t>(k),
                            qt.get());
                    y.set_since_last();
                }
                x.inject(t_dftB, (m+n)*(m+n));

                double t_conv=0;
                for(bt y = bt(t_conv, mc); !y.done(); ++y) {
                    fft_prepare(fti, tC.get());
                    fft_zero(fti, tC.get());
                    fft_addcompose(fti, tC.get(), tA.get(), tB.get(), tt.get(), qt.get());
                    y.set_since_last();
                }
                x.inject(t_conv, m*(m+n)*(m+n));

                double const t_iftC=0;
                for(bt y = bt(t_conv, mc); !y.done(); ++y) {
                    fft_prepare(fti, tC.get());
                    fft_ift(fti, (mp_limb_t*)C, 
                            runtime_numeric_cast<mp_size_t>(input_length),
                            tC.get(), qt.get());
                    y.set_since_last();
                }
                x.inject(t_iftC, m*(m+n));
            }
        }

        ab->free(A);
        ab->free(B);
        ab->free(C);
        
        cout << input_length
            << " " << finder.summarize_for_this_length(k)
            << extra_info.str()
            << "\n";
    }
    hup_caught = 0;

    auto table = finder.export_best_table();

    for(auto & x : table)
        x.second = nadjs-1-x.second;

    cout << "/* FFT depth adjustments for "
                << (m)<<"*"<<(m+n)
                <<" times "
                << (m+n)<<"*"<<(m+n)<<" middle-products */\n";
    cout << "#define MP_FTI_DEPTH_ADJ_" <<m<<"_"<<(m+n)<<"_"<<(m+n)
        << " " << decltype(finder)::print_result(table) << "\n";

    if (cl_out) {
        *cl_out = (cutoff_list) malloc((table.size()+1)*sizeof(**cl_out));
        for(unsigned int i = 0 ; i < table.size(); i++) {
            (*cl_out)[i].k = table[i].first;
            (*cl_out)[i].choice = table[i].second;
        }
        (*cl_out)[table.size()].k = SIZE_MAX;
        (*cl_out)[table.size()].choice = INT_MAX;
    }
}/*}}}*/


template<bool is_binary>
static void lingen_tune_mul(typename matpoly<is_binary>::arith_hard * ab, unsigned int m, unsigned int n, cutoff_list cl MAYBE_UNUSED)/*{{{*/
{
    using fft_type = fft_transform_info;
    cxx_gmp_randstate rstate;

    gmp_randseed_ui(rstate, 1);
#define TUNE_MUL_FINDER_NMETHODS 4

    using timer_t = timer_rusage;
    cutoff_finder<timer_rusage> finder(TUNE_MUL_FINDER_NMETHODS);
    finder.set_method_name(0, "polymat-basecase");
    finder.set_method_name(1, "polymat-karatsuba");
    finder.set_method_name(2, "matpoly-kronecker-schönhage");
    finder.set_method_name(3, "matpoly-ft-kronecker-schönhage-caching");

    polymat_cutoff_info always_basecase[1];
    polymat_cutoff_info improved[1];

    cout << "# Tuning "
                << (m+n)<<"*"<<(m+n)
                <<" times "
                << (m+n)<<"*"<<(m+n)<<" products\n";

    cout << "# inputlength, ncoeffs";
    for(unsigned int i = 0 ; i < finder.ntests ; i++) {
        cout << ", " << finder.method_name(i);
    }
    cout << "\n";
    cout << "# Note: for input length k within lingen,"
        << " we use ncoeffs=k*m/(m+n) = "<<(double)m/(m+n)<<"*k\n";
    /* input length k means we consider the pi polynomial which is
     * created by k successive steps. So this mulitplication, in effect,
     * occurs at the end of a recursive procedure of length 2k, which
     * collects the pi matrices of its two child calls.
     */

    /* make sure we don't stop the bench absurdly early. We set the bench
     * as input_length \approx 2000 */
    size_t const min_bench = bench_atleast_uptothis * m / (m+n);

    /* Beware, k is the length of piL, not a degree. Hence length 1
     * clearly makes no sense */
    for(size_t k = 2 ; !hup_caught && (k < min_bench || !finder.done()) ; k=finder.next_length(k)) {
        size_t const input_length = (m+n) * k / m;
        polymat piL  (ab, m+n, m+n, k);
        polymat piR  (ab, m+n, m+n, k);
        polymat piref;
        piL.fill_random(k, rstate);
        piR.fill_random(k, rstate);
        
        ostringstream const extra_info;

        if (finder.still_meaningful_to_test(0)) {
            polymat pi;
            /* disable kara for a minute */
            polymat_set_mul_kara_cutoff(always_basecase, nullptr);
            for(small_bench<timer_t> x = finder.micro_bench(0); !x.done(); ++x) {
                pi.mul(piL, piR);
                x.set_since_last();
            }
            if (piref.get_size() == 0) {
                piref = std::move(pi);
                // fprintf(stderr, "BASIS0\n");
            } else if (pi.cmp(piref) != 0) {
                fprintf(stderr, "MISMATCH0!\n");
            }
        }

        if (finder.still_meaningful_to_test(1)) {
            polymat pi;
            /* This temporarily sets the cutoff table to enable karatsuba for
             * length >=k (hence for this test) at least, possibly using
             * karatsuba one or more times in the recursive calls depending
             * on what has been measured as best so far.
             */
            finder.export_kara_cutoff_data_force_kara_now(improved, k);
            polymat_set_mul_kara_cutoff(improved, nullptr);
            for(small_bench<timer_t> x = finder.micro_bench(1); !x.done(); ++x) {
                pi.mul(piL, piR);
                x.set_since_last();
            }
            if (piref.get_size() == 0) {
                piref = std::move(pi);
                // fprintf(stderr, "BASIS1\n");
            } else if (pi.cmp(piref) != 0) {
                fprintf(stderr, "MISMATCH1!\n");
            }
        }

        /* don't exaggerate our memory requirements */
        matpoly<is_binary> xpiL;
        xpiL.set_polymat(piL);
        piL = polymat();

        matpoly<is_binary> xpiR;
        xpiR.set_polymat(piR);
        piR = polymat();

        matpoly<is_binary> xpiref;
        if (piref.get_size()) {
            xpiref.set_polymat(piref);
            piref = polymat();
        }

        /* we should make the effort of converting polymat to matpoly,
         * right ? */

        if (finder.still_meaningful_to_test(2)) {
            matpoly<is_binary> xpi;
            /* The matpoly layer is just completetly different -- and gets
             * faster quite early on... */
            for(small_bench<timer_t> x = finder.micro_bench(2); !x.done(); ++x) {
                xpi = matpoly<is_binary>::mul(xpiL, xpiR);
                x.set_since_last();
            }
            if (xpiref.get_size() == 0) {
                xpiref = std::move(xpi);
                // fprintf(stderr, "BASIS2\n");
            } else if (xpi.cmp(xpiref) != 0) {
                fprintf(stderr, "MISMATCH2!\n");
            }
        }

        if (finder.still_meaningful_to_test(3)) {
            matpoly<is_binary> xpi;
            int const adj = cutoff_list_get(cl, k);
            {
                ostringstream o;
                o << "matpoly-ft-kronecker-schönhage-caching@adj" << adj;
                finder.set_method_name(3, o.str());
            }
#if 0
            matpoly_ft tpi, tpiL, tpiR;
            mpz_t p;
            mpz_init(p);
            abfield_characteristic(ab, p);
            struct fft_transform_info fti[1];
            fft_get_transform_info_fppol(fti, p, xpiL->size, xpiR->size, xpiL->n);
            int s = 0;
            matpoly_clear(ab, xpi);
            matpoly_init(ab, xpi, m+n, m+n, xpiL->size + xpiR->size - 1);
            for(small_bench<timer_t> x = finder.micro_bench(3); !x.done(); ++x) {
                matpoly_ft_init(fti, ab, tpiL, xpiL->m, xpiL->n);
                matpoly_ft_init(fti, ab, tpiR, xpiR->m, xpiR->n);
                matpoly_ft_init(fti, ab, tpi, xpiL->m, xpiR->n);
                matpoly_ft_dft(fti, ab, tpiL, xpiL);
                matpoly_ft_dft(fti, ab, tpiR, xpiR);
                matpoly_ft_mul(fti, ab, tpi, tpiL, tpiR);
                xpi->size = xpiL->size + xpiR->size - 1;
                ASSERT_ALWAYS(xpi->size <= xpi->alloc);
                matpoly_ft_ift(fti, ab, xpi, tpi);
                matpoly_ft_clear(fti, ab, tpiL);
                matpoly_ft_clear(fti, ab, tpiR);
                matpoly_ft_clear(ab, tpi,  fti);
                x.set_since_last();
                s++;
            }
            mpz_clear(p);
#else
            /* TODO: matpoly_mul_caching_adj, being the back-end of
             * matpoly_mul_caching, is likely to gain a stats argument
             * someday (not if deemed useless, though). If this happens,
             * we'll have to add the "stats" argument here only so that the
             * whole thing compiles. Alas this won't work out of the box, the
             * stats structure must receive information regarding the
             * planned substeps, prior to the call. So either we fix the
             * code here, either we loosen our API complexity to this
             * regard, or we ditch this tune-cutoffs file entirely.
             * (anyway the more important plingen_tuning code also needs
             * some update so that this change can happen).
             */
            tree_stats stats;
            for(small_bench<timer_t> x = finder.micro_bench(3); !x.done(); ++x) {
                xpi = matpoly_ft<fft_type>::mul_caching_adj(stats, xpiL, xpiR, adj, nullptr);
                x.set_since_last();
            }
#endif
            if (xpiref.get_size() == 0) {
                xpiref = std::move(xpi);
            } else if (xpi.cmp(xpiref) != 0) {
                fprintf(stderr, "MISMATCH3!\n");
            }
        }

        cout << input_length
            << " " << finder.summarize_for_this_length(k)
            << extra_info.str()
            << "\n";
    }
    hup_caught = 0;


    auto const table = finder.export_kara_cutoff_data(improved);
    polymat_set_mul_kara_cutoff(improved, nullptr);

    cout << "/* Cutoffs (0=basecase, 1=kara) for "
                << (m+n)<<"*"<<(m+n)
                <<" times "
                << (m+n)<<"*"<<(m+n)<<" products: */\n";
    cout << "#define MUL_CUTOFFS_" <<(m+n)<<"_"<<(m+n)<<"_"<<(m+n)
        << " " << finder.print_result(table) << "\n";
}/*}}}*/

template<bool is_binary>
static void lingen_tune_mp(typename matpoly<is_binary>::arith_hard * ab, unsigned int m, unsigned int n, cutoff_list cl MAYBE_UNUSED)/*{{{*/
{
    using fft_type = fft_transform_info;
    cxx_gmp_randstate rstate;

    gmp_randseed_ui(rstate, 1);
#define TUNE_MP_FINDER_NMETHODS 4

    using timer_t = timer_rusage;
    cutoff_finder<timer_rusage> finder(TUNE_MP_FINDER_NMETHODS);
    finder.set_method_name(0, "polymat-basecase");
    finder.set_method_name(1, "polymat-karatsuba");
    finder.set_method_name(2, "matpoly-kronecker-schönhage");
    finder.set_method_name(3, "matpoly-ft-kronecker-schönhage-caching");

    polymat_cutoff_info always_basecase[1];
    polymat_cutoff_info improved[1];

    cout << "# Tuning "
                << (m)<<"*"<<(m+n)
                <<" times "
                << (m+n)<<"*"<<(m+n)<<" middle-products\n";

    cout << "# inputlength, ncoeffs";
    for(unsigned int i = 0 ; i < finder.ntests ; i++) {
        cout << ", " << finder.method_name(i);
    }
    cout << "\n";

    cout << "# Note: for input length k within lingen, "
         << " we use MP((1+m/(m+n))k,m/(m+n)k)->k"
         << " = MP("    <<1+(double)m/(m+n)<<"*k"
                        <<", "
                        <<(double)m/(m+n)<<"*k"
                        << ")\n";
    /* we use the same semantics as for testing the product. input length
     * k means we want to count how it takes to do the required MP with
     * the pi matrix which has been computed after k successive steps.
     * IOW, this is what happens at the middle of a recursive procedure
     * on length 2k. In terms of degrees, this will be:
     *     2k times m/(m+n)k --> k  *BUT*...
     * we want coefficients [k..2k[ ; but then, the contribution from the
     * first (k-m/(m+n)k) coefficients is useless. Therefore, we rewrite
     * the polynomial lengths as:
     *      (1+m/(m+n))k times m/(m+n)k --> k
     */

    /* make sure we don't stop the bench absurdly early. We set the bench
     * as input_length \approx 2000 */
    size_t const min_bench = bench_atleast_uptothis * m / (m+n);

    /* Beware, k is a length of piL, not a degree. Hence length 1 clearly
     * makes no sense */
    for(size_t k = 2 ; !hup_caught && (k < min_bench || !finder.done()) ; k=finder.next_length(k)) {
        size_t const input_length = (m+n) * k / m;
        /* MP(degree a, degree b) -> degree b-a
         *    length la=a+1, length lb=b+1 -> length b-a+1 = lb-la+1
         * we do want lc=input_length, we have set la=k, whence lb =
         * k+input_length-1
         */
        size_t const E_length = k + input_length - 1;
        polymat E(ab, m, m+n, E_length);
        polymat piL(ab, m+n, m+n, k);
        polymat ERref;
        E.fill_random(E_length, rstate);
        piL.fill_random(k, rstate);

        ostringstream const extra_info;

        if (finder.still_meaningful_to_test(0)) {
            polymat ER;
            /* disable kara for a minute */
            polymat_set_mp_kara_cutoff(always_basecase, nullptr);
            for(small_bench<timer_t> x = finder.micro_bench(0); !x.done(); ++x) {
                ER.mp(E, piL);
                x.set_since_last();
            }
            if (ERref.get_size() == 0) {
                ERref = std::move(ER);
                // fprintf(stderr, "BASIS0\n");
            } else if (ER.cmp(ERref) != 0) {
                fprintf(stderr, "MISMATCH0!\n");
            }
        }

        if (finder.still_meaningful_to_test(1)) {
            polymat ER;
            /* This temporarily sets the cutoff table to enable karatsuba for
             * length >=k (hence for this test) at least, possibly using
             * karatsuba one or more times in the recursive calls depending
             * on what has been measured as best so far.
             */
            finder.export_kara_cutoff_data_force_kara_now(improved, k);
            polymat_set_mp_kara_cutoff(improved, nullptr);
            for(small_bench<timer_t> x = finder.micro_bench(1); !x.done(); ++x) {
                ER.mp(E, piL);
                x.set_since_last();
            }
            if (ERref.get_size() == 0) {
                ERref = std::move(ER);
                // fprintf(stderr, "BASIS1\n");
            } else if (ER.cmp(ERref) != 0) {
                fprintf(stderr, "MISMATCH1!\n");
            }
        }

        /* don't exaggerate our memory requirements */
        matpoly<is_binary> xpiL(ab, m+n, m+n, k);
        xpiL.set_polymat(piL);
        piL = polymat();

        matpoly<is_binary> xE(ab, m, m+n, E_length);
        xE.set_polymat(E);
        E = polymat();

        matpoly<is_binary> xERref;
        if (ERref.get_size()) {
            xERref.set_polymat(ERref);
        }


        /* we should make the effort of converting polymat to matpoly,
         * right ? */

        if (finder.still_meaningful_to_test(2)) {
            matpoly<is_binary> xER;
            /* The matpoly layer is just completetly different -- and gets
             * faster quite early on... */
            for(small_bench<timer_t> x = finder.micro_bench(2); !x.done(); ++x) {
                xER = matpoly<is_binary>::mp(xE, xpiL);
                x.set_since_last();
            }
            if (xERref.get_size() == 0) {
                xERref = std::move(xER);
                // fprintf(stderr, "BASIS2\n");
            } else if (xER.cmp(xERref) != 0) {
                fprintf(stderr, "MISMATCH2!\n");
            }
        }

        if (finder.still_meaningful_to_test(3)) {
            matpoly<is_binary> xER;
            int adj = UINT_MAX;
            if (cl) {
                adj = cutoff_list_get(cl, k);
                ostringstream o;
                o << "matpoly-ft-kronecker-schönhage-caching@adj" << adj;
                finder.set_method_name(3, o.str());
            }
#if 0
            matpoly_ft tER, tpiL, tE;
            mpz_t p;
            mpz_init(p);
            abfield_characteristic(ab, p);
            struct fft_transform_info fti[1];
            fft_get_transform_info_fppol_mp(fti, p, xpiL->size, xE->size, xpiL->m);
            int s = 0;
            matpoly_clear(ab, xER);
            matpoly_init(ab, xER, m, m+n, input_length);
            for(small_bench<timer_t> x = finder.micro_bench(3); !x.done(); ++x) {
                matpoly_ft_init(fti, ab, tpiL, xpiL->m, xpiL->n);
                matpoly_ft_init(fti, ab, tE, xE->m, xE->n);
                matpoly_ft_init(fti, ab, tER, xE->m, xpiL->n);
                matpoly_ft_dft(fti, ab, tE, xE);
                matpoly_ft_dft(fti, ab, tpiL, xpiL);
                /* length E_length * length k ==> length input_length
                 * with E_length = input_length + k - 1
                 *
                 * we'll shift by k-1 coefficients, because k is the
                 * smallest length
                 */
                matpoly_ft_mul(fti, ab, tER, tE, tpiL);
                xER->size = input_length;
                ASSERT_ALWAYS(xER->size <= xER->alloc);
                matpoly_ft_ift_mp(fti, ab, xER, tER, k-1);
                matpoly_ft_clear(fti, ab, tpiL);
                matpoly_ft_clear(fti, ab, tE);
                matpoly_ft_clear(ab, tER,  fti);
                x.set_since_last();
                s++;
            }
            mpz_clear(p);
#endif
            /* see remark above about matpoly_mul_caching gaining a stats
             * argument someday */
            tree_stats stats;
            for(small_bench<timer_t> x = finder.micro_bench(3); !x.done(); ++x) {
                xER = matpoly_ft<fft_type>::mp_caching_adj(stats, xE, xpiL, adj, nullptr);
                x.set_since_last();
            }
            if (xERref.get_size() == 0) {
                xERref = std::move(xER);
            } else if (xER.cmp(xERref) != 0) {
                fprintf(stderr, "MISMATCH3!\n");
            }
        }

        cout << input_length
            << " " << finder.summarize_for_this_length(k)
            << extra_info.str()
            << "\n";
    }
    hup_caught = 0;

    const auto table = finder.export_kara_cutoff_data(improved);
    polymat_set_mp_kara_cutoff(improved, nullptr);

    cout << "/* Cutoffs (0=basecase, 1=kara) for "
                << (m)<<"*"<<(m+n)
                <<" times "
                << (m+n)<<"*"<<(m+n)<<" middle-products */\n";
    cout << "#define MP_CUTOFFS_" <<m<<"_"<<(m+n)<<"_"<<(m+n)
        << " " << finder.print_result(table) << "\n";

}/*}}}*/

#if 0
/* 20150826: bigpolymat deleted */
void lingen_tune_bigmul(abdst_field ab, unsigned int m, unsigned int n, unsigned int m1, unsigned int n1, MPI_Comm comm)/*{{{*/
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    cxx_gmp_randstate rstate;
    gmp_randseed_ui(rstate, 1);
    /* arguments to the ctor:
     * 2 : we are benching 2 methods
     * 1 : size makes sense only for >=1
     */
    cutoff_finder<> finder(2, 1);

    bigpolymat model;
    bigpolymat_init_model(model, comm, m1, n1);

    /* Beware, k is a length, not a degree. Hence length 1 clearly makes
     * no sense */
    for(unsigned int k = 2 ; !finder.done() ; k=finder.next_length(k)) {
        /* Note: we are benching degree k, but the degree we are
         * interested in for pi is k*m/(m+n) */
        polymat pi, piL, piR;
        bigpolymat bpi, bpiL, bpiR;
        polymat_init(ab, piL, m+n, m+n, k);
        polymat_init(ab, piR, m+n, m+n, k);
        polymat_init(ab, pi, m+n, m+n, 2*k-1);
        bigpolymat_init(ab, bpiL, model, m+n, m+n, k);
        bigpolymat_init(ab, bpiR, model, m+n, m+n, k);
        bigpolymat_init(ab, bpi, model, m+n, m+n, 2*k-1);
        polymat_fill_random(ab, piL, k, rstate);
        polymat_fill_random(ab, piR, k, rstate);
        polymat_fill_random(ab, bigpolymat_my_cell(bpiL), k, rstate);
        polymat_fill_random(ab, bigpolymat_my_cell(bpiR), k, rstate);

        double ttloc;
        if (rank == 0) {
            for(small_bench<timer_t> x = finder.micro_bench<timer_wct>(ttloc); !x.done(); ++x) {
                polymat_mul(ab, pi, piL, piR);
                x.set_since_last();
            }
        }
        MPI_Bcast(&ttloc, 1, MPI_DOUBLE, 0, comm);

        double ttmpi;
        for(auto x = finder.micro_bench<timer_wct_synchronized>(ttmpi); !x.done(); ++x) {
            bigpolymat_mul(ab, bpi, bpiL, bpiR);
            x.set_since_last();
        }
        if (rank == 0)
            printf("%d %1.6f %1.6f %1.1f\n", k, ttloc, ttmpi, ttmpi/ttloc);
        finder.new_winner(k, ttmpi < ttloc); /* < : mpi wins: 1 */

        polymat_clear(ab, piL);
        polymat_clear(ab, piR);
        polymat_clear(ab, pi);
        bigpolymat_clear(ab, bpiL);
        bigpolymat_clear(ab, bpiR);
        bigpolymat_clear(ab, bpi);
    }

    bigpolymat_clear_model(model);

    if (rank == 0) {
        cout << "/* Cutoffs for "<<(m+n)<<"*"<<(m+n)<<"*"<<(m+n)<<" MPI products: */\n";
        cout << "#define MUL_MPI_CUTOFFS_" <<(m+n)<<"_"<<(m+n)<<"_"<<(m+n)
            << " " << finder.result() << "\n";
    }
}/*}}}*/
#endif

template<bool is_binary>
void lingen_tune_cutoffs(bw_dimensions<is_binary> & d, MPI_Comm comm MAYBE_UNUSED, cxx_param_list & pl)
{
    cxx_gmp_randstate rstate;

    int catchsig=0;

    setvbuf(stdout, nullptr, _IONBF, 0);
    setvbuf(stderr, nullptr, _IONBF, 0);

    param_list_parse(pl, "B", bench_atleast_uptothis);
    param_list_parse(pl, "catchsig", catchsig);

    typename matpoly<is_binary>::arith_hard * ab = & d.ab;
    unsigned int const m = d.m;
    unsigned int const n = d.n;
    cxx_mpz p(ab->characteristic());

    gmp_randseed_ui(rstate, 1);

    if (catchsig) {
#ifdef  HAVE_SIGHUP
        catch_control_signals();
#else
        fprintf(stderr, "Warning: support for signal trapping not enabled at compile time\n");
#endif
    }

    printf("# Tuning for m=%d n=%d p=[%zu %u-bits words]\n",
            m, n, mpz_size(p), GMP_LIMB_BITS);

    /* XXX BUG: with depth adjustment==0 early on, we get check failures.
     * Must investigate */

    cutoff_list cl_mp = nullptr;
    // lingen_tune_mp_fti_depth(ab, m, n, &cl_mp);
    lingen_tune_mp<is_binary>(ab, m, n, cl_mp);

    cutoff_list cl_mul = nullptr;
    // lingen_tune_mul_fti_depth(ab, m, n, &cl_mul);
    lingen_tune_mul<is_binary>(ab, m, n, cl_mul);

    // int tune_bm_basecase = 1;
    int const tune_mp = 1;
    /* record the list of cutoffs, and the method which wins from there
     */

    if (tune_mp) {
        /* Now for benching mp and plain mul */
        size_t const maxtune = 10000 / (m*n);
        /* Bench the products which would come together with a k-steps
         * basecase algorithm. IOW a 2k, one-level recursive call incurs
         * twice the k-steps basecase, plus once the timings counted here
         * (presently, this has Karatsuba complexity)
         */
        for(size_t k = 10 ; k < maxtune ; k+= k/10) {
            size_t const sE = k*(m+2*n)/(m+n);
            size_t const spi = k*m/(m+n);
            polymat E(ab, m, m+n, sE);
            polymat piL(ab, m+n, m+n, spi);
            polymat piR(ab, m+n, m+n, spi);
            polymat pi(ab, m+n, m+n, spi*2);
            polymat Er(ab, m, m+n, sE-spi+1);
            E.set_size(sE);
            for(size_t v = 0 ; v < E.m * E.n * E.get_size() ; v++) {
                ab->set_random(ab->vec_item(E.x, v), rstate);
            }
            piL.set_size(spi);
            piR.set_size(spi);
            for(size_t v = 0 ; v < piL.m * piL.n * piL.get_size() ; v++) {
                ab->set_random(ab->vec_item(piL.x, v), rstate);
                ab->set_random(ab->vec_item(piR.x, v), rstate);
            }
            double ttmp = 0, ttmul = 0;
            ttmp -= seconds();
            Er.mp(E, piL);
            ttmp += seconds();
            ttmul -= seconds();
            pi.mul(piL, piR);
            ttmul += seconds();
            double const ttmpq = ttmp / (k*k);
            double const ttmulq = ttmul / (k*k);
            double const ttmpk = ttmp / pow(k, 1.58);
            double const ttmulk = ttmul / pow(k, 1.58);
            printf("%zu [%.2e+%.2e = %.2e] [%.2e+%.2e = %.2e]\n",
                    k,
                    ttmpq, ttmulq, ttmpq + ttmulq,
                    ttmpk, ttmulk, ttmpk + ttmulk
                    );
            // (seconds()-tt) / (k*k)); // ((sE-spi) * spi) / (m*(m+n)*(m+n)));
            // printf("%zu %.2e\n", E.get_size(), (seconds()-tt) / (k*k)); // (spi * spi) / ((m+n)*(m+n)*(m+n)));
        }
    }

#if 0   /* for basecase with polymat */
    if (tune_bm_basecase) {
        bmstatus bm;
        bmstatus_init(bm, m, n);
        abfield_specify(bm->d.ab, MPFQ_PRIME_MPZ, (mpz_srcptr) p);
        unsigned int t0 = bm->t;
        unsigned int maxtune = 10000 / (m * n);
        for(unsigned int k = 10 ; k < maxtune ; k += k/10) {
            unsigned int * delta = new unsigned int[m + n];
            polymat E, pi;
            polymat_init(ab, pi, 0, 0, 0);
            polymat_init(ab, E, m, m+n, maxtune);
            E->size = k;
            for(unsigned int v = 0 ; v < E->m * E->n * E->size ; v++) {
                abrandom(ab, E->x[v], rstate);
            }
            double tt = seconds();
            for(unsigned int j = 0 ; j < m + n ; delta[j++]=1);
            bm->t = 1;
            bw_lingen_basecase(bm, pi, E, delta);
            printf("%zu %.2e\n", E->size, (seconds()-tt) / (k * k));
            polymat_clear(ab, pi);
            polymat_clear(ab, E);
            delete[] delta;
        }
        bmstatus_clear(bm);
    }
#endif

#if 0
    /* This should normally be reasonably quick, and running it every
     * time can be considered as an option */
    if (rank == 0) {
        lingen_tune_mul(ab, m, n);
        lingen_tune_mp(ab, m, n);
    }
    extern polymat_cutoff_info polymat_mul_kara_cutoff;
    extern polymat_cutoff_info polymat_mp_kara_cutoff;
    bigpolymat_bcast_polymat_cutoff(&polymat_mul_kara_cutoff, 0, comm);
    bigpolymat_bcast_polymat_cutoff(&polymat_mp_kara_cutoff, 0, comm);

    lingen_tune_bigmul(ab, m, n, mpi[0]*thr[0], mpi[1]*thr[1], comm);
#endif
}

template void lingen_tune_cutoffs<false>(bw_dimensions<false> & d, MPI_Comm comm, cxx_param_list & pl);
