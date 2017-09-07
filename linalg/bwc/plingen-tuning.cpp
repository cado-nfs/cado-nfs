#include "cado.h"

#include <cstddef>      /* see https://gcc.gnu.org/gcc-4.9/porting_to.html */
#include <sys/time.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>
#ifdef  HAVE_SIGHUP
#include <signal.h>
#endif

#include "portability.h"
#include "macros.h"
#include "utils.h"
#include "mpfq_layer.h"
#include "lingen-polymat.h"
#include "lingen-matpoly.h"
// #include "lingen-bigpolymat.h" // 20150826: deleted.
#ifdef  HAVE_MPIR
#include "lingen-matpoly-ft.h"
#endif
#include "plingen.h"
#include "plingen-tuning.h"

#include <vector>
#include <utility>
#include <map>
#include <string>
#include <ostream>
#include <iostream>
#include <sstream>

void plingen_tuning_decl_usage(param_list_ptr pl)
{
    param_list_decl_usage(pl, "B",
            "minimum end of bench span window");
    param_list_decl_usage(pl, "catchsig",
            "enable intercepting ^C");
}
void plingen_tuning_lookup_parameters(param_list_ptr pl)
{
    param_list_lookup_string(pl, "B");
    param_list_lookup_string(pl, "catchsig");
}

/* interface to C programs for list of cutoffs we compute *//*{{{*/
/* This is the simple C-type cutoff table which is used down the line by
 * C programs */

struct cutoff_list_s {
    unsigned int k;     /* this is the number of coefficients (length) of
                         * pi.  Admittedly,  it is slightly bogus to have
                         * this as a main variable.  Length [k] for pi
                         * corresponds to an "input length" of [k*(m+n)/m].
                         * In the MP case, this means that E has length
                         * [input_length+k-1], and the output E_right has
                         * length [input_length]. */
    unsigned int choice;        /* semantics vary. For FFT-based stuff,
                                 * this is the depth reduction from the
                                 * flit default.
                                 */
};

typedef cutoff_list_s * cutoff_list;

/* cutoff list *ALWAYS* finish with UINT_MAX, UINT_MAX */
unsigned int cutoff_list_get(cutoff_list cl, unsigned int k)
{
    if (!cl) return 0;
    if (k < cl->k) return 0;
    for( ; k >= cl->k ; cl++);
    return (--cl)->choice;
}
/*}}}*/

/* The -B argument may be used to request printing of timing results at
 * least up to this input length */
unsigned int bench_atleast_uptothis = 0;

using namespace std;

/* timer backends *//*{{{*/
#ifdef  HAVE_GCC_STYLE_AMD64_INLINE_ASM
struct timer_rdtsc {
    inline double operator()() const {
        uint64_t c = cputicks();
        return c/1.0e9;
    }
    static const char * timer_unit() { return "Gccyc"; }
};
#endif


struct timer_rusage {
    inline double operator()() const { return seconds(); }
    static const char * timer_unit() { return "seconds (rusage)"; }
};
struct timer_wct {
    inline double operator()() const { return wct_seconds(); }
    static const char * timer_unit() { return "seconds (wct)"; }
};
/* It is important, when doing MPI benches, that we use this algorithm */
struct timer_wct_synchronized {
    inline double operator()() const { 
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
    int enough_repeats;
    double minimum_time;
    double enough_time;
    measurement_choice() {
        /* defaults */
        enough_repeats = 1000;
        minimum_time = 0.5;
        enough_time = 1.0;
    }
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
    double tt;

    int n;
    explicit small_bench(double& t, measurement_choice mc = measurement_choice()) : mc(mc), tfinal(t) {tt=n=0;reset();}
    /* are we done ? when we are, we set tfinal to the average value. */
    inline int done() {
        if (tt < mc.enough_time) {
            if (tt < mc.minimum_time && n < mc.enough_repeats)
                return 0;
        }
        // fprintf(stderr, "%.2f %d\n", tt, n);
        tfinal = tt / n;
        return 1;
    }
    inline small_bench& operator++() { n++; return *this; }
    inline void reset() { last_starting_tt = T()(); }
    /* adds to tt the time since the last timer reset, and resets the
     * timer.  */
    inline void add_since_last(int weight = 1) { 
        double now = T()();
        tt += weight * (now - last_starting_tt);
        last_starting_tt = now;
    }
    /* *sets* tt to the time since the last timer reset. Do not reset the
     * timer.  */
    inline void set_since_last(int weight = 1) { 
        double now = T()();
        tt = weight * (now - last_starting_tt);
    }
    inline void inject(double t, int weight = 1) { tt += weight*t; }
};

/* {{{ cutoff_finder and its child small_bench */
template<typename Timer_backend = timer_rusage>
struct cutoff_finder {
    measurement_choice mc;
    double scale;
    int stable_cutoff_break;
    unsigned int ntests;
    int stable;
    vector<double> benches;
    vector<bool> meaningful;
    vector<pair<unsigned int, pair<vector<double>, int> > > all_results;
    map<int, string> method_names;
        double slowness_ratio;
        unsigned int age_slow_discard;

    cutoff_finder(unsigned int ntests, measurement_choice mc = measurement_choice())
        :
        mc(mc),
        ntests(ntests),
        benches(ntests),
        meaningful(ntests, true)
    {
        /* Fill defaults */
        scale = 1.1;
        stable_cutoff_break = 4;
        stable = 0;
        slowness_ratio = 2;
        age_slow_discard = 5;
    }

    inline void set_method_name(int k, string const& s) {
        method_names.insert(make_pair(k, s));
    }
    inline string method_name(int i) const {
        ostringstream v;
        map<int, string>::const_iterator z = method_names.find(i);
        if (z == method_names.end()) {
            v << "method " << i;
            return v.str();
        }
        return z->second;
    }

    inline unsigned int next_length(unsigned int k) const {
        return MAX(k+1, scale*k);
    }

    inline int done() const {
        int nactive=0;
        for(unsigned int i = 0 ; i < ntests ; i++) {
            nactive += meaningful[i];
        }
        if (nactive > 1) return 0;
        return !(stable < stable_cutoff_break);
    }

    inline bool still_meaningful_to_test(int i) { return meaningful[i]; }

    string summarize_for_this_length(unsigned int k)
    {
        std::ostringstream comments;
        int best = -1;
        for(unsigned int i = 0 ; i < ntests ; i++) {
            if (meaningful[i] && (best < 0 || benches[i] < benches[best])) best = i;
        }
        ASSERT_ALWAYS(best >= 0);
        all_results.push_back(make_pair(k,
                    make_pair(benches, best)));

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
                int best_there = all_results[all_results.size()-1-age_slow].second.second;
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

    /* This is really limited to karatsuba-like cuttofs. It's ugly */
    vector<pair<unsigned int, int> > export_best_table()
    {
        vector<pair<unsigned int, int> > steps;
        steps.push_back(make_pair(1,0));

        typedef vector<pair<unsigned int, pair<vector<double>, int> > >::const_iterator it_t;
        for(it_t it = all_results.begin() ; it != all_results.end() ; ++it) {
            unsigned int size = it->first;
            int best = it->second.second;
            if (steps.empty() || best != steps.back().second) {
                steps.push_back(make_pair(size, best));
            }
        }
        return steps;
    }
    vector<pair<unsigned int, int> > export_kara_cutoff_data(struct polymat_cutoff_info * dst)
    {
        /* This size will eventually feed the ->subdivide field for the
         * cutoff info */
        unsigned int first_kara_size = UINT_MAX;

        /* This size will eventually feed the ->cut field for the
         * cutoff info. For sizes above this, we will _always_ use
         * karatsuba. */
        unsigned int first_alwayskara_size = UINT_MAX;

        vector<pair<unsigned int, int> > steps;
        steps.push_back(make_pair(1,0));

        typedef vector<pair<unsigned int, pair<vector<double>, int> > >::const_iterator it_t;
        for(it_t it = all_results.begin() ; it != all_results.end() ; ++it) {
            unsigned int size = it->first;
            vector<double> const& benches(it->second.first);
            int best = it->second.second;
            /* In case fft wins, we invent something which will use
             * karatsuba still */
            if (best > 1) best = benches[1] < benches[0];
            if (steps.empty() || best != steps.back().second) {
                steps.push_back(make_pair(size, best));
                if (best == 1 && size < first_kara_size)
                    first_kara_size = size;
                /* assign it multiple times */
                first_alwayskara_size = size;
            }
        }
        dst->cut = first_alwayskara_size;
        dst->subdivide = first_kara_size;
        dst->table_size = steps.size();
        dst->table = (unsigned int (*)[2]) realloc(dst->table,
                dst->table_size * sizeof(unsigned int[2]));
        for(unsigned int v = 0 ; v < dst->table_size ; v++) {
            dst->table[v][0] = steps[v].first;
            dst->table[v][1] = steps[v].second;
        };
        return steps;
    }
    vector<pair<unsigned int, int> > export_kara_cutoff_data_force_kara_now(struct polymat_cutoff_info * dst, unsigned int size)
    {
        vector<double> allz(ntests);
        all_results.push_back(make_pair(size, make_pair(allz, 1)));
        vector<pair<unsigned int, int> > x = export_kara_cutoff_data(dst);
        all_results.pop_back();
        return x;
    }
    static string print_result(vector<pair<unsigned int, int> > const& tab) {
        ostringstream s;
        s << "{";
        typedef vector<pair<unsigned int, int> >::const_iterator it_t;
        for(it_t y = tab.begin() ; y != tab.end() ; y++) {
            s << " { " << y->first << ", " << y->second << " },";
        }
        s << " }";
        return s.str();
    }

};
/* }}} */

#ifdef  HAVE_SIGHUP
double last_hup = 0;
int hup_caught = 0;

void sighandler(int sig MAYBE_UNUSED)
{
    double t = wct_seconds();
    if (t < last_hup + 0.5) {
        fprintf(stderr, "Interrupt twice in half a second; exiting\n");
        exit(1);
    }
    last_hup = wct_seconds();
    hup_caught++;
}


void catch_control_signals()
{
    struct sigaction sa[1];
    memset(sa, 0, sizeof(sa));
    sa->sa_handler = sighandler;
    sigaction(SIGHUP, sa, NULL);
    sigaction(SIGINT, sa, NULL);

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

#ifdef HAVE_MPIR
void plingen_tune_mul_fti_depth(abdst_field ab, unsigned int m, unsigned int n, cutoff_list *cl_out)/*{{{*/
{
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 1);

    typedef timer_rdtsc timer_t;

    int nadjs=7;
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
    cout << "# Note: for input length k within plingen,"
        << " we use ncoeffs=k*m/(m+n) = "<<(double)m/(m+n)<<"*k\n";

    /* This is for forcing the bench to run until a large length. This is
     * unnecessary here */
    unsigned int min_bench = 0;

    /* Beware, k is the length of piL, not a degree. Hence length 1
     * clearly makes no sense */
    for(unsigned int k = 2 ; !hup_caught && (k < min_bench || !finder.done()) ; k=finder.next_length(k)) {
        unsigned int input_length = (m+n) * k / m;

        abvec   A,   B,   C;
        void  *tA, *tB, *tC;

        mpz_t p;
        mpz_init(p);
        abfield_characteristic(ab, p);

        abvec_init(ab, &A, k);
        abvec_init(ab, &B, k);
        abvec_init(ab, &C, 2*k-1);
        
        abvec_random(ab, A, k, rstate);
        abvec_random(ab, B, k, rstate);
        abvec_set_zero(ab, C, 2*k-1);

        ostringstream extra_info;

        for(int index = 0 ; index < nadjs ; index++) {
            struct fft_transform_info fti[1];
            fft_get_transform_info_fppol(fti, p, k, k, m+n);
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
            fft_get_transform_allocs(fft_alloc_sizes, fti);
            void * tt = malloc(fft_alloc_sizes[2]);
            void * qt = malloc(fft_alloc_sizes[1]);

            tA = malloc(fft_alloc_sizes[0]);
            tB = malloc(fft_alloc_sizes[0]);
            tC = malloc(fft_alloc_sizes[0]);

            typedef small_bench<timer_t> bt;
            for(bt x = finder.micro_bench(index); !x.done(); ++x) {

                double t_dftA=0;
                for(bt y = bt(t_dftA, mc); !y.done(); ++y) {
                    fft_transform_prepare(tA, fti);
                    fft_do_dft_fppol(tA, (mp_limb_t*)A, k, qt, fti, p);
                    y.set_since_last();
                }
                x.inject(t_dftA, (m+n)*(m+n));

                double t_dftB=0;
                for(bt y = bt(t_dftB, mc); !y.done(); ++y) {
                    fft_transform_prepare(tB, fti);
                    fft_do_dft_fppol(tB, (mp_limb_t*)B, k, qt, fti, p);
                    y.set_since_last();
                }
                x.inject(t_dftB, (m+n)*(m+n));

                double t_conv=0;
                for(bt y = bt(t_conv, mc); !y.done(); ++y) {
                    fft_transform_prepare(tC, fti);
                    fft_zero(tC, fti);
                    fft_addmul(tC, tA, tB, tt, qt, fti);
                    y.set_since_last();
                }
                x.inject(t_conv, (m+n)*(m+n)*(m+n));

                double t_iftC=0;
                for(bt y = bt(t_conv, mc); !y.done(); ++y) {
                    fft_transform_prepare(tC, fti);
                    fft_do_ift_fppol((mp_limb_t*)C, 2*k-1, tC, qt, fti, p);
                    y.set_since_last();
                }
                x.inject(t_iftC, (m+n)*(m+n));
            }
            free(tA);
            free(tB);
            free(tC);
        }

        mpz_clear(p);

        cout << input_length
            << " " << finder.summarize_for_this_length(k)
            << extra_info.str()
            << "\n";
    }
    hup_caught = 0;

    vector<pair<unsigned int, int> > table = finder.export_best_table();

    typedef vector<pair<unsigned int, int> >::iterator it_t;
    for(it_t x = table.begin() ; x != table.end() ; ++x)
        x->second = nadjs-1-x->second;

    cout << "/* FFT depth adjustments for "
                << (m)<<"*"<<(m+n)
                <<" times "
                << (m+n)<<"*"<<(m+n)<<" products */\n";
    cout << "#define MUL_FTI_DEPTH_ADJ_" <<m+n<<"_"<<(m+n)<<"_"<<(m+n)
        << " " << finder.print_result(table) << endl;

    if (cl_out) {
        *cl_out = (cutoff_list) malloc((table.size()+1)*sizeof(**cl_out));
        for(unsigned int i = 0 ; i < table.size(); i++) {
            (*cl_out)[i].k = table[i].first;
            (*cl_out)[i].choice = table[i].second;
        }
        (*cl_out)[table.size()].k      = UINT_MAX;
        (*cl_out)[table.size()].choice = UINT_MAX;
    }

    gmp_randclear(rstate);
}/*}}}*/
void plingen_tune_mp_fti_depth(abdst_field ab, unsigned int m, unsigned int n, cutoff_list * cl_out)/*{{{*/
{
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 1);

    typedef timer_rdtsc timer_t;

    int nadjs=7;
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
    cout << "# Note: for input length k within plingen,"
        << " we use ncoeffs=k*m/(m+n) = "<<(double)m/(m+n)<<"*k\n";

    /* This is for forcing the bench to run until a large length. This is
     * unnecessary here */
    unsigned int min_bench = 0;

    /* Beware, k is the length of piL, not a degree. Hence length 1
     * clearly makes no sense */
    for(unsigned int k = 2 ; !hup_caught && (k < min_bench || !finder.done()) ; k=finder.next_length(k)) {
        unsigned int input_length = (m+n) * k / m;

        unsigned int E_length = k + input_length - 1;

        abvec   A,   B,   C;
        void  *tA, *tB, *tC;

        mpz_t p;
        mpz_init(p);
        abfield_characteristic(ab, p);

        abvec_init(ab, &A, E_length);
        abvec_init(ab, &B, k);
        abvec_init(ab, &C, input_length);
        
        abvec_random(ab, A, k, rstate);
        abvec_random(ab, B, k, rstate);
        abvec_set_zero(ab, C, k);

        ostringstream extra_info;

        for(int index = 0 ; index < nadjs ; index++) {
            struct fft_transform_info fti[1];
            fft_get_transform_info_fppol_mp(fti, p, k, E_length, m+n);
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
            fft_get_transform_allocs(fft_alloc_sizes, fti);
            void * tt = malloc(fft_alloc_sizes[2]);
            void * qt = malloc(fft_alloc_sizes[1]);

            tA = malloc(fft_alloc_sizes[0]);
            tB = malloc(fft_alloc_sizes[0]);
            tC = malloc(fft_alloc_sizes[0]);

            typedef small_bench<timer_t> bt;
            for(bt x = finder.micro_bench(index); !x.done(); ++x) {

                double t_dftA=0;
                for(bt y = bt(t_dftA, mc); !y.done(); ++y) {
                    fft_transform_prepare(tA, fti);
                    fft_do_dft_fppol(tA, (mp_limb_t*)A, E_length, qt, fti, p);
                    y.set_since_last();
                }
                x.inject(t_dftA, m*(m+n));

                double t_dftB=0;
                for(bt y = bt(t_dftB, mc); !y.done(); ++y) {
                    fft_transform_prepare(tB, fti);
                    fft_do_dft_fppol(tB, (mp_limb_t*)B, k, qt, fti, p);
                    y.set_since_last();
                }
                x.inject(t_dftB, (m+n)*(m+n));

                double t_conv=0;
                for(bt y = bt(t_conv, mc); !y.done(); ++y) {
                    fft_transform_prepare(tC, fti);
                    fft_zero(tC, fti);
                    fft_addmul(tC, tA, tB, tt, qt, fti);
                    y.set_since_last();
                }
                x.inject(t_conv, m*(m+n)*(m+n));

                double t_iftC=0;
                for(bt y = bt(t_conv, mc); !y.done(); ++y) {
                    fft_transform_prepare(tC, fti);
                    fft_do_ift_fppol_mp((mp_limb_t*)C, input_length, tC, qt, fti, p, k-1);
                    y.set_since_last();
                }
                x.inject(t_iftC, m*(m+n));
            }
            free(tA);
            free(tB);
            free(tC);
        }

        mpz_clear(p);

        cout << input_length
            << " " << finder.summarize_for_this_length(k)
            << extra_info.str()
            << "\n";
    }
    hup_caught = 0;

    vector<pair<unsigned int, int> > table = finder.export_best_table();

    typedef vector<pair<unsigned int, int> >::iterator it_t;
    for(it_t x = table.begin() ; x != table.end() ; ++x)
        x->second = nadjs-1-x->second;

    cout << "/* FFT depth adjustments for "
                << (m)<<"*"<<(m+n)
                <<" times "
                << (m+n)<<"*"<<(m+n)<<" middle-products */\n";
    cout << "#define MP_FTI_DEPTH_ADJ_" <<m<<"_"<<(m+n)<<"_"<<(m+n)
        << " " << finder.print_result(table) << endl;

    if (cl_out) {
        *cl_out = (cutoff_list) malloc((table.size()+1)*sizeof(**cl_out));
        for(unsigned int i = 0 ; i < table.size(); i++) {
            (*cl_out)[i].k = table[i].first;
            (*cl_out)[i].choice = table[i].second;
        }
        (*cl_out)[table.size()].k = UINT_MAX;
        (*cl_out)[table.size()].choice = UINT_MAX;
    }

    gmp_randclear(rstate);
}/*}}}*/
#endif  /* HAVE_MPIR */


void plingen_tune_mul(abdst_field ab, unsigned int m, unsigned int n, cutoff_list cl MAYBE_UNUSED)/*{{{*/
{
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 1);
#ifdef HAVE_MPIR
#define TUNE_MUL_FINDER_NMETHODS 4
#else  /* HAVE_MPIR */
#define TUNE_MUL_FINDER_NMETHODS 3
#endif  /* HAVE_MPIR */

    typedef timer_rusage timer_t;
    cutoff_finder<timer_rusage> finder(TUNE_MUL_FINDER_NMETHODS);
    finder.set_method_name(0, "polymat-basecase");
    finder.set_method_name(1, "polymat-karatsuba");
    finder.set_method_name(2, "matpoly-kronecker-schönhage");
#ifdef  HAVE_MPIR
    finder.set_method_name(3, "matpoly-ft-kronecker-schönhage-caching");
#endif

    polymat_cutoff_info always_basecase[1];
    polymat_cutoff_info improved[1];

    polymat_cutoff_info_init(always_basecase);
    polymat_cutoff_info_init(improved);

    cout << "# Tuning "
                << (m+n)<<"*"<<(m+n)
                <<" times "
                << (m+n)<<"*"<<(m+n)<<" products\n";

    cout << "# inputlength, ncoeffs";
    for(unsigned int i = 0 ; i < finder.ntests ; i++) {
        cout << ", " << finder.method_name(i);
    }
    cout << "\n";
    cout << "# Note: for input length k within plingen,"
        << " we use ncoeffs=k*m/(m+n) = "<<(double)m/(m+n)<<"*k\n";
    /* input length k means we consider the pi polynomial which is
     * created by k successive steps. So this mulitplication, in effect,
     * occurs at the end of a recursive procedure of length 2k, which
     * collects the pi matrices of its two child calls.
     */

    /* make sure we don't stop the bench absurdly early. We set the bench
     * as input_length \approx 2000 */
    unsigned int min_bench = bench_atleast_uptothis * m / (m+n);

    /* Beware, k is the length of piL, not a degree. Hence length 1
     * clearly makes no sense */
    for(unsigned int k = 2 ; !hup_caught && (k < min_bench || !finder.done()) ; k=finder.next_length(k)) {
        unsigned int input_length = (m+n) * k / m;
        polymat pi, piref, piL, piR;
        matpoly xpi, xpiref, xpiL, xpiR;
        polymat_init(ab, piL, m+n, m+n, k);
        polymat_init(ab, piR, m+n, m+n, k);
        polymat_init(ab, pi, m+n, m+n, 2*k-1);
        polymat_init(ab, piref, m+n, m+n, 2*k-1);
        polymat_fill_random(ab, piL, k, rstate);
        polymat_fill_random(ab, piR, k, rstate);
        
        /* of course we know the allocation needs, but we'll lazy-allocate
         * here, to save ram */
        matpoly_init(ab, xpiL, 0, 0, 0);
        matpoly_init(ab, xpiR, 0, 0, 0);
        matpoly_init(ab, xpi, 0, 0, 0);
        matpoly_init(ab, xpiref, 0, 0, 0);

        ostringstream extra_info;

        if (finder.still_meaningful_to_test(0)) {
            /* disable kara for a minute */
            polymat_set_mul_kara_cutoff(always_basecase, NULL);
            for(small_bench<timer_t> x = finder.micro_bench(0); !x.done(); ++x) {
                polymat_mul(ab, pi, piL, piR);
                x.set_since_last();
            }
            if (piref->size == 0) {
                polymat_swap(pi, piref);
                // fprintf(stderr, "BASIS0\n");
            } else if (polymat_cmp(ab, pi, piref) != 0) {
                fprintf(stderr, "MISMATCH0!\n");
            }
        }

        if (finder.still_meaningful_to_test(1)) {
            /* This temporarily sets the cutoff table to enable karatsuba for
             * length >=k (hence for this test) at least, possibly using
             * karatsuba one or more times in the recursive calls depending
             * on what has been measured as best so far.
             */
            finder.export_kara_cutoff_data_force_kara_now(improved, k);
            polymat_set_mul_kara_cutoff(improved, NULL);
            for(small_bench<timer_t> x = finder.micro_bench(1); !x.done(); ++x) {
                polymat_mul(ab, pi, piL, piR);
                x.set_since_last();
            }
            if (piref->size == 0) {
                polymat_swap(pi, piref);
                // fprintf(stderr, "BASIS1\n");
            } else if (polymat_cmp(ab, pi, piref) != 0) {
                fprintf(stderr, "MISMATCH1!\n");
            }
        }

        /* don't exaggerate our memory requirements */
        matpoly_set_polymat(ab, xpiL, piL);
        polymat_clear(ab, piL);
        matpoly_set_polymat(ab, xpiR, piR);
        polymat_clear(ab, piR);
        if (piref->size) {
            matpoly_set_polymat(ab, xpiref, piref);
        }
        polymat_clear(ab, piref);
        matpoly_init(ab, xpi, m+n, m+n, 2*k-1);
        polymat_clear(ab, pi);

        /* we should make the effort of converting polymat to matpoly,
         * right ? */

        if (finder.still_meaningful_to_test(2)) {
            /* The matpoly layer is just completetly different -- and gets
             * faster quite early on... */
            for(small_bench<timer_t> x = finder.micro_bench(2); !x.done(); ++x) {
                matpoly_mul(ab, xpi, xpiL, xpiR);
                x.set_since_last();
            }
            if (xpiref->size == 0) {
                matpoly_swap(xpi, xpiref);
                // fprintf(stderr, "BASIS2\n");
            } else if (matpoly_cmp(ab, xpi, xpiref) != 0) {
                fprintf(stderr, "MISMATCH2!\n");
            }
        }

#ifdef HAVE_MPIR
        if (finder.still_meaningful_to_test(3)) {
            unsigned int adj = cutoff_list_get(cl, k);
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
                matpoly_ft_init(ab, tpiL, xpiL->m, xpiL->n, fti);
                matpoly_ft_init(ab, tpiR, xpiR->m, xpiR->n, fti);
                matpoly_ft_init(ab, tpi, xpiL->m, xpiR->n, fti);
                matpoly_ft_dft(ab, tpiL, xpiL, fti);
                matpoly_ft_dft(ab, tpiR, xpiR, fti);
                matpoly_ft_mul(ab, tpi, tpiL, tpiR, fti);
                xpi->size = xpiL->size + xpiR->size - 1;
                ASSERT_ALWAYS(xpi->size <= xpi->alloc);
                matpoly_ft_ift(ab, xpi, tpi, fti);
                matpoly_ft_clear(ab, tpiL, fti);
                matpoly_ft_clear(ab, tpiR, fti);
                matpoly_ft_clear(ab, tpi,  fti);
                x.set_since_last();
                s++;
            }
            mpz_clear(p);
#else
            for(small_bench<timer_t> x = finder.micro_bench(3); !x.done(); ++x) {
                matpoly_mul_caching_adj(ab, xpi, xpiL, xpiR, adj);
                x.set_since_last();
            }
#endif
            if (xpiref->size == 0) {
                matpoly_swap(xpi, xpiref);
            } else if (matpoly_cmp(ab, xpi, xpiref) != 0) {
                fprintf(stderr, "MISMATCH3!\n");
            }
#if 0
#ifdef  TIME_FFT
            char msg[256];
            snprintf(msg, sizeof(msg),
                    " [depth %d, %d*%.2f dft %d*%.2f conv %d*%.2f ift]",
                    (int) fti->depth,
                    fti->dft.n/s, 1000*fti->dft.t/fti->dft.n,
                    fti->conv.n/s, 1000*fti->conv.t/fti->conv.n,
                    fti->ift.n/s, 1000*fti->ift.t/fti->ift.n);
            extra_info << msg;
            /* the following one-liner may be used to grab only
             * fft-related data */
    // perl -ne '/^\d+\s(\d+).*depth (\d+).*\*([\d\.]+)\sdft.*\*([\d\.]+)\sconv.*\*([\d\.]+)\sift/ && print "$1 $2 $3 $4 $5\n";' 
#endif
#endif
        }
#endif

        cout << input_length
            << " " << finder.summarize_for_this_length(k)
            << extra_info.str()
            << "\n";

        polymat_clear(ab, piL);
        polymat_clear(ab, piR);
        polymat_clear(ab, pi);
        polymat_clear(ab, piref);
        matpoly_clear(ab, xpiL);
        matpoly_clear(ab, xpiR);
        matpoly_clear(ab, xpi);
        matpoly_clear(ab, xpiref);
    }
    hup_caught = 0;


    vector<pair<unsigned int, int> > table = finder.export_kara_cutoff_data(improved);
    polymat_set_mul_kara_cutoff(improved, NULL);

    cout << "/* Cutoffs (0=basecase, 1=kara) for "
                << (m+n)<<"*"<<(m+n)
                <<" times "
                << (m+n)<<"*"<<(m+n)<<" products: */\n";
    cout << "#define MUL_CUTOFFS_" <<(m+n)<<"_"<<(m+n)<<"_"<<(m+n)
        << " " << finder.print_result(table) << endl;
    gmp_randclear(rstate);

    polymat_cutoff_info_clear(always_basecase);
    polymat_cutoff_info_clear(improved);
}/*}}}*/

void plingen_tune_mp(abdst_field ab, unsigned int m, unsigned int n, cutoff_list cl MAYBE_UNUSED)/*{{{*/
{
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    gmp_randseed_ui(rstate, 1);
#ifdef HAVE_MPIR
#define TUNE_MP_FINDER_NMETHODS 4
#else  /* HAVE_MPIR */
#define TUNE_MP_FINDER_NMETHODS 3
#endif  /* HAVE_MPIR */

    typedef timer_rusage timer_t;
    cutoff_finder<timer_rusage> finder(TUNE_MP_FINDER_NMETHODS);
    finder.set_method_name(0, "polymat-basecase");
    finder.set_method_name(1, "polymat-karatsuba");
    finder.set_method_name(2, "matpoly-kronecker-schönhage");
#ifdef  HAVE_MPIR
    finder.set_method_name(3, "matpoly-ft-kronecker-schönhage-caching");
#endif

    polymat_cutoff_info always_basecase[1];
    polymat_cutoff_info improved[1];

    polymat_cutoff_info_init(always_basecase);
    polymat_cutoff_info_init(improved);

    cout << "# Tuning "
                << (m)<<"*"<<(m+n)
                <<" times "
                << (m+n)<<"*"<<(m+n)<<" middle-products\n";

    cout << "# inputlength, ncoeffs";
    for(unsigned int i = 0 ; i < finder.ntests ; i++) {
        cout << ", " << finder.method_name(i);
    }
    cout << "\n";

    cout << "# Note: for input length k within plingen, "
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
    unsigned int min_bench = bench_atleast_uptothis * m / (m+n);

    /* Beware, k is a length of piL, not a degree. Hence length 1 clearly
     * makes no sense */
    for(unsigned int k = 2 ; !hup_caught && (k < min_bench || !finder.done()) ; k=finder.next_length(k)) {
        unsigned int input_length = (m+n) * k / m;
        /* MP(degree a, degree b) -> degree b-a
         *    length la=a+1, length lb=b+1 -> length b-a+1 = lb-la+1
         * we do want lc=input_length, we have set la=k, whence lb =
         * k+input_length-1
         */
        unsigned int E_length = k + input_length - 1;
        polymat ER, E, piL, ERref;
        matpoly xER, xE, xpiL, xERref;
        polymat_init(ab, E, m, m+n, E_length);
        polymat_init(ab, piL, m+n, m+n, k);
        polymat_init(ab, ER, m, m+n, input_length);
        polymat_init(ab, ERref, m, m+n, input_length);
        polymat_fill_random(ab, E, E_length, rstate);
        polymat_fill_random(ab, piL, k, rstate);

        matpoly_init(ab, xE, m, m+n, E_length);
        matpoly_init(ab, xpiL, m+n, m+n, k);
        matpoly_init(ab, xER, m, m+n, input_length);
        matpoly_init(ab, xERref, m, m+n, input_length);

        ostringstream extra_info;

        if (finder.still_meaningful_to_test(0)) {
            /* disable kara for a minute */
            polymat_set_mp_kara_cutoff(always_basecase, NULL);
            for(small_bench<timer_t> x = finder.micro_bench(0); !x.done(); ++x) {
                polymat_mp(ab, ER, E, piL);
                x.set_since_last();
            }
            if (ERref->size == 0) {
                polymat_swap(ER, ERref);
                // fprintf(stderr, "BASIS0\n");
            } else if (polymat_cmp(ab, ER, ERref) != 0) {
                fprintf(stderr, "MISMATCH0!\n");
            }
        }

        if (finder.still_meaningful_to_test(1)) {
            /* This temporarily sets the cutoff table to enable karatsuba for
             * length >=k (hence for this test) at least, possibly using
             * karatsuba one or more times in the recursive calls depending
             * on what has been measured as best so far.
             */
            finder.export_kara_cutoff_data_force_kara_now(improved, k);
            polymat_set_mp_kara_cutoff(improved, NULL);
            for(small_bench<timer_t> x = finder.micro_bench(1); !x.done(); ++x) {
                polymat_mp(ab, ER, E, piL);
                x.set_since_last();
            }
            if (ERref->size == 0) {
                polymat_swap(ER, ERref);
                // fprintf(stderr, "BASIS1\n");
            } else if (polymat_cmp(ab, ER, ERref) != 0) {
                fprintf(stderr, "MISMATCH1!\n");
            }
        }

        /* don't exaggerate our memory requirements */
        matpoly_set_polymat(ab, xpiL, piL);
        polymat_clear(ab, piL);
        matpoly_set_polymat(ab, xE, E);
        polymat_clear(ab, E);
        if (ERref->size) {
            matpoly_set_polymat(ab, xERref, ERref);
        }
        polymat_clear(ab, ERref);
        matpoly_init(ab, xER, m, m+n, input_length);
        polymat_clear(ab, ER);

        /* we should make the effort of converting polymat to matpoly,
         * right ? */

        if (finder.still_meaningful_to_test(2)) {
            /* The matpoly layer is just completetly different -- and gets
             * faster quite early on... */
            for(small_bench<timer_t> x = finder.micro_bench(2); !x.done(); ++x) {
                matpoly_mp(ab, xER, xE, xpiL);
                x.set_since_last();
            }
            if (xERref->size == 0) {
                matpoly_swap(xER, xERref);
                // fprintf(stderr, "BASIS2\n");
            } else if (matpoly_cmp(ab, xER, xERref) != 0) {
                fprintf(stderr, "MISMATCH2!\n");
            }
        }

#ifdef HAVE_MPIR
        if (finder.still_meaningful_to_test(3)) {
            unsigned int adj = UINT_MAX;
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
                matpoly_ft_init(ab, tpiL, xpiL->m, xpiL->n, fti);
                matpoly_ft_init(ab, tE, xE->m, xE->n, fti);
                matpoly_ft_init(ab, tER, xE->m, xpiL->n, fti);
                matpoly_ft_dft(ab, tE, xE, fti);
                matpoly_ft_dft(ab, tpiL, xpiL, fti);
                /* length E_length * length k ==> length input_length
                 * with E_length = input_length + k - 1
                 *
                 * we'll shift by k-1 coefficients, because k is the
                 * smallest length
                 */
                matpoly_ft_mul(ab, tER, tE, tpiL, fti);
                xER->size = input_length;
                ASSERT_ALWAYS(xER->size <= xER->alloc);
                matpoly_ft_ift_mp(ab, xER, tER, k-1, fti);
                matpoly_ft_clear(ab, tpiL, fti);
                matpoly_ft_clear(ab, tE, fti);
                matpoly_ft_clear(ab, tER,  fti);
                x.set_since_last();
                s++;
            }
            mpz_clear(p);
#endif
            for(small_bench<timer_t> x = finder.micro_bench(3); !x.done(); ++x) {
                matpoly_mp_caching_adj(ab, xER, xE, xpiL, adj);
                x.set_since_last();
            }
            if (xERref->size == 0) {
                matpoly_swap(xER, xERref);
            } else if (matpoly_cmp(ab, xER, xERref) != 0) {
                fprintf(stderr, "MISMATCH3!\n");
            }
#if 0
#ifdef  TIME_FFT
            char msg[256];
            snprintf(msg, sizeof(msg),
                    " [depth %d, %d*%.2f dft %d*%.2f conv %d*%.2f ift]",
                    (int) fti->depth,
                    fti->dft.n/s, 1000*fti->dft.t/fti->dft.n,
                    fti->conv.n/s, 1000*fti->conv.t/fti->conv.n,
                    fti->ift.n/s, 1000*fti->ift.t/fti->ift.n);
            extra_info << msg;
            /* the following one-liner may be used to grab only
             * fft-related data */
    // perl -ne '/^\d+\s(\d+).*depth (\d+).*\*([\d\.]+)\sdft.*\*([\d\.]+)\sconv.*\*([\d\.]+)\sift/ && print "$1 $2 $3 $4 $5\n";' 
#endif
#endif
        }
#endif

        cout << input_length
            << " " << finder.summarize_for_this_length(k)
            << extra_info.str()
            << "\n";

        polymat_clear(ab, E);
        polymat_clear(ab, piL);
        polymat_clear(ab, ER);
        polymat_clear(ab, ERref);
        matpoly_clear(ab, xE);
        matpoly_clear(ab, xpiL);
        matpoly_clear(ab, xER);
        matpoly_clear(ab, xERref);
    }
    hup_caught = 0;

    vector<pair<unsigned int, int> > table = finder.export_kara_cutoff_data(improved);
    polymat_set_mp_kara_cutoff(improved, NULL);

    cout << "/* Cutoffs (0=basecase, 1=kara) for "
                << (m)<<"*"<<(m+n)
                <<" times "
                << (m+n)<<"*"<<(m+n)<<" middle-products */\n";
    cout << "#define MP_CUTOFFS_" <<m<<"_"<<(m+n)<<"_"<<(m+n)
        << " " << finder.print_result(table) << endl;
    gmp_randclear(rstate);

    polymat_cutoff_info_clear(always_basecase);
    polymat_cutoff_info_clear(improved);
}/*}}}*/

#if 0
/* 20150826: bigpolymat deleted */
void plingen_tune_bigmul(abdst_field ab, unsigned int m, unsigned int n, unsigned int m1, unsigned int n1, MPI_Comm comm)/*{{{*/
{
    int rank;
    MPI_Comm_rank(comm, &rank);
    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
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
            << " " << finder.result() << endl;
    }
    gmp_randclear(rstate);
}/*}}}*/
#endif

void plingen_tuning(abdst_field ab, unsigned int m, unsigned int n, MPI_Comm comm, param_list_ptr pl)
{
    mpz_t p;
    int thr[2] = {1,1};
    int mpi[2] = {1,1};
    int rank;
    int catchsig=0;

    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    param_list_parse_uint(pl, "B", &bench_atleast_uptothis);
    param_list_parse_int(pl, "catchsig", &catchsig);

    MPI_Comm_rank(comm, &rank);

    mpz_init(p);
    abfield_characteristic(ab, p);
    param_list_parse_intxint(pl, "mpi", mpi);
    param_list_parse_intxint(pl, "thr", thr);
    if (rank == 0) {
        printf("# Tuning for m=%d n=%d p=[%zu %u-bits words]"
                " mpi=%dx%d thr=%dx%d\n",
                m, n, mpz_size(p), GMP_LIMB_BITS,
                mpi[0], mpi[1], thr[0], thr[1]);
    }
    mpz_clear(p);

    if (catchsig) {
#ifdef  HAVE_SIGHUP
        catch_control_signals();
#else
        fprintf(stderr, "Warning: support for signal trapping not enabled at compile time\n");
#endif
    }

    /* XXX BUG: with depth adjustment==0 early on, we get check failures.
     * Must investigate */
    cutoff_list cl_mul = NULL;
    cutoff_list cl_mp = NULL;

    // plingen_tune_mp_fti_depth(ab, m, n, &cl_mp);
    // plingen_tune_mul_fti_depth(ab, m, n, &cl_mul);
    plingen_tune_mp(ab, m, n, cl_mp);
    plingen_tune_mul(ab, m, n, cl_mul);

#if 0
    /* This should normally be reasonably quick, and running it every
     * time can be considered as an option */
    if (rank == 0) {
        plingen_tune_mul(ab, m, n);
        plingen_tune_mp(ab, m, n);
    }
    extern polymat_cutoff_info polymat_mul_kara_cutoff;
    extern polymat_cutoff_info polymat_mp_kara_cutoff;
    bigpolymat_bcast_polymat_cutoff(&polymat_mul_kara_cutoff, 0, comm);
    bigpolymat_bcast_polymat_cutoff(&polymat_mp_kara_cutoff, 0, comm);

    plingen_tune_bigmul(ab, m, n, mpi[0]*thr[0], mpi[1]*thr[1], comm);
#endif

#if 0
    int tune_bm_basecase = 1;
    int tune_mp = 1;
    /* record the list of cutoffs, and the method which wins from there
     */

    if (tune_mp) {
        /* Now for benching mp and plain mul */
        unsigned int maxtune = 10000 / (m*n);
        /* Bench the products which would come together with a k-steps
         * basecase algorithm. IOW a 2k, one-level recursive call incurs
         * twice the k-steps basecase, plus once the timings counted here
         * (presently, this has Karatsuba complexity)
         */
        for(unsigned int k = 10 ; k < maxtune ; k+= k/10) {
            polymat E, piL, piR, pi, Er;
            unsigned int sE = k*(m+2*n)/(m+n);
            unsigned int spi = k*m/(m+n);
            polymat_init(E, m, m+n, sE);
            polymat_init(piL, m+n, m+n, spi);
            polymat_init(piR, m+n, m+n, spi);
            polymat_init(pi, m+n, m+n, spi*2);
            polymat_init(Er, m, m+n, sE-spi+1);
            E->size = sE;
            for(unsigned int v = 0 ; v < E->m * E->n * E->size ; v++) {
                abrandom(ab, E->x[v], rstate);
            }
            piL->size = spi;
            piR->size = spi;
            for(unsigned int v = 0 ; v < piL->m * piL->n * piL->size ; v++) {
                abrandom(ab, piL->x[v], rstate);
                abrandom(ab, piR->x[v], rstate);
            }
            double ttmp = 0, ttmul = 0;
            ttmp -= seconds();
            polymat_mp(ab, Er, E, piL);
            ttmp += seconds();
            ttmul -= seconds();
            polymat_mul(ab, pi, piL, piR);
            ttmul += seconds();
            double ttmpq = ttmp / (k*k);
            double ttmulq = ttmul / (k*k);
            double ttmpk = ttmp / pow(k, 1.58);
            double ttmulk = ttmul / pow(k, 1.58);
            printf("%u [%.2e+%.2e = %.2e] [%.2e+%.2e = %.2e]\n",
                    k,
                    ttmpq, ttmulq, ttmpq + ttmulq,
                    ttmpk, ttmulk, ttmpk + ttmulk
                    );
            // (seconds()-tt) / (k*k)); // ((sE-spi) * spi) / (m*(m+n)*(m+n)));
            // printf("%zu %.2e\n", E->size, (seconds()-tt) / (k*k)); // (spi * spi) / ((m+n)*(m+n)*(m+n)));
            polymat_clear(E);
            polymat_clear(piL);
            polymat_clear(piR);
            polymat_clear(pi);
            polymat_clear(Er);
        }
    }
    if (tune_bm_basecase) {
        unsigned int maxtune = 10000 / (m * n);
        for(unsigned int k = 10 ; k < maxtune ; k += k/10) {
            unsigned int * delta = malloc((m + n) * sizeof(unsigned int));
            polymat E, pi;
            polymat_init(pi, 0, 0, 0);
            polymat_init(E, m, m+n, maxtune);
            E->size = k;
            for(unsigned int v = 0 ; v < E->m * E->n * E->size ; v++) {
                abrandom(ab, E->x[v], rstate);
            }
            double tt = seconds();
            for(unsigned int j = 0 ; j < m + n ; delta[j++]=1);
            bm->t = 1;
            bw_lingen_basecase(bm, pi, E, delta);
            printf("%zu %.2e\n", E->size, (seconds()-tt) / (k * k));
            polymat_clear(pi);
            polymat_clear(E);
            free(delta);
        }
    }


#endif

    return;
}

