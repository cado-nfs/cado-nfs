#ifndef CADO_LAS_SIEVER_CONFIG_HPP
#define CADO_LAS_SIEVER_CONFIG_HPP

#include <cstring>

#include <algorithm>
#include <map>
#include <vector>
#include <tuple>
#include <utility>

#include "fb.hpp"
#include "params.h"
#include "las-side-config.hpp"

#include "las-special-q-task-tree.hpp"
#include "las-special-q-task-simple.hpp"

#include "sieve-methods.hpp"

struct special_q; // IWYU pragma: keep

/* siever_config */
 
/* The following structure lists the fields with an impact on the siever.
 * Different values for these fields will correspond to different siever
 * structures.
 */
struct siever_config {
    /* The bit size of the special-q. Counting in bits is no necessity,
     * we could imagine being more accurate. This is set by
     * siever_config_pool::get_config_for_q  */
    // unsigned int bitsize __attribute__((deprecated));  /* bitsize == 0 indicates end of table */
    // int side __attribute__((deprecated));              /* special-q side */

    int logA;
    int logI;   /* see below. logI is initialized late in the game */

    int adjust_strategy = 0;

    /* This does not really belong here. I'd rather have it at the
     * las_info level. However for obscure reasons,
     * las.get_strategies(siever_config&)
     * wants it.
     *
     * FIXME: er. now that get_strategies lives below las, we can do this
     * move, right ?
     */
    unsigned int sublat_bound;


    /* These four parameters are as they are provided in the command
     * line. In truth, the ones that really matter are the ones in the
     * fb_factorbase::key_type object that is stored within the
     * nfs_work structure (in sides[side].fbK).
     *
     * In particular, bucket_thresh and bucket_thresh1 below have a
     * default value that is dependent on:
     *  - logI -- which is not set here. Well, it exists here, but it is
     *    set late in the game, after all other fields
     *    (sieve_range_adjust does that).  And we don't want to go the
     *    route of making the default value for some of the fields in
     *    this struct be dynamic depending on a setter function for logI
     *  - the side, as well.
     */
    private:
    /* access should rather be
     * via ws.sides[side].fbK.thresholds[0,1]
     * and ws.sides[side].fbK.{td_thresh, skipped}
     */
    unsigned long bucket_thresh = 0;  // bucket sieve primes >= bucket_thresh
#if MAX_TOPLEVEL >= 2
    unsigned long bucket_thresh1 = 0; // primes above are 2-level bucket-sieved
#endif
#if MAX_TOPLEVEL >= 3
    unsigned long bucket_thresh2 = 0; // primes above are 3-level bucket-sieved
#endif
    static_assert(MAX_TOPLEVEL == 3);

    unsigned int td_thresh = 1024;
    unsigned int skipped = 1;         // don't sieve below this

    public:
    /* the only way to access the four fields above */
    fb_factorbase::key_type instantiate_thresholds(int side) const;

    public:

    /* unsieve threshold is not related to the factor base. */
    unsigned int unsieve_thresh = 100;

    std::vector<siever_side_config> sides;

    /*
     * The logic behind skipping the resieving is as follows. Briefly
     * put, it should be a two-sides thing only.
     *
     * When there are two sides, and only side s is sieved (lim>0), then
     * side 1-s is handled by cofactorization only (probably batch). For
     * this to make the slightest bit of sense, it has to be that only a
     * tiny fraction of the (a,b) pairs survive sieving on side
     * s. So once we've identified these survivors on side s, we move on
     * to cofactoring on side 1-s (which will keep a fraction of its
     * input), and a priori later on finish with cofactoring on side s if
     * there are cofactors to be found. However when we reach the latter
     * step, we're speaking of a fraction of a tiny fraction, and sieving
     * or resieving are not worth the trouble (we'll compute norms that
     * have numerous factors below lim, but we'll find them the hard way,
     * it's not that hard).
     *
     * When we have only one side, we no longer have this cumulative
     * "fraction of a tiny fraction" effect. BUT we still have some
     * interest in removing known primes from the norms, as we do in the
     * "normal" (= nowhere batch) case.
     */
    bool needs_resieving() const {
        for(auto const & s : sides)
            if (s.lim == 0)
                return false;
        return true;
    }

    template<sieve_method Algo>
    static void declare_usage(cxx_param_list & pl);
    template<sieve_method Algo>
    static bool parse_default(siever_config & sc, cxx_param_list & pl, int);

    /*{{{ has_same_config */
    bool operator==(siever_config const & o) const { return memcmp(this, &o, sizeof(*this)) == 0; }

    struct has_same_config {
        siever_config const & sc;
        has_same_config(siever_config const & sc) : sc(sc) {}
        bool operator()(siever_config const& o) const { return o == sc; }
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
    };
    has_same_config same_config() const { return has_same_config(*this); }
    /*}}}*/
#if 0
    /*{{{ has_same_config_q */
    struct has_same_config_q {
        siever_config const & sc;
        has_same_config_q(siever_config const & sc) : sc(sc) {}
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
        bool operator()(siever_config const& o) const {
            return sc.side == o.side && sc.bitsize == o.bitsize;
        }
    };
    has_same_config_q same_config_q() const {
        return has_same_config_q(*this);
    }
    /*}}}*/
#endif
    /* {{{ has_same_fb_parameters */
    struct has_same_fb_parameters {
        siever_config const & sc;
        has_same_fb_parameters(siever_config const & sc) : sc(sc) {}
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
        bool operator()(siever_config const& o) const {
            if (sc.sides.size() != o.sides.size())
                return false;
            bool ok = true;
            for(unsigned int side = 0 ; side < sc.sides.size(); side++) {
                ok = ok && sc.sides[side].lim == o.sides[side].lim;
                ok = ok && sc.sides[side].powlim == o.sides[side].powlim;
            }
            return ok;
        }
    };
    has_same_fb_parameters same_fb_parameters() const { return { *this }; }
    /*}}}*/
#if 0
    /*{{{ has_same_sieving -- currently duplicates has_same_fb_parameters */
    struct has_same_sieving {
        siever_config const & sc;
        has_same_sieving(siever_config const & sc) : sc(sc) {}
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
        bool operator()(siever_config const& o) const {
            return has_same_fb_parameters(sc)(o);
        }
    };
    has_same_sieving same_sieving() const { return has_same_sieving(*this); }
    /*}}}*/
#endif
    /*{{{ has_same_cofactoring */
    struct has_same_cofactoring {
        siever_config const & sc;
        has_same_cofactoring(siever_config const & sc) : sc(sc) {}
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
        bool operator()(siever_config const& o) const { 
            if (sc.sides.size() != o.sides.size())
                return false;
            bool ok = true;
            for(unsigned int side = 0 ; side < sc.sides.size(); side++) {
                ok = ok && sc.sides[side].lambda == o.sides[side].lambda;
                ok = ok && sc.sides[side].lpb == o.sides[side].lpb;
                ok = ok && sc.sides[side].mfb == o.sides[side].mfb;
                ok = ok && sc.sides[side].ncurves == o.sides[side].ncurves;
            }
            return ok;
        }
        struct comparison {
            bool operator()(siever_config const& a, siever_config const& b) const {
                unsigned int min = std::min(a.sides.size(), b.sides.size());
                for(unsigned int s = 0 ; s < min; s++) {
                    auto as = std::tie(a.sides[s].lambda, a.sides[s].lpb,
                                       a.sides[s].mfb, a.sides[s].ncurves);
                    auto bs = std::tie(b.sides[s].lambda, b.sides[s].lpb,
                                       b.sides[s].mfb, b.sides[s].ncurves);
                    if (as == bs)
                        continue;
                    return as < bs;
                }
                return a.sides.size() < b.sides.size();
            }
        };

    };
    has_same_cofactoring same_cofactoring() const { return { *this }; }
    /*}}}*/
};


/* {{{ descent_hint
 *
 * This is used for the descent. For each factor size, we provide a
 * reasonable siever_config value
 *
 * We also provide, based on experience, info relative to how long it
 * takes to finish the smoothing process for a prime factor of this size.
 */
/* }}} */

struct siever_config_pool {
    typedef std::pair<int, unsigned int> key_type;

    struct descent_hint : public siever_config {
        double expected_time;
        double expected_success;
    };

    typedef std::map<key_type, descent_hint> hint_table_t;
    hint_table_t hints;

    static constexpr int max_increase_lpb_default = 0;
    static constexpr int max_increase_logA_default = 4;

    int max_increase_lpb = max_increase_lpb_default;
    int max_increase_logA = max_increase_logA_default;

    int max_descent_attempts_allowed() const {
        return (base.adjust_strategy != 2) + max_increase_logA + max_increase_lpb;
    }

    descent_hint const * get_hint(int side, unsigned int bitsize) const {
        auto it = hints.find(key_type(side, bitsize));
        if (it == hints.end())
            return nullptr;
        else
            return &it->second;
    }

    /* The siever_config in [base] needs not be complete. The
     * default_config_ptr field points here if it is complete. If not,
     * the fields here are just used as a base for initializing the other
     * configurations.
     *
     * Note that while the "base for initializing" functionality is
     * likely to stay, the notion of a "default config" seems to be
     * screwed altogether, and we would rather like to see it disappear
     * someday. Currently it is used for displaying memory usage, setting
     * defaults for dupqmin, and getting the lim abd lpb parameters
     * before going to the batch step.
     */

    siever_config const * default_config_ptr = nullptr;
    siever_config base;

    siever_config get_config_for_q(special_q_task const & doing) const;

    template<sieve_method Algo>
    siever_config_pool(cxx_param_list& pl, int nb_polys, Algo);

    private:
    void parse_hints_file(const char * filename);
    public:

    double hint_expected_time(key_type const &K) const {
        if (hints.find(K) == hints.end())
            return -1;
        return hints.at(K).expected_time;
    }
    double hint_expected_success(key_type const &K) const {
        if (hints.find(K) == hints.end())
            return -1;
        return hints.at(K).expected_success;
    }

    static void declare_usage(cxx_param_list & pl);
};

#endif	/* CADO_LAS_SIEVER_CONFIG_HPP */
