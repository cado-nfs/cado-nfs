#ifndef LAS_SIEVER_CONFIG_HPP_
#define LAS_SIEVER_CONFIG_HPP_

#include "las-base.hpp"
#include "las-todo-entry.hpp"
#include "fb-types.h"
#include <string.h>
#include <string>
#include <map>
#include "params.h"

/* siever_config */
 
/* The following structure lists the fields with an impact on the siever.
 * Different values for these fields will correspond to different siever
 * structures.
 */
struct siever_config : public _padded_pod<siever_config> {
    /* The bit size of the special-q. Counting in bits is no necessity,
     * we could imagine being more accurate */
    unsigned int bitsize;  /* bitsize == 0 indicates end of table */
    int side;              /* special-q side */
    int logA;


    /* For a given logA, we may trigger configurations for various logI
     * values. logI_adjusted is a sieving-only parameter. */
    int logI_adjusted;
    sublat_t sublat;
    unsigned long bucket_thresh;    // bucket sieve primes >= bucket_thresh
    unsigned long bucket_thresh1;   // primes above are 2-level bucket-sieved
    unsigned int td_thresh;
    unsigned int skipped;           // don't sieve below this
    unsigned int unsieve_thresh;
    struct side_config {
        unsigned long lim;    /* factor base bound */
        unsigned long powlim; /* bound on powers in the factor base */
        int lpb;              /* large prime bound is 2^lpb */
        int mfb;              /* bound for residuals is 2^mfb */
        int ncurves;          /* number of cofactorization curves */
        double lambda;        /* lambda sieve parameter */
        unsigned long qmin;   /* smallest q sieved on this side, for dupsup */
        unsigned long qmax;   /* largest q sieved on this side, for dupsup */
    };
    side_config sides[2];

    void display() const;

    static void declare_usage(param_list_ptr pl);
    static bool parse_default(siever_config & sc, param_list_ptr pl);

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
    /* {{{ has_same_fb_parameters */
    struct has_same_fb_parameters {
        siever_config const & sc;
        has_same_fb_parameters(siever_config const & sc) : sc(sc) {}
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
        bool operator()(siever_config const& o) const {
            bool ok = true;
            // ok = ok && sc.logI_adjusted == o.logI_adjusted;
            ok = ok && sc.bucket_thresh == o.bucket_thresh;
            ok = ok && sc.bucket_thresh1 == o.bucket_thresh1;
            ok = ok && sc.td_thresh == o.td_thresh;
            ok = ok && sc.skipped == o.skipped;
            // ok = ok && sc.bk_multiplier == o.bk_multiplier;
            ok = ok && sc.unsieve_thresh == o.unsieve_thresh;
            for(int side = 0 ; side < 2 ; side++) {
                ok = ok && sc.sides[side].lim == o.sides[side].lim;
                ok = ok && sc.sides[side].powlim == o.sides[side].powlim;
            }
            return ok;
        }
    };
    has_same_fb_parameters same_fb_parameters() const { return has_same_fb_parameters(*this); }
    /*}}}*/
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
    /*{{{ has_same_cofactoring */
    struct has_same_cofactoring {
        siever_config const & sc;
        has_same_cofactoring(siever_config const & sc) : sc(sc) {}
        template<typename OtherLasType>
        bool operator()(OtherLasType const& o) const { return (*this)(o.conf); }
        bool operator()(siever_config const& o) const { 
            bool ok = true;
            for(int side = 0 ; side < 2 ; side++) {
                ok = ok && sc.sides[side].lambda == o.sides[side].lambda;
                ok = ok && sc.sides[side].lpb == o.sides[side].lpb;
                ok = ok && sc.sides[side].mfb == o.sides[side].mfb;
                ok = ok && sc.sides[side].ncurves == o.sides[side].ncurves;
            }
            return ok;
        }
    };
    has_same_cofactoring same_cofactoring() const { return has_same_cofactoring(*this); }
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

    descent_hint const * get_hint(int side, unsigned int bitsize) const {
        hint_table_t::const_iterator it = hints.find(key_type(side, bitsize));
        if (it == hints.end())
            return NULL;
        else
            return &it->second;
    }

    siever_config const * default_config_ptr;

    /* This needs not be complete. The default_config_ptr field points
     * here if it is complete. If not, the fields here are just used as a
     * base for initializing the other configurations */
    siever_config base;

    siever_config get_config_for_q(las_todo_entry const& doing) const;

    siever_config_pool(cxx_param_list& pl);

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
};

#endif	/* LAS_SIEVER_CONFIG_HPP_ */
