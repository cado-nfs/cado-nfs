/*****************************************************************
*                Functions for the factor base                  *
*****************************************************************/

#ifndef FB_HPP_
#define FB_HPP_

#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>
#include <gmp.h>
#include <array>
#include <vector>
#include <deque>
#include <list>
#include <algorithm>
#include <map>
#include "multityped_array.hpp"
#include "las-config.h"
#include "cado_poly.h" // for MAX_DEGREE
#include "fb-types.h"
#include "las-qlattice.hpp"
#include "las-plattice.hpp"
#include "las-base.hpp"

#ifdef HAVE_GLIBC_VECTOR_INTERNALS
#include "mmappable_vector.hpp"
#else
/* yes, it's a hack */
#define mmappable_vector std::vector
#endif

/* The *factor base* is made of *entries*. We have several vectors of
 * entries, each with primes splitting in the same number of roots.
 *
 * For each set of *thresholds* and each *scale* value, we define a
 * *slicing*. A slicing is a division of the factor base into *parts*.
 * Then each part contains several vectors of *slices* (in the same
 * manner as we have several vectors of entries in the factor base). A
 * slice is formed by basically two pointers into the vector of entries
 * in the factor base that have this number of roots.
 *
 * The factor base also contains auxiliary data attached to each vector
 * of entries, namely the cumulative distribution function of the weight.
 * This is used to subdivide slices into pieces of equal weight when
 * needed.
 */

/* {{{ fb entries: sets of prime ideals above one single rational prime p. */
/* Forward declaration so fb_entry_general can use it in constructors */
template <int Nr_roots>
class fb_entry_x_roots;

/* A root modulo a prime power q. q is specified externally */
struct fb_general_root {
  /* exp and oldexp are maximal such that:
     If not projective and a == br (mod p^k), then p^exp | F(a,b)
     -"-                   a == br (mod p^(k-1)), then p^oldexp | F(a,b)
     If projective and ar == b  -"- */
  fbroot_t r;
  bool proj;
  unsigned char exp, oldexp;

  fb_general_root (){}
  fb_general_root (const fbroot_t r, const unsigned char nexp=1,
                   const unsigned char oldexp=0, const bool proj=false) :
                   r(r), proj(proj), exp(nexp), oldexp(oldexp) {}
  /* Create a root from a linear polynomial */
  fb_general_root (fbprime_t q, cxx_mpz_poly const & poly, const unsigned char nexp=1,
                   const unsigned char oldexp=0);

  /* Constructor from the old format of storing projective roots, which has q
     added to the root if the root is projective */
  fb_general_root (const unsigned long long old_r, const fbprime_t q,
                   const unsigned char nexp=1, const unsigned char oldexp=0) :
                   r((old_r >= q) ? (old_r - q) : old_r),
                   proj(old_r >= q), exp(nexp), oldexp(oldexp) {}

  /* A root is simple if it is not projective and the exp goes from 0 to 1 */
  bool is_simple() const {return exp == 1 && oldexp == 0 && !proj;}

  /* Convert a root to the old format of storing projective roots with q added */
  unsigned long long to_old_format(const fbprime_t q) const {
    return (unsigned long long) r + (proj ? q : 0);
  }

  /* Print one root. Projective roots are printed as r+q */
  void fprint(FILE *out, const fbprime_t q) const {
    fprintf(out, "%llu", to_old_format(q));
    if (oldexp != 0 || this->exp != 1)
      fprintf(out, ":%hhu:%hhu", oldexp, this->exp);
  }

  void transform(fb_general_root &result, const fbprime_t q,
                 const redc_invp_t invq,
                 const qlattice_basis &basis) const {
    unsigned long long t = to_old_format(q);
    t = fb_root_in_qlattice(q, t, invq, basis);
    result = fb_general_root(t, q, exp, oldexp);
  }
};

/* General entries are anything that needs auxiliary information:
   Prime powers, projective roots, ramified primes where exp != oldexp + 1,
   etc. They could, of course, also store the simple cases, but for those we
   use the simple struct to conserve memory and to decide algorithms (batch
   inversion, etc.) statically. */
class fb_entry_general : public _padded_pod<fb_entry_general> {
  void read_roots (const char *, unsigned char, unsigned char, unsigned long);
public:
  typedef fb_entry_general transformed_entry_t;
  fbprime_t q, p; /* q = p^k */
  redc_invp_t invq; /* invq = -1/q (mod 2^32), or (mod 2^64), depending on
		       the size of redc_invp_t */
  fb_general_root roots[MAX_DEGREE];
  unsigned char k, nr_roots;
  /* Static class members to allow fb_vector<> to distinguish between and
     operate on both kind of entries */
  static const bool is_general_type = true;
  static const unsigned char fixed_nr_roots = 0;
  inline int get_nr_roots() const { return nr_roots; }

  fb_entry_general() {}
  template <int Nr_roots>
  fb_entry_general (const fb_entry_x_roots<Nr_roots> &e);
  fbprime_t get_q() const {return q;}
  fbroot_t get_r(const size_t i) const {return roots[i].r;};
  fbroot_t get_proj(const size_t i) const {return roots[i].proj;};
  void parse_line (const char *line, unsigned long linenr);
  void merge (const fb_entry_general &);
  void fprint(FILE *out) const;
  bool is_simple() const;
  void transform_roots(transformed_entry_t &, const qlattice_basis &) const;
  double weight() const {return 1./q * nr_roots;}
  /* Allow sorting by p */
  bool operator<(const fb_entry_general &other) const {return this->p < other.p;}
  bool operator>(const fb_entry_general &other) const {return this->p > other.p;}
  struct sort_byq {
    bool operator()(fb_entry_general const & a, fb_entry_general const & b) const { return a.get_q() < b.get_q(); };
  };
};

template <int Nr_roots>
class fb_transformed_entry_x_roots {
public:
  fbprime_t p;
  fbroot_t roots[Nr_roots];
  bool proj[Nr_roots];
  static const unsigned char k = 1, nr_roots = Nr_roots;
  /* Static class members to allow fb_vector<> to distinguish between and
     operate on both kind of entries */
  static const bool is_general_type = false;
  static const unsigned char fixed_nr_roots = Nr_roots;
  fbprime_t get_q() const {return p;}
  fbroot_t get_r(const size_t i) const {return roots[i];};
  fbroot_t get_proj(const size_t i) const {return proj[i];};
};

/* "Simple" factor base entries. We imply q=p, k=1, oldexp=0, exp=1,
   and projective=false for all roots. */
template <int Nr_roots>
class fb_entry_x_roots : public _padded_pod<fb_entry_x_roots<Nr_roots> > {
public:
  typedef fb_transformed_entry_x_roots<Nr_roots> transformed_entry_t;
  fbprime_t p;
  redc_invp_t invq; /* invq = -1/q (mod 2^32), or (mod 2^64), depending on
		       the size of redc_invp_t */
  fbroot_t roots[Nr_roots];
  /* Static class members to allow fb_vector<> to distinguish between and
     operate on both kind of entries */
  static const unsigned char k = 1, nr_roots = Nr_roots;
  static const bool is_general_type = false;
  static const unsigned char fixed_nr_roots = Nr_roots;
  inline int get_nr_roots() const { return Nr_roots; }
  fb_entry_x_roots() {};
  /* Allow assignment-construction from general entries */
  fb_entry_x_roots(const fb_entry_general &e) : p(e.p), invq(e.invq) {
    ASSERT_ALWAYS(Nr_roots == e.nr_roots);
    ASSERT(e.is_simple());
    for (int i = 0; i < Nr_roots; i++)
      roots[i] = e.roots[i].r;
  }
  fbprime_t get_q() const {return p;}
  fbroot_t get_r(const size_t i) const {return roots[i];};
  fbroot_t get_proj(const size_t i MAYBE_UNUSED) const {return false;};
  double weight() const {return 1./p * Nr_roots;}
  /* Allow sorting by p */
  bool operator<(const fb_entry_x_roots<Nr_roots> &other) const {return this->p < other.p;}
  bool operator>(const fb_entry_x_roots<Nr_roots> &other) const {return this->p > other.p;}
  void fprint(FILE *) const;
  void transform_roots(transformed_entry_t &, const qlattice_basis &) const;
};

/* }}} */

class fb_slice_interface {
    public:
        fb_slice_interface() = default;
        virtual ~fb_slice_interface(){}
        virtual size_t size() const = 0;
        virtual unsigned char get_logp() const = 0;
        virtual fbprime_t get_prime(slice_offset_t offset) const = 0;
        virtual unsigned char get_k(slice_offset_t offset) const = 0;
        virtual slice_index_t get_index() const = 0;
        virtual double get_weight() const = 0;
        virtual bool is_general() const = 0;
        virtual int get_nr_roots() const = 0;
};

/* This is one of the function objects in fb.cpp that needs to tinker
 * with the internals of fb_slice.
 */
struct helper_functor_subdivide_slices;

/* see fb_slice_weight.hpp */
template<typename FB_ENTRY_TYPE>
class fb_slice_weight_estimator;

template<typename FB_ENTRY_TYPE>
class fb_slice : public fb_slice_interface {
    friend class fb_slice_weight_estimator<FB_ENTRY_TYPE>;
    typedef mmappable_vector<FB_ENTRY_TYPE> fb_entry_vector;
    typename fb_entry_vector::const_iterator _begin, _end;
    unsigned char logp;
    slice_index_t index;
    double weight;
    friend struct helper_functor_subdivide_slices;
    fb_slice(typename fb_entry_vector::const_iterator it, unsigned char logp) : _begin(it), _end(it), logp(logp), index(0), weight(0) {}
    fb_slice(typename fb_entry_vector::const_iterator it, typename fb_entry_vector::const_iterator jt, unsigned char logp) : _begin(it), _end(jt), logp(logp), index(0), weight(0) {}
    public:
    typedef FB_ENTRY_TYPE entry_t;
    inline typename fb_entry_vector::const_iterator begin() const { return _begin; }
    inline typename fb_entry_vector::const_iterator end() const { return _end; }
    inline size_t size() const override { return _end - _begin; }
    unsigned char get_logp() const override { return logp; }
    fbprime_t get_prime(slice_offset_t offset) const override {
        /* While it may make sense to manipulate slices that exceed the
         * expected max size during the work to construct them, it should
         * be pretty clear that we'll never ever try to use get_prime,
         * which is limited in its type, on a slice that has more prime
         * ideals that this function can address */
        ASSERT(size() <= std::numeric_limits<slice_offset_t>::max());
        return _begin[offset].p;
    }
    unsigned char get_k(slice_offset_t offset) const override {
        return _begin[offset].k;
        /* power. Well, most often it's a constant ! We need to
         * access it from the virtual base though. This is way we're
         * not folding it to a template access.  */
    }
    slice_index_t get_index() const override {return index;}
    double get_weight() const override {return weight;}
    bool is_general() const override { return FB_ENTRY_TYPE::is_general_type; }
    /* get_nr_roots() on a fb_slice returns zero for slices of general
     * type ! */
    int get_nr_roots() const override { return FB_ENTRY_TYPE::fixed_nr_roots;}
};

/* entries and general_entries: we declare general entries, as well
 * as entries for all numbers of roots between 1 and MAX_ROOTS
 *
 * (notationally, general_entries corresponds to -1 as a number of
 * roots).
 * */

template<typename T> struct entries_and_cdf {
    typedef mmappable_vector<T> container_type;
    typedef mmappable_vector<double> weight_container_type;
    struct type : public container_type {
        typedef typename entries_and_cdf<T>::container_type container_type;
        typedef typename entries_and_cdf<T>::weight_container_type weight_container_type;
        /* cumulative distribution function. This is set up by
         * helper_functor_append. We use it to split into slices.
         * weight_cdf[i] is \sum_{j < i} super[j].weight
         *
         * Thus we always need a lone 0 to initialize. See
         * helper_functor_put_first_0, used in the fb_factorbase ctor.
         */
        /* We *might* want to consider the cdf only for several entries
         * at a time (say, 16), so as to minimize the cost of finding the
         * split points */
        weight_container_type weight_cdf;
        inline weight_container_type::const_iterator weight_begin() const { 
            return weight_cdf.begin();
        }
        inline weight_container_type::const_iterator weight_end() const { 
            return weight_cdf.end();
        }
        inline weight_container_type::size_type weight_size() const { 
            return weight_cdf.size();
        }
        double weight_delta(size_t a, size_t b) const {
            return weight_cdf[b] - weight_cdf[a];
        }
        double weight_delta(typename container_type::const_iterator a, typename container_type::const_iterator b) const {
            return weight_cdf[b-container_type::begin()] - weight_cdf[a-container_type::begin()];
        }
        double weight_cdf_at(size_t a) const {
            return weight_cdf[a];
        }
        double weight_cdf_at(typename container_type::const_iterator a) const {
            return weight_cdf[a-container_type::begin()];
        }
    };
};
template<int n> struct fb_entries_factory {
    typedef typename entries_and_cdf<fb_entry_x_roots<n>>::type type;
};
template<> struct fb_entries_factory<-1> {
    typedef typename entries_and_cdf<fb_entry_general>::type type;
};
template<int n> struct fb_slices_factory {
    typedef std::vector<fb_slice<fb_entry_x_roots<n>>> type;
};
template<> struct fb_slices_factory<-1> {
    typedef std::vector<fb_slice<fb_entry_general>> type;
};

class fb_factorbase {
    public:
        /* Has to be <= MAX_DEGREE */
    static const int MAX_ROOTS = MAX_DEGREE;

    private:

    cxx_mpz_poly f;
    int side;
    unsigned long lim;
    unsigned long powlim;

    typedef multityped_array<fb_entries_factory, -1, MAX_ROOTS+1> entries_t;
    entries_t entries;

    public:
    typedef std::array<size_t, MAX_ROOTS+2> threshold_pos;
    threshold_pos get_threshold_pos(fbprime_t) const;

    private:
    /* The threshold position cache is used to accelerate the creation of
     * slicings. */
    mutable std::map<fbprime_t, threshold_pos> threshold_pos_cache;

    /* {{{ append. This inserts a pool of entries to the factor base. We
     * do it with many entries in a row so as to avoid looping MAX_ROOTS
     * times for each root, and so that we don't feel sorry to use
     * multityped_array_foreach (even though in truth, it's quite
     * probably fine)
     */
    private:
    struct helper_functor_append;
    public:
    void append(std::list<fb_entry_general> &pool);
    /* }}} */

    /* {{{ Various ways to count primes, prime ideals, and hit ratio
     * ("weight") of entries in the whole factor base, or in subranges
     */
    private:
    struct helper_functor_count_primes;
    struct helper_functor_count_prime_ideals;
    struct helper_functor_count_weight;
    struct helper_functor_count_combined;
    struct helper_functor_count_primes_interval;
    struct helper_functor_count_prime_ideals_interval;
    struct helper_functor_count_weight_interval;
    struct helper_functor_count_combined_interval;
    public:
    size_t count_primes() const;
    size_t count_prime_ideals() const;
    size_t count_weight() const;
    /* }}} */


    /* the key type lists the things with respect to which we're willing
     * to divide the views on our factor base.
     */
        struct key_type {
            fbprime_t thresholds[FB_MAX_PARTS];
            fbprime_t td_thresh;
            fbprime_t skipped;
            double scale;
            /* This might seem non obvious, but this parameters controls
             * the size of the slices, because we want to enforce some
             * relatively-even division. It's not entirely clear that we
             * want it here, but we definitely want it somewhere. */
            unsigned int nb_threads;

            bool operator<(key_type const& x) const {
                int r;
                for(int i = 0; i < FB_MAX_PARTS; ++i) {
                    r = (thresholds[i] > x.thresholds[i]) - (x.thresholds[i] > thresholds[i]);
                    if (r) return r < 0;
                }
                r = (td_thresh > x.td_thresh) - (x.td_thresh > td_thresh);
                if (r) return r < 0;
                r = (skipped > x.skipped) - (x.skipped > skipped);
                if (r) return r < 0;
                return scale < x.scale;
            }
        };

    public:
    class slicing {
        public:
        struct stats_type {
            /* explicit-initializing as below forces zero-initialization
             * of members */
            std::array<size_t, FB_MAX_PARTS> primes {{}};
            std::array<size_t, FB_MAX_PARTS> ideals {{}};
            std::array<double, FB_MAX_PARTS> weight {{}};
        };
        stats_type stats;

        class part {
            /* slices and general_slices: actually general_slices is
             * slices for number of roots -1, and that is 
             * always a vector of length 1. And we have slices (vectors
             * of fb_slice objects) for all numbers of roots between 1 and
             * MAX_ROOTS.
             */
            typedef multityped_array<fb_slices_factory, -1, MAX_ROOTS+1> slices_t;
            friend struct helper_functor_subdivide_slices;
            slices_t slices;
            slice_index_t first_slice_index = 0;
            public:
            template<int n>
                typename fb_slices_factory<n>::type& get_slices_vector_for_nroots() {
                    return slices.get<n>();
                }
            private:
            
            /* the general vector (the one with index -1 in the
             * multityped array) is made of "anything that is not
             * simple" (as per the logic that used to be in
             * fb_part::append and fb_entry_general::is_simple). That means:
             *  - all powers go to the general vector
             *  - everything below bkthresh, too (but see below).
             *  - primes with multiple roots as well (but to be honest,
             *    it's not entirely clear to me why -- after all, there's
             *    no special treatment to speak of).
             */

            friend class slicing;

            /* We use this to access our arrays called slices and
             * general_slices.
             *
             * 0 <= i <slicesG.size()   -->  return &(slicesG[i])
             * slicesG.size() <= i < slices0.size() --> return &(slices0[i-slicesG.size()])
             * slices0.size() <= i < slices1.size() --> return &(slices1[i-slicesG.size()-slices0.size()])
             * and so on.
             */
            struct helper_functor_get {
                typedef fb_slice_interface const * type;
                typedef slice_index_t key_type;
                template<typename T>
                    type operator()(T const & x, slice_index_t & k) {
                        if ((size_t) k < x.size()) {
                            return &(x[k]);
                        } else {
                            k -= x.size();
                            return NULL;
                        }
                    }
            };

            fb_slice_interface const * get(slice_index_t & index) const {
                /* used to be in fb_vector::get_slice
                 *
                 * and in fb_part::get_slice for the lookup of the
                 * vector.
                 *
                 * TODO: dichotomy, perhaps ? Can we do the dichotomy
                 * elegantly ?
                 */
                if (index < first_slice_index) return NULL;
                slice_index_t idx = index - first_slice_index;
                if (idx >= nslices()) {
                    index = idx - nslices();
                    return NULL;
                }
                const fb_slice_interface * res = multityped_array_locate<helper_functor_get>()(slices, idx);
                ASSERT_ALWAYS(res);
                ASSERT_ALWAYS(res->get_index() == index);
                if (res) return res;
                return NULL;
            }

            /* {{{ use caching for the number of slices, as it's a handy
             * thing to query */
            mutable slice_index_t _nslices = std::numeric_limits<slice_index_t>::max();
            struct helper_functor_nslices {
                template<typename T>
                    slice_index_t operator()(slice_index_t t, T const & x) const {
                        return t + x.size();
                    }
            };
            slice_index_t nslices() const {
                if (_nslices != std::numeric_limits<slice_index_t>::max())
                    return _nslices;
                helper_functor_nslices N;
                return _nslices = multityped_array_fold(N, 0, slices);
            }
            /* }}} */

            template<typename F>
            struct foreach_slice_s {
                F & f;
                template<typename T>
                void operator()(T & x) { for(auto & a : x) f(a); }
                template<typename T>
                void operator()(T const & x) { for(auto const & a : x) f(a); }
            };
            public:
            template<typename F>
            void foreach_slice(F & f) {
                foreach_slice_s<F> FF { f };
                multityped_array_foreach(FF, slices);
            }
            template<typename F>
            void foreach_slice(F & f) const {
                foreach_slice_s<F> FF { f };
                multityped_array_foreach(FF, slices);
            }
            /*
             * old g++ seems to have difficulties with this variant, and
             * is puzzled by the apparent ambiguity -- newer g++ groks it
             * correctly, as does clang
            template<typename F>
            void foreach_slice(F && f) {
                multityped_array_foreach(foreach_slice_s<F> { f }, slices);
            }
            template<typename F>
            void foreach_slice(F && f) const {
                multityped_array_foreach(foreach_slice_s<F> { f }, slices);
            }
            */

            fb_slice_interface const & operator[](slice_index_t index) const {
                /* This bombs out at runtime if get returns NULL, but
                 * then it should be an indication of a programmer
                 * mistake.
                 */
                return *get(index);
            }
        };

        private:

        /* We have "parts" that correspond to the various layers of the
         * bucket sieving
         *
         *  ==> THERE IS NO "part" FOR THE SMALL SIEVE <==
         *
         * This means that parts[0] will always be a meaningless
         * placeholder. Indeed, the small sieve does not care about
         * slices!
         */
        std::array<part, FB_MAX_PARTS> parts;

        /* toplevel is set by the ctor */
        int toplevel = 0;

        fb_slice_interface const * get(slice_index_t index) const {
            for (auto const & p : parts) {
                if (index < p.first_slice_index + p.nslices()) 
                    return p.get(index);
            }
            return NULL;
        }
        public:

        part const & get_part(int i) const { return parts[i]; }

        inline int get_toplevel() const { return toplevel; }

        /* This accesses the *fb_slice* with this index. Not the vector of
         * slices ! */
        fb_slice_interface const & operator[](slice_index_t index) const {
            return *get(index);
        }

        public:

        /* Here, when computing the slices, we stow away the small
         * primes, and arrange them according to the internal logic of
         * the small sieve. In particular, we want to keep track of where
         * the resieved primes are.
         *
         * Note that the data that we prepare here is still not attached
         * to any special-q. We're replacing what used to be done in
         * init_fb_smallsieved, but not in small_sieve_init.
         *
         * Primes below the bucket-sieve threshold are small-sieved. Some
         * will be treated specially when cofactoring:
         *
         *  - primes such that nb_roots / p >= 1/tdhresh will be
         *    trial-divided.
         *  - primes (ideals) that are not trial divided and that are not
         *    powers are resieved.
         *
         * And of course, depending on the special-q, small-sieved prime
         * ideals become projective, and therefore escape the general
         * treatment.
         *
         * Note that the condition on trial division is not clear-cut
         * w.r.t p. The small sieve code is somewhat dependent on the
         * number of hits per row, which is I/p. And it matters to have
         * small-sieved primes sorted according to this value. Hence, as
         * p increases, we might have:
         *  small p's that are trial-divided because 1/p >=
         *  1/td_thresh.
         *  some p's with p<=2*td_thresh, among which some are trial
         *   divided if 2 roots or more are above, otherwise they're
         *   resieved.
         *  some p's with p<=3*td_thresh, again with a distinction...
         *  and so on up to p>d*tdshresh, (d being the polynomial
         *  degree). Above this bound we're sure that we no longer see
         *  any trial-divided prime.
         *
         * For each slicing, we elect to compute two copies of the lists
         * of prime ideals below bkthresh:
         *  - first, prime ideals that we know will be resieved. Those
         *    form the largest part of the list.
         *  - second, prime ideals above primes that will be trial-divided.
         *
         */

        /* This contains all factor base prime ideals below the bucket
         * threshold.  */
        struct small_sieve_entries_t {
            /* general ones. Not powers, not trial-divided ones. */
            std::vector<fb_entry_general> resieved;
            /* the rest. some are powers, others are trial-divided */
            std::vector<fb_entry_general> rest;
            /* from "rest" above, we can infer the list of trial-divided
             * primes by merely restricting to entries with k==1 */

            /* this is sort of a trash can. We won't use these primes for
             * the small sieve, but they do matter for trial division */
            std::vector<unsigned long> skipped;
        };
        small_sieve_entries_t small_sieve_entries;
        /* TODO: I doubt that the struct above will stay. Seems awkward.
         * We'd better have a single vector, and the position of the
         * "end-of-resieved" thing.
         */

        template<typename T>
        void foreach_slice(T & f) {
            for (auto & p : parts) {
                p.foreach_slice(f);
            }
        }
        template<typename T>
        void foreach_slice(T && f) {
            for (auto & p : parts) {
                p.foreach_slice(f);
            }
        }

        slicing() = default;
        slicing(fb_factorbase const & fb, key_type const & K);
    };

    private:
        std::map<key_type, slicing> cache;
        int read(const char * const filename);

    public:
        /* accessors.
         * As in the std::map case, the non-const accessor is allowed to
         * create stuff. */
    
        inline slicing & operator[](key_type const& K) { 
            auto it = cache.find(K);
            if (it != cache.end())
                return it->second;

            /* create a new slot */
            return cache[K] = slicing(*this, K);
        }

        inline slicing const & operator[](key_type const& K) const {
            return cache.at(K);
        }

    private:
        void make_linear();
        void make_linear_threadpool (unsigned int nb_threads);

    public:
        fb_factorbase(cxx_cado_poly const & cpoly, int side, unsigned long lim, unsigned long powlim, cxx_param_list & pl, const char * fbc_filename);

    private:
        struct sorter {
            template<typename T>
            void operator()(T & x) {
                /* not entirely clear to me. Shall we sort by q or by p ?
                 */
                typedef typename T::value_type X;
                auto by_q = [](X const & a, X const & b) { return a.get_q() < b.get_q(); };
                /*
                auto by_p_then_q = [](X const & a, X const & b) {
                    return
                        a.get_p() < b.get_p() ||
                        a.get_p() == b.get_p() &&
                        a.get_q() < b.get_q();
                };
                */
                std::sort(x.begin(), x.end(), by_q);
            }
        };

    public:
        void finalize() {
            sorter S;
            multityped_array_foreach(S, entries);
        }
};

std::ostream& operator<<(std::ostream& o, fb_factorbase::key_type const &);

unsigned char   fb_log (double x, double y, double z);
unsigned char	fb_log_delta (fbprime_t, unsigned long, unsigned long, double);
fbprime_t       fb_is_power (fbprime_t, unsigned long *);

/* This is defined in fb-slice-weight.cpp */
void print_slice_weight_estimator_stats();

#endif  /* FB_HPP_ */
