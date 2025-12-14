#ifndef CADO_FB_HPP
#define CADO_FB_HPP

#include <cstddef>
#include <cstdio>

#include <algorithm>
#include <array>
#include <iosfwd>
#include <limits>
#include <list>
#include <map>
#include <mutex>
#include <type_traits>
#include <utility>
#include <vector>

#include "cado_poly.h"
#include "macros.h"
#include "mpz_poly.h"

#include "fb-types.hpp"
#include "las-config.hpp"
#include "lock_guarded_container.hpp"
#include "mmappable_vector.hpp"
#include "multityped_array.hpp"

struct qlattice_basis;
struct cxx_param_list;

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
template <int Nr_roots> class fb_entry_x_roots;

/* A root modulo a prime power q=p^k. q is specified externally */
struct fb_general_root {
    /* exp and oldexp are maximal such that:
       If not projective and a == br (mod p^k), then p^exp | F(a,b)
       -"-                   a == br (mod p^(k-1)), then p^oldexp | F(a,b)
       If projective and ar == b  -"- */

    /* Note that we must arrange so that this structure has no hidden
     * padding, and has all its fields default-initialized to zero, so that
     * we can write to the factor base cache file without transferring
     * potentially uninitialized bytes.
     */
    fbroot_t r = 0;
    bool proj = false;
    unsigned char exp = 0, oldexp = 0;

  private:
    unsigned char dummy_padding_byte MAYBE_UNUSED_PRIVATE_DATA_MEMBER = 0;

  public:
    fb_general_root() = default;
    fb_general_root(fbroot_t const r, unsigned char const nexp,
                    unsigned char const oldexp, bool const proj = false)
        : r(r)
        , proj(proj)
        , exp(nexp)
        , oldexp(oldexp)
    {
    }
    explicit fb_general_root(fbroot_t r)
        : r(r)
        , exp(1)
    {
    }
    /* Create a root from a linear polynomial */
    static fb_general_root fb_linear_root(fbprime_t q,
                                          cxx_mpz_poly const & poly,
                                          unsigned char nexp,
                                          unsigned char oldexp);

  private:
    friend class fb_entry_general;
    /* Constructor from the old format of storing projective roots, which has q
       added to the root if the root is projective */
    explicit fb_general_root(fb_root_p1 const R, unsigned char const nexp = 1,
                    unsigned char const oldexp = 0)
        : r(R.r)
        , proj(R.proj)
        , exp(nexp)
        , oldexp(oldexp)
    {
    }

    /* A root is simple if it is not projective and the exp goes from 0 to 1 */
    bool is_simple() const { return exp == 1 && oldexp == 0 && !proj; }

    /* Convert a root to the old format of storing projective roots with q added
     */
    unsigned long long to_old_format(fbprime_t const q) const
    {
        return (unsigned long long)r + (proj ? q : 0);
    }

    /* Print one root. Projective roots are printed as r+q */
    void fprint(FILE * out, fbprime_t const q) const
    {
        fprintf(out, "%llu", to_old_format(q));
        if (oldexp != 0 || exp != 1)
            fprintf(out, ":%hhu:%hhu", oldexp, exp);
    }

    void transform(fb_general_root & result, fbprime_t q,
                   redc_invp_t invq, qlattice_basis const & basis) const;
};
/* std::has_unique_object_representations should be almost what we want,
 * but unfortunately it does not seem to work as it should (tested on
 * g++13 and clang-{14,15,16})
 */
static_assert(std::has_unique_object_representations_v<fb_general_root>,
              "fb_general_root must have no padding");
static_assert(sizeof(fb_general_root) == 8,
              "fb_general_root must have no padding");

/* General entries are anything that needs auxiliary information:
   Prime powers, projective roots, ramified primes where exp != oldexp + 1,
   etc. They could, of course, also store the simple cases, but for those we
   use the simple struct to conserve memory and to decide algorithms (batch
   inversion, etc.) statically. */
class fb_entry_general
{
    void read_roots(char const *, unsigned char, unsigned char, unsigned long);

  public:
    using transformed_entry_t = fb_entry_general;
    fbprime_t q = 0, p = 0;   /* q = p^k */
    redc_invp_t invq = 0; /* invq = -1/q (mod 2^32), or (mod 2^64), depending on
                         the size of redc_invp_t */
    unsigned char k = 0, nr_roots = 0;

  private:
    unsigned char dummy_padding_byte1 MAYBE_UNUSED_PRIVATE_DATA_MEMBER = 0;
    unsigned char dummy_padding_byte2 MAYBE_UNUSED_PRIVATE_DATA_MEMBER = 0;

  public:
    fb_general_root roots[MAX_DEGREE];
    /* Static class members to allow fb_vector<> to distinguish between and
       operate on both kind of entries */
    static bool const is_general_type = true;
    static unsigned char const fixed_nr_roots = 0;
    int get_nr_roots() const { return nr_roots; }

    fb_entry_general() = default;
    template <int Nr_roots>
    explicit fb_entry_general(fb_entry_x_roots<Nr_roots> const & e);
    fbprime_t get_q() const { return q; }
    fbroot_t get_r(size_t const i) const { return roots[i].r; };
    fbroot_t get_proj(size_t const i) const { return roots[i].proj; };
    void parse_line(char const * line, unsigned long linenr);
    bool can_merge(fb_entry_general const &) const;
    void merge(fb_entry_general const &);
    void fprint(FILE * out) const;
    bool is_simple() const;
    void transform_roots(transformed_entry_t &, qlattice_basis const &) const;
    double weight() const { return 1. / q * nr_roots; }
    /* Allow sorting by p */
    bool operator<(fb_entry_general const & other) const
    {
        return this->p < other.p;
    }
    bool operator>(fb_entry_general const & other) const
    {
        return this->p > other.p;
    }
    struct sort_byq {
        bool operator()(fb_entry_general const & a,
                        fb_entry_general const & b) const
        {
            return a.get_q() < b.get_q();
        };
    };
};
/* std::has_unique_object_representations should be almost what we want,
 * but unfortunately it does not seem to work as it should (tested on
 * g++13 and clang-{14,15,16})
 */
static_assert(std::has_unique_object_representations_v<fb_entry_general>,
              "fb_entry_general must have no padding");
static_assert(sizeof(fb_entry_general) == 4 * 4 + MAX_DEGREE * 8,
              "fb_entry_general must have no padding");

template <int Nr_roots> class fb_transformed_entry_x_roots
{
  public:
    fbprime_t p;
    std::array<fbroot_t, Nr_roots> roots;
    std::array<bool, Nr_roots> proj;
    static unsigned char const k = 1, nr_roots = Nr_roots;
    /* Static class members to allow fb_vector<> to distinguish between and
       operate on both kind of entries */
    static bool const is_general_type = false;
    static unsigned char const fixed_nr_roots = Nr_roots;
    fbprime_t get_q() const { return p; }
    fbroot_t get_r(size_t const i) const { return roots[i]; };
    fbroot_t get_proj(size_t const i) const { return proj[i]; };
};

/* "Simple" factor base entries. We imply q=p, k=1, oldexp=0, exp=1,
   and projective=false for all roots. */
template <int Nr_roots> class fb_entry_x_roots
{
  public:
    using transformed_entry_t = fb_transformed_entry_x_roots<Nr_roots>;
    fbprime_t p;
    redc_invp_t invq; /* invq = -1/q (mod 2^32), or (mod 2^64), depending on
                         the size of redc_invp_t */
    std::array<fbroot_t, Nr_roots> roots;
    /* Static class members to allow fb_vector<> to distinguish between and
       operate on both kind of entries */
    static unsigned char const k = 1, nr_roots = Nr_roots;
    static bool const is_general_type = false;
    static unsigned char const fixed_nr_roots = Nr_roots;
    int get_nr_roots() const { return Nr_roots; }
    // fb_entry_x_roots() {};
    fb_entry_x_roots(fbprime_t p, redc_invp_t invq, const fbroot_t * roots)
        : p(p)
        , invq(invq)
    {
        std::copy_n(roots, Nr_roots, this->roots.begin());
    }
    /* Allow assignment-construction from general entries */
    explicit fb_entry_x_roots(fb_entry_general const & e)
        : p(e.p)
        , invq(e.invq)
    {
        ASSERT_ALWAYS(Nr_roots == e.nr_roots);
        ASSERT(e.is_simple());
        for (int i = 0; i < Nr_roots; i++)
            roots[i] = e.roots[i].r;
    }
    fbprime_t get_q() const { return p; }
    fbroot_t get_r(size_t const i) const { return roots[i]; }
    fbroot_t get_proj(size_t const i MAYBE_UNUSED) const { return false; }
    double weight() const { return 1. / p * Nr_roots; }
    /* Allow sorting by p */
    bool operator<(fb_entry_x_roots<Nr_roots> const & other) const
    {
        return this->p < other.p;
    }
    bool operator>(fb_entry_x_roots<Nr_roots> const & other) const
    {
        return this->p > other.p;
    }
    void fprint(FILE *) const;
    void transform_roots(transformed_entry_t &, qlattice_basis const &) const;
};

static_assert(
    std::has_unique_object_representations_v<fb_entry_x_roots<1>>,
    "fb_entry_x_roots<1> must not have padding");
static_assert(
    std::has_unique_object_representations_v<fb_entry_x_roots<2>>,
    "fb_entry_x_roots<2> must not have padding");

/* }}} */

class fb_slice_interface
{
  public:
    fb_slice_interface() = default;
    virtual ~fb_slice_interface() = default;
    virtual size_t size() const = 0;
    virtual unsigned char get_logp() const = 0;
    virtual fbprime_t get_prime(slice_offset_t offset) const = 0;
    virtual unsigned char get_k(slice_offset_t offset) const = 0;
    /* global index across all fb parts */
    virtual slice_index_t get_index() const = 0;
    virtual double get_weight() const = 0;
    virtual bool is_general() const = 0;
    virtual int get_nr_roots() const = 0;
};

/* see fb_slice_weight.hpp */
template <typename FB_ENTRY_TYPE> class fb_slice_weight_estimator;

template <typename FB_ENTRY_TYPE> class fb_slice : public fb_slice_interface
{
    friend class fb_slice_weight_estimator<FB_ENTRY_TYPE>;
    using fb_entry_vector = mmappable_vector<FB_ENTRY_TYPE>;
    typename fb_entry_vector::const_iterator _begin, _end;
    unsigned char logp;
    slice_index_t index; /* global index across all fb parts */
    double weight;
    friend struct helper_functor_subdivide_slices;
    fb_slice(typename fb_entry_vector::const_iterator it, unsigned char logp)
        : _begin(it)
        , _end(it)
        , logp(logp)
        , index(0)
        , weight(0)
    {
    }
    fb_slice(typename fb_entry_vector::const_iterator it,
             typename fb_entry_vector::const_iterator jt, unsigned char logp)
        : _begin(it)
        , _end(jt)
        , logp(logp)
        , index(0)
        , weight(0)
    {
    }

  public:
    using entry_t = FB_ENTRY_TYPE;
    typename fb_entry_vector::const_iterator begin() const
    {
        return _begin;
    }
    typename fb_entry_vector::const_iterator end() const { return _end; }
    size_t size() const override { return _end - _begin; }
    unsigned char get_logp() const override { return logp; }
    fbprime_t get_prime(slice_offset_t offset) const override
    {
        /* While it may make sense to manipulate slices that exceed the
         * expected max size during the work to construct them, it should
         * be pretty clear that we'll never ever try to use get_prime,
         * which is limited in its type, on a slice that has more prime
         * ideals that this function can address */
        ASSERT(size() <= std::numeric_limits<slice_offset_t>::max());
        return _begin[offset].p;
    }
    unsigned char get_k(slice_offset_t offset) const override
    {
        return _begin[offset].k;
        /* power. Well, most often it's a constant ! We need to
         * access it from the virtual base though. This is way we're
         * not folding it to a template access.  */
    }
    /* global index across all fb parts */
    slice_index_t get_index() const override { return index; }
    double get_weight() const override { return weight; }
    bool is_general() const override { return FB_ENTRY_TYPE::is_general_type; }
    /* get_nr_roots() on a fb_slice returns zero for slices of general
     * type ! */
    int get_nr_roots() const override { return FB_ENTRY_TYPE::fixed_nr_roots; }
};

/* entries and general_entries: we declare general entries, as well
 * as entries for all numbers of roots between 1 and MAX_ROOTS
 *
 * (notationally, general_entries corresponds to -1 as a number of
 * roots).
 * */

template <typename T> struct entries_and_cdf {
    using container_type = mmappable_vector<T>;
    using weight_container_type = mmappable_vector<double>;
    struct type : public container_type {
        using container_type = entries_and_cdf<T>::container_type;
        using weight_container_type = entries_and_cdf<T>::weight_container_type;
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
        weight_container_type::const_iterator weight_begin() const
        {
            return weight_cdf.begin();
        }
        weight_container_type::const_iterator weight_end() const
        {
            return weight_cdf.end();
        }
        weight_container_type::size_type weight_size() const
        {
            return weight_cdf.size();
        }
        double weight_delta(size_t a, size_t b) const
        {
            return weight_cdf[b] - weight_cdf[a];
        }
        double weight_delta(typename container_type::const_iterator a,
                            typename container_type::const_iterator b) const
        {
            return weight_cdf[b - container_type::begin()] -
                   weight_cdf[a - container_type::begin()];
        }
        double weight_cdf_at(size_t a) const { return weight_cdf[a]; }
        double weight_cdf_at(typename container_type::const_iterator a) const
        {
            return weight_cdf[a - container_type::begin()];
        }
    };
};
template <typename T>
using entries_and_cdf_t = entries_and_cdf<T>::type;

template <int n>
struct works_with_mmappable_vector<fb_entry_x_roots<n>>
    : public std::true_type {
};
template <>
struct works_with_mmappable_vector<fb_entry_general> : public std::true_type {
};

template <int n> struct fb_entries_factory {
    using type = entries_and_cdf_t<fb_entry_x_roots<n>>;
};
template <> struct fb_entries_factory<-1> {
    using type = entries_and_cdf_t<fb_entry_general>;
};
template <int n> struct fb_slices_factory {
    using type = std::vector<fb_slice<fb_entry_x_roots<n>>>;
};
template <> struct fb_slices_factory<-1> {
    using type = std::vector<fb_slice<fb_entry_general>>;
};
template <typename T, typename U> struct std::common_type<fb_slice<T>, fb_slice<U>> {
    using type = fb_slice_interface;
};

class fb_factorbase
{
  public:
    /* Has to be <= MAX_DEGREE */
    static int const MAX_ROOTS = MAX_DEGREE;

  private:
    cxx_mpz_poly f;
    int side = -1;
    unsigned long lim = 0;
    unsigned long powlim = 0;

  public:
    bool empty() const { return lim == 0; }

  private:
    using entries_t = cado::multityped_array<fb_entries_factory, -1, MAX_ROOTS + 1>;
    entries_t entries;

  public:
    using threshold_pos = std::array<size_t, MAX_ROOTS + 2>;
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
    void append(std::list<fb_entry_general> & pool);
    /* }}} */

    /* {{{ Various ways to count primes, prime ideals, and hit ratio
     * ("weight") of entries in the whole factor base, or in subranges
     */
  private:
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
        std::array<fbprime_t, FB_MAX_PARTS> thresholds;
        fbprime_t td_thresh;
        fbprime_t skipped = 0;
        double scale = 0;
        /* This might seem non obvious, but this parameters controls
         * the size of the slices, because we want to enforce some
         * relatively-even division. It's not entirely clear that we
         * want it here, but we definitely want it somewhere. */
        unsigned int nb_threads = 0;

        bool operator<(key_type const & x) const
        {
            int r;
            for (int i = 0; i < FB_MAX_PARTS; ++i) {
                r = (thresholds[i] > x.thresholds[i]) -
                    (x.thresholds[i] > thresholds[i]);
                if (r)
                    return r < 0;
            }
            r = (td_thresh > x.td_thresh) - (x.td_thresh > td_thresh);
            if (r)
                return r < 0;
            r = (skipped > x.skipped) - (x.skipped > skipped);
            if (r)
                return r < 0;
            return scale < x.scale;
        }
    };

  public:
    class slicing
    {
      public:
        struct stats_type {
            /* explicit-initializing as below forces zero-initialization
             * of members */
            std::array<size_t, FB_MAX_PARTS> primes {{}};
            std::array<size_t, FB_MAX_PARTS> ideals {{}};
            std::array<double, FB_MAX_PARTS> weight {{}};
        };
        stats_type stats;

        class part
        {
            /* slices and general_slices: actually general_slices is
             * slices for number of roots -1, and that is
             * always a vector of length 1. And we have slices (vectors
             * of fb_slice objects) for all numbers of roots between 1 and
             * MAX_ROOTS.
             */
            using slices_t = cado::multityped_array<fb_slices_factory, -1, MAX_ROOTS + 1>;
            friend struct helper_functor_subdivide_slices;

          public:
            slices_t slices;

            slice_index_t first_slice_index = 0;
            template <int n>
            typename fb_slices_factory<n>::type & get_slices_vector_for_nroots()
            {
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
             * slicesG.size() <= i < slices0.size() --> return
             * &(slices0[i-slicesG.size()]) slices0.size() <= i < slices1.size()
             * --> return &(slices1[i-slicesG.size()-slices0.size()]) and so on.
             */
            struct helper_functor_get {
                using type = fb_slice_interface const *;
                using key_type = slice_index_t;
                template <typename T>
                type operator()(T const & x, slice_index_t & k)
                {
                    if ((size_t)k < x.size()) {
                        return &(x[k]);
                    } else {
                        k -= x.size();
                        return nullptr;
                    }
                }
            };

            /* index is the global index across all fb parts */
            fb_slice_interface const * get(slice_index_t index) const
            {
                /* used to be in fb_vector::get_slice
                 *
                 * and in fb_part::get_slice for the lookup of the
                 * vector.
                 *
                 * TODO: dichotomy, perhaps ? Can we do the dichotomy
                 * elegantly ?
                 */


                if (index < first_slice_index)
                    return nullptr;

                /*
                const slice_index_t idx = index - first_slice_index;
                if (idx >= nslices()) {
                    // index = idx - nslices();
                    return nullptr;
                }
                */

                auto nth = [](auto const & B, size_t v) {
                    auto locate = [&](auto const & x) -> fb_slice_interface const * {
                        if (v < x.size())
                            return &(x[v]);
                        v -= x.size();
                        return nullptr;
                    };
                    return B.find(locate);
                };

                auto const * res = nth(slices, index - first_slice_index);
                ASSERT_ALWAYS(res);
                ASSERT_ALWAYS(res->get_index() == index);
                return res;
            }

            /* {{{ use caching for the number of slices, as it's a handy
             * thing to query */
            mutable slice_index_t _nslices =
                std::numeric_limits<slice_index_t>::max();

          public:
            slice_index_t nslices() const
            {
                if (_nslices != std::numeric_limits<slice_index_t>::max())
                    return _nslices;
                auto totalsize = [](auto const & B) {
                    size_t s = 0;
                    B.foreach([&](auto const & x) { s += x.size(); });
                    return s;
                };
                return _nslices = totalsize(slices);
            }
            /* }}} */

          public:
            template <typename F> void foreach_slice(F && f)
            {
                slices.foreach([&](auto & x) {
                    for (auto & a: x)
                        std::forward<F>(f)(a);
                });
            }
            template <typename F> void foreach_slice(F && f) const
            {
                slices.foreach([&](auto const & x) {
                    for (auto const & a: x)
                        std::forward<F>(f)(a);
                });
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

            /* index: global index across all fb parts */
            fb_slice_interface const & operator[](slice_index_t index) const
            {
                /* This bombs out at runtime if get returns nullptr, but
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

      public:
        size_t nparts() const { return parts.size(); }

        /* index: global index across all fb parts */
        fb_slice_interface const * get(slice_index_t index, size_t i0, size_t i1) const
        {
            for (size_t i = i0 ; i < i1 ; i++) {
                auto const & p = parts[i];
                if (index < p.first_slice_index + p.nslices())
                    return p.get(index);
            }
            return nullptr;
        }
        fb_slice_interface const * get(slice_index_t index, size_t i0) const
        {
            return get(index, i0, parts.size());
        }
        fb_slice_interface const * get(slice_index_t index) const
        {
            return get(index, 0, parts.size());
        }

        part const & get_part(int i) const { return parts[i]; }

        int get_toplevel() const { return toplevel; }

        /* This accesses the *fb_slice* with this index. Not the vector of
         * slices !
         * index: global index across all fb parts */
        fb_slice_interface const & operator[](slice_index_t index) const
        {
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

        template <typename T> void foreach_slice(T & f)
        {
            for (auto & p: parts) {
                p.foreach_slice(f);
            }
        }

        slicing() = default;
        slicing(fb_factorbase const & fb, key_type const & K);
    };

  private:
    lock_guarded_container<std::map<key_type, slicing>> cache;
    int read(char const * filename);

  public:
    /* accessors.
     * As in the std::map case, the non-const accessor is allowed to
     * create stuff. */

    slicing & operator[](key_type const & K)
    {
        const std::lock_guard<std::mutex> foo(cache.mutex());
        auto it = cache.find(K);
        if (it != cache.end())
            return it->second;

        /* create a new slot */
        slicing & res(cache[K] = slicing(*this, K));
        return res;
    }

    slicing const & operator[](key_type const & K) const
    {
        return cache.at(K);
    }

  private:
    void make_linear();
    void make_linear_threadpool(unsigned int nb_threads);

  public:
    /* Note that the intent here is to read the factor base once and
     * for all. In the descent context, where we have several
     * different fb lim parameters, we should get away with different
     * slicings, instead of trying to redo the factor base
     * initialization.
     */
    fb_factorbase(cxx_cado_poly const & cpoly, int side, cxx_param_list & pl,
                  char const * fbc_filename, int nthreads = 1);
    fb_factorbase() = default;
    /*
    fb_factorbase(fb_factorbase &&) = default;
    fb_factorbase & operator=(fb_factorbase &&) = default;
    */

  public:
    void finalize()
    {
        entries.foreach([](auto & x) {
            /* not entirely clear to me. Shall we sort by q or by p ?
             */
            using X = std::remove_reference_t<decltype(x)>::value_type;
            auto by_q = [](X const & a, X const & b) {
                return a.get_q() < b.get_q();
            };
            /*
            auto by_p_then_q = [](X const & a, X const & b) {
                return
                    a.get_p() < b.get_p() ||
                    a.get_p() == b.get_p() &&
                    a.get_q() < b.get_q();
            };
            */
            std::ranges::sort(x, by_q);
            });
    }
};

std::ostream & operator<<(std::ostream & o, fb_factorbase::key_type const &);

namespace fmt {
    template<> struct formatter<fb_factorbase::key_type> : ostream_formatter {};
} /* namespace fmt */

unsigned char fb_log(double x, double y, double z);
unsigned char fb_log_delta(fbprime_t, unsigned long, unsigned long, double);
fbprime_t fb_is_power(fbprime_t, unsigned long *);

#endif /* CADO_FB_HPP */
