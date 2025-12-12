#ifndef CADO_RENUMBER_HPP
#define CADO_RENUMBER_HPP

#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdio>

#include <algorithm>
#include <array>
#include <istream>
#include <iterator>
#include <map>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <type_traits>

#include <sys/types.h>

#include "badideals.hpp"
#include "cado_poly.h"
#include "cxx_mpz.hpp"
#include "macros.h"
#include "mpz_poly.h"
#include "typedefs.h"

struct cxx_param_list; // IWYU pragma: keep

/* To build a renumber table in memory in the simplest way, the
 * process goes as follows

 renumber_t renumber_table(cpoly);
 renumber_table.set_lpb(lpb);
 renumber_table.build();

 (there are various ways to control the build process, including a way to
 stream the computed table to a file, and compute free relations as well.
 This is done in freerel.cpp ; but the heavylifting really happens in
 renumber proper.)

 Note that by default, the renumber tables include the bad ideals
 information. There is currently no way to turn it off. It should
 probably be done.

 * To read a renumber table from a file, this goes as:

 renumber_t renumber_table(cpoly);
 renumber_table.read_from_file(renumberfilename);

*/

#define RENUMBER_MAX_LOG_CACHED 20

struct renumber_t {
    struct corrupted_table : public std::runtime_error {/*{{{*/
        explicit corrupted_table(std::string const &);
    };/*}}}*/
    struct p_r_side {/*{{{*/
        p_r_values_t p;
        p_r_values_t r;
        int side;
        bool operator==(p_r_side const & x) const { return p == x.p && r == x.r  && side == x.side; }
        bool same_p(p_r_side const & x) const { return p == x.p && side == x.side; }
        bool operator<(p_r_side const & x) const { 
            int s;
            if ((s = (side > x.side) - (x.side > side)) != 0) return s < 0;
            if ((s = (p > x.p) - (x.p > p)) != 0) return s < 0;
            if ((s = (r > x.r) - (x.r > r)) != 0) return s < 0;
            return false;
        }
    };/*}}}*/

    /* As of 20220311, I'm killing the old and transitional formats.
     *
     * 20220411 is the same as 20200515, but with some hex/bin mess
     * cleaned up
     *
     */
    static constexpr const int format_flat = 20220411;

private: /*{{{ internal data fields*/

    int format = format_flat;

    /* all the (p,r,side) description of the bad ideals */
    std::vector<std::pair<p_r_side, badideal>> bad_ideals;

    p_r_values_t bad_ideals_max_p = 0;

    cxx_cado_poly cpoly;

    /* Only for format_flat. */
    std::vector<std::array<p_r_values_t, 2>> flat_data;

    std::vector<unsigned int> lpb;
    std::vector<index_t> index_from_p_cache;

    /* This is computed on the fly when the table is built. It's cheap
     * enough anyway. We rarely use it, and we're a bit too lazy to
     * change the file format to include these.
     *
     * note that the ramified primes are a subset of these.
     */
    std::vector<std::vector<std::pair<cxx_mpz,int>>> small_primes;

    mutable std::map<p_r_side, int> exceptional_inertia;

    /*
     * [0..above_add): additional columns
     * [above_add..above_bad): bad ideals
     * [above_bad..above_cache): ideals in the big data table, p cached 
     * [above_cache..above_all): rest
     *
     * These are the outer indices only. The internal table size depends
     * on the renumber format that is is use (see renumber.cpp).
     *
     * format_flat: flat_data.size() == above_all - above_bad
     */
    index_t above_add = 0;
    index_t above_bad = 0;
    index_t above_cache = 0;
    index_t above_all = 0;
/*}}}*/

public:
    /* various accessors {{{*/
    int get_format() const { return format; }
    unsigned int get_lpb(int i) const { return lpb[i]; }
    unsigned int get_max_lpb() const { return *std::max_element(lpb.begin(), lpb.end()); }
    unsigned int get_min_lpb() const { return *std::min_element(lpb.begin(), lpb.end()); }
    size_t size() const { return above_all; }
    int get_nb_polys() const { return cpoly->nb_polys; }
    mpz_poly_srcptr get_poly(int side) const { return cpoly->pols[side]; }
    int get_poly_deg(int side) const { return get_poly(side)->deg; }
    int get_rational_side() const {
        for(int side = 0 ; side < get_nb_polys() ; side++) {
            if (get_poly_deg(side) == 1) return side;
        }
        return -1;
    }
    index_t get_max_index() const { return above_all; }
    index_t get_max_cached_index() const { return above_cache; }
    index_t number_of_additional_columns() const { return above_add; }
    /* sides where the J ideal is non trivial (even in the degree 1
     * case). So it's really the list of sides where f is not monic
     */
    std::vector<int> get_sides_of_additional_columns() const;
    index_t number_of_bad_ideals() const { return above_bad - above_add; }
    size_t get_memory_size() const {
        return flat_data.size() * sizeof(decltype(flat_data)::value_type);
    }
/*}}}*/

    /*{{{ default ctors */
    renumber_t() = default;
    ~renumber_t() = default;
    renumber_t(renumber_t const &) = delete;
    renumber_t& operator=(renumber_t const &) = delete;
    renumber_t(renumber_t &&) noexcept = default;
    renumber_t& operator=(renumber_t &&) = default;
    /*}}}*/

    explicit renumber_t(cxx_cado_poly const & cpoly) : cpoly(cpoly), lpb(cpoly->nb_polys, 0) {}
    renumber_t(cxx_cado_poly const & cpoly,
            std::string const & filename,
            bool for_dl)
        : renumber_t(cpoly)
    {
        read_from_file(filename.c_str(), for_dl);
    }

    /*{{{ configuration when creating the table */
    void set_lpb(std::vector<unsigned int> const & x) {
        ASSERT_ALWAYS(x.size() == lpb.size());
        lpb = x;
    }
    /*}}}*/

    /*{{{ reading the table */
    void read_from_file(const char * filename, bool for_dl);
    void recompute_debug_number_theoretic_stuff();
    /*}}}*/

    /*{{{ most important outer-visible routines: lookups */

    /* special lookups:
     *  a lookup of an additional column on side s returns {0, 0, s}
     *  a lookup of a bad ideal (given by (p,r)) returns the index of the
     *  _first_ ideal above this (p,r)
     * except for the non-injectivity of the mapping index->(p,r) for bad
     * ideals, the map is "almost" a bijection.
     */

    /* return the number of bad ideals above x (and therefore zero if
     * x is not bad) ; likewise for index h.
     * If the ideal is bad, put in the reference [first] the
     * first index that corresponds to the bad ideals.
     */
    int is_bad(p_r_side) const;
    int is_bad(index_t &, p_r_side) const;
    int is_bad(index_t & first_index, index_t h) const;

    /* two convenience shortcuts, to avoid curlies */
    int is_bad(p_r_values_t p, p_r_values_t r, int side) const {
        return is_bad({p, r, side});
    }
    int is_bad(index_t & index, p_r_values_t p, p_r_values_t r, int side) const {
        return is_bad(index, {p, r, side});
    }

    bool is_bad (index_t h) const {
        return h >= above_add && h < above_bad;
    }
    bool is_additional_column (index_t h) const {
        return h < above_add;
    }

    bool has_merged_additional_column() const {
        return get_nb_polys() == 2 && get_sides_of_additional_columns().size() == 2;
    }

    index_t index_from_p_r (p_r_side) const;
    index_t index_from_p_r (p_r_values_t p, p_r_values_t r, int side) const {
        return index_from_p_r({p, r, side});
    }
    p_r_side p_r_from_index (index_t) const;

    int inertia_from_p_r(p_r_side) const;
    int inertia_from_p_r(p_r_values_t p, p_r_values_t r, int side) const {
        return inertia_from_p_r({p, r, side});
    }


    class const_iterator;

    index_t index_from_p(p_r_values_t p0) const;
    const_iterator iterator_from_p(p_r_values_t p0) const;
    index_t index_from_p(p_r_values_t p0, int side) const;
    const_iterator iterator_from_p(p_r_values_t p0, int side) const;

    /* This second interface works for bad ideals as well. */
    std::pair<index_t, std::vector<int>> indices_from_p_a_b(p_r_side x, int e, int64_t a, uint64_t b) const;
    /*}}}*/

    /* {{{ build() functionality */
    /* To build a renumber table in memory in the simplest way, the
     * process goes as follows
       renumber_t renumber_table(cpoly);
       renumber_table.set_lpb(lpb);
       renumber_table.build();
     */

    struct cooked {
        std::vector<int> nroots;
        std::vector<std::array<p_r_values_t, 2>> flat;
        std::string text;
        bool empty() const { return flat.empty(); }
    };

    struct hook {
        virtual void operator()(renumber_t & R, p_r_values_t p, index_t idx, renumber_t::cooked const & C) = 0;
        virtual ~hook() = default;
    };

    static void builder_declare_usage(cxx_param_list &);
    static void builder_lookup_parameters(cxx_param_list &);
    index_t build(cxx_param_list &, bool for_dl, hook * = nullptr);
    index_t build(bool for_dl, hook * = nullptr);
    /* }}} */

    /*{{{ debugging aids*/
    std::string debug_data(index_t i) const;
    std::string debug_data_sagemath(index_t i) const;
    std::string debug_data_machine_description(index_t i) const;
    void info(std::ostream & os) const;
    void more_info(std::ostream & os) const;
    /*}}}*/

private:/*{{{ more implementation-level stuff. */
    void read_header(std::istream& os);
    void read_bad_ideals(std::istream& is);
    /* there's no write_table, because writing the table is done by
     * the build() function (called from freerel) */
    void read_table(std::istream& is);
    void compute_bad_ideals();
    void compute_bad_ideals_from_dot_badideals_hint(std::istream&, unsigned int = UINT_MAX);
    void compute_ramified_primes();
    void write_header(std::ostream& os) const;
    void write_bad_ideals(std::ostream& os) const;
    /* these two could be made public, I believe. The public way to do
     * the same is to use the param_list argument to build()
     */
    void use_additional_columns_for_dl();
    void set_format(int);

    unsigned int needed_bits() const;
    /* this returns an index i such that data[i - above_bad] points to
     * the beginning of data for p.
     */
    const_iterator get_first_iterator_from_p(p_r_values_t p) const;

    p_r_values_t compute_vr_from_p_r_side (p_r_side x) const;
    p_r_side compute_p_r_side_from_p_vr (p_r_values_t p, p_r_values_t vr) const;
    p_r_values_t compute_vp_from_p (p_r_values_t p) const;
    p_r_values_t compute_p_from_vp (p_r_values_t vp) const;

    /* The "cook" function can be used asynchronously to prepare the
     * fragments of the renumber table in parallel. use_cooked must use
     * the same data, but synchronously -- and stores it to the table, of
     * course. use_cooked_nostore does the same, except that it is made
     * for the situation where we have no interest in keeping track of
     * the renumber table itself. The only thing that matters is keeping
     * track of the above_all index, which is done by the input and
     * output index_t values.
     */
    cooked cook(unsigned long p, std::vector<std::vector<unsigned long>> &) const;
    void use_cooked(p_r_values_t p, cooked const & C);
    index_t use_cooked_nostore(index_t n0, p_r_values_t p, cooked const & C);

    struct builder; // IWYU pragma: keep
    friend struct builder;
/*}}}*/

public:

    friend class const_iterator;

    class const_iterator
    {
        friend struct renumber_t;
        private:
            renumber_t const * table;
            /* these are outer indices when below above_bad, and then we have
             * above_bad + the inner index.  Subtract table.above_bad to get
             * inner table indices.
             */
            index_t i;
        public:
            typedef p_r_side                value_type;
            typedef ptrdiff_t               difference_type;
            typedef p_r_side const *        const_pointer;
            typedef p_r_side const &        const_reference;
            typedef std::input_iterator_tag iterator_category;

            const_iterator(renumber_t const & table,
                    index_t i)
                : table(&table)
                , i(i)
            {}

            p_r_side operator*() const;
            std::array<p_r_values_t, 2> raw() const;
            bool operator==(const const_iterator& other) const { return i == other.i; }
            bool operator!=(const const_iterator& other) const { return !(*this == other); }
            const_iterator operator++(int);
            const_iterator& operator++();
    };

    const_iterator begin() const;
    const_iterator end() const;
};
#endif /* CADO_RENUMBER_HPP */

