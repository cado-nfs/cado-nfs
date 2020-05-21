#ifndef RENUMBER_HPP_
#define RENUMBER_HPP_

#include <cstdint>     /* AIX wants it first (it's a bug) */
#include <algorithm>
#include <cstdio>
#include <cstring>
#include <string>
#include <cerrno>
#include <array>
#include <vector>
#include <stdexcept>
#include <type_traits>
#include "typedefs.h"
#include "cado_poly.h"
#include "badideals.hpp"

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
    struct corrupted_table : public std::runtime_error {
        corrupted_table(std::string const &);
    };
    struct p_r_side {
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
    };
private:
    /* all the (p,r,side) description of the bad ideals */
    std::vector<std::pair<p_r_side, badideal> > bad_ideals;

    p_r_values_t bad_ideals_max_p = 0;

    cxx_cado_poly cpoly;

    /* What actually goes in the data[] table is implementation
     * dependent. Currently we have several (3) choices. See in
     * renumber.cpp
     */
    std::vector<p_r_values_t> traditional_data;
    std::vector<std::array<p_r_values_t, 2>> flat_data;

    std::vector<unsigned int> lpb;

    std::vector<index_t> index_from_p_cache;

    /*
     * [0..above_add): additional columns
     * [above_add..above_bad): bad ideals
     * [above_bad..above_cache): ideals in the big data table, p cached 
     * [above_cache..above_all): rest
     *
     * These are the outer indices only. The internal table size depends
     * on the renumber format that is is use (see renumber.cpp).
     *
     * traditional: traditional_data.size() == above_all - above_bad
     * variant: traditional_data.size() == above_all - above_bad + nprimes
     * flat: flat_data.size() == above_all - above_bad
     */
    index_t above_add = 0;
    index_t above_bad = 0;
    index_t above_cache = 0;
    index_t above_all = 0;

public:
    inline unsigned int get_lpb(unsigned int i) const { return lpb[i]; }
    inline unsigned int get_max_lpb() const { return *std::max_element(lpb.begin(), lpb.end()); }
    inline unsigned int get_min_lpb() const { return *std::min_element(lpb.begin(), lpb.end()); }
    inline uint64_t get_size() const { return above_all; }
    inline unsigned int get_nb_polys() const { return cpoly->nb_polys; }
    inline p_r_values_t get_max_bad_p() const { return bad_ideals_max_p; }
    inline mpz_poly_srcptr get_poly(int side) const { return cpoly->pols[side]; }
    inline int get_poly_deg(int side) const { return get_poly(side)->deg; }
    inline int get_rational_side() const {
        for(unsigned int side = 0 ; side < get_nb_polys() ; side++) {
            if (get_poly_deg(side) == 1) return side;
        }
        return -1;
    }
    inline index_t get_max_index() const { return above_all; }
    inline index_t get_max_cached_index() const { return above_cache; }
    inline p_r_values_t get_cached_primes_bound() const { return index_from_p_cache.size(); }
    inline index_t number_of_additional_columns() const { return above_add; }
    std::vector<int> get_sides_of_additional_columns() const;
    inline index_t number_of_bad_ideals() const { return above_bad - above_add; }
    inline size_t get_memory_size() const {
        return traditional_data.size() * sizeof(decltype(traditional_data)::value_type)
            + flat_data.size() * sizeof(decltype(flat_data)::value_type);
    }
    renumber_t() = default;
    // ~renumber_t() = default;
    renumber_t(renumber_t const &) = delete;
    renumber_t& operator=(renumber_t const &) = delete;
    renumber_t(renumber_t &&) = default;
    renumber_t& operator=(renumber_t &&) = default;

    renumber_t(cxx_cado_poly const & cpoly) : cpoly(cpoly), lpb(cpoly->nb_polys, 0) {}

    void set_lpb(std::vector<unsigned int> const & x) {
        ASSERT_ALWAYS(x.size() == lpb.size());
        lpb = x;
    }
    void use_additional_columns_for_dl();
    void compute_bad_ideals();
    void compute_bad_ideals_from_dot_badideals_hint(std::istream&, unsigned int = UINT_MAX);

    void read_header(std::istream& os);
    void write_header(std::ostream& os) const;
    void read_bad_ideals(std::istream& is);
    void write_bad_ideals(std::ostream& os) const;
    /* transitory, for interface that produces renumber file in legacy
     * format.
     */
    void read_bad_ideals_info(std::istream & is);
    /* there's no write_table, because writing the table is done by
     * the build() function (called from freerel)
     */
    void read_table(std::istream& is);

    /* first version uses only the renumber table. Second version is for
     * the old format, and we need the badidealinfo file as well (not the
     * badideals file, which does not contain enough info.
     */
    void read_from_file(const char * filename);
    void read_from_file(const char *, const char *);
    // renumber_t(const char *);
    // renumber_t(const char *, const char *);

    /* return the number of bad ideals above x (and therefore zero if
     * x is not bad). If the ideal is bad, put in the reference [first] the
     * first index that corresponds to the bad ideals.
     */
    int is_bad(p_r_side) const;
    int is_bad(index_t &, p_r_side) const;
    /* two convenience shortcuts, to avoid curlies */
    inline int is_bad(p_r_values_t p, p_r_values_t r, int side) const {
        return is_bad({p, r, side});
    }
    inline int is_bad(index_t & index, p_r_values_t p, p_r_values_t r, int side) const {
        return is_bad(index, {p, r, side});
    }

    bool is_bad (index_t h) const {
        return h >= above_add && h < above_bad;
    }
    bool is_additional_column (index_t h) const {
        return h < above_add;
    }

    /* None of these access interfaces work for additional columns. */
    index_t index_from_p_r (p_r_side) const;
    inline index_t index_from_p_r (p_r_values_t p, p_r_values_t r, int side) const {
        return index_from_p_r({p, r, side});
    }
    p_r_side p_r_from_index (index_t) const;

    /* This second interface works for bad ideals as well. */
    std::pair<index_t, std::vector<int>> indices_from_p_a_b(p_r_side x, int e, int64_t a, uint64_t b) const;

    struct cooked {
        std::vector<unsigned int> nroots;
        std::vector<p_r_values_t> traditional;
        std::vector<std::array<p_r_values_t, 2>> flat;
        std::string text;
    };
    /* This is purely descriptive, used for debugging */
    std::string debug_data(index_t i) const;

    struct hook {
        virtual void operator()(renumber_t & R, p_r_values_t p, index_t idx, renumber_t::cooked const & C) = 0;
        virtual ~hook() = default;
    };

    /* To build a renumber table in memory in the simplest way, the
     * process goes as follows
       renumber_t renumber_table(cpoly);
       renumber_table.set_lpb(lpb);
       renumber_table.build();
     */

    static void builder_configure_switches(cxx_param_list &);
    static void builder_declare_usage(cxx_param_list &);
    static void builder_lookup_parameters(cxx_param_list &);
    index_t build(cxx_param_list &, hook * = nullptr);
    index_t build(hook * = nullptr);


private:
    unsigned int needed_bits() const;
    /* this returns an index i such that data[i - above_bad] points to
     * the beginning of data for p. Note that in the
     * renumber_format_variant case, this does not mean that i is the
     * index of the first prime above p !
     */
    index_t get_first_index_from_p(p_r_side x) const;

    p_r_values_t compute_vr_from_p_r_side (p_r_side x) const;
    p_r_side compute_p_r_side_from_p_vr (p_r_values_t p, p_r_values_t vr) const;
    p_r_values_t compute_vp_from_p (p_r_values_t p) const;
    p_r_values_t compute_p_from_vp (p_r_values_t vp) const;

    bool traditional_get_largest_nonbad_root_mod_p (p_r_side & x) const;
    index_t traditional_backtrack_until_vp(index_t i, index_t min) const;
    bool traditional_is_vp_marker(index_t i) const;


    /* The "cook" function can be used asynchronously to prepare the
     * fragments of the renumber table in parallel. use_cooked must use
     * the same data, but synchronously -- and stores it to the table, of
     * course. The use_cooked_nostore does the same, except that it is
     * made for the situation where we have no interest in keeping track
     * of the renumber table itself. The only thing that matters is
     * keeping track of the above_all index, which is done by the input
     * and output index_t values. */
    cooked cook(unsigned long p, std::vector<std::vector<unsigned long>> &) const;
    void use_cooked(p_r_values_t p, cooked const & C);
    index_t use_cooked_nostore(index_t n0, p_r_values_t p, cooked const & C);

    struct builder;
    friend struct builder;
};
#endif /* RENUMBER_HPP_ */
