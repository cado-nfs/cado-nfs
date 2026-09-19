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
#include "cado_poly.hpp"
#include "cxx_mpz.hpp"
#include "macros.h"
#include "mmappable_vector.hpp"
#include "mpz_poly.h"
#include "typedefs.h"

namespace cado::params {
struct cxx_param_list;
}
using cxx_param_list = cado::params::cxx_param_list;

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

    /* format_binary is the same table, but the (p, vr) pairs are stored
     * as a binary blob at a known offset in the file, so that reading
     * the table is an mmap() and not a parse. The header stays text,
     * and carries the offset of the blob, the size in bytes of one
     * p_r_values_t on the machine that wrote it, and the byte order.
     * See renumber.cpp for the layout.
     *
     * The header is the one of format_flat, plus a last line giving
     * that offset, the entry width, and the byte order. That line is
     * padded with spaces so that it ends where the data begins, which
     * is how a reader that cannot seek gets there.
     *
     * The offset is a fixed multiple of 4096. That is a property of the
     * format, not of the page size of either machine: mmap() wants a
     * page-aligned offset, but mmapped_file::mapping obtains one by
     * rounding down at run time, so a table written where pages are 4kB
     * reads fine where they are 16kB, and conversely. What the offset
     * does have to respect is the alignment of the entries.
     */
    static constexpr const int format_binary = 20250919;

private: /*{{{ internal data fields*/

    int format = format_flat;

    /* all the (p,r,side) description of the bad ideals */
    std::vector<std::pair<p_r_side, badideal>> bad_ideals;

    p_r_values_t bad_ideals_max_p = 0;

    cxx_cado_poly cpoly;

    /* This is an mmappable_vector because we want to be able to map it
     * straight from the file when the format allows it.
     */
    mmappable_vector<std::array<p_r_values_t, 2>> flat_data;

    /* Only meaningful for format_binary, and only once the header has
     * been read: where the blob starts in the file, and how wide its
     * entries are.
     */
    size_t binary_data_offset = 0;
    unsigned int binary_element_size = 0;

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
    int get_nb_polys() const { return cpoly.nsides(); }
    mpz_poly_srcptr get_poly(int side) const { return cpoly[side]; }
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

    explicit renumber_t(cxx_cado_poly const & cpoly) : cpoly(cpoly), lpb(cpoly.nsides(), 0) {}
    renumber_t(cxx_cado_poly const & cpoly,
            std::string const & filename,
            bool for_dl)
        : renumber_t(cpoly)
    {
        read_from_file(filename, for_dl);
    }

    /*{{{ configuration when creating the table */
    void set_lpb(std::vector<unsigned int> const & x) {
        ASSERT_ALWAYS(x.size() == lpb.size());
        lpb = x;
    }
    /*}}}*/

    /*{{{ reading the table */
    void read_from_file(std::string const & filename, bool for_dl);
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
    std::pair<index_t, std::vector<int>> indices_from_p_a_b(
            p_r_side x,
            int e,
            mpz_srcptr a,
            mpz_srcptr b) const;
    /*}}}*/

    /* {{{ build() functionality */
    /* To build a renumber table in memory in the simplest way, the
     * process goes as follows
       renumber_t renumber_table(cpoly);
       renumber_table.set_lpb(lpb);
       renumber_table.build();
     */

    /* What one thread computes for one interval of primes: the entries
     * of the table for those primes, and, for each prime that has at
     * least one ideal above it and in increasing order, the number of
     * roots on each side. The prime itself needs not be stored, since
     * it is the first coordinate of the entries.
     *
     * The index of the first entry of the fragment in the whole table
     * is only known once the previous fragments are complete, which is
     * why anything that needs it (the hook, below) runs in a second
     * pass.
     */
    struct fragment {
        std::vector<std::array<p_r_values_t, 2>> flat;
        std::vector<uint8_t> nroots;    /* nsides per prime */
        std::string text;               /* only for format_flat */
        std::string hook_text;
        index_t base = 0;
        uint64_t nprimes_seen = 0;
        bool empty() const { return flat.empty(); }
        void clear() {
            flat.clear();
            nroots.clear();
            text.clear();
            hook_text.clear();
            base = 0;
            nprimes_seen = 0;
        }
    };

    struct hook {
        /* Called from several threads at once, on distinct fragments,
         * and hence forbidden to touch any shared state: the output
         * goes to the per-fragment string.
         */
        virtual void operator()(renumber_t const & R, p_r_values_t p,
                index_t idx, uint8_t const * nroots, std::string & out) = 0;
        /* Called single-threaded, with the fragments in increasing
         * order of the primes they cover.
         */
        virtual void flush(std::string const & out) = 0;
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
    void read_table_binary(std::istream& is, std::string const & filename,
            bool may_mmap);
    /* fills index_from_p_cache and above_cache, once flat_data is
     * there. Cheap: it only looks at the primes below 2^20.
     */
    void compute_index_from_p_cache();
    /* header (+ bad ideals) as a string, padded so that the binary blob
     * that follows starts on a page boundary
     */
    std::string header_string_with_padding() const;
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

    /* Append to the fragment what the table holds for the prime p,
     * given its roots on each side. Called from several threads at
     * once, on distinct fragments.
     */
    void cook_into(unsigned long p,
            std::vector<std::vector<unsigned long>> & roots,
            fragment & F) const;

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

