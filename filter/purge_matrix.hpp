#ifndef CADO_PURGE_MATRIX_HPP
#define CADO_PURGE_MATRIX_HPP

#include <cstddef>
#include <cstdint>
#include <cstdio>

#include <compare>
#include <limits>
#include <mutex>
#include <type_traits>
#include <utility>
#include <vector>

#include "memalloc.hpp"
#include "typedefs.h"

/* XXX implementation for this structure is split in purge_matrix.cpp
 * (most stuff, including singleton removal) and clique_removal.cpp
 */
template <typename T>
concept container_with_h = requires {
    T().begin();
    T().end();
    T()[0].h;
};

template <typename F>
concept weight_function = requires { float(F()(3)); };

struct purge_matrix : simple_minded_chunk_allocator<index_t> {
    static constexpr index_t END_OF_ROW = std::numeric_limits<index_t>::max();
    static constexpr weight_t OVERWEIGHT = std::numeric_limits<weight_t>::max();
    std::vector<index_t *> rows;
    std::vector<weight_t> column_weights;
    size_t remaining_rows = 0;
    size_t remaining_columns = 0;
    private:
    mutable std::mutex m;
    public:

    std::make_signed_t<size_t> excess() const
    {
        return remaining_rows - remaining_columns;
    }

    bool is_active(size_t i) const { return rows[i] != nullptr; }

    /* Insert row i into the matrix, containing the primes in the given
     * list that are within the column index range [col0, col1).
     *
     * remaining_rows increases by one
     * remaining_columns increases if needed.
     *
     * the matrix is reallocated if needed, with lock protection.
     */
    template <container_with_h C>
    void new_row(size_t i, size_t col0, size_t col1, C const & primes)
    {
        size_t active_weight = 0;
        for (auto & c: primes)
            active_weight += (c.h >= col0 && c.h < col1);

        auto * p = alloc(active_weight + 1);
        p[active_weight] = END_OF_ROW;
        auto * q = p;
        size_t newcols = 0;
        for (auto & c: primes) {
            auto h = c.h;
            if (h >= column_weights.size()) {
                std::scoped_lock const dummy(m);
                const size_t z = column_weights.size();
                if (h >= z)
                    column_weights.insert(column_weights.end(), h + 1 - z, 0);
            }
            auto & w = column_weights[h];
            if (w == 0)
                newcols++;
            if (w < OVERWEIGHT)
                w++;
            if (h >= col0 && h < col1)
                *q++ = h;
        }
        std::scoped_lock const dummy(m);
        remaining_rows++;
        remaining_columns += newcols;
        if (i >= rows.size())
            rows.insert(rows.end(), i + 1 - rows.size(), nullptr);
        rows[i] = p;
    }

    /* for the purposes of clique removal, this computes the addition
     * i1+i2 of the row indices of the non-zero coefficients,
     * specifically in the columns of weight 2.
     */
    using sum2_type = std::vector<size_t>;

    void print_stats_column_weights(FILE *, int verbose) const;
    void print_stats_row_weights(FILE *, int verbose) const;
    void print_stats_on_cliques (FILE *, sum2_type const &, int verbose) const;

    int64_t singleton_removal(size_t nthreads, int verbose);

    /* everything about clique removal is in clique_removal.cpp */
    void clique_removal(int64_t target_excess, size_t nthreads, int verbose);

    struct connected_component {
        size_t i = 0;     /* smallest row */
        float w = 0;      /* weight */

        /* compare weights, first. */
        std::partial_ordering operator<=>(connected_component const & o) const {
            if (auto r = w <=> o.w; r != 0) return r;
            /* same weight, use first row to be deterministic */
            return i <=> o.i;
        }
        auto operator==(connected_component const & o) const {
            return i == o.i;
        }
    };

    /* in clique_removal.cpp */
    static void print_clique_removal_weight_function();

    private:

    /* compute_connected_component() has a hidden quadratic behavior.
     * It's fine if all connected components remain tiny, but it can
     * become a pain quite soon if they don't. The set version has better
     * complexity. As a stopgap measure, we keep both for the moment, and
     * we'll see if there's a point in having only the set version in
     * longer term
     */
    static constexpr size_t ccc_quadratic_abort = 10;

    template <weight_function F>
    std::pair<size_t, connected_component>
    compute_connected_component(size_t i0, std::vector<size_t> const & sum2,
                                    std::vector<size_t> & buf,
                                    F const & f) const;

    template <weight_function F>
    std::pair<size_t, connected_component>
    compute_connected_component_with_set(size_t i0,
            std::vector<size_t> const & sum2,
            F const & f) const;

    void delete_connected_component(size_t i0,
            std::vector<size_t> const & sum2,
            std::vector<size_t> & buf);

    size_t clique_removal_core (sum2_type const &,
            int64_t target_excess, size_t max_nb_comp,
            size_t nthreads, int verbose);

    size_t one_singleton_removal(size_t nthreads);

    /* this function needs to touch the matrix data atomically. For this
     * reason, we prefer to limit its exposition. It is used concurrently
     * within singleton_removal and clique_removal
     */
    void delete_row(size_t i);

    sum2_type compute_sum2(size_t nthreads = 1) const;

    static void print_stats_generic(FILE * out, std::vector<size_t> const & w,
                        char const * name, char const * unit, int verbose);
};

#endif /* CADO_PURGE_MATRIX_HPP */
