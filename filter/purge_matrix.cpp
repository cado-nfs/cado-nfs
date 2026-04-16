#include "cado.h" // IWYU pragma: keep

#include <cmath>
#include <cstdint>
#include <cstdio>

#include <algorithm>
#include <atomic>
#include <limits>
#include <string>
#include <thread>
#include <vector>

#include "fmt/base.h"
#include "fmt/format.h"

#include "macros.h"
#include "purge_matrix.hpp"
#include "subdivision.hpp"
#include "timing.h"
#include "typedefs.h"
#include "utils_cxx.hpp"

void purge_matrix::print_stats_generic(FILE * out, std::vector<size_t> const & w,
                        char const * name, char const * unit, int verbose)
{
    size_t av = 0, std = 0, nzero = 0;
    size_t min = std::numeric_limits<size_t>::max();
    size_t max = 0;
    for (auto c: w) {
        if (c) {
            nzero++;
            min = std::min(min, c);
            max = std::max(max, c);
            av += c;
            std += c * c;
        }
    }
    double const av_f = double_ratio(av, nzero + !nzero);
    double const std_f = sqrt(double_ratio(std, nzero + !nzero) - av_f * av_f);

    fmt::print(out, "# STATS on {}: #{} = {}\n", name, name, nzero);
    fmt::print(out, "# STATS on {}: min {} = {}\n", name, unit, min);
    fmt::print(out, "# STATS on {}: max {} = {}\n", name, unit, max);
    fmt::print(out, "# STATS on {}: av {} = {:.2f}\n", name, unit, av_f);
    fmt::print(out, "# STATS on {}: std {} = {:.2f}\n", name, unit, std_f);

    if (verbose > 1) {
        std::vector<size_t> dist(max - min + 1, 0);
        for (auto c: w)
            if (c)
                dist[c - min]++;

        for (size_t i = min; i <= max; i++)
            if (dist[i - min] > 0)
                fmt::print(out, "# STATS on {}: #{} of {} {} : {}\n", name,
                           name, unit, i, dist[i - min]);
    }
    fflush(out);
}

void purge_matrix::print_stats_column_weights(FILE * out, int verbose) const
{
    std::vector<size_t> w(column_weights.size(), 0);

    for (auto const * p: rows) {
        if (p)
            for (auto const * q = p; *q != END_OF_ROW; q++)
                w[*q]++;
    }

    print_stats_generic(out, w, "cols", "weight", verbose);
}

void purge_matrix::print_stats_row_weights(FILE * out, int verbose) const
{
    std::vector<size_t> w(rows.size(), 0);

    for (size_t i = 0; i < rows.size(); i++) {
        if (auto const * p = rows[i]) {
            auto const * q = p;
            for (; *q != END_OF_ROW; q++)
                ;
            w[i] = q - p;
        }
    }
    print_stats_generic(out, w, "rows", "weight", verbose);
}

template <typename F>
static void apply_fragment_multithread(size_t nrows, size_t nthreads, F const & f)
{
    std::vector<std::thread> T;
    T.reserve(nthreads);
    auto const S = subdivision(nrows, nthreads);
    for (size_t i = 0; i < nthreads; i++) {
        auto [i0, i1] = S.nth_block(i);
        T.emplace_back(f, i0, i1);
    }
    for (auto & t: T)
        t.join();
}

auto purge_matrix::compute_sum2(size_t nthreads) const -> sum2_type
{
    sum2_type res(column_weights.size(), 0);
    auto fragment = [&](size_t i0, size_t i1) {
        for (size_t i = i0; i < i1; i++) {
            if (auto const * p = rows[i]) {
                for (; *p != END_OF_ROW; p++) {
                    if (LIKELY(column_weights[*p] == 2))
                        std::atomic_ref(res[*p]) += i;
                }
            }
        }
    };
    apply_fragment_multithread(rows.size(), nthreads, fragment);
    return res;
}

void purge_matrix::delete_row(size_t i)
{
    ASSERT_ALWAYS(is_active(i));
    auto const * p = rows[i];
    size_t killed = 0;
    for (; *p != std::numeric_limits<index_t>::max(); p++) {
        auto w = std::atomic_ref(column_weights[*p]);
        ASSERT(w.load());
        /* Decrease only if not equal to the maximum value */
        /* If weight becomes 0, we just remove a column */
        if (w.load() < OVERWEIGHT)
            if (!--w)
                killed++;
    }
    rows[i] = nullptr;
    if (killed)
        std::atomic_ref(remaining_columns) -= killed;
    std::atomic_ref(remaining_rows)--;
}

/* returns the number of rows that are removed */
size_t purge_matrix::one_singleton_removal(size_t nthreads)
{
    const size_t rows_before = remaining_rows;
    auto fragment = [&](size_t i0, size_t i1) {
        for (size_t i = i0; i < i1; i++) {
            auto const * p = rows[i];
            if (p) {
                for (auto const * q = p; *q != END_OF_ROW; q++) {
                    if (std::atomic_ref(column_weights[*q]).load() == 1) {
                        delete_row(i);
                        break;
                    }
                }
            }
        }
    };
    apply_fragment_multithread(rows.size(), nthreads, fragment);
    return rows_before - remaining_rows;
}

/* Perform a complete singleton removal step:  call singleton_removal_oneiter_*
 * until there is no more singleton.
 * Return the excess at the end *
 */
int64_t purge_matrix::singleton_removal(size_t nthreads, int verbose)
{
    int64_t excess;
    unsigned int iter = 0;
    for (;;) {
        excess = remaining_rows - remaining_columns;
        if (verbose >= 0) {
            std::string const s =
                iter ? fmt::format("  iter {:03}", iter) : "begin with";
            fmt::print("Sing. rem.: {}:"
                       " nrows={} ncols={} excess={} at {:2.2f}\n",
                       s, remaining_rows, remaining_columns, excess, seconds());
            fflush(stdout);
        }

        if (iter++ ; !one_singleton_removal(nthreads))
            break;
    }

    if (verbose >= 0)
        fmt::print("Sing. rem.:   iter {:03}:"
                   " No more singletons, finished at {:2.2f}\n",
                   iter, seconds());

    return excess;
}
