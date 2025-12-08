#include "cado.h" // IWYU pragma: keep

#include <cstdio>

#include <map>
#include <string>
#include <mutex>
#include <ranges>
#include <utility>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "las-auxiliary-data.hpp"

#include "macros.h"

#include "las-config.hpp"
#include "las-info.hpp"
#include "tdict.hpp"
#include "arith/ularith.h"
#include "verbose.h"

#include "fmt/ranges.h"

void
sieve_checksum::update(const unsigned int other)
{
    unsigned long r;
    ularith_addmod_ul_ul(&r, checksum, other, checksum_prime);
    checksum = r;
}

void
sieve_checksum::update(const unsigned char *data, const size_t len)
{
    cxx_mpz mb;
    unsigned int new_checksum;

    mpz_import(mb, len, -1, sizeof(unsigned char), -1, 0, data);
    new_checksum = mpz_tdiv_ui(mb, checksum_prime);

    this->update(new_checksum);
}

nfs_aux::~nfs_aux()
{
    ASSERT_ALWAYS_NOTHROW(!rt.timer.running());

    ASSERT_ALWAYS_NOTHROW(dest_rt);

    if (!complete) {
        std::lock_guard<std::mutex> const lock(rt.mm);
        dest_rt->rep.accumulate_and_clear(std::move(rt.rep));
        dest_rt->timer += rt.timer;
        return;
    }

    /* locking is not needed here, since rt.rep is our private timer and all
     * threads completed their subtasks.
     */
    for (auto & T : th) {
        rt.rep.accumulate_and_clear(std::move(T.rep));
        rt.timer += T.timer;
        for (int side = 0; side < (int) T.checksum_post_sieve.size(); side++)
            checksum_post_sieve[side].update(T.checksum_post_sieve[side]);
    }

    verbose_output_start_batch();

    if (tdict::global_enable >= 2) {
        verbose_fmt_print (0, 1, "{}", rt.timer.display());

        timetree_t::timer_data_type t = 0;
        for(auto const &c : rt.timer.filter_by_category()) {
            verbose_fmt_print (0, 1, "# {}: {:.2f}\n",
                    coarse_las_timers::explain(c.first), c.second);
            t += c.second;
        }
        verbose_fmt_print (0, 1, "# total counted time: {:.2f}\n", t);
    }

    // we're in a dtor, exceptions can turn your computer into a coconut.
    // Well, we do have ASSERT_ALWAYS down below...
    // coverity[fun_call_w_exception]
    rt.rep.display_survivor_counters();

    auto p = [](auto const & ck) { return ck.get_checksum(); };
    verbose_fmt_print(0, 2,
            "# Checksums over sieve region: after all sieving: {}\n",
            fmt::join(checksum_post_sieve | std::views::transform(p), " "));

    verbose_fmt_print(0, 1, "# {} {}\n",
            rt.rep.reports, las.batch ? "survivor(s) saved" : "relation(s)");

    if (las.suppress_duplicates)
        verbose_fmt_print(0, 1, "# number of eliminated duplicates: {}\n", rt.rep.duplicates);

    wct_qt0 = wct_seconds() - wct_qt0;

    if (rt.rep.survivors.after_sieve != rt.rep.survivors.not_both_even) {
        verbose_fmt_print(0, 1, "# Warning: found {} hits with i,j both even (not a bug, but should be very rare)\n", rt.rep.survivors.after_sieve - rt.rep.survivors.not_both_even);
    }

    auto D = rt.timer.filter_by_category();
    timetree_t::timer_data_type const qtcpu = rt.timer.total_counted_time();

    verbose_fmt_print (0, 1, "# Time for {}: ({:.1f} elapsed)",
            doing.sq(), wct_qt0);
    if (las_production_mode) {
        verbose_fmt_print (0, 1, " [-production mode, no timings]");
    } else {
        verbose_fmt_print (0, 1, " {:1.2f}s", qtcpu);
        verbose_fmt_print (0, 2,
                " [norm {:1.2f}+{:1.1f}, sieving {:1.1f}"
                " ({:1.1f}+{:1.1f} + {:1.1f}),"
                " factor {:1.1f} ({:1.1f}+{:1.1f} + {:1.1f}),"
                " rest {:1.1f}]",
                D[coarse_las_timers::norms(0)],
                D[coarse_las_timers::norms(1)],

                D[coarse_las_timers::sieving(0)]+
                D[coarse_las_timers::sieving(1)]+
                D[coarse_las_timers::search_survivors()]+
                D[coarse_las_timers::sieving_mixed()],

                D[coarse_las_timers::sieving(0)],
                D[coarse_las_timers::sieving(1)],
                D[coarse_las_timers::search_survivors()]+
                D[coarse_las_timers::sieving_mixed()],

                D[coarse_las_timers::cofactoring(0)]+
                D[coarse_las_timers::cofactoring(1)]+
                D[coarse_las_timers::cofactoring_mixed()],

                D[coarse_las_timers::cofactoring(0)],
                D[coarse_las_timers::cofactoring(1)],
                D[coarse_las_timers::cofactoring_mixed()],

                D[coarse_las_timers::bookkeeping()]
                    );
    }
    verbose_fmt_print (0, 1, "\n");

    verbose_output_end_batch();

    {
        std::lock_guard<std::mutex> const lock(dest_rt->mm);
        dest_rt->rep.accumulate_and_clear(std::move(rt.rep));
        dest_rt->timer += rt.timer;
    }
}

#ifndef DISABLE_TIMINGS
/* we really wish to have a single timing slot for all the instantiations
 * of fill_in_buckets_toplevel_wrapper */
tdict::slot tdict_slot_for_fibt("fill_in_buckets_toplevel");
tdict::slot tdict_slot_for_alloc_buckets("allocate_buckets");
tdict::slot tdict_slot_for_threads("multithreaded tasks");
tdict::slot_parametric tdict_slot_for_side("side ", "");
#endif

