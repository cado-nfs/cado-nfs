#include "cado.h" // IWYU pragma: keep

#include <cmath>

#include <ostream>
#include <string>

#include <gmp.h>
#include "fmt/format.h"

#include "las-globals.hpp"
#include "las-info.hpp"
#include "las-descent-candidate-relation.hpp"
#include "las-special-q-task-tree.hpp"
#include "las-siever-config.hpp"
#include "relation.hpp"
#include "timing.h"
#include "typedefs.h"
#include "verbose.h"

static std::string shortcode(special_q_task_tree::status_code const & s)
{
    switch(s) {
        case special_q_task_tree::PENDING:   return "?";
        case special_q_task_tree::DONE:      return ".";
        case special_q_task_tree::ABANDONED: return "#";
        case special_q_task_tree::IN_PROGRESS:  return "/";
        case special_q_task_tree::IN_RECURSION: return "*";
    }
    /* unreachable, but some compilers are stupid. */
    return {};
}

std::ostream& operator<<(std::ostream& os, special_q_task_tree::prefixed const & T)
{
    special_q_task_tree const * t = T.t;
    std::string pre = T.prefix;

    if (t->sq()) {
        pre += shortcode(t->status);
        std::string comment;
        if (t->try_again)
            comment += fmt::format(" try{}", t->try_again);
        os << fmt::format("{} {} {:1.4f}{}\t{}\n",
                pre,
                t->shortname(), t->spent, comment, t->sq());
    } else {
        os << pre << "[ROOT NODE]\n";
    }

    for(auto const & c : t->children_by_status) {
        for(auto const * t : c.second)
            os << special_q_task_tree::prefixed { t, pre };
    };

    return os;
}

std::ostream& operator<<(std::ostream& os, special_q_task_tree::status_code const & s)
{
    switch(s) {
        case special_q_task_tree::PENDING: os << "PENDING"; break;
        case special_q_task_tree::DONE: os << "DONE"; break;
        case special_q_task_tree::ABANDONED: os << "ABANDONED"; break;
        case special_q_task_tree::IN_PROGRESS: os << "IN_PROGRESS"; break;
        case special_q_task_tree::IN_RECURSION: os << "IN_RECURSION"; break;
    }
    return os;
}

bool special_q_task_tree::new_candidate_relation(las_info const & las, relation & rel, std::mutex & mm)
{
    /* This returns true only if this descent node is now done, either based
     * on the new relation we have registered, or because the previous
     * relation is better anyway */

    if (las.tree->must_avoid(rel)) {
        verbose_fmt_print(0, 1, "# [descent] Warning: we have already used this relation, avoiding\n");
        return true;
    }

    /* compute rho for all primes, even on the rational side */
    rel.fixup_r(true);

    descent_candidate_relation newcomer;
    newcomer.rel = rel;
    double time_left = 0;

    for(int side = 0 ; side < las.cpoly->nb_polys ; side++) {
        for(unsigned int i = 0 ; i < rel.sides[side].size() ; i++) {
            special_q v(side, rel.sides[side][i]);
            if (mpz_cmp(p, v.p) == 0)
                continue;

            p_r_values_t const p = mpz_get_ui(v.p);
            if (v.p.fits<p_r_values_t>()) {
                p_r_values_t const r = mpz_get_ui(v.r);
                if (las.dlog_base.is_known(side, p, r))
                    continue;
            }

            siever_config_pool::key_type const K(side, mpz_sizeinbase(v.p, 2));
            double const e = las.config_pool.hint_expected_time(K);
            if (e < 0) {
                /* This is not worrisome per se. We just do
                 * not have the info in the descent hint table,
                 * period.
                 */
                verbose_fmt_print (0, 1,
                            "# [descent] Warning: cannot estimate"
                            " refactoring time for relation involving {} ({})\n",
                            v.shortname(), v);
                time_left = INFINITY;
            } else {
                if (std::isfinite(time_left))
                    time_left += e;
            }
            newcomer.outstanding.emplace_back(v);
        }
    }
    verbose_fmt_print(0, 1, "# [descent] This relation entails an additional time of {:.2f} for the smoothing process ({} children)\n",
            time_left, newcomer.outstanding.size());

    /* when we're re-examining this special-q because of a previous
     * failure, there's absolutely no reason to hurry up on a relation */
    newcomer.set_time_left(time_left, try_again ? INFINITY : general_grace_time_ratio);

    {
        const std::lock_guard<std::mutex> lock(mm);
        if (newcomer < contender) {
            if (newcomer.outstanding.empty()) {
                verbose_fmt_print(0, 1, "# [descent] Yiippee, splitting done\n");
            } else if (std::isfinite(contender.deadline)) {
                // This implies that newcomer.deadline is also finite 
                const double delta = contender.time_left-newcomer.time_left;
                verbose_fmt_print(0, 1, "# [descent] Improved ETA by {:.2f}\n", delta);
            } else if (contender) {
                // This implies that we have fewer outstanding special-q's
                verbose_fmt_print(0, 1, "# [descent] Improved number of children to split from {} to {}\n",
                        contender.outstanding.size(),
                        newcomer.outstanding.size());
            }
            contender = newcomer;
            if (!contender.outstanding.empty()) {
                verbose_fmt_print(0, 1, "# [descent] still searching for {:.2f}\n", contender.deadline - seconds());
            }
        }
    }
    return must_take_decision();
}

