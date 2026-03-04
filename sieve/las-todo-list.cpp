#include "cado.h" // IWYU pragma: keep

#include <cctype>
#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <thread>
#include <vector>
#include <mutex>
#include <string>
#include <stdexcept>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "cado_poly.h"
#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "las-galois.hpp"
#include "special-q.hpp"
#include "las-todo-list.hpp"
#include "macros.h"
#include "mpz_poly.h"
#include "params.h"
#include "rootfinder.h"
#include "verbose.h"

void todo_list_base::configure_switches(cxx_param_list & pl)
{
    param_list_configure_switch(pl, "-allow-compsq", nullptr);
    param_list_configure_switch(pl, "-print-todo-list", nullptr);
}

void todo_list_base::declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "sqside", "put special-q on this side");
    param_list_decl_usage(pl, "seed", "Use this seed for random state seeding (currently used only by --random-sample)");
    param_list_decl_usage(pl, "random-sample", "Sample this number of special-q's at random, within the range [q0,q1[");
    param_list_decl_usage(pl, "nq", "Process this number of special-q's and stop");
    param_list_decl_usage(pl, "todo", "provide file with a list of special-q to sieve instead of qrange");
    param_list_decl_usage(pl, "allow-compsq", "allows composite special-q");
    param_list_decl_usage(pl, "qfac-min", "factors of q must be at least that");
    param_list_decl_usage(pl, "qfac-max", "factors of q must be at most that");
    param_list_decl_usage(pl, "print-todo-list", "only print the special-q's to be sieved");
}

todo_list_base::todo_list_base(cxx_cado_poly const & cpoly, cxx_param_list & pl)
    : cpoly(cpoly), galois(param_list_lookup_string(pl, "galois"))
{
    unsigned long seed = 0;
    if (param_list_parse_ulong(pl, "seed", &seed))
        gmp_randseed_ui(rstate, seed);

    if (param_list_parse(pl, "random-sample", nq_max)) {
        random_sampling = 1;
        if (param_list_parse(pl, "nq", nq_max)) {
            fmt::print(stderr, "# Warning: both options -nq and -random-sample "
                               "found. Limiting the number of generated primes "
                               "to {}\n", nq_max);
        }
    } else {
        param_list_parse(pl, "nq", nq_max);
    }

    sqside = cpoly->nb_polys == 1 ? 0 : 1;
    if (!param_list_parse_int(pl, "sqside", &sqside) && cpoly->nb_polys > 1) {
        verbose_fmt_print(0, 1, "# Warning: sqside not given, "
                "assuming side 1 for backward compatibility.\n");
    }
    ASSERT_ALWAYS(sqside >= 0 && sqside < cpoly->nb_polys);

    /* Init and parse info regarding work to be done by the siever */
    /* Actual parsing of the command-line fragments is done within
     * las_todo_feed, but this is an admittedly contrived way to work */
    const char * filename = param_list_lookup_string(pl, "todo");
    if (filename) {
        todo_list_fd = std::make_unique<std::ifstream>(filename);
        if (todo_list_fd->fail()) {
            fprintf(stderr, "%s: %s\n", filename, strerror(errno));
            /* There's no point in proceeding, since it would really change
             * the behaviour of the program to do so */
            exit(EXIT_FAILURE);
        }
    }

    /* composite special-q ? */
    if ((allow_composite_q = param_list_parse_switch(pl, "-allow-compsq"))) {
        /* defaults are set in the class description */
        param_list_parse_uint64(pl, "qfac-min", &qfac_min);
        param_list_parse_uint64(pl, "qfac-max", &qfac_max);
    }

    print_todo_list_flag = param_list_parse_switch(pl, "-print-todo-list");
}


/* {{{ Populating the todo list */

/* Format of a file with a list of special-q (-todo option):
 *   - Comments are allowed (start line with #)
 *   - Blank lines are ignored
 *   - Each valid line must have the form
 *       s q r
 *     where s is the side (0 or 1) of the special q, and q and r are as usual.
 */
bool todo_list_base::feed_qlist()
{
    /* XXX this is called with a lock held on this->mm, so care must be
     * taken to call this->push_unlocked() and not this->push() !!
     */
    if (!super::empty())
        return true;

    if (created == nq_max)
        return false;

    std::string line;
    auto is_comment = [](std::string const & s) {
        for(auto c : s) {
            if (c == '#') return true;
            if (!isspace(c)) return false;
        }
        return true;
    };

    for( ; std::getline(*todo_list_fd, line) && is_comment(line) ; ) ;

    if (!*todo_list_fd)
        return false;

    /* We have a new entry to parse */
    cxx_mpz p, r = -1;
    int side = -1;
    std::istringstream is(line);
    std::string tail;
    is >> side >> p;
    if (!is.eof()) {
        is >> r;
    }
    if (!is || (!is.eof() && (is >> tail, !is_comment(tail))))
        throw std::runtime_error(
                fmt::format("parse error in todo file while reading {}", line));
    auto const & f = cpoly->pols[side];
    /* specifying the rational root as <0
     * means that it must be recomputed. Putting 0 does not have this
     * effect, since it is a legitimate value after all.
     */
    if (r < 0) {
        // For rational side, we can compute the root easily.
        ASSERT_ALWAYS(f->deg == 1);
        cxx_gmp_randstate rstate;
        std::vector<cxx_mpz> roots = mpz_poly_roots(f, p, rstate);
        ASSERT_ALWAYS(roots.size() == 1);
        r = roots[0];
    }

    ASSERT_ALWAYS(p > 0);
    ASSERT_ALWAYS(r >= 0);

    push_unlocked(special_q(p, r, side));
    return true;
}


/* This exists because of the race condition between feed() and pop()
 */
special_q todo_list_base::feed_and_pop()
{
    const std::lock_guard<std::mutex> foo(mm);

    if (super::empty()) {
        if (todo_list_fd)
            feed_qlist();
        else
            feed_qrange(rstate);
    }
    if (super::empty())
        return {};
    auto ret = super::top();
    // if (bool(ret)) pulled++;
    super::pop();
    return ret;
}
/* }}} */

void todo_list_base::print_todo_list(cxx_param_list const & pl, int nthreads) const
{

    if (random_sampling) {
        /* Then we cannot print the todo list in a multithreaded way
         * without breaking the randomness predictability. It's a bit
         * annoying, yes. On the other hand, it's highly likely in this
         * case that the number of q's to pick is small enough anyway !
         */
        nthreads = 1;
    }

    std::vector<std::vector<special_q>> lists(nthreads);

    auto segment = [&, this](unsigned int i) {
        auto todo2 = this->create_sub_todo_list(i, nthreads, pl);
        for(;;) {
            auto doing = todo2->feed_and_pop();
            if (!doing) break;
            if (nthreads == 1) {
                verbose_fmt_print(0, 1, "{} {} {}\n",
                        doing.side, doing.p, doing.r);
            } else {
                lists[i].push_back(doing);
            }
        }
    };

    if (nthreads > 1) {
        verbose_fmt_print(0, 1, "# Collecting the todo list in memory from {} "
                                "using {} threads\n", *this, nthreads);
    }

    std::vector<std::thread> subjobs;
    subjobs.reserve(nthreads);
    for(int subjob = 0 ; subjob < nthreads ; ++subjob)
        subjobs.emplace_back(segment, subjob);
    for(auto & t : subjobs) t.join();
    for(auto const & v : lists) {
        for(auto const & doing : v) {
            verbose_fmt_print(0, 1, "{} {} {}\n",
                    doing.side, doing.p, doing.r);
        }
    }
}


void las_todo_list::declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "q0",   "left bound of special-q range");
    param_list_decl_usage(pl, "q1",   "right bound of special-q range");
    param_list_decl_usage(pl, "rho",  "sieve only root r mod q0");
    todo_list_base::declare_usage(pl);
}

las_todo_list::las_todo_list(cxx_cado_poly const & cpoly, cxx_param_list & pl)
    : todo_list_base(cpoly, pl)
{
    if (allow_composite_q && galois) {
        fprintf(stderr, "-galois and -allow-compsq are incompatible options "
                        "at the moment\n");
        exit(EXIT_FAILURE);
    }

    if (nq_max != SIZE_MAX) {
        if (param_list_lookup_string(pl, "rho")) {
            fprintf(stderr, "Error: argument -nq is incompatible with -rho\n");
            exit(EXIT_FAILURE);
        }
        if (param_list_lookup_string(pl, "q1"))
            verbose_fmt_print(0, 1, "# Warning: arguments nq and q1 will both "
                                    "limit the q range\n");
    }

    /* It's not forbidden to miss -q0 */
    param_list_parse_mpz(pl, "q0", q0);
    param_list_parse_mpz(pl, "q1", q1);

    if (mpz_cmp_ui(q0, 0) == 0) {
        if (!todo_list_fd) {
            fprintf(stderr, "Error: Need either -todo or -q0\n");
            exit(EXIT_FAILURE);
        }
        return;
    }

    cxx_gmp_randstate rstate;

    if (mpz_cmp_ui(q1, 0) != 0) {
        next_legitimate_specialq(q0, q0, 0);
    } else {
        /* We don't have -q1. If we have -rho, we sieve only <q0,
         * rho>. */
        cxx_mpz rho;
        if (param_list_parse(pl, "rho", rho)) {
            cxx_mpz q0_cmdline = q0;
            auto fac_q = next_legitimate_specialq(q0, q0, 0);
            if (mpz_cmp(q0, q0_cmdline) != 0) {
                fprintf(stderr, "Error: q0 is not a legitimate special-q\n");
                exit(EXIT_FAILURE);
            }
            std::vector<cxx_mpz> roots = mpz_poly_roots(cpoly->pols[sqside], q0, fac_q, rstate);
            if (std::ranges::find(roots, rho) == roots.end()) {
                fprintf(stderr, "Error: rho is not a root modulo q0\n");
                exit(EXIT_FAILURE);
            }
            push_unlocked(special_q(q0, rho, sqside));
            /* Set empty interval [q0 + 1, q0] as special-q interval */
            mpz_set(q1, q0);
            mpz_add_ui (q0, q0, 1);
        } else {
            /* If we don't have -rho, we sieve only q0, but all roots of it.
               If -q0 does not give a legitimate special-q value, advance to the
               next legitimate one. */
            mpz_set(rho, q0);
            next_legitimate_specialq(q0, q0, 0);
            mpz_set(q1, q0);
        }
    }

    if (random_sampling) {
        if (mpz_cmp_ui(q0, 0) == 0 || mpz_cmp_ui(q1, 0) == 0) {
            fprintf(stderr, "Error: --random-sample requires -q0 and -q1\n");
            exit(EXIT_FAILURE);
        }
        /* For random sampling, it's important that for all integers in
         * the range [q0, q1[, their nextprime() is within the range, and
         * that at least one such has roots mod f. Make sure that
         * this is the case.
         */
        cxx_mpz q, q1_orig = q1;
        /* we need to know the limit of the q range */
        for(unsigned long i = 1 ; ; i++) {
            mpz_sub_ui(q, q1, i);
            const std::vector<uint64_t> fac_q = next_legitimate_specialq(q, q, 0);
            if (mpz_cmp(q, q1) >= 0)
                continue;
            if (!mpz_poly_roots(cpoly->pols[sqside], q, fac_q, rstate).empty())
                break;
            /* small optimization: avoid redoing root finding
             * several times */
            mpz_set (q1, q);
            i = 1;
        }
        /* now q is the largest prime < q1 with f having roots mod q */
        mpz_add_ui (q1, q, 1);

        /* so now if we pick x an integer in [q0, q1[, then nextprime(x-1)
         * will be in [q0, q1_orig[, which is what we look for,
         * really.
         */
        if (mpz_cmp(q0, q1) > 0) {
            fmt::print(stderr, 
                    "Error: range [{},{}[ contains no prime with roots mod f\n",
                    q0,
                    q1_orig);
            exit(EXIT_FAILURE);
        }
    }
}

/* Put in r the smallest legitimate special-q value that it at least
 * s + diff (note that if s+diff is already legitimate, then r = s+diff
 * will result.
 * In case of composite sq, also returns the factorization of r
 */
std::vector<uint64_t> las_todo_list::next_legitimate_specialq(cxx_mpz & r, cxx_mpz const & s, const unsigned long diff) const
{
    std::vector<uint64_t> fac;
    if (allow_composite_q) {
        unsigned long tfac[64];
        int const nf = next_mpz_with_factor_constraints(r, tfac,
                s, diff, qfac_min, qfac_max);
        fac.assign(tfac, tfac + nf);
    } else {
        mpz_add_ui(r, s, diff);
        /* mpz_nextprime() returns a prime *greater than* its input argument,
           which we don't always want, so we subtract 1 first. */
        mpz_sub_ui(r, r, 1);
        mpz_nextprime(r, r);
    }
    return fac;
}

/* See below in main() for documentation about the q-range and q-list
 * modes */
/* These functions return non-zero if the todo list is not empty.
 * Note: contrary to the qlist mode, here the q-range will be pushed at
 * once (but the caller doesn't need to know that).
 *
 * Note: the random state is used by both the random sampler **AND** the
 * rootfinder.
 */
bool las_todo_list::feed_qrange(gmp_randstate_t rstate)
{
    /* XXX this is called with a lock held on this->mm, so care must be
     * taken to call this->push_unlocked() and not this->push() !!
     */

    /* If we still have entries in the stack, don't add more now */
    if (!super::empty())
        return true;

    mpz_poly_ptr f = cpoly->pols[sqside];

    if (!random_sampling) {
        /* We're going to process the sq's and put them into the list
           The loop processes all special-q in [q0, q1]. On loop entry,
           the value in q0 is known to be a legitimate special-q. Its
           factorization is lost, so we recompute it. */

        /* handy aliases */
        cxx_mpz & q = q0;
        auto fac_q = next_legitimate_specialq(q, q, 0);

        std::vector<special_q> my_list;

        int nb_no_roots = 0;
        int nb_rootfinding = 0;
        /* If nq_max is specified, then q1 has no effect, even though it
         * has been set equal to q */
        for ( ; (nq_max < SIZE_MAX || mpz_cmp(q, q1) < 0) &&
                created + my_list.size() < nq_max ; )
        {
            std::vector<cxx_mpz> roots = mpz_poly_roots(f, q, fac_q, rstate);

            nb_rootfinding++;
            if (roots.empty()) nb_no_roots++;

            if (galois) {
                skip_galois_roots(q, roots, galois);
            }

            for (auto const & r : roots)
                my_list.emplace_back(q, r, sqside);

            fac_q = next_legitimate_specialq(q, q, 1);
        }

        if (nb_no_roots) {
            verbose_fmt_print(0, 1, "# polynomial has no roots for {} of the {} primes that were tried\n", nb_no_roots, nb_rootfinding);
        }

        // Truncate to nq_max if necessary and push the sq in reverse
        // order, because they are processed via a stack (required for
        // the descent).
        size_t push_here = my_list.size();
        if (nq_max < SIZE_MAX)
            push_here = std::min(push_here, nq_max - created);
        for(size_t i = 0 ; i < push_here ; i++) {
            size_t const ind = push_here-i-1;
            push_unlocked(my_list[ind]);
        }
    } else { /* random sampling case */
        /* we care about being uniform here */
        cxx_mpz q;
        cxx_mpz diff;
        mpz_sub(diff, q1, q0);
        ASSERT_ALWAYS(created == 0 || created == nq_max);
        unsigned int spin = 0;
        for ( ; created < nq_max ; ) {
            /* try in [q0 + k * (q1-q0) / n, q0 + (k+1) * (q1-q0) / n[ */
            cxx_mpz q0l, q1l;
	    /* we use k = n-1-created instead of k=created so that
	       special-q's are sieved in increasing order */
	    size_t const k = nq_max - 1 - created;
            mpz_mul_ui(q0l, diff, k);
            mpz_mul_ui(q1l, diff, k + 1);
            mpz_fdiv_q_ui(q0l, q0l, nq_max);
            mpz_fdiv_q_ui(q1l, q1l, nq_max);
            mpz_add(q0l, q0, q0l);
            mpz_add(q1l, q0, q1l);

            mpz_sub(q, q1l, q0l);
            mpz_urandomm(q, rstate, q);
            mpz_add(q, q, q0l);
            auto const fac_q = next_legitimate_specialq(q, q, 0);
            std::vector<cxx_mpz> roots = mpz_poly_roots(f, q, fac_q, rstate);
            if (roots.empty()) {
                spin++;
                if (spin >= 1000) {
                    fprintf(stderr, "Error: cannot find primes with roots in one of the sub-ranges for random sampling\n");
                    exit(EXIT_FAILURE);
                }
                continue;
            }
            spin = 0;
            if (galois) {
                skip_galois_roots(q, roots, galois);
            }
            push_unlocked(special_q(q,
                        roots[gmp_urandomm_ui(rstate, roots.size())],
                        sqside));
        }
    }

    return !super::empty();
}

std::unique_ptr<todo_list_base> las_todo_list::create_sub_todo_list(
    unsigned int i, int nthreads, cxx_param_list const & pl) const
{
    cxx_param_list pl2 = pl;
    cxx_mpz tmp;
    mpz_sub(tmp, q1, q0);
    mpz_mul_ui(tmp, tmp, i);
    mpz_fdiv_q_ui(tmp, tmp, nthreads);
    mpz_add(tmp, q0, tmp);
    {
        std::ostringstream os;
        os << tmp;
        param_list_add_key(pl2, "q0", os.str().c_str(), PARAMETER_FROM_CMDLINE);
    }

    mpz_sub(tmp, q1, q0);
    mpz_mul_ui(tmp, tmp, i + 1);
    mpz_fdiv_q_ui(tmp, tmp, nthreads);
    mpz_add(tmp, q0, tmp);
    {
        std::ostringstream os;
        os << tmp;
        param_list_add_key(pl2, "q1", os.str().c_str(), PARAMETER_FROM_CMDLINE);
    }

    return todo_list_base::create<las_todo_list>(cpoly, pl2);
}

std::string las_todo_list::string_repr() const
{
    return fmt::format("q0={} to q1={}", q0, q1);
}

/******************************************************************************/
/* SIQS ***********************************************************************/
/******************************************************************************/
void siqs_todo_list::declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "qidx0",   "left bound of special-q range");
    param_list_decl_usage(pl, "qidx1",   "right bound of special-q range");
    param_list_decl_usage(pl, "qfac-nfac", "number of factors of q");
    todo_list_base::declare_usage(pl);
}

siqs_todo_list::siqs_todo_list(cxx_cado_poly const & cpoly, cxx_param_list & pl)
    : todo_list_base(cpoly, pl), qidx0(-1), qidx1(-1)
{
    if (!allow_composite_q) {
        verbose_fmt_print(0, 1, "# Warning, -allow-compsq will be assumed\n");
        allow_composite_q = true;
    }

    /* Needed here, because was not parsed if allow_composite_q was false. */
    param_list_parse_uint64(pl, "qfac-min", &qfac_min);
    param_list_parse_uint64(pl, "qfac-max", &qfac_max);

    if (qfac_min <= 2) {
        fmt::print(stderr, "# Error, -qfac-min must be > 2\n");
        exit(EXIT_FAILURE);
    }

    if (!param_list_parse_uint64(pl, "qfac-nfac", &qfac_nfac)) {
        fmt::print(stderr, "# Error, -qfac-nfac is mandatory\n");
        exit(EXIT_FAILURE);
    }

    /* It's not forbidden to miss -q0 */
    param_list_parse_mpz(pl, "qidx0", qidx0);
    param_list_parse_mpz(pl, "qidx1", qidx1);

    if (mpz_cmp_si(qidx0, -1) == 0) {
        if (!todo_list_fd) {
            fprintf(stderr, "Error: Need either -todo or -qidx0\n");
            exit(EXIT_FAILURE);
        }
        return;
    }

    if (param_list_lookup_string(pl, "qidx1") != nullptr
            && param_list_lookup_string(pl, "nq") != nullptr) {
        /* if both nq and qidx1 were given, we should have qidx1=qidx0+nq */
        cxx_mpz tmp;
        mpz_add_ui(tmp, qidx0, nq_max);
        if (tmp != qidx1) {
            fmt::print(stderr, "Error: incompatible '-nq {}' / '-qidx1 {}' for "
                               "qidx0={}\n", nq_max, qidx1, qidx0);
            exit(EXIT_FAILURE);
        }
    }

    cxx_gmp_randstate rstate;

    if (mpz_cmp_si(qidx1, -1) == 0) {
        if (nq_max < SIZE_MAX) {
            mpz_add_ui(qidx1, qidx0, nq_max);
        } else {
            /* We don't have -q1, we sieve only q0 */
            mpz_add_ui(qidx1, qidx0, 1U);
        }
    }

    if (random_sampling) {
        if (mpz_cmp_si(qidx0, -1) == 0 || mpz_cmp_si(qidx1, -1) == 0) {
            fprintf(stderr, "Error: --random-sample requires -qidx0 and -qidx1\n");
            exit(EXIT_FAILURE);
        }
    }
}

bool siqs_todo_list::feed_qrange(gmp_randstate_t rstate)
{
    /* If we still have entries in the stack, don't add more now */
    if (!super::empty())
        return true;

    std::vector<uint64_t> scratch_vec;
    cxx_mpz scratch_mpz;

    if (!random_sampling) {
        cxx_mpz idx;
        unsigned long nadded = 0;
        // Push the sq in reverse order, because they are processed via a stack
        if (nq_max < UINT_MAX) {
            idx = qidx0 + ((int)nq_max - (int)created - 1);
        } else {
            idx = qidx1 - 1U;
        }
        for( ; idx >= qidx0; idx-=1U, ++nadded) {
            push_unlocked(special_q_from_index(idx, scratch_mpz, scratch_vec));
        }
        qidx0 += nadded;
    } else { /* random sampling case */
        /* we care about being uniform here */
        cxx_mpz qidx;
        cxx_mpz diff;
        mpz_sub(diff, qidx1, qidx0);
        ASSERT_ALWAYS(created == 0 || created == nq_max);
        unsigned long const n = nq_max;
        for ( ; created < n ; ) {
            /* try in [qidx0 + k * (qidx1-qidx0) / n, qidx0 + (k+1) * (qidx1-qidx0) / n[ */
            cxx_mpz qidx0l, qidx1l;
            /* we use k = n-1-nq_pushed instead of k=nq_pushed so that
             * special-q's are sieved in increasing order.
             */
            unsigned long const k = n - 1 - created;
            mpz_mul_ui(qidx0l, diff, k);
            mpz_mul_ui(qidx1l, diff, k + 1);
            mpz_fdiv_q_ui(qidx0l, qidx0l, n);
            mpz_fdiv_q_ui(qidx1l, qidx1l, n);
            mpz_add(qidx0l, qidx0, qidx0l);
            mpz_add(qidx1l, qidx0, qidx1l);

            mpz_sub(qidx, qidx1l, qidx0l);
            mpz_urandomm(qidx, rstate, qidx);
            mpz_add(qidx, qidx, qidx0l);
            push_unlocked(special_q_from_index(qidx, scratch_mpz, scratch_vec));
        }
    }

    return !super::empty();
}

void siqs_todo_list::k_combination_from_index(
        std::vector<uint64_t>::reverse_iterator it, cxx_mpz & idx, uint64_t k)
{
    if (k == 0) {
        ASSERT_ALWAYS(idx == 0U);
    } else if (idx == 0U) {
        *it = k-1;
        k_combination_from_index(++it, idx, k-1);
    } else {
        cxx_mpz t;
        do {
            ++(*it);
            mpz_bin_uiui(t, *it, k);
        } while (t <= idx);
        --(*it);
        mpz_bin_uiui(t, *it, k);
        idx -= t;
        k_combination_from_index(++it, idx, k-1);
    }
}

special_q siqs_todo_list::special_q_from_index(cxx_mpz const & idx,
                                          cxx_mpz & scratch_q,
                                          std::vector<uint64_t> & scratch_vec) {
    scratch_vec.assign(qfac_nfac, 0U);
    scratch_q = idx;
    k_combination_from_index(scratch_vec.rbegin(), scratch_q, qfac_nfac);
    /* scratch_vec now contains the index of the factors, now the actual factor
     */
    scratch_q = 1U;
    for (auto & f: scratch_vec) {
        f = get_qfac_prime(f);
        mpz_mul_ui(scratch_q, scratch_q, f);
    }
    return { scratch_q, -1, sqside, scratch_vec };
}

uint64_t siqs_todo_list::get_qfac_prime(uint64_t i) {
    cxx_mpz D;
    mpz_neg(D, mpz_poly_coeff_const(cpoly->pols[0], 0));
    while (i >= qfac_primes.size()) {
        cxx_mpz p;
        if (qfac_primes.empty()) {
            p = qfac_min-1U; /* -1 in case qfac_min is prime */
        } else {
            p = qfac_primes.back();
        }
        do {
            mpz_nextprime(p, p);
        } while (mpz_legendre(D, p) != 1);
        ASSERT_ALWAYS(p <= qfac_max);
        qfac_primes.emplace_back(std::move(p));
    }
    return mpz_get_ui(qfac_primes[i]);
}

std::unique_ptr<todo_list_base> siqs_todo_list::create_sub_todo_list(
    unsigned int i, int nthreads, cxx_param_list const & pl) const
{
    cxx_param_list pl2 = pl;
    cxx_mpz tmp;
    mpz_sub(tmp, qidx1, qidx0);
    mpz_mul_ui(tmp, tmp, i);
    mpz_fdiv_q_ui(tmp, tmp, nthreads);
    mpz_add(tmp, qidx0, tmp);
    {
        std::ostringstream os;
        os << tmp;
        param_list_add_key(pl2, "qidx0", os.str().c_str(), PARAMETER_FROM_CMDLINE);
    }

    mpz_sub(tmp, qidx1, qidx0);
    mpz_mul_ui(tmp, tmp, i + 1);
    mpz_fdiv_q_ui(tmp, tmp, nthreads);
    mpz_add(tmp, qidx0, tmp);
    {
        std::ostringstream os;
        os << tmp;
        param_list_add_key(pl2, "qidx1", os.str().c_str(), PARAMETER_FROM_CMDLINE);
    }
    param_list_remove_key(pl2, "nq");

    return todo_list_base::create<siqs_todo_list>(cpoly, pl2);
}

std::string siqs_todo_list::string_repr() const
{
    return fmt::format("{}[{},{}[ with {} factor(s) in [{}, {}[",
                       random_sampling ? "random sampling from " : "",
                       qidx0, qidx1, qfac_nfac, qfac_min, qfac_max);
}
