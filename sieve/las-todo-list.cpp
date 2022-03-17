#include "cado.h" // IWYU pragma: keep

#include <cctype>          // for isspace, isdigit
#include <cerrno>          // for errno
#include <climits>         // for UINT_MAX
#include <cstdlib>         // for exit, EXIT_FAILURE, strtoul
#include <cstring>         // for strerror
#include <cstdarg>         // IWYU pragma: keep
#include <algorithm>       // for find, min
#include <memory>          // for allocator_traits<>::value_type
#include <vector>          // for vector, vector<>::iterator
#include <sstream>
#include <thread>                         // for thread
#include <gmp.h>           // for mpz_set, mpz_cmp_ui, mpz_t, mpz_cmp, gmp_r...
#include "gmp_aux.h"
#include "las-todo-list.hpp"
#include "las-galois.hpp"  // for skip_galois_roots
#include "macros.h"        // for ASSERT_ALWAYS
#include "mpz_poly.h"
#include "rootfinder.h" // mpz_poly_roots_ulong
#include "verbose.h"    // verbose_output_print
#include "params.h"

/* Put in r the smallest legitimate special-q value that it at least
 * s + diff (note that if s+diff is already legitimate, then r = s+diff
 * will result.
 * In case of composite sq, also store the factorization of r in fac_r
 */

static void next_legitimate_specialq(cxx_mpz & r, std::vector<uint64_t> & fac_r, cxx_mpz const & s, const unsigned long diff, las_todo_list const & L)
{
    if (L.allow_composite_q) {
        unsigned long tfac[64];
        int nf = next_mpz_with_factor_constraints(r, tfac,
                s, diff, L.qfac_min, L.qfac_max);
        fac_r.assign(tfac, tfac + nf);
    } else {
        mpz_add_ui(r, s, diff);
        /* mpz_nextprime() returns a prime *greater than* its input argument,
           which we don't always want, so we subtract 1 first. */
        mpz_sub_ui(r, r, 1);
        mpz_nextprime(r, r);
    }
}

static void next_legitimate_specialq(cxx_mpz & r, cxx_mpz const & s, const unsigned long diff, las_todo_list const & L)
{
    std::vector<uint64_t> t;
    next_legitimate_specialq(r, t, s, diff, L);
}

void las_todo_list::configure_switches(cxx_param_list & pl)
{
    param_list_configure_switch(pl, "-allow-compsq", NULL);
    param_list_configure_switch(pl, "-print-todo-list", NULL);
}

void las_todo_list::declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "sqside", "put special-q on this side");
    param_list_decl_usage(pl, "q0",   "left bound of special-q range");
    param_list_decl_usage(pl, "q1",   "right bound of special-q range");
    param_list_decl_usage(pl, "rho",  "sieve only root r mod q0");
    param_list_decl_usage(pl, "random-sample", "Sample this number of special-q's at random, within the range [q0,q1]");
    param_list_decl_usage(pl, "nq", "Process this number of special-q's and stop");
    param_list_decl_usage(pl, "todo", "provide file with a list of special-q to sieve instead of qrange");
    param_list_decl_usage(pl, "allow-compsq", "allows composite special-q");
    param_list_decl_usage(pl, "qfac-min", "factors of q must be at least that");
    param_list_decl_usage(pl, "qfac-max", "factors of q must be at most that");
    param_list_decl_usage(pl, "print-todo-list", "only print the special-q's to be sieved");
}

las_todo_list::~las_todo_list()
{
    if (todo_list_fd) {
        fclose(todo_list_fd);
        todo_list_fd = NULL;
    }
}

las_todo_list::las_todo_list(cxx_cado_poly const & cpoly, cxx_param_list & pl)
    : cpoly(cpoly)
{
    nq_pushed = 0;
    nq_max = UINT_MAX;
    random_sampling = 0;
    if (param_list_parse_uint(pl, "random-sample", &nq_max)) {
        random_sampling = 1;
    } else if (param_list_parse_uint(pl, "nq", &nq_max)) {
        if (param_list_lookup_string(pl, "rho")) {
            fprintf(stderr, "Error: argument -nq is incompatible with -rho\n");
            exit(EXIT_FAILURE);
        }
        if (param_list_lookup_string(pl, "q1"))
            verbose_output_print(0, 1, "# Warning: argument -nq takes priority over -q1 ; -q1 ignored\n");
    }

    sqside = 1;
    if (!param_list_parse_int(pl, "sqside", &sqside)) {
        verbose_output_print(0, 1, "# Warning: sqside not given, "
                "assuming side 1 for backward compatibility.\n");
    }

    /* Init and parse info regarding work to be done by the siever */
    /* Actual parsing of the command-line fragments is done within
     * las_todo_feed, but this is an admittedly contrived way to work */
    const char * filename = param_list_lookup_string(pl, "todo");
    if (filename) {
        todo_list_fd = fopen(filename, "r");
        if (todo_list_fd == NULL) {
            fprintf(stderr, "%s: %s\n", filename, strerror(errno));
            /* There's no point in proceeding, since it would really change
             * the behaviour of the program to do so */
            exit(EXIT_FAILURE);
        }
    } else {
        todo_list_fd = NULL;
    }

    /* composite special-q ? Note: this block is present both in
     * las-todo-list.cpp and las-info.cpp */
    if ((allow_composite_q = param_list_parse_switch(pl, "-allow-compsq"))) {
        /* defaults are set in the class description */
        param_list_parse_uint64(pl, "qfac-min", &qfac_min);
        param_list_parse_uint64(pl, "qfac-max", &qfac_max);
    }

    galois = param_list_lookup_string(pl, "galois");

    if (allow_composite_q && galois) {
        fprintf(stderr, "-galois and -allow-compsq are incompatible options at the moment");
        exit(EXIT_FAILURE);
    }

    print_todo_list_flag = param_list_parse_switch(pl, "-print-todo-list");

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

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

    if (mpz_cmp_ui(q1, 0) != 0) {
        next_legitimate_specialq(q0, q0, 0, *this);
    } else {
        /* We don't have -q1. If we have -rho, we sieve only <q0,
         * rho>. */
        cxx_mpz t;
        if (param_list_parse_mpz(pl, "rho", (mpz_ptr) t)) {
            cxx_mpz q0_cmdline = q0;
            std::vector<uint64_t> fac_q;
            next_legitimate_specialq(q0, fac_q, q0, 0, *this);
            if (mpz_cmp(q0, q0_cmdline) != 0) {
                fprintf(stderr, "Error: q0 is not a legitimate special-q\n");
                exit(EXIT_FAILURE);
            }
            std::vector<cxx_mpz> roots = mpz_poly_roots(cpoly->pols[sqside], q0, fac_q, rstate);
            if (std::find(roots.begin(), roots.end(), t) == roots.end()) {
                fprintf(stderr, "Error: rho is not a root modulo q0\n");
                exit(EXIT_FAILURE);
            }
            push_unlocked(q0, t, sqside);
            /* Set empty interval [q0 + 1, q0] as special-q interval */
            mpz_set(q1, q0);
            mpz_add_ui (q0, q0, 1);
        } else {
            /* If we don't have -rho, we sieve only q0, but all roots of it.
               If -q0 does not give a legitimate special-q value, advance to the
               next legitimate one. */
            mpz_set(t, q0);
            next_legitimate_specialq(q0, q0, 0, *this);
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
            std::vector<uint64_t> fac_q;
            next_legitimate_specialq(q, fac_q, q, 0, *this);
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
            gmp_fprintf(stderr, "Error: range [%Zd,%Zd[ contains no prime with roots mod f\n",
                    (mpz_srcptr) q0,
                    (mpz_srcptr) q1_orig);
            exit(EXIT_FAILURE);
        }
    }

    gmp_randclear(rstate);
}

/* {{{ Populating the todo list */
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
    /* If we still have entries in the stack, don't add more now */
    if (!super::empty())
        return true;

    mpz_poly_ptr f = cpoly->pols[sqside];

    std::vector<uint64_t> fac_q;

    if (!random_sampling) {
        /* We're going to process the sq's and put them into the list
           The loop processes all special-q in [q0, q1]. On loop entry,
           the value in q0 is known to be a legitimate special-q. Its
           factorization is lost, so we recompute it. */

        /* handy aliases */
        cxx_mpz & q = q0;
        next_legitimate_specialq(q, fac_q, q, 0, *this);

        struct q_r_pair {
            cxx_mpz q;
            cxx_mpz r;
            q_r_pair(const mpz_t _q, const mpz_t _r) {
                mpz_set(q, _q);
                mpz_set(r, _r);
            }
        };

        std::vector<q_r_pair> my_list;

        int nb_no_roots = 0;
        int nb_rootfinding = 0;
        /* If nq_max is specified, then q1 has no effect, even though it
         * has been set equal to q */
        for ( ; (nq_max < UINT_MAX || mpz_cmp(q, q1) < 0) &&
                nq_pushed + my_list.size() < nq_max ; )
        {
            std::vector<cxx_mpz> roots = mpz_poly_roots(f, q, fac_q, rstate);

            nb_rootfinding++;
            if (roots.empty()) nb_no_roots++;

            if (galois) {
                size_t nroots = skip_galois_roots(roots.size(), q, (mpz_t*)roots.data(), galois);
                roots.erase(roots.begin() + nroots, roots.end());
            }

            for (auto const & r : roots)
                my_list.push_back({q, r});

            next_legitimate_specialq(q, fac_q, q, 1, *this);
        }

        if (nb_no_roots) {
            verbose_output_vfprint(0, 1, gmp_vfprintf, "# polynomial has no roots for %d of the %d primes that were tried\n", nb_no_roots, nb_rootfinding);
        }

        // Truncate to nq_max if necessary and push the sq in reverse
        // order, because they are processed via a stack (required for
        // the descent).
        int push_here = my_list.size();
        if (nq_max < UINT_MAX)
            push_here = std::min(push_here, int(nq_max - nq_pushed));
        for(int i = 0 ; i < push_here ; i++) {
            nq_pushed++;
            int ind = push_here-i-1;
            push_unlocked(my_list[ind].q, my_list[ind].r, sqside);
        }
    } else { /* random sampling case */
        /* we care about being uniform here */
        cxx_mpz q;
        cxx_mpz diff;
        mpz_sub(diff, q1, q0);
        ASSERT_ALWAYS(nq_pushed == 0 || nq_pushed == nq_max);
	unsigned long n = nq_max;
        unsigned int spin = 0;
        for ( ; nq_pushed < n ; ) {
            /* try in [q0 + k * (q1-q0) / n, q0 + (k+1) * (q1-q0) / n[ */
            cxx_mpz q0l, q1l;
	    /* we use k = n-1-nq_pushed instead of k=nq_pushed so that
	       special-q's are sieved in increasing order */
	    unsigned long k = n - 1 - nq_pushed;
            mpz_mul_ui(q0l, diff, k);
            mpz_mul_ui(q1l, diff, k + 1);
            mpz_fdiv_q_ui(q0l, q0l, n);
            mpz_fdiv_q_ui(q1l, q1l, n);
            mpz_add(q0l, q0, q0l);
            mpz_add(q1l, q0, q1l);

            mpz_sub(q, q1l, q0l);
            mpz_urandomm(q, rstate, q);
            mpz_add(q, q, q0l);
            next_legitimate_specialq(q, fac_q, q, 0, *this);
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
                size_t nroots = skip_galois_roots(roots.size(), q, (mpz_t*)roots.data(), galois);
                roots.erase(roots.begin() + nroots, roots.end());
            }
            nq_pushed++;
            push_unlocked(q, roots[gmp_urandomm_ui(rstate, roots.size())], sqside);
        }
    }

    return !super::empty();
}

/* Format of a file with a list of special-q (-todo option):
 *   - Comments are allowed (start line with #)
 *   - Blank lines are ignored
 *   - Each valid line must have the form
 *       s q r
 *     where s is the side (0 or 1) of the special q, and q and r are as usual.
 */
bool las_todo_list::feed_qlist()
{
    if (!super::empty())
        return true;

    char line[1024];
    char * x;
    for( ; ; ) {
        /* We should consider the blocking case as well. (e.g. we're
         * reading from a pipe.)  */
        x = fgets(line, sizeof(line), todo_list_fd);
        /* Tolerate comments and blank lines */
        if (x == NULL) return 0;
        if (*x == '#') continue;
        for( ; *x && isspace(*x) ; x++) ;
        if (!*x) continue;
        break;
    }

    /* We have a new entry to parse */
    cxx_mpz p, r;
    int side = -1;
    int rc;
    switch(*x++) {
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
                   x--;
                   errno = 0;
                   side = strtoul(x, &x, 0);
                   ASSERT_ALWAYS(!errno);
                   ASSERT_ALWAYS(side < cpoly->nb_polys);
                   break;
        default:
                   fprintf(stderr,
                           "parse error in todo file, while reading: %s\n",
                           line);
                   /* We may as well default on the command-line switch */
                   exit(EXIT_FAILURE);
    }

    int nread1 = 0;
    int nread2 = 0;

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);
    mpz_set_ui(r, 0);
    for( ; *x && !isdigit(*x) ; x++) ;
    rc = gmp_sscanf(x, "%Zi%n %Zi%n", (mpz_ptr) p, &nread1, (mpz_ptr) r, &nread2);
    ASSERT_ALWAYS(rc == 1 || rc == 2); /* %n does not count */
    x += (rc==1) ? nread1 : nread2;
    {
        mpz_poly_ptr f = cpoly->pols[side];
        /* specifying the rational root as <0
         * means that it must be recomputed. Putting 0 does not have this
         * effect, since it is a legitimate value after all.
         */
        if (rc < 2 || (f->deg == 1 && rc == 2 && mpz_cmp_ui(r, 0) < 0)) {
            // For rational side, we can compute the root easily.
            ASSERT_ALWAYS(f->deg == 1);
            /* ugly cast, yes */
            int nroots = mpz_poly_roots ((mpz_t*) &r, f, p, rstate);
            ASSERT_ALWAYS(nroots == 1);
        }
    }
    gmp_randclear(rstate);

    for( ; *x ; x++) ASSERT_ALWAYS(isspace(*x));
    push_unlocked(p, r, side);
    return true;
}


bool las_todo_list::feed(gmp_randstate_t rstate)
{
    std::lock_guard<std::mutex> foo(mm);
    if (!super::empty())
        return true;
    if (todo_list_fd)
        return feed_qlist();
    else
        return feed_qrange(rstate);
}

/* This exists because of the race condition between feed() and pop()
 */
las_todo_entry * las_todo_list::feed_and_pop(gmp_randstate_t rstate)
{
    std::lock_guard<std::mutex> foo(mm);
    if (super::empty()) {
        if (todo_list_fd)
            feed_qlist();
        else
            feed_qrange(rstate);
    }
    if (super::empty())
        return nullptr;
    las_todo_entry doing = super::top();
    super::pop();
    history.push_back(doing);
    return &history.back();
}
/* }}} */

void las_todo_list::print_todo_list(cxx_param_list & pl, gmp_randstate_ptr rstate, int nthreads) const
{

    if (random_sampling) {
        /* Then we cannot print the todo list in a multithreaded way
         * without breaking the randomness predictability. It's a bit
         * annoying, yes. On the other hand, it's highly likely in this
         * case that the number of q's to pick is small enough anyway !
         */
        nthreads = 1;
    }

    std::vector<std::vector<las_todo_entry>> lists(nthreads);
    std::vector<std::thread> subjobs;

    auto segment = [&, this](unsigned int i) {
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
        las_todo_list todo2(cpoly, pl2);
        for(;;) {
            las_todo_entry * doing_p = todo2.feed_and_pop(rstate);
            if (!doing_p) break;
            las_todo_entry& doing(*doing_p);
            if (nthreads == 1) {
                verbose_output_vfprint(0, 1, gmp_vfprintf,
                        "%d %Zd %Zd\n",
                        doing.side,
                        (mpz_srcptr) doing.p,
                        (mpz_srcptr) doing.r);
            } else {
                lists[i].push_back(doing);
            }
        }
    };

    if (nthreads > 1) {
        verbose_output_vfprint(0, 1, gmp_vfprintf,
                "# Collecting the todo list in memory from q0=%Zd to q1=%Zd using %d threads\n",
                (mpz_srcptr) q0, (mpz_srcptr) q1, nthreads);
    }

    for(int subjob = 0 ; subjob < nthreads ; ++subjob)
        subjobs.push_back(std::thread(segment, subjob));
    for(auto & t : subjobs) t.join();
    for(auto const & v : lists) {
        for(auto const & doing : v) {
            verbose_output_vfprint(0, 1, gmp_vfprintf,
                    "%d %Zd %Zd\n",
                    doing.side,
                    (mpz_srcptr) doing.p,
                    (mpz_srcptr) doing.r);
        }
    }
}
