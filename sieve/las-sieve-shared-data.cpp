#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include <mutex>
#include <memory>
#include <utility>

#include "cado_poly.h"
#include "ecm/facul_strategies.hpp"
#include "gmp_aux.h"
#include "las-cofactor.hpp"
#include "las-sieve-shared-data.hpp"
#include "las-unsieve.hpp"
#include "macros.h"
#include "memusage.h"
#include "misc.h"
#include "timing.h"
#include "verbose.h"
#include "params.h"
#include "las-siever-config.hpp"
#include "las-side-config.hpp"


void sieve_shared_data::declare_usage(cxx_param_list & pl)
{
    cxx_cado_poly::declare_usage(pl);
    siever_side_config::declare_usage(pl);
    param_list_decl_usage(pl, "fbc",  "factor base cache file");
}

void sieve_shared_data::lookup_parameters(cxx_param_list & pl, int nsides)
{
    /* We don't expect that cxx_cado_poly can be looked up late, so
     * there's no reason to thaw it.
     */
    siever_side_config::lookup_parameters(pl, nsides);
    param_list_lookup_string(pl, "fbc");
}

sieve_shared_data::side_data::side_data(int side,
                                        cxx_cado_poly const& cpoly,
                                        cxx_param_list& pl,
                                        int nthreads)
  : f(cpoly->pols[side])
  , fb(cpoly, side, pl, param_list_lookup_string(pl, "fbc"), nthreads)
{}


/* FIXME: get_trialdiv_data is currently in las-trialdiv.cpp ; it belongs
 * here (that is, the bulk should stay there, but the interaction with
 * the cache mechanism should be be here instead).
 */

/* The fb_factorbase ctor parses the lim[01] and powlim[01] directly from
 * the command line. These get interpreted as the "max" bounds, and we
 * compute the factor base up to that limit.
 */
sieve_shared_data::sieve_shared_data( /*{{{*/
        cxx_cado_poly const & cpoly,
        cxx_param_list & pl)
    : cpoly(cpoly)
    , sides{(size_t) cpoly->nb_polys}
    , cofactfilename { param_list_lookup_string (pl, "file-cofact") }
{
}
void sieve_shared_data::load_factor_base(cxx_param_list & pl, int nthreads) /*{{{*/
{
    for (int i = 0; i < cpoly->nb_polys; i ++)
        sides[i] = side_data {i, cpoly, pl, nthreads};
}
/*}}}*/
/*}}}*/
unsieve_data const * sieve_shared_data::get_unsieve_data(siever_config const & conf) /* {{{ */
{
    std::pair<int, int> const p(conf.logI, conf.logA);
    std::lock_guard<std::mutex> const dummy(us_cache.mutex());
    auto it = us_cache.find(p);
    if (it != us_cache.end()) {
        return &it->second;
    }
    auto itb = us_cache.insert(std::make_pair(p, unsieve_data(p)));
    ASSERT(itb.second);
    return &(*itb.first).second;
}/*}}}*/
j_divisibility_helper const * sieve_shared_data::get_j_divisibility_helper(int J) /* {{{ */
{
    ASSERT_ALWAYS(J);
    /* Round to the next power of two */
    unsigned int const Jround = 1 << nbits(J-1);
    std::lock_guard<std::mutex> const dummy(jdiv_cache.mutex());
    auto it = jdiv_cache.find(Jround);
    if (it != jdiv_cache.end()) {
        return &it->second;
    }
    auto itb = jdiv_cache.insert(std::pair<unsigned int, j_divisibility_helper>(Jround, Jround));
    ASSERT(itb.second);
    return &(*itb.first).second;
}/*}}}*/
facul_strategies const * sieve_shared_data::get_strategies(siever_config const & conf) /* {{{ */
{
    std::lock_guard<std::mutex> const dummy(facul_strategies_cache.mutex());
    auto it = facul_strategies_cache.find(conf);
    if (it != facul_strategies_cache.end()) {
        return it->second.get();
    }

    double const time_strat = seconds();

    std::unique_ptr<FILE, delete_FILE> file;

    if (cofactfilename)
        file.reset(fopen (cofactfilename, "r"));

    auto itb = facul_strategies_cache.insert(std::make_pair(conf,
            std::shared_ptr<facul_strategies>(
                facul_make_strategies (conf, file.get(), 0))));

    ASSERT_ALWAYS(itb.second);
    verbose_fmt_print(0, 1, "# Building/reading strategies took {:1.1f}s\n",
            seconds() - time_strat);

    if (!(*itb.first).second.get()) {
        fprintf (stderr, "impossible to read %s\n", cofactfilename);
        abort ();
    }

    return (*itb.first).second.get();
}/*}}}*/

sieve_shared_data::~sieve_shared_data()
{
    verbose_fmt_print(0, 2,
            "# Getting rid of sieve_shared_data structure [rss={}]\n",
            size_disp_fine(1024UL * Memusage2(), 10000.0));
}
