#include "cado.h" // IWYU pragma: keep
#include <cstdio>            // for fclose, fopen, fprintf, NULL, FILE, stderr
#include <cstdlib>           // for abort
#include <mutex>             // for lock_guard, mutex
#include "ecm/facul.hpp"         // for facul_clear_strategies, facul_strategies_t
#include "gmp_aux.h"    // nbits
#include "las-cofactor.hpp"  // for facul_make_strategies
#include "las-sieve-shared-data.hpp"
#include "las-unsieve.hpp"   // for j_divisibility_helper, unsieve_data
#include "macros.h"          // for ASSERT, ASSERT_ALWAYS
#include "memusage.h"   // Memusage2
#include "misc.h"          // size_disp_fine
#include "timing.h"             // for seconds
#include "verbose.h"             // verbose_output_print
#include "params.h"
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
{
    cofactfilename = param_list_lookup_string (pl, "file-cofact");
}
void sieve_shared_data::load_factor_base(cxx_param_list & pl, int nthreads) /*{{{*/
{
    sides[0] = side_data {0, cpoly, pl, nthreads};
    sides[1] = side_data {1, cpoly, pl, nthreads};
}
/*}}}*/
/*}}}*/
unsieve_data const * sieve_shared_data::get_unsieve_data(siever_config const & conf) /* {{{ */
{
    std::pair<int, int> p(conf.logI, conf.logA);
    std::lock_guard<std::mutex> dummy(us_cache.mutex());
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
    unsigned int Jround = 1 << nbits(J-1);
    std::lock_guard<std::mutex> dummy(jdiv_cache.mutex());
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
    std::lock_guard<std::mutex> dummy(facul_strategies_cache.mutex());
    auto it = facul_strategies_cache.find(conf);
    if (it != facul_strategies_cache.end()) {
        return it->second.get();
    }

    double time_strat = seconds();

    FILE* file = NULL;
    if (cofactfilename != NULL) /* a file was given */
        file = fopen (cofactfilename, "r");

    auto itb = facul_strategies_cache.insert(std::make_pair(conf,
            std::shared_ptr<facul_strategies>(
                facul_make_strategies (conf, file, 0))));

    if (file)
        fclose (file);

    ASSERT_ALWAYS(itb.second);
    verbose_output_print(0, 1, "# Building/reading strategies took %1.1fs\n",
            seconds() - time_strat);

    if (!(*itb.first).second.get()) {
        fprintf (stderr, "impossible to read %s\n", cofactfilename);
        abort ();
    }

    return (*itb.first).second.get();
}/*}}}*/

sieve_shared_data::~sieve_shared_data()
{
    char buf1[16];
    verbose_output_print(0, 2, "# Getting rid of sieve_shared_data structure [rss=%s]\n",
            size_disp_fine(1024UL * Memusage2(), buf1, 10000.0));
}
