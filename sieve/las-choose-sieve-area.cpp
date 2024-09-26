#include "cado.h" // IWYU pragma: keep
#include <cinttypes>               // for PRId64
#include <cstdio>              // IWYU pragma: keep
#include <cstdarg>             // IWYU pragma: keep
#include <gmp.h>                   // for mpz_srcptr, gmp_vfprintf
#include "cado_poly.h"   // cxx_cado_poly
#include "cxx_mpz.hpp"   // for cxx_mpz
#include "fb-types.h"              // for sublat_t
#include "las-auxiliary-data.hpp"  // for nfs_aux, report_and_timer
#include "las-choose-sieve-area.hpp"
#include "las-info.hpp"            // for las_info, HILIGHT_END, HILIGHT_START
#include "las-multiobj-globals.hpp"
#include "las-norms.hpp"           // for sieve_range_adjust
#include "las-qlattice.hpp"        // for qlattice_basis, qlattice_basis::to...
#include "las-siever-config.hpp"   // for siever_config, siever_config_pool
#include "las-todo-entry.hpp"      // for las_todo_entry
#include "macros.h"                // for MAYBE_UNUSED
#include "tdict.hpp"               // for timetree_t
#include "verbose.h"    // verbose_output_vfprint

int never_discard = 0;      /* only enabled for las_descent */

static bool choose_sieve_area(las_info const & las,
        timetree_t * ptimer MAYBE_UNUSED,
        las_todo_entry const & doing,
        siever_config & conf,
        qlattice_basis & Q,
        uint32_t & J)
{
    sieve_range_adjust Adj;

    int adjust_strategy = las.adjust_strategy;
    {

    /* Our business: find an appropriate siever_config, that is
     * appropriate for this special-q. Different special-q's may lead to
     * different siever_config's, it is allowed.
     *
     * This process completes when we also set conf.logI.
     */
    conf = las.config_pool.get_config_for_q(doing);

    /* The whole business about sieve_range_adjust is to compute the
     * optimal logI and J for this special-q. We have several
     * strategies for that.
     */
    try {
        Adj = sieve_range_adjust(doing, las.cpoly, conf);
    } catch (qlattice_basis::too_skewed const & x) {
        verbose_output_vfprint(0, 1, gmp_vfprintf,
                "# "
                HILIGHT_START
                "Discarding side-%d q=%Zd; rho=%Zd (q-lattice basis does not fit)\n"
                HILIGHT_END,
                doing.side,
                (mpz_srcptr) doing.p,
                (mpz_srcptr) doing.r);
        return false;
    }

    if (!support_large_q && !Adj.Q.fits_31bits()) { // for fb_root_in_qlattice_31bits
        verbose_output_print(2, 1,
                "# Warning, special-q basis is too skewed,"
                " skipping this special-q."
                " Define SUPPORT_LARGE_Q to proceed anyway.\n");
        return false;
    }

    }

    /* Must be done before any sort of lognorm computation, even before the
     * adjustment. (only sublat_bound matters, here -- which actual
     * sublattice we will sieve is irrelevant at this point, as what we
     * do will be shared for all sublattices).
     */
    Adj.Q.sublat.m = conf.sublat_bound;

    {

    /* Try strategies for adopting the sieving range */

    int should_discard = !Adj.sieve_info_adjust_IJ();

    if (should_discard) {
        if (never_discard) {
            Adj.set_minimum_J_anyway();
        } else {
            verbose_output_vfprint(0, 1, gmp_vfprintf,
                    "# "
                    HILIGHT_START
                    "Discarding side-%d q=%Zd; rho=%Zd;"
                    HILIGHT_END
                    " a0=%" PRId64
                    "; b0=%" PRId64
                    "; a1=%" PRId64
                    "; b1=%" PRId64
                    "; raw_J=%u;\n", 
                    doing.side,
                    (mpz_srcptr) doing.p,
                    (mpz_srcptr) doing.r,
                    Adj.Q.a0, Adj.Q.b0, Adj.Q.a1, Adj.Q.b1, Adj.J);
            return false;
        }
    }

    /* With adjust_strategy == 2, we want to display the other
     * values, too. Also, strategy 0 wants strategy 1 to run first.
     */
    if (adjust_strategy != 1)
        Adj.sieve_info_update_norm_data_Jmax();

    if (adjust_strategy >= 2)
        Adj.adjust_with_estimated_yield();

    if (adjust_strategy >= 3) {
        /* Let's change that again. We tell the code to keep logI as
         * it is currently. */
        Adj.sieve_info_update_norm_data_Jmax(true);
    }

    /* check whether J is too small after the adjustments */
    if (Adj.J < Adj.get_minimum_J())
    {
        if (never_discard) {
            Adj.set_minimum_J_anyway();
        } else {
            verbose_output_vfprint(0, 1, gmp_vfprintf,
                    "# "
                    HILIGHT_START
                    "Discarding side-%d q=%Zd; rho=%Zd;"
                    HILIGHT_END,
                    doing.side,
                    (mpz_srcptr) doing.p,
                    (mpz_srcptr) doing.r);
            verbose_output_print(0, 1,
                    " a0=%" PRId64
                    "; b0=%" PRId64
                    "; a1=%" PRId64
                    "; b1=%" PRId64
                    "; raw_J=%u;\n",
                    Adj.Q.a0, Adj.Q.b0, Adj.Q.a1, Adj.Q.b1, Adj.J);
            return false;
        }
    }

    /* last check: we don't want the homography induced by the qlattice
     * to trip over a potential small rational root of the polynomial.
     * Of course it is highly unexpected. It can only
     * happen for a linear polynomial, of course, and when it happens, it
     * means that the bivariate fij would only on one of its variables,
     * so we're not going to do anything very interesting.
     *
     */
    {
        int64_t H[4] = { Adj.Q.a0, Adj.Q.b0, Adj.Q.a1, Adj.Q.b1 };
        cxx_mpz_poly fij;
    
        mpz_poly_homography (fij, las.cpoly->pols[doing.side], H);

        if (fij->deg < las.cpoly->pols[doing.side]->deg) {
            verbose_output_vfprint(0, 1, gmp_vfprintf,
                    "# "
                    HILIGHT_START
                    "Discarding side-%d q=%Zd; rho=%Zd;"
                    HILIGHT_END,
                    doing.side,
                    (mpz_srcptr) doing.p,
                    (mpz_srcptr) doing.r);
            verbose_output_print(0, 1,
                    " a0=%" PRId64
                    "; b0=%" PRId64
                    "; a1=%" PRId64
                    "; b1=%" PRId64
                    "; raw_J=%u; // explanation: tripped over rational root.\n",
                    Adj.Q.a0, Adj.Q.b0, Adj.Q.a1, Adj.Q.b1, Adj.J);
            return false;
        }
    }


    /* At this point we've decided on a new configuration for the
     * siever.
     */
    conf.logI = Adj.logI;
    Q = Adj.Q;
    J = Adj.J;

    }

    return true;
}

bool choose_sieve_area(las_info const & las,
        std::shared_ptr<nfs_aux> aux_p,
        las_todo_entry const & doing,
        siever_config & conf,
        qlattice_basis & Q,
        uint32_t & J)
{
    timetree_t & timer(aux_p->rt.timer);
    return choose_sieve_area(las, &timer, doing, conf, Q, J);
}

bool choose_sieve_area(las_info const & las,
        las_todo_entry const & doing,
        siever_config & conf,
        qlattice_basis & Q,
        uint32_t & J)
{
    return choose_sieve_area(las, nullptr, doing, conf, Q, J);
}

