#include "cado.h" // IWYU pragma: keep

#include <cstdint>

#include <memory>

#include <gmp.h>

#include "cado_poly.h"
#include "fb-types.hpp"
#include "las-auxiliary-data.hpp"
#include "las-choose-sieve-area.hpp"
#include "las-info.hpp"
#include "las-multiobj-globals.hpp"
#include "las-norms.hpp"
#include "las-qlattice.hpp"
#include "las-siever-config.hpp"
#include "las-special-q-task.hpp"
#include "macros.h"
#include "mpz_poly.h"
#include "tdict.hpp"
#include "verbose.h"

int never_discard = 0;      /* only enabled for las_descent */

static bool choose_sieve_area(las_info const & las,
        timetree_t * ptimer MAYBE_UNUSED,
        special_q_task const & doing,
        siever_config & conf,
        qlattice_basis & Q,
        uint32_t & J)
{
    /* Q is output-only, it's overwritten by the result of the
     * sieve_range_adjust ctor */

    std::unique_ptr<sieve_range_adjust> Adj;

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
        Adj = std::make_unique<sieve_range_adjust>(doing, las.cpoly, conf);
    } catch (qlattice_basis::too_skewed const & x) {
        verbose_fmt_print(0, 1,
                "# Discarding {} (q-lattice basis does not fit)\n",
                doing.sq());
        return false;
    }

    if (!support_large_q && !Adj->Q.fits_31bits()) { // for fb_root_in_qlattice_31bits
        verbose_fmt_print(2, 1,
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
    Adj->Q.sublat.m = conf.sublat_bound;

    {

    /* Try strategies for adopting the sieving range */

    int const should_discard = !Adj->sieve_info_adjust_IJ();

    if (should_discard) {
        if (never_discard) {
            Adj->set_minimum_J_anyway();
        } else {
            verbose_fmt_print(0, 1,
                    "# Discarding {}; raw_J=%u;\n", Adj->Q, Adj->J);
            return false;
        }
    }

    int const adjust_strategy = conf.adjust_strategy;

    /* With adjust_strategy == 2, we want to display the other
     * values, too. Also, strategy 0 wants strategy 1 to run first.
     */
    if (adjust_strategy != 1)
        Adj->sieve_info_update_norm_data_Jmax();

    if (adjust_strategy >= 2)
        Adj->adjust_with_estimated_yield();

    if (adjust_strategy >= 3) {
        /* Let's change that again. We tell the code to keep logI as
         * it is currently. */
        Adj->sieve_info_update_norm_data_Jmax(true);
    }

    /* check whether J is too small after the adjustments */
    if (Adj->J < Adj->get_minimum_J())
    {
        if (never_discard) {
            Adj->set_minimum_J_anyway();
        } else {
            verbose_fmt_print(0, 1,
                    "# Discarding {}; raw_J={};\n", Adj->Q, Adj->J);
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
        cxx_mpz_poly fij = las.cpoly[doing.side].homography(
                { Adj->Q.a0, Adj->Q.b0, Adj->Q.a1, Adj->Q.b1 });

        if (fij->deg < las.cpoly->pols[doing.side]->deg) {
            verbose_fmt_print(0, 1,
                    "# Discarding {}; raw_J={};"
                    " // explanation: tripped over rational root."
                    "\n", Adj->Q, Adj->J);
            return false;
        }
    }

    /* At this point we've decided on a new configuration for the
     * siever.
     */
    conf.logI = Adj->logI;
    Q = Adj->Q;
    J = Adj->J;

    }

    return true;
}

static bool choose_sieve_area(las_info const & las,
        timetree_t * ptimer MAYBE_UNUSED,
        special_q_task const & doing,
        siever_config & conf,
        siqs_special_q_data & Q,
        uint32_t & J)
{
    conf = las.config_pool.get_config_for_q(doing);
    auto n = doing.sq().prime_factors.size();
    conf.logI = conf.logA-(n-1);

    /* -n/2*q - I/2*q <= a = rj + i*q < n/2*q + I/2*q */
    cxx_mpz a_ub = doing.p << (conf.logI-1);
    mpz_addmul_ui(a_ub, doing.p, n/2);
    if (!support_large_q && a_ub.bits() > 63) {
        verbose_output_print(2, 1,
                "# Warning, special-q and I are too large to fit in 63 bits,"
                " skipping this special-q."
                " Define SUPPORT_LARGE_Q to proceed anyway.\n");
        return false;
    }

    Q = siqs_special_q_data(doing.sq(), las.cpoly);
    J = 1U << (n-1);
    return true;
}

template<>
bool choose_sieve_area(las_info const & las,
        std::shared_ptr<nfs_aux> const & aux_p,
        special_q_task const & doing,
        siever_config & conf,
        qlattice_basis & Q,
        uint32_t & J)
{
    timetree_t & timer(aux_p->rt.timer);
    return choose_sieve_area(las, &timer, doing, conf, Q, J);
}

template<>
bool choose_sieve_area(las_info const & las,
        std::shared_ptr<nfs_aux> const & aux_p,
        special_q_task const & doing,
        siever_config & conf,
        siqs_special_q_data & Q,
        uint32_t & J)
{
    timetree_t & timer(aux_p->rt.timer);
    return choose_sieve_area(las, &timer, doing, conf, Q, J);
}

template<>
bool choose_sieve_area(las_info const & las,
        special_q_task const & doing,
        siever_config & conf,
        qlattice_basis & Q,
        uint32_t & J)
{
    return choose_sieve_area(las, nullptr, doing, conf, Q, J);
}

template<>
bool choose_sieve_area(las_info const & las,
        special_q_task const & doing,
        siever_config & conf,
        siqs_special_q_data & Q,
        uint32_t & J)
{
    return choose_sieve_area(las, nullptr, doing, conf, Q, J);
}
