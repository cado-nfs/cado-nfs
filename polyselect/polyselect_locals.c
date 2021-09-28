#include "cado.h"
#include "polyselect_locals.h"

void polyselect_thread_locals_init(polyselect_thread_locals_ptr loc, polyselect_main_data_ptr main, int idx)
{
    mpz_init(loc->ad);
    loc->main = main;
    mpz_set_ui(loc->ad, main->incr);
    mpz_mul_si(loc->ad, loc->ad, idx);
    mpz_add(loc->ad, loc->ad, main->admin);

    polyselect_poly_header_init(loc->header, main->N, main->d, loc->ad);
    polyselect_proots_init(loc->R, main->lenPrimes);

    polyselect_stats_init(loc->stats);
    loc->rstate = loc->stats->rstate;
}

void polyselect_thread_locals_clear(polyselect_thread_locals_ptr loc)
{
    polyselect_stats_clear(loc->stats);
    polyselect_poly_header_clear(loc->header);
    polyselect_proots_clear(loc->R);
    mpz_clear(loc->ad);

}

