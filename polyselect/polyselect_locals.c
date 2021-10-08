#include "cado.h"
#include "polyselect_locals.h"
#include "polyselect_match.h"
#include "dllist.h"

void polyselect_thread_locals_init(polyselect_thread_locals_ptr loc, polyselect_main_data_srcptr main)
{
    mpz_init(loc->ad);
    loc->main = main;

    polyselect_poly_header_init(loc->header, main->N, main->d);

    polyselect_proots_init(loc->R, main->d, main->lenPrimes);

    polyselect_stats_init(loc->stats, main->keep);
    loc->rstate = loc->stats->rstate;

    dllist_init_head(&loc->async_jobs);
    dllist_init_head(&loc->empty_job_slots);
}

void polyselect_thread_locals_set_idx(polyselect_thread_locals_ptr loc, int idx)
{
    mpz_set_ui(loc->ad, loc->main->incr);
    mpz_mul_si(loc->ad, loc->ad, idx);
    mpz_add(loc->ad, loc->ad, loc->main->admin);

    polyselect_poly_header_set_ad(loc->header, loc->ad);
}

void polyselect_thread_locals_provision_job_slots(polyselect_thread_locals_ptr loc, size_t n)
{
    for( ; n-- ; ) {
        polyselect_match_info_ptr job = malloc(sizeof(polyselect_match_info_t));
        polyselect_match_info_init_trivial(job);
        dllist_push_back(&loc->empty_job_slots, &job->queue);
    }
}


void polyselect_thread_locals_unprovision_job_slots(polyselect_thread_locals_ptr loc)
{
    dllist_for_each_safe(&loc->empty_job_slots, ptr) {
        dllist_pop(ptr);
        polyselect_match_info_ptr job = dllist_entry(ptr, struct polyselect_match_info_s, queue);
        polyselect_match_info_clear(job);
        free(job);
    }
}

void polyselect_thread_locals_clear(polyselect_thread_locals_ptr loc)
{
    ASSERT_ALWAYS(dllist_is_empty(&loc->async_jobs));

    polyselect_thread_locals_unprovision_job_slots(loc);

    polyselect_stats_clear(loc->stats);
    polyselect_poly_header_clear(loc->header);
    polyselect_proots_clear(loc->R);
    mpz_clear(loc->ad);
}

