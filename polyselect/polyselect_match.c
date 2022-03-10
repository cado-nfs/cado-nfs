#include "cado.h"
#include "polyselect_match.h"
#include "polyselect_thread.h"

void polyselect_match_info_clear(polyselect_match_info_ptr task)
{
    /* we have to support trivial instances. */
    if (mpz_cmp_ui(task->header->N, 0) == 0) return;
    polyselect_poly_header_clear(task->header);
    mpz_clear(task->rq);
}

void polyselect_match_info_init_trivial(polyselect_match_info_ptr task)
{
    memset(task, 0, sizeof(polyselect_match_info_t));
    mpz_init(task->rq);
    polyselect_poly_header_init(task->header);
}

void polyselect_match_info_set(polyselect_match_info_ptr job, unsigned long p1, unsigned long p2, const int64_t i, uint64_t q, mpz_srcptr rq, polyselect_thread_srcptr loc)
{
    job->p1 = p1;
    job->p2 = p2;
    job->i = i;
    job->q = q;
    if (rq) {
        mpz_set(job->rq, rq);
    } else {
        mpz_set_ui(job->rq, 0);
        job->q = 1;
    }
    polyselect_poly_header_set(job->header, loc->team->header);
    dllist_init_node(&job->queue);
}

void polyselect_match_info_init(polyselect_match_info_ptr job, unsigned long p1, unsigned long p2, const int64_t i, uint64_t q, mpz_srcptr rq, polyselect_thread_srcptr loc)
{
    polyselect_match_info_init_trivial(job);
    polyselect_match_info_set(job, p1, p2, i, q, rq, loc);
}

