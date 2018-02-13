#include "cado.h"
#include <stdlib.h>
#include <string.h>
#include "las-report-stats.hpp"

void las_report_init(las_report_ptr p)
{
    memset(p, 0, sizeof(las_report));
    p->report_sizes = (long unsigned int (*)[256]) malloc(sizeof(unsigned long) << 16);
    p->survivor_sizes = (long unsigned int (*)[256]) malloc(sizeof(unsigned long) << 16);
}

void las_report_clear(las_report_ptr p)
{
    free(p->report_sizes);
    free(p->survivor_sizes);
    memset(p, 0, sizeof(las_report));
}

void las_report_copy(las_report_ptr p, las_report_ptr q)
{
    unsigned long (*ss)[256] = p->survivor_sizes;
    unsigned long (*rs)[256] = p->report_sizes;
    memcpy(p, q, sizeof(las_report));
    p->survivor_sizes = ss;
    p->report_sizes = rs;
    memcpy(p->survivor_sizes, q->survivor_sizes, sizeof(unsigned long) << 16);
    memcpy(p->report_sizes, q->report_sizes, sizeof(unsigned long) << 16);
}

void las_report_accumulate_and_clear(las_report_ptr p, las_report_ptr q)
{
    unsigned long (*ss)[256] = q->survivor_sizes;
    unsigned long (*rs)[256] = q->report_sizes;
    p->reports += q->reports;
    unsigned long * ps = (unsigned long*) & p->survivors;
    unsigned long * qs = (unsigned long*) & q->survivors;
    for(size_t i = 0 ; i < sizeof(p->survivors) / sizeof(unsigned long) ; i++) {
        ps[i] += qs[i];
    }
    p->ttbuckets_fill  += q->ttbuckets_fill;
    p->ttbuckets_apply += q->ttbuckets_apply;
    p->ttf     += q->ttf;
    p->ttcof     += q->ttcof;
    for(int side = 0 ; side < 2 ; side++) {
        p->tn[side]  += q->tn[side];
    }
    for (int i1 = 0; i1 < 256; i1++) {
        for (int i2 = 0; i2 < 256; i2++) {
            p->survivor_sizes[i1][i2] += q->survivor_sizes[i1][i2];
            p->report_sizes[i1][i2] += q->report_sizes[i1][i2];
        }
    }
    memset(q, 0, sizeof(las_report));
    q->survivor_sizes = ss;
    q->report_sizes = rs;
}
