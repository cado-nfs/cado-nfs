#ifndef STRATEGY_H
#define STRATEGY_H

#include "tab_fm.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct strategy {
    tabular_fm_t *tab_fm;
    double proba;
    double time;
    int *side;
    int len_side;
    //In practice, we use side only one time. So the real and physical
    //size are the same, and it's not necessary to allocate it for now.

} strategy_t;

strategy_t *strategy_create();

void strategy_free(strategy_t * t);

tabular_fm_t *strategy_get_tab_fm(strategy_t * t);

double strategy_get_proba(strategy_t * t);

double strategy_get_time(strategy_t * t);

void strategy_set_proba(strategy_t * t, double proba);

void strategy_set_time(strategy_t * t, double time);

void strategy_add_fm(strategy_t * t, fm_t * elem);

void strategy_add_fm_side(strategy_t * t, fm_t * elem, int side);

strategy_t *strategy_copy(strategy_t * t);

int strategy_fprint(FILE * file, strategy_t * t);

int strategy_print(strategy_t * t);

#ifdef __cplusplus
}
#endif

#endif				/* STRATEGY_H */
