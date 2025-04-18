#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include "fm.hpp"
#include "macros.h"
#include "strategy.hpp"
#include "tab_fm.hpp"
#include "utils_cxx.hpp"

strategy_t * strategy_create()
{
    strategy_t * t = (strategy_t *)malloc(sizeof(strategy_t));
    ASSERT_ALWAYS(t != nullptr);
    t->tab_fm = tabular_fm_create();
    t->proba = 0;
    t->time = 0;
    t->side = nullptr;
    t->len_side = 0;
    return t;
}

void strategy_free(strategy_t * t)
{
    if (t != nullptr) {
        tabular_fm_free(t->tab_fm);
        if (t->side != nullptr)
            free(t->side);
        free(t);
    }
}

tabular_fm_t * strategy_get_tab_fm(strategy_t const * t)
{
    return t->tab_fm;
}

double strategy_get_proba(strategy_t const * t)
{
    return t->proba;
}

double strategy_get_time(strategy_t const * t)
{
    return t->time;
}

void strategy_set_proba(strategy_t * t, double proba)
{
    t->proba = proba;
}

void strategy_set_time(strategy_t * t, double time)
{
    t->time = time;
}

void strategy_add_fm(strategy_t * t, fm_t * elem)
{
    tabular_fm_add_fm(t->tab_fm, elem);
}

void strategy_add_fm_side(strategy_t * t, fm_t * elem, int side)
{
    tabular_fm_add_fm(t->tab_fm, elem);
    if (t->side == nullptr) {
        t->len_side = t->tab_fm->size;
        t->side = (int *)calloc(t->len_side, sizeof(int));
        ASSERT(t->side != nullptr);
    } else {
        t->len_side++;
        checked_realloc(t->side, t->len_side);
    }
    t->side[t->len_side - 1] = side;
}

strategy_t * strategy_copy(strategy_t * t)
{
    strategy_t * elem = strategy_create();
    tabular_fm_concat(elem->tab_fm, t->tab_fm);
    elem->proba = t->proba;
    elem->time = t->time;
    // side
    if (t->side != nullptr) {
        elem->len_side = t->len_side;
        elem->side = (int *)malloc(sizeof(int) * (elem->len_side));
        for (unsigned int i = 0; i < elem->len_side; i++)
            elem->side[i] = t->side[i];
    }
    return elem;
}

int strategy_fprint(FILE * output_file, strategy_t const * t)
{
    if (output_file == nullptr)
        return -1;

    tabular_fm_t const * tmp = t->tab_fm;
    // test if the varaible side is used!
    int is_alloced_side = false;
    if (t->side != nullptr && t->len_side == tmp->size)
        is_alloced_side = true;

    for (unsigned int i = 0; i < tmp->size; i++) {
        fm_t * fm = tmp->tab[i];
        for (unsigned int j = 0; j < fm->len_method; j++)
            fprintf(output_file, "%lu ", fm->method[j]);
        if (is_alloced_side)
            fprintf(output_file, "%d", t->side[i]);
        else
            fprintf(output_file, "%d", 0); // default side!

        fprintf(output_file, "\n");
    }
    fprintf(output_file, "Probability: %1.10lf\n", t->proba);
    fprintf(output_file, "Time: %lf\n", t->time);
    return 0;
}

int strategy_print(strategy_t const * t)
{
    return strategy_fprint(stdout, t);
}
