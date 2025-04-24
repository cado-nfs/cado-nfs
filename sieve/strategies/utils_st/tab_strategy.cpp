#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include "fm.hpp"
#include "macros.h"
#include "strategy.hpp"
#include "tab_strategy.hpp"
#include "utils_cxx.hpp"
#include "macros.h"

tabular_strategy_t * tabular_strategy_create(void)
{
    tabular_strategy_t * t =
        (tabular_strategy_t *)malloc(sizeof(tabular_strategy_t));
    ASSERT_ALWAYS(t != NULL);

    t->size = 0;
    t->alloc = 2;

    t->tab = (strategy_t **)malloc(t->alloc * sizeof(strategy_t *));
    ASSERT_ALWAYS(t->tab != NULL);

    return t;
}

void tabular_strategy_free(tabular_strategy_t * t)
{
    if (t != NULL) {
        for (unsigned int i = 0; i < t->size; i++)
            strategy_free(t->tab[i]);
        free(t->tab);
        free(t);
    }
}

void tabular_strategy_realloc(tabular_strategy_t * t)
{
    checked_realloc(t->tab, t->alloc * 2);
    t->alloc *= 2;
}

tabular_strategy_t * tabular_strategy_copy(tabular_strategy_t * t)
{
    tabular_strategy_t * res = tabular_strategy_create();
    unsigned int len = t->size;
    for (unsigned int i = 0; i < len; i++) {
        tabular_strategy_add_strategy(res, t->tab[i]);
    }
    return res;
}

unsigned int tabular_strategy_get_size(tabular_strategy_t const * t)
{
    return t->size;
}

void tabular_strategy_add_strategy(tabular_strategy_t * t,
                                   strategy_t * strategy)
{
    if (t->size >= t->alloc)
        tabular_strategy_realloc(t);
    strategy_t * elem = strategy_copy(strategy);
    t->tab[t->size] = elem;
    t->size++;
}

void tabular_strategy_concat(tabular_strategy_t * t1, tabular_strategy_t * t2)
{
    unsigned int len = t2->size;
    for (unsigned int i = 0; i < len; i++)
        tabular_strategy_add_strategy(t1, t2->tab[i]);
}

tabular_strategy_t * tabular_strategy_concat_st(tabular_strategy_t * t1,
                                                tabular_strategy_t * t2)
{
    tabular_strategy_t * t = tabular_strategy_create();
    tabular_strategy_concat(t, t1);
    tabular_strategy_concat(t, t2);
    return t;
}

/************************************************************************/
/*           PRINT AND SCAN OUR FILES OF FACTORING METHODS              */
/************************************************************************/

int tabular_strategy_fprint(FILE * output_file, tabular_strategy_t * t)
{
    for (unsigned int i = 0; i < t->size; i++)
        if (strategy_fprint(output_file, t->tab[i]) == -1)
            return -1;
    return 0;
}

int tabular_strategy_print(tabular_strategy_t * t)
{
    return tabular_strategy_fprint(stdout, t);
}

static int is_number(char const c)
{
    return (c >= 48 && c <= 57);
}

static void next_number(FILE * file, int * current_char)
{
    // end the current number
    while (is_number(*current_char))
        *current_char = fgetc(file);

    // find the next number
    while (*current_char != EOF && !is_number(*current_char)) {
        *current_char = fgetc(file);
    }
    // WTF? ungetc maybe?
    int rc = fseek(file, -1, SEEK_CUR);
    DIE_ERRNO_DIAG(rc < 0, "rewind(%s)", "strategy file");
}

tabular_strategy_t * tabular_strategy_fscan(FILE * file)
{
    if (file == NULL)
        return NULL;

    // allocate tabular_strategy
    tabular_strategy_t * tab = tabular_strategy_create();

    // collect data
    int current_char = fgetc(file);
    int side;
    int rc;

    while (current_char != EOF) {
        //{{collect strat
        strategy_t * strat = strategy_create();
        //{{{collect fm
        while (is_number(current_char)) {
            fm_t * elem = fm_create();
            fseek(file, -1, SEEK_CUR);
            rc = fscanf(file, "%lu", &elem->method[0]);
            ASSERT_ALWAYS(rc == 1);
            next_number(file, &current_char);

            rc = fscanf(file, "%lu", &elem->method[1]);
            ASSERT_ALWAYS(rc == 1);
            next_number(file, &current_char);

            rc = fscanf(file, "%lu", &elem->method[2]);
            ASSERT_ALWAYS(rc == 1);
            next_number(file, &current_char);

            rc = fscanf(file, "%lu", &elem->method[3]);
            ASSERT_ALWAYS(rc == 1);
            next_number(file, &current_char);

            rc = fscanf(file, "%d", &side);
            ASSERT_ALWAYS(rc == 1);
            // go to end of line: 10 = '\t'
            while (current_char != 10)
                current_char = fgetc(file);
            current_char = fgetc(file);

            strategy_add_fm_side(strat, elem, side);
            fm_free(elem);
        }
        next_number(file, &current_char);
        rc = fscanf(file, "%lf", &strat->proba);
        ASSERT_ALWAYS(rc == 1);

        next_number(file, &current_char);
        rc = fscanf(file, "%lf", &strat->time);
        ASSERT_ALWAYS(rc == 1);

        next_number(file, &current_char);
        //}}
        // add to the tabular
        tabular_strategy_add_strategy(tab, strat);
        strategy_free(strat);
    }
    return tab;
}
