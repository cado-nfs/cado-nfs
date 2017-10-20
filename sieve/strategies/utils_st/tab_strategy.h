#ifndef TAB_STRATEGY_H
#define TAB_STRATEGY_H

#include "strategy.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct tabular_strategy {
    strategy_t **tab;
    int index;
    int size;
} tabular_strategy_t;

tabular_strategy_t *tabular_strategy_create(void);

void tabular_strategy_free(tabular_strategy_t * t);

void tabular_strategy_realloc(tabular_strategy_t * t);

tabular_strategy_t *tabular_strategy_copy(tabular_strategy_t * t);

int tabular_strategy_get_index(tabular_strategy_t * t);

strategy_t *tabular_strategy_get_strategy(tabular_strategy_t * t, int index);

void
tabular_strategy_add_strategy(tabular_strategy_t * t, strategy_t * strategy);

void tabular_strategy_concat(tabular_strategy_t * t1, tabular_strategy_t * t2);

tabular_strategy_t *tabular_strategy_concat_st(tabular_strategy_t * t1,
					       tabular_strategy_t * t2);

int tabular_strategy_fprint (FILE* file, tabular_strategy_t * t);

int tabular_strategy_print(tabular_strategy_t * t);

tabular_strategy_t* tabular_strategy_fscan (FILE* file);

#ifdef __cplusplus
}
#endif

#endif				/* TAB_STRATEGY_H */
