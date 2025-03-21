#ifndef TAB_DECOMP_HPP
#define TAB_DECOMP_HPP

#include <cstdio>
#include "decomp.hpp"

typedef struct tabular_decomp {
    decomp_t **tab;
    int index;
    int alloc;
} tabular_decomp_t;

tabular_decomp_t *tabular_decomp_create();

void tabular_decomp_free(tabular_decomp_t * t);

void tabular_decomp_realloc(tabular_decomp_t * t);

void tabular_decomp_add_decomp(tabular_decomp_t * t, decomp_t * decomp);

void
tabular_decomp_add(tabular_decomp_t * t, unsigned int len, unsigned int *tab, double nb_elem);

decomp_t *tabular_decomp_get_decomp(tabular_decomp_t * t, int index);

void tabular_decomp_concat(tabular_decomp_t * t1, tabular_decomp_t * t2);

int tabular_decomp_fprint(FILE * output_file, tabular_decomp_t * t);

int tabular_decomp_print(tabular_decomp_t * t);

tabular_decomp_t *tabular_decomp_fscan(FILE * file);

#endif				/* TAB_DECOMP_HPP */
