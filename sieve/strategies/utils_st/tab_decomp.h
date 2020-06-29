#ifndef TAB_DECOMP
#define TAB_DECOMP

#include <stdio.h>
#include "decomp.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct tabular_decomp {
    decomp_t **tab;
    int index;
    int size;
} tabular_decomp_t;

tabular_decomp_t *tabular_decomp_create();

void tabular_decomp_free(tabular_decomp_t * t);

void tabular_decomp_realloc(tabular_decomp_t * t);

void tabular_decomp_add_decomp(tabular_decomp_t * t, decomp_t * decomp);

void
tabular_decomp_add(tabular_decomp_t * t, int len, int *tab, double nb_elem);

decomp_t *tabular_decomp_get_decomp(tabular_decomp_t * t, int index);

void tabular_decomp_concat(tabular_decomp_t * t1, tabular_decomp_t * t2);

int tabular_decomp_fprint(FILE * output_file, tabular_decomp_t * t);

int tabular_decomp_print(tabular_decomp_t * t);

tabular_decomp_t *tabular_decomp_fscan(FILE * file);

#ifdef __cplusplus
}
#endif

#endif				/* TAB_DECOMP */
