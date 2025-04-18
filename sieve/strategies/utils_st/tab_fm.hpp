#ifndef TAB_FM_HPP
#define TAB_FM_HPP

#include <cstdio>

#include "fm.hpp"

typedef struct tabular_fm {
    fm_t ** tab;
    unsigned int size;
    unsigned int alloc;
} tabular_fm_t;

tabular_fm_t * tabular_fm_create();

void tabular_fm_free(tabular_fm_t * t);

void tabular_fm_realloc(tabular_fm_t * t);

unsigned int tabular_fm_get_size(tabular_fm_t const * t);

fm_t * tabular_fm_get_fm_rw(tabular_fm_t * t, unsigned int index);
fm_t const * tabular_fm_get_fm(tabular_fm_t const * t, unsigned int index);

void tabular_fm_add_fm(tabular_fm_t * t, fm_t const * fm);

void tabular_fm_add(tabular_fm_t * t,
        unsigned long const * method, unsigned int len_method,
        double const * proba, unsigned int len_proba,
        double const * time, unsigned int len_time,
        unsigned int len_p_min);

void tabular_fm_set_fm_index(tabular_fm_t * t, fm_t const * fm, unsigned int ind);

void tabular_fm_set_index(tabular_fm_t * t,
        unsigned long const * method, unsigned int len_method,
        double const * proba, unsigned int len_proba,
        double const * time, unsigned int len_time,
        unsigned int len_p_min, unsigned int ind);

/* concatenate t2 at the end of t1 */
void tabular_fm_concat(tabular_fm_t * t1, tabular_fm_t * t2);

void tabular_fm_put_zero(tabular_fm_t * t, unsigned int index);

bool tabular_fm_is_zero(tabular_fm_t const * t, unsigned int index);

tabular_fm_t * extract_fm_method(tabular_fm_t const * t, int method, int curve);

int tabular_fm_fprint(FILE * output_file, tabular_fm_t const * t);

int tabular_fm_print(tabular_fm_t const * t);

tabular_fm_t * tabular_fm_fscan(FILE * file);

/************************************************************************/
/*                      SORT_TAB_FM                                     */
/************************************************************************/
/* the comparative is according to the probabilities! */
int fm_cmp(fm_t * el1, fm_t * el2);

void fm_swap(tabular_fm_t * t, unsigned int index1, unsigned int index2);

// according to the probabilities!
void tabular_fm_sort(tabular_fm_t * t);

#endif /* TAB_FM_HPP */
