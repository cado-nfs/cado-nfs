#ifndef TAB_FM
#define TAB_FM

#include <stdio.h>
#ifndef __cplusplus
#include <stdbool.h>
#endif

#include "fm.h"

typedef struct tabular_fm {
    fm_t **tab;
    int index;
    int size;
} tabular_fm_t;

#ifdef __cplusplus
extern "C" {
#endif


tabular_fm_t *tabular_fm_create(void);

void tabular_fm_free(tabular_fm_t * t);

void tabular_fm_realloc(tabular_fm_t * t);

int tabular_fm_get_index(tabular_fm_t * t);

fm_t *tabular_fm_get_fm(tabular_fm_t * t, int index);

void tabular_fm_add_fm(tabular_fm_t * t, fm_t * fm);

void
tabular_fm_add(tabular_fm_t * t, unsigned long *method, int len_method,
	       double *proba, int len_proba, double *time, int len_time,
	       int len_p_min);

void tabular_fm_set_fm_index(tabular_fm_t * t, fm_t * fm, int ind);

void
tabular_fm_set_index(tabular_fm_t * t, unsigned long *method, int len_method,
		     double *proba, int len_proba, double *time, int len_time,
		     int len_p_min, int ind);

/* concatenate t2 at the end of t1 */
void tabular_fm_concat(tabular_fm_t * t1, tabular_fm_t * t2);

void tabular_fm_put_zero(tabular_fm_t * t, int index);

bool tabular_fm_is_zero(tabular_fm_t * t, int index);

tabular_fm_t *extract_fm_method(tabular_fm_t * t, int method, int curve);

int tabular_fm_fprint(FILE * output_file, tabular_fm_t * t);

int tabular_fm_print(tabular_fm_t * t);

 tabular_fm_t * tabular_fm_fscan(FILE * file);

/************************************************************************/
/*                      SORT_TAB_FM                                     */
/************************************************************************/
/* the comparative is according to the probabilities! */
int fm_cmp(fm_t * el1, fm_t * el2);

void fm_swap(tabular_fm_t * t, int index1, int index2);

//according to the probabilities!
void tabular_fm_sort(tabular_fm_t * t);

#ifdef __cplusplus
}
#endif
#endif				/* TAB_FM */
