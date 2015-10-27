#include "cado.h"
#include "portability.h"
#include "utils.h"
#include "tab_fm.h"

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "facul.h"

static const double EPSILON_DBL = 0.000001;

tabular_fm_t *tabular_fm_create(void)
{
    tabular_fm_t *t = malloc(sizeof(*t));
    ASSERT(t != NULL);

    t->index = 0;
    t->size = 2;

    t->tab = malloc(t->size * sizeof(fm_t *));
    ASSERT(t->tab != NULL);

    return t;
}

void tabular_fm_free(tabular_fm_t * t)
{
    if (t != NULL)
	{
	    for (int i = 0; i < t->index; i++)	//size
		fm_free(t->tab[i]);
	    free(t->tab);
	    free(t);
	}
}

void tabular_fm_realloc(tabular_fm_t * t)
{
    t->tab = realloc(t->tab, t->size * 2 * (sizeof(fm_t *)));
    ASSERT(t->tab != NULL);
    t->size *= 2;
}

int tabular_fm_get_index(tabular_fm_t * t)
{
    return t->index;
}

void tabular_fm_add_fm(tabular_fm_t * t, fm_t * fm)
{
    tabular_fm_set_fm_index(t, fm, t->index);
}

void
tabular_fm_add(tabular_fm_t * t, unsigned long *method, int len_method,
	       double *proba, int len_proba, double *time, int len_time,
	       int len_p_min)
{
    tabular_fm_set_index(t, method, len_method, proba, len_proba, time,
			 len_time, len_p_min, t->index);
}

void tabular_fm_set_fm_index(tabular_fm_t * t, fm_t * fm, int ind)
{
    tabular_fm_set_index(t, fm->method, fm->len_method, fm->proba,
			 fm->len_proba, fm->time, fm->len_time, fm->len_p_min,
			 ind);
}

void
tabular_fm_set_index(tabular_fm_t * t, unsigned long *method, int len_method,
		     double *proba, int len_proba, double *time, int len_time,
		     int len_p_min, int ind)
{
    if (ind >= t->size)
	tabular_fm_realloc(t);

    if (ind >= t->index) {
	t->tab[ind] = fm_create();
	ASSERT(t->tab[ind] != NULL);
    }

    fm_set_method(t->tab[ind], method, len_method);
    fm_set_proba(t->tab[ind], proba, len_proba, len_p_min);
    fm_set_time(t->tab[ind], time, len_time);
    if (ind >= t->index)
	t->index++;
}

fm_t *tabular_fm_get_fm(tabular_fm_t * t, int index)
{
    if (index >= t->index)
	return NULL;
    return t->tab[index];
}

void tabular_fm_concat(tabular_fm_t * t1, tabular_fm_t * t2)
{
    int len = t2->index;
    for (int i = 0; i < len; i++)
	tabular_fm_add_fm(t1, t2->tab[i]);
}

void tabular_fm_put_zero(tabular_fm_t * t, int index)
{
    if (index < t->index)
	fm_put_zero(t->tab[index]);
}

bool tabular_fm_is_zero(tabular_fm_t * t, int index)
{
    if(index >= t->index)
	return 0;
    return fm_is_zero(t->tab[index]);
}

tabular_fm_t *extract_fm_method(tabular_fm_t * t, int method, int curve)
{
    tabular_fm_t *res = tabular_fm_create();
    int len = t->index;
    for (int i = 0; i < len; i++) {
	fm_t *el = t->tab[i];
	if ((int)el->method[0] == method) {
	    if (method == EC_METHOD) {
		if ((int)el->method[1] == curve)
		    tabular_fm_add_fm(res, el);
	    } else
		tabular_fm_add_fm(res, el);
	}
    }
    return res;
}

/************************************************************************/
/*           PRINT AND SCAN OUR FILES OF FACTORING METHODS              */
/************************************************************************/

int tabular_fm_print(tabular_fm_t * t)
{
    return tabular_fm_fprint(stdout, t);
}

int tabular_fm_fprint(FILE * file, tabular_fm_t * t)
{
    int len = t->index;
    for (int i = 0; i < len; i++) {
	fm_t *elem = tabular_fm_get_fm(t, i);
	if (fm_fprint(file, elem) < 0)
	    return -1;
    }
    return 0;
}

static int is_elem(const char c)
{
    return (c >= 48 && c <= 57) || c == '|';
}

static void next_elem(FILE * file, int *current_char)
{
    //end the current element
    while (*current_char != EOF && is_elem(*current_char))
	*current_char = fgetc(file);

    //find the next element
    while (*current_char != EOF && !is_elem(*current_char)) {
	*current_char = fgetc(file);
    }
    fseek(file, -1, SEEK_CUR);
}

static fm_t *sub_routine_fm_fscanf(FILE * file, int *current_char)
{
    fm_t *fm = fm_create();

    int len_method = 100;
    int len_proba = 100;
    int len_time = 100;

    unsigned long method[len_method];
    double proba[len_proba];
    double time[len_time];

    int ind = 0;
    int rc;
    while (ind < len_method && *current_char != '|') {
	rc = fscanf(file, "%lu", &method[ind++]);
        ASSERT_ALWAYS(rc == 1);
	next_elem(file, current_char);
    }
    next_elem(file, current_char);

    fm_set_method(fm, method, ind);

    rc = fscanf(file, "%d", &fm->len_p_min);
    ASSERT_ALWAYS(rc == 1);
    next_elem(file, current_char);
    ind = 0;
    while (ind < len_proba && *current_char != '|') {
	rc = fscanf(file, "%lf", &proba[ind++]);
        ASSERT_ALWAYS(rc == 1);
	next_elem(file, current_char);
    }
    next_elem(file, current_char);

    fm_set_proba(fm, proba, ind, fm->len_p_min);

    ind = 0;
    while (ind < len_time && *current_char != '|') {
	rc = fscanf(file, "%lf", &time[ind++]);
        ASSERT_ALWAYS(rc == 1);
	next_elem(file, current_char);
    }
    next_elem(file, current_char);

    fm_set_time(fm, time, ind);

    return fm;
}

tabular_fm_t* tabular_fm_fscan(FILE * file)
{
    if (file == NULL)
	return NULL;
    tabular_fm_t * res = tabular_fm_create ();
    int current_char = fgetc(file);
    fseek(file, -1, SEEK_CUR);
    while (current_char != EOF) {
	fm_t *fm = sub_routine_fm_fscanf(file, &current_char);
	tabular_fm_add_fm(res, fm);
	fm_free(fm);
    }

    return res;
}

/************************************************************************/
/*                      SORT_TAB_FM                                     */
/************************************************************************/

//return a positive value if el1 is greater than el2 and a negative
//value otherwise.
int fm_cmp(fm_t * el1, fm_t * el2)
{
    if (fm_is_zero(el1))
	return -1;
    else if (fm_is_zero(el2))
	return 1;
    /*assume that the variable len_p_min is the same for the both
      fm_t.*/
    //compare the probabilities!
    int len1 = el1->len_proba;
    int len2 = el2->len_proba;

    int len = (len1 < len2) ? len1 : len2;
    double diff_proba = 0;
    for (int i = 0; i < len; i++) 
	if (el1->proba[i] > EPSILON_DBL && el2->proba[i] > EPSILON_DBL) 
	    diff_proba += el1->proba[i] - el2->proba[i];

    return (diff_proba>EPSILON_DBL)? 1: -1;
}

void fm_swap(tabular_fm_t * t, int index1, int index2)
{
    fm_t *c = t->tab[index1];
    t->tab[index1] = t->tab[index2];
    t->tab[index2] = c;
}

static void tabular_fm_sort_rec(tabular_fm_t * t, int begin, int end)
{
    int index_max = begin;
    for (int i = begin; i < end; i++) {
	if (fm_cmp(t->tab[i], t->tab[index_max]) > 0) {
	    index_max = i;
	}
    }
    fm_swap(t, end - 1, index_max);
}

void tabular_fm_sort(tabular_fm_t * t)
{
    int max = t->index;
    while (max > 0) {
	tabular_fm_sort_rec(t, 0, max);
	max--;
    }
}

