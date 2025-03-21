#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include "macros.h"
#include "decomp.h"

decomp_t *decomp_create(unsigned int len, const unsigned int *tab, double nb_elem)
{
    decomp_t *t = malloc(sizeof(*t));
    ASSERT_ALWAYS(t != NULL);
    t->len = len;
    t->tab = malloc(t->len * sizeof(unsigned int));
    ASSERT_ALWAYS(t->tab != NULL);
    for (unsigned int i = 0; i < t->len; i++)
	t->tab[i] = tab[i];
    t->nb_elem = nb_elem;
    return t;
}

void decomp_free(decomp_t * t)
{
    if (t != NULL)
	{
	    free(t->tab);
	    free(t);
	}
}

double decomp_get_nb_elem(decomp_t * t)
{
    return t->nb_elem;
}

const unsigned int *decomp_get_tab(decomp_t * t)
{
    return t->tab;
}

unsigned int decomp_get_len(decomp_t * t)
{
    return t->len;
}

void decomp_set_decomp(decomp_t * t, const unsigned int *tab, unsigned int len)
{
    t->len = len;
    CHECKED_REALLOC(t->tab, t->len, unsigned int);
    for (unsigned int i = 0; i < t->len; i++)
	t->tab[i] = tab[i];
}

void decomp_set_nb_elem(decomp_t * t, double nb_elem)
{
    t->nb_elem = nb_elem;
}

int decomp_fprint(FILE * output_file, decomp_t * t)
{
    if (output_file == NULL)
	return -1;
    if (t == NULL)
	return -2;
    fprintf(output_file, "[ ");
    for (unsigned int l = 0; l < t->len - 1; l++)
	fprintf(output_file, "%d, ", t->tab[l]);
    fprintf(output_file, "%d ] : %f\n", t->tab[t->len - 1], t->nb_elem);
    return 0;
}

int decomp_print(decomp_t * t)
{
    return decomp_fprint(stdout, t);
}
