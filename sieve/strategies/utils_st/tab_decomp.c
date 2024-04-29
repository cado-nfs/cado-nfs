#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <regex.h>
#include "tab_decomp.h"
#include "macros.h"

tabular_decomp_t *tabular_decomp_create(void)
{
    tabular_decomp_t *t = malloc(sizeof(*t));
    ASSERT_ALWAYS(t != NULL);

    t->index = 0;
    t->size = 2;

    t->tab = malloc(t->size * sizeof(decomp_t *));
    ASSERT_ALWAYS(t->tab != NULL);

    return t;
}

void tabular_decomp_free(tabular_decomp_t * t)
{
    if (t != NULL)
	{
	    for (int i = 0; i < t->index; i++)
		decomp_free(t->tab[i]);
	    free(t->tab);
	    free(t);
	}
}

void tabular_decomp_realloc(tabular_decomp_t * t)
{
    t->tab = realloc(t->tab, t->size * 2 * (sizeof(decomp_t *)));
    ASSERT_ALWAYS(t->tab != NULL);
    t->size *= 2;
}

void tabular_decomp_add(tabular_decomp_t * t, int len, int *tab, double nb_elem)
{
    if (t->index >= t->size)
	tabular_decomp_realloc(t);
    t->tab[t->index] = decomp_create(len, tab, nb_elem);
    t->index++;
}

void tabular_decomp_add_decomp(tabular_decomp_t * t, decomp_t * decomp)
{
    tabular_decomp_add(t, decomp->len, decomp->tab, decomp->nb_elem);
}

void tabular_decomp_concat(tabular_decomp_t * t1, tabular_decomp_t * t2)
{
    int len = t2->index;
    for (int i = 0; i < len; i++)
	tabular_decomp_add_decomp(t1, t2->tab[i]);
}

decomp_t *tabular_decomp_get_decomp(tabular_decomp_t * t, int index)
{
    return t->tab[index];
}

int tabular_decomp_fprint(FILE * output_file, tabular_decomp_t * t)
{
    for (int i = 0; i < t->index; i++)
	{
	    int err = decomp_fprint(output_file, t->tab[i]);
	    if (err < 0)
		return err;
	}
    return 0;
}

int tabular_decomp_print(tabular_decomp_t * t)
{
    return tabular_decomp_fprint(stdout, t);
}

static decomp_t *analyse_line(char *line)
{
    regex_t preg_decomp;
    const char *str_preg_decomp = "([[:digit:]]+.[[:digit:]]*)";
    regcomp(&preg_decomp, str_preg_decomp, REG_ICASE | REG_EXTENDED);

    //process the ligne
    const char *str_process = &line[0];
    const int len_max = 10;
    char **res = malloc(sizeof(*res) * len_max);
    int ind_res = -1;
    while (str_process[0] != '\0') {
	//printf ("line-->'%s'\n", str_process);
	/*TEST REGULAR EXPRESSION  'preg_decomp */
	size_t nmatch = 2;
	regmatch_t *pmatch = calloc(nmatch, sizeof(*pmatch));
	regexec(&preg_decomp, str_process, nmatch, pmatch, 0);
	if (pmatch[0].rm_so != pmatch[0].rm_eo) {
	    int start = pmatch[1].rm_so;
	    int end = pmatch[1].rm_eo;
	    if (start != -1) {
		int size = end - start;
		char *el = malloc(sizeof(*el) * (size + 1));
		ASSERT_ALWAYS(el != NULL);
		strncpy(el, &str_process[start], size);
		el[size] = '\0';
		ind_res++;
		res[ind_res] = el;	//strtoul (el, NULL, 10);
		//printf ("el = %s\n", el);
	    }
	} else {
	    free(pmatch);
	    break;
	}
	str_process = &str_process[pmatch[0].rm_eo];
	free(pmatch);
    }
    decomp_t *dec = NULL;
    //create decomp_t* dec
    if (ind_res != -1) {
	int tab[ind_res];
	for (int i = 0; i < ind_res; i++) {
	    tab[i] = atoi(res[i]);
	    free(res[i]);
	}
	double nb_elem = strtod(res[ind_res], NULL);
	free(res[ind_res]);
	dec = decomp_create(ind_res, tab, nb_elem);
    }
    //free
    regfree(&preg_decomp);
    free(res);
    return dec;
}

/*
  This function extracts all decompositions of a cofactor from a file.
 */
tabular_decomp_t *tabular_decomp_fscan(FILE * file)
{
    if (file == NULL)
	return NULL;

    tabular_decomp_t *res = tabular_decomp_create();

    const int len_line = 1000;
    char line[len_line];
    while (fgets(line, len_line, file) != 0) {
	decomp_t *dec = analyse_line(line);
	if (dec != NULL) {
	    tabular_decomp_add_decomp(res, dec);
	    decomp_free(dec);
	}
    }
    return res;
}
