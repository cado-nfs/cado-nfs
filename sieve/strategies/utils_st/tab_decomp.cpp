#include "cado.h" // IWYU pragma: keep

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <istream>
#include <ostream>

#include "decomp.hpp"
#include "tab_decomp.hpp"

std::ostream & operator<<(std::ostream & os, tabular_decomp const & t)
{
    for (auto const & D: t)
        os << D << "\n";
    return os;
}

std::istream & operator>>(std::istream & is, tabular_decomp & t)
{
    for(;;) {
        is >> std::ws;
        if (is.eof())
            break;
        decomp D;
        is >> D;
        t.push_back(std::move(D));
    }
    return is;
}

#if 0
static decomp_t *analyse_line(char *line)
{
    regex_t preg_decomp;
    const char *str_preg_decomp = "([[:digit:]]+.[[:digit:]]*)";
    regcomp(&preg_decomp, str_preg_decomp, REG_ICASE | REG_EXTENDED);

    //process the ligne
    const char *str_process = &line[0];
    const int len_max = 10;
    char **res = (char **) malloc(sizeof(*res) * len_max);
    unsigned int ind_res = UINT_MAX;
    while (str_process[0] != '\0') {
	//printf ("line-->'%s'\n", str_process);
	/*TEST REGULAR EXPRESSION  'preg_decomp */
	size_t nmatch = 2;
	regmatch_t *pmatch = (regmatch_t *) calloc(nmatch, sizeof(*pmatch));
	regexec(&preg_decomp, str_process, nmatch, pmatch, 0);
	if (pmatch[0].rm_so != pmatch[0].rm_eo) {
	    int start = pmatch[1].rm_so;
	    int end = pmatch[1].rm_eo;
	    if (start != -1) {
		int size = end - start;
		char *el = (char *) malloc(sizeof(*el) * (size + 1));
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
    if (ind_res != UINT_MAX) {
	unsigned int * tab = (unsigned int *) malloc(ind_res * sizeof(unsigned int));
	for (unsigned int i = 0; i < ind_res; i++) {
            int z = atoi(res[i]);
            ASSERT_ALWAYS(z >= 0);
	    tab[i] = (unsigned int) z;
	    free(res[i]);
	}
	double nb_elem = strtod(res[ind_res], NULL);
	free(res[ind_res]);
	dec = decomp_create(ind_res, tab, nb_elem);
        free(tab);
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
#endif
