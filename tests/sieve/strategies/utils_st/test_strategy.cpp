#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include "fm.h"            // for fm_set_method, fm_set_proba, fm_set_time
#include "strategy.h"      // for strategy_t, strategy_add_fm, strategy_free
#include "tab_fm.h"        // for tabular_fm_get_fm, tabular_fm_t
#include "tab_strategy.h"
#include "macros.h"

//test equality between two fm!
int fm_are_equals (fm_t* t1, fm_t* t2){
    //test just method because fscan and fprint strategy don't take in
    //consideration the probabilities and the time for each factoring
    //method.
    int len_met1 = fm_get_len_method (t1);
    int len_met2 = fm_get_len_method (t2);
    if (len_met2 != len_met1)
	return 0;
    unsigned long* tab_t1 = fm_get_method (t1);
    unsigned long* tab_t2 = fm_get_method (t2);
    for (int i = 0; i < len_met1; i++)
	if (tab_t1[i] != tab_t2[i])
	    return 0;
    return 1;
}


//test equality between two tab_fm!
int tab_fm_are_equals (tabular_fm_t* t1, tabular_fm_t* t2){
    int len1 = t1->index;
    int len2 = t2->index;
    if (len1 != len2)
	return 0;
    
    for (int i = 0; i < len1; i++)
	if (!fm_are_equals (tabular_fm_get_fm (t1, i),
			    tabular_fm_get_fm (t2, i)))
	    return 0;
    return 1;
}
//test equality between two strategies!
int strategies_are_equals (strategy_t* t1, strategy_t* t2)
{
    //equality between two tab_fm
    if ( !tab_fm_are_equals (strategy_get_tab_fm (t1),
			     strategy_get_tab_fm (t2)))
	return 0;

    //equality probab. and time
    double p1 = strategy_get_proba (t1);
    double p2 = strategy_get_proba (t2);
    double prec = 0.0001;
    if (p1-p2 > prec || p2-p1 > prec)
	return 0;
    return 1;
}

//test equality between two tabular_strategy_t!
int tabular_strategies_are_equals (tabular_strategy_t* t1,
				   tabular_strategy_t* t2)
{
    //equality between two tab_strategy
    int len1 = t1->index;
    int len2 = t2->index;
    if (len1 != len2)
	return 0;
    for (int i = 0; i < len1; i++)
	{
	    if (!strategies_are_equals (t1->tab[i], t2->tab[i]))
		return 0;
	}
    return 1;
}


int main()
{
    //test strategy
    strategy_t* t1 = strategy_create ();
    fm_t* elem = fm_create ();
    unsigned long value[4] = {1,2,3,4};
    double p = 0.5;
    double a = 2;
    fm_set_method (elem, value,4);
    fm_set_proba (elem, &p, 1,0);
    fm_set_time (elem, &a, 1);
    strategy_add_fm (t1, elem);

    p = 0.5;
    a = 21;
    fm_set_method (elem, value,4);
    fm_set_proba (elem, &p, 1,0);
    fm_set_time (elem, &a, 1);
    strategy_add_fm (t1, elem);

    unsigned long value2[6] = {1,0,10,103};
    p = 0.7;
    a = 221;
    fm_set_method (elem, value2,4);
    fm_set_proba (elem, &p, 1,0);
    fm_set_time (elem, &a, 1);

    strategy_add_fm (t1, elem);
    strategy_set_proba(t1, 0.5);
    strategy_set_time(t1, 23.458);

    //test strategy_copy()
    strategy_t* t2 = strategy_copy(t1);

    //t1 equal to t2
    if (strategies_are_equals (t1, t2) != 1)
	{
	    fprintf (stderr, "error with the test(1)!!!\n");
	    return EXIT_FAILURE;
	}
    strategy_set_proba(t2, 0.6);
    if (strategies_are_equals (t1, t2) != 0)
	{
	    fprintf (stderr, "error with the test(2)!!!\n");
	    return EXIT_FAILURE;
	}

    //test side
    strategy_add_fm_side(t2, elem, 1);
    if (t2->side[t2->len_side-1]!=1)
	{
	    fprintf (stderr, "error with the test(3)!!!\n");
	    return EXIT_FAILURE;
	}

    //Test tabular_strategy
    //add t and t2 in tab
    tabular_strategy_t *tab = tabular_strategy_create();
    tabular_strategy_add_strategy (tab, t1);
    if(strategies_are_equals (t1, tab->tab[0]) != 1)
	{
	    fprintf (stderr, "error with the test(4)!!!\n");
	    return EXIT_FAILURE;
	}

    tabular_strategy_add_strategy (tab, t2);
    if (strategies_are_equals (t1, tab->tab[1]) != 0)
	{
	    fprintf (stderr, "error with the test(5)!!!\n");
	    return EXIT_FAILURE;
	}

    //test print and scan

    // coverity complains about insecure temp files. For tests, I don't
    // think it's a problem, really.
    // coverity[secure_temp]
    FILE* file = tmpfile();
    DIE_ERRNO_DIAG(file == NULL, "tmpfile(%s)", "");
    int errf = (tabular_strategy_fprint (file, tab) == -1);
    if (errf)
	{
            fprintf (stderr, "write error on temp file\n");
	    exit (EXIT_FAILURE);
	}
    fseek(file, 0, SEEK_SET);
    tabular_strategy_t* tab2 = tabular_strategy_fscan (file);
    if (tab2 == NULL)
	{
            fprintf (stderr, "read error on temp file\n");
	    exit (EXIT_FAILURE);
	}
    fclose (file);

    if (tabular_strategies_are_equals (tab2, tab) != 1)
	{
	    fprintf (stderr, "error with the test(6)!!!\n");
	    return EXIT_FAILURE;
	}
    //free 
    fm_free (elem);
    strategy_free (t1);
    strategy_free (t2);
    tabular_strategy_free (tab);
    tabular_strategy_free (tab2);
    return EXIT_SUCCESS;
}

