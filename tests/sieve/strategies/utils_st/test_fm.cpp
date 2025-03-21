#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>              /* for _O_BINARY */
#include "fm.h"      // for fm_set_proba, fm_set_method, fm_set_time, fm_free
#include "macros.h"
#include "tab_fm.h"

//test equality between two fm!  check all parameters in the structure
//fm, while fm_is_equal just see the parameter method!

int fm_are_equals (fm_t* t1, fm_t* t2){
    //test method
    int len_met1 = fm_get_len_method (t1);
    int len_met2 = fm_get_len_method (t2);
    if (len_met2 != len_met1)
	return 0;
    unsigned long* tab_t1 = fm_get_method (t1);
    unsigned long* tab_t2 = fm_get_method (t2);
    for (int i = 0; i < len_met1; i++)
	if (tab_t1[i] != tab_t2[i])
	    return 0;
    //test proba
    int len_prob1 = fm_get_len_proba (t1);
    int len_prob2 = fm_get_len_proba (t2);
    if (len_prob2 != len_prob1)
	return 0;
    double* tabp_t1 = fm_get_proba (t1);
    double* tabp_t2 = fm_get_proba (t2);
    double prec = 0.0001;
    for (int i = 0; i < len_prob1; i++)
	if (tabp_t1[i]- tabp_t2[i]> prec || tabp_t2[i]- tabp_t1[i]> prec)
	    return 0;
    //test time
    int len_time1 = fm_get_len_time (t1);
    int len_time2 = fm_get_len_time (t2);
    if (len_time2 != len_time1)
	return 0;
    double* tabt_t1 = fm_get_time (t1);
    double* tabt_t2 = fm_get_time (t2);
    for (int i = 0; i < len_time1; i++)
	if (tabt_t1[i]- tabt_t2[i]> prec || tabt_t2[i]- tabt_t1[i]> prec)
	    return 0;
    return 1;
}


//test equality between two tab_fm!
int tabular_fm_are_equals (tabular_fm_t* t1, tabular_fm_t* t2){
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


//check if the tab_fm is really sorted!
int check_sort (tabular_fm_t* tab)
{
    int len = tab->index;
    if (len == 1)
	return 1;
    for (int i = 0; i < len-1; i++)
	{
	    if (fm_cmp (tab->tab[i], tab->tab[i+1]) > 0)
		return 0;
	}
    return 1;
}



int main ()
{
#ifdef HAVE_MINGW
    _fmode = _O_BINARY;		/* Binary open for all files */
#endif

    //create a first fm!
    fm_t* t = fm_create();
    unsigned long elem[4];
    for (int i=0; i < 4; i++)
	elem[i] = i;
    fm_set_method (t, elem, 4);

    double elem2[10];
    for (int i=0; i < 10; i++)
	elem2[i] = i*0.01;
    fm_set_proba (t, elem2, 10, 0);
    
    double elem3[10];
    for (int i=0; i < 10; i++)
	elem3[i] = i*100;
    fm_set_time (t, elem3, 10);
    //test tab
    //test add
    tabular_fm_t* tab1 = tabular_fm_create();
    tabular_fm_add_fm (tab1, t);
    //test fm_set_fm_index
    tabular_fm_set_fm_index(tab1, t, 1);
    if (!fm_are_equals(tabular_fm_get_fm(tab1, 0), tabular_fm_get_fm(tab1, 1)))
	{
	    fprintf (stderr, "error with the test(1)!!!\n");
	    return EXIT_FAILURE;
	}
    //test concat
    double elem4[5] = {0.2,0.3,0.9,0.99,0.999};
    fm_set_proba (t, elem4, 5, 0);
    tabular_fm_t* tab2 = tabular_fm_create();
    tabular_fm_add_fm (tab2, t);
    tabular_fm_concat (tab1, tab2);

    if (!fm_are_equals(tabular_fm_get_fm(tab2, 0),
		       tabular_fm_get_fm(tab1, tab1->index - tab2->index )))
	{
	    fprintf (stderr, "error with the test(2)!!!\n");
	    return EXIT_FAILURE;
	}
    //test copy
    fm_t* t2 = fm_copy(tabular_fm_get_fm(tab1, 0));
    if (!fm_are_equals(tabular_fm_get_fm(tab1, 0), t2))
	{
	    fprintf (stderr, "error with the test(3)!!!\n");
	    return EXIT_FAILURE;
	}  
    if (fm_are_equals(tabular_fm_get_fm(tab2, 0), t2))
	{
	    fprintf (stderr, "error with the test(4)!!!\n");
	    return EXIT_FAILURE;
	}  
    //test fprint fscan
 
    // coverity complains about insecure temp files. For tests, I don't
    // think it's a problem, really.
    // coverity[secure_temp]
    FILE* file = tmpfile();
    DIE_ERRNO_DIAG(file == NULL, "tmpfile(%s)", "");
    int errf = (tabular_fm_fprint (file, tab1) == -1);
    if (errf)
	{
            fprintf (stderr, "write error on temp file\n");
	    exit (EXIT_FAILURE);
	}
    fseek(file, 0, SEEK_SET);
    tabular_fm_t* tab3 = tabular_fm_fscan (file);
    if (tab3 == NULL)
	{
            fprintf (stderr, "read error on temp file\n");
	    exit (EXIT_FAILURE);
	}
    fclose (file);
    if (tabular_fm_are_equals (tab3, tab1) != 1)
	{
	    fprintf (stderr, "error with the test(5)!!!\n");
	    return EXIT_FAILURE;
	}

    //sort (and swap);
    //random test
    tabular_fm_t* tab4 = tabular_fm_create();
    fm_set_method (t, elem, 4);
    fm_set_time (t, elem3, 10);
    for (int i = 0; i < 10; i++)
    	{
    	    double elem_p[5];
    	    for (int j = 0; j < 5; j++)
    		elem_p[j] = (rand()%1000)/1000;
    	    fm_set_proba (t, elem_p, 5,0);
    	    tabular_fm_add_fm (tab4, t);
    	}
    tabular_fm_sort (tab4);
    if (!check_sort (tab4))
    	{
    	    fprintf (stderr, "error with the test(6)!!!\n");
    	    return EXIT_FAILURE;
    	}

    //constant test
    tabular_fm_t* tab5 = tabular_fm_create();
    fm_set_method (t, elem, 4);
    fm_set_time (t, elem3, 10);
    double elem_p1[5] = {0.9, 0.8,0.7,0.6,0.3};
    fm_set_proba (t, elem_p1, 5,0);
    tabular_fm_add_fm (tab5, t);
    double elem_p2[5] = {0.8999, 0.8001,0.6999,0.55,0.2999};
    fm_set_proba (t, elem_p2, 5,0);
    tabular_fm_add_fm (tab5, t);
    double elem_p3[5] = {0.5, 0.399,0.19,0.0006,0.00003};
    fm_set_proba (t, elem_p3, 5,0);
    tabular_fm_add_fm (tab5, t);
    double elem_p4[5] = {0.5, 0.399,0.19,0.0006,0.00003};
    fm_set_proba (t, elem_p4, 5,0);
    tabular_fm_add_fm (tab5, t);
    double elem_p5[5] = {0.2, 0.18,0.07,0.0006,0.0000003};
    fm_set_proba (t, elem_p5, 5,0);
    tabular_fm_add_fm (tab5, t);
    double elem_p6[5] = {0.6, 0.5,0.4,0.3,0.2};
    fm_set_proba (t, elem_p6, 5,0);
    tabular_fm_add_fm (tab5, t);
    tabular_fm_sort (tab5);
    if (!check_sort (tab5))
    	{
    	    fprintf (stderr, "error with the test(7)!!!\n");
    	    return EXIT_FAILURE;
    	}
    fm_swap(tab5, 3, 1);
    if (check_sort (tab5))
    	{
    	    fprintf (stderr, "error with the test(8)!!!\n");
    	    return EXIT_FAILURE;
    	}
    //free
    tabular_fm_free (tab1);
    tabular_fm_free (tab2);
    tabular_fm_free (tab3);
    tabular_fm_free (tab4);
    tabular_fm_free (tab5);
    fm_free (t);
    fm_free (t2);
    
	
    return EXIT_SUCCESS;
}



