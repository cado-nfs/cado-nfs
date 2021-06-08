#include "cado.h" // IWYU pragma: keep

#include <stdio.h>
#include <stdlib.h>
#include "decomp.h"      // for decomp_t, decomp_create, decomp_free, decomp...
#include "macros.h"
#include "tab_decomp.h"

int decomp_are_equal (decomp_t* t1, decomp_t* t2)
{
    int len1 = t1->len;
    int len2 = t2->len;
    if (len1 != len2)
	return 0;

    for (int i = 0; i < len1; i++)
	if (t1->tab[i]!=t2->tab[i])
	    return 0;
    double prec = 0.1;
    if (t1->nb_elem - t2->nb_elem > prec || t2->nb_elem - t1->nb_elem > prec)
	return 0;
    return 1;
}

int tabular_decomp_are_equal (tabular_decomp_t* t1, tabular_decomp_t* t2)
{
    int len1 = t1->index;
    int len2 = t2->index;
    if (len1 != len2)
	return 0;
    for (int i = 0; i < len1; i++)
	if (!decomp_are_equal(t1->tab[i], t2->tab[i]))
	    return 0;
    return 1;
}

int main()
{
    //create decomp_t!
    int tab[3] = {1,2,3};
    decomp_t* el1 = decomp_create (3, tab, 1000);
    tabular_decomp_t* t = tabular_decomp_create ();
    tabular_decomp_add_decomp (t, el1);
    if (!decomp_are_equal(tabular_decomp_get_decomp(t,0), el1))
	{
	    fprintf (stderr, "error with the test(1)!!!\n");
	    return EXIT_FAILURE;
	}
    tabular_decomp_add_decomp (t, el1);
    //test realloc()
    int tab2[5] = {11,10,9,8,7};
    decomp_t* el2 = decomp_create (5, tab2, 10000);
    tabular_decomp_add_decomp (t, el2);
    if (!decomp_are_equal(tabular_decomp_get_decomp(t,2), el2))
	{
	    fprintf (stderr, "error with the test(2)!!!\n");
	    return EXIT_FAILURE;
	}
    if (decomp_are_equal(tabular_decomp_get_decomp(t,2), el1))
	{
	    fprintf (stderr, "error with the test(3)!!!\n");
	    return EXIT_FAILURE;
	}
    //set and get
    decomp_set_decomp (el1, decomp_get_tab(el2), decomp_get_len(el2));
    decomp_set_nb_elem (el1, decomp_get_nb_elem(el2));
    if (!decomp_are_equal(el2, el1))
	{
	    fprintf (stderr, "error with the test(4)!!!\n");
	    return EXIT_FAILURE;
	}  
    //fprint fscan
 
    // coverity complains about insecure temp files. For tests, I don't
    // think it's a problem, really.
    // coverity[secure_temp]
    FILE* file = tmpfile();
    DIE_ERRNO_DIAG(file == NULL, "tmpfile(%s)", "");
    int errf = tabular_decomp_fprint (file, t) ;
    if (errf < 0)
	{
	    if (errf == -1)
		fprintf (stderr, "error with your 'tab_decomp_t'\n");
	    else //==-2
		fprintf (stderr, "write error on temp file\n");
	    exit (EXIT_FAILURE);
	}
    fseek(file, 0, SEEK_SET);
    tabular_decomp_t* t2 = tabular_decomp_fscan (file);
    if (t2 == NULL)
	{
	    fprintf (stderr, "read error on temp file\n");
	    exit (EXIT_FAILURE);
	}
    fclose (file);

    if (tabular_decomp_are_equal (t2, t) != 1)
	{
	    fprintf (stderr, "error with the test(5)!!!\n");
	    return EXIT_FAILURE;
	}
    decomp_set_nb_elem (t2->tab[1], 21);
    if (tabular_decomp_are_equal (t2, t) != 0)
	{
	    fprintf (stderr, "error with the test(6)!!!\n");
	    return EXIT_FAILURE;
	}
    //free
    tabular_decomp_free (t);
    tabular_decomp_free (t2);
    decomp_free (el1);
    decomp_free (el2);
    
    return EXIT_SUCCESS;
}
