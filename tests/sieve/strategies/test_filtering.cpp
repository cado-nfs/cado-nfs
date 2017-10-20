#include "cado.h"
#include "generate_factoring_method.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//check if the spread of remaining methods is homogeneous.
int check_filt (tabular_fm_t* res, int init_nb_method)
{
    int len = res->index;
    if (len <= 2)//because less than two points isn't really representative!!!
	return 1;
    double aver_gap_prat = 0;
    double aver_gap_theo = init_nb_method/(double)len;
    for (int i = 0; i < len-1; i++)
	aver_gap_prat += (res->tab[i+1]->method[0] - res->tab[i]->method[0]);

    aver_gap_prat /= (len-1);
    double perc_gap = (aver_gap_prat - aver_gap_theo)/aver_gap_theo;
    if (perc_gap < 0)
	perc_gap *= -1;
    /* printf ("per = %lf\n", perc_gap); */
    /* printf ("aver_gap_theo = %lf\n", aver_gap_theo); */
    /* printf ("aver_gap_prat = %lf\n", aver_gap_prat); */
    if (perc_gap > 0.4)
	return 0;
    return 1;
}

int main ()
{
    unsigned int nb_fm = 10;
    int final_nb_fm = 4;

    tabular_fm_t* tab = tabular_fm_create ();

    fm_t* fm = fm_create();
    for (unsigned int i = 0; i < nb_fm; i++)
	{
	    unsigned long elem[4] = {i, 0, i*(1+rand()%10),  i*(1+rand()%10)};
	    fm_set_method (fm, elem, 4);
	    //proba
	    int len_proba = 4;
	    double proba[len_proba];
	    for (int j = 0; j < len_proba; j++)
		proba[j] = i/((double)nb_fm+1)+0.01*j*i;
	    fm_set_proba (fm, proba, len_proba, 0);
	    //time
	    double time[4];
	    for (int j = 0; j < 4; j++)
		time[j] = (i*i*i)*(pow(10,j));
	    fm_set_time (fm, time, 4);
	    //add
	    tabular_fm_add_fm (tab, fm);
	}
    fm_free (fm);
    tabular_fm_t* res = filtering (tab, final_nb_fm);
    //tabular_fm_print (res);
    int err = 0;
    if (!check_filt (res, nb_fm))
	err = 1;
    //free
    tabular_fm_free (tab);
    tabular_fm_free (res);

    return err?EXIT_FAILURE:EXIT_SUCCESS;
}
