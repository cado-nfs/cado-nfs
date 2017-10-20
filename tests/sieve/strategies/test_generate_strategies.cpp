#include "cado.h"

#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "portability.h"
#include "utils.h"
#include "facul.hpp"
#include "pm1.h"
#include "pp1.h"
#include "facul_ecm.h"

#include "generate_strategies.h"
#include "generate_factoring_method.hpp"

static double EPSILON_DBL = LDBL_EPSILON;

//convert type: from strategy_t to facul_strategy_t.
facul_strategy_t* convert_strategy_to_facul_strategy (strategy_t* t)
{
    tabular_fm_t* tab_fm = strategy_get_tab_fm (t);
    int nb_methods = tab_fm->index;

    facul_strategy_t *strategy;
    strategy = (facul_strategy_t*) malloc(sizeof(facul_strategy_t));
    strategy->methods = (facul_method_t*) malloc((nb_methods+1) * sizeof(facul_method_t));
    /*
      without this second method, the function
      facul_clear_strategy (strategy) failed!
    */
    strategy->methods[1].method = 0;
    strategy->methods[1].plan = NULL;

    strategy->lpb = ULONG_MAX;
    strategy->assume_prime_thresh = 0;
    strategy->BBB = 0;

    for (int i = 0; i < nb_methods; i++)
	{
	    fm_t* fm = tab_fm->tab[i];
	    int method = (int)fm->method[0];
	    int curve = (int)fm->method[1];
	    unsigned long B1 = fm->method[2];
	    unsigned long B2 = fm->method[3];
	    strategy->methods[i].method = method;

	    if (method == PM1_METHOD) {
		strategy->methods[i].plan = malloc(sizeof(pm1_plan_t));
		ASSERT(strategy->methods[i].plan != NULL);
		pm1_make_plan((pm1_plan_t*) strategy->methods[i].plan, B1, B2, 0);
	    } else if (method == PP1_27_METHOD || method == PP1_65_METHOD) {
		strategy->methods[i].plan = malloc(sizeof(pp1_plan_t));
		ASSERT(strategy->methods[i].plan != NULL);
		pp1_make_plan((pp1_plan_t*) strategy->methods[i].plan, B1, B2, 0);
	    } else if (method == EC_METHOD) {
		long sigma;
		if (curve == MONTY16) {
		    sigma = 1;
		} else
		    sigma = 2 + rand();

		strategy->methods[i].plan = malloc(sizeof(ecm_plan_t));
		ASSERT((ecm_plan_t*) strategy->methods[i].plan != NULL);
		ecm_make_plan((ecm_plan_t*) strategy->methods[i].plan, B1, B2, curve,
			      labs(sigma), 0, 0);
	    } else {
		exit(EXIT_FAILURE);
	    }
	}
    strategy->methods[nb_methods].method = 0;
    strategy->methods[nb_methods].plan = NULL;
    
    return strategy;
}

int
select_random_index_dec(double sum_nb_elem, tabular_decomp_t* t)
{
    //100000 to consider approximation of distribution
    double alea = rand() % (unsigned long)sum_nb_elem;
    int i = 0;
    double bound = t->tab[0]->nb_elem;
    int len = t->index;
    while (i < (len-1) && (alea - bound) >= 1) {
	i++;
	bound += t->tab[i]->nb_elem;
    }
    return i;
}



/*
  This function compute directly the probabity and the time to find a
  factor in a good decomposition in a cofactor of length r.
 */
weighted_success
bench_proba_time_st(gmp_randstate_t state, facul_strategy_t* strategy,
		    tabular_decomp_t* init_tab, int r, int lpb)
{
    int nb_succes_max = 10000, nb_succes = 0;
    int nb_test = 0;
    int nb_test_max = 10 * nb_succes_max; //20 for five percent
    double time = 0;
    
    double sum_dec = 0;
    for (int i = 0; i < init_tab->index; i++)
	sum_dec += init_tab->tab[i]->nb_elem;
    //struct timeval st_test, end_test;

    std::vector<cxx_mpz> f;

    while (nb_succes < nb_succes_max && (nb_test<nb_test_max))
	{
	    int index = select_random_index_dec(sum_dec, init_tab);
	    int len_p = init_tab->tab[index]->tab[0];

	    cxx_mpz N = generate_composite_integer(state, len_p, r);
            f.clear();

	    double starttime = microseconds();
	    //gettimeofday (&st_test, NULL);
	    int nfound = facul(f, N, strategy);
	    //gettimeofday (&end_test, NULL);
	    double endtime = microseconds();

	    /* time += (end_test.tv_sec - st_test.tv_sec)*1000000L */
	    /* 	+ end_test.tv_usec - st_test.tv_usec; */
      	    time += endtime - starttime;
	    if (nfound != 0) {
                int len_factor = mpz_sizeinbase(f[0], 2);
                if (len_factor <= lpb && (r-len_factor) <= lpb)
                    nb_succes++;
            }
	    nb_test++;
	}

    return weighted_success(nb_succes, time, nb_test);
}


int
get_nb_word (int r)
{
    /*
      We add 0.5 to the length of one word, because for our times the
      lenght is inclusive. For example, if MODREDCUL_MAXBITS = 64
      bits, a cofactor is in one word if is lenght is less OR equal to
      64 bits. So, if you don't add 0.5 to MODREDCUL_MAXBITS, you
      lost the equal and thus insert an error in your maths. 
    */
    double half_word = (MODREDCUL_MAXBITS+0.5)/2.0;
    int number_half_wd = floor(r /half_word);
    int ind = (number_half_wd <2)? 0: number_half_wd - 1;
    return ind;
}

/*
  This function is like bench_time() but only compute the time for one
  length. In fact, i remove the unnecessary computation to reduce the
  cost of this binary.
*/
void bench_time_mini(gmp_randstate_t state, tabular_fm_t * fm, int r)
{
    unsigned long *param;
    fm_t *elem;
    int len = fm->index;	//number of methods!
    //precompute 4 arrays for our bench_time!
    //{{
    int len_n;
    if (r <= MODREDCUL_MAXBITS)
	len_n = MODREDCUL_MAXBITS;
    else if (r <= MODREDC15UL_MAXBITS)
	len_n = MODREDC15UL_MAXBITS;
    else if (r <= MODREDC2UL2_MAXBITS)
	len_n = MODREDC2UL2_MAXBITS;
    else
	len_n = MODREDC2UL2_MAXBITS +30;
    int nb_test = 100000;
    std::vector<cxx_mpz> N;
    N.reserve(nb_test);
    for (int i = 0; i < nb_test; i++) {
	N.push_back(generate_prime_factor(state, len_n));
    }
    //}}
    int ind = get_nb_word(r);
    for (int i = 0; i < len; i++) {
	elem = tabular_fm_get_fm(fm, i);
	param = fm_get_method(elem);
	unsigned long method = param[0];
	unsigned long curve = param[1];
	unsigned long B1 = param[2];
	unsigned long B2 = param[3];
	if (B1 != 0 || B2 != 0) {
	    facul_strategy_t *st = generate_fm(method, curve, B1, B2);

	    double res[4];
            ASSERT_ALWAYS(ind < 4);
	    res[ind] = bench_time_fm_onelength(st, N, nb_test);
	    fm_set_time(elem, res, 4);
	    facul_clear_strategy(st);
	} else {
	    double time[4] = { 0, 0, 0, 0 };
	    fm_set_time(elem, time, 4);
	}
	
    }
}
/*
  This function is like bench_proba() but only compute the necessary
  probabilities. In fact, i remove the unnecessary computation to
  reduce the cost of this binary.
*/

void bench_proba_mini(gmp_randstate_t state, tabular_fm_t * fm, int* val_p,
		      int len_val_p, int len_p_min)
{
    int len = fm->index;	//number of methods!
    int p_max = 100;
    double *proba = (double*) calloc(p_max, sizeof(double));
    ASSERT(proba != NULL);

    unsigned long *param;
    fm_t *elem;

    //{{Will contain the precomputes of composite integer!
    std::vector<std::vector<cxx_mpz>> N(len_val_p);
    int nb_test_max = 10000;

    //}}

    for (int i = 0; i < len; i++) {
	elem = tabular_fm_get_fm(fm, i);
	param = fm_get_method(elem);
	unsigned long method = param[0];
	unsigned long curve = param[1];
	unsigned long B1 = param[2];
	unsigned long B2 = param[3];

	facul_strategy_t *st = generate_fm(method, curve, B1, B2);

	int max_index = 0;
	for (int j = 0; j < len_val_p; j++)
	    {
		int ind_proba = val_p[j] - len_p_min;
		if (B1 == 0 && B2 == 0)
		    proba[ind_proba] = 0;
		else {
		    int len_p = len_p_min + ind_proba;
		    int len_n = 60 + len_p;
		    proba[ind_proba] = bench_proba_fm(st, state, len_p,
						      len_n, N[j],
						      nb_test_max);
		}
		if (ind_proba > max_index)
		    max_index = ind_proba;
	    }
	fm_set_proba(elem, proba, max_index+1, len_p_min);
	//free
	facul_clear_strategy(st);
    }
    free(proba);
}


/************************************************************************/
/*     MAIN                                                             */
/************************************************************************/


int main()
{
    gmp_randstate_t state;
    gmp_randinit_default(state);

    int fbb = 15;//25;
    int lpb = 19;//29;
    int r = 35;//55;

    //input decomp!
    tabular_decomp_t* init_tab = tabular_decomp_create ();
    int tmp1[2] = {17,18};
    decomp_t* el1 = decomp_create (2, tmp1, 100);
    tabular_decomp_add_decomp (init_tab, el1);
    decomp_free (el1);
    int tmp2[2] = {15,20};
    decomp_t* el2 = decomp_create (2, tmp2, 200);
    tabular_decomp_add_decomp (init_tab, el2);
    decomp_free (el2);

    int len_val_factor = 4;
    int val_factor[4] = {15, 17, 18, 20};

    //fm
    fm_t* pm1 = fm_create ();
    unsigned long elem1[4] = {PM1_METHOD, 0, 50, 500};
    fm_set_method (pm1, elem1, 4);
    fm_t* pp1 = fm_create ();
    unsigned long elem2[4] = {PP1_65_METHOD, 0, 70, 700};
    fm_set_method (pp1, elem2, 4);
    fm_t* ecm = fm_create ();
    unsigned long elem3[4] = {EC_METHOD, BRENT12, 80, 1000};
    fm_set_method (ecm, elem3, 4);
    tabular_fm_t* tab = tabular_fm_create();
    tabular_fm_add_fm (tab, pm1);
    tabular_fm_add_fm (tab, pp1);
    tabular_fm_add_fm (tab, ecm);
    //bench our method
    bench_proba_mini(state, tab, val_factor, len_val_factor, fbb);
    //bench_time_mini (state, tab, r);
    /* printf ("mini\n"); */
    /* tabular_fm_print (tab); */

    /* bench_proba (state, tab, fbb); */
    /* bench_time(state, tab); */
    /* printf ("\n all \n"); */
    /* tabular_fm_print (tab); */
    //generate some examples of strategies!
    strategy_t* strat1 =strategy_create ();
    strategy_add_fm (strat1, tab->tab[0]);//pm1
    strategy_add_fm (strat1, tab->tab[1]);//pp1
    strategy_add_fm (strat1, tab->tab[2]);//ecm
    double prob1 = compute_proba_strategy(init_tab, strat1, fbb, lpb);
    //double time1 = compute_time_strategy(init_tab, strat1, r);

    //Create our strategy to use facul().
    facul_strategy_t* st = convert_strategy_to_facul_strategy (strat1);
    //bench our strategies!
    weighted_success res = bench_proba_time_st(state, st, init_tab, r, lpb);
    double prob2 = res.prob;
    //double time2 = res[1];
    
    printf ("prob1 = %lf, prob2 = %lf\n", prob1, prob2);
    // printf ("time1 = %lf, time2 = %lf\n", time1, time2);
    
    double precision_p = 0.05;
    if ((prob1 - prob2) > precision_p || (prob2 - prob1) > precision_p)
    	{
    	    fprintf (stderr, "error with the test(1)\n");
    	    return EXIT_FAILURE;
    	}

    //{{test the function: bench_proba_time_pset()
    int c = tab->tab[0]->method[3]/tab->tab[0]->method[2];
    int param[6] = {(int) tab->tab[0]->method[2], (int) tab->tab[0]->method[2], 1, c,c,1};

    tabular_fm_t* tmp =
    	bench_proba_time_pset (tab->tab[0]->method[0], tab->tab[0]->method[1], state,
    			       17, 20, 18*3, param);

    /*check this probability: the probability to find a prime number
    of length 17 or 18 bits must be more than this to find 18 bits and
    less than this for 18 bits.*/
    fm_t* pm1_bis = tmp->tab[1];
    /* fm_print (tab->tab[0]); */
    /* fm_print (tmp->tab[1]); */
    if (pm1_bis->proba[0] < tab->tab[0]->proba[20-fbb] ||
	pm1_bis->proba[0] > tab->tab[0]->proba[17-fbb])
	{
    	    fprintf (stderr, "error with the test(2)\n");
    	    return EXIT_FAILURE;
	}
    tabular_fm_free (tmp);
    
    //}}
    
    //test the generation of our strategy for one pair of cofactor!!!
    //{{
    strategy_set_proba(strat1, 0.4);
    strategy_set_time(strat1, 40);    
    tabular_strategy_t * strat_r0 = tabular_strategy_create ();
    tabular_strategy_t * strat_r1 = tabular_strategy_create ();

    //firstly, create the zero strategy!
    strategy_t* zero_st = strategy_create();
    unsigned long tab0[4] = {PM1_METHOD, 0, 0, 0};
    fm_t* zero_fm = fm_create();
    fm_set_method (zero_fm, tab0, 4);
    strategy_add_fm (zero_st, zero_fm);
    strategy_set_proba(zero_st, 0.0);    

    tabular_strategy_add_strategy (strat_r0, zero_st);
    tabular_strategy_add_strategy (strat_r0, strat1);
    tabular_strategy_add_strategy (strat_r1, zero_st);
    tabular_strategy_t * res2 = generate_strategy_r0_r1(strat_r0, strat_r1);

    //check proba + time!
    if (!(res2->index == 1 && res2->tab[0]->proba < EPSILON_DBL
	  && res2->tab[0]->time < EPSILON_DBL))
	{
	    fprintf (stderr, "error with the test(3)\n");
	    return EXIT_FAILURE;
	}
    strategy_set_proba(strat_r1->tab[0], 1.0);
    tabular_strategy_t * res3 = generate_strategy_r0_r1(strat_r0, strat_r1);

    //check proba + time!
    if (!(res3->index == 2 && res3->tab[0]->proba < EPSILON_DBL
	  && res3->tab[0]->time < EPSILON_DBL &&
	  (res3->tab[1]->proba - strat1->proba) < EPSILON_DBL &&
	  (res3->tab[1]->time - strat1->time) < EPSILON_DBL))
	{
	    fprintf (stderr, "error with the test(3)\n");
	    return EXIT_FAILURE;
	}
    
    //}}
    //free
    facul_clear_strategy(st);
    strategy_free (strat1);
    tabular_fm_free (tab);
    fm_free (pm1);
    fm_free (pp1);
    fm_free (ecm);
    gmp_randclear (state);
    tabular_decomp_free (init_tab);

    fm_free (zero_fm);
    strategy_free (zero_st);
    tabular_strategy_free (res2);
    tabular_strategy_free (res3);
    tabular_strategy_free (strat_r0);
    tabular_strategy_free (strat_r1);


    return EXIT_SUCCESS;
}
