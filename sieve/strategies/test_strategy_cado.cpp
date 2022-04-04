#include "cado.h" // IWYU pragma: keep
#include <array>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cstring>
#include <gmp.h>
#include "pm1.h"
#include "pp1.h"
#include "facul_ecm.h"
#include "facul.hpp"
#include "facul_fwd.hpp"
#include "fm.h" // fm_t fm_set_method
#include "finding_good_strategy.h"
#include "generate_factoring_method.hpp"
#include "generate_strategies.h"
#include "tab_strategy.h"
#include "tab_fm.h"
#include "decomp.h"     // decomp_t
#include "tab_decomp.h"
#include "cxx_mpz.hpp"  // cxx_mpz
#include "timing.h"  // microseconds
#include "params.h"     // param_list
#include "macros.h"
#include "modredc_ul.h" // MODREDCUL_MAXBITS
#include "stage2.h" // stage2_plan_t
#include "strategy.h"                     // for strategy_t, strategy_add_fm
#include "misc.h"   // u64_random  // IWYU: keep

int CONST_TEST_R = 55;

//#define CADO_INTERLEAVING
/*
  This binary allows to test our procedure choosing optimal
  strategies. In fact, using the strategy of CADO, we can use it to
  approximate the "theorical" number of relations found per second found
  with certains parameters. Then, comparing the "theorical" value with the
  real value, we could detect if a problem exists in our procedure.
*/

/* There's some dead code in here, which is never called. If we ever come
 * back to this strategy business, there might be something to look at in
 * here, but the odds are low.
 */
#define COMPILE_DEAD_CODE

#ifdef COMPILE_DEAD_CODE
#include "facul.hpp"

static MAYBE_UNUSED tabular_fm_t *generate_methods_cado(const int lpb)
{
    /* we set mfb = 3*lpb to avoid the special case of 2 large primes */
    int n = nb_curves (lpb, 3 * lpb);
    tabular_fm_t *res = tabular_fm_create();
    fm_t *fm = fm_create();
    unsigned long method[4];

    /* run one P-1 curve with B1=315 and B2=2205 */
    method[0] = PM1_METHOD;	//method
    method[1] = 0;		//curve
    method[2] = 315;		//B1
    method[3] = 2205;		//B2
    fm_set_method(fm, method, 4);
    tabular_fm_add_fm(res, fm);

    /* run one P+1 curve with B1=525 and B2=3255 */
    method[0] = PP1_27_METHOD;	//method
    method[1] = 0;		//curve
    method[2] = 525;		//B1
    method[3] = 3255;		//B2
    fm_set_method(fm, method, 4);
    tabular_fm_add_fm(res, fm);

    /* run one ECM curve with Montgomery parametrization, B1=105, B2=3255 */
    method[0] = EC_METHOD;	//method
    method[1] = MONTY12;	//curve
    method[2] = 105;		//B1
    method[3] = 3255;		//B2
    fm_set_method(fm, method, 4);
    tabular_fm_add_fm(res, fm);

    if (n > 0) {
	method[0] = EC_METHOD;	//method
	method[1] = BRENT12;	//curve
	method[2] = 315;	//B1
	method[3] = 5355;	//B2
	fm_set_method(fm, method, 4);
	tabular_fm_add_fm(res, fm);
    }

    /* heuristic strategy where B1 is increased by sqrt(B1) at each curve */
    double B1 = 105.0;
    for (int i = 4; i < n + 3; i++) {
	double B2;
	unsigned int k;
	B1 += sqrt(B1);//TEst
	B2 = 17.0 * B1;
	/* we round B2 to (2k+1)*105, thus k is the integer nearest to
	   B2/210-0.5 */
	k = B2 / 210.0;
	method[0] = EC_METHOD;	//method
	method[1] = MONTY12;	//curve
	method[2] = B1;
	method[3] = (2 * k + 1) * 105;	//B2
	//printf ("B1 = %lf, %d\n", B1, (2 * k + 1) * 105);
	//getchar();
	fm_set_method(fm, method, 4);
	tabular_fm_add_fm(res, fm);
    }
    fm_free(fm);
    return res;
}
#endif  /* COMPILE_DEAD_CODE */

/*
This function generates the strategy of cado and computes the
  probability and the time to find a prime divisor in a cofactor of
  'r' bits with the bound fbb and lpb.  This strategy is the
  concatenation of all methods in 'methods'.
*/
static tabular_strategy_t *generate_strategy_cado(tabular_fm_t * methods,
					   tabular_decomp_t * tab_dec, int fbb,
					   int lpb, int r)
{
    tabular_strategy_t *tab_strat = tabular_strategy_create();
    strategy_t *strat = strategy_create();

    int lim = 2 * fbb - 1;

    ASSERT_ALWAYS((tab_dec == NULL) == (r < lim));

    if (r < lim) {
	fm_t *zero = fm_create();
	unsigned long method[4] = { PM1_METHOD, 0, 0, 0 };
	fm_set_method(zero, method, 4);
	strategy_add_fm(strat, zero);
	strategy_set_time(strat, 0.0);

	if (r != 1 && (r < fbb || r > lpb))
	    strategy_set_proba(strat, 0.0);
	else			// r==1 or fbb<= r0 <= lpb
	    strategy_set_proba(strat, 1.0);
	fm_free(zero);
    } else {
        /* we set mfb = 3*lpb to avoid the special case of 2 large primes */
        int len = 3 + nb_curves (lpb, 3 * lpb);
        //printf ("len  = %d\n", len);
        ASSERT(len <= methods->index);
        for (int i = 0; i < len; i++)
            strategy_add_fm(strat, methods->tab[i]);

        //eval
        double p = compute_proba_strategy(tab_dec, strat, fbb, lpb);

        double t = compute_time_strategy(tab_dec, strat, r);
        /* if (r == CONST_TEST_R) */
        /*   { */
        /*     printf ("p =  %lf, t = %lf\n", p, t); */
        /*     getchar(); */
        /*   } */
        strategy_set_proba(strat, p);
        strategy_set_time(strat, t);
    }

    tabular_strategy_add_strategy(tab_strat, strat);
    strategy_free(strat);

    return tab_strat;
}



/*
  generate the matrix with the strategy of CADO.
*/

static tabular_strategy_t ***generate_matrix_cado(const char *name_directory_decomp,
					   tabular_fm_t * methods,
					   unsigned long lim0, int lpb0, int mfb0,
					   unsigned long lim1, int lpb1, int mfb1)
{
    /*
       allocates the matrix which contains all optimal strategies for
       each couple (r0, r1): r0 is the lenght of the cofactor in
       the first side, and r1 for the second side.
     */
    tabular_strategy_t ***matrix = (tabular_strategy_t***) malloc(sizeof(*matrix) * (mfb0 + 1));
    ASSERT(matrix != NULL);

    for (int r0 = 0; r0 <= mfb0; r0++) {
	matrix[r0] = (tabular_strategy_t**) malloc(sizeof(tabular_strategy_t *) * (mfb1 + 1));
	ASSERT(matrix[r0] != NULL);
    }

    int fbb0 = ceil (log2 ((double) (lim0 + 1)));
    int fbb1 = ceil (log2 ((double) (lim1 + 1)));

    /*
       For each r0 in [fbb0..mfb0], we precompute and strore the data
       for all strategies.  Whereas, for each r1, the values will be
       computed sequentially. 
     */
    fm_t *zero = fm_create();
    unsigned long method_zero[4] = { 0, 0, 0, 0 };
    fm_set_method(zero, method_zero, 4);

    tabular_strategy_t **data_rat = (tabular_strategy_t**) malloc(sizeof(tabular_strategy_t *) * (mfb0 + 1));
    ASSERT(data_rat);

    int lim = 2 * fbb0 - 1;
    for (int r0 = 0; r0 <= mfb0; r0++) {
	tabular_decomp_t *tab_decomp = NULL;
	if (r0 >= lim) {
	    char name_file[200];
            snprintf(name_file, sizeof(name_file),
		    "%s/decomp_%lu_%d", name_directory_decomp, lim0, r0);
	    FILE *file = fopen(name_file, "r");

	    tab_decomp = tabular_decomp_fscan(file);

	    if (tab_decomp == NULL) {
		fprintf(stderr, "impossible to read '%s'\n", name_file);
		exit(EXIT_FAILURE);
	    }
	    fclose(file);
	}
	data_rat[r0] = generate_strategy_cado(methods, tab_decomp,
					      fbb0, lpb0, r0);
	tabular_decomp_free(tab_decomp);
    }

    /*
       read good elements for r_2 in the array data_r1 and compute the
       data for each r_1. So :
     */
    lim = 2 * fbb1 - 1;
    for (int r1 = 0; r1 <= mfb1; r1++) {
	tabular_decomp_t *tab_decomp = NULL;
	if (r1 >= lim) {
	    char name_file[200];
	    snprintf(name_file, sizeof(name_file),
		    "%s/decomp_%lu_%d", name_directory_decomp, lim1, r1);
	    FILE *file = fopen(name_file, "r");

	    tab_decomp = tabular_decomp_fscan(file);

	    if (tab_decomp == NULL) {
		fprintf(stderr, "impossible to read '%s'\n", name_file);
		exit(EXIT_FAILURE);
	    }
	    fclose(file);
	}

	tabular_strategy_t *strat_r1 =
	    generate_strategy_cado(methods, tab_decomp, fbb1, lpb1, r1);

	tabular_decomp_free(tab_decomp);

	for (int r0 = 0; r0 <= mfb0; r0++) {
	  tabular_strategy_t *res =
	    generate_strategy_r0_r1(data_rat[r0], strat_r1);
	  matrix[r0][r1] = res;
	}
	tabular_strategy_free(strat_r1);
    }

    //free
    for (int r0 = 0; r0 <= mfb0; r0++)
	tabular_strategy_free(data_rat[r0]);
    free(data_rat);

    fm_free(zero);
    return matrix;
}


/************************************************************************/
/*                  To interleave our strategies                        */
/************************************************************************/
static MAYBE_UNUSED int is_good_decomp(decomp_t * dec, int len_p_min, int len_p_max)
{
    int len = dec->len;
    for (int i = 0; i < len; i++)
	if (dec->tab[i] > len_p_max || dec->tab[i] < len_p_min)
	    return false;
    return true;
}

//problem if tab_edc == NULL;
static double compute_time_strategy_ileav(tabular_decomp_t ** init_tab, strategy_t * strat,
				   int* fbb, int* lpb, int* r)
{
    int nb_fm = tabular_fm_get_index(strat->tab_fm);
    tabular_fm_t *tab_fm = strat->tab_fm;

    //{{
    int last_method_side[2] = {0,0};
    for (int index_fm = 0; index_fm < nb_fm; index_fm++) {
      int side = strat->side[index_fm];
      last_method_side[side] = index_fm;
    }
    int end_of_one_side = (last_method_side[0] < last_method_side[1])?
      last_method_side[0]: last_method_side[1];
    
    //}}
    
    //{{
    /*
      We add 0.5 to the length of one word, because for our times the
      lenght is inclusive. For example, if MODREDCUL_MAXBITS = 64
      bits, a cofactor is in one word if is lenght is less OR equal to
      64 bits. So, if you don't add 0.5 to MODREDCUL_MAXBITS, you
      lost the equal and thus insert an error in your maths.
    */
    int ind_time[2];
    for (int side = 0; side < 2; side++)
      {
	double half_word = (MODREDCUL_MAXBITS+0.5)/2.0;
	int number_half_wd = floor(r[side] /half_word);
	//the next computation is necessary in the relation with the
	//benchmark in gfm!
	ind_time[side] = (number_half_wd <2)? 0: number_half_wd - 1;
      }
    //}}
    //double prob = 0;
    double time_average = 0;
    //store the number of elements in the different decompositions!
    double all = 0.0;
    int nb_decomp0 = init_tab[0]->index;
    int nb_decomp1 = init_tab[1]->index;
    for (int index_decomp0 = 0; index_decomp0 < nb_decomp0; index_decomp0++)
      for (int index_decomp1 = 0; index_decomp1 < nb_decomp1; index_decomp1++)
	{
	  decomp_t *dec[2];
	  dec[0] = init_tab[0]->tab[index_decomp0];
	  dec[1] = init_tab[1]->tab[index_decomp1];
	  double time_dec = 0;
	  double proba_fail_side[2] = {1, 1};//didn't find a non-trivial factor!

	  double proba_run_next_fm = 1;

	  //to know if a side should be consider or not!
	  int is_bad_dec[2] = {!is_good_decomp (dec[0], fbb[0], lpb[0]),
			       !is_good_decomp (dec[1], fbb[1], lpb[1])};
	  
	  double time_method = 0;

	  //compute the time of each decomposition
	  for (int index_fm = 0; index_fm < nb_fm; index_fm++) {	    
	    fm_t* elem = tabular_fm_get_fm(tab_fm, index_fm);
	    int side = strat->side[index_fm];

	    //probability to run the next method!
	    /*the next method is used if:
	      - the previous method of the same side did not find a factor!
	      - you didn't have a failure in the other side 
	      (i.e.: find a bad decomposition, or find anything 
	      and it was our last method for this side!).
	    */

	    if (is_bad_dec[(1+side)%2] == 1)
	      proba_run_next_fm = proba_fail_side[side] *
		(proba_fail_side[(1+side)%2]);
	    //the side will continue if the other side is smooth
	    //(proba=1-proba echec)!
	    else if (index_fm > last_method_side[(1+side)%2]
		     && is_bad_dec[(1+side)%2] == 0)
	      proba_run_next_fm = proba_fail_side[side] *
		(1-proba_fail_side[(1+side)%2]);
	    
	    else 
	      proba_run_next_fm = proba_fail_side[side];
	    
	    int len_time = fm_get_len_time (elem);
	    if (ind_time[side] >= len_time)
	      time_method = elem->time[len_time-1];
	    else
	      time_method = elem->time[ind_time[side]];
	    
	    time_dec += time_method * proba_run_next_fm;

	    double proba_fail_method =
	      compute_proba_method_one_decomp (dec[side], elem);
	    
	    proba_fail_side[side] *= proba_fail_method;

	    
	    if (elem->method[0] == PM1_METHOD ||
		elem->method[0] == PP1_27_METHOD ||
		elem->method[0] == PP1_65_METHOD)
	      {
		//because if you chain PP1||PM1 to PM1||PP1-->they are
		//not independant.
		proba_fail_side[side] = (proba_fail_side[side] + proba_fail_method) / 2;
	      }
	    if (index_fm == end_of_one_side)
	      {
		//this side is a bad decomposition in all cases!
		if (is_bad_dec[side] == 1)
		  proba_fail_side[side] = 0; 
	      }
	  }	  
	  double nb_elem = dec[1]->nb_elem*dec[0]->nb_elem;
	  time_average += time_dec * nb_elem;
	  /* int is_good[2] = {is_good_decomp (dec[0], fbb[0], lpb[0]), */
	  /* 		    is_good_decomp (dec[1], fbb[1], lpb[1])};	   */
	  /* prob += nb_elem * (1-proba_fail_side[0]) * is_good[0]* */
	  /*   (1-proba_fail_side[1]) * is_good[1]; */
	  
	  all += nb_elem;
	}
    //printf ("prob = %lf\n", prob/all);
    if (all < 0.00000001) //all==0
      {
	return 0;
      }
    return time_average / all;
}


static strategy_t*
gen_strat_r0_r1_ileav_st_rec (strategy_t * strat_r0, int index_r0,
			      strategy_t * strat_r1, int index_r1,
			      tabular_decomp_t ** init_tab,
			      int* fbb, int* lpb, int* r,
			      strategy_t* current_st, int current_index)
{
  int len_r0 = strat_r0->tab_fm->index;
  int len_r1 = strat_r1->tab_fm->index;
  int max_len = len_r0 + len_r1;

  /* if (index_r0 >= len_r0 && */
  /*     index_r1 >= len_r1) */
  /*   return NULL; */
  /* printf ("%d (max= %d), %d (max=%d)\n", index_r0, len_r0, */
  /* 	  index_r1, len_r1); */
  
  if (current_st == NULL)
    {
      current_st = strategy_create ();
      current_index = 0;
      current_st->len_side = max_len;
      current_st->side = (int*) malloc(sizeof(int) * (current_st->len_side));
    }
    
  strategy_t *tmp0 = NULL, *tmp1 = NULL;

  if (index_r0 < len_r0)
    {
      if (index_r0 < 5)
	{
	  tabular_fm_set_fm_index (current_st->tab_fm, strat_r0->tab_fm->tab[index_r0],
				   current_index);
	  current_st->side[current_index] = 0;
	  tmp0 = gen_strat_r0_r1_ileav_st_rec (strat_r0, index_r0 + 1,
					       strat_r1, index_r1,
					       init_tab, fbb, lpb, r,
					       current_st, current_index+1);
	}
      else //concat the stay of our method!!!!
	{
	  int len = len_r0 - index_r0;
	  for (int i = 0; i < len; i++)
	    {
	      tabular_fm_set_fm_index (current_st->tab_fm,
				       strat_r0->tab_fm->tab[index_r0 + i],
				       current_index + i);
	      current_st->side[current_index + i] = 0;
	    }
	  tmp0 = gen_strat_r0_r1_ileav_st_rec (strat_r0, len_r0,
					       strat_r1, index_r1,
					       init_tab, fbb, lpb, r,
					       current_st, current_index + len);

	}
      //strategy_print (tmp0);
    }

  if (index_r1 < len_r1)
    {
      if (index_r1 < 5)
	{
	  tabular_fm_set_fm_index (current_st->tab_fm, strat_r1->tab_fm->tab[index_r1],
				   current_index);
	  current_st->side[current_index] = 1;
	  tmp1 = gen_strat_r0_r1_ileav_st_rec (strat_r0, index_r0,
					       strat_r1, index_r1 + 1,
					       init_tab, fbb, lpb, r,
					       current_st, current_index+1);
	}
      else //concat the stay of our method!!!!
	{
	  int len = len_r1 - index_r1;
	  for (int i = 0; i < len; i++)
	    {
	      tabular_fm_set_fm_index (current_st->tab_fm,
				       strat_r1->tab_fm->tab[index_r1 + i],
				       current_index + i);
	      current_st->side[current_index + i] = 1;
	    }
	  tmp1 = gen_strat_r0_r1_ileav_st_rec (strat_r0, index_r0,
					       strat_r1, len_r1,
					       init_tab, fbb, lpb, r,
					       current_st, current_index + len);

	}

      //strategy_print (tmp1);
    }
  
  if (current_index == 0) //it's the first round of our recursion!
    strategy_free (current_st);

  //end of the recursion!
  
  if (tmp0 == NULL && tmp1 == NULL)
    {
      strategy_t* final_st = strategy_create ();
      final_st->len_side = max_len;
      final_st->side = (int*) malloc(sizeof(int) * (final_st->len_side));
      //copy the current strategy!
      for (int i = 0; i < current_index; i++)
	{
	  strategy_add_fm(final_st, current_st->tab_fm->tab[i]);
	  final_st->side[i] = current_st->side[i];
	}
      double prob = strat_r0->proba * strat_r1->proba;
      strategy_set_proba (final_st, prob);
      
      double time = compute_time_strategy_ileav (init_tab, final_st, fbb, lpb, r);
      strategy_set_time(final_st, time);
      /* printf ("final_st\n"); */
      /* strategy_print (final_st); */
      /* getchar (); */
      /* printf ("temps = %lf, index_current = %d\n", time, current_index); */
      return final_st;
    }
  
  if (tmp0 == NULL)
    return tmp1;
  else if (tmp1 == NULL)
    return tmp0;
  else if (tmp0->time < tmp1->time)
    {
      strategy_free (tmp1);
      return tmp0;
    }
  else //(tmp0->time > tmp1->time)
    {
      strategy_free (tmp0);
      return tmp1;
    }
}
			      
static MAYBE_UNUSED strategy_t*
gen_strat_r0_r1_ileav_st(strategy_t * strat_r0,
			 strategy_t * strat_r1,
			 tabular_decomp_t ** init_tab,
			 int* fbb, int* lpb, int* r)
{
  int len_r0 = strat_r0->tab_fm->index;
  int len_r1 = strat_r1->tab_fm->index;
  int max_len = len_r0 + len_r1;

  strategy_t* st = strategy_create ();
  st->len_side = max_len;
  st->side = (int*) malloc(sizeof(int) * (st->len_side));
  int index_r0 = 0, index_r1 = 0;
  int i = 0;
  int sequence = 2;
  while (i < max_len)
    {
      int test = sequence%2;
      if (index_r0 < len_r0 && (test==0 || index_r1 >= len_r1))
	{
	  strategy_add_fm(st, strat_r0->tab_fm->tab[index_r0++]);
	  st->side[i++] = 0;
	}
      else if (index_r1 < len_r1)
	{
	  strategy_add_fm(st, strat_r1->tab_fm->tab[index_r1++]);
	  st->side[i++] = 1;
	}
      else if (index_r1 == len_r1 && index_r0 == len_r0)
	break;
      sequence = (sequence -test)/2;
    }
  double prob = strat_r0->proba * strat_r1->proba;
  strategy_set_proba (st, prob);
      
  double time = compute_time_strategy_ileav (init_tab, st, fbb, lpb, r);
  strategy_set_time(st, time);

  return st;
}

static tabular_strategy_t *gen_strat_r0_r1_ileav(tabular_strategy_t * strat_r0,
					  tabular_strategy_t * strat_r1,
					  tabular_decomp_t ** init_tab,
					  int* fbb, int* lpb, int* r)
{
  tabular_strategy_t* res = tabular_strategy_create();//generate_strategy_r0_r1(strat_r0, strat_r1);
  int len0 = strat_r0->index;
  int len1 = strat_r1->index;
  //printf ("r0 = %d, r1=%d\n", r[0], r[1]);
  //printf ("CLASSIC proba=%lf, time=%lf\n", res->tab[0]->proba, res->tab[0]->time);
  //  strategy_print (res->tab[0]);
  for (int r0 = 0; r0 < len0; r0++)
    for (int r1 = 0; r1 < len1; r1++)
      {
	/* strategy_t* st1 = gen_strat_r0_r1_ileav_st(strat_r0->tab[r0], */
	/* 					   strat_r1->tab[r1], */
	/* 					   init_tab, fbb, lpb, r); */
	strategy_t* st2 = gen_strat_r0_r1_ileav_st_rec(strat_r0->tab[r0], 0,
						       strat_r1->tab[r1], 0,
						       init_tab, fbb, lpb, r,
						       NULL, 0);
	//printf ("INTERL proba=%lf, time=%lf\n", st->proba, st->time);
	//strategy_print (st);
	tabular_strategy_add_strategy (res, st2);
	//Test entrelacement:
	/* int current_side = st2->side[0]; */
	/* int nb_chg = 0; */
	/* for (int i = 1; i < st2->len_side; i++) */
	/*   { */
	/*     if (st2->side[i] != current_side) */
	/*       nb_chg++; */
	/*     current_side = st2->side[i]; */
	/*   } */
	/* if (nb_chg > 1) */
	/*   { */
	/*     printf ("r0 = %d, r1=%d\n", r[0], r[1]); */
	/* printf ("CLASSI proba=%lf, time=%lf\n", st1->proba, st1->time); */
	/* printf ("INTERL proba=%lf, time=%lf\n", st2->proba, st2->time); */
	    //getchar();
	//}
	//compare probabilities: 
	//strategy_free (st1);
	strategy_free (st2);
      }
  return res;
}

#ifdef COMPILE_DEAD_CODE
/*
  Test an interleaving with cado!!
 */
static MAYBE_UNUSED tabular_strategy_t *** generate_matrix_cado_ileav(
        const char *name_directory_decomp,
        tabular_fm_t * methods,
        unsigned long lim0, int lpb0, int mfb0,
        unsigned long lim1, int lpb1, int mfb1)
{
    /*
       allocates the matrix which contains all optimal strategies for
       each couple (r0, r1): r0 is the lenght of the cofactor in
       the first side, and r1 for the second side.
     */
    tabular_strategy_t ***matrix = (tabular_strategy_t***) malloc(sizeof(*matrix) * (mfb0 + 1));
    ASSERT(matrix != NULL);

    for (int r0 = 0; r0 <= mfb0; r0++) {
	matrix[r0] = (tabular_strategy_t**) malloc(sizeof(tabular_strategy_t *) * (mfb1 + 1));
	ASSERT(matrix[r0] != NULL);
    }
    
    int fbb0 = ceil (log2 ((double) (lim0 + 1)));
    int fbb1 = ceil (log2 ((double) (lim1 + 1)));

    /*
       For each r0 in [fbb0..mfb0], we precompute and strore the data
       for all strategies.  Whereas, for each r1, the values will be
       computed sequentially. 
     */
    fm_t *zero = fm_create();
    unsigned long method_zero[4] = { 0, 0, 0, 0 };
    fm_set_method(zero, method_zero, 4);

    tabular_strategy_t **data_rat = (tabular_strategy_t**) malloc(sizeof(tabular_strategy_t *) * (mfb0 + 1));
    ASSERT(data_rat);

    int lim = 2 * fbb0 - 1;
    for (int r0 = 0; r0 <= mfb0; r0++) {
	tabular_decomp_t *tab_decomp = NULL;
	if (r0 >= lim) {
	    char name_file[200];
	    snprintf(name_file, sizeof(name_file),
		    "%s/decomp_%lu_%d", name_directory_decomp, lim0, r0);
	    FILE *file = fopen(name_file, "r");

	    tab_decomp = tabular_decomp_fscan(file);

	    if (tab_decomp == NULL) {
		fprintf(stderr, "impossible to read '%s'\n", name_file);
		exit(EXIT_FAILURE);
	    }
	    fclose(file);
	}
	data_rat[r0] = generate_strategy_cado(methods, tab_decomp,
					      fbb0, lpb0, r0);
	tabular_decomp_free(tab_decomp);
    }

    /*
       read good elements for r_2 in the array data_r1 and compute the
       data for each r_1. So :
     */
    lim = 2 * fbb1 - 1;
    for (int r1 = 0; r1 <= mfb1; r1++) {
      printf ("r1 = %d\n", r1);
	tabular_decomp_t *tab_decomp = NULL;
	char name_file[200];
	if (r1 >= lim) {
	    snprintf(name_file, sizeof(name_file),
		    "%s/decomp_%lu_%d", name_directory_decomp, lim1, r1);
	    FILE *file = fopen(name_file, "r");

	    tab_decomp = tabular_decomp_fscan(file);

	    if (tab_decomp == NULL) {
		fprintf(stderr, "impossible to read '%s'\n", name_file);
		exit(EXIT_FAILURE);
	    }
	    fclose(file);
	}

	tabular_strategy_t *strat_r1 =
	    generate_strategy_cado(methods, tab_decomp, fbb1, lpb1, r1);


	for (int r0 = 0; r0 <= mfb0; r0++) {
	  if (r0 < 2*fbb0-1 || r1 < 2*fbb1-1)
	    {
	      matrix[r0][r1] =
		generate_strategy_r0_r1(data_rat[r0], strat_r1);
	    }
	  else
	    {
	      snprintf(name_file, sizeof(name_file),
		      "%s/decomp_%lu_%d", name_directory_decomp, lim0, r0);
	      FILE *file = fopen(name_file, "r");

	      tabular_decomp_t* tab_decomp_r0 = tabular_decomp_fscan(file);
	      
	      if (tab_decomp_r0 == NULL) {
		fprintf(stderr, "impossible to read '%s'\n", name_file);
		exit(EXIT_FAILURE);
	      }
	      fclose(file);

	      
	      tabular_decomp_t* init_tab[2] = {tab_decomp_r0, tab_decomp}; //todo!!
	      int fbb[2] = {fbb0, fbb1};
	      int lpb[2] = {lpb0, lpb1};
	      int r[2] = {r0, r1};
	      matrix[r0][r1] = gen_strat_r0_r1_ileav(data_rat[r0],
						     strat_r1, init_tab,
						     fbb, lpb, r);
	    }
	}
	tabular_strategy_free(strat_r1);
	tabular_decomp_free(tab_decomp);
    }

    //free
    for (int r0 = 0; r0 <= mfb0; r0++)
	tabular_strategy_free(data_rat[r0]);
    free(data_rat);

    fm_free(zero);
    return matrix;
}
#endif  /* COMPILE_DEAD_CODE */

static tabular_strategy_t ***generate_matrix_ileav(const char *name_directory_decomp,
					    const char *name_directory_str,
					    unsigned long lim0, int lpb0,
					    int mfb0,
					    unsigned long lim1, int lpb1,
					    int mfb1)
{
    /*
       allocates the matrix which contains all optimal strategies for
       each couple (r0, r1): r0 is the lenght of the cofactor in
       the first side, and r1 for the second side.
     */
    tabular_strategy_t ***matrix = (tabular_strategy_t***) malloc(sizeof(*matrix) * (mfb0 + 1));
    ASSERT(matrix != NULL);

    for (int r0 = 0; r0 <= mfb0; r0++) {
	matrix[r0] = (tabular_strategy_t**) malloc(sizeof(tabular_strategy_t *) * (mfb1 + 1));
	ASSERT(matrix[r0] != NULL);
    }
    
    int fbb0 = ceil (log2 ((double) (lim0 + 1)));
    int fbb1 = ceil (log2 ((double) (lim1 + 1)));

    /*
       For each r0 in [fbb0..mfb0], we precompute and strore the data
       for all strategies.  Whereas, for each r1, the values will be
       computed sequentially. 
     */
    fm_t *zero = fm_create();
    unsigned long method_zero[4] = { 0, 0, 0, 0 };
    fm_set_method(zero, method_zero, 4);

    tabular_strategy_t **data_rat = (tabular_strategy_t**) malloc(sizeof(tabular_strategy_t *) * (mfb0 + 1));
    ASSERT(data_rat);

    int lim = 2 * fbb0 - 1;
    for (int r0 = 0; r0 <= mfb0; r0++) {
	char name_file_in[strlen(name_directory_str) + 64];
	FILE * file_in;
	//get back the best strategies for r0!
	snprintf(name_file_in, sizeof(name_file_in), "%s/strategies%lu_%d",
		name_directory_str, lim0, r0);
	file_in = fopen(name_file_in, "r");
	data_rat[r0] = tabular_strategy_fscan(file_in);
	if (data_rat[r0] == NULL) {
	  fprintf(stderr,
		  "Parser error: can't read the file '%s'\n",
		  name_file_in);
	  exit(EXIT_FAILURE);
	}
	fclose(file_in);
    }

    /*
       read good elements for r_2 in the array data_r1 and compute the
       data for each r_1. So :
     */
    lim = 2 * fbb1 - 1;
    for (int r1 = 0; r1 <= mfb1; r1++) {
      printf ("r1 = %d\n", r1);
	tabular_decomp_t *tab_decomp = NULL;
	char name_file[200];
	if (r1 >= lim) {
            snprintf(name_file, sizeof(name_file),
		    "%s/decomp_%lu_%d", name_directory_decomp, lim1, r1);
	    FILE *file = fopen(name_file, "r");

	    tab_decomp = tabular_decomp_fscan(file);

	    if (tab_decomp == NULL) {
		fprintf(stderr, "impossible to read '%s'\n", name_file);
		exit(EXIT_FAILURE);
	    }
	    fclose(file);
	}
	char name_file_in[strlen(name_directory_str) + 64];
	FILE * file_in;
	//get back the best strategies for r0!
	snprintf(name_file_in, sizeof(name_file_in), 
                "%s/strategies%lu_%d",
		name_directory_str, lim1, r1);
	file_in = fopen(name_file_in, "r");
	tabular_strategy_t *strat_r1 = tabular_strategy_fscan(file_in);
	if (strat_r1 == NULL) {
	  fprintf(stderr,
		  "Parser error: can't read the file '%s'\n",
		  name_file_in);
	  exit(EXIT_FAILURE);
	}
	fclose(file_in);

	for (int r0 = 0; r0 <= mfb0; r0++) {
	  if (r0 < 2*fbb0-1 || r0 > 3*fbb0-2 ||
	      r1 < 2*fbb1-1 || r1 > 3*fbb1-2)
	    {
	      matrix[r0][r1] =
		generate_strategy_r0_r1(data_rat[r0], strat_r1);
	    }
	  else
	    {
                snprintf(name_file, sizeof(name_file),
		      "%s/decomp_%lu_%d", name_directory_decomp, lim0, r0);
	      FILE *file = fopen(name_file, "r");

	      tabular_decomp_t* tab_decomp_r0 = tabular_decomp_fscan(file);
	      
	      if (tab_decomp_r0 == NULL) {
		fprintf(stderr, "impossible to read '%s'\n", name_file);
		exit(EXIT_FAILURE);
	      }
	      fclose(file);

	      
	      tabular_decomp_t* init_tab[2] = {tab_decomp_r0, tab_decomp}; //todo!!
	      int fbb[2] = {fbb0, fbb1};
	      int lpb[2] = {lpb0, lpb1};
	      int r[2] = {r0, r1};
	      matrix[r0][r1] =
		gen_strat_r0_r1_ileav(data_rat[r0], strat_r1,
					      init_tab, fbb, lpb, r);
	    }
	}
	tabular_strategy_free(strat_r1);
	tabular_decomp_free(tab_decomp);
    }

    //free
    for (int r0 = 0; r0 <= mfb0; r0++)
	tabular_strategy_free(data_rat[r0]);
    free(data_rat);

    fm_free(zero);
    return matrix;
}



/************************************************************************/
/*                  To bench our strategies                             */
/************************************************************************/

#ifdef COMPILE_DEAD_CODE
static int
select_random_index_dec(double sum_nb_elem, tabular_decomp_t* t, gmp_randstate_ptr state)
{
    //100000 to consider approximation of distribution
    double alea = gmp_urandomm_ui(state, 10000);
    int i = 0;
    double bound = (t->tab[0]->nb_elem/sum_nb_elem) * 10000;
    int len = t->index;
    while (i < (len-1) && (alea - bound) >= 1) {
	i++;
	bound += (t->tab[i]->nb_elem/sum_nb_elem) * 10000;
    }
    return i;
}
#endif  /* COMPILE_DEAD_CODE */


/*
  This function compute directly the probabity and the time to find a
  factor in a good decomposition in a cofactor of length r.

TODO: a function of the same name, yet behaving somewhat differently, exists in tests/sieve/strategies/test_generate_strategies.cpp
 */
#ifdef COMPILE_DEAD_CODE
static MAYBE_UNUSED weighted_success
bench_proba_time_st(gmp_randstate_t state, facul_strategy_oneside const & strategy,
		    tabular_decomp_t* init_tab, int r, MAYBE_UNUSED int lpb)
{
    int nb_test = 0, nb_success = 0;
    int nb_test_max = 100000;
    double time = 0;
    
    double sum_dec = 0;
    for (int i = 0; i < init_tab->index; i++)
	sum_dec += init_tab->tab[i]->nb_elem;

    std::vector<cxx_mpz> f;

    while (nb_test < nb_test_max)
	{
	    /*
	       f will contain the prime factor of N that the strategy
	       found.  Note that N is composed by two prime factors by the
	       previous function.
	    */
            int index = select_random_index_dec(sum_dec, init_tab, state);
            int len_p = init_tab->tab[index]->tab[0];

            cxx_mpz N = generate_composite_integer(state, len_p, r);

            f.clear();
            time -= microseconds();
	    int nfound = facul(f, N, strategy);
            time += microseconds();
            nb_success += (nfound != 0);
	    nb_test++;
	}
    
    return weighted_success(nb_success, time, nb_test);
}
#endif  /* COMPILE_DEAD_CODE */

#ifdef COMPILE_DEAD_CODE
//convert type: from strategy_t to facul_strategy_t.
static facul_strategy_oneside convert_strategy_to_facul_strategy (strategy_t* t, unsigned long lim,
						      int lpb, int side)
{
    tabular_fm_t* tab_fm = strategy_get_tab_fm (t);

    std::vector<facul_method::parameters> mps;

    for (int i = 0; i < tab_fm->index; i++) {
	  if (t->side[i] != side)
              continue;

          fm_t* fm = tab_fm->tab[i];
          int method = (int)fm->method[0];
          ec_parameterization_t curve = (ec_parameterization_t) fm->method[1];
          unsigned long B1 = fm->method[2];
          unsigned long B2 = fm->method[3];

          mps.emplace_back(
                  method, B1, B2,
                  curve,
                  curve == MONTY16 ? 1 : 4,
                  0     // extra_primes. It's 1 almost everywhere else, wtf?
            );
    }

    return facul_strategy_oneside(lim, lpb, 4 * lpb, mps, 0);
}
#endif  /* COMPILE_DEAD_CODE */

#ifdef COMPILE_DEAD_CODE
facul_strategies convert_strategy_to_facul_strategies (strategy_t* t,
        const unsigned int * r,
        const unsigned long * fbb,
        const unsigned int * lpb, 
        const unsigned int * mfb,
        gmp_randstate_ptr rstate)
{
    /* This is dead code, and it probably doesn't make sense to think
     * about strategies for >2 sides anyway.
     */
    std::vector<unsigned long> B(2);
    std::vector<unsigned int> lpb_(2);
    std::vector<unsigned int> mfb_(2);
    std::array<unsigned int, 2> r_;
    auto pB = B.begin();
    auto plpb = lpb_.begin();
    auto pmfb = mfb_.begin();
    auto pr = r_.begin();
    for(int side = 0 ; side < 2 ; side++) {
        *pB++ = *fbb++;
        *plpb++ = *lpb++;
        *pmfb++ = *mfb++;
        *pr++ = *r++;
    }

    std::vector<facul_method::parameters_with_side> mps;

    tabular_fm_t* tab_fm = strategy_get_tab_fm (t);
    for (int i = 0; i < tab_fm->index; i++)
    {
        fm_t* fm = tab_fm->tab[i];
        int method = (int)fm->method[0];
        ec_parameterization_t curve = (ec_parameterization_t) fm->method[1];
        unsigned long B1 = fm->method[2];
        unsigned long B2 = fm->method[3];
        int side = (t->side != NULL)?t->side[i]:0;
        unsigned long parameter = 1;
        if (method == EC_METHOD && curve != MONTY16) {
            for( ; (parameter = u64_random(rstate)) < 2 ; );
        }
        mps.emplace_back(side, method, B1, B2, curve, parameter, 1);
    }

    facul_strategies::strategy_file S;
    S[r_] = mps;

    return facul_strategies(B, lpb_, mfb_, true, S, 0);
}
#endif  /* COMPILE_DEAD_CODE */

#ifdef COMPILE_DEAD_CODE
static MAYBE_UNUSED weighted_success
bench_proba_time_st_both(gmp_randstate_t state,
                         strategy_t * t,
			 tabular_decomp_t ** init_tab,
			 const unsigned int * r,
                         const unsigned long * fbb,
                         const unsigned int * lpb,
                         const unsigned int * mfb)
{
    /* This is dead code, and it probably doesn't make sense to think
     * about strategies for >2 sides anyway.
     */
    unsigned long lim[2] = {1UL << (fbb[0]-1), 1UL << (fbb[1]-1) };

    int nb_test = 0, nb_success = 0;
    int nb_test_max = 10000;
    double time = 0;
    
    gmp_randstate_t state_copy;
    gmp_randinit_set(state_copy, state);
  
    double sum_dec[2] = {0,0};
    for (int side = 0; side < 2; side++)
      for (int i = 0; i < init_tab[side]->index; i++)
	sum_dec[side] += init_tab[side]->tab[i]->nb_elem;

#if 1
    // Classic: bench without interleaving

    facul_strategy_oneside facul_st_s0 = convert_strategy_to_facul_strategy (t,lim[0], lpb[0], 0);
    facul_strategy_oneside facul_st_s1 = convert_strategy_to_facul_strategy (t,lim[1], lpb[1], 1);
    printf ("classic\n");
    nb_success = 0;
    nb_test = 0;
    time = 0;
    {
        std::vector<cxx_mpz> f;
        while (nb_test < nb_test_max)
        {
            cxx_mpz N[2];
            for(int side = 0 ; side < 2 ; side++) {
                int index = select_random_index_dec(sum_dec[side], init_tab[side], state);
                int len_p = init_tab[side]->tab[index]->tab[1];
                N[side] = generate_composite_integer(state, len_p, r[side]);
            }
            /*
               f will contain the prime factor of N that the strategy
               found.  Note that N is composed by two prime factors by the
               previous function.
               */
            f.clear();
            time -= microseconds();
            int nfound = facul(f, N[0], facul_st_s0);
            if (nfound) {
                f.clear();
                nfound = facul(f, N[1], facul_st_s1);
                nb_success += (nfound != 0);
            }
            time += microseconds();
            nb_test++;
            //getchar ();
        }
    }
    weighted_success res(nb_success, time, nb_test);
    printf ("classic: prob = %lf, temps = %lf\n", res.prob, res.time);
#endif

#if 1
    printf ("interl\n");
    facul_strategies facul_st = convert_strategy_to_facul_strategies (t,r,lim,lpb, mfb, state);

    nb_success = 0;
    nb_test = 0;
    time = 0;

    {
        std::vector<std::vector<cxx_mpz>> f(2);
        while (nb_test < nb_test_max)
        {
            std::vector<cxx_mpz> N(2);
            for(int side = 0 ; side < 2 ; side++) {
                int index = select_random_index_dec(sum_dec[side], init_tab[side], state);
                int len_p = init_tab[side]->tab[index]->tab[1];
                N[side] = generate_composite_integer(state_copy, len_p, r[side]);
            }
            /*
               f will contain the prime factor of N that the strategy
               found.  Note that N is composed by two prime factors by the
               previous function.
               */
            f[0].clear();
            f[1].clear();
            time -= microseconds();
            int is_smooth[2] = {0,0};
            facul_both(f, N, facul_st, is_smooth);
            if (is_smooth[0]==1 && is_smooth[1]==1)
                nb_success++;
            time += microseconds();
            nb_test++;
            /* printf ("[%d, %d]\n", (int)mpz_sizeinbase (N[nb_test][0], 2), */
            /* 	  (int)mpz_sizeinbase (N[nb_test][1], 2)); */
            /* printf ("nb_success = %d\n", nb_success); */
            //getchar();
        }
    }
    weighted_success res2(nb_success, time, nb_test);
    printf ("interL: proba = %lf, time = %lf\n", res2.prob, res2.time);
#endif

    gmp_randclear(state_copy);

    return res;
}
#endif  /* COMPILE_DEAD_CODE */


/************************************************************************/
/*                            USAGE                                     */
/************************************************************************/

static void declare_usage(param_list pl)
{
    param_list_usage_header(pl,
			    "This binary allows to test the strategy of cado,"
			    "and especially compute the theorical number of "
			    "relations found per second by this strategy.\n");

    param_list_decl_usage(pl, "lim0",
			  "set rationnal factor base bound to lim0\n");
    param_list_decl_usage(pl, "lim1",
			  "set algebraic factor base bound to lim1\n");
    param_list_decl_usage(pl, "lpb0",
			  "set rational large prime bound to 2^lpb0");
    param_list_decl_usage(pl, "lpb1",
			  "set algebraic large prime bound to 2^lpb1");
    param_list_decl_usage(pl, "mfb0", "set the first cofactor bound to 2^mfb0");
    param_list_decl_usage(pl, "mfb1",
			  "set the second cofactor bound to 2^mfb1");
    param_list_decl_usage(pl, "decomp",
			  "to locate the file or the directory , according to\n"
			  "\t \t if you need one or several files,\n"
			  "\t \t which contain(s) the file(s) of cofactors decompositions.");
    param_list_decl_usage(pl, "dist",
			  "the pathname of our file which contains the distribution\n"
			  "\t\t of our pairs of cofactors.");
    param_list_decl_usage(pl, "t",
			  "specify the time (seconds) to optain cofactors in the file\n"
			  "\t\t given by the option 'dist'.");
    param_list_decl_usage(pl, "out",
			  "the output file which contain our strategies\n");
}


/************************************************************************/
/*     MAIN                                                             */
/************************************************************************/

// coverity[root_function]
int main(int argc, char *argv[])
{
    param_list pl;
    param_list_init(pl);
    declare_usage(pl);

    if (argc <= 1) {
	param_list_print_usage(pl, argv[0], stderr);
	exit(EXIT_FAILURE);
    }

    argv++, argc--;
    for (; argc;) {
	if (param_list_update_cmdline(pl, &argc, &argv)) {
	    continue;
	}
	/* Could also be a file */
	FILE *f;
	if ((f = fopen(argv[0], "r")) != NULL) {
	    param_list_read_stream(pl, f, 0);
	    fclose(f);
	    argv++, argc--;
	    continue;
	}
	fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
	param_list_print_usage(pl, argv[0], stderr);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }

    //option parser
    unsigned long lim0 = 0;
    int lpb0 = -1;
    unsigned long lim1 = 0;
    int lpb1 = -1;
    int mfb0 = -1;
    int mfb1 = -1;
    double C0 = -1;

    param_list_parse_ulong(pl, "lim0", &lim0);
    param_list_parse_int(pl, "lpb1", &lpb1);
    param_list_parse_ulong(pl, "lim1", &lim1);
    param_list_parse_int(pl, "mfb1", &mfb1);
    param_list_parse_int(pl, "lpb0", &lpb0);
    param_list_parse_int(pl, "mfb0", &mfb0);
    param_list_parse_double(pl, "t", &C0);

    if (lim0 == 0 || lpb0 == -1 || mfb0 == -1 ||
	lim1 == 0 || lpb1 == -1 || mfb1 == -1 || C0 == -1) {
	fputs("ALL parameters are mandatory!\n", stderr);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }

    //{{just to obtain our matrix
    /* tabular_fm_t* methods = generate_methods_cado(lpb0); */
    /* tabular_strategy_t* tab = tabular_strategy_create (); */
    /* strategy_t* st = strategy_create (); */
    /* int len = 3+nb_curves(lpb0); */
    /* for (int i = 0; i < len; i++) */
    /*   strategy_add_fm(st, methods->tab[i]); */
    /* tabular_strategy_add_strategy (tab, st); */
    /* tabular_strategy_t* res=  generate_strategy_r0_r1 (tab, tab); */
    /* FILE* filee = fopen("strategy_check_28_33_66_99","w"); */
    /* for (int r0 = 0; r0 <= mfb0; r0++) */
    /*   for (int r1 = 0; r1 <= mfb1; r1++) */
    /* 	{ */
    /* 	  fprintf(filee, */
    /* 		  "[r0=%d, r1=%d] : (p = %lf, t = %lf)\n", */
    /* 		  r0, r1, 0.0,0.0); */
    /* 	  strategy_fprint_design(filee, res->tab[0]); */
    /* 	} */
    /* tabular_strategy_free (res); */
    /* tabular_strategy_free (tab); */
    /* strategy_free (st); */
    /* tabular_fm_free (methods); */
    /* fclose (filee); */
    /* exit(1); */
    /* //}} */

    //int fbb0 = ceil (log2 ((double) (lim0+1)));
    //int fbb1 = ceil (log2 ((double) (lim1+1)));
    //int fbb = (fbb0 < fbb1) ? fbb0 : fbb1;
    //    int lpb = (lpb0 > lpb1) ? lpb0 : lpb1;
    //convert the time in micro-s. because all previous binaries
    //compute their times in micro-s.
    C0 *= 1000000;		//s-->micro-s 

    //option: tab decomp
    const char *name_directory_decomp;
    //  "/localdisk/trichard/results/decomp_cofactor/decomp_tmp";
    if ((name_directory_decomp =
	 param_list_lookup_string(pl, "decomp")) == NULL) {
	fputs("Parser error: Please re-run with the option "
	      "-decomp and a valid directory name.\n", stderr);
	param_list_clear(pl);
	exit(EXIT_FAILURE);
    }
    //option: distribution cofactors 
    const char *name_file_cofactor;
    //"/localdisk/trichard/cado768/cofactors";
    if ((name_file_cofactor = param_list_lookup_string(pl, "dist")) == NULL) {
    	fputs("Parser error: Please re-run with the option -dist "
    	      "followed by the pathname of the file which stores the "
    	      "distribution of our cofactors.\n", stderr);
    	param_list_clear(pl);
    	exit(EXIT_FAILURE);
    }

    FILE *file_C = fopen(name_file_cofactor, "r");
    unsigned long **distrib_C = extract_matrix_C(file_C, mfb0 + 1, mfb1 + 1);
    if (distrib_C == NULL) {
    	fprintf(stderr, "Error while reading file %s\n", name_file_cofactor);
    	param_list_clear(pl);
    	exit(EXIT_FAILURE);
    }
    fclose(file_C);

    gmp_randstate_t state;
    gmp_randinit_default(state);

    //select our methods
    //tabular_fm_t *methods = generate_methods_cado(lpb);
    //benchmark
    //bench_proba(state, methods, fbb, 0, 0);
    //bench_time(state, methods, 0);
    //FILE* filee = fopen ("bench_data_cado","w");
    //tabular_fm_fprint (filee, methods);
    //fclose (filee);
    FILE* file_in = fopen ("/localdisk/trichard/cadoRSA155/data_fm_25", "r");
    if (file_in == NULL)
      {
    	fprintf (stderr,
    		 "impossible to read: /localdisk/trichard/cadoRSA155/data_fm_25\n");
    	  exit(1);
      }
    tabular_fm_t *methods = tabular_fm_fscan (file_in);
    /* we set mfb = 3*lpb0 to avoid the special-case of 2 large primes */
    printf ("len  = %d, (%d)\n", methods->index, 3 + nb_curves (lpb0, 3*lpb0));
    
    //test computation of probabilities
    //{tab_init
    /* int fbb0 = ceil (log2 ((double) (lim0+1))); */
    /* int fbb1 = ceil (log2 ((double) (lim1+1))); */

    /* int r0 = 62, r1=63; */
    /* assert (r0 <= mfb0 && r1 <= mfb1); */
    /* printf ("r0 = %d, r1= %d\n", r0, r1); */
    /* char name_file[200]; */
    /* snprintf(name_file, sizeof(name_file), */
    /* 	    "%s/decomp_%lu_%d", name_directory_decomp, lim0, r0);//modify it!! */
    /* printf ("%s\n", name_file); */
    /* fflush(stdout); */
    /* FILE *file = fopen(name_file, "r"); */
    /* tabular_decomp_t* init_tab0 = tabular_decomp_fscan (file); */
    /* fclose (file); */
    /* snprintf(name_file, sizeof(name_file), */
    /* 	    "%s/decomp_%lu_%d", name_directory_decomp, lim1, r1);//modify it!! */
    /* printf ("%s\n", name_file); */
    /* file = fopen(name_file, "r"); */
    /* tabular_decomp_t* init_tab1 = tabular_decomp_fscan (file); */
    /* fclose (file); */
    /* //Test{{{ */
    /* /\* init_tab0->index = 1; *\/ */
    /* /\* init_tab1->index = 1; *\/ */
    /* /\* tabular_decomp_print (init_tab0); *\/ */
    /* //}}} */
    /* //  } */
    /* strategy_t* strat1 = strategy_create(); */
    /* for (int i = 0; i < methods->index; i++) */
    /*   strategy_add_fm (strat1, methods->tab[i]); */
    /* /\* strategy_add_fm (strat1, methods->tab[0]); *\/ */
    /* /\* strategy_add_fm (strat1, methods->tab[1]); *\/ */
    /* /\* strategy_add_fm (strat1, methods->tab[2]); *\/ */
    /* /\* strategy_add_fm (strat1, methods->tab[3]); *\/ */
    /* /\* strategy_add_fm (strat1, methods->tab[4]); *\/ */
    /* /\* strategy_add_fm (strat1, methods->tab[5]); *\/ */
    /* /\* strategy_add_fm (strat1, methods->tab[6]); *\/ */
    /* strategy_t* strat2 = strategy_create(); */
    /* for (int i = 0; i < methods->index; i++) */
    /*   strategy_add_fm (strat2, methods->tab[i]); */

    /* /\* strategy_add_fm (strat2, methods->tab[0]); *\/ */
    /* /\* strategy_add_fm (strat2, methods->tab[1]); *\/ */
    /* /\* strategy_add_fm (strat2, methods->tab[2]); *\/ */
    /* /\* strategy_add_fm (strat2, methods->tab[3]); *\/ */
    /* /\* strategy_add_fm (strat2, methods->tab[4]); *\/ */
    /* /\* strategy_add_fm (strat2, methods->tab[5]); *\/ */
    /* /\* strategy_add_fm (strat2, methods->tab[6]); *\/ */
    /* /\* strategy_print (strat1); *\/ */
    /* /\* strategy_print (strat2); *\/ */
    /* double proba1 = compute_proba_strategy (init_tab0, strat1, fbb0, lpb0); */
    /* double proba2 = compute_proba_strategy (init_tab1, strat2, fbb1, lpb1); */
    /* double temps1 = compute_time_strategy (init_tab0, strat1, r0); */
    /* double temps2 = compute_time_strategy (init_tab1, strat2, r1); */
    /* strategy_set_proba (strat1, proba1); */
    /* strategy_set_proba (strat2, proba2); */
    /* strategy_set_time (strat1, temps1); */
    /* strategy_set_time (strat2, temps2); */
    /* printf ("proba1 = %lf, proba2 = %lf, combo = %lf\n", */
    /* 	    proba1, proba2, proba1*proba2); */
    /* printf ("temps1 = %lf, temps2 = %lf, combo = %lf, %lf\n", */
    /* 	    temps1, temps2, temps1 + proba1*temps2, temps2 + proba2*temps1); */
    /* int fbb[2] = {fbb0, fbb1}; */
    /* int lpb[2] ={lpb0, lpb1}; */
    /* int r[2] = {r0, r1}; */
    /* tabular_decomp_t* init_tab[2] = {init_tab0, init_tab1}; */
    /* //{ */
    /* strat1->len_side = 3; */
    /* strat1->side = malloc (sizeof (int) * strat1->len_side); */
    /* strat1->side[0] = strat1->side[1] = strat1->side[2] =0; */
    /* //  } */
    /* /\* double temps = compute_time_strategy_ileav(init_tab, strat1, fbb, lpb, r); *\/ */
    /* /\* printf ("temps %lf\n", temps); *\/ */
    /* /\* strategy_t *res = gen_strat_r0_r1_ileav_st(strat1, strat2, init_tab, *\/ */
    /* /\* 					       fbb, lpb, r); *\/ */
    /* strategy_t* res = gen_strat_r0_r1_ileav_st_rec(strat1, 0, */
    /* 						   strat2, 0, */
    /* 						   init_tab, fbb, lpb, r, */
    /* 						   NULL, 0); */
    /* printf ("resutat!!!!\n"); */
    /* printf ("theo: proba = %lf, tps = %lf\n", res->proba, res->time); */
    /* //strategy_print (res); */
    /* //bench time with our interleaving methods! */
    /* int mfb[2] = {mfb0, mfb1}; */
    /* printf ("bench proba time\n"); */
    /* double* tmp = bench_proba_time_st_both(state, res, init_tab, r, fbb, lpb, mfb); */
    /* printf ("bench both interL: proba = %lf, time = %lf\n", tmp[0], tmp[1]); */
    /* exit(1); */
    //}}
    
    //generate our strategies
    //remark: for each pair (r0, r1), we have only one strategy!!
#ifndef CADO_INTERLEAVING
    //clear me!!!!
    //compute our strategy
    tabular_strategy_t ***matrix_strat;
    if (0) //my classic matrix!
      {
    	char pathname_st[200] = "../res_matrix";///localdisk/trichard/cado704_new/res_matrix/";
    	matrix_strat =
    	  extract_matrix_strat(pathname_st, mfb0 + 1, mfb1 + 1);
      }
    else //compute interleaving strategies!
      {
    	const char name_directory_str[100] =
    	  "/tmp/res_precompt_st";
    	matrix_strat = generate_matrix_ileav(name_directory_decomp,
    					     name_directory_str,
    					     lim0, lpb0, mfb0,
    					     lim1, lpb1, mfb1);
      }
    printf ("our strategy_file!\n");
    strategy_t ***matrix_strat_res =
      compute_best_strategy(matrix_strat, distrib_C,  mfb0 + 1, mfb1 + 1,
    			    C0);

    /* exit(1); */
    //}}
    
    tabular_strategy_t ***matrix =
      generate_matrix_cado(name_directory_decomp, methods,
    			   lim0, lpb0, mfb0,
    			   lim1, lpb1, mfb1);
    
    //eval our strategy!
    double Y = 0, T = C0;
    for (int r0 = 0; r0 <= mfb0; r0++) {
    	for (int r1 = 0; r1 <= mfb1; r1++) {
    	    Y += distrib_C[r0][r1] * matrix[r0][r1]->tab[0]->proba;
    	    T += distrib_C[r0][r1] * matrix[r0][r1]->tab[0]->time;
	    //test: Add test with our strategy!
	    if (0)//matrix_strat_res[r0][r1] != NULL)
	      {
	    	printf ("cado r0=%d, r1=%d, p =%lf, t = %lf\n", r0, r1,
	    		matrix[r0][r1]->tab[0]->proba,
	    		matrix[r0][r1]->tab[0]->time);
	    	printf ("file r0=%d, r1=%d, p =%lf, t = %lf\n", r0, r1,
	    		matrix_strat_res[r0][r1]->proba,
	    		matrix_strat_res[r0][r1]->time);
	    	printf ("number of pair: %lu\n", distrib_C[r0][r1]);
	    /* 	strategy_print (matrix_strat_res[r0][r1]); */
	    /* //getchar(); */
	    }
	    //end of test
	}
    }
    //print the result!
    printf ("without interleaving with cado!\n");
    printf(" Y = %lf relations, T = %lf s., yt = %1.10lf s/rel\n", Y,
    	   T / 1000000,  T / (Y * 1000000));

    const char *pathname_output;
    pathname_output = param_list_lookup_string(pl, "out");
    FILE *file_output = fopen(pathname_output, "w");
    fprint_final_strategy(file_output, matrix_strat_res,
			  mfb0 + 1, mfb1 + 1);
    fclose (file_output);
    /* if (pathname_output != NULL) { */
    /* 	FILE *file_output = fopen(pathname_output, "w"); */
    /* 	for (int r0 = 0; r0 <= mfb0; r0++) */
    /* 	    for (int r1 = 0; r1 <= mfb1; r1++) */
    /* 		if (matrix[r0][r1]->tab[0] != NULL) { */
    /* 		    fprintf(file_output, */
    /* 			    "[r0=%d, r1=%d] : (p = %lf, t = %lf)\n", */
    /* 			    r0, r1, matrix[r0][r1]->tab[0]->proba, */
    /* 			    matrix[r0][r1]->tab[0]->time); */
    /* 		    strategy_fprint_design(file_output, matrix[r0][r1]->tab[0]); */
    /* 		} */
    /* 	fclose(file_output); */
    /* } */
#else  //interleaving!
    
    tabular_strategy_t ***matrix =
      generate_matrix_cado_ileav(name_directory_decomp, methods,
				 lim0, lpb0, mfb0,
				 lim1, lpb1, mfb1);
    
    strategy_t ***matrix_strat_res =
      compute_best_strategy(matrix, distrib_C, mfb0 + 1, mfb1 + 1,
			    C0);

    const char *pathname_output;
    pathname_output = param_list_lookup_string(pl, "out");

    if (pathname_output != NULL)
      {
	FILE *file_output = fopen(pathname_output, "w");
	fprint_final_strategy(file_output, matrix_strat_res,
			      mfb0 + 1, mfb1 + 1);
	fclose (file_output);
      }
    //free matrix_strat_res
    for (int r0 = 0; r0 <= mfb0; r0++){
      for (int r1 = 0; r1 <= mfb1; r1++)
	strategy_free(matrix_strat_res[r0][r1]);
      free(matrix_strat_res[r0]);
    }
    free(matrix_strat_res);
#endif

    //{{test bench strategy!
    /* int r0 = CONST_TEST_R; */
    /* int r1 = CONST_TEST_R; */
    /* printf ("r0 = %d, r1 = %d, lpb0 = %d, lpb1 = %d\n", r0, r1, lpb0, lpb1); */
    /* //printf ("nb call = %lu\n", distrib_C[r0][r1]); */
    /* printf ("proba found: %lf, %lf\n", matrix[r0][r1]->tab[0]->proba, matrix[r0][r1]->tab[0]->time); */
    /* char name_file[200]; */
    /* snprintf(name_file, sizeof(name_file), */
    /* 	    "%s/decomp_%d_%d", name_directory_decomp, fbb0, r0); */
    /* FILE *file = fopen(name_file, "r"); */
    
    /* tabular_decomp_t* tab_decomp = tabular_decomp_fscan(file); */
    
    /* tabular_decomp_print (tab_decomp); */
    /* fclose (file); */
    
    /* facul_strategy_t* facul_st = facul_make_strategy (fbb0, lpb0, 0, 0); */

    /* double* res = bench_proba_time_st(state, facul_st, tab_decomp, r0, lpb0); */
    /* double p0 = res[0], t0 = res[1]; */
    /* printf ("side = 0, proba = %lf, time = %lf\n", res[0], res[1]); */
    /* tabular_decomp_free (tab_decomp); */
    /* facul_clear_strategy(facul_st); */
    /* //r1 */
    /* snprintf(name_file, sizeof(name_file), */
    /* 	    "%s/decomp_%d_%d", name_directory_decomp, fbb1, r1); */
    /* file = fopen(name_file, "r"); */
    
    /* tab_decomp = tabular_decomp_fscan(file); */
    /* fclose (file); */
    
    /* facul_st = facul_make_strategy (fbb1, lpb1, 0, 0); */
    /* res = bench_proba_time_st(state, facul_st, tab_decomp, r1, lpb1); */
    /* printf ("side = 1, proba = %lf, time = %lf\n", res[0], res[1]); */
    /* printf ("two sides, proba = %lf, time = %lf\n", res[0]*p0, t0+p0*res[1]); */
    /* tabular_decomp_free (tab_decomp); */
    /* //Free facil_st */
    /* facul_clear_strategy(facul_st); */
    //}}
    
    //free
    for (int r0 = 0; r0 <= mfb0; r0++) {
	for (int r1 = 0; r1 <= mfb1; r1++)
	    tabular_strategy_free(matrix[r0][r1]);
	free(distrib_C[r0]);
	free(matrix[r0]);
    }
    free(distrib_C);
    free(matrix);
    gmp_randclear(state);
    param_list_clear(pl);

    return EXIT_SUCCESS;
    
}
