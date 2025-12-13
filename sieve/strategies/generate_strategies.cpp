#include "cado.h" // IWYU pragma: keep

#include <cmath>
#include <cstdlib>
#include <cfloat>

#include <fstream>

#include "arith/modredc_ul.h"
#include "convex_hull.hpp"
#include "decomp.hpp"
#include "facul_ecm.h"
#include "facul_method.hpp"
#include "fm.hpp"
#include "generate_strategies.hpp"
#include "macros.h"
#include "point.hpp"
#include "strategy.hpp"
#include "tab_decomp.hpp"
#include "tab_fm.hpp"
#include "tab_point.hpp"
#include "tab_strategy.hpp"
#include "utils_cxx.hpp"



/************************************************************************/
/*                      COLLECT DATA FOR ONLY ONE COFACTOR              */
/************************************************************************/
/*
  Compute the probability to doesn't find a non-trivial factor when we
  use this method 'fm'.
 */
double
compute_proba_method_one_decomp (decomp const & D, fm_t const * fm)
{
    const double *proba_suc = fm_get_proba(fm);
    unsigned int const len_proba = fm_get_len_proba(fm);
    double proba_fail = 1; 
    for (auto j : D) {
	j -= fm->len_p_min;
	if (j < len_proba)
	    proba_fail *= 1 - proba_suc[j];
	//else //the probability seems closest to 0.
    }   
    return proba_fail;
}

/*
  Compute the probability to find a non-trivial factor in a good decomposition.
 */
double
compute_proba_strategy(tabular_decomp const & init_tab, strategy_t * strat,
		       unsigned int len_p_min, unsigned int len_p_max)
{
    double all = 0.0;
    double nb_found_elem = 0;
    unsigned int const nb_fm = tabular_fm_get_size(strat->tab_fm);
    tabular_fm_t *tab_fm = strat->tab_fm;

    for (auto const & D : init_tab) {
	if (is_good_decomp(D, len_p_min, len_p_max)) {
	    //the probability to doesn't find a non trivial factor
	    //with all methods in our strategy.
	    double p_fail_all = 1;
	    for (unsigned int index_fm = 0; index_fm < nb_fm; index_fm++) {
		fm_t const * elem = tabular_fm_get_fm(tab_fm, index_fm);
		double const p_fail_one = compute_proba_method_one_decomp (D, elem);
		p_fail_all *= p_fail_one;
		if (elem->method[0] == PM1_METHOD ||
		    elem->method[0] == PP1_27_METHOD ||
		    elem->method[0] == PP1_65_METHOD)
		    //because if you chain PP1||PM1 to PM1||PP1-->they are
		    //not independant.
		    p_fail_all = (p_fail_one + p_fail_all) / 2;
	    }
	    nb_found_elem += (1 - p_fail_all) * D.nb_elem;
	}
	all += D.nb_elem;
    }
    if (all < (double)LDBL_EPSILON) // all == 0.0 -> it exists any decomposition!
      return 0;
    return nb_found_elem / all;
}

/*
  Compute the average time when we apply our strategy 'strat' in a
  cofactor of r bits!
*/
double compute_time_strategy(tabular_decomp const & init_tab, strategy_t * strat, unsigned int r)
{
    unsigned int const nb_fm = tabular_fm_get_size(strat->tab_fm);
    tabular_fm_t *tab_fm = strat->tab_fm;
    //{{
    /*
      We add 0.5 to the length of one word, because for our times the
      lenght is inclusive. For example, if MODREDCUL_MAXBITS = 64
      bits, a cofactor is in one word if is lenght is less OR equal to
      64 bits. So, if you don't add 0.5 to MODREDCUL_MAXBITS, you
      lost the equal and thus insert an error in your maths. 
    */
    double const half_word = (MODREDCUL_MAXBITS+0.5)/2.0;
    unsigned int const number_half_wd = floor(r /half_word);
    //the next computation is necessary in the relation with the
    //benchmark in gfm!
    unsigned int const ind_time = (number_half_wd <2)? 0: number_half_wd - 1;
    //}}

    double time_average = 0;
    //store the number of elements in the different decompositions!
    double all = 0.0;
    for (auto const & D : init_tab) {
	double time_dec = 0;
	double proba_fail_all = 1;
	double time_method = 0;
	//compute the time of each decomposition
	for (unsigned int index_fm = 0; index_fm < nb_fm; index_fm++) {
	    fm_t const * elem = tabular_fm_get_fm(tab_fm, index_fm);
	    unsigned int const len_time = fm_get_len_time (elem);
	    if (ind_time >= len_time)
	      time_method = elem->time[len_time-1];
	    else
	      time_method = elem->time[ind_time];
	    time_dec += time_method * proba_fail_all;
	    
	    double const proba_fail_method = 
	      compute_proba_method_one_decomp (D, elem);
	    proba_fail_all *= proba_fail_method;
	    if (elem->method[0] == PM1_METHOD ||
		elem->method[0] == PP1_27_METHOD ||
		elem->method[0] == PP1_65_METHOD)
	      //because if you chain PP1||PM1 to PM1||PP1-->they are
	      //not independant.
	      proba_fail_all = (proba_fail_all + proba_fail_method) / 2;
	}

	time_average += time_dec * D.nb_elem;
	all += D.nb_elem;
    }
    if (all < (double)LDBL_EPSILON) // all == 0.0 -> it exists any decomposition!
      return 0;
    
    return time_average / all;
}

/*
  As their name suggests, this function adds one strategy to our array
  't' without the zero methods.
*/
static void
tabular_strategy_add_strategy_without_zero(tabular_strategy_t * t,
					   strategy_t * strategy)
{
    if (t->size >= t->alloc)
	tabular_strategy_realloc(t);

    strategy_t *elem = strategy_create();
    unsigned int const len = strategy->tab_fm->size;
    int strat_is_zero = true;
    for (unsigned int i = 0; i < len; i++) {
	if (!tabular_fm_is_zero(strategy->tab_fm, i)) {
	    strat_is_zero = false;
	    tabular_fm_add_fm(elem->tab_fm, strategy->tab_fm->tab[i]);
	}
    }
    if (strat_is_zero)
	tabular_fm_add_fm(elem->tab_fm, strategy->tab_fm->tab[0]);

    elem->proba = strategy->proba;
    elem->time = strategy->time;
    t->tab[t->size] = elem;
    t->size++;
}

/************************************************************************/
/*                   GENERATE MATRIX                                    */
/************************************************************************/
/*
  This function add differents chains of ecm to the strategy 'strat'.
  Note that the variable 'lbucket' allows to add ecm by bucket of
  length 'lbucket'.
 */
static void
generate_collect_iter_ecm(fm_t * zero, tabular_fm_t * ecm,
			  unsigned int ind_ecm, strategy_t * strat,
                          unsigned int ind_tab,
			  unsigned int index_iter, unsigned int len_iteration,
                          unsigned int lbucket,
			  tabular_decomp const & init_tab,
			  tabular_strategy_t**all_strat_ptr,
			  unsigned int fbb, 
                          unsigned int lpb,
                          unsigned int r, int is_already_used_B12M16)
{
    tabular_strategy_t * all_strat = *all_strat_ptr;
    if (index_iter >= len_iteration) {
	unsigned int const nb_strat = all_strat->size;
	tabular_strategy_add_strategy_without_zero(all_strat, strat);
	double const proba =
	    compute_proba_strategy(init_tab, all_strat->tab[nb_strat], fbb,lpb);
	double const time = 
	    compute_time_strategy(init_tab, all_strat->tab[nb_strat], r);

	strategy_set_proba(all_strat->tab[nb_strat], proba);
	strategy_set_time(all_strat->tab[nb_strat], time);
    } else {
	for (unsigned int i = ind_ecm; i < ecm->size; i++) {
	  /* The curve BRENT12 and MONTY16 are curves with only one
	  sigma. So use it only one time.*/
	  if (ecm->tab[i]->method[1] == MONTY16 || //MONTY16
	      ecm->tab[i]->method[1] == BRENT12) //BRENT12
	    {
	      if (is_already_used_B12M16)
		continue;
	      tabular_fm_set_fm_index(strat->tab_fm, ecm->tab[i], ind_tab);
	      generate_collect_iter_ecm(zero, ecm, i+1, strat,ind_tab+1,
					index_iter + 1, len_iteration,
					lbucket, init_tab, all_strat_ptr,
					fbb, lpb, r, true);			
	    }
	  else //MONTY12
	    {
	      for (unsigned int j = 0; j < lbucket && ind_tab+j < strat->tab_fm->size; j++)
		{
		  tabular_fm_set_fm_index(strat->tab_fm, ecm->tab[i], ind_tab+j);
		}
	      generate_collect_iter_ecm(zero, ecm, i, strat,ind_tab+lbucket,
					index_iter + lbucket, len_iteration,
					lbucket+1, init_tab, all_strat_ptr,
					fbb, lpb,r, is_already_used_B12M16);
	    }	
	}
        all_strat = *all_strat_ptr; // might have changed during recursive call
    	/* to protect the ram, we reduce the number of strategies by
	   the convex hull when this number become too big.*/
	if (all_strat->size < 100000) {
	    tabular_strategy_t* ch = convex_hull_strategy(all_strat);
	    //clear previous collect and start a new collect.
	    tabular_strategy_free(all_strat);
	    *all_strat_ptr = tabular_strategy_create();
	    tabular_strategy_concat(*all_strat_ptr, ch);
	    tabular_strategy_free(ch);
	}
    }
}

/*
  This function allows to generate strategies, computes their data,
  and selects the convex hull of these strategies for one size of
  cofactor 'r'.  It allows to avoid to full the RAM when we generated
  strategies!  
  Moreover, this generator chains our factoring methods
  such that: 
 PM1 (0/1) + PP1 (0/1) + ECM-M12(0/1/2...)+ ECM-M16/B12(0/1)+ ECM-M12(0/1/...)
*/
tabular_strategy_t *generate_strategies_oneside(tabular_decomp const & init_tab,
						fm_t * zero, tabular_fm_t * pm1,
						tabular_fm_t * pp1,
						tabular_fm_t * ecm,
						int ncurves, unsigned long lim,
						unsigned int lpb, unsigned int r)
{
    //contains the final result
    tabular_strategy_t *res;

    //check the cases where r is trivial!!
    //{{
    unsigned int const fbb = ceil (log2 ((double) (lim + 1)));
    unsigned int const lim_is_prime = 2 * fbb - 1;

    ASSERT_ALWAYS(!init_tab.empty() == (r >= lim_is_prime));

    /*
      In this case, r is already a prime number!
      We define two zero strategies to manage:
      - the case where r is a good prime number (fbb<r<lpb)
      |-> the probability is equal to 1.
      - the case where r is a prime number too big (lpb < r <
      fbb^2) or because the lenght of r is impossible (r<fbb).
      - |-> the probability is equal to 0.
    */
    if (r < lim_is_prime) {
	strategy_t *st_zero = strategy_create();
	strategy_add_fm(st_zero, zero);
	strategy_set_time(st_zero, 0.0);

	if (r != 1 && (r < fbb || r > lpb))
	    strategy_set_proba(st_zero, 0);
	else	// r==1 or fbb<= r0 <= lpb
	    strategy_set_proba(st_zero, 1.0);
	res = tabular_strategy_create();
	tabular_strategy_add_strategy(res, st_zero);
	strategy_free(st_zero);
	return res;
    }
    //}}
    //contains strategies which will be processed.
    tabular_strategy_t *all_strat = tabular_strategy_create();

    unsigned int const len_pm1 = pm1->size;
    unsigned int const len_pp1 = pp1->size;

    //init strat
    int const len_strat = 2 + ncurves;
    //printf ("len_strat = %d, ncurve = %d\n",len_strat, ncurves);
    strategy_t *strat = strategy_create();
    for (int i = 0; i < len_strat; i++)
	strategy_add_fm(strat, zero);

    tabular_fm_t *tab_strat = strat->tab_fm;
    //PM1
    for (unsigned int ind_pm1 = 0; ind_pm1 < len_pm1; ind_pm1++) {
      double const current_proba_pm1 = pm1->tab[ind_pm1]->proba[0];
      tabular_fm_set_fm_index(tab_strat, pm1->tab[ind_pm1], 0);
      //PP1      
      for (unsigned int ind_pp1 = 0; ind_pp1 < len_pp1; ind_pp1++) {
        if ( pp1->tab[ind_pp1]->method[2] != 0 //B1==0
	     && pp1->tab[ind_pp1]->proba[0] < current_proba_pm1)
	  continue;
	tabular_fm_set_fm_index(tab_strat, pp1->tab[ind_pp1], 1);

	//ECM (M12-B12-M16)
	generate_collect_iter_ecm(zero, ecm, 0, strat, 2,
				  0, ncurves, 0, init_tab,
				  &all_strat, fbb, lpb, r, false);
      }
    } 
    //compute the final convex hull.
    res = convex_hull_strategy(all_strat);
    tabular_strategy_free(all_strat);
    strategy_free(strat);
    return res;
}

/*
  As their name suggests, this function concatenates two strategies but
  take into account the side of each strategy: st1 will be the
  first_side and st2 the other side.
*/

static strategy_t *concat_strategies(strategy_t * st1, strategy_t * st2,
				     int first_side)
{
    strategy_t *st = strategy_create();
    int const len1 = st1->tab_fm->size;
    int const len2 = st2->tab_fm->size;
    st->len_side = len1 + len2;
    if (st->side == nullptr)
	st->side = (int*) malloc(sizeof(int) * (st->len_side));
    else
        checked_realloc(st->side, st->len_side);
    int side = first_side;
    for (int i = 0; i < len1; i++) {
	strategy_add_fm(st, st1->tab_fm->tab[i]);
	st->side[i] = side;
    }
    side = first_side ? 0 : 1;
    for (int i = 0; i < len2; i++) {
	strategy_add_fm(st, st2->tab_fm->tab[i]);
	st->side[len1 + i] = side;
    }
    return st;
}

/*
  returns the best strategies to factor a couple (r0, r1), from a set
  of optimal strategies for each side. Note that, the probability and
  the time to find a non-trivial for each side must be previously
  computed!
 */
tabular_strategy_t *generate_strategy_r0_r1(tabular_strategy_t * strat_r0,
					    tabular_strategy_t * strat_r1)
{
    unsigned int const len_r0 = strat_r0->size;
    unsigned int const len_r1 = strat_r1->size;

    tabular_strategy_t *strat_r0_r1 = tabular_strategy_create();
    tabular_strategy_t *ch = tabular_strategy_create();
    unsigned long nb_strat = 0;

    /*
       for each array of strategies, the first one is the zero strategy.
       There are two cases : 
       -first, we have a success at 0%.
       -Secondly, we have a success at 100%, because the cofactor is already
       prime and not too big.  The four lines below allow to avoid to
       have several unless methods with a zero probability.
     */

    strategy_t *st;
    for (unsigned int r = 0; r < len_r0; r++)	//first side
    {
	double const p0 = strat_r0->tab[r]->proba;
	double const c0 = strat_r0->tab[r]->time;
	for (unsigned int a = 0; a < len_r1; a++)	//second side
	{
	    nb_strat++;
	    //compute success ans cost:
	    double const p1 = strat_r1->tab[a]->proba;
	    double const c1 = strat_r1->tab[a]->time;
	    double const proba = p0 * p1;
	    double const tps0 = c0 + p0 * c1;
	    //time when we begin with side 0
	    double const tps1 = c1 + p1 * c0;
            //time when we begin with side 1
	    if (tps0 < tps1) {
		st = concat_strategies(strat_r0->tab[r], strat_r1->tab[a], 0);
		strategy_set_proba(st, proba);
		strategy_set_time(st, tps0);
	    } else {
		st = concat_strategies(strat_r1->tab[a], strat_r0->tab[r], 1);
		strategy_set_proba(st, proba);
		strategy_set_time(st, tps1);
	    }
	    tabular_strategy_add_strategy(strat_r0_r1, st);
	    strategy_free(st);
	}

	//process data to avoid to full the RAM!!!
	if (nb_strat > 100000 || r == (len_r0 - 1)) {
	    //add strategies of the old convexhull and recompute the new
	    //convex hull

	    tabular_strategy_concat(strat_r0_r1, ch);
	    tabular_strategy_free(ch);

	    ch = convex_hull_strategy(strat_r0_r1);
	    //clear previous collect and start a new collect.
	    tabular_strategy_free(strat_r0_r1);
	    strat_r0_r1 = tabular_strategy_create();
	    nb_strat = 0;
	}
    }

    //free
    tabular_strategy_free(strat_r0_r1);

    return ch;
}

/*
  returns the best strategies for each couple of cofactors of lenght
  (r0, r1), from a set of factoring methods. Note that this function
  use the previous functions, and need all probabilities and
  times for each method must be previously computed. (to do that, you
  could use the binary gfm).
 */

tabular_strategy_t ***generate_matrix(const char *name_directory_decomp,
				      tabular_fm_t* pm1, tabular_fm_t* pp1,
				      tabular_fm_t*ecm, int ncurves,
				      unsigned long lim0, unsigned int lpb0, unsigned int mfb0,
				      unsigned long lim1, unsigned int lpb1, unsigned int mfb1)
{

    unsigned int const fbb0 = ceil (log2 ((double) (lim0 + 1)));
    unsigned int const fbb1 = ceil (log2 ((double) (lim1 + 1)));

    /*
       allocates the matrix which contains all optimal strategies for
       each couple (r0, r1): r0 is the lenght of the cofactor in
       the first side, and r1 for the second side.
     */
    auto * matrix = new tabular_strategy_t **[mfb0 + 1];

    for (unsigned int r0 = 0; r0 <= mfb0; r0++)
	matrix[r0] = new tabular_strategy_t*[mfb1 + 1];

    /*
       For each r0 in [fbb0..mfb0], we precompute and strore the data
       for all strategies.  Whereas, for each r1, the values will be
       computed sequentially. 
     */
    fm_t *zero = fm_create();
    unsigned long method_zero[4] = { 0, 0, 0, 0 };
    fm_set_method(zero, method_zero, 4);

    auto *data_rat = new tabular_strategy_t*[mfb0 + 1];
    ASSERT (data_rat);

    unsigned int lim_is_prime = 2 * fbb0 - 1;
    for (unsigned int r0 = 0; r0 <= mfb0; r0++) {
	tabular_decomp tab_decomp;
	if (r0 >= lim_is_prime) {
            auto filename = fmt::format("{}/decomp_{}_{}", name_directory_decomp, lim0, r0);
            std::ifstream is(filename);
            if (!(is >> tab_decomp)) {
                fmt::print(stderr, "Cannot read {}\n", filename);
		exit(EXIT_FAILURE);
	    }
	}
	data_rat[r0] =
	    generate_strategies_oneside(tab_decomp, zero, pm1, pp1,
					ecm, ncurves, lim0, lpb0, r0);
    }

    /*
       read good elements for r_2 in the array data_r1 and compute the
       data for each r_1. So :
     */
    lim_is_prime = 2 * fbb1 - 1;
    for (unsigned int r1 = 0; r1 <= mfb1; r1++) {
	tabular_decomp tab_decomp;
	if (r1 >= lim_is_prime) {
            auto filename = fmt::format("{}/decomp_{}_{}", name_directory_decomp, lim1, r1);
            std::ifstream is(filename);
            if (!(is >> tab_decomp)) {
                fmt::print(stderr, "Cannot read {}\n", filename);
		exit(EXIT_FAILURE);
	    }
	}

	tabular_strategy_t *strat_r1 =
	  generate_strategies_oneside(tab_decomp, zero, pm1, pp1,
				      ecm, ncurves, lim1, lpb1, r1);

	for (unsigned int r0 = 0; r0 <= mfb0; r0++) {
	    matrix[r0][r1] = generate_strategy_r0_r1(data_rat[r0], strat_r1);
	}
	tabular_strategy_free(strat_r1);
    }

    //free
    for (unsigned int r0 = 0; r0 <= mfb0; r0++)
	tabular_strategy_free(data_rat[r0]);
    delete[] data_rat;

    fm_free(zero);
    return matrix;
}

/************************************************************************/
/*                      CONVEX_HULL_FM                                  */
/************************************************************************/

/*
  These following functions allow to use the module convex_hull, to
  compute the convex hull of a set of strategies.
 */

tabular_point convert_tab_point_to_tab_strategy(tabular_strategy_t * t)
{
    tabular_point res;
    strategy_t *elem;
    for (unsigned int i = 0; i < t->size; i++) {
	elem = t->tab[i];
	res.emplace_back(point { .number=i, .x=elem->proba, .y=elem->time });
    }
    return res;
}

tabular_strategy_t *convert_tab_strategy_to_tab_point(tabular_point const & t,
						      tabular_strategy_t * init)
{
    tabular_strategy_t *res = tabular_strategy_create();
    for (auto const & p : t) {
	tabular_strategy_add_strategy(res, init->tab[p.number]);
    }
    return res;
}

tabular_strategy_t *convex_hull_strategy(tabular_strategy_t * t)
{
    tabular_point res = convex_hull(convert_tab_point_to_tab_strategy(t));
    tabular_strategy_t *res_strat = convert_tab_strategy_to_tab_point(res, t);
    return res_strat;
}
