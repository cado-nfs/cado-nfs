#include "cado.h"

#include <float.h>
#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>


#include "portability.h"
#include "utils.h"
#include "facul.h"
#include "pm1.h"
#include "pp1.h"
#include "ecm.h"
#include "getprime.h"

#include "generate_factoring_method.h"
#include "convex_hull.h"

/*
  BOUND_SIGMA is used when you generate a random value of sigma.
*/
const int BOUND_SIGMA = 15;

/* 
   The test fails with probabily equals to 1/4^15 (approx 0) (MILLER RABIN).
*/
const int NB_PRIMALITY_TEST = 15;

const double EPSILON_DBL = LDBL_EPSILON;

/*
  If the probability to find a prime number of p bits is less than 3%, 
  we can suppose that the probability to find p+i bits (i>0) is equal to 0.
*/
const double BENCH_MIN_PROBA = 0.03;

/************************************************************************/
/*                      To generate numbers n=pq                        */
/************************************************************************/

/*
  this function allows to compute an array with the probability of 
  scoring of each size of prime number in [min, ..., max].
 */
double *distribution_prime_number(int min, int max)
{
    int len = max - min +1;
    ASSERT(len >= 0);

    double *tab = malloc(len * sizeof(double));
    ASSERT(tab != NULL);

    unsigned long elem = pow(2, min);
    int log_elem = min;
    for (int i = 0; i < len; i++) {
	tab[i] = (elem - 1) / log_elem;
	elem *= 2;
	log_elem++;
    }

    double proba = 1;
    double sum = 0;
    for (int i = 0; i < len; i++) {
	tab[i] *= proba;
	proba /= 2;
	sum += tab[i];
    }

    ASSERT(sum != 0);
    for (int i = 0; i < len; i++)
	tab[i] /= sum;

    return tab;
}

/*
  This function generates in 'res' a random prime number of 'lenFact' bits.
*/

void generate_prime_factor(mpz_t res, gmp_randstate_t state, int lenFact)
{
    mpz_t min;
    mpz_init(min);
    mpz_ui_pow_ui(min, 2, lenFact - 1);
    do {
	mpz_urandomm(res, state, min);
	mpz_add(res, res, min);
    }
    while (mpz_probab_prime_p(res, NB_PRIMALITY_TEST) == 0);

    //clear
    mpz_clear(min);
}

/*
  This function generates in 'res' a random composite number of 'lenFactall' 
  bits with a prime factor of 'lenFact1' bits.
*/
void
generate_composite_integer(mpz_t res, gmp_randstate_t state,
			   int lenFact1, int lenFactall)
{
    mpz_t bound_max, bound_min;
    mpz_init(bound_max);
    mpz_init(bound_min);
    mpz_ui_pow_ui(bound_min, 2, lenFactall - 1);	//2^(lenfactall)
    mpz_mul_ui(bound_max, bound_min, 2);
    mpz_sub_ui(bound_max, bound_max, 1);	//2^(lenfactall+1)-1

    mpz_t p;
    mpz_init(p);

    mpz_t q;
    mpz_init(q);
    do {
	/*p */
	generate_prime_factor(p, state, lenFact1);
	/*q */
	generate_prime_factor(q, state, lenFactall - lenFact1);

	/* compute n */
	mpz_mul(res, q, p);

    } while (mpz_cmp(bound_min, res) > 0 || mpz_cmp(bound_max, res) < 0);

    // clear
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(bound_max);
    mpz_clear(bound_min);
}


int
select_random_index_according_dist(double *dist, int len)
{
    int alea = rand(); /* 0 <= alea <= RAND_MAX */
    int i = 0;
    int bound = (int) (dist[0] * (double) RAND_MAX);
    while (i < (len-1) && alea > bound) {
	i++;
	bound += (int) (dist[i] * (double) RAND_MAX);
    }
    return i;
}

/*
  This function generates in 'res' a random composite number of 'lenFactall' 
  bits with a prime factor where this size depends to the distrubution 
  of prime numbers (in 'dist')!
*/
int
generate_composite_integer_interval(mpz_t res, gmp_randstate_t state,
				    double *dist, int lenFact1_min,
				    int lenFact1_max, int lenFactall)
{
    int len = lenFact1_max -lenFact1_min +1;
    int index = select_random_index_according_dist(dist, len);
    generate_composite_integer(res, state, lenFact1_min + index, lenFactall);
    return lenFact1_min + index;
}

/************************************************************************/
/*                To model our factoring methods                        */
/************************************************************************/

/* 
   This function allows to create a type facul_strategy_t which contains all
   informations for our factoring method.
   This is necessary, when you would use the function facul () to realize 
   our bench.
*/

facul_strategy_t *generate_fm(int method, int curve, unsigned long B1,
			      unsigned long B2)
{
    facul_strategy_t *strategy;
    strategy = malloc(sizeof(facul_strategy_t));
    strategy->methods = malloc(2 * sizeof(facul_method_t));
    /*
      without this second method, the function
      facul_clear_strategy (strategy) failed!
    */
    strategy->methods[1].method = 0;
    strategy->methods[1].plan = NULL;

    strategy->lpb = ULONG_MAX;
    strategy->assume_prime_thresh = 0;
    strategy->BBB = 0;
    strategy->methods[0].method = method;

    if (method == PM1_METHOD) {
	strategy->methods[0].plan = malloc(sizeof(pm1_plan_t));
	ASSERT(strategy->methods[0].plan != NULL);
	pm1_make_plan(strategy->methods[0].plan, B1, B2, 0);
    } else if (method == PP1_METHOD) {
	strategy->methods[0].plan = malloc(sizeof(pp1_plan_t));
	ASSERT(strategy->methods[0].plan != NULL);
	pp1_make_plan(strategy->methods[0].plan, B1, B2, curve, 0);
    } else if (method == EC_METHOD) {
	long sigma;
	if (curve == MONTY16)
	  sigma = 1;
	else if (curve == BRENT12)
	  sigma = 11;
	else 
	  sigma = 2 + rand()%BOUND_SIGMA;

	strategy->methods[0].plan = malloc(sizeof(ecm_plan_t));
	ASSERT(strategy->methods[0].plan != NULL);
	ecm_make_plan(strategy->methods[0].plan, B1, B2, curve,
		      labs(sigma), 1, 0);
    } else {
      exit(EXIT_FAILURE);
    }
    return strategy;
}

/************************************************************************/
/*            BENCH THE PROBABILITY AND TIME OF OUR METHODS             */
/************************************************************************/

/*
  This function collects the probability of success to find a prime number 
  (of 'len_p' bits) in a composite number (of 'len_n' bits) 
  with the strategy.
*/

double
bench_proba_fm(facul_strategy_t * strategy, gmp_randstate_t state,
	       unsigned long len_p, unsigned long len_n, mpz_t* N,
	       int nb_test_max)
{
    int nb_success_max = 1000, nb_success = 0;
    int nb_test = 0;

    while (nb_success < nb_success_max && (nb_test < nb_test_max)) {
	/* 
	   f will contain the prime factor of N that the strategy
	   found.  Note that N is composed by two prime factors by the
	   previous function.
	*/
	unsigned long f[2];
	if (mpz_cmp_ui (N[nb_test], 0) == 0) {
	    generate_composite_integer(N[nb_test], state, len_p, len_n);
	    mpz_init (N[nb_test+1]);
	    mpz_set_ui (N[nb_test+1], 0);
	}
	f[0] = 0;
	facul(f, N[nb_test], strategy);

	if (f[0] != 0)
	    nb_success++;
	nb_test++;
    }

    return nb_success / ((double)nb_test);
}

double
bench_time_fm_onelength(facul_strategy_t * method, mpz_t* N, int nb_test)
{
    double tps = 0;
    unsigned long f[2];

    uint64_t starttime, endtime;
    starttime = microseconds();

    for (int i = 0; i < nb_test; i++)
	facul(f, N[i], method);
    
    endtime = microseconds();
    tps += endtime - starttime;

    return tps / ((double)nb_test);

}

/*
  This function allows to collect for each factoring methods, 
  the probability to find a prime number of a certain size, 
  from 'len_p_min' until the probability become null. 
*/

void bench_proba(gmp_randstate_t state, tabular_fm_t * fm, int len_p_min)
{
    int len = fm->index;	//number of methods!
    int p_max = 100;
    double *proba = malloc(p_max * sizeof(double));
    ASSERT(proba != NULL);

    unsigned long *param;
    fm_t *elem;
    //{{Will contain the our composite integers!
    int nb_test_max = 10000;
    mpz_t*N[p_max];
    for (int i = 0; i < p_max; i++)
	{
	    N[i] = malloc(sizeof (mpz_t) * (nb_test_max+1));
	    ASSERT (N[i] != NULL);
	    mpz_init(N[i][0]);
	    mpz_set_ui (N[i][0], 0);
	}
    //}}
    
    for (int i = 0; i < len; i++) {
	elem = tabular_fm_get_fm(fm, i);
	param = fm_get_method(elem);
	unsigned long method = param[0];
	unsigned long curve = param[1];
	unsigned long B1 = param[2];
	unsigned long B2 = param[3];

	facul_strategy_t *st = generate_fm(method, curve, B1, B2);

	int ind_proba = 0;
	do {
	    if (B1 == 0 && B2 == 0)
		proba[ind_proba] = 0;
	    else {
		int len_p = len_p_min + ind_proba;
		int len_n = 60 + len_p;
		proba[ind_proba] = bench_proba_fm(st, state, len_p,
						  len_n, N[ind_proba],
						  nb_test_max);
	    }
	    ind_proba++;
	} while ((proba[ind_proba - 1] - BENCH_MIN_PROBA) > EPSILON_DBL
		 && ind_proba < p_max);

	fm_set_proba(elem, proba, ind_proba, len_p_min);
	//free
	facul_clear_strategy(st);
    }
    //free
    for (int i = 0; i < p_max; i++) {
	if (N[i] != NULL) {
	    int j = 0;
	    while(mpz_cmp_ui(N[i][j], 0) != 0)
		{
		    mpz_clear(N[i][j]);
		    j++;
		}
	    mpz_clear(N[i][j]);
	}
	free(N[i]);
    }
    free(proba);
}

/*
  This function allows to collect for each factoring methods, the
  differents time to find a prime number in a interger of different
  bits size (MODREDCUL_MAXBITS, MODREDC15UL_MAXBITS,
  MODREDC2UL2_MAXBITS and MODREDC2UL2_MAXBITS+30).
*/
void bench_time(gmp_randstate_t state, tabular_fm_t * fm)
{
    unsigned long *param;
    fm_t *elem;
    int len = fm->index;	//number of methods!
    //precompute 4 arrays for our bench_time!
    //{{
    int nb_test = 100000;
    mpz_t* N1 = malloc(sizeof(mpz_t) * nb_test);
    mpz_t* N2 = malloc(sizeof(mpz_t) * nb_test);
    mpz_t* N3 = malloc(sizeof(mpz_t) * nb_test);
    mpz_t* N4 = malloc(sizeof(mpz_t) * nb_test);
    for (int i = 0; i < nb_test; i++) {
	mpz_init(N1[i]);
	mpz_init(N2[i]);
	mpz_init(N3[i]);
	mpz_init(N4[i]);
	generate_prime_factor(N1[i], state, MODREDCUL_MAXBITS);
	generate_prime_factor(N2[i], state, MODREDC15UL_MAXBITS);
	generate_prime_factor(N3[i], state, MODREDC2UL2_MAXBITS);
	generate_prime_factor(N4[i], state, MODREDC2UL2_MAXBITS +30);
    }
    //}}
    for (int i = 0; i < len; i++) {
	elem = tabular_fm_get_fm(fm, i);
	param = fm_get_method(elem);
	unsigned long method = param[0];
	unsigned long curve = param[1];
	unsigned long B1 = param[2];
	unsigned long B2 = param[3];
	if (B1 != 0 || B2 != 0) {
	    facul_strategy_t *st = generate_fm(method, curve, B1, B2);
	    double time[4];
	    time[0]= bench_time_fm_onelength(st, N1, nb_test);
	    time[1]= bench_time_fm_onelength(st, N2, nb_test);
	    time[2]= bench_time_fm_onelength(st, N3, nb_test);
	    time[3]= bench_time_fm_onelength(st, N4, nb_test);
	    fm_set_time(elem, time, 4);
	    //free
	    facul_clear_strategy(st);
	} else {
	    double time[4] = { 0, 0, 0, 0 };
	    fm_set_time(elem, time, 4);
	}
    }
    //free
    for (int i = 0; i < nb_test; i++) {
	mpz_clear (N1[i]);
	mpz_clear (N2[i]);
	mpz_clear (N3[i]);
	mpz_clear (N4[i]);
    }
    free (N1);
    free (N2);
    free (N3);
    free (N4);
}

/************************************************************************/
/*       BENCH according to an interval of size of prime numbers        */
/************************************************************************/

/* 
   This function returns default parameters for the collection of 
   factoring methods and stores their in 'res' such that: 
   res[0]=b1_min 
   res[1]=b1_max
   res[2]=b1_step 
   res[3]=c_min 
   res[4]=c_max 
   res[5]=c_step 
   TODO: these parameters will be improved :)
*/
int *choice_parameters(int method, int len_p_min)
{
    int b1_min, b1_max, b1_step, c_min, c_max, c_step;
    if (method == PM1_METHOD ||
	method == PP1_METHOD) {
	if (len_p_min < 20) {
	    b1_min = 20;
	    b1_max = 1000;
	    b1_step = 10;
	    c_min = 10;
	    c_max = 50;
	    c_step = 5;
	} else if (len_p_min <= 25) {
	    b1_min = 20;
	    b1_max = 4500;
	    b1_step = 20;
	    c_min = 40;
	    c_max = 80;
	    c_step = 5;
	} else if (len_p_min < 30) {
	    b1_min = 100;
	    b1_max = 10000;
	    b1_step = 50;
	    c_min = 20;
	    c_max = 200;
	    c_step = 10;
	} else {
	    b1_min = 100;
	    b1_max = 30000;
	    b1_step = 100;
	    c_min = 20;
	    c_max = 500;
	    c_step = 10;
	}

    } else if (method == EC_METHOD) {
	if (len_p_min < 20) {
	    b1_min = 20;
	    b1_max = 5000;
	    b1_step = 10;
	    c_min = 10;
	    c_max = 100;
	    c_step = 10;
	} else if (len_p_min <= 25) {
	    b1_min = 100;
	    b1_max = 10000;
	    b1_step = 10;
	    c_min = 10;
	    c_max = 200;
	    c_step = 40;
	} else if (len_p_min < 30) {
	    b1_min = 100;
	    b1_max = 13000;
	    b1_step = 100;
	    c_min = 10;
	    c_max = 500;
	    c_step = 20;
	} else {
	    b1_min = 200;
	    b1_max = 30000;
	    b1_step = 100;
	    c_min = 10;
	    c_max = 1000;
	    c_step = 20;
	}

    } else {
	return NULL;
    }
    int *result = malloc(6 * sizeof(int));
    ASSERT(result != NULL);
    result[0] = b1_min;
    result[1] = b1_max;
    result[2] = b1_step;
    result[3] = c_min;
    result[4] = c_max;
    result[5] = c_step;
    return result;
}
/*
  This function allows to compute the probability and the time of a strategy
  to find a prime number in an interval [2**len_p_min, 2**len_p_max].
*/
static double *bench_proba_time_pset_onefm(facul_strategy_t *strategy,
					   mpz_t* N, int nb_test_max)
{
    double *res = malloc(2 * sizeof(double));
    ASSERT(res != NULL);
    int nb_succes_lim = 1000, nb_succes = 0;
    int nb_test = 0;

    double starttime, endtime;
    starttime = microseconds();

    while ((nb_succes < nb_succes_lim) && (nb_test < nb_test_max)) {

      //computes the the time of execution
	/*
	   f will contain the prime factor of N that the strategy
	   found.  Note that N is composed by two prime factors by the
	   previous function.
	 */
	unsigned long f[2];
	f[0] = 0;

	facul(f, N[nb_test], strategy);

	if (f[0] != 0)
	    nb_succes++;
	nb_test++;
    }
    endtime = microseconds();
    double tps = endtime - starttime;

    //clear

    res[0] = nb_succes / ((double)nb_test);
    res[1] = tps / ((double)nb_test);
    return res;
}

/*
  This function returns an array which contains, for each method, 
  the probability and the cost to find a prime number for each size between 
  'len_p_min' and 'len_p_max'. The type of these methods is specified by 
  two parameters : 'method' and 'curve', and the parameters B1 and B2 are 
  given by 'param_region' or, if it's equal to NULL, use the default sieve.
*/

tabular_fm_t *bench_proba_time_pset(int method, int curve,
				    gmp_randstate_t state,
				    int len_p_min, int len_p_max,
				    int len_n, int *param_region)
{
    //define the sieve region
    int c_min, c_max;
    int b1_min, b1_max;
    int c_pas, b1_pas;
    if (param_region != NULL) {
	b1_min = param_region[0];
	b1_max = param_region[1];
	b1_pas = param_region[2];

	c_min = param_region[3];
	c_max = param_region[4];
	c_pas = param_region[5];
    } else {
	//default parameters for the sieve region.
	int *param = choice_parameters(method, len_p_min);
	b1_min = param[0];
	b1_max = param[1];
	b1_pas = param[2];

	c_min = param[3];
	c_max = param[4];
	c_pas = param[5];
	free(param);
    }
    
    //{{Will contain the our composite integers!
    int nb_test_max = 10000;
    mpz_t*N = malloc(sizeof (mpz_t) * (nb_test_max));
    ASSERT (N != NULL);
    double *disp = distribution_prime_number(len_p_min, len_p_max);
    for (int i = 0; i < nb_test_max; i++)
	{
	    /* generation of the integers N[i] (which will be
	       factoring). To avoid to waste time, we precompute these
	       values only one time 
	    */
	    mpz_init(N[i]);
	    generate_composite_integer_interval(N[i], state, disp, len_p_min,
						len_p_max, len_n);
	}
    //}}

    tabular_fm_t *tab_fusion = tabular_fm_create();

    //add zero method
    unsigned long tmp_method[4] = { method, 0, 0, 0 };
    double zero = 0;
    tabular_fm_add(tab_fusion, tmp_method, 4, &zero, 1, &zero, 1, len_p_min);

    for (int c = c_min; c <= c_max; c += c_pas) {
	int B1;
	int B2;

	B1 = b1_min;
	B2 = B1 * c;

	tabular_fm_t *tab = tabular_fm_create();
	unsigned long elem[4];
	double proba = 0;
	double tps = 0;
	double max_proba = 0.9;

	while (B1 <= b1_max && proba < max_proba) {
	    facul_strategy_t *fm = generate_fm(method, curve, B1, B2);
	    double* res = bench_proba_time_pset_onefm(fm, N, nb_test_max);
	    proba = res[0];
	    tps = res[1];
	    free (res);
	    
	    elem[0] = method;
	    elem[1] = curve;
	    elem[2] = B1;
	    elem[3] = B2;

	    tabular_fm_add(tab, elem, 4, &proba, 1, &tps, 1, len_p_min);

	    facul_clear_strategy(fm);

	    B1 = B1 + b1_pas;
	    B2 = B1 * c;

	}
	//merge arrays
	tabular_fm_concat(tab_fusion, tab);
	tabular_fm_free(tab);
    }
    //free
    free (disp);
    for (int i = 0; i < nb_test_max; i++)
	mpz_clear(N[i]);
    free(N);
   
    return tab_fusion;
}

/************************************************************************/
/*                     GENERATE FACTORING METHODS                       */
/************************************************************************/

/*
  This function returns for each sub-interval of prime number, the probability
  and the time to find a prime number in N (an interger of len_n bits).
  With this function, we can specify the method, the curve and the sieve 
  region of parameters for B1 and B2. Note that an option 'opt_ch' 
  allows to apply the selection by convex hull.
*/
tabular_fm_t *generate_factoring_methods_mc(gmp_randstate_t state,
					    int len_p_min, int len_p_max,
					    int len_n, int method, int curve,
					    int opt_ch, int *param_sieve)
{
    ASSERT(len_p_min <= len_p_max);

    tabular_fm_t *gfm;

    tabular_fm_t *collect = bench_proba_time_pset
	(method, curve, state, len_p_min, len_p_max, len_n, param_sieve);

    if (opt_ch) {
	//apply the convex hull
	gfm = convex_hull_fm(collect);
	tabular_fm_free(collect);
    } else 
	gfm = collect;

    return gfm;
}

/*
  This function calls the previous function (generate_factoring_methods_mc ())
  for each one of our functions.
*/

tabular_fm_t *generate_factoring_methods(gmp_randstate_t state, int len_p_min,
					 int len_p_max, int len_n, int opt_ch,
					 int *param_sieve)
{

    tabular_fm_t *gfm = tabular_fm_create();

    //we begin by the first method:
    int ind_method = 0;
    int ind_curve = 0;
    int curve[3] = {MONTY12, MONTY16, BRENT12};
    int method[4] = {PM1_METHOD, PP1_METHOD, EC_METHOD};
    while (ind_method <= 3 && ind_curve < 3) {

	printf("method = %d, curve = %d\n",
	       method[ind_method], curve[ind_curve]);
	tabular_fm_t *res = generate_factoring_methods_mc
	    (state, len_p_min, len_p_max, len_n, method[ind_method],
	     curve[ind_curve], opt_ch, param_sieve);

	tabular_fm_concat(gfm, res);

	//free
	tabular_fm_free(res);

	//index
	if (method[ind_method] != EC_METHOD)
	    ind_method++;
	else
	    ind_curve++;
    }

    return gfm;
}

/************************************************************************/
/*                    MERGE and compute the convex hull                 */
/************************************************************************/
/*
  This function collects factoring methods from the file 'file_in' 
  and make a selection by convex hull and prints them in file_out.
 */
tabular_fm_t *convex_hull_from_file(FILE * file_in, FILE * file_out)
{
    tabular_fm_t *all_st = tabular_fm_fscan(file_in);
    if (all_st == NULL)
	return NULL;

    tabular_fm_t *res = convex_hull_fm(all_st);

    tabular_fm_free(all_st);

    int err = tabular_fm_fprint(file_out, res);
    if (err < 0)
	return NULL;

    return res;
}

/************************************************************************/
/*                      FILTERING                                       */
/************************************************************************/

static int
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
filtering: most homegenous method that allows to keep
methods of differents probabilities. With a classic version, the
remaining methods tend to have very high probabilities, and it's not
necessarily that we want!
*/
tabular_fm_t *filtering(tabular_fm_t * tab, int final_nb_methods)
{
    //create the matrix with the average dist between a pair of methods!
    int nb_methods = tab->index;
    double **dist = malloc(nb_methods * sizeof(double *));
    ASSERT(dist != NULL);
    for (int i = 0; i < nb_methods; i++) {
	dist[i] = malloc((nb_methods) * sizeof(double));
	ASSERT(dist[i] != NULL);
    }
    for (int i = 0; i < nb_methods; i++) {
	dist[i][i] = 0;
	//trade off for i.
	fm_t *eli = tabular_fm_get_fm(tab, i);
	double compromis_i[eli->len_proba];
	for (int p = 0; p < eli->len_proba; p++) {
	    if (eli->proba[p] > EPSILON_DBL)
		compromis_i[p] = eli->time[get_nb_word(p)] / eli->proba[p];
	    else //if (eli->time[p] < EPSILON_DBL): fm zero!
		compromis_i[p] = 0;

	}
	for (int j = i + 1; j < nb_methods; j++) {
	    //trade off for j.
	    fm_t *elj = tabular_fm_get_fm(tab, j);
	    double compromis_j[elj->len_proba];
	    for (int p = 0; p < elj->len_proba; p++) {
		if (elj->proba[p] > EPSILON_DBL)
		    compromis_j[p] = elj->time[get_nb_word(p)] / elj->proba[p];
		else //if (elj->time[p] < EPSILON_DBL): fm zero!
		    compromis_j[p] = 0;
	    }

	    //compute dist
	    double moy_dist = 0;
	    int nb_elem = (elj->len_proba < eli->len_proba)?
		elj->len_proba:eli->len_proba;
	    for (int p = 0; p < nb_elem; p++) {
		double tmp = compromis_i[p] - compromis_j[p];
		tmp *=tmp;
		moy_dist += tmp;
	    }
	    moy_dist = sqrt(moy_dist)/nb_elem;

	    dist[i][j] = moy_dist;
	    dist[j][i] = moy_dist;
	}
    }

    //sort the pairs of methods according to the dist!
    int nb_pair = (nb_methods-1)*(nb_methods)/2;
    int sort_dist[nb_pair][2];
    //todo: improve this method! For now, it's a naive method.
    int k = 0;
    while (k < nb_pair) {
	int i_min = -1;
	int j_min = -1;
	double ratio = INFINITY;
	for (int i = 0; i < nb_methods; i++) {
	    for (int j = i + 1; j < nb_methods; j++) {
		if (dist[i][j] < ratio)
		    {
			i_min = i;
			j_min = j;
			ratio = dist[i][j];
		    }
	    }
	}
	if (i_min == -1 || j_min == -1)
	    break;
	sort_dist[k][0] = i_min;
	sort_dist[k][1] = j_min;
	dist[i_min][j_min] = INFINITY;
	k++;
    }
    nb_pair = k;

    //clear method until you have the good numbers of methods.
    int* tab_fm_is_removed = calloc (nb_methods, sizeof (int));
    int nb_rem_methods = nb_methods;

    while(nb_rem_methods > final_nb_methods){
	int ind = 0;
	while (ind < nb_pair && nb_rem_methods > final_nb_methods)
	    {
		if (tab_fm_is_removed[sort_dist[ind][0]] == 0 &&
		    tab_fm_is_removed[sort_dist[ind][1]] == 0)
		    {
			fm_t* el0 = tabular_fm_get_fm(tab, sort_dist[ind][0]);
			fm_t* el1 = tabular_fm_get_fm(tab, sort_dist[ind][1]);
			double ratio0 = (el0->proba[0] < EPSILON_DBL)?
			    0:el0->time[0]/el0->proba[0];
			double ratio1 = (el1->proba[0] < EPSILON_DBL)?
			    0:el1->time[0]/el1->proba[0];

			if ( ratio0 > ratio1)
			    {
				tab_fm_is_removed[sort_dist[ind][0]] = -1;
				tab_fm_is_removed[sort_dist[ind][1]] = 1;
			    }
			else
			    {
				tab_fm_is_removed[sort_dist[ind][1]] = -1;
				tab_fm_is_removed[sort_dist[ind][0]] = 1;
			    }
			nb_rem_methods--;
		    }
		ind++;
	    }
	for (int i = 0; i < nb_methods; i++)
	    if (tab_fm_is_removed[i] == 1)
		tab_fm_is_removed[i] = 0;
    }
    //build the final tab_fm with the remaining factoring methods!
    tabular_fm_t *res = tabular_fm_create();
    for (int i = 0; i < nb_methods; i++)
	{
	    if (tab_fm_is_removed[i] != -1)
		tabular_fm_add_fm (res, tabular_fm_get_fm(tab, i));
	}
    //free
    free (tab_fm_is_removed);
    for (int i = 0; i < nb_methods; i++)
	free (dist[i]);
    free (dist);

    tabular_fm_sort (res);

    return res;
}










/************************************************************************/
/*                      CONVEX_HULL_FM                                  */
/************************************************************************/

/*
  These differents functions allow to use the module convex_hull, 
  to compute the convex hull of a set of factoring methods.
 */

tabular_point_t *convert_tab_point_to_tab_fm(tabular_fm_t * t)
{
    tabular_point_t *res = tabular_point_create();
    int len = tabular_fm_get_index(t);
    fm_t *elem;
    for (int i = 0; i < len; i++) {
	elem = tabular_fm_get_fm(t, i);
	tabular_point_add(res, i, elem->proba[0], elem->time[0]);
    }
    return res;
}

tabular_fm_t *convert_tab_fm_to_tab_point(tabular_point_t * t,
					  tabular_fm_t * init)
{
    tabular_fm_t *res = tabular_fm_create();
    int len = tabular_point_get_index(t);
    for (int i = 0; i < len; i++) {
	int index = point_get_number(tabular_point_get_point(t, i));
	tabular_fm_add_fm(res, init->tab[index]);
    }
    return res;
}

tabular_fm_t *convex_hull_fm(tabular_fm_t * t)
{
    tabular_point_t *tmp = convert_tab_point_to_tab_fm(t);
    tabular_point_t *res = convex_hull(tmp);
    tabular_fm_t *res_fm = convert_tab_fm_to_tab_point(res, t);
    tabular_point_free(tmp);
    tabular_point_free(res);
    return res_fm;
}
