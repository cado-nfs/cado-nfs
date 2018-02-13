#include "cado.h"

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "facul_ecm.h"
#include "facul.hpp"
#include "portability.h"
#include "utils.h"
#include "finding_good_strategy.h"

//#define STATS 
/* if define then print in a file 'result_strat' the different values
of s that are tested. */

static double EPSILON_DBL = LDBL_EPSILON;

/************************************************************************/
/*                    EXTRACT MATRIX STRATREGIES                        */
/************************************************************************/

/*
  This function read the directory 'file_name' to extract the
  strategies for each pair (r0,r1); Note that 'len_abs' and 'len_ord'
  are the size of the matrix.
*/

tabular_strategy_t ***extract_matrix_strat(const char *pathname_st,
					   const int len_abs, const int len_ord)
{

    //allocate matrix
    tabular_strategy_t ***matrix = (tabular_strategy_t***) malloc(sizeof(*matrix) * len_abs);
    ASSERT(matrix != NULL);
    for (int r1 = 0; r1 < len_abs; r1++) {
	matrix[r1] = (tabular_strategy_t**) malloc(sizeof(*matrix[r1]) * len_ord);
	assert(matrix[r1] != NULL);
    }

    for (int r0 = 0; r0 < len_abs; r0++) {
	for (int r1 = 0; r1 < len_ord; r1++) {
	    char name_file[strlen(pathname_st) + 64];
	    snprintf(name_file, sizeof(name_file), "%s/strategies_%d_%d", pathname_st, r0, r1);
	    FILE *file = fopen(name_file, "r");
	    matrix[r0][r1] = tabular_strategy_fscan(file);
	    if (matrix[r0][r1] == NULL) {
		fprintf(stderr,
			"Impossible to read the file '%s'\n", name_file);
		return NULL;
	    }
	    fclose(file);
	}
    }
    return matrix;
}

/************************************************************************/
/*                         EXTRACT THE SET C                            */
/************************************************************************/

/*
  This function extracts the distribution of the size of cofactors
  given by the file 'file' ans store the result in a matrix
  [len_abs][len_ord]. Note that if the size of cofactors are greater
  than len_abs or len_ord, then they aren't stored in the matrix.
*/
unsigned long **extract_matrix_C(FILE * file, int len_abs, int len_ord)
{
    if (file == NULL)
	return NULL;
    //create matrix
    unsigned long **matrix_call = (unsigned long**) malloc(sizeof(*matrix_call) * len_abs);
    assert(matrix_call != NULL);

    for (int i = 0; i < len_abs; i++) {
	matrix_call[i] = (unsigned long *) calloc(sizeof(*matrix_call[i]), len_ord);
	assert(matrix_call[i]);
    }

    //collect data
    int i, j;
    unsigned long c, unused_s;
    while (!feof(file)) {
	if (fscanf(file, "%d %d %lu %lu\n", &i, &j, &c, &unused_s) != 4)
	    return NULL;
	if (i < len_abs && j < len_ord)
	    matrix_call[i][j] = c;
    }

    return matrix_call;
}

/************************************************************************/
/*          SEARCH THE BEST STRATEGIES TO MAXIMIZE Y/T                  */
/************************************************************************/

static int
subroutine_dicho(const tabular_strategy_t * tab_strat, double s, int ind_min,
		 int ind_max, double slope_min, double slope_max)
{
    if (ind_max - ind_min <= 1)
	return ind_min;

    int middle = floor((ind_max + ind_min) / 2);
    double slope_middle;
    if (middle == 0)
	slope_middle = INFINITY;
    else
	{
	    strategy_t *elem1 = tab_strat->tab[middle-1];
	    strategy_t *elem2 = tab_strat->tab[middle];
	    slope_middle = 
		(elem2->proba - elem1->proba) / (elem2->time - elem1->time)*1000000;
	}

    if (slope_middle < s)
	return subroutine_dicho(tab_strat, s, ind_min, middle, slope_min,
				slope_middle);
    else
	return subroutine_dicho(tab_strat, s, middle, ind_max, slope_middle,
				slope_max);
}

static int
subroutine_compute_slope_yt_dicho(const tabular_strategy_t * tab_strat,
				  double s)
{
    int len = tab_strat->index;
    if (len == 1)
	return 0;
    else
	return subroutine_dicho(tab_strat, s, 0, len, INFINITY, 0);
}


static double
compute_slope_yt(tabular_strategy_t *** matrix_strat, unsigned long **distrib_C,
		 int len_abs, int len_ord, double s, double C0)
{
    double Y = 0, T = C0;

    for (int r1 = 0; r1 < len_abs; r1++) {
	for (int r2 = 0; r2 < len_ord; r2++) {
	    if (distrib_C[r1][r2] > EPSILON_DBL) {
		int index =
		    subroutine_compute_slope_yt_dicho(matrix_strat[r1][r2], s);
		Y += distrib_C[r1][r2] *
		    matrix_strat[r1][r2]->tab[index]->proba;
		T += distrib_C[r1][r2] * matrix_strat[r1][r2]->tab[index]->time;
	    }
	}
    }
    return Y / T;
}


static double
sampling_function_interval(tabular_strategy_t *** matrix_strat,
			   unsigned long **distrib_C, int len_abs, int len_ord,
			   double C0, double init_s, double maxi_s, double pas)
{
#ifdef STATS
    //print in this file the différent value of s and the value of yt
    //found for each s.
    FILE *result_file = fopen("result_strat", "w+");
#endif
    double max_s = 0, max_yt = 0;
    for (double s = init_s; s < maxi_s; s += pas) {
	double yt =
	    compute_slope_yt(matrix_strat, distrib_C, len_abs, len_ord, s, C0);
	if (yt > max_yt) {
	    max_s = s;
	    max_yt = yt;
	}
#ifdef STATS
	fprintf(result_file, "%lf \t %1.10lf\n", s, yt);
#endif
    }
    return max_s;
#ifdef STATS
    fclose(result_file);
#endif
}


static double
sampling_function(tabular_strategy_t *** matrix_strat,
		  unsigned long **distrib_C, int len_abs, int len_ord,
		  double C0, double init_s, double pas)
{
#ifdef STATS
    //print in this file the différent value of s and the ratio found
    //for each s.
    FILE *result_file = fopen("result_strat", "w+");
#endif
    double max_s = 0, max_yt = 0;
    double s = init_s;
    int chronos = 0;
    while (chronos < 100) {
	double yt =
	    compute_slope_yt(matrix_strat, distrib_C, len_abs, len_ord, s, C0);
	if (yt > max_yt) {
	    chronos = 0;
	    max_s = s;
	    max_yt = yt;
	} else
	    chronos++;
	s += pas;
#ifdef STATS
	fprintf(result_file, "%lf \t %1.10lf\n", s, yt);
#endif
    }
#ifdef STATS
    fclose(result_file);
#endif

    return max_s;
}

/*
  This function looks for the best choice of s!  
  And Our goal is to choice a good value for s that optimates the 
  number of relations per seconde in the sieving step.
  (cofactor + precomputation sieve).
*/
strategy_t ***compute_best_strategy(tabular_strategy_t *** matrix_strat,
				    unsigned long **distrib_C, int len_abs,
				    int len_ord, double C0)
{
    double step_s = 1;
    //find a interesting interval to search optimal value for s.
    double max_s = sampling_function(matrix_strat, distrib_C,
				     len_abs, len_ord, C0, 0, step_s);

    //search the optimal value of s.
    double min = (max_s - step_s > 0) ? max_s - step_s : 0;

    max_s = sampling_function_interval(matrix_strat, distrib_C,
				       len_abs, len_ord,
				       C0, min, min + 2*step_s, 0.00010);

    //printf("Our choice for the slope s: %lf (rel/s)\n", max_s);
    double Y = 0, T = C0;
    //build the matrix with the optimal strategies.
    strategy_t ***matrix_res = (strategy_t***) malloc(sizeof(*matrix_res) * len_abs);
    assert(matrix_res != NULL);
    for (int r1 = 0; r1 < len_abs; r1++) {
	matrix_res[r1] = (strategy_t**) malloc(sizeof(*matrix_res[r1]) * len_ord);
	assert(matrix_res[r1] != NULL);
	for (int r2 = 0; r2 < len_ord; r2++) {
	    if (distrib_C[r1][r2] < EPSILON_DBL)
		matrix_res[r1][r2] = NULL;
	    else {
		int index =
		    subroutine_compute_slope_yt_dicho(matrix_strat[r1][r2],
						      max_s);
		matrix_res[r1][r2] =
		    strategy_copy(matrix_strat[r1][r2]->tab[index]);

		Y += distrib_C[r1][r2] *
		    matrix_strat[r1][r2]->tab[index]->proba;

		T += distrib_C[r1][r2] * matrix_strat[r1][r2]->tab[index]->time;
	    }
	}
    }
    printf(" Y = %lf relations, T = %lf s., yt = %1.10lf s/rel.\n", Y,
	   T / 1000000, T /(Y * 1000000));

    return matrix_res;
}

/************************************************************************/
/*                  DESIGN OUR RESULT                                   */
/************************************************************************/
/*
  These two functions allow to print in a file our final strategies
  for each pair (r0,r1) where at least one call has been done by CADO
  in the distribution of cofactors.
 */
void strategy_fprint_design(FILE * output_file, const strategy_t * t)
{
    tabular_fm_t *tmp = t->tab_fm;
    for (int i = 0; i < tmp->index; i++) {
	fprintf(output_file, "[S%d: ", t->side[i]);

	fm_t *fm = tmp->tab[i];
	//method: type, curve
	if (fm->method[0] == PP1_27_METHOD)
	    fprintf(output_file, "PP1-27, %lu, %lu ] ", fm->method[2],
		    fm->method[3]);
	else if (fm->method[0] == PP1_65_METHOD)
	    fprintf(output_file, "PP1-65, %lu, %lu ] ", fm->method[2],
		    fm->method[3]);
	else if (fm->method[0] == PM1_METHOD)
	    fprintf(output_file, "PM1, %lu, %lu ] ", fm->method[2],
		    fm->method[3]);
	else			//==EC_METHOD
	{
	    if (fm->method[1] == BRENT12)
		fprintf(output_file, "ECM-B12, %lu, %lu ] ", fm->method[2],
			fm->method[3]);
	    else if (fm->method[1] == MONTY12)
		fprintf(output_file, "ECM-M12, %lu, %lu ] ", fm->method[2],
			fm->method[3]);
	    else		//==MONTY16
		fprintf(output_file, "ECM-M16, %lu, %lu ] ", fm->method[2],
			fm->method[3]);
	}
    }
    fprintf(output_file, "\n");
}

int fprint_final_strategy(FILE * file, strategy_t *** matrix_strat_res,
			  const int len_abs, const int len_ord)
{
    if (file == NULL)
	return -1;
    for (int r0 = 0; r0 < len_abs; r0++)
	for (int r1 = 0; r1 < len_ord; r1++)
	    if (matrix_strat_res[r0][r1] != NULL) {
		fprintf(file,
			"[r0=%d, r1=%d] : (p = %lf, t = %lf)\n",
			r0, r1, matrix_strat_res[r0][r1]->proba,
			matrix_strat_res[r0][r1]->time);
		strategy_fprint_design(file, matrix_strat_res[r0][r1]);
	    }
    return 0;
}
