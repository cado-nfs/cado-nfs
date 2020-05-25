#ifndef FINDING_GOOD_STRATEGY
#define FINDING_GOOD_STRATEGY

#include <stdio.h>
#include "strategy.h" // strategy_t
#include "tab_strategy.h"

#ifdef __cplusplus
extern "C" {
#endif


tabular_strategy_t ***extract_matrix_strat(const char *pathname_st,
					   const int len_abs,
					   const int len_ord);

unsigned long **extract_matrix_C(FILE * file, int len_abs, int len_ord);

strategy_t ***compute_best_strategy(tabular_strategy_t *** matrix_strat,
				    unsigned long **distrib_C,
				    int len_abs, int len_ord, double C0);


//to print our final strategies!
void strategy_fprint_design(FILE * output_file, const strategy_t * t);

int
fprint_final_strategy(FILE * file, strategy_t *** matrix_strat_res,
		      const int len_abs, const int len_ord);

#ifdef __cplusplus
}
#endif

#endif				/* FINDING_GOOD_STRATEGY */
