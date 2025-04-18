#ifndef FINDING_GOOD_STRATEGY_HPP
#define FINDING_GOOD_STRATEGY_HPP

#include <cstdio>
#include "strategy.hpp"
#include "tab_strategy.hpp"

tabular_strategy_t ***extract_matrix_strat(const char *pathname_st,
					   unsigned int len_abs,
					   unsigned int len_ord);

unsigned long **extract_matrix_C(FILE * file, unsigned int len_abs, unsigned int len_ord);

strategy_t ***compute_best_strategy(tabular_strategy_t *** matrix_strat,
				    unsigned long **distrib_C,
				    unsigned int len_abs, unsigned int len_ord, double C0);


//to print our final strategies!
void strategy_fprint_design(FILE * output_file, const strategy_t * t);

int
fprint_final_strategy(FILE * file, strategy_t *** matrix_strat_res,
		      unsigned int len_abs, unsigned int len_ord);

#endif				/* FINDING_GOOD_STRATEGY_HPP */
