#ifndef GENERATE_STRATEGIES_HPP
#define GENERATE_STRATEGIES_HPP

#include "tab_decomp.hpp"
#include "tab_strategy.hpp"
#include "tab_point.hpp"
#include "decomp.hpp"
#include "fm.hpp"
#include "tab_fm.hpp"
#include "strategy.hpp"

/************************************************************************/
/*                      COLLECT DATA FOR ONLY ONE COFACTOR              */
/************************************************************************/

double compute_proba_method_one_decomp (decomp const & dec, fm_t const * fm);

double compute_proba_strategy(tabular_decomp const & init_tab, strategy_t * strat,
			      unsigned int len_p_min, unsigned int len_p_max);

double compute_time_strategy(tabular_decomp const & init_tab, strategy_t * strat, 
			     unsigned int r);

/************************************************************************/
/*                   GENERATE MATRIX                                    */
/************************************************************************/

tabular_strategy_t *generate_strategies_oneside(tabular_decomp const & tab_decomp,
						fm_t * zero,
						tabular_fm_t * pm1,
						tabular_fm_t * pp1,
						tabular_fm_t * ecm,
						int nb_curve, unsigned long lim,
						unsigned int lpb, unsigned int r);

tabular_strategy_t *generate_strategy_r0_r1(tabular_strategy_t * strat_r0,
					    tabular_strategy_t * strat_r1);

tabular_strategy_t ***generate_matrix(const char *name_directory_decomp,
				      tabular_fm_t * pm1, tabular_fm_t * pp1,
				      tabular_fm_t * ecm, int nb_curve, 
				      unsigned long lim0, unsigned int lpb0, unsigned int mfb0,
				      unsigned long lim1, unsigned int lpb1, unsigned int mfb1);

/************************************************************************/
/*                      CONVEX_HULL_ST                                  */
/************************************************************************/

tabular_point convert_tab_point_to_tab_strategy(tabular_strategy_t * t);

tabular_strategy_t *convert_tab_strategy_to_tab_point(tabular_point const & t,
						      tabular_strategy_t * init);

tabular_strategy_t *convex_hull_strategy(tabular_strategy_t * t);

#endif				/* GENERATE_STRATEGIES_HPP */
