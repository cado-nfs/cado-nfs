#ifndef GENERATE_STRATEGIES
#define GENERATE_STRATEGIES

#include "tab_decomp.h"
#include "tab_strategy.h"
#include "tab_point.h"
#include "decomp.h"
#include "fm.h"
#include "tab_fm.h"
#include "strategy.h"

#ifdef __cplusplus
extern "C" {
#endif

/************************************************************************/
/*                      COLLECT DATA FOR ONLY ONE COFACTOR              */
/************************************************************************/

double
compute_proba_method_one_decomp (decomp_t* dec, fm_t* fm);

double compute_proba_strategy(tabular_decomp_t * init_tab, strategy_t * strat,
			      int len_p_min, int len_p_max);

double compute_time_strategy(tabular_decomp_t * init_tab, strategy_t * strat, 
			     int r);

/************************************************************************/
/*                   GENERATE MATRIX                                    */
/************************************************************************/

tabular_strategy_t *generate_strategies_oneside(tabular_decomp_t * tab_decomp,
						fm_t * zero,
						tabular_fm_t * pm1,
						tabular_fm_t * pp1,
						tabular_fm_t * ecm,
						int nb_curve, unsigned long lim,
						int ub, int r);

tabular_strategy_t *generate_strategy_r0_r1(tabular_strategy_t * strat_r0,
					    tabular_strategy_t * strat_r1);

tabular_strategy_t ***generate_matrix(const char *name_directory_decomp,
				      tabular_fm_t * pm1, tabular_fm_t * pp1,
				      tabular_fm_t * ecm, int nb_curve, 
				      unsigned long lim0, int lpb0, int mfb0,
				      unsigned long lim1, int lpb1, int mfb1);

/************************************************************************/
/*                      CONVEX_HULL_ST                                  */
/************************************************************************/

tabular_point_t *convert_tab_point_to_tab_strategy(tabular_strategy_t * t);

tabular_strategy_t *convert_tab_strategy_to_tab_point(tabular_point_t * t,
						      tabular_strategy_t *
						      init);

tabular_strategy_t *convex_hull_strategy(tabular_strategy_t * t);

#ifdef __cplusplus
}
#endif

#endif				/* GENERATE_STRATEGIES */
