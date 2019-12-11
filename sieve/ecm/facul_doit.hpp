#ifndef FACUL_DOIT_H
#define FACUL_DOIT_H

#include "modredc_ul.h"
#include "modredc_15ul.h"
#include "modredc_2ul2.h"
#include "mod_mpz.h"
#include "facul.hpp"

int facul_doit_ul (std::vector<cxx_mpz> &, const modulusredcul_t, 
		   const facul_strategy_t *, const int);
int facul_doit_15ul (std::vector<cxx_mpz> &, const modulusredc15ul_t, 
		     const facul_strategy_t *, const int);
int facul_doit_2ul2 (std::vector<cxx_mpz> &, const modulusredc2ul2_t, 
		     const facul_strategy_t *, const int);
int facul_doit_mpz (std::vector<cxx_mpz> &, const modulusmpz_t, 
		    const facul_strategy_t *, const int);

int
facul_doit_onefm_ul (std::vector<cxx_mpz> &, const modulusredcul_t,
		     const facul_method_t, const FaculModulusBase * &,
		     const FaculModulusBase * &, unsigned long, double, double);
int
facul_doit_onefm_15ul (std::vector<cxx_mpz> &, const modulusredc15ul_t,
		       const facul_method_t, const FaculModulusBase * &,
		       const FaculModulusBase * &, unsigned long, double, double);
int
facul_doit_onefm_2ul2 (std::vector<cxx_mpz> &, const modulusredc2ul2_t,
		       const facul_method_t, const FaculModulusBase * &,
		       const FaculModulusBase * &, unsigned long, double, double);
int
facul_doit_onefm_mpz (std::vector<cxx_mpz> &, const modulusmpz_t,
		      const facul_method_t, const FaculModulusBase * &,
		      const FaculModulusBase * &, unsigned long, double, double);

/* int* */
/* facul_both (unsigned long**, mpz_t* , */
/* 	    const facul_strategies_t *); */

#endif /* FACUL_DOIT_H */
