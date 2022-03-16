#ifndef FACUL_DOIT_H
#define FACUL_DOIT_H

#include "modredc_ul.h"
#include "modredc_15ul.h"
#include "modredc_2ul2.h"
#include "mod_mpz.h"
#include "facul.hpp"
#include "modset.hpp"

int facul_doit (std::vector<cxx_mpz> &, const modulusredcul_t, 
		facul_strategy_oneside const &, const int);
int facul_doit (std::vector<cxx_mpz> &, const modulusredc15ul_t, 
		facul_strategy_oneside const &, const int);
int facul_doit (std::vector<cxx_mpz> &, const modulusredc2ul2_t, 
		facul_strategy_oneside const &, const int);
int facul_doit (std::vector<cxx_mpz> &, const modulusmpz_t, 
		facul_strategy_oneside const &, const int);

int
facul_doit_onefm (std::vector<cxx_mpz> &, const modulusredcul_t,
		  facul_method const &, const FaculModulusBase * &,
		  const FaculModulusBase * &, unsigned long, double, double);
int
facul_doit_onefm (std::vector<cxx_mpz> &, const modulusredc15ul_t,
		  facul_method const &, const FaculModulusBase * &,
		  const FaculModulusBase * &, unsigned long, double, double);
int
facul_doit_onefm (std::vector<cxx_mpz> &, const modulusredc2ul2_t,
		  facul_method const &, const FaculModulusBase * &,
		  const FaculModulusBase * &, unsigned long, double, double);
int
facul_doit_onefm (std::vector<cxx_mpz> &, const modulusmpz_t,
		  facul_method const &, const FaculModulusBase * &,
		  const FaculModulusBase * &, unsigned long, double, double);

/* int* */
/* facul_both (unsigned long**, mpz_t* , */
/* 	    const facul_strategies_t *); */

#endif /* FACUL_DOIT_H */
