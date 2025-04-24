#ifndef FACUL_DOIT_HPP
#define FACUL_DOIT_HPP

#include <memory>
#include <vector>

#include "arith/modredc_ul.h"
#include "arith/modredc_15ul.h"
#include "arith/modredc_2ul2.h"
#include "arith/mod_mpz.h"
#include "facul.hpp"
#include "facul_method.hpp"
#include "modset.hpp"

facul_status
facul_doit_onefm (std::vector<cxx_mpz> &, const modulusredcul_t,
		  facul_method const &,
                  std::vector<std::unique_ptr<FaculModulusBase>> &,
		  unsigned long, double, double);
facul_status
facul_doit_onefm (std::vector<cxx_mpz> &, const modulusredc15ul_t,
		  facul_method const &,
                  std::vector<std::unique_ptr<FaculModulusBase>> &,
		  unsigned long, double, double);
facul_status
facul_doit_onefm (std::vector<cxx_mpz> &, const modulusredc2ul2_t,
		  facul_method const &,
                  std::vector<std::unique_ptr<FaculModulusBase>> &,
		  unsigned long, double, double);
facul_status
facul_doit_onefm (std::vector<cxx_mpz> &, const modulusmpz_t,
		  facul_method const &,
                  std::vector<std::unique_ptr<FaculModulusBase>> &,
		  unsigned long, double, double);

/* int* */
/* facul_both (unsigned long**, mpz_t* , */
/* 	    const facul_strategies_t *); */

#endif /* FACUL_DOIT_HPP */
