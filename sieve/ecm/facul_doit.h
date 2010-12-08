#ifndef FACUL_DOIT_H
#define FACUL_DOIT_H

#include "modredc_ul.h"
#include "modredc_15ul.h"
#include "modredc_2ul2.h"
#include "facul.h"

int primetest_ul (const modulusredcul_t m);
int primetest_15ul (const modulusredc15ul_t m);
int primetest_2ul2 (const modulusredc2ul2_t m);
int facul_doit_ul (unsigned long *, const modulusredcul_t, 
		   const facul_strategy_t *, const int);
int facul_doit_15ul (unsigned long *, const modulusredc15ul_t, 
		     const facul_strategy_t *, const int);
int facul_doit_2ul2 (unsigned long *, const modulusredc2ul2_t, 
		     const facul_strategy_t *, const int);

#endif /* FACUL_DOIT_H */
