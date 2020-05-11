#ifndef LAS_GALOIS_HPP_
#define LAS_GALOIS_HPP_

#include <ostream>

#include <gmp.h>    // for mpz_t
struct relation;

int skip_galois_roots(const int orig_nroots, const mpz_t q, mpz_t *roots,
		  const char *galois_autom);

void add_relations_with_galois(const char *galois, std::ostream& os,
				      const char *comment, unsigned long *cpt,
				      relation &rel);

#endif	/* LAS_GALOIS_HPP_ */
