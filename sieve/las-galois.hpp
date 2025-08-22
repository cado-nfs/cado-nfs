#ifndef CADO_LAS_GALOIS_HPP
#define CADO_LAS_GALOIS_HPP

#include <ostream>
#include <vector>

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "galois_action.hpp"
#include "relation.hpp"

void skip_galois_roots(cxx_mpz const & q, std::vector<cxx_mpz> & roots,
                       galois_action const & gal_action);

void add_relations_with_galois(const char *galois, std::ostream& os,
				      const char *comment, unsigned long *cpt,
				      relation &rel);

#endif	/* CADO_LAS_GALOIS_HPP */
