#ifndef CADO_CLIQUE_REMOVAL_HPP
#define CADO_CLIQUE_REMOVAL_HPP

/* A clique is a connected component of the graph where the nodes are the rows
   and the edges are the columns of weight 2.
   /!\ It is not a clique is the sense of graph theory.
   We will try to use the "connected component" terminology instead.
*/

#include <stdint.h>
#include "purge_matrix.hpp"

void comp_print_info_weight_function ();

/********************** main functions ****************************************/
void cliques_removal (purge_matrix_ptr, int64_t, unsigned int, int);

#endif /* CADO_CLIQUE_REMOVAL_HPP */

