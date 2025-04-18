#ifndef CADO_CLIQUE_REMOVAL_H
#define CADO_CLIQUE_REMOVAL_H

/* A clique is a connected component of the graph where the nodes are the rows
   and the edges are the columns of weight 2.
   /!\ It is not a clique is the sense of graph theory.
   We will try to use the "connected component" terminology instead.
*/

#include <stdint.h>
#include "purge_matrix.h"

#ifdef __cplusplus
extern "C" {
#endif

void comp_print_info_weight_function ();

/********************** main functions ****************************************/
void cliques_removal (purge_matrix_ptr, int64_t, unsigned int, int);

#ifdef __cplusplus
}
#endif

#endif /* CADO_CLIQUE_REMOVAL_H */

