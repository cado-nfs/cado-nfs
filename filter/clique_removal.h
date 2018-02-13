#ifndef CLIQUE_REMOVAL_H_
#define CLIQUE_REMOVAL_H_

/* A clique is a connected component of the graph where the nodes are the rows
   and the edges are the columns of weight 2.
   /!\ It is not a clique is the sense of graph theory.
   We will try to use the "connected component" terminology instead.
*/

#ifdef __cplusplus
extern "C" {
#endif

void comp_print_info_weight_function ();

/********************** main functions ****************************************/
void cliques_removal (purge_matrix_ptr, int64_t, unsigned int, int);

#ifdef __cplusplus
}
#endif

#endif /* CLIQUE_REMOVAL_H_ */

