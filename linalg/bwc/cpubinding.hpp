#ifndef CADO_CPUBINDING_HPP
#define CADO_CPUBINDING_HPP

#include "params.h"     // param_list

void cpubinding_decl_usage(cxx_param_list &);
void cpubinding_lookup_parameters(cxx_param_list & pl);

/* This returns an opaque pointer to data which will be used to perform
 * the actual cpu binding. This function must be called in
 * single-threaded context.
 * 
 * This returns NULL if cpubinding failed.
 *
 * If messages is not NULL, it is set to point to a newly allocated
 * string indicating all messages from the cpubinding engine. This is
 * meant to collect messages for various nodes in an MPI context, and
 * print only the unique ones (see parallelizing_info.c)
 */
void * cpubinding_get_info(char ** messages, cxx_param_list & pl, unsigned int, unsigned int);

/* perform the actual pinning. This must be called for each thread */
void cpubinding_do_pinning(void * pinning_info, int i, int j);

/* free the opaque pointer */
void cpubinding_free_info(void * pinning_info, unsigned int, unsigned int);

#endif	/* CPUBINDING_HPP_ */
