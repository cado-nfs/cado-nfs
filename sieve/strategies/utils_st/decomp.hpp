#ifndef DECOMP_HPP
#define DECOMP_HPP

#include <cstdio>

typedef struct decomp {
    unsigned int *tab;
    unsigned int len;		//length of the decomposition
    double nb_elem;		//number of elements which satisfy
				//this decomposition
} decomp_t;

decomp_t *decomp_create(unsigned int len, const unsigned int *tab, double nb_elem);

void decomp_free(decomp_t * t);

double decomp_get_nb_elem(decomp_t * t);

const unsigned int *decomp_get_tab(decomp_t * t);

unsigned int decomp_get_len(decomp_t * t);

void decomp_set_decomp(decomp_t * t, const unsigned int *tab, unsigned int len);

void decomp_set_nb_elem(decomp_t * t, double nb_elem);

int decomp_fprint(FILE * output_file, decomp_t * t);

int decomp_print(decomp_t * t);

#endif				/* DECOMP_HPP */
