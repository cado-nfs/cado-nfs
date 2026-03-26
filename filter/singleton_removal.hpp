#ifndef CADO_SINGLETON_REMOVAL_HPP
#define CADO_SINGLETON_REMOVAL_HPP

#include "purge_matrix.hpp"
#include <cstdint>

void singleton_removal_oneiter_mono (purge_matrix_ptr mat);
void singleton_removal_oneiter_mt (purge_matrix_ptr, unsigned int);
int64_t singleton_removal (purge_matrix_ptr, unsigned int, int);

#endif /* CADO_SINGLETON_REMOVAL_HPP */
