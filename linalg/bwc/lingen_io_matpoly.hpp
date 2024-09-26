#ifndef LINGEN_IO_MATPOLY_HPP_
#define LINGEN_IO_MATPOLY_HPP_
// The selection is done within lingen_matpoly_select.hpp only
// IWYU pragma: no_include "mpfq_layer.h"
// IWYU pragma: no_include "mpfq_fake.hpp"
#include <iosfwd>
#include <cstdio>       // FILE
#include <vector>
#include "lingen_matpoly_select.hpp" // IWYU pragma: keep
struct cxx_param_list;

extern unsigned int io_matpoly_block_size;

void lingen_io_matpoly_decl_usage(cxx_param_list & pl);
void lingen_io_matpoly_lookup_parameters(cxx_param_list & pl);
void lingen_io_matpoly_interpret_parameters(cxx_param_list & pl);


int matpoly_read(matpoly::arith_hard * ab, FILE * f, matpoly & M, unsigned int k0, unsigned int k1, int ascii, int transpose);
int matpoly_write(matpoly::arith_hard * ab, std::ostream& os, matpoly const & M, unsigned int k0, unsigned int k1, int ascii, int transpose);
int matpoly_write_split(matpoly::arith_hard * ab, std::vector<std::ofstream> & fw, matpoly const & M, unsigned int k0, unsigned int k1, int ascii);

#endif	/* LINGEN_IO_MATPOLY_HPP_ */
