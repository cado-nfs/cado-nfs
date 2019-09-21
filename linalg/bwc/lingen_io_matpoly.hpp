#ifndef LINGEN_IO_MATPOLY_HPP_
#define LINGEN_IO_MATPOLY_HPP_

#include "lingen.hpp"
#include <ostream>
#include <fstream>
#include <cstdio>
#include <vector>
#include "params.h"

#ifdef SELECT_MPFQ_LAYER_u64k1
#include "lingen_matpoly_binary.hpp"
#include "lingen_matpoly_ft.hpp"
#else
/* lingen-matpoly is the default code. */
#include "lingen_matpoly.hpp"
#if 0
#ifdef ENABLE_MPI_LINGEN        /* in lingen.hpp */
#include "lingen_bigmatpoly.hpp"
#include "lingen_bigmatpoly_ft.hpp"
#endif
#endif
#endif

extern unsigned int io_matpoly_block_size;

void lingen_io_matpoly_decl_usage(cxx_param_list & pl);
void lingen_io_matpoly_lookup_parameters(cxx_param_list & pl);
void lingen_io_matpoly_interpret_parameters(cxx_param_list & pl);


int matpoly_read(abdst_field ab, FILE * f, matpoly & M, unsigned int k0, unsigned int k1, int ascii, int transpose);
int matpoly_write(abdst_field ab, std::ostream& os, matpoly const & M, unsigned int k0, unsigned int k1, int ascii, int transpose);
int matpoly_write_split(abdst_field ab, std::vector<std::ofstream> & fw, matpoly const & M, unsigned int k0, unsigned int k1, int ascii);

#endif	/* LINGEN_IO_MATPOLY_HPP_ */
