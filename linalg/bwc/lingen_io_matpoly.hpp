#ifndef CADO_LINGEN_IO_MATPOLY_HPP
#define CADO_LINGEN_IO_MATPOLY_HPP

#include <cstdio>

#include <iosfwd>
#include <vector>

#include "lingen_matpoly_select.hpp"

struct cxx_param_list;

template<bool is_binary>
struct lingen_io_matpoly {
    static unsigned int block_size;
    
    static void decl_usage(cxx_param_list & pl);
    static void lookup_parameters(cxx_param_list & pl);
    static void interpret_parameters(cxx_param_list & pl);

    typedef matpoly<is_binary> matpoly_type;
    typedef typename matpoly_type::arith_hard arith_hard_t;

    static int read(arith_hard_t * ab, FILE * f, matpoly_type & M, unsigned int k0, unsigned int k1, int ascii, int transpose);
    static int write(arith_hard_t * ab, std::ostream& os, matpoly_type const & M, unsigned int k0, unsigned int k1, int ascii, int transpose);
    static int write_split(arith_hard_t * ab, std::vector<std::ofstream> & fw, matpoly_type const & M, unsigned int k0, unsigned int k1, int ascii);
};

/* both lingen_io_matpoly<true> and lingen_io_matpoly<false> have a
 * complete set of full specializations defined. It's admittedly not very
 * satisfactory.
 */
#ifdef LINGEN_BINARY
template<>
int lingen_io_matpoly<true>::write(arith_hard_t * ab, std::ostream& os, matpoly_type const & M, unsigned int k0, unsigned int k1, int ascii, int transpose);
template<>
int lingen_io_matpoly<true>::write_split(arith_hard_t * ab, std::vector<std::ofstream> & fw, matpoly_type const & M, unsigned int k0, unsigned int k1, int ascii);
template<>
int lingen_io_matpoly<true>::read(arith_hard_t * ab, FILE * f, matpoly_type & M, unsigned int k0, unsigned int k1, int ascii, int transpose);
extern template struct lingen_io_matpoly<true>;
#else
template<>
int lingen_io_matpoly<false>::write(arith_hard_t * ab, std::ostream& os, matpoly_type const & M, unsigned int k0, unsigned int k1, int ascii, int transpose);
template<>
int lingen_io_matpoly<false>::write_split(matpoly_type::arith_hard * ab, std::vector<std::ofstream> & fw, matpoly_type const & M, unsigned int k0, unsigned int k1, int ascii);
template<>
int lingen_io_matpoly<false>::read(matpoly_type::arith_hard * ab, FILE * f, matpoly_type & M, unsigned int k0, unsigned int k1, int ascii, int transpose);
extern template struct lingen_io_matpoly<false>;
#endif


#if 0
extern unsigned int io_matpoly_block_size;

void lingen_io_matpoly_decl_usage(cxx_param_list & pl);
void lingen_io_matpoly_lookup_parameters(cxx_param_list & pl);
void lingen_io_matpoly_interpret_parameters(cxx_param_list & pl);


template<bool is_binary>
int matpoly_read(typename matpoly_type::arith_hard * ab, FILE * f, matpoly_type & M, unsigned int k0, unsigned int k1, int ascii, int transpose);
template<bool is_binary>
int matpoly_write(typename matpoly_type::arith_hard * ab, std::ostream& os, matpoly_type const & M, unsigned int k0, unsigned int k1, int ascii, int transpose);
template<bool is_binary>
int matpoly_write_split(typename matpoly_type::arith_hard * ab, std::vector<std::ofstream> & fw, matpoly_type const & M, unsigned int k0, unsigned int k1, int ascii);

#endif

#endif	/* LINGEN_IO_MATPOLY_HPP_ */
