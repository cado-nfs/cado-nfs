#ifndef LINGEN_CHECKPOINTS_HPP_
#define LINGEN_CHECKPOINTS_HPP_

#include <stddef.h>               // for size_t
#include <stdexcept>              // for runtime_error
#include <string>                 // for string
#include "lingen_bigmatpoly.hpp"  // for bigmatpoly
class matpoly;
struct bmstatus;
struct cxx_param_list;

enum cp_which {
    LINGEN_CHECKPOINT_E,
    LINGEN_CHECKPOINT_PI,
};

template<typename matpoly_type>
int load_checkpoint_file(bmstatus & bm, cp_which which, matpoly_type &, unsigned int t0, unsigned int t1);
template<typename matpoly_type>
int save_checkpoint_file(bmstatus & bm, cp_which which, matpoly_type const &, unsigned int t0, unsigned int t1);

template<>
int load_checkpoint_file<matpoly>(bmstatus & bm, cp_which which, matpoly &, unsigned int t0, unsigned int t1);
template<>
int save_checkpoint_file<matpoly>(bmstatus & bm, cp_which which, matpoly const &, unsigned int t0, unsigned int t1);

template<>
int load_checkpoint_file<bigmatpoly>(bmstatus & bm, cp_which which, bigmatpoly &, unsigned int t0, unsigned int t1);
template<>
int save_checkpoint_file<bigmatpoly>(bmstatus & bm, cp_which which, bigmatpoly const &, unsigned int t0, unsigned int t1);

struct lingen_checkpoint {
    bmstatus & bm;
    int level;
    unsigned int t0;
    unsigned int t1;
    unsigned int target_t;
    int mpi;
    int rank;
    std::string auxfile;
    std::string sdatafile;
    std::string gdatafile;
    std::string datafile; /* a copy of either sdatafile or gdatafile */
    /* be sure to change when needed */
    static constexpr unsigned long format = 4;
    static std::string get_cp_basename(bmstatus & bm, cp_which which, unsigned int t0, unsigned int t1);
    lingen_checkpoint(bmstatus & bm, unsigned int t0, unsigned int t1, int mpi, std::string base);
    lingen_checkpoint(bmstatus & bm, cp_which which, unsigned int t0, unsigned int t1, int mpi) : lingen_checkpoint(bm, t0, t1, mpi, get_cp_basename(bm, which, t0, t1)) {}
    bool save_aux_file(size_t Xsize) const;
    bool load_aux_file(size_t & Xsize);
    bool checkpoint_already_present() const;
    int load_data_file(matpoly & X);
    int save_data_file(matpoly const & X);
    struct invalid_aux_file : public std::runtime_error {
        invalid_aux_file(std::string const & s) : std::runtime_error(s) {}
    };
    static void decl_usage(cxx_param_list & pl);
    static void lookup_parameters(cxx_param_list & pl);
    static void interpret_parameters(cxx_param_list & pl);
    static const char * directory;
    static unsigned int threshold;
    static int save_gathered;
};

#endif	/* LINGEN_CHECKPOINTS_HPP_ */
