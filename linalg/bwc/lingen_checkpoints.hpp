#ifndef CADO_LINGEN_CHECKPOINTS_HPP
#define CADO_LINGEN_CHECKPOINTS_HPP

#include <cstddef>

#include <stdexcept>
#include <vector>
#include <istream>
#include <string>

#include "cxx_mpz.hpp"
#include "lingen_bigmatpoly.hpp"

template<bool is_binary> class matpoly;
template<bool is_binary> struct bmstatus;
struct cxx_param_list;

enum cp_which {
    LINGEN_CHECKPOINT_E,
    LINGEN_CHECKPOINT_PI,
};

template<bool is_binary>
int load_checkpoint_file(bmstatus<is_binary> & bm, cp_which which, matpoly<is_binary> &, unsigned int t0, unsigned int t1);
template<bool is_binary>
int save_checkpoint_file(bmstatus<is_binary> & bm, cp_which which, matpoly<is_binary> const &, unsigned int t0, unsigned int t1);

template<bool is_binary>
int load_checkpoint_file(bmstatus<is_binary> & bm, cp_which which, bigmatpoly<is_binary> &, unsigned int t0, unsigned int t1);
template<bool is_binary>
int save_checkpoint_file(bmstatus<is_binary> & bm, cp_which which, bigmatpoly<is_binary> const &, unsigned int t0, unsigned int t1);

template<bool is_binary_arg>
struct lingen_checkpoint {
    static constexpr bool is_binary = is_binary_arg;
    bmstatus<is_binary> & bm;
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
    static std::string get_cp_basename(bmstatus<is_binary> & bm, cp_which which, unsigned
int t0, unsigned int t1);
    lingen_checkpoint(bmstatus<is_binary> & bm, unsigned int t0, unsigned int t1, int mpi, std::string base);
    lingen_checkpoint(bmstatus<is_binary> & bm, cp_which which, unsigned int t0, unsigned int t1, int mpi) : lingen_checkpoint(bm, t0, t1, mpi, get_cp_basename(bm, which, t0, t1)) {}
    bool save_aux_file(size_t Xsize) const;
    bool load_aux_file(size_t & Xsize);
    bool checkpoint_already_present() const;
    int load_data_file(matpoly<is_binary> & X);
    int save_data_file(matpoly<is_binary> const & X);
    struct invalid_aux_file : public std::runtime_error {
        explicit invalid_aux_file(std::string const & s)
            : std::runtime_error(s)
        {}
    };
    static void decl_usage(cxx_param_list & pl);
    static void lookup_parameters(cxx_param_list & pl);
    static void interpret_parameters(cxx_param_list & pl);
    static std::string default_directory;
    static unsigned int threshold;
    static int save_gathered;

    struct header_info
    {
        static constexpr bool is_binary = is_binary_arg;
        unsigned int m = 0;
        unsigned int n = 0;
        int level = 0;
        unsigned int t0 = 0;
        unsigned int t1 = 0;
        unsigned int t = 0;
        unsigned int ncoeffs = 0;
        cxx_mpz p;
        std::vector<unsigned int> delta;
        std::vector<int> lucky;
        int done = 0;

        std::istream& input(std::istream&);
    };
};

#ifdef LINGEN_BINARY
inline std::istream& operator>>(std::istream & is, lingen_checkpoint<true>::header_info & u) { return u.input(is); }
#else
inline std::istream& operator>>(std::istream & is, lingen_checkpoint<false>::header_info & u) { return u.input(is); }
#endif
#ifdef LINGEN_BINARY
extern template struct lingen_checkpoint<true>;
#else
extern template struct lingen_checkpoint<false>;
#endif

#endif	/* LINGEN_CHECKPOINTS_HPP_ */
