#ifndef CADO_MF_BAL_HPP
#define CADO_MF_BAL_HPP

#include <string>

#include "params.h"      // param_list_ptr

struct mf_bal_args {
    std::string rwfile;
    std::string cwfile;
    std::string mfile;
    std::string bfile;
    int quiet = 0;
    int nh = 0;
    int nv = 0;
    int withcoeffs = 0;
    int rectangular = 0;
    int skip_decorrelating_permutation = 0;
    enum { MF_BAL_PERM_YES, MF_BAL_PERM_NO, MF_BAL_PERM_AUTO } do_perm[2] { MF_BAL_PERM_AUTO, MF_BAL_PERM_AUTO };
};

void mf_bal(struct mf_bal_args * mba);
void mf_bal_adjust_from_option_string(struct mf_bal_args * mba, const char * opts);
void mf_bal_decl_usage(cxx_param_list& pl);
void mf_bal_configure_switches(cxx_param_list & pl, struct mf_bal_args * mba);
void mf_bal_parse_cmdline(struct mf_bal_args * mba, cxx_param_list & pl, int * p_argc, char const *** p_argv);
void mf_bal_interpret_parameters(struct mf_bal_args * mba, cxx_param_list & pl);

#endif	/* MF_BAL_HPP_ */
