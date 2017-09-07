#ifndef MF_BAL_H_
#define MF_BAL_H_

#ifdef __cplusplus
extern "C" {
#endif

struct mf_bal_args {
    const char * rwfile;
    const char * cwfile;
    const char * mfile;
    const char * bfile;
    int quiet;
    int nh;
    int nv;
    int withcoeffs;
    int rectangular;
    int skip_decorrelating_permutation;
    enum { MF_BAL_PERM_YES, MF_BAL_PERM_NO, MF_BAL_PERM_AUTO } do_perm[2];
};

void mf_bal(struct mf_bal_args * mba);
void mf_bal_adjust_from_option_string(struct mf_bal_args * mba, const char * opts);
void mf_bal_decl_usage(param_list_ptr pl);
void mf_bal_configure_switches(param_list_ptr pl, struct mf_bal_args * mba);
void mf_bal_parse_cmdline(struct mf_bal_args * mba, param_list_ptr pl, int * p_argc, char *** p_argv);
void mf_bal_interpret_parameters(struct mf_bal_args * mba, param_list_ptr pl);

#ifdef __cplusplus
}
#endif

#endif	/* MF_BAL_H_ */
