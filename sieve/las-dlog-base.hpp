#ifndef LAS_DLOG_BASE_HPP_
#define LAS_DLOG_BASE_HPP_

#include <vector>    // for vector
#include "cado_poly.h"   // cxx_cado_poly
#include "renumber.hpp"   // for renumber_t
#include "typedefs.h"
struct cxx_param_list;

struct las_dlog_base {
    cxx_cado_poly const & cpoly;
private:
    char * renumberfilename;
    char * logfilename;

    renumber_t renumber_table;
    std::vector<bool> known_logs;
    unsigned long lpb[2];

    void read();
public:
    bool is_known(int side, p_r_values_t p, p_r_values_t r) const;
    las_dlog_base(cxx_cado_poly const &, cxx_param_list & pl);
    ~las_dlog_base();
    static void declare_usage(cxx_param_list & pl);
    static void lookup_parameters(cxx_param_list & pl);

};

#endif	/* LAS_DLOG_BASE_HPP_ */
