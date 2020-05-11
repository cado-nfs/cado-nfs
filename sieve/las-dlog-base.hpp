#ifndef LAS_DLOG_BASE_HPP_
#define LAS_DLOG_BASE_HPP_

#include <cstdint>  // for uint64_t
#include <vector>    // for vector
#include "utils.h"   // for renumber_t

struct las_dlog_base {

private:
    char * renumberfilename;
    char * logfilename;

    renumber_t renumber_table;
    std::vector<bool> known_logs;
    unsigned long lpb[2];

    void read();
public:
    bool is_known(int side, uint64_t p, uint64_t r) const;
    las_dlog_base(cxx_param_list & pl);
    ~las_dlog_base();
    static void declare_usage(cxx_param_list & pl);
    static void lookup_parameters(cxx_param_list & pl);

};

#endif	/* LAS_DLOG_BASE_HPP_ */
