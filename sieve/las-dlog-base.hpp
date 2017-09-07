#ifndef LAS_DLOG_BASE_HPP_
#define LAS_DLOG_BASE_HPP_

#include "utils.h"
#include "params.h"

#include <vector>

struct las_dlog_base {
private:
    char * renumberfilename;
    char * logfilename;

    renumber_t renumber_table;
    std::vector<bool> known_logs;
    unsigned long lpb[2];

    void lookup_parameters(param_list pl);
    void read();
public:
    bool is_known(int side, uint64_t p, uint64_t r) const;
    las_dlog_base(param_list_ptr pl);
    ~las_dlog_base();

    static void declare_parameter_usage(param_list pl);
};

#endif	/* LAS_DLOG_BASE_HPP_ */
