#include "cado.h"

#include <algorithm>

#include "lingen_bmstatus.hpp"

lingen_call_companion & bmstatus::companion(int depth, size_t L)/*{{{*/
{
    lingen_hints::key_type K { depth, L };

    if (hints.find(K) != hints.end())
        return hints[K];

    fprintf(stderr, "# No tuned configuration for"
            " depth=%d input_length=%zu\n",
            depth, L);

    struct same_depth {
        int d;
        bool operator()(lingen_hints::value_type const & kv) const {
            return kv.first.depth == d;
        }
    };

    auto it = std::find_if(hints.begin(), hints.end(), same_depth { depth } );
    if (it == hints.end()) {
        fprintf(stderr, "# No tuned configuration for"
                " depth=%d !!!\n", depth);
        exit(EXIT_FAILURE);
    }
    fprintf(stderr, "# Using nearby configuration for"
            " depth=%d input_length=%zu\n",
            it->first.depth, it->first.L);

    hints[K]=it->second;

    return it->second;
}/*}}}*/


