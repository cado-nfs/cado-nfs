#include "cado.h"
#include "params.h"
#include <string>
#include <vector>
#include "portability.h"

/* This is just a helper function that is used form matmul_top.c
 *
 * We do this in order to avoid a super-weird interface in params.h
 *
 */

/* returns an allocated string holding the name of the midx-th submatrix */
char * matrix_list_get_item(param_list_ptr pl, const char * key, int midx)
{
    std::vector<std::string> m;
    int rc = param_list_parse(pl, key, m);
    if (rc == 0)
        return NULL;
    ASSERT_ALWAYS(midx < (int) m.size());
    return strdup(m[midx].c_str());
}

