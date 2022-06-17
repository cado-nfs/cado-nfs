#include "cado.h" // IWYU pragma: keep
#include "bwc_config.h" // BUILD_DYNAMICALLY_LINKABLE_BWC // IWYU pragma: keep
#include "matmul_facade.hpp"

void MATMUL_NAME(rebind)(matmul_ptr mm)
{
    REBIND_ALL(mm);
}

#ifdef  BUILD_DYNAMICALLY_LINKABLE_BWC
extern "C" void matmul_solib_do_rebinding(matmul_ptr mm)
{
    REBIND_ALL(mm);
}
#endif

