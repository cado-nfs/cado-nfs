#include "cado.h" // IWYU pragma: keep
#include "bwc_config.h" // BUILD_DYNAMICALLY_LINKABLE_BWC // IWYU pragma: keep
#include "matmul_facade.hpp"

extern matmul_interface_ctor_t CADO_CONCATENATE4(new_matmul_, ARITH_LAYER, _, MM_IMPL);

#ifdef  BUILD_DYNAMICALLY_LINKABLE_BWC
extern "C" matmul_interface_ctor_t * matmul_solib_reach_ctor()
{
    return &CADO_CONCATENATE4(new_matmul_, ARITH_LAYER, _, MM_IMPL);
}
#endif
