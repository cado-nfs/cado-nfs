#include "cado.h" // IWYU pragma: keep

#include <cstdlib> // for abort

#include "modredc126.hpp"

#include "arithxx_common.hpp"

/* Only the .cpp source files that emit the non-inline symbols will
 * include this impl header file. So even though it does not look like
 * we're using it, in fact we are!  */
#include "arithxx_api_impl.hpp"      // IWYU pragma: keep
#include "arithxx_api128_impl.hpp"  // IWYU pragma: keep
#include "arithxx_redc_impl.hpp"  // IWYU pragma: keep
#include "arithxx_redc128_impl.hpp"  // IWYU pragma: keep
#include "arithxx_batch_Q_to_Fp_impl.hpp"  // IWYU pragma: keep

#include "arithxx_redc128.hpp"


// scan-headers: stop here

template struct arithxx_details::batch_Q_to_Fp_context<arithxx_modredc126>;
template struct arithxx_details::redc<arithxx_modredc126>;
template struct arithxx_details::redc128<arithxx_modredc126>;
template struct arithxx_details::api<arithxx_modredc126>;
template struct arithxx_details::api_bysize<arithxx_modredc126>;
