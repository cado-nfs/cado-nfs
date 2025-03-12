#include "cado.h" // IWYU pragma: keep
#include <cstdint>

#include <vector>

#include "arith-generic.hpp"
#include "matmul_top_vec.hpp"
#include "parallelizing_info.hpp" // for parallelizing_info_s, pi_comm, seria...
#include "xdotprod.hpp"

void x_dotprod(arith_generic::elt * dst, std::vector<uint32_t> const & xv,
               unsigned int j0, unsigned int j1, unsigned int nx,
               mmt_vec const & v, int sign)
{
    /* We're reading from the shared right vector data -- this area is
     * written to by the other threads in the column. Some of them might
     * be lingering in reduce operations, so we have to wait for them
     */
    if (mmt_vec_is_shared(v)) {
        serialize_threads(v.pi->wr[v.d]);
    } else {
        // I'd presume that no locking is needed here. But it's unchecked
        // ASSERT_ALWAYS(0);
    }

    for (unsigned int j = j0; j < j1; j++) {
        arith_generic::elt & where = v.abase->vec_item(dst, j - j0);
        for (unsigned int t = 0; t < nx; t++) {
            uint32_t const i = xv[j * nx + t];
            unsigned int const vi0 = v.i0 + mmt_my_own_offset_in_items(v);
            unsigned int const vi1 = vi0 + mmt_my_own_size_in_items(v);
            if (i < vi0 || i >= vi1)
                continue;
            arith_generic::elt const & coeff = v.abase->vec_item(v.v, i - v.i0);
            if (sign > 0) {
                v.abase->add_and_reduce(where, coeff);
            } else {
                v.abase->sub_and_reduce(where, coeff);
            }
        }
    }
}
