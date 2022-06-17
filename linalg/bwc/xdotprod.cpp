#include "cado.h" // IWYU pragma: keep
#include "arith-generic.hpp"
#include "parallelizing_info.hpp"  // for parallelizing_info_s, pi_comm, seria...
#include "xdotprod.hpp"

void x_dotprod(arith_generic::elt * dst, uint32_t * xv, unsigned int m, unsigned int nx, mmt_vec_ptr v, int sign)
{
    /* We're reading from the shared right vector data -- this area is
     * written to by the other threads in the column. Some of them might
     * be lingering in reduce operations, so we have to wait for them
     */
    if (!v->siblings) {
        serialize_threads(v->pi->wr[v->d]);
    } else {
        // I'd presume that no locking is needed here. But it's unchecked
        // ASSERT_ALWAYS(0);
    }

    for(unsigned int j = 0 ; j < m ; j++) {
        arith_generic::elt & where = v->abase->vec_item(dst, j);
        for(unsigned int t = 0 ; t < nx ; t++) {
            uint32_t i = xv[j*nx+t];
            unsigned int vi0 = v->i0 + mmt_my_own_offset_in_items(v);
            unsigned int vi1 = vi0 + mmt_my_own_size_in_items(v);
            if (i < vi0 || i >= vi1)
                continue;
            arith_generic::elt & coeff = v->abase->vec_item(v->v, i - v->i0);
            if (sign > 0) {
                v->abase->add_and_reduce(where, coeff);
            } else {
                v->abase->sub_and_reduce(where, coeff);
            }
        }
    }
}

