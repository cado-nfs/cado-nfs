#include "cado.h" // IWYU pragma: keep
#include "mpfq/mpfq_vbase.h"     // for mpfq_vbase_s
#include "parallelizing_info.h"  // for parallelizing_info_s, pi_comm, seria...
#include "xdotprod.h"

void x_dotprod(void * dst, uint32_t * xv, unsigned int m, unsigned int nx, mmt_vec_ptr v, int sign)
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
        void * where = v->abase->vec_subvec(v->abase, dst, j);
        for(unsigned int t = 0 ; t < nx ; t++) {
            uint32_t i = xv[j*nx+t];
            unsigned int vi0 = v->i0 + mmt_my_own_offset_in_items(v);
            unsigned int vi1 = vi0 + mmt_my_own_size_in_items(v);
            if (i < vi0 || i >= vi1)
                continue;
            void * coeff = v->abase->vec_subvec(v->abase, v->v, i - v->i0);
            if (sign > 0) {
                v->abase->add(v->abase, where, where, coeff);
            } else {
                v->abase->sub(v->abase, where, where, coeff);
            }
        }
    }
}

