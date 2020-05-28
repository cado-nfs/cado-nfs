#include "cado.h" // IWYU pragma: keep
#include <cstring>
#include <cstdlib>                      // for free, malloc
#include "macros.h"                      // for ASSERT_ALWAYS
#include "bblas_mat64.hpp"  // for mat64
#include "bblas_perm_matrix.hpp"
#include "misc.h"      // cado_ctz64

/* perm_matrix stuff */
void perm_matrix_init(perm_matrix_ptr x, int n)
{
    x->n=n;
    x->v=(int*) malloc(n*sizeof(int));
    for(int i = 0 ; i < n ; i++) x->v[i]=i;
}
void perm_matrix_clear(perm_matrix_ptr x)
{
    free(x->v); x->v = NULL; x->n = 0;
}
void perm_matrix_transpose(perm_matrix_ptr x, perm_matrix_srcptr y)
{
    if (x == y) {
        perm_matrix t;
        perm_matrix_init(t, y->n);
        perm_matrix_transpose(x, t);
        perm_matrix_clear(t);
        return;
    }
    ASSERT_ALWAYS(x->n == y->n);
    for(int i = 0 ; i < x->n ; i++)
        x->v[i]=-1;
    for(int i = 0 ; i < x->n ; i++) {
        if (y->v[i] >= 0)
            x->v[y->v[i]]=i;
    }
}

void perm_matrix_get_matrix(mat64 * qm, perm_matrix_ptr qp)
{
    int * phi = qp->v;
    int n = qp->n;
    ASSERT_ALWAYS((n%64)==0);
    int nb = n/64;
    memset((void *) qm, 0, nb*nb*sizeof(mat64));
    uint64_t * qq = (uint64_t*) qm;
    for(int k = 0 ; k < n ; ) {
        for(int jq = 0 ; jq < 64 ; jq++, k++, qq++) {
            int v = phi[k];
            ASSERT_ALWAYS(v >= 0);
            qq[v&~63]=((uint64_t)1)<<(v%64);
        }
        qq+=(nb-1)*64;
    }
}


/* phi_p must point to a zone where at least as many limbs as the number
 * of set bits in the bits[] array
 */
void perm_matrixtab_complete(int * phi, uint64_t * bits, int nbits)
{
    ASSERT_ALWAYS(nbits % 64 == 0);
    for(int offset=0 ; offset < nbits ; offset+=64) {
        for(uint64_t w = bits[offset/64], z ; w ; w^=z) {
            z = w^(w&(w-1));
            uint64_t j = cado_ctz64(w);
            *phi++ = offset + j;
        }
    }
}

/* Given phi as in PLUQ64_inner below, return two permutation matrices
 * such that p*u*transpose(q) is a diagonal matrix.
 */
void pqperms_from_phi(perm_matrix_ptr p, perm_matrix_ptr q, int * phi, int m, int n)
{
    int ip = 0;ASSERT_ALWAYS((m%64)==0);ASSERT_ALWAYS(p->n==m);
    int jq = 0;ASSERT_ALWAYS((n%64)==0);ASSERT_ALWAYS(q->n==n);
    uint64_t ms[m/64]; for(int i = 0 ; i < m/64 ; i++) ms[i]=~((uint64_t)0);
    uint64_t ns[n/64]; for(int j = 0 ; j < n/64 ; j++) ns[j]=~((uint64_t)0);
    uint64_t w;
    for(int k = 0 ; k < m ; k++) {
        int v = phi[k];
        if (v < 0) continue;
        p->v[ip++] = k; w=((uint64_t)1)<<(k%64); ms[k/64]^=w;
        q->v[jq++] = v; w=((uint64_t)1)<<(v%64); ns[v/64]^=w;
    }
    perm_matrixtab_complete(p->v + ip, ms, m);
    perm_matrixtab_complete(q->v + jq, ns, n);
}


