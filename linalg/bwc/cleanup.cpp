#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <gmp.h>
#include "blockmatrix.hpp"
#include "bblas_gauss.h"

#include "params.h"     // param_list
#include "macros.h"

void usage()
{
    fprintf(stderr, "Usage: ./cleanup -ncols <N> -out <file> <file.0> <file.1> ...\n");
}
// coverity[root_function]
int main(int argc, char **argv)
{
    param_list pl;
    param_list_init(pl);
    /* Vectors must be *in order* !!! */
    argv++,argc--;
    unsigned int ncols = 0;
    const char * outfile = NULL;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        /* might also be a kernel file */
        break;
    }
    outfile = param_list_lookup_string(pl, "out");
    param_list_parse_uint(pl, "ncols", &ncols);
    if (param_list_warn_unused(pl) || ncols == 0 || outfile == 0) {
        usage();
        fflush (stderr);
        exit(EXIT_FAILURE);
    }

        ASSERT_ALWAYS(ncols % 64 == 0);

    blockmatrix S(ncols, ncols);
    blockmatrix ST(ncols, ncols);
    blockmatrix T(ncols, ncols);
    S.set_identity();
    uint64_t * kzone = (uint64_t *) malloc(FLAT_BYTES_WITH_READAHEAD(ncols, ncols));
    int limbs_per_row = iceildiv(ncols, 64);

    unsigned int common_nrows = 0;
    unsigned int nrows;


    if (argc == 0) {
        fprintf(stderr, "Error: no input files\n");
        usage();
        exit(EXIT_FAILURE);
    }

    /* read very first file, deduce some important info about sizes */
    {
        struct stat sbuf[1];
        int rc = stat(argv[0], sbuf);
        if (rc < 0) { perror(argv[0]); exit(EXIT_FAILURE); }
        ASSERT_ALWAYS(sbuf->st_size % (ncols/8) == 0);
        nrows = sbuf->st_size / (ncols/8);
        common_nrows = nrows;
    }

    blockmatrix k(nrows, ncols);
    blockmatrix kprev(nrows, ncols);
    blockmatrix kfinal(nrows, ncols);
    kfinal.set_zero();
    blockmatrix kS(nrows, ncols);
    uint64_t * zone = (uint64_t *) malloc(FLAT_BYTES_WITH_READAHEAD(ncols, nrows));
    int limbs_per_col = iceildiv(nrows, 64);
    int prevrank = ncols;
    int rank0 = 0;


    for(int i = 0 ; i < argc ; i++) {
        struct stat sbuf[1];
        int rc = stat(argv[i], sbuf);
        if (rc < 0) { perror(argv[i]); exit(EXIT_FAILURE); }
        ASSERT_ALWAYS(sbuf->st_size % (ncols/8) == 0);
        nrows = sbuf->st_size / (ncols/8);
        fprintf(stderr, "%s: %u x %u\n", argv[i], nrows, ncols);
        ASSERT_ALWAYS(common_nrows == nrows);

        k.read_from_flat_file(0, 0, argv[i], nrows, ncols);

        /* we would like to have an in-place multiply. Trivial to do for
         * ncols==64, harder to get it right as well for >1 column
         * blocks. Therefore, we stick to simple and stupid code.
         */
        blockmatrix::mul_smallb(k, S);

        blockmatrix::copy_transpose_to_flat(zone, limbs_per_col, k);
        // coverity[swapped_arguments]
        int rank = spanned_basis(
                (mp_limb_t *) kzone,
                (mp_limb_t *) zone,
                ncols,
                nrows,
                sizeof(uint64_t) / sizeof(mp_limb_t) * limbs_per_col,
                sizeof(uint64_t) / sizeof(mp_limb_t) * limbs_per_row,
                NULL);
        // kzone*transpose(kS) is reduced
        // kS*transpose(kzone) is reduced (equivalent formulation)
        blockmatrix::copy_transpose_from_flat(T, kzone, limbs_per_row);
        // blockmatrix_reverse_columns(T, T);

        /* multiply kprev, k, and S by T */
        /* same comment as above applies, btw. */
        blockmatrix::mul_smallb(S, T);
        blockmatrix::mul_smallb(k, T);
        /* only columns [0..rank-1] in T are non-zero */

        if (i) {
            blockmatrix::mul_smallb(kprev, T);
        }

        if (i) {
            ASSERT_ALWAYS(rank <= prevrank);
            if (rank < prevrank) {
                printf("%s: rank drops from %d to %d\n", argv[i], prevrank, rank);
                /* bits [rank..prevrank[ of kprev are kernel vectors.
                 * They can be added to bits [rank..prevrank[ of kfinal
                 */
                kfinal.copy_colrange(kprev, rank, prevrank);
            }
        } else {
            printf("%s: rank %d\n", argv[i], rank);
            rank0 = rank;
        }

        //printf("%s: rank %d\n", argv[i], r);

        prevrank = rank;
        k.swap(kprev);
    }
    // finish assuming rank 0.
    kfinal.copy_colrange(kprev, 0, prevrank);

    /* Oh, now we need to check the combined rank of all these */
    blockmatrix::copy_transpose_to_flat(zone, limbs_per_col, kfinal);
    int rankf = spanned_basis(
            (mp_limb_t *) kzone,
            (mp_limb_t *) zone,
            ncols,
            nrows,
            sizeof(uint64_t) / sizeof(mp_limb_t) * limbs_per_col,
            sizeof(uint64_t) / sizeof(mp_limb_t) * limbs_per_row,
            NULL);
    blockmatrix::copy_transpose_from_flat(T, kzone, limbs_per_row);
    blockmatrix::mul_smallb(kfinal, T);
    if (rankf < rank0) {
        printf("final adjustment: rank drops from %d to %d\n", rank0, rankf);
    }

    kfinal.write_to_flat_file(outfile, 0, 0, common_nrows, ncols);
    printf("%s: written %d kernel vectors\n", outfile, rankf);
    free(zone);
    free(kzone);
    param_list_clear(pl);

    return 0;
}

// N:=20;i:=1;part:=RandomPartition(N);M:=Matrix(GF(2),N,N,[]);for s in part do for k in [1..s-1] do M[i+k-1,i+k]:=1; end for; i+:=s; end for;M:=Transpose(M);k:=Random([1..N]);V:=VectorSpace(GF(2),N);vs:=[Random(V):i in [1..k]]; k-Dimension(sub<V|[v*M^i:i in [0..N],v in vs]> meet Nullspace(M));
