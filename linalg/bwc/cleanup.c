#include "cado.h" // IWYU pragma: keep

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "blockmatrix.h"
#include "gauss.h"
#include "params.h"     // param_list
#include "macros.h"

void usage()
{
    fprintf(stderr, "Usage: ./cleanup -ncols <N> -out <file> <file.0> <file.1> ...\n");
}
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

    blockmatrix S = NULL;
    blockmatrix T = NULL;
    blockmatrix ST = NULL;
    S = blockmatrix_alloc(ncols, ncols);
    ST = blockmatrix_alloc(ncols, ncols);
    T = blockmatrix_alloc(ncols, ncols);
    blockmatrix_set_identity(S);
    uint64_t * kzone = malloc(FLAT_BYTES_WITH_READAHEAD(ncols, ncols));
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

    blockmatrix k = blockmatrix_alloc(nrows, ncols);
    blockmatrix kprev = blockmatrix_alloc(nrows, ncols);
    blockmatrix kfinal = blockmatrix_alloc(nrows, ncols);
    blockmatrix_set_zero(kfinal);
    blockmatrix kS = blockmatrix_alloc(nrows, ncols);
    uint64_t * zone = malloc(FLAT_BYTES_WITH_READAHEAD(ncols, nrows));
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

        blockmatrix_read_from_flat_file(k, 0, 0, argv[i], nrows, ncols);

        /* we would like to have an in-place multiply. Trivial to do for
         * ncols==64, harder to get it right as well for >1 column
         * blocks. Therefore, we stick to simple and stupid code.
         */
        blockmatrix_mul_smallb(k, k, S);

        blockmatrix_copy_transpose_to_flat(zone, limbs_per_col, 0, 0, k);
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
        blockmatrix_copy_transpose_from_flat(T, kzone, limbs_per_row, 0, 0);
        // blockmatrix_reverse_columns(T, T);

        /* multiply kprev, k, and S by T */
        /* same comment as above applies, btw. */
        blockmatrix_mul_smallb(S, S, T);
        blockmatrix_mul_smallb(k, k, T);
        /* only columns [0..rank-1] in T are non-zero */

        if (i) {
            blockmatrix_mul_smallb(kprev, kprev, T);
        }

        if (i) {
            ASSERT_ALWAYS(rank <= prevrank);
            if (rank < prevrank) {
                printf("%s: rank drops from %d to %d\n", argv[i], prevrank, rank);
                /* bits [rank..prevrank[ of kprev are kernel vectors.
                 * They can be added to bits [rank..prevrank[ of kfinal
                 */
                blockmatrix_copy_colrange(kfinal, kprev, rank, prevrank);
            }
        } else {
            printf("%s: rank %d\n", argv[i], rank);
            rank0 = rank;
        }

        //printf("%s: rank %d\n", argv[i], r);

        prevrank = rank;
        blockmatrix_swap(kprev, k);
    }
    // finish assuming rank 0.
    blockmatrix_copy_colrange(kfinal, kprev, 0, prevrank);

    /* Oh, now we need to check the combined rank of all these */
    blockmatrix_copy_transpose_to_flat(zone, limbs_per_col, 0, 0, kfinal);
    int rankf = spanned_basis(
            (mp_limb_t *) kzone,
            (mp_limb_t *) zone,
            ncols,
            nrows,
            sizeof(uint64_t) / sizeof(mp_limb_t) * limbs_per_col,
            sizeof(uint64_t) / sizeof(mp_limb_t) * limbs_per_row,
            NULL);
    blockmatrix_copy_transpose_from_flat(T, kzone, limbs_per_row, 0, 0);
    blockmatrix_mul_smallb(kfinal, kfinal, T);
    if (rankf < rank0) {
        printf("final adjustment: rank drops from %d to %d\n", rank0, rankf);
    }

    blockmatrix_write_to_flat_file(outfile, kfinal, 0, 0, common_nrows, ncols);
    printf("%s: written %d kernel vectors\n", outfile, rankf);
    blockmatrix_free(k);
    blockmatrix_free(kprev);
    blockmatrix_free(kfinal);
    blockmatrix_free(kS);
    blockmatrix_free(S);
    blockmatrix_free(ST);
    blockmatrix_free(T);
    free(zone);
    free(kzone);
    param_list_clear(pl);

    return 0;
}

// N:=20;i:=1;part:=RandomPartition(N);M:=Matrix(GF(2),N,N,[]);for s in part do for k in [1..s-1] do M[i+k-1,i+k]:=1; end for; i+:=s; end for;M:=Transpose(M);k:=Random([1..N]);V:=VectorSpace(GF(2),N);vs:=[Random(V):i in [1..k]]; k-Dimension(sub<V|[v*M^i:i in [0..N],v in vs]> meet Nullspace(M));
