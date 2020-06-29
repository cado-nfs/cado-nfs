#include "cado.h" // IWYU pragma: keep

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "matmul.h"
#include "macros.h"
#include "mpfq/mpfq.h"
#include "mpfq/mpfq_vbase.h"
#include "matmul-mf.h"
#include "portability.h" // asprintf // IWYU pragma: keep
#include "params.h"
#include "raw_matrix_u32.h"   // for matrix_u32


void declare_usage(param_list_ptr pl)
{
    param_list_decl_usage(pl, "matrix-file", "matrix file to work with");
    param_list_decl_usage(pl, "prime", "characteristic of the base field [default=2]");
    param_list_decl_usage(pl, "direction", "direction of the product, left for v*M, right for M*v [default=left for p=2, right otherwise]");
    param_list_decl_usage(pl, "withcoeffs", "whether we have coefficients in the matrix (i.e. not only 1's). Defaults to 1 (true) for p>2");
    param_list_decl_usage(pl, "impl", "name of the implementation backend. Defaults to bucket for p==2, basicp for p>2");
    param_list_decl_usage(pl, "groupsize", "number of vectors to consider together (defaults to 64 for p==2, 1 for p>2)");
    param_list_decl_usage(pl, "tmpdir", "directory where matrix cache file is saved (defaults to /tmp)\n");
}

int main(int argc, char * argv[])
{
    mpfq_vbase xx;
    matmul_t mm;
    char *argv0 = argv[0];

    mpz_t prime;
    int withcoeffs = 0;
    int direction = 0;
    int groupsize;
    mpz_init_set_ui(prime, 2);
    const char * impl = NULL;
    const char * matrixfile = NULL;
    const char * tmp = NULL;
    const char * tmpdir = "/tmp";

    setvbuf(stdout,NULL,_IONBF,0);
    setvbuf(stderr,NULL,_IONBF,0);

    param_list pl;
    param_list_init(pl);

    declare_usage(pl);


    for(argv++, argc-- ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    param_list_parse_mpz(pl, "prime", prime);
    direction = mpz_cmp_ui(prime, 2) != 0;      /* 0 == left */
    withcoeffs = mpz_cmp_ui(prime, 2) != 0;     /* 0 == no coeffs */
    groupsize = mpz_cmp_ui(prime, 2) == 0 ? 64 : 1;

    param_list_parse_int(pl, "withcoeffs", &withcoeffs);
    
    
    if ((tmp = param_list_lookup_string(pl, "direction")) != NULL) {
        if (strcmp(tmp, "left") == 0) {
            direction = 0;
        } else {
            direction = 1;
        }
    }

    if ((tmp = param_list_lookup_string(pl, "tmpdir")) != NULL) {
        tmpdir = tmp;
    }

    if ((impl = param_list_lookup_string(pl, "impl")) == NULL) {
        impl = mpz_cmp_ui(prime, 2) == 0 ? "bucket" : "basicp";
    }

    if ((matrixfile = param_list_lookup_string(pl, "matrix-file")) == NULL) {
        fprintf(stderr, "Error: argument matrix-file is mandatory\n");
        exit(EXIT_FAILURE);
    }
    param_list_warn_unused(pl);

    mpfq_vbase_oo_field_init_byfeatures(xx,
                MPFQ_PRIME_MPZ, prime,
                MPFQ_SIMD_GROUPSIZE, groupsize,
                MPFQ_DONE);

    if (direction == 1) {
        fprintf(stderr, "Saving cache for matrix-times-vector\n");
    } else {
        fprintf(stderr, "Saving cache for vector-times-matrix\n");
    }

    /* build a file name for the cache file */
    char * locfile;
    {
        char * matrixfile_copy = strdup(matrixfile);
        char * basename = strrchr(matrixfile_copy, '/');
        char * tmp;
        if (basename == NULL) {
            basename = matrixfile_copy;
        }
        if ((tmp = strstr(basename, ".bin")) != NULL) {
            *tmp='\0';
        }
        asprintf(&locfile, "%s/%s", tmpdir, basename);
        free(matrixfile_copy);
    }

    mm = matmul_init(xx, 0, 0, locfile, impl, pl, direction);

    /* uh ? */
    ASSERT_ALWAYS(mm->store_transposed == !direction);
    matrix_u32 m;
    memset(m, 0, sizeof(matrix_u32));
    m->mfile = matrixfile;
    /* The bfile here makes very little sense -- we're working with the
     * local file anyway */
    m->bfile = NULL;

    mf_prepare_matrix_u32(mm, m, matrixfile, withcoeffs);

    matmul_build_cache(mm, m);
    matmul_save_cache(mm);
    matmul_clear(mm);

    /* done here just because with have some lookup_string's into the pl
     * struct here and there.
     */
    param_list_clear(pl);
    free(locfile);
}
