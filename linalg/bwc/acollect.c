#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/stat.h>
#include <unistd.h>
#include <limits.h>
#include <dirent.h>
#include <gmp.h>
#include "macros.h"
#include "bw-common.h"
#include "portability.h" // asprintf // IWYU pragma: keep
#include "params.h"

/* This program is rather standalone. It checks the current directory
 * for files matching the pattern A%u-%u.%u-%u, and concatenates them.
 *
 * It needs to know m, though (auto-detection would be easy).
 */

void usage()
{
    fprintf(stderr, "Usage: acollect [m=<m>] [wdir=<path>] [--remove-old]\n");
    exit(EXIT_FAILURE);
}

int remove_old = 0;
int bits_per_coeff = 1;

struct afile_s {
    unsigned int n0, n1, j0, j1;
};
typedef struct afile_s afile[1];
typedef struct afile_s * afile_ptr;


typedef int (*sortfunc_t) (const void *, const void *);

int afile_cmp(afile_ptr a, afile_ptr b)
{
    int dj0 = a->j0 - b->j0;
    int dj1 = a->j1 - b->j1;
    int dn0 = a->n0 - b->n0;
    int dn1 = a->n1 - b->n1;

    if (dn0) return dn0;
    if (dn1) return dn1;
    if (dj0) return dj0;
    if (dj1) return dj1;
    return 0;
}

struct afile_list {
    afile * a;
    unsigned int n;
    unsigned int alloc;
};

int read_afiles(struct afile_list * a)
{
    a->n = 0;
    DIR * dir = opendir(".");
    struct dirent * de;

    for( ; (de = readdir(dir)) != NULL ; ) {
        if (a->n >= a->alloc) {
            a->alloc += 32 + a->alloc / 4;
            a->a = realloc(a->a, a->alloc * sizeof(afile));
        }
        afile_ptr A = a->a[a->n];
        int k;
        int rc = sscanf(de->d_name, "A%u-%u.%u-%u%n", &A->n0, &A->n1, &A->j0, &A->j1, &k);
        /* rc is expected to be 4 or 5 depending on our reading of the
         * standard */
        if (rc < 4 || k != (int) strlen(de->d_name)) {
            // fprintf(stderr, "skipped %s\n", de->d_name);
            continue;
        }
        if ((A->n1 * bits_per_coeff) % CHAR_BIT || (A->n0 * bits_per_coeff) % CHAR_BIT) {
            fprintf(stderr, "%s has bad boundaries\n",
                    de->d_name);
            exit(EXIT_FAILURE);
        }
        struct stat sbuf[1];
        rc = stat(de->d_name, sbuf);
        if (rc < 0) {
            fprintf(stderr, "stat(%s): %s\n", de->d_name, strerror(errno));
            exit(EXIT_FAILURE);
        }
        ssize_t expected = bits_per_coeff;
        expected *= A->n1-A->n0;
        expected /= CHAR_BIT;
        expected *= A->j1-A->j0;
        expected *= bw->m;

        if (sbuf->st_size != expected) {
            fprintf(stderr, "%s does not have expected size %zu\n",
                    de->d_name, expected);
            exit(EXIT_FAILURE);
        }

        a->n++;
    }
    closedir(dir);

    if (a->n == 0) {
        return 0;
    }

    if (a->n == 1) {
        char * tmp;
        int rc = asprintf(&tmp, "A%u-%u.%u-%u",
                a->a[0]->n0,a->a[0]->n1,a->a[0]->j0,a->a[0]->j1);
        ASSERT_ALWAYS(rc >= 0);
        printf("%s\n", tmp);
        free(tmp);
        free(a->a);
        return 0;
    }

    qsort(a->a, a->n, sizeof(afile), (sortfunc_t) &afile_cmp);

    return a->n;
}

int main(int argc, char * argv[])
{
    param_list pl;

    bw_common_init(bw, &argc, &argv);
    param_list_init(pl);

    bw_common_decl_usage(pl);
    /* {{{ declare local parameters and switches */
    param_list_decl_usage(pl, "remove-old",
            "discard original A file once the concatenated file has been successfully written");
    param_list_configure_switch(pl, "--remove-old", &remove_old);
    /* }}} */

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);

    /* {{{ interpret our parameters */
    if (mpz_cmp_ui(bw->p, 2) > 0) {
        bits_per_coeff = ULONG_BITS * iceildiv(mpz_sizeinbase(bw->p, 2), ULONG_BITS);
    } else {
        bits_per_coeff = 1;
    }
    /* }}} */

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    param_list_clear(pl);

    int rc = 0;
    struct afile_list a[1];
    memset(a,0,sizeof(a));
    if (read_afiles(a) == 0) {
        rc = 0;
        goto paradise;
    }



    /* First merge all files we find with similar [n0..n1[ range. Since
     * these files are expected to come from different sites, it's normal
     * to have different [n0..n1[ ranges progress at different speeds.
     */
    int did_merge = 0;
    for(unsigned int k0=0, k1 ; k0 < a->n ; k0 = k1) {
        unsigned int j0 = a->a[k0]->j0;
        unsigned int j1 = a->a[k0]->j1;
        k1 = k0;
        unsigned int n0 = a->a[k0]->n0;
        unsigned int n1 = a->a[k0]->n1;
        unsigned int j = j0;
        for( ; k1 < a->n && a->a[k1]->n0 == n0 ; k1++) {
            if (a->a[k1]->n1 != n1 || a->a[k1]->j0 != j) {
                fprintf(stderr, "Found inconsistent files A%u-%u.%u-%u and A%u-%u.%u-%u\n",
                        a->a[k0]->n0, a->a[k0]->n1, a->a[k0]->j0, a->a[k0]->j1,
                        a->a[k1]->n0, a->a[k1]->n1, a->a[k1]->j0, a->a[k1]->j1);
                exit(EXIT_FAILURE);
            }
            j = a->a[k1]->j1;
        }
        j1 = j;
        
        /* creating file [n0..n1[.., [j0..j1[ ; since we're just catting
         * the files, no trouble.
         */
        if (k1-k0 == 1)
            continue;
        did_merge++;
        FILE * f = fopen("A.temp", "wb");
        for(unsigned int k = k0 ; k < k1 ; k++) {
            char * tmp;
            int rc = asprintf(&tmp, "A%u-%u.%u-%u",
                    a->a[k]->n0,a->a[k]->n1,a->a[k]->j0,a->a[k]->j1);
            ASSERT_ALWAYS(rc >= 0);
            FILE * g = fopen(tmp, "rb");
            char buf[BUFSIZ];
            for( ; ; ) {
                int nr = fread(buf, 1, BUFSIZ, g);
                if (nr < BUFSIZ && ferror(g)) {
                    fprintf(stderr, "%s: %s\n", tmp, strerror(errno));
                    exit(EXIT_FAILURE);
                }
                if (nr == 0)
                    break;
                int nw = fwrite(buf, 1, nr, f);
                if (nr < nw) {
                    fprintf(stderr, "copying %s: %s\n", tmp, strerror(errno));
                    exit(EXIT_FAILURE);
                }
                if (nr < BUFSIZ)
                    break;
            }
            fclose(g);

            if (remove_old) {
                if (unlink(tmp) < 0) {
                    fprintf(stderr, "unlink(%s): %s\n", tmp, strerror(errno));
                    exit(EXIT_FAILURE);
                }
            }
            free(tmp);
        }
        fclose(f);
        {
            char * tmp;
            int r = asprintf(&tmp, "A%u-%u.%u-%u", n0,n1,j0,j1);
            ASSERT_ALWAYS(r >= 0);
            r = rename("A.temp", tmp);
            if (r < 0) {
                fprintf(stderr, "rename(A.temp, %s): %s\n",
                        tmp, strerror(errno));
                exit(EXIT_FAILURE);
            }
            free(tmp);
        }
    }

    if (did_merge && !remove_old) {
        fprintf(stderr, "Done some merges, but cannot continue unless --remove-old is specified\n");
        exit(EXIT_FAILURE);
    }

    /* Good. Now merge the other way around. Not clear it's really
     * something we want to do like this, though.
     */

    afile final;
    /* start with unset values */
    final->n0 = UINT_MAX;
    final->n1 = UINT_MAX;
    final->j0 = UINT_MAX;
    final->j1 = UINT_MAX;

    if (read_afiles(a) == 0) {
        rc = 0;
        goto paradise;
    }

    FILE * f = fopen("A.temp", "wb");
    for(unsigned int k0 = 0, k1 ; k0 < a->n ; k0 = k1) {
        unsigned int j0 = a->a[k0]->j0;
        unsigned int j1 = a->a[k0]->j1;
        k1 = k0;
        unsigned int n0 = a->a[k0]->n0;
        unsigned int n1 = a->a[k0]->n1;
        unsigned int n = n0;
        for( ; k1 < a->n && a->a[k1]->j0 == j0 ; k1++) {
            if (a->a[k1]->j1 != j1 || a->a[k1]->n0 != n) {
                fprintf(stderr, "Found inconsistent files A%u-%u.%u-%u and A%u-%u.%u-%u\n",
                        a->a[k0]->n0, a->a[k0]->n1, a->a[k0]->j0, a->a[k0]->j1,
                        a->a[k1]->n0, a->a[k1]->n1, a->a[k1]->j0, a->a[k1]->j1);
                exit(EXIT_FAILURE);
            }
            n = a->a[k1]->n1;
        }
        n1 = n;
        
        /* Files [k0] to [k1] can be collected */

        if (    (final->n0 != UINT_MAX && n0 != final->n0) ||
                (final->n1 != UINT_MAX && n1 != final->n1) ||
                (final->j1 != UINT_MAX && j0 != final->j1))
        {
            fprintf(stderr, "Cannot append A%u-%u.%u-%u to A%u-%u.%u-%u\n",
                    n0,n1,j0,j1,
                    final->n0,final->n1,final->j0,final->j1);
            rc = 1;
            break;
        }

        FILE ** rs = malloc((k1-k0) * sizeof(FILE *));
        for(unsigned int k = k0 ; k < k1 ; k++) {
            char * tmp;
            int rc = asprintf(&tmp, "A%u-%u.%u-%u",
                    a->a[k]->n0,a->a[k]->n1,a->a[k]->j0,a->a[k]->j1);
            ASSERT_ALWAYS(rc >= 0);
            rs[k - k0] = fopen(tmp, "rb");
            if (rs[k-k0] == NULL) {
                fprintf(stderr, "fopen(%s): %s\n", tmp, strerror(errno));
                exit(EXIT_FAILURE);
            }
            free(tmp);
        }

        char * buf = malloc((n1-n0)*bits_per_coeff/CHAR_BIT);

        final->j1 = j0;

        for(unsigned int j = j0 ; j < j1 ; j++) {
            for(int i = 0 ; i < bw->m ; i++) {
                char * ptr = buf;
                size_t rz;
                size_t sz;
                for(unsigned int k = k0 ; k < k1 ; k++) {
                    sz = (a->a[k]->n1 - a->a[k]->n0) * bits_per_coeff/ CHAR_BIT;
                    rz = fread(ptr, 1, sz, rs[k-k0]);
                    if (rz < sz) {
                        rc = 2;
                        fprintf(stderr, "fwrite: short read\n");
                        exit(EXIT_FAILURE);
                    }
                    ptr += sz;
                }
                sz = (n1-n0)* bits_per_coeff/CHAR_BIT;
                rz = fwrite(buf, 1, sz, f);
                if (rz != sz) {
                    fprintf(stderr, "fwrite: short write\n");
                    exit(EXIT_FAILURE);
                }
            }
            final->j1++;
        }
        free(buf);

        if (final->j0 == UINT_MAX) { final->j0 = j0; }
        if (final->n0 == UINT_MAX) { final->n0 = n0; }
        if (final->n1 == UINT_MAX) { final->n1 = n1; }

        if (fflush(f) != 0) {
            // we're in trouble
            fprintf(stderr, "fflush(): %s\n", strerror(errno));
            exit(EXIT_FAILURE);
        }

        final->j1 = j1;

        for(unsigned int k = k0 ; k < k1 ; k++) {
            fclose(rs[k-k0]);
            if (!remove_old) continue;
            char * tmp;
            int rc = asprintf(&tmp, "A%u-%u.%u-%u",
                    a->a[k]->n0,a->a[k]->n1,a->a[k]->j0,a->a[k]->j1);
            ASSERT_ALWAYS(rc >= 0);
            if (unlink(tmp) < 0) {
                fprintf(stderr, "unlink(%s): %s\n", tmp, strerror(errno));
                exit(EXIT_FAILURE);
            }
            free(tmp);
        }
        free(rs);
    }
    fclose(f);
    if (final->j0 != UINT_MAX) {
        char * tmp;
        int r;
        r = asprintf(&tmp, "A%u-%u.%u-%u", final->n0,final->n1,final->j0,final->j1);
        r = rename("A.temp", tmp);
        if (r < 0) {
            fprintf(stderr, "rename(A.temp, %s): %s\n", tmp, strerror(errno));
            exit(EXIT_FAILURE);
        }
        printf("%s\n",tmp);
        free(tmp);
    }
    free(a->a);

paradise:
    bw_common_clear(bw);

    return rc;
}
