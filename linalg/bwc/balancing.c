#include "cado.h" // IWYU pragma: keep
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdint.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdlib.h>

#include "balancing.h"
#include "cheating_vec_init.h"
#include "portability.h" // asprintf // IWYU pragma: keep
#include "fix-endianness.h" // fread32_little
#include "crc.h"        // cado_crc_lfsr
#include "misc.h"       // has_suffix
#include "macros.h"

void balancing_set_row_col_count(balancing_ptr bal)
{
    unsigned int s = bal->h->nh * bal->h->nv;
    unsigned int b = iceildiv(bal->h->nrows, s);
    for( ; b % (FORCED_ALIGNMENT_ON_MPFQ_VEC_TYPES / MINIMUM_ITEM_SIZE_OF_MPFQ_VEC_TYPES) ; b++);
    bal->trows = s * b;
    b = iceildiv(bal->h->ncols, s);
    for( ; b % (FORCED_ALIGNMENT_ON_MPFQ_VEC_TYPES / MINIMUM_ITEM_SIZE_OF_MPFQ_VEC_TYPES) ; b++);
    bal->tcols = s * b;
    if (bal->h->flags & FLAG_REPLICATE) {
        bal->tcols = bal->trows = MAX(bal->trows, bal->tcols);
    }
}

void balancing_finalize(balancing_ptr bal)
{
    cado_crc_lfsr l;
    cado_crc_lfsr_init(l);
    uint32_t w = 0;
    balancing_set_row_col_count(bal);

    /* It does not make sense to say that we want to replicate
     * permutations if we have both a row and a column permutation in the
     * file, right?
     */
    uint32_t c = bal->h->flags & FLAG_ROWPERM;
    uint32_t r = bal->h->flags & FLAG_COLPERM;
    uint32_t a = bal->h->flags & FLAG_REPLICATE;
    ASSERT_ALWAYS(!(c && r && a));

    if (bal->h->flags & FLAG_ROWPERM) {
        w = cado_crc_lfsr_turn32_little(l, bal->rowperm, bal->trows * sizeof(uint32_t));
    }
    if (bal->h->flags & FLAG_COLPERM) {
        w = cado_crc_lfsr_turn32_little(l, bal->colperm, bal->tcols * sizeof(uint32_t));
    }
    cado_crc_lfsr_clear(l);
    bal->h->checksum = w;
    if (bal->h->flags & FLAG_REPLICATE) {
        // a trick to identify conjugated perms.
        bal->h->checksum &= ~0xff;
    }
}

void balancing_write_inner(balancing_ptr bal, const char * filename)
{
    FILE * pfile;
    fprintf(stderr, "Writing balancing data to %s\n", filename);
    pfile = fopen(filename, "wb");
    if (pfile == NULL) {
        perror(filename);
        abort();
    }
    int rc = 0;
    /* Any change to the balancing_header structure must propagate here */
    ASSERT_ALWAYS(sizeof(struct balancing_header_s) == 64);
    rc += fwrite32_little(&bal->h->zero, 1, pfile);
    rc += fwrite32_little(&bal->h->magic, 1, pfile);
    rc += fwrite32_little(&bal->h->nh, 1, pfile);
    rc += fwrite32_little(&bal->h->nv, 1, pfile);
    rc += fwrite32_little(&bal->h->nrows, 1, pfile);
    rc += fwrite32_little(&bal->h->ncols, 1, pfile);
    rc += fwrite32_little(&bal->h->nzrows, 1, pfile);
    rc += fwrite32_little(&bal->h->nzcols, 1, pfile);
    rc += fwrite64_little(&bal->h->ncoeffs, 1, pfile);
    rc += fwrite32_little(&bal->h->checksum, 1, pfile);
    rc += fwrite32_little(&bal->h->flags, 1, pfile);
    rc += fwrite32_little(bal->h->pshuf, 2, pfile);
    rc += fwrite32_little(bal->h->pshuf_inv, 2, pfile);
    ASSERT_ALWAYS(rc == 15);
    if (bal->h->flags & FLAG_ROWPERM) {
        rc = fwrite32_little(bal->rowperm, bal->trows, pfile);
    }
    if (bal->h->flags & FLAG_COLPERM) {
        rc = fwrite32_little(bal->colperm, bal->tcols, pfile);
    }
    fclose(pfile);
}

void balancing_write(balancing_ptr bal, const char * mfile, const char * suggest)
{
    /* the semantics of -out for this program are farily weird. If it's
     * a file, then we'll use that as an output name (this is the call to
     * balancing_write_inner() early on below). If it's a
     * directory, we'll place the balancing file named the standard way
     * there, as done by the asprintf naming later on.
     */
    int suggestion_is_directory = 0;
    if (suggest && strlen(suggest)) {
        struct stat sb[1];
        suggestion_is_directory = (stat(suggest, sb) == 0 && S_ISDIR(sb->st_mode));
    }

    if (suggest && strlen(suggest) && !suggestion_is_directory) {
        balancing_write_inner(bal, suggest);
        /* If we received "-out", don't store the balancing file with the
         * default name -- that would be rather odd.
         */
        return;
    }

    char * dup_prefix=strdup(mfile);
    char * ext[2] = { ".bin", ".txt" };
    for(int j = 0 ; j < 2 ; j++) {
        if (has_suffix(dup_prefix, ext[j])) {
            dup_prefix[strlen(dup_prefix)-strlen(ext[j])]='\0';
            break;
        }
    }

    char * filename;
    int rc;
    if (suggestion_is_directory) {
        char * q = strrchr(dup_prefix, '/');
        if (q) { q++; } else { q = dup_prefix; }
        rc = asprintf(&filename, "%s/%s.%dx%d.%08" PRIx32 ".bin",
                suggest, q, bal->h->nh, bal->h->nv, bal->h->checksum);
    } else {
        rc = asprintf(&filename, "%s.%dx%d.%08" PRIx32 ".bin",
                dup_prefix, bal->h->nh, bal->h->nv, bal->h->checksum);
    }
    ASSERT_ALWAYS(rc >= 0);
    free(dup_prefix);
    balancing_write_inner(bal, filename);
    free(filename);
}

void balancing_read_header_inner(balancing_ptr bal, FILE * pfile)
{
    int rc = 0;
    ASSERT_ALWAYS(pfile);
    ASSERT_ALWAYS(sizeof(struct balancing_header_s) == 64);
    rc += fread32_little(&bal->h->zero, 1, pfile);
    rc += fread32_little(&bal->h->magic, 1, pfile);
    rc += fread32_little(&bal->h->nh, 1, pfile);
    rc += fread32_little(&bal->h->nv, 1, pfile);
    rc += fread32_little(&bal->h->nrows, 1, pfile);
    rc += fread32_little(&bal->h->ncols, 1, pfile);
    rc += fread32_little(&bal->h->nzrows, 1, pfile);
    rc += fread32_little(&bal->h->nzcols, 1, pfile);
    rc += fread64_little(&bal->h->ncoeffs, 1, pfile);
    rc += fread32_little(&bal->h->checksum, 1, pfile);
    rc += fread32_little(&bal->h->flags, 1, pfile);
    rc += fread32_little(bal->h->pshuf, 2, pfile);
    rc += fread32_little(bal->h->pshuf_inv, 2, pfile);
    ASSERT_ALWAYS(rc == 15);
    if (bal->h->zero != 0 || bal->h->magic != BALANCING_MAGIC) {
        fprintf(stderr, "Incompatible balancing file\n");
        exit(EXIT_FAILURE);
    }
    /* It does not make sense to say that we want to replicate
     * permutations if we have both a row and a column permutation in the
     * file, right?
     */
    uint32_t c = bal->h->flags & FLAG_ROWPERM;
    uint32_t r = bal->h->flags & FLAG_COLPERM;
    uint32_t a = bal->h->flags & FLAG_REPLICATE;
    ASSERT_ALWAYS(!(c && r && a));

}

void balancing_read_header(balancing_ptr bal, const char * filename)
{
    FILE * pfile;
    char * derived = derived_filename(filename, "hdr", ".bin");
    ASSERT_ALWAYS(filename);
    pfile = fopen(filename, "rb");
    if (pfile == NULL) {
        pfile = fopen(derived, "rb");
        if (pfile == NULL) {
            fprintf(stderr, "Cannot read %s nor %s: %s\n",
                    filename, derived, strerror(errno));
            abort();
        }
    }
    balancing_read_header_inner(bal, pfile);
    balancing_set_row_col_count(bal);
    fclose(pfile);
    free(derived);
}

void balancing_read(balancing_ptr bal, const char * filename)
{
    FILE * pfile;

    ASSERT_ALWAYS(filename);
    pfile = fopen (filename, "rb");
    if (pfile == NULL) {
        perror(filename);
        abort();
    }
    balancing_read_header_inner(bal, pfile);
    balancing_set_row_col_count(bal);
    if (bal->h->flags & FLAG_ROWPERM) {
        bal->rowperm = malloc(bal->trows * sizeof(uint32_t));
        int rc = fread32_little(bal->rowperm, bal->trows, pfile);
        ASSERT_ALWAYS(rc == (int) bal->trows);
    }
    if (bal->h->flags & FLAG_COLPERM) {
        bal->colperm = malloc(bal->tcols * sizeof(uint32_t));
        int rc = fread32_little(bal->colperm, bal->tcols, pfile);
        ASSERT_ALWAYS(rc == (int) bal->tcols);
    }
    fclose(pfile);
}

void balancing_init(balancing_ptr bal)
{
    memset(bal, 0, sizeof(balancing));
    bal->h->magic = BALANCING_MAGIC;
}
void balancing_clear(balancing_ptr bal)
{
    if (bal->colperm) free(bal->colperm);
    if (bal->rowperm) free(bal->rowperm);
    memset(bal, 0, sizeof(balancing));
}
