#include "cado.h" // IWYU pragma: keep
#include <string.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdint.h>
#include <sys/stat.h>
#include <errno.h>
#include <stdlib.h>

#include "balancing.hpp"
#include "portability.h" // asprintf // IWYU pragma: keep
#include "fix-endianness.h" // fread32_little
#include "crc.h"        // cado_crc_lfsr
#include "misc.h"       // has_suffix
#include "macros.h"

void balancing_set_row_col_count(balancing & bal)
{
    unsigned int s = bal.nh * bal.nv;
    unsigned int b = iceildiv(bal.nrows, s);
    for( ; b % MINIMUM_ITEMS_IN_BWC_CHUNKS ; b++);
    bal.trows = s * b;
    b = iceildiv(bal.ncols, s);
    for( ; b % MINIMUM_ITEMS_IN_BWC_CHUNKS ; b++);
    bal.tcols = s * b;
    if (bal.flags & FLAG_REPLICATE) {
        bal.tcols = bal.trows = MAX(bal.trows, bal.tcols);
    }
}

void balancing_finalize(balancing & bal)
{
    cado_crc_lfsr l;
    cado_crc_lfsr_init(l);
    uint32_t w = 0;
    balancing_set_row_col_count(bal);

    /* It does not make sense to say that we want to replicate
     * permutations if we have both a row and a column permutation in the
     * file, right?
     */
    uint32_t c = bal.flags & FLAG_ROWPERM;
    uint32_t r = bal.flags & FLAG_COLPERM;
    uint32_t a = bal.flags & FLAG_REPLICATE;
    ASSERT_ALWAYS(!(c && r && a));

    if (bal.flags & FLAG_ROWPERM) {
        w = cado_crc_lfsr_turn32_little(l, bal.rowperm, bal.trows * sizeof(uint32_t));
    }
    if (bal.flags & FLAG_COLPERM) {
        w = cado_crc_lfsr_turn32_little(l, bal.colperm, bal.tcols * sizeof(uint32_t));
    }
    cado_crc_lfsr_clear(l);
    bal.checksum = w;
    if (bal.flags & FLAG_REPLICATE) {
        // a trick to identify conjugated perms.
        bal.checksum &= ~0xff;
    }
}

void balancing_write_inner(balancing const & bal, const char * filename)
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
    ASSERT_ALWAYS(sizeof(balancing_header) == 64);
    rc += fwrite32_little(&bal.zero, 1, pfile);
    rc += fwrite32_little(&bal.magic, 1, pfile);
    rc += fwrite32_little(&bal.nh, 1, pfile);
    rc += fwrite32_little(&bal.nv, 1, pfile);
    rc += fwrite32_little(&bal.nrows, 1, pfile);
    rc += fwrite32_little(&bal.ncols, 1, pfile);
    rc += fwrite32_little(&bal.nzrows, 1, pfile);
    rc += fwrite32_little(&bal.nzcols, 1, pfile);
    rc += fwrite64_little(&bal.ncoeffs, 1, pfile);
    rc += fwrite32_little(&bal.checksum, 1, pfile);
    rc += fwrite32_little(&bal.flags, 1, pfile);
    rc += fwrite32_little(bal.pshuf, 2, pfile);
    rc += fwrite32_little(bal.pshuf_inv, 2, pfile);
    ASSERT_ALWAYS(rc == 15);
    if (bal.flags & FLAG_ROWPERM) {
        rc = fwrite32_little(bal.rowperm, bal.trows, pfile);
        ASSERT_ALWAYS(rc == (int) bal.trows);
    }
    if (bal.flags & FLAG_COLPERM) {
        rc = fwrite32_little(bal.colperm, bal.tcols, pfile);
        ASSERT_ALWAYS(rc == (int) bal.tcols);
    }
    fclose(pfile);
}

void balancing_write(balancing const & bal, const char * mfile, const char * suggest)
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
    char const * ext[2] = { ".bin", ".txt" };
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
                suggest, q, bal.nh, bal.nv, bal.checksum);
    } else {
        rc = asprintf(&filename, "%s.%dx%d.%08" PRIx32 ".bin",
                dup_prefix, bal.nh, bal.nv, bal.checksum);
    }
    ASSERT_ALWAYS(rc >= 0);
    free(dup_prefix);
    balancing_write_inner(bal, filename);
    free(filename);
}

void balancing_read_header_inner(balancing & bal, FILE * pfile)
{
    int rc = 0;
    ASSERT_ALWAYS(pfile);
    ASSERT_ALWAYS(sizeof(balancing_header) == 64);
    rc += fread32_little(&bal.zero, 1, pfile);
    rc += fread32_little(&bal.magic, 1, pfile);
    rc += fread32_little(&bal.nh, 1, pfile);
    rc += fread32_little(&bal.nv, 1, pfile);
    rc += fread32_little(&bal.nrows, 1, pfile);
    rc += fread32_little(&bal.ncols, 1, pfile);
    rc += fread32_little(&bal.nzrows, 1, pfile);
    rc += fread32_little(&bal.nzcols, 1, pfile);
    rc += fread64_little(&bal.ncoeffs, 1, pfile);
    rc += fread32_little(&bal.checksum, 1, pfile);
    rc += fread32_little(&bal.flags, 1, pfile);
    rc += fread32_little(bal.pshuf, 2, pfile);
    rc += fread32_little(bal.pshuf_inv, 2, pfile);
    ASSERT_ALWAYS(rc == 15);
    if (bal.zero != 0 || bal.magic != BALANCING_MAGIC) {
        fprintf(stderr, "Incompatible balancing file\n");
        exit(EXIT_FAILURE);
    }
    /* It does not make sense to say that we want to replicate
     * permutations if we have both a row and a column permutation in the
     * file, right?
     */
    uint32_t c = bal.flags & FLAG_ROWPERM;
    uint32_t r = bal.flags & FLAG_COLPERM;
    uint32_t a = bal.flags & FLAG_REPLICATE;
    ASSERT_ALWAYS(!(c && r && a));

}

void balancing_read_header(balancing & bal, std::string const & filename)
{
    FILE * pfile;
    ASSERT_ALWAYS(!filename.empty());
    char * derived = derived_filename(filename.c_str(), "hdr", ".bin");
    pfile = fopen(filename.c_str(), "rb");
    if (pfile == NULL) {
        pfile = fopen(derived, "rb");
        if (pfile == NULL) {
            fprintf(stderr, "Cannot read %s nor %s: %s\n",
                    filename.c_str(), derived, strerror(errno));
            abort();
        }
    }
    balancing_read_header_inner(bal, pfile);
    balancing_set_row_col_count(bal);
    fclose(pfile);
    free(derived);
}

void balancing_read(balancing & bal, std::string const & filename)
{
    FILE * pfile;

    ASSERT_ALWAYS(!filename.empty());
    pfile = fopen (filename.c_str(), "rb");
    if (pfile == NULL) {
        perror(filename.c_str());
        abort();
    }
    balancing_read_header_inner(bal, pfile);
    balancing_set_row_col_count(bal);
    if (bal.flags & FLAG_ROWPERM) {
        bal.rowperm = (uint32_t *) malloc(bal.trows * sizeof(uint32_t));
        int rc = fread32_little(bal.rowperm, bal.trows, pfile);
        ASSERT_ALWAYS(rc == (int) bal.trows);
    }
    if (bal.flags & FLAG_COLPERM) {
        bal.colperm = (uint32_t *) malloc(bal.tcols * sizeof(uint32_t));
        int rc = fread32_little(bal.colperm, bal.tcols, pfile);
        ASSERT_ALWAYS(rc == (int) bal.tcols);
    }
    fclose(pfile);
}

void balancing_init(balancing & bal)
{
    bal.magic = BALANCING_MAGIC;
}
void balancing_clear(balancing & bal)
{
    if (bal.colperm) free(bal.colperm);
    if (bal.rowperm) free(bal.rowperm);
}
