#include "cado.h" // IWYU pragma: keep

#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <memory>
#include <string>

#include <sys/stat.h>

#include "fmt/base.h"
#include "fmt/format.h"

#include "balancing.hpp"
#include "portability.h"
#include "fix-endianness.h"
#include "crc.h"
#include "misc.h"
#include "macros.h"
#include "utils_cxx.hpp"

void balancing_set_row_col_count(balancing & bal)
{
    unsigned int const s = bal.nh * bal.nv;
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
    uint32_t const c = bal.flags & FLAG_ROWPERM;
    uint32_t const r = bal.flags & FLAG_COLPERM;
    uint32_t const a = bal.flags & FLAG_REPLICATE;
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

void balancing_write_inner(balancing const & bal, std::string const & filename)
{
    fmt::print(stderr, "Writing balancing data to {}\n", filename);
    std::unique_ptr<FILE, delete_FILE> const pfile(fopen(filename.c_str(), "wb"));
    if (!pfile) {
        perror(filename.c_str());
        abort();
    }
    size_t rc = 0;
    /* Any change to the balancing_header structure must propagate here */
    ASSERT_ALWAYS(sizeof(balancing_header) == 64);
    rc += fwrite32_little(&bal.zero, 1, pfile.get());
    rc += fwrite32_little(&bal.magic, 1, pfile.get());
    rc += fwrite32_little(&bal.nh, 1, pfile.get());
    rc += fwrite32_little(&bal.nv, 1, pfile.get());
    rc += fwrite32_little(&bal.nrows, 1, pfile.get());
    rc += fwrite32_little(&bal.ncols, 1, pfile.get());
    rc += fwrite32_little(&bal.nzrows, 1, pfile.get());
    rc += fwrite32_little(&bal.nzcols, 1, pfile.get());
    rc += fwrite64_little(&bal.ncoeffs, 1, pfile.get());
    rc += fwrite32_little(&bal.checksum, 1, pfile.get());
    rc += fwrite32_little(&bal.flags, 1, pfile.get());
    rc += fwrite32_little(bal.pshuf, 2, pfile.get());
    rc += fwrite32_little(bal.pshuf_inv, 2, pfile.get());
    ASSERT_ALWAYS(rc == 15);
    if (bal.flags & FLAG_ROWPERM) {
        rc = fwrite32_little(bal.rowperm, bal.trows, pfile.get());
        ASSERT_ALWAYS(rc == bal.trows);
    }
    if (bal.flags & FLAG_COLPERM) {
        rc = fwrite32_little(bal.colperm, bal.tcols, pfile.get());
        ASSERT_ALWAYS(rc == bal.tcols);
    }
}

void balancing_write(balancing const & bal, std::string const & mfile, std::string const & suggest)
{
    /* the semantics of -out for this program are farily weird. If it's
     * a file, then we'll use that as an output name (this is the call to
     * balancing_write_inner() early on below). If it's a
     * directory, we'll place the balancing file named the standard way
     * there, as done by the asprintf naming later on.
     */
    bool suggestion_is_directory = false;
    if (!suggest.empty()) {
        struct stat sb[1];
        suggestion_is_directory = (stat(suggest.c_str(), sb) == 0 && S_ISDIR(sb->st_mode));
    }

    if (!suggest.empty() && !suggestion_is_directory) {
        balancing_write_inner(bal, suggest);
        /* If we received "-out", don't store the balancing file with the
         * default name -- that would be rather odd.
         */
        return;
    }

    /* TODO: refactor at least part of this with code in build_matcache? */
    std::string locfile;
    {
        if (suggestion_is_directory) {
            auto it = mfile.rfind('/');
            it = (it == std::string::npos) ? 0 : (it + 1);
            locfile = fmt::format("{}/{}", suggest, mfile.substr(it));
        } else {
            locfile = mfile;
        }
        auto it = locfile.rfind(".bin");
        if (it != std::string::npos)
            locfile.erase(it, locfile.size());
    }

    locfile = fmt::format("{}.{}{}.{}.bin",
                locfile, bal.nh, bal.nv, bal.checksum);

    balancing_write_inner(bal, locfile);
}

static void balancing_read_header_inner(balancing & bal, FILE * pfile)
{
    size_t rc = 0;
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
    uint32_t const c = bal.flags & FLAG_ROWPERM;
    uint32_t const r = bal.flags & FLAG_COLPERM;
    uint32_t const a = bal.flags & FLAG_REPLICATE;
    ASSERT_ALWAYS(!(c && r && a));

    ASSERT_ALWAYS(bal.nh);
    ASSERT_ALWAYS(bal.nv);
}

void balancing_read_header(balancing & bal, std::string const & filename)
{
    ASSERT_ALWAYS(!filename.empty());
    std::string hdrname = std::unique_ptr<char, free_delete<char>>(
            derived_filename(filename.c_str(), "hdr", ".bin")).get();
    std::unique_ptr<FILE, delete_FILE> pfile(fopen(filename.c_str(), "rb"));
    if (!pfile) {
        pfile.reset(fopen(hdrname.c_str(), "rb"));
        if (!pfile) {
            fmt::print(stderr, "Cannot read {} nor {}: {}\n",
                    filename, hdrname, strerror(errno));
            abort();
        }
    }
    balancing_read_header_inner(bal, pfile.get());
    balancing_set_row_col_count(bal);
}

void balancing_read(balancing & bal, std::string const & filename)
{
    ASSERT_ALWAYS(!filename.empty());
    const std::unique_ptr<FILE, delete_FILE> pfile (fopen (filename.c_str(), "rb"));
    if (!pfile) {
        perror(filename.c_str());
        abort();
    }
    balancing_read_header_inner(bal, pfile.get());
    balancing_set_row_col_count(bal);
    if (bal.flags & FLAG_ROWPERM) {
        bal.rowperm = (uint32_t *) malloc(bal.trows * sizeof(uint32_t));
        size_t const rc = fread32_little(bal.rowperm, bal.trows, pfile.get());
        ASSERT_ALWAYS(rc == bal.trows);
    }
    if (bal.flags & FLAG_COLPERM) {
        bal.colperm = (uint32_t *) malloc(bal.tcols * sizeof(uint32_t));
        size_t const rc = fread32_little(bal.colperm, bal.tcols, pfile.get());
        ASSERT_ALWAYS(rc == bal.tcols);
    }
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
