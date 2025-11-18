#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cinttypes>

#include <sys/stat.h>

#include "select_mpi.h"
#include "abfiles.hpp"
#include "macros.h"
#include "mpi_proxies.hpp"
#include "portability.h"        // strdup

using cado_mpi::mpi_data_agrees;

// interface for reading the list of (a,b)'s, with sort of a random
// access (for the starting point only).

#define ABFILE_MAX_LINE_LENGTH  256

void ab_source_init(ab_source_ptr ab, const char * fname, int rank, int root, MPI_Comm comm)
{
    memset(ab, 0, sizeof(ab_source));
    ab->fname0 = fname;
    const char * magic;
    if ((magic = strstr(fname, ".prep.")) != NULL) {
        // then assume kleinjung format.
        ab->prefix = (char*) malloc(magic - fname + 1);
        strncpy(ab->prefix, fname, magic-fname);
        ab->prefix[magic-fname]='\0';
        magic++;
        int fnum;
        if (sscanf(magic, "prep.%d.rel.%d", &ab->depnum, &fnum) == 2) {
            ab->nfiles = -1;  // to be determined later on.
        } else {
            FATAL_ERROR_CHECK(1, "error in parsing filename");
        }
    } else if ((magic = strstr(fname, ".dep.alg.")) != NULL) {
        // assume cado format (means only one file, so we don't need to
        // parse, really.
        ab->prefix = (char*) malloc(magic - fname + 1);
        strncpy(ab->prefix, fname, magic-fname);
        ab->prefix[magic-fname]='\0';
        magic++;
        if (sscanf(magic, "dep.alg.%d", &ab->depnum) == 1) {
            ab->nfiles = 0;
        } else {
            FATAL_ERROR_CHECK(1, "error in parsing filename");
        }
    } else if ((magic = strstr(fname, ".dep.side1.")) != NULL) {
        // assume cado format (means only one file, so we don't need to
        // parse, really.
        ab->prefix = (char*) malloc(magic - fname + 1);
        strncpy(ab->prefix, fname, magic-fname);
        ab->prefix[magic-fname]='\0';
        magic++;
        if (sscanf(magic, "dep.side1.%d", &ab->depnum) == 1) {
            ab->nfiles = 0;
        } else {
            FATAL_ERROR_CHECK(1, "error in parsing filename");
        }
    } else if ((magic = strstr(fname, ".dep.")) != NULL) {
        // assume cado format (means only one file, so we don't need to
        // parse, really.
        ab->prefix = (char*) malloc(magic - fname + 1);
        strncpy(ab->prefix, fname, magic-fname);
        ab->prefix[magic-fname]='\0';
        magic++;
        if (sscanf(magic, "dep.%d", &ab->depnum) == 1) {
            ab->nfiles = 0;
        } else {
            FATAL_ERROR_CHECK(1, "error in parsing filename");
        }
    } else {
        FATAL_ERROR_CHECK(1, "error in parsing filename");
    }

    // do some size estimations;
    size_t tsize = 0;
    struct stat sbuf[1];
    int rc;
    if (ab->nfiles == 0) {
        rc = stat(fname, sbuf);
        ASSERT_ALWAYS(mpi_data_agrees(rc, comm));
        ASSERT_ALWAYS(rc == 0);
        tsize=sbuf->st_size;
        // we have 2.5 non-digit bytes per file line. However we
        // don't know the line count, so we can't subtract. As a
        // guess, we read the first 8kb, and count the number of
        // lines in there.
        if (rank == root) {
            char buf[8192];
            FILE * f = fopen(fname, "r");
            rc = fread(buf, 1, sizeof(buf), f);
            ASSERT_ALWAYS(rc == sizeof(buf));
            fclose(f);
            int nrows_16k = 0;
            for(unsigned int i = 0 ; i < sizeof(buf) ; i++) {
                nrows_16k += buf[i] == '\n';
            }
            ab->nab_estim = (double) tsize * nrows_16k / sizeof(buf);
            // ab->digitbytes_estim = tsize - 2.5 * ab->nab_estim;
        }
        // XXX OK, this requires endianness consistency.
        MPI_Bcast(&ab->nab_estim, 1, CADO_MPI_SIZE_T, root, comm);
        ab->file_bases = (size_t *) malloc(2 * sizeof(size_t));
        ab->file_bases[0] = 0;
        ab->file_bases[1] = tsize;
        ab->totalsize = tsize;
    } else {
        size_t dummy;
        size_t hdrbytes;
        char line[ABFILE_MAX_LINE_LENGTH];
        if (rank == root) {
            FILE * f = fopen(fname, "r");
            char * xx = fgets(line, sizeof(line), f);
            DIE_ERRNO_DIAG(xx == NULL, "fgets(%s)", fname);
            rc = sscanf(line, "AB %zu %zu", &ab->nab, &dummy);
            DIE_ERRNO_DIAG(rc != 2, "parse(%s)", fname);
            hdrbytes = ftell(f);
            fclose(f);
        }
        MPI_Bcast(&hdrbytes, 1, CADO_MPI_SIZE_T, root, comm);

        ab->sname_len = strlen(fname) + 10;
        ab->sname = (char *) malloc(ab->sname_len);
        for(ab->nfiles = 0 ; ; ab->nfiles++) {
            snprintf(ab->sname, ab->sname_len, "%s.prep.%d.rel.%d",
                    ab->prefix, ab->depnum, ab->nfiles);
            rc = stat(ab->sname, sbuf);
            ASSERT_ALWAYS(rc == 0 || errno == ENOENT);
            ASSERT_ALWAYS(mpi_data_agrees(rc, comm));
            if (rc < 0) break;
            tsize += sbuf->st_size - hdrbytes;
        }
        ASSERT_ALWAYS(ab->nfiles > 0);
        ab->nab_estim = ab->nab;
        // ab->digitbytes_estim = tsize - 5 * ab->nab_estim;
        ab->file_bases = (size_t *) malloc((ab->nfiles+1) * sizeof(size_t));
        ab->file_bases[0] = 0;
        for(int i = 0 ; i < ab->nfiles ; i++) {
            snprintf(ab->sname, ab->sname_len, "%s.prep.%d.rel.%d",
                    ab->prefix, ab->depnum, i);
            rc = stat(ab->sname, sbuf);
            ASSERT_ALWAYS(rc == 0);
            ASSERT_ALWAYS(mpi_data_agrees(rc, comm));
            ab->file_bases[i+1]=ab->file_bases[i] + sbuf->st_size;
        }
        ab->totalsize = ab->file_bases[ab->nfiles];
    }
}

void ab_source_rewind(ab_source_ptr ab)
{
    if (ab->f) fclose(ab->f);
    ab->f = NULL;
    ab->c = 0;
    ab->nab = 0;
    ab->cpos = 0;
    ab->tpos = 0;
}

void ab_source_init_set(ab_source_ptr ab, ab_source_ptr ab0)
{
    memcpy(ab, ab0, sizeof(struct ab_source_s));
    ab->prefix = strdup(ab0->prefix);
    ab->sname = (char *) malloc(ab->sname_len);
    ab->file_bases = (size_t *) malloc((ab->nfiles+1) * sizeof(size_t));
    memcpy(ab->file_bases, ab0->file_bases, (ab->nfiles+1) * sizeof(size_t));
    ab->f = NULL;
    ab_source_rewind(ab);
}

void ab_source_clear(ab_source_ptr ab)
{
    ab_source_rewind(ab);
    free(ab->prefix);
    free(ab->sname);
    free(ab->file_bases);
    memset(ab, 0, sizeof(ab_source));
}

static int ab_openfile_internal(ab_source_ptr ab)
{
    const char * s;
    if (ab->nfiles == 0) {
        ab->f = fopen(s=ab->fname0, "r");
        ab->tpos = 0;
    } else {
        snprintf(ab->sname, ab->sname_len, "%s.prep.%d.rel.%d",
                ab->prefix, ab->depnum, ab->c);
        ab->f = fopen(s=ab->sname, "r");
        if (ab->f == NULL && errno == ENOENT)
            return 0;
        ab->tpos = ab->file_bases[ab->c];
        ab->cpos = 0;
        char header[80];
        char * rp = fgets(header, sizeof(header), ab->f);
        ASSERT_ALWAYS(rp);
    }
    ab->cpos = ftell(ab->f);
    ab->tpos += ab->cpos;
    DIE_ERRNO_DIAG(ab->f == NULL, "fopen(%s)", s);
    return 1;
}

int ab_source_next(ab_source_ptr ab, int64_t * a, uint64_t * b)
{
    if (ab->f) {
        int rc MAYBE_UNUSED;
        char line[ABFILE_MAX_LINE_LENGTH];
        char * xx = fgets(line, sizeof(line), ab->f);
        size_t cpos = ftell(ab->f);
        if (xx) {
            if (ab->nfiles == 0) {
                rc = sscanf(line, "%" SCNd64 " %" SCNu64, a, b);
                ASSERT(rc == 2);
            } else {
                int dummy;
                rc = sscanf(line, "%d %" SCNd64 " %" SCNu64, &dummy, a, b);
                ASSERT(rc == 3);
            }
            ab->tpos += cpos - ab->cpos;
            ab->cpos = cpos;
            ab->nab++;
            return 1;
        }
        fclose(ab->f); ab->f = NULL;
        ab->tpos += cpos - ab->tpos;
        // don't update cpos, it is not defined in this situation.
        if (ab->nfiles == 0)
            return 0;
        ab->c++;
    }
    if (ab_openfile_internal(ab) == 0)
        return 0;
    return ab_source_next(ab, a, b);
}

void ab_source_move_afterpos(ab_source_ptr ab, size_t offset)
{
    /* move the file pointer to the earliest non-header line starting at
     * an offset which is greater than or equal to offset.
     *
     * IF the current file position is already >= offset, do nothing.
     *
     * IF the current file position is < offset, then the returned offset
     * is >= offset. This is achieved by seeking to some pre_offset <
     * offset, and
     * advancing to the next line starting at >= offset.
     *
     * This is used to read a collection of files in chunks whose size is
     * governed by the file size.
     */

    /* note that the test below always succeeds for offset == 0 */
    if (ab->tpos >= offset)
        return;
    // otherwise it's really, really a can of worms.
    ASSERT_ALWAYS(ab->f == NULL);

    // which file ?
    FATAL_ERROR_CHECK(offset >= ab->totalsize,
            "attempt to seek beyond end of files");

    size_t pre_offset = MAX(offset, 10) - 10;
    ASSERT_ALWAYS(pre_offset < offset);

    if (ab->nfiles == 0) {
        // well, does not really make a lot of sense here, but anyway.
        ab_openfile_internal(ab);
        fseek(ab->f, pre_offset, SEEK_SET);
    } else {
        for( ; ab->file_bases[ab->c+1] <= pre_offset ; ab->c++) ;
        ab_openfile_internal(ab);
        // so we know that
        // ab->file_bases[ab->c] <= pre_offset < ab->file_bases[ab->c+1]
        fseek(ab->f, pre_offset - ab->file_bases[ab->c], SEEK_SET);
        // note that it is possible that tpos, as obtained after
        // openfile, is >= offset. but we know that we'll be able to seek
        // to pre_offset, and that one is appropriately < offset, so no
        // special case.
    }
    char line[ABFILE_MAX_LINE_LENGTH];
    char * xx = fgets(line, sizeof(line), ab->f);
    DIE_ERRNO_DIAG(xx == NULL, "fgets(%s)", ab->nfiles ? ab->sname : ab->fname0);
    size_t cpos = ftell(ab->f);
    ab->tpos += cpos - ab->cpos;
    ab->cpos = cpos;
    for(int n_adjust = 0 ; ab->tpos < offset ; n_adjust++) {
        FATAL_ERROR_CHECK(n_adjust > 10, "adjustment on the runaway");
        int64_t a;
        uint64_t b;
        int r = ab_source_next(ab, &a, &b);
        FATAL_ERROR_CHECK(r == 0, "adjustment failed");
    }
}
