#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cinttypes>
#include <cerrno>

#include <sys/stat.h>

#include "select_mpi.h"
#include "ab_source.hpp"
#include "macros.h"
#include "mpi_proxies.hpp"
#include "portability.h"        // strdup

using cado_mpi::mpi_data_agrees;

// interface for reading the list of (a,b)'s, with sort of a random
// access (for the starting point only).

#define ABFILE_MAX_LINE_LENGTH  256

ab_source::ab_source(std::string const & fname, int rank, int root, MPI_Comm comm)
{
    init(fname, rank, root, comm);
}

ab_source::ab_source(ab_source const & o)
    : fname0(o.fname0)
    , sname(o.sname)
    , sname_len(o.sname_len)
    , nab(o.nab)
    , prefix(o.prefix)
    , depnum(o.depnum)
    , nab_estim(o.nab_estim)
    , file_bases(o.file_bases)
    , nfiles(o.nfiles)
    , totalsize(o.totalsize)
{
}

ab_source & ab_source::operator=(ab_source const & o)
{
    if (this != &o) {
        rewind();
        fname0 = o.fname0;
        sname = o.sname;
        sname_len = o.sname_len;
        nab = o.nab;
        prefix = o.prefix;
        depnum = o.depnum;
        nab_estim = o.nab_estim;
        file_bases = o.file_bases;
        nfiles = o.nfiles;
        totalsize = o.totalsize;
        f = nullptr;
        c = 0;
        cpos = 0;
        tpos = 0;
    }
    return *this;
}

ab_source::ab_source(ab_source && o) noexcept
    : fname0(std::move(o.fname0))
    , sname(std::move(o.sname))
    , sname_len(o.sname_len)
    , nab(o.nab)
    , prefix(std::move(o.prefix))
    , depnum(o.depnum)
    , nab_estim(o.nab_estim)
    , file_bases(std::move(o.file_bases))
    , nfiles(o.nfiles)
    , totalsize(o.totalsize)
    , f(o.f)
    , c(o.c)
    , cpos(o.cpos)
    , tpos(o.tpos)
{
    o.f = nullptr;
    o.c = 0;
    o.cpos = 0;
    o.tpos = 0;
}

ab_source & ab_source::operator=(ab_source && o) noexcept
{
    if (this != &o) {
        rewind();
        fname0 = std::move(o.fname0);
        sname = std::move(o.sname);
        sname_len = o.sname_len;
        nab = o.nab;
        prefix = std::move(o.prefix);
        depnum = o.depnum;
        nab_estim = o.nab_estim;
        file_bases = std::move(o.file_bases);
        nfiles = o.nfiles;
        totalsize = o.totalsize;
        f = o.f;
        c = o.c;
        cpos = o.cpos;
        tpos = o.tpos;
        o.f = nullptr;
        o.c = 0;
        o.cpos = 0;
        o.tpos = 0;
    }
    return *this;
}

ab_source::~ab_source()
{
    rewind();
}

void ab_source::rewind()
{
    if (f) {
        fclose(f);
        f = nullptr;
    }
    c = 0;
    nab = 0;
    cpos = 0;
    tpos = 0;
}

void ab_source::init(std::string const & fname, int rank, int root, MPI_Comm comm)
{
    rewind();
    fname0 = fname;

    size_t pos;
    if ((pos = fname.find(".prep.")) != std::string::npos) {
        // then assume kleinjung format.
        prefix = fname.substr(0, pos);
        int fnum;
        if (sscanf(fname.c_str() + pos + 1, "prep.%d.rel.%d", &depnum, &fnum) == 2) {
            nfiles = -1;  // to be determined later on.
        } else {
            FATAL_ERROR_CHECK(1, "error in parsing filename");
        }
    } else if ((pos = fname.find(".dep.alg.")) != std::string::npos) {
        // assume cado format (means only one file, so we don't need to parse, really.
        prefix = fname.substr(0, pos);
        if (sscanf(fname.c_str() + pos + 1, "dep.alg.%d", &depnum) == 1) {
            nfiles = 0;
        } else {
            FATAL_ERROR_CHECK(1, "error in parsing filename");
        }
    } else if ((pos = fname.find(".dep.side1.")) != std::string::npos) {
        // assume cado format
        prefix = fname.substr(0, pos);
        if (sscanf(fname.c_str() + pos + 1, "dep.side1.%d", &depnum) == 1) {
            nfiles = 0;
        } else {
            FATAL_ERROR_CHECK(1, "error in parsing filename");
        }
    } else if ((pos = fname.find(".dep.")) != std::string::npos) {
        // assume cado format
        prefix = fname.substr(0, pos);
        if (sscanf(fname.c_str() + pos + 1, "dep.%d", &depnum) == 1) {
            nfiles = 0;
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
    if (nfiles == 0) {
        rc = stat(fname.c_str(), sbuf);
        ASSERT_ALWAYS(mpi_data_agrees(rc, comm));
        ASSERT_ALWAYS(rc == 0);
        tsize = sbuf->st_size;
        if (rank == root) {
            char buf[8192];
            FILE * file = fopen(fname.c_str(), "r");
            rc = fread(buf, 1, sizeof(buf), file);
            ASSERT_ALWAYS(rc == sizeof(buf));
            fclose(file);
            int nrows_16k = 0;
            for(unsigned int i = 0 ; i < sizeof(buf) ; i++) {
                nrows_16k += buf[i] == '\n';
            }
            nab_estim = (double) tsize * nrows_16k / sizeof(buf);
        }
        MPI_Bcast(&nab_estim, 1, CADO_MPI_SIZE_T, root, comm);
        file_bases = {0, tsize};
        totalsize = tsize;
    } else {
        size_t dummy;
        size_t hdrbytes;
        char line[ABFILE_MAX_LINE_LENGTH];
        if (rank == root) {
            FILE * file = fopen(fname.c_str(), "r");
            char * xx = fgets(line, sizeof(line), file);
            DIE_ERRNO_DIAG(xx == NULL, "fgets(%s)", fname.c_str());
            rc = sscanf(line, "AB %zu %zu", &nab, &dummy);
            DIE_ERRNO_DIAG(rc != 2, "parse(%s)", fname.c_str());
            hdrbytes = ftell(file);
            fclose(file);
        }
        MPI_Bcast(&hdrbytes, 1, CADO_MPI_SIZE_T, root, comm);

        sname_len = fname.size() + 20;
        for(nfiles = 0 ; ; nfiles++) {
            sname = prefix + ".prep." + std::to_string(depnum) + ".rel." + std::to_string(nfiles);
            rc = stat(sname.c_str(), sbuf);
            ASSERT_ALWAYS(rc == 0 || errno == ENOENT);
            ASSERT_ALWAYS(mpi_data_agrees(rc, comm));
            if (rc < 0) break;
            tsize += sbuf->st_size - hdrbytes;
        }
        ASSERT_ALWAYS(nfiles > 0);
        nab_estim = nab;
        file_bases.assign(nfiles + 1, 0);
        file_bases[0] = 0;
        for(int i = 0 ; i < nfiles ; i++) {
            sname = prefix + ".prep." + std::to_string(depnum) + ".rel." + std::to_string(i);
            rc = stat(sname.c_str(), sbuf);
            ASSERT_ALWAYS(rc == 0);
            ASSERT_ALWAYS(mpi_data_agrees(rc, comm));
            file_bases[i+1] = file_bases[i] + sbuf->st_size;
        }
        totalsize = file_bases[nfiles];
    }
}

int ab_source::openfile_internal()
{
    const char * s;
    if (nfiles == 0) {
        f = fopen(s = fname0.c_str(), "r");
        tpos = 0;
    } else {
        sname = prefix + ".prep." + std::to_string(depnum) + ".rel." + std::to_string(c);
        f = fopen(s = sname.c_str(), "r");
        if (f == NULL && errno == ENOENT)
            return 0;
        tpos = file_bases[c];
        cpos = 0;
        char header[80];
        char * rp = fgets(header, sizeof(header), f);
        ASSERT_ALWAYS(rp);
    }
    DIE_ERRNO_DIAG(f == NULL, "fopen(%s)", s);
    cpos = ftell(f);
    tpos += cpos;
    return 1;
}

int ab_source::next(int64_t & a, uint64_t & b)
{
    if (f) {
        int rc MAYBE_UNUSED;
        char line[ABFILE_MAX_LINE_LENGTH];
        char * xx = fgets(line, sizeof(line), f);
        size_t cpos_new = ftell(f);
        if (xx) {
            if (nfiles == 0) {
                rc = sscanf(line, "%" SCNd64 " %" SCNu64, &a, &b);
                ASSERT(rc == 2);
            } else {
                int dummy;
                rc = sscanf(line, "%d %" SCNd64 " %" SCNu64, &dummy, &a, &b);
                ASSERT(rc == 3);
            }
            tpos += cpos_new - cpos;
            cpos = cpos_new;
            nab++;
            return 1;
        }
        fclose(f);
        f = nullptr;
        tpos += cpos_new - tpos;
        if (nfiles == 0)
            return 0;
        c++;
    }
    if (openfile_internal() == 0)
        return 0;
    return next(a, b);
}

void ab_source::move_afterpos(size_t offset)
{
    if (tpos >= offset)
        return;
    ASSERT_ALWAYS(f == nullptr);

    FATAL_ERROR_CHECK(offset >= totalsize,
            "attempt to seek beyond end of files");

    size_t pre_offset = MAX(offset, 10) - 10;
    ASSERT_ALWAYS(pre_offset < offset);

    if (nfiles == 0) {
        openfile_internal();
        fseek(f, pre_offset, SEEK_SET);
    } else {
        for( ; file_bases[c+1] <= pre_offset ; c++) ;
        openfile_internal();
        fseek(f, pre_offset - file_bases[c], SEEK_SET);
    }
    char line[ABFILE_MAX_LINE_LENGTH];
    char * xx = fgets(line, sizeof(line), f);
    DIE_ERRNO_DIAG(xx == nullptr, "fgets(%s)", nfiles ? sname.c_str() : fname0.c_str());
    size_t cpos_new = ftell(f);
    tpos += cpos_new - cpos;
    cpos = cpos_new;
    for(int n_adjust = 0 ; tpos < offset ; n_adjust++) {
        FATAL_ERROR_CHECK(n_adjust > 10, "adjustment on the runaway");
        int64_t a;
        uint64_t b;
        int r = next(a, b);
        FATAL_ERROR_CHECK(r == 0, "adjustment failed");
    }
}
