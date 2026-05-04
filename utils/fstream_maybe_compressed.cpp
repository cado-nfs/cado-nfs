#include "cado.h" // IWYU pragma: keep

#include <cstdio>

#include <ios>
#include <stdexcept>
#include <string>
#include <memory>
#include <fstream>

#include <unistd.h>

#include "fmt/base.h"
#include "fmt/format.h"

#include "macros.h"
#include "cado_pipe_streambuf.hpp"
#include "gzip.h"
#include "fstream_maybe_compressed.hpp"

streambase_maybe_compressed::streambase_maybe_compressed(std::string const & name, std::ios_base::openmode mode)
{
    open(name, mode);
}

void streambase_maybe_compressed::open(std::string const & name_arg, std::ios_base::openmode mode)
{
    std::string name = name_arg;
    orig_name = name;
    if (mode & std::ios_base::out) {
        // fmtlib's fmt::format oddly mentions that it can throw a format
        // error, while its constexpr nature should be able to mark it as
        // impossible.
        // coverity[exception_thrown]
        tempname = fmt::format("{}.tmp.{}", name, getpid());
        name = tempname;
        /* make sure that the tempname is writable. This way, we can
         * point to /dev/stdout, for instance. It's still a bit fragile
         * because of TOCTOU, but well, it's better than nothing.
         *
         * Note that it's something we can't do with the C version
         * fopen_maybe_compressed.
         */
        std::ofstream dummy(tempname, std::ios::out);
        if (dummy.is_open()) {
            dummy.close();
            unlink(tempname.c_str());
        } else {
            name = orig_name;
            tempname.clear();
        }
    }

    pbuf.reset();
    fbuf.reset();

    if (mode & std::ios_base::in && access(name.c_str(), R_OK) != 0)
        throw std::runtime_error("cannot open file for reading");
    /* creating is ok, of course
    if (mode & std::ios_base::out && access(name, W_OK) != 0)
        throw std::runtime_error("cannot open file for writing");
     */
    for(auto const & r : cado::details::supported_compression_formats) {
        if (!orig_name.ends_with(r.suffix)) continue;
        std::string command;
        if (mode & std::ios_base::in && *r.pfmt_in)
            command = fmt::format(fmt::runtime(r.pfmt_in), name);
        else if (mode & std::ios_base::out && *r.pfmt_out)
            command = fmt::format(fmt::runtime(r.pfmt_out), name);

        if (!command.empty()) {
            /* apparently popen() under Linux does not accept the 'b' modifier */
            pbuf = std::make_unique<cado_pipe_streambuf>(command.c_str(), mode);
            init(pbuf.get());
        } else {
            fbuf = std::make_unique<std::filebuf>();
            fbuf->open(name, mode);
            init(fbuf.get());
        }
        /* hmmm */
        return;
    }
    /* if exceptions() includes std::ios_base::badbit, then this will throw.
     * it might actually be a reasonable thing to do, but the truth is
     * that the exceptions mask is set only as a side-effect of something
     * (perhaps fmt::ostream_formatter but I'm not sure). In any case, we
     * do not want rdbuf() to point to a previously set buffer, and
     * having badbit is an _expected_ and correct post-condition of
     * closing.
     */
    exceptions(std::ios_base::goodbit);
    rdbuf(nullptr);
    this->setstate(ios_base::failbit);
};

void streambase_maybe_compressed::close()
{
    /* we need to close explicitly, otherwise the dtor won't flush.
     */
    if (pbuf) pbuf->close();
    if (fbuf) fbuf->close();
    /* see above */
    exceptions(std::ios_base::goodbit);
    rdbuf(nullptr);
    pbuf.reset();
    fbuf.reset();
    if (!tempname.empty()) {
        int const rc = rename(tempname.c_str(), orig_name.c_str());
        ASSERT_ALWAYS(rc == 0);
        tempname.clear();
    }
}

// we're in a dtor, exceptions can turn your computer into a coconut.
// yet we have an ASSERT_ALWAYS in close()
// coverity[exn_spec_violation]
streambase_maybe_compressed::~streambase_maybe_compressed()
{
    close();
}
