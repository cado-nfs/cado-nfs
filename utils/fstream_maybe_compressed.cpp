#include "cado.h" // IWYU pragma: keep

#include <ios>
#include <stdexcept>
#include <string>

#include <fcntl.h>
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
    init(buf);
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
    }

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
            buf = pbuf.get();
            pipe = true;
        } else {
            fbuf = std::make_unique<std::filebuf>();
            fbuf->open(name, mode);
            buf = fbuf.get();
            pipe = false;
        }
        /* hmmm */
        return;
    }
};

void streambase_maybe_compressed::close()
{
    if (pipe) pbuf->close();
    else fbuf->close();
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
