#include "cado.h" // IWYU pragma: keep

#include <stdexcept>
#include <ios>

#include "fd_streambuf.hpp"
#include "cado_pipe_streambuf.hpp"
#include "cado_popen.h"
#include "macros.h"

namespace {
const char*
ios_mode_to_stdio_mode(std::ios_base::openmode mode)
{
    /* See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=81665
     *
     * This is corrected in g++-15, not before.
     */
#if GNUC_VERSION_ATMOST(14,99,0)
    auto F = [](auto x) { return int(x); };
#else
    auto F = [](auto x) { return x; };
#endif

    switch (F(mode)) {
        case F(std::ios_base::in):
            return "r";
        case F(std::ios_base::out):
            return "w";
        case F(std::ios_base::in | std::ios_base::binary):
            return "rb";
        case F(std::ios_base::out | std::ios_base::binary):
            return "wb";
        default:
            throw std::invalid_argument("open mode not understood");
    }
    return nullptr;
}
} // namespace

cado_pipe_streambuf::cado_pipe_streambuf(const char* command, const char* mode)
  : fd_streambuf(cado_fd_popen(command, mode))
{}

cado_pipe_streambuf::cado_pipe_streambuf(const char* command, std::ios_base::openmode mode)
  : fd_streambuf(cado_fd_popen(command, ios_mode_to_stdio_mode(mode)))
{}

cado_pipe_streambuf::~cado_pipe_streambuf()
{
    if (fd_ < 0)
        return;
    cado_fd_pclose(fd_);
    fd_ = -1;
}

void
cado_pipe_streambuf::close()
{
    sync();
    cado_fd_pclose(fd_);
    fd_ = -1;
};
