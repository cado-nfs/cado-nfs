#include "cado.h" // IWYU pragma: keep
#include <stdexcept>
#include <ios>
#include <cstddef>      // for NULL
#include "cado_pipe_streambuf.hpp"
#include "cado_popen.h"

static const char*
ios_mode_to_stdio_mode(std::ios_base::openmode mode)
{
    switch (mode) {
        case std::ios_base::in:
            return "r";
        case std::ios_base::out:
            return "w";
            /* does not work. why ?
        case std::ios_base::in | std::ios_base::binary:
            return "rb";
        case std::ios_base::out | std::ios_base::binary:
            return "wb";
            */
        default:
            throw std::invalid_argument("open mode not understood");
    }
    return NULL;
}

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
