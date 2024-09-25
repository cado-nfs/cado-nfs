#include "cado.h" // IWYU pragma: keep
#include <streambuf>
#include <new>
#include <unistd.h>
#include "fd_streambuf.hpp"
#include "macros.h"

// constexpr const size_t fd_streambuf::bufferSize;

fd_streambuf::fd_streambuf(int d)
    : fd_(d)
    , readBuf(nullptr)
    , writeBuf(nullptr)
{
}

fd_streambuf::~fd_streambuf()
{
    if (writeBuf) { delete[] writeBuf; writeBuf = nullptr; }
    if (readBuf) { delete[] readBuf; readBuf = nullptr; }
    if (fd_ < 0) return;
    close();
}


void fd_streambuf::close()
{
    /* sync() is virtual. We are only assuming that we're resolving at
     * the current level, and we want to make that explicit.
     */
    fd_streambuf::sync();
    ::close(fd_);

    setg(nullptr, nullptr, nullptr);
    setp(nullptr, nullptr);

    delete [] readBuf;
    readBuf = nullptr;

    delete [] writeBuf;
    writeBuf = nullptr;

    fd_ = -1;
}

fd_streambuf::int_type fd_streambuf::underflow()
{
    ASSERT(gptr() == egptr());

    if (fd_ < 0) {
        return traits_type::eof();
    }

    if (readBuf) {
        *readBuf = *(gptr() - 1);
    } else {
        readBuf = new (std::nothrow) char[bufferSize];
        if (!readBuf) {
            return traits_type::eof();
        }
    }

    auto nRead = read(fd_, readBuf + 1, bufferSize - 1);
    if (nRead <= 0) {
        return traits_type::eof();
    }
    setg(readBuf, readBuf + 1, readBuf + 1 + nRead);
    return traits_type::to_int_type(*gptr());
}

fd_streambuf::int_type fd_streambuf::overflow(int_type c)
{
    ASSERT(pptr() == epptr());

    if (fd_ < 0) {
        return traits_type::eof();
    }

    if (!writeBuf) {
        writeBuf = new (std::nothrow) char[bufferSize];
        if (!writeBuf) {
            return traits_type::eof();
        }
    }

    if (c == traits_type::eof() || sync() == -1) {
        return traits_type::eof();
    }

    *pptr() = traits_type::to_char_type(c);
    pbump(1);
    return c;
}

int fd_streambuf::sync()
{
    if (fd_ < 0 || !writeBuf) {
        return 0;
    }

    char *p = pbase();
    while (p < pptr()) {
        auto nWritten = write(fd_, p, pptr() - p);
        if (nWritten <= 0) {
            return -1;
        }
        p += nWritten;
    }

    setp(writeBuf, writeBuf + bufferSize);
    return 0;
}
