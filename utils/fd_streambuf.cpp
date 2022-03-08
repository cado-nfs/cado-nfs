#include "cado.h" // IWYU pragma: keep
#include <streambuf>
#include <new>
#include <unistd.h>
#include "fd_streambuf.hpp"
#include "macros.h"

constexpr const size_t fd_streambuf::bufferSize;

fd_streambuf::fd_streambuf(int d) : fd_(d), readBuf(0), writeBuf(0)
{
}

fd_streambuf::~fd_streambuf()
{
    if (writeBuf) { delete[] writeBuf; writeBuf = 0; }
    if (readBuf) { delete[] readBuf; readBuf = 0; }
    if (fd_ < 0) return;
    close();
}

/*
   void fd_streambuf::fd(int v)
   {
   if (fd_ == v) {
   return;
   }

   if (fd_ >= 0) {
   sync();
   ::close(fd_);
   }

   setg(0, 0, 0);
   setp(0, 0);

   delete [] readBuf;
   readBuf = 0;

   delete [] writeBuf;
   writeBuf = 0;

   fd_ = v;
   }
   */

void fd_streambuf::close()
{
    /* sync() is virtual. We are only assuming that we're resolving at
     * the current level, and we want to make that explicit.
     */
    fd_streambuf::sync();
    ::close(fd_);

    setg(0, 0, 0);
    setp(0, 0);

    delete [] readBuf;
    readBuf = 0;

    delete [] writeBuf;
    writeBuf = 0;

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

    int nRead = read(fd_, readBuf + 1, bufferSize - 1);
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
        int nWritten = write(fd_, p, pptr() - p);
        if (nWritten <= 0) {
            return -1;
        }
        p += nWritten;
    }

    setp(writeBuf, writeBuf + bufferSize);
    return 0;
}
