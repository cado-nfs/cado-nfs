#ifndef CADO_UTILS_FD_STREAMBUF_HPP
#define CADO_UTILS_FD_STREAMBUF_HPP

#include <cstddef>

#include <streambuf>
#include <vector>
#include <unistd.h>

#include <sys/types.h>

class fd_streambuf : public std::streambuf {
public:
    explicit fd_streambuf(int fd, size_t bufferSize = 4096)
        : fd_(fd)
        , bufferSize_(bufferSize)
    {
        readBuf.resize(bufferSize_);
        writeBuf.resize(bufferSize_);

        setp(writeBuf.data(), writeBuf.data() + writeBuf.size());
        setg(nullptr, nullptr, nullptr);
    }

    ~fd_streambuf() override { do_flush(); }
    fd_streambuf(const fd_streambuf&) = delete;
    fd_streambuf& operator=(const fd_streambuf&) = delete;
    fd_streambuf(fd_streambuf&& other) noexcept = default;
    fd_streambuf& operator=(fd_streambuf&& other) noexcept = default;

    bool flush() {
        return do_flush() == 0;
    }

protected:
    int sync() override { return do_flush(); }

    int_type underflow() override {
        if (fd_ < 0) return traits_type::eof();

        const ssize_t n = ::read(fd_, readBuf.data(), bufferSize_);
        if (n <= 0) return traits_type::eof();

        setg(readBuf.data(), readBuf.data(), readBuf.data() + n);
        return traits_type::to_int_type(*gptr());
    }

    int_type overflow(int_type c) override {
        if (sync() == -1) return traits_type::eof();

        if (c != traits_type::eof()) {
            *pptr() = traits_type::to_char_type(c);
            pbump(1);
        }
        return traits_type::not_eof(c);
    }

private:
    int do_flush() {
        if (fd_ < 0 || pptr() == pbase()) return 0;

        const size_t n = pptr() - pbase();
        if (::write(fd_, pbase(), n) != static_cast<ssize_t>(n)) return -1;

        setp(writeBuf.data(), writeBuf.data() + writeBuf.size());
        return 0;
    }
    int fd_;    // we don't own this fd!
    size_t bufferSize_;
    std::vector<char> readBuf;
    std::vector<char> writeBuf;
};

#endif /* CADO_UTILS_FD_STREAMBUF_HPP */
