#ifndef CADO_UTILS_FD_STREAMBUF_HPP
#define CADO_UTILS_FD_STREAMBUF_HPP

#include <streambuf>
#include <cstdio>       // reuse BUFSIZ from stdio.

class fd_streambuf : public std::streambuf
{
    public:
    static constexpr const size_t bufferSize = BUFSIZ;
    explicit fd_streambuf(int d = -1);
    ~fd_streambuf() override;
    int fd() const { return fd_; }
    // void fd(int v);
    void close(); // { fd(-1); }

    protected:
    int_type underflow() override;
    int_type overflow(int_type c) override;
    int sync() override;

    private:
    int fd_;
    char *readBuf, *writeBuf;
    friend class cado_pipe_streambuf;

    public:
    fd_streambuf(fd_streambuf const &) = delete;
    fd_streambuf(fd_streambuf &&) = delete;
    fd_streambuf& operator=(fd_streambuf const &) = delete;
    fd_streambuf& operator=(fd_streambuf &&) = delete;
};

#endif /* CADO_UTILS_FD_STREAMBUF_HPP */
