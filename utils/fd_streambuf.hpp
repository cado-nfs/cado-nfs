#ifndef CXX_FDSTREAMBUF_HPP_
#define CXX_FDSTREAMBUF_HPP_

#include <streambuf>
#include <cstdio>       // reuse BUFSIZ from stdio.

class fd_streambuf : public std::streambuf
{
    public:
    static constexpr const size_t bufferSize = BUFSIZ;
    fd_streambuf(int d = -1);
    ~fd_streambuf();
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
};

#endif /* CXX_FDSTREAMBUF_HPP_ */
