#ifndef CADO_PIPE_STREAMBUF_HPP
#define CADO_PIPE_STREAMBUF_HPP

#include <ios>
#include "fd_streambuf.hpp"

class cado_pipe_streambuf : public fd_streambuf {
    private:
    int fd_ = -1;   /* we own this. fd_streambuf is a copy */
    cado_pipe_streambuf(int fd, const char* command, const char* mode);
    public:
    cado_pipe_streambuf(const char * command, const char * mode);
    cado_pipe_streambuf(const char * command, std::ios_base::openmode mode);
    ~cado_pipe_streambuf() override;
    // void close();
    bool is_open() const { return fd_ >= 0; }

//     protected:
//     int_type underflow() override { return fd_streambuf::underflow(); }
//     int_type overflow(int_type c) override { return fd_streambuf::overflow(c); }
//     int sync() override { return fd_streambuf::sync(); }
// 
};


#endif	/* CADO_PIPE_STREAMBUF_HPP */
