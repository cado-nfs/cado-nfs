#ifndef CADO_PIPE_STREAMBUF_HPP_
#define CADO_PIPE_STREAMBUF_HPP_

#include <ios>
#include "fd_streambuf.hpp"

class cado_pipe_streambuf : public fd_streambuf {
    public:
    cado_pipe_streambuf(const char * command, const char * mode);
    cado_pipe_streambuf(const char * command, std::ios_base::openmode mode);
    ~cado_pipe_streambuf();
    void close();

//     protected:
//     int_type underflow() override { return fd_streambuf::underflow(); }
//     int_type overflow(int_type c) override { return fd_streambuf::overflow(c); }
//     int sync() override { return fd_streambuf::sync(); }
// 
};


#endif	/* CADO_PIPE_STREAMBUF_HPP_ */
