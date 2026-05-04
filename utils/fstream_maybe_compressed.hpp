#ifndef UTILS_FSTREAM_COMPRESSED_HPP_
#define UTILS_FSTREAM_COMPRESSED_HPP_

#include <istream>
#include <ostream>
#include <memory>
#include <fstream>
#include <string>
#include <ios>

#include "cado_pipe_streambuf.hpp"

class streambase_maybe_compressed : virtual public std::ios {
    protected:
    std::unique_ptr<cado_pipe_streambuf> pbuf;
    std::unique_ptr<std::filebuf> fbuf;
    std::string orig_name;
    std::string tempname;
    public:
    bool is_open() const { return rdbuf(); }
    /* ofstream_maybe_compressed and ifstream_maybe_compressed must be
     * default-constructible */
    streambase_maybe_compressed() = default;
    streambase_maybe_compressed(std::string const & name, std::ios_base::openmode mode);
    /* Note that in output mode, the file will first be created with a
     * temp name, and eventually only the dtor will move it from that
     * temp name to the final location.
     * (this behaviour might be system-dependent).
     */
    ~streambase_maybe_compressed() override;
    void open(std::string const &, std::ios_base::openmode mode);
    void close();
    bool is_pipe() const { return pbuf.get(); }
    streambase_maybe_compressed(streambase_maybe_compressed const &) = delete;
    streambase_maybe_compressed(streambase_maybe_compressed &&) = delete;
    streambase_maybe_compressed& operator=(streambase_maybe_compressed const &) = delete;
    streambase_maybe_compressed& operator=(streambase_maybe_compressed &&) = delete;
};

template <class charT, class Traits = std::char_traits<charT> >
class basic_ifstream_maybe_compressed : public streambase_maybe_compressed, public std::basic_istream<charT, Traits> {
public:
    basic_ifstream_maybe_compressed(std::string const & name,
            std::ios_base::openmode mode)
        : streambase_maybe_compressed(name, mode)
        , std::basic_istream<charT, Traits>(rdbuf())
    {}
    basic_ifstream_maybe_compressed() = default;
    explicit basic_ifstream_maybe_compressed(std::string const & name)
        : basic_ifstream_maybe_compressed(name, std::ios_base::in)
    {}
    void open(std::string const & name, std::ios_base::openmode mode = std::ios_base::in) {
        streambase_maybe_compressed::open(name, mode);
    }
};

template <class charT, class Traits = std::char_traits<charT> >
class basic_ofstream_maybe_compressed : public streambase_maybe_compressed, public std::basic_ostream<charT, Traits> {
public:
    basic_ofstream_maybe_compressed(std::string const & name,
            std::ios_base::openmode mode)
        : streambase_maybe_compressed(name, mode)
        , std::basic_ostream<charT, Traits>(rdbuf())
    {}
    basic_ofstream_maybe_compressed() = default;
    explicit basic_ofstream_maybe_compressed(std::string const & name)
        : basic_ofstream_maybe_compressed(name, std::ios_base::out)
    {}
    void open(std::string const & name, std::ios_base::openmode mode = std::ios_base::out) {
        streambase_maybe_compressed::open(name, mode);
    }
};

using ifstream_maybe_compressed = basic_ifstream_maybe_compressed<char>;
using ofstream_maybe_compressed = basic_ofstream_maybe_compressed<char>;

#endif	/* UTILS_FSTREAM_COMPRESSED_HPP_ */
