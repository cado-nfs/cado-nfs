#ifndef UTILS_FSTREAM_COMPRESSED_HPP_
#define UTILS_FSTREAM_COMPRESSED_HPP_


#include <istream>
#include <ostream>
#include <memory>
#include <fstream>      // std::filebuf // IWYU pragma: keep
#include <string>
#include <ios>
#include <streambuf>

class cado_pipe_streambuf;

class streambase_maybe_compressed : virtual public std::ios {
    bool pipe = false;
    protected:
    std::unique_ptr<cado_pipe_streambuf> pbuf;
    std::unique_ptr<std::filebuf> fbuf;
    std::streambuf * buf = nullptr;
    std::string orig_name;
    std::string tempname;
    public:
    /* I don't think that we need a default ctor, do we ? */
    streambase_maybe_compressed(std::string const & name, std::ios_base::openmode mode);
    /* Note that in output mode, the file will first be created with a
     * temp name, and eventually only the dtor will move it from that
     * temp name to the final location.
     * (this behaviour might be system-dependent).
     */
    ~streambase_maybe_compressed() override;
    void open(std::string const &, std::ios_base::openmode mode);
    void close();
    bool is_pipe() const { return pipe; }
    streambase_maybe_compressed(streambase_maybe_compressed const &) = delete;
    streambase_maybe_compressed(streambase_maybe_compressed &&) = delete;
    streambase_maybe_compressed& operator=(streambase_maybe_compressed const &) = delete;
    streambase_maybe_compressed& operator=(streambase_maybe_compressed &&) = delete;
};

template <class charT, class Traits = std::char_traits<charT> >
class basic_ifstream_maybe_compressed : public streambase_maybe_compressed, public std::basic_istream<charT, Traits> {
public:
    explicit basic_ifstream_maybe_compressed(std::string const & name)
        : streambase_maybe_compressed(name, std::ios::in)
        , std::basic_istream<charT, Traits>(buf)
    {}
    void open(std::string const & name) {
        streambase_maybe_compressed::open(name, std::ios::in);
    }
};

template <class charT, class Traits = std::char_traits<charT> >
class basic_ofstream_maybe_compressed : public streambase_maybe_compressed, public std::basic_ostream<charT, Traits> {
public:
    explicit basic_ofstream_maybe_compressed(std::string const & name)
        : streambase_maybe_compressed(name, std::ios::out)
        , std::basic_ostream<charT, Traits>(buf)
    {}
    void open(std::string const & name) {
        streambase_maybe_compressed::open(name, std::ios::out);
    }
};

using ifstream_maybe_compressed = basic_ifstream_maybe_compressed<char>;
using ofstream_maybe_compressed = basic_ofstream_maybe_compressed<char>;

#endif	/* UTILS_FSTREAM_COMPRESSED_HPP_ */
