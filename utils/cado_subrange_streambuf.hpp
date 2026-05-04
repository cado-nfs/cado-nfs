#ifndef UTILS_CADO_SUBRANGE_STREAMBUF_HPP_
#define UTILS_CADO_SUBRANGE_STREAMBUF_HPP_

#include <iosfwd>
#include <ios>
#include <algorithm>
#include <limits>

namespace cado {
class subrange_streambuf : public std::streambuf {
    static constexpr off_type eofpos = -1;
public:
    subrange_streambuf(std::streambuf* source, std::streampos start, std::streampos end)
        : source_(source)
        , start_(start)
        , end_(end)
        , length_(end - start)
    {
        // Ensure the range is valid
        length_ = std::max(length_, static_cast<off_type>(0));
    }

protected:
    // Called when the internal buffer is exhausted
    int_type underflow() override {
        if (current_pos_ >= length_) {
            return traits_type::eof();
        }

        // Seek the underlying source buffer to the correct absolute position
        if (source_->pubseekpos(start_ + current_pos_, std::ios_base::in) == std::streampos(eofpos)) {
            return traits_type::eof();
        }

        return source_->sgetc();
    }

    // Move to the next character
    int_type uflow() override {
        const int_type ch = underflow();
        if (ch != traits_type::eof()) {
            current_pos_ += 1;
            source_->sbumpc(); // Advance the source buffer
        }
        return ch;
    }

    // Handles seeking within the fragment
    pos_type seekoff(off_type off, std::ios_base::seekdir dir, std::ios_base::openmode) override {
        off_type new_pos;
        if (dir == std::ios_base::beg)      new_pos = off;
        else if (dir == std::ios_base::cur) new_pos = current_pos_ + off;
        else if (dir == std::ios_base::end) new_pos = length_ + off;
        else return {eofpos};

        if (new_pos < 0 || new_pos > length_) return {eofpos};
        
        current_pos_ = new_pos;
        return {current_pos_};
    }

    pos_type seekpos(pos_type sp, std::ios_base::openmode which) override {
        return seekoff(sp, std::ios_base::beg, which);
    }

    // Optimization: Bulk read
    std::streamsize xsgetn(char* s, std::streamsize n) override {
        const std::streamsize remaining = length_ - current_pos_;
        const std::streamsize to_read = std::min(n, remaining);

        if (to_read <= 0) return 0;

        source_->pubseekpos(start_ + current_pos_, std::ios_base::in);
        const std::streamsize actually_read = source_->sgetn(s, to_read);
        
        current_pos_ += actually_read;
        return actually_read;
    }

private:
    std::streambuf* source_;
    std::streampos start_;
    std::streampos end_;
    off_type current_pos_ = 0;
    off_type length_ = std::numeric_limits<off_type>::max();
};
} /* namespace cado */

#endif	/* UTILS_CADO_SUBRANGE_STREAMBUF_HPP_ */
