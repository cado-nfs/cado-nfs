#ifndef CADO_LINGEN_SUBSTEP_SCHEDULE_HPP
#define CADO_LINGEN_SUBSTEP_SCHEDULE_HPP

#include <array>
#include <istream>
#include <ostream>
#include <sstream>

#include "lingen_fft_select.hpp"

struct lingen_substep_schedule {
    /* output characteristics -- the ones we have to choose */

    /* batch, shrink0, shrink2. Used to be static parameters, now they're
     * dynamic.
     *
     * batch: with batch=b, schedule b times more simultaneous
     *  transforms, so that we can have more parallelism. Costs b times
     *  more RAM
     * shrink0: divide the number of transforms on dimension 0 (m/n for
     *  MP, (m+n)/r for MUL) so as to save memory
     * shrink2: divide the number of transforms on dimension 2 (m+n)/r
     *  for both MP and MUL) so as to save memory
     */

    enum fft_type_t { FFT_NONE, FFT_FLINT, FFT_CANTOR, FFT_TERNARY } fft_type;

    /* shrink0 and shrink2 are important. Here, we restrict our operation
     * to shrink0*shrink2 multiplications of matrices of size
     * (n0/shrink0) times (n2/shrink2).
     */
    unsigned int shrink0;
    unsigned int shrink2;

    /* We may batch in all three dimensions. At most, batch is {n0, n1,
     * n2}. More generally, if it is {b0,b1,b2}, then the amount of RAM
     * spent for transforms is r * b1 * (b0 + b2) transforms. This brings
     * the benefit of computing b0*b1 transforms simultaneously on A, and
     * b1*b2 transforms simultaneously on B.
     *
     * The outer loop, of size n1/r, is processed b1 steps at a time.
     */
    std::array<unsigned int, 3> batch;
    private:
    static constexpr const char * io_token_shrink = "shrink";
    static constexpr const char * io_token_batch = "batch";
    static constexpr const char * io_token_fft_none = "fft=none";
    static constexpr const char * io_token_fft_flint = "fft=flint";
    static constexpr const char * io_token_fft_cantor = "fft=cantor";
    static constexpr const char * io_token_fft_ternary = "fft=ternary";
    public:

    lingen_substep_schedule() : fft_type(FFT_NONE), shrink0(1), shrink2(1), batch {{1,1,1}}  {}
    /*
     * well, defaults are fine, it's probably best to say nothing
     *
    lingen_substep_schedule(lingen_substep_schedule const&) = default;
    lingen_substep_schedule(lingen_substep_schedule &&) = default;
    lingen_substep_schedule& operator=(lingen_substep_schedule const&) = default;
    lingen_substep_schedule& operator=(lingen_substep_schedule &&) = default;
    */

    static const char * fft_name(fft_type_t fft_type) {
        switch(fft_type) {
            case FFT_NONE:   return io_token_fft_none;
            case FFT_FLINT:  return io_token_fft_flint;
            case FFT_CANTOR: return io_token_fft_cantor;
            case FFT_TERNARY:return io_token_fft_ternary;
            default:
              throw std::runtime_error("invalid data (fft_type)");
        }
    }
    const char * fft_name() const { return fft_name(fft_type); }

    static std::ostream& fft_type_serialize(std::ostream& os, fft_type_t fft_type)
    {
        return os << fft_name(fft_type);
    }

    std::ostream& serialize(std::ostream& os) const
    {
        os << " ";
        fft_type_serialize(os, fft_type);
        os << " " << io_token_shrink << " " << shrink0 << " " << shrink2;
        os << " " << io_token_batch << " " << batch[0] << " " << batch[1] << " " << batch[2];
        return os;
    }

    std::string serialize() const {
        std::stringstream ss;
        serialize(ss);
        return ss.str();
    }

    static std::istream& fft_type_unserialize(std::istream& is, fft_type_t & fft_type)
    {
        std::string s;
        is >> s;
        if (s == io_token_fft_none) {
            fft_type = FFT_NONE;
        } else if (s == io_token_fft_flint) {
            fft_type = FFT_FLINT;
        } else if (s == io_token_fft_cantor) {
            fft_type = FFT_CANTOR;
        } else if (s == io_token_fft_ternary) {
            fft_type = FFT_TERNARY;
        } else {
            is.setstate(std::ios::failbit);
            return is;
        }
        return is;
    }
    std::istream& unserialize(std::istream& is)
    {
        if (!fft_type_unserialize(is, fft_type))
            return is;
        std::string s;
        is >> s;
        if (s != io_token_shrink) {
            is.setstate(std::ios::failbit);
            return is;
        }
        is >> shrink0 >> shrink2;
        is >> s;
        if (s != io_token_batch) {
            is.setstate(std::ios::failbit);
            return is;
        }
        is >> batch[0] >> batch[1] >> batch[2];
        return is;
    }

    auto operator<=>(lingen_substep_schedule const & o) const {
        if (auto c = fft_type <=> o.fft_type; c != 0) return c;
        if (auto c = shrink0 <=> o.shrink0; c != 0) return c;
        if (auto c = shrink2 <=> o.shrink2; c != 0) return c;
        if (auto c = batch <=> o.batch; c != 0) return c;
        return std::strong_ordering::equal;
    }
    bool operator==(lingen_substep_schedule const & o) const {
        return fft_type == o.fft_type
            && shrink0 == o.shrink0
            && shrink2 == o.shrink2
            && batch == o.batch;
    }
    bool check() const {
        return shrink0 && shrink2 && batch[0] && batch[1] && batch[2];
    }
};

template<typename T> struct encode_fft_type_details {
    static constexpr lingen_substep_schedule::fft_type_t value = lingen_substep_schedule::FFT_NONE;
};
#ifndef LINGEN_BINARY
template<> struct encode_fft_type_details<fft_transform_info> {
    static constexpr lingen_substep_schedule::fft_type_t value = lingen_substep_schedule::FFT_FLINT;
};
#else
template<> struct encode_fft_type_details<gf2x_cantor_fft_info> {
    static constexpr lingen_substep_schedule::fft_type_t value = lingen_substep_schedule::FFT_CANTOR;
};
template<> struct encode_fft_type_details<gf2x_ternary_fft_info> {
    static constexpr lingen_substep_schedule::fft_type_t value = lingen_substep_schedule::FFT_TERNARY;
};
#endif
template<typename T>
inline constexpr lingen_substep_schedule::fft_type_t encode_fft_type = encode_fft_type_details<T>::value;


#endif	/* LINGEN_SUBSTEP_SCHEDULE_HPP_ */
