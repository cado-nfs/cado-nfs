#ifndef LINGEN_SUBSTEP_SCHEDULE_HPP_
#define LINGEN_SUBSTEP_SCHEDULE_HPP_

#include <array>
#include <istream>

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

    lingen_substep_schedule() : shrink0(1), shrink2(1), batch {{1,1,1}}  {}
    lingen_substep_schedule(lingen_substep_schedule const&) = default;

    std::ostream& serialize(std::ostream& os) const
    {
        switch(fft_type) {
            case FFT_NONE: os << " " << io_token_fft_none; break;
            case FFT_FLINT: os << " " << io_token_fft_flint; break;
            case FFT_CANTOR: os << " " << io_token_fft_cantor; break;
            case FFT_TERNARY: os << " " << io_token_fft_ternary; break;
        }
        os << " " << io_token_shrink << " " << shrink0 << " " << shrink2;
        os << " " << io_token_batch << " " << batch[0] << " " << batch[1] << " " << batch[2];
        return os;
    }

    std::istream& unserialize(std::istream& is)
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

    bool operator<(lingen_substep_schedule const & o) const {
        if (fft_type < o.fft_type) return true;
        if (fft_type > o.fft_type) return false;
        if (shrink0 < o.shrink0) return true;
        if (shrink0 > o.shrink0) return false;
        if (shrink2 < o.shrink2) return true;
        if (shrink2 > o.shrink2) return false;
        if (batch < o.batch) return true;
        if (batch > o.batch) return false;
        return false;
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

#endif	/* LINGEN_SUBSTEP_SCHEDULE_HPP_ */
