#ifndef CADO_LINGEN_IO_WRAPPERS_HPP
#define CADO_LINGEN_IO_WRAPPERS_HPP

#include "cado_config.h"

#include <cstdio>

#include <fstream>
#include <vector>
#include <memory>
#include <array>
#include <string>
#include <tuple>

#include <sys/types.h>

#include "gmp_aux.h"
#include "arith-hard.hpp"
#include "lingen_bigmatpoly.hpp"
#include "lingen_bw_dimensions.hpp"
#include "lingen_matpoly_select.hpp"
#include "macros.h"
#include "select_mpi.h"
#include "sha1.h"

template<bool is_binary>
struct bmstatus;

/* This layer intends to absorb some 1400 lines of lingen.cpp, in a
 * "pipes and fittings" approach. This should be more usable than the
 * very intricated bm_io layer that we currently have. The price that we
 * accept to pay in exchange is some caching matpoly's at some places
 * (those will also bring the benefit of maximizing the performance for
 * the binary case, and might even do good in the prime field case either).
 */
template<bool is_binary>
struct lingen_io_wrapper_base
{
    typename matpoly<is_binary>::arith_hard * ab;
    unsigned int nrows;
    unsigned int ncols;

    /* average_matsize is the size in bytes of the final storage target
     * of the stream. Typically that means the bytes on disk. The base
     * class has a method that determines this average matrix size from
     * the ab layer. This corresponds to the in-memory size.
     * However, streams that are tied to an actual file may
     * also possess the information of whether the file is in ascii or
     * not, in which case this function is overridden to represent the
     * number of bytes of the on-disk data.
     */
    virtual double average_matsize() const;

    /* This is an indicative size of the read batches that this layer
     * feels is appropriate for I/O. Currently the only means to control
     * that preferred window size is via the io_matpoly_block_size
     * parameter.
     *
     * A priori this is only dependent on the average matrix size, and we
     * can provide a default implementation in the base class.
     */
    virtual unsigned int preferred_window() const;

    lingen_io_wrapper_base(typename matpoly<is_binary>::arith_hard * ab, unsigned int nrows, unsigned int ncols)
      : ab(ab)
      , nrows(nrows)
      , ncols(ncols)
    {}

    protected:
    ~lingen_io_wrapper_base() {}
};

template<bool is_binary>
struct lingen_input_wrapper_base : public lingen_io_wrapper_base<is_binary>
{
    using lingen_io_wrapper_base<is_binary>::nrows;
    using lingen_io_wrapper_base<is_binary>::ncols;
    using lingen_io_wrapper_base<is_binary>::ab;
    using lingen_io_wrapper_base<is_binary>::average_matsize;
    lingen_input_wrapper_base(typename matpoly<is_binary>::arith_hard * ab,
                              unsigned int nrows,
                              unsigned int ncols)
      : lingen_io_wrapper_base<is_binary>(ab, nrows, ncols)
    {}
    /* The input source is regarded as a stream, with no possibility of
     * random access. Here, k0 and k1 refer to positions within the
     * destination matpoly. As far as the source is concerned, only
     * the k1-k0 "next" coefficient matrices will be fetched.
     *
     * The concrete implementation of this class _may_ keep track of a
     * position with respect to the concrete source, but that bit of
     * information is not exposed in the interface.
     */
    virtual ssize_t read_to_matpoly(matpoly<is_binary>& dst,
                                    unsigned int k0,
                                    unsigned int k1) = 0;
    virtual size_t guessed_length() const = 0;

    virtual ~lingen_input_wrapper_base() {}
};

template<bool is_binary>
struct lingen_output_wrapper_base : public lingen_io_wrapper_base<is_binary>
{
    using lingen_io_wrapper_base<is_binary>::nrows;
    using lingen_io_wrapper_base<is_binary>::ncols;
    using lingen_io_wrapper_base<is_binary>::ab;

    lingen_output_wrapper_base(typename matpoly<is_binary>::arith_hard * ab,
                               unsigned int nrows,
                               unsigned int ncols)
      : lingen_io_wrapper_base<is_binary>(ab, nrows, ncols)
    {}
    /* The output source is regarded as a stream, with no possibility of
     * random access. Here, k0 and k1 refer to positions within the
     * source matpoly. As far as the destination is concerned, we will
     * only store k1-k0 coefficient matrices "after" the last one stored.
     *
     * The concrete implementation of this class _may_ keep track of a
     * position with respect to the concrete destination, but that bit of
     * information is not exposed in the interface.
     */
    virtual ssize_t write_from_matpoly(matpoly<is_binary> const& src,
                                       unsigned int k0,
                                       unsigned int k1) = 0;

    virtual ~lingen_output_wrapper_base() {}
};

template<bool is_binary>
class lingen_file_input : public lingen_input_wrapper_base<is_binary>
{
    FILE* f;
    std::string filename;
    bool ascii;
    unsigned int length_hint;
    void open_file();
    void close_file();

    public:
    lingen_file_input(typename matpoly<is_binary>::arith_hard * ab,
                      unsigned int nrows,
                      unsigned int ncols,
                      std::string filename,
                      bool ascii,
                      unsigned int length_hint)
      : lingen_input_wrapper_base<is_binary>(ab, nrows, ncols)
      , filename(std::move(filename))
      , ascii(ascii)
      , length_hint(length_hint)
    {
        open_file();
    }
    unsigned int preferred_window() const override;
    ~lingen_file_input() override { close_file(); }
    lingen_file_input(lingen_file_input const&) = delete;
    double average_matsize() const override;
    size_t guessed_length() const override;
    ssize_t read_to_matpoly(matpoly<is_binary>& dst,
                            unsigned int k0,
                            unsigned int k1) override;
};

template<bool is_binary>
struct lingen_random_input : public lingen_input_wrapper_base<is_binary>
{
    using lingen_input_wrapper_base<is_binary>::average_matsize;
    cxx_gmp_randstate & rstate;
    size_t next_src_k = 0;
    size_t length;
    lingen_random_input(typename matpoly<is_binary>::arith_hard * ab,
                        unsigned int nrows,
                        unsigned int ncols,
                        cxx_gmp_randstate & rstate,
                        size_t length)
      : lingen_input_wrapper_base<is_binary>(ab, nrows, ncols)
      , rstate(rstate)
      , length(length)
    {}

    unsigned int preferred_window() const override;
    size_t guessed_length() const override { return length; }

    ssize_t read_to_matpoly(matpoly<is_binary>& dst,
                            unsigned int k0,
                            unsigned int k1) override;
};

template<bool is_binary>
struct lingen_F0 : protected bw_dimensions<is_binary>
{
    public:
        using bw_dimensions<is_binary>::m;
        using bw_dimensions<is_binary>::n;
        using bw_dimensions<is_binary>::nrhs;

    // pairs are (exponent, column number)
    std::vector<std::array<unsigned int, 2>> fdesc;
    /* This is initialized by lingen_E_from_A::initial_read
     * The default value that we put here is only for static analysis
     * whinings, really.
     */
    unsigned int t0 = UINT_MAX;
    void share(int root, MPI_Comm comm);
    lingen_F0(bw_dimensions<is_binary> const & d)
        : bw_dimensions<is_binary>(d)
    {}
    /* This method is valid only once F0 is completely filled
     *
     * this returns a pair (degree kA, column number jA) meaning that
     * column jE of E is actually X^{-kA}*column jA of A.
     */
    std::tuple<unsigned int, unsigned int> column_data_from_A(unsigned int jE) const;
    /* This one does the same, but with respect to the matrix with column
     * j shifted right if j >= nrhs
     */
    std::tuple<unsigned int, unsigned int> column_data_from_Aprime(unsigned int jE) const;
};

template<bool is_binary>
class lingen_E_from_A
  : public lingen_F0<is_binary>
  , public lingen_input_wrapper_base<is_binary>
{
    void initial_read();
    lingen_input_wrapper_base<is_binary>& A;

    /* cache contains degrees [cache_k0..cache_k1[ of A.
     * the span cache_k1 - cache_k0 is not fixed a priori.
     */
    matpoly<is_binary> cache;
    matpoly<is_binary> tail;
    unsigned int cache_k0 = 0;
    unsigned int cache_k1 = 0;
    // unsigned int next_src_k = 0; // same as cache_k1, in fact !


    void share(int root, MPI_Comm comm);
    void refresh_cache_upto(unsigned int k);

    static constexpr unsigned int simd = is_binary ? ULONG_BITS : 1;

    public:

    using lingen_io_wrapper_base<is_binary>::nrows;
    using lingen_io_wrapper_base<is_binary>::ncols;
    using lingen_F0<is_binary>::t0;
    using lingen_F0<is_binary>::m;
    using lingen_F0<is_binary>::n;
    using lingen_F0<is_binary>::nrhs;
    using lingen_F0<is_binary>::fdesc;

    lingen_E_from_A(bw_dimensions<is_binary> const & d,
            lingen_input_wrapper_base<is_binary> & A)
      : lingen_F0<is_binary>(d)
      , lingen_input_wrapper_base<is_binary>(A.ab, d.m, d.m + d.n)
      , A(A)
      , cache(A.ab, A.nrows, A.ncols, 0)
      , tail(A.ab, nrows, ncols, 0)
    {
        ASSERT_ALWAYS(A.nrows == d.m);
        ASSERT_ALWAYS(A.ncols == d.n);
        initial_read();
    }
    unsigned int preferred_window() const override {
        return A.preferred_window();
    }
    size_t guessed_length() const override {
        return A.guessed_length() - t0;
    }
    ssize_t read_to_matpoly(matpoly<is_binary>& dst,
                            unsigned int k0,
                            unsigned int k1) override;
};

template<typename matpoly_type>
class lingen_scatter : public lingen_output_wrapper_base<matpoly_type::is_binary>
{
    static constexpr bool is_binary = matpoly_type::is_binary;
    matpoly_type& E;
    unsigned int next_dst_k = 0;

    public:
    using lingen_output_wrapper_base<is_binary>::nrows;
    using lingen_output_wrapper_base<is_binary>::ncols;
    using lingen_output_wrapper_base<is_binary>::ab;
    explicit lingen_scatter(matpoly_type& E)
      : lingen_output_wrapper_base<is_binary>(E.ab, E.m, E.n)
      , E(E)
    {}
    ssize_t write_from_matpoly(matpoly<is_binary> const& src,
                               unsigned int k0,
                               unsigned int k1) override;
};

#if 0
/* yes, we must insist on extern template.
 * https://github.com/OpenKinect/libfreenect2/issues/157
 */
#ifdef LINGEN_BINARY
extern template class lingen_scatter<matpoly<true>>;
extern template class lingen_scatter<bigmatpoly<true>>;
#else
extern template class lingen_scatter<matpoly<false>>;
extern template class lingen_scatter<bigmatpoly<false>>;
#endif
#endif

#if 0
template<typename matpoly_type>
class shared_or_common_size {}

template<>
class shared_or_common_size<matpoly> {
    size_t shared;
public:
    shared_or_common_size(matpoly const &);
    inline size_t get_size(matpoly const &) const { return shared; }
};
template<>
class shared_or_common_size<bigmatpoly> {
    shared_or_common_size(bigmatpoly const &) {}
    inline size_t get_size(bigmatpoly const & pi) const { return pi.get_size(); }
};
#endif

template<typename matpoly_type>
class lingen_gather : public lingen_input_wrapper_base<matpoly_type::is_binary>
                      // , private shared_or_common_size<matpoly_type>
{
    // typedef shared_or_common_size<matpoly_type> size_accessor;
    static constexpr bool is_binary = matpoly_type::is_binary;
    matpoly_type& pi;
    unsigned int next_src_k = 0;

    public:
    using lingen_input_wrapper_base<is_binary>::nrows;
    using lingen_input_wrapper_base<is_binary>::ncols;
    using lingen_input_wrapper_base<is_binary>::ab;

    lingen_gather(matpoly_type& pi)
      : lingen_input_wrapper_base<is_binary>(pi.ab, pi.m, pi.n)
      // , size_accessor(pi)
      , pi(pi)
    {}
    size_t guessed_length() const override {
        // return size_accessor::get_size(pi);
        return pi.get_size();
    }
    ssize_t read_to_matpoly(matpoly<is_binary>& dst,
                            unsigned int k0,
                            unsigned int k1) override;
};

template<typename matpoly_type>
class lingen_gather_reverse : public lingen_input_wrapper_base<matpoly_type::is_binary>
                              // , private shared_or_common_size<matpoly_type>
{
    static constexpr bool is_binary = matpoly_type::is_binary;
    // typedef shared_or_common_size<matpoly_type> size_accessor;
    matpoly_type& pi;
    /* Since the source is written in reverse order, it's a bit of a
     * misnomer, really. We're really referring to the count which is 0
     * when we're about to read the leading coefficient.
     */
    unsigned int next_src_k = 0;

    public:
    lingen_gather_reverse(matpoly_type& pi)
      : lingen_input_wrapper_base<is_binary>(pi.ab, pi.m, pi.n)
      // , size_accessor(pi)
      , pi(pi)
    {}
    inline size_t guessed_length() const override {
        // return size_accessor::get_size(pi);
        return pi.get_size();
    }
    ssize_t read_to_matpoly(matpoly<is_binary>& dst,
                            unsigned int k0,
                            unsigned int k1) override;
};

template<bool is_binary>
class lingen_F_from_PI
  : public lingen_F0<is_binary>
  , public lingen_input_wrapper_base<is_binary>
{
    public:
    using lingen_F0<is_binary>::nrhs;
    using lingen_F0<is_binary>::m;
    using lingen_F0<is_binary>::n;
    using lingen_io_wrapper_base<is_binary>::nrows;
    using lingen_io_wrapper_base<is_binary>::ncols;

    private:
    static constexpr unsigned int simd = is_binary ? ULONG_BITS : 1;
    lingen_input_wrapper_base<is_binary> & pi;
    matpoly<is_binary> cache;
    matpoly<is_binary> tail;
    unsigned int cache_k0 = 0;
    unsigned int cache_k1 = 0;
    // unsigned int next_src_k = 0; // same as cache_k1, in fact !

    matpoly<is_binary> rhs;
    struct sol_desc
    {
        unsigned int j;
        unsigned int shift;
    };
    std::vector<sol_desc> sols;
    matpoly<is_binary> recompute_rhs();
    void reorder_solutions();
    /* This returns (iF, s), such that the reversal of pi_{ipi, jpi}
     * contributes to entry (iF, jF), once shifted right by s.
     */
    std::tuple<unsigned int, unsigned int> get_shift_ij(unsigned int ipi, unsigned jF) const;
    public:
    lingen_F_from_PI(bmstatus<is_binary> const &, lingen_input_wrapper_base<is_binary>& pi, lingen_F0<is_binary> const& F0);
    using lingen_F0<is_binary>::t0;
    size_t guessed_length() const override {
        return pi.guessed_length() + t0;
    }
    ssize_t read_to_matpoly(matpoly<is_binary>& dst,
                            unsigned int k0,
                            unsigned int k1) override;
    void write_rhs(lingen_output_wrapper_base<is_binary> & Srhs);
};

template<bool is_binary>
class lingen_output_to_singlefile : public lingen_output_wrapper_base<is_binary>
{
    std::string filename;
    std::unique_ptr<std::ofstream> os;
    bool ascii;
    void open_file();
    bool done_open = false;

    public:
    lingen_output_to_singlefile(typename matpoly<is_binary>::arith_hard * ab,
                                unsigned int nrows,
                                unsigned int ncols,
                                std::string filename,
                                bool ascii = false)
      : lingen_output_wrapper_base<is_binary>(ab, nrows, ncols)
      , filename(std::move(filename))
      , ascii(ascii)
    {
    }

    lingen_output_to_singlefile(lingen_output_to_singlefile const&) = delete;

    ssize_t write_from_matpoly(matpoly<is_binary> const& src,
                               unsigned int k0,
                               unsigned int k1) override;
};

template<bool is_binary>
class lingen_output_to_splitfile : public lingen_output_wrapper_base<is_binary>
{
    std::string pattern;
    std::vector<std::ofstream> fw;
    bool ascii;
    void open_file();
    bool done_open = false;

    public:
    using lingen_output_wrapper_base<is_binary>::nrows;
    using lingen_output_wrapper_base<is_binary>::ncols;

    lingen_output_to_splitfile(typename matpoly<is_binary>::arith_hard * ab,
                               unsigned int nrows,
                               unsigned int ncols,
                               std::string pattern,
                               bool ascii = false);

    ssize_t write_from_matpoly(matpoly<is_binary> const& src,
                               unsigned int k0,
                               unsigned int k1) override;
};

/* This just prints the checksum to stdout on the dtor */
template<bool is_binary>
class lingen_output_to_sha1sum : public lingen_output_wrapper_base<is_binary>
{
    sha1_checksumming_stream f;
    std::string who;
    size_t written = 0;

    public:
    lingen_output_to_sha1sum(typename matpoly<is_binary>::arith_hard * ab,
                             unsigned int nrows,
                             unsigned int ncols,
                             std::string who)
      : lingen_output_wrapper_base<is_binary>(ab, nrows, ncols)
      , who(std::move(who))
    {}
    ~lingen_output_to_sha1sum() override;
    ssize_t write_from_matpoly(matpoly<is_binary> const& src,
                               unsigned int k0,
                               unsigned int k1) override;
};

template<bool is_binary>
void
pipe(lingen_input_wrapper_base<is_binary> & in,
     lingen_output_wrapper_base<is_binary> & out,
     const char * action, bool skip_trailing_zeros = false);


#endif /* LINGEN_IO_WRAPPERS_HPP_ */
