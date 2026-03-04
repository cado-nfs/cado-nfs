#ifndef CADO_LAS_NORMS_HPP
#define CADO_LAS_NORMS_HPP

#include <cstdint>
#include <cstring>
#include <ostream>
#include <string>
#include <vector>
#include <array>

#include <gmp.h>
#include "cado_poly.h"
#include "las-config.hpp"
#include "las-qlattice.hpp"
#include "las-siever-config.hpp"
#include "logapprox.hpp"
#include "polynomial.hpp"
#include "macros.h"
#include "mpz_poly.h"
#include "sieve-methods.hpp"

struct special_q; // IWYU pragma: keep

/* Only relevant with --adjust-strategy 2 */
#define ADJUST_STRATEGY2_MIN_SQUEEZE 0
#define ADJUST_STRATEGY2_MAX_SQUEEZE 3

double get_maxnorm_rectangular (polynomial<double> const & src_poly, double X, double Y);

struct lognorm_base {/*{{{*/
    int logI;
    uint32_t J;

    unsigned char bound; /* A sieve array entry is a sieve survivor if it is
                            at most "bound" on each side */
    protected:
    double maxlog2;      /* Bound on the log in base 2. This is
                            intermediary data, really. */
    public:
    double get_maxlog2() const { return maxlog2; }
    cxx_mpz_poly fij;  /* coefficients of F(a0*i+a1*j, b0*i+b1*j)
                        * (divided by q on the special-q side) */

    polynomial<double> fijd;   /* coefficients of F_q (divided by q
                                * on the special q side) */
    polynomial<long double> fijld; /* same with long double */

    double scale;      /* scale used for logarithms for fb and norm.
                        * must be of form (int)x * 0.1 */

    private:
    void ctor_common(
            siever_config const & sc,
            cxx_cado_poly const & cpoly,
            int side);

    public:
    lognorm_base() = default;
    lognorm_base(lognorm_base const &) = default;
    lognorm_base& operator=(lognorm_base const &) = default;
    lognorm_base(
            siever_config const & sc,
            cxx_cado_poly const & cpoly,
            int side,
            special_q_data_class auto const & Q,
            int logI,
            uint32_t J);
    virtual ~lognorm_base() = default;

    void norm(
            cxx_mpz & x,
            int i,
            unsigned int j,
            special_q_data_class auto const & Q) const;

    unsigned char lognorm(
            int i,
            unsigned int j,
            special_q_data_class auto const & Q) const;

    virtual void fill(unsigned char * S, unsigned int N MAYBE_UNUSED) const {
        /* Whether we put something or not here is not really important.
         * A no-op would do as well. */
        memset(S, 255, 1U << LOG_BUCKET_REGION);
    }
};

/*}}}*/
struct lognorm_reference : public lognorm_base {/*{{{*/
    /* See init_degree_X_norms_bucket_region_referencecode for the
     * explanation of this table. */

    /* Number of bits used to estimate the norms with the old reference code.
     * Unused otherwise.
     * This should be large enough: it must be such that all norms are
     * smaller than 2^(2^NORM_BITS)
     * This imposes NORM_BITS >= 8, or even >= 9 for large factorizations. */
    static const int NORM_BITS = 10;

    std::array<unsigned char, 1 << NORM_BITS> lognorm_table {};

    lognorm_reference() = default;
    lognorm_reference(lognorm_reference const &) = default;
    lognorm_reference& operator=(lognorm_reference const &) = default;
    lognorm_reference(siever_config const & sc, cxx_cado_poly const & cpoly, int side, qlattice_basis const & Q, int logI, uint32_t J);
    ~lognorm_reference() override = default;
    void fill(unsigned char * S, unsigned int N) const override;
    private:
    void fill_alg(unsigned char * S, uint32_t N) const;
    void fill_rat(unsigned char * S, uint32_t N) const;
    void fill_siqs(unsigned char * S, uint32_t N) const;
};

/*}}}*/
struct lognorm_smart : public lognorm_base {/*{{{*/
    /* This table depends on the scale of the logarithm, so clearly it
     * can't be shared between sides.
     */
    std::array<double, 257> cexp2 {};
    private:
    static std::array<double, 257> cexp2_init(double scale);
    public:
    /* For degree>1 only: a piecewise linear approximation of the
     * polynomial, which is within an multiplicative factor of the
     * original one on the segment [-I,I]x{1}.
     */
    piecewise_linear_function G;
    lognorm_smart() = default;
    lognorm_smart(lognorm_smart const &) = default;
    lognorm_smart& operator=(lognorm_smart const &) = default;
    lognorm_smart(lognorm_smart &&) = default;
    lognorm_smart & operator=(lognorm_smart &&) = default;
    lognorm_smart(
            siever_config const & sc,
            cxx_cado_poly const & cpoly,
            int side,
            special_q_data_class auto const & Q,
            int logI,
            uint32_t J);
    ~lognorm_smart() override = default;
    void fill(unsigned char * S, unsigned int N) const override;
    private:
    void fill_rat_inner (unsigned char *S, int i0, int i1, unsigned int j0, unsigned int j1, polynomial<double> const & fijd) const;
    void fill_alg(unsigned char * S, uint32_t N) const;
    void fill_rat(unsigned char * S, uint32_t N) const;
    void fill_siqs(unsigned char * S, uint32_t N) const;
};

/*}}}*/
struct sieve_range_adjust {/*{{{*/
    qlattice_basis Q;
private:
    siever_config conf;         /* This "conf" field is only used for a
                                 * few fields:
                                 *      logA
                                 *      lpb
                                 * We're specifically *not* using the
                                 * sieving fields, since by design these
                                 * can be decided *after* the adjustment.
                                 */
    cxx_cado_poly const & cpoly;
    // int nb_threads;  // no longer needed.
    std::vector<polynomial<double>> fijd;
    int logA;
public:
    int logI = 0;
    uint32_t J = 0;

#if 0
    sieve_range_adjust(special_q const & doing, las_info const & las)
        : doing(doing), cpoly(las.cpoly), nb_threads(las.nb_threads)
    {
        /* See whether for this size of special-q, we have predefined
         * parameters (note: we're copying the default config, and then
         * we replace by an adjusted one if needed). */
        conf = las.config_pool.get_config_for_q(doing);
        /* These two will be adjusted in the process */
        logA = conf.logA;
        logI = J = 0;
    }
#endif

    /* This is only for desperate cases. In las-duplicates, for the
     * moment it seems that we're lacking the las_info structure...
     *
     * Note that the ctor for qlattice_basis calls SkewGauss
     */
    sieve_range_adjust(special_q const & doing, cxx_cado_poly const & cpoly, siever_config const & conf)
        : Q(doing, cpoly->skew)
        , conf(conf)
        , cpoly(cpoly)
        , logA(conf.logA)
    {
    }

    /* There are three strategies to do a post-SkewGauss adjustment of
     * the q-lattice basis.  */

    /* implementation is in las-norms.cpp */
    // all these functions return 0 if they feel that the special-q
    // should be discarded.
    int sieve_info_adjust_IJ();    // "raw" J.
    int sieve_info_update_norm_data_Jmax(bool keep_logI = false);
    int adjust_with_estimated_yield();

    // a fall-back measure for desperate cases.
    // XXX when estimated_yield() wins, this will probably no longer be
    // necessary.
    uint32_t get_minimum_J() const;
    void set_minimum_J_anyway();

    siever_config const& config() const { return conf; }
private:
    template<typename T> struct mat {
        std::array<T, 4> x;
        T const& operator()(int i, int j) const { return x[2*i+j]; }
        T & operator()(int i, int j) { return x[2*i+j]; }
        mat(T a, T b, T c, T d) : x { a, b, c, d } {}
        explicit mat(T y[4]) : x { y[0], y[1], y[2], y[3] } {}
        std::ostream& print_me(std::ostream& o) const {
            o << "["
                << x[0] << ", "
                << x[1] << ", "
                << x[2] << ", "
                << x[3] << "]";
            return o;
        }
    };
    template<typename T> struct vec {
        std::array<T, 2> x;
        vec(T a, T b) : x { a, b } {}
        explicit vec(T y[2]) : x { y[0], y[1] } {}
        T const& operator[](int i) const { return x[i]; }
        T & operator[](int i) { return x[i]; }
        T const& operator()(int i) const { return x[i]; }
        T & operator()(int i) { return x[i]; }
        std::ostream& print_me(std::ostream& o) const {
            o << "["
                << x[0] << ", "
                << x[1] << "]";
            return o;
        }
    };
    friend sieve_range_adjust::vec<double> operator*(sieve_range_adjust::vec<double> const& a, sieve_range_adjust::mat<int> const& m) ;
    friend qlattice_basis operator*(sieve_range_adjust::mat<int> const& m, qlattice_basis const& Q) ;
    void prepare_fijd();
    int round_to_full_bucket_regions(const char *, std::string const & s = std::string());
    double estimate_yield_in_sieve_area(mat<int> const& shuffle, int squeeze, unsigned int N);
};/*}}}*/

extern sieve_range_adjust::vec<double> operator*(sieve_range_adjust::vec<double> const& a, sieve_range_adjust::mat<int> const& m) ;
extern qlattice_basis operator*(sieve_range_adjust::mat<int> const& m, qlattice_basis const& Q) ;

#endif	/* CADO_LAS_NORMS_HPP */
