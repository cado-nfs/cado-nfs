#ifndef LAS_NORMS_HPP_
#define LAS_NORMS_HPP_

#include <stdint.h>
#include "las-siever-config.hpp"
#include "las-qlattice.hpp"
#include "mpz_poly.h"
#include "double_poly.h"
#include "cado_poly.h"
#include "logapprox.hpp"

double get_maxnorm_rectangular (double_poly_srcptr src_poly, const double X, const double Y);

struct lognorm_base {/*{{{*/
    int logI, J;

    unsigned char bound; /* A sieve array entry is a sieve survivor if it is
                            at most "bound" on each side */
    protected:
    double maxlog2;      /* Bound on the log in base 2. This is
                            intermediary data, really. */
    public:
    cxx_mpz_poly fij;  /* coefficients of F(a0*i+a1*j, b0*i+b1*j)
                        * (divided by q on the special-q side) */

    cxx_double_poly fijd;      /* coefficients of F_q (divided by q
                                * on the special q side) */

    double scale;      /* scale used for logarithms for fb and norm.
                        * must be of form (int)x * 0.1 */

    lognorm_base(siever_config const & sc, cxx_cado_poly const & cpoly, int side, qlattice_basis const & Q, int logI, int J);

    void norm(mpz_ptr x, int i, unsigned int j) const;
    unsigned char lognorm(int i, unsigned int j) const;

    virtual void fill(unsigned char * S, int N MAYBE_UNUSED) const {
        /* Whether we put something or not here is not really important.
         * A no-op would do as well. */
        memset(S, 255, 1U << LOG_BUCKET_REGION);
    }
};

/*}}}*/
struct lognorm_reference : public lognorm_base {/*{{{*/
    /* See init_degree_X_norms_bucket_region_referencecode for the
     * explanation of this table. */
    unsigned char lognorm_table[1 << NORM_BITS];

    lognorm_reference(siever_config const & sc, cxx_cado_poly const & cpoly, int side, qlattice_basis const & Q, int logI, int J);
    virtual void fill(unsigned char * S, int N) const;
};

/*}}}*/
struct lognorm_smart : public lognorm_base {/*{{{*/
    /* This table depends on the scale of the logarithm, so clearly it
     * can't be shared between sides.
     */
    double cexp2[257];
    /* For degree>1 only: a piecewise linear approximation of the
     * polynomial, which is within an multiplicative factor of the
     * original one on the segment [-I,I]x{1}.
     */
    piecewise_linear_function G;
    lognorm_smart(siever_config const & sc, cxx_cado_poly const & cpoly, int side, qlattice_basis const & Q, int logI, int J);
    virtual void fill(unsigned char * S, int N) const;
};

/*}}}*/
struct sieve_range_adjust {/*{{{*/
    friend struct sieve_info;
private:
    las_todo_entry doing;
    siever_config conf;         /* This "conf" field is only used for a
                                 * few fields:
                                 *      logA
                                 *      sublat
                                 *      lpb
                                 * We're specifically *not* using the
                                 * sieving fields, since by design these
                                 * can be decided *after* the adjustment.
                                 */
    cado_poly_srcptr cpoly;
    int nb_threads;
    cxx_double_poly fijd[2];
    int logA;
public:
    int logI;
    uint32_t J;
    qlattice_basis Q;

#if 0
    sieve_range_adjust(las_todo_entry const & doing, las_info const & las)
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
     * moment it seems that we're lacking the las_info structure... */
    sieve_range_adjust(las_todo_entry const & doing, cado_poly_srcptr cpoly, siever_config const & conf, int nb_threads = 1)
        : doing(doing), conf(conf), cpoly(cpoly), nb_threads(nb_threads)
    {
        logA = conf.logA;
        logI = J = 0;
    }
    sieve_range_adjust() = default;
    sieve_range_adjust(sieve_range_adjust&&) = default;
    sieve_range_adjust& operator=(sieve_range_adjust&&) = default;


    int SkewGauss() {
        int ret = ::SkewGauss(Q, doing.p, doing.r, cpoly->skew);
        Q.set_q(doing.p, doing.prime_sq);
        if (!Q.prime_sq)
            Q.prime_factors = Q.prime_factors;
        return ret;
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
    uint32_t get_minimum_J();
    void set_minimum_J_anyway();

    siever_config const& config() const { return conf; }
private:
    template<typename T> struct mat {
        T x[4];
        T const& operator()(int i, int j) const { return x[2*i+j]; }
        T & operator()(int i, int j) { return x[2*i+j]; }
        mat(T a, T b, T c, T d) { x[0]=a; x[1]=b; x[2]=c; x[3]=d; }
        mat(T y[4]) { x[0]=y[0]; x[1]=y[1]; x[2]=y[2]; x[3]=y[3]; }
    };
    template<typename T> struct vec {
        T x[2];
        vec(T a, T b) { x[0] = a; x[1] = b; }
        vec(T y[2]) { x[0] = y[0]; x[1] = y[1]; }
        T const& operator[](int i) const { return x[i]; }
        T & operator[](int i) { return x[i]; }
        T const& operator()(int i) const { return x[i]; }
        T & operator()(int i) { return x[i]; }
    };
    friend sieve_range_adjust::vec<double> operator*(sieve_range_adjust::vec<double> const& a, sieve_range_adjust::mat<int> const& m) ;
    friend qlattice_basis operator*(sieve_range_adjust::mat<int> const& m, qlattice_basis const& Q) ;
    void prepare_fijd();
    int adapt_threads(const char *);
    double estimate_yield_in_sieve_area(mat<int> const& shuffle, int squeeze, int N);
};/*}}}*/

#endif	/* LAS_NORMS_HPP_ */
