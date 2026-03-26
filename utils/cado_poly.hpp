#ifndef CADO_POLY_HPP
#define CADO_POLY_HPP

#include <cstdio>
#include <cstdlib>

#include <vector>
#include <ostream>

#include <gmp.h>

#include "params.h"
#include "mpz_poly.h"
#include "cxx_mpz.hpp"

/* The maximum degree of polynomials supported. Used for statically 
   allocating storage (i.e. "mpz_t poly[MAX_DEGREE]") */
#define MAX_DEGREE 10

/* Same idea as for cxx_mpz and friends */
struct cxx_cado_poly : public std::vector<cxx_mpz_poly>
{
    using super = std::vector<cxx_mpz_poly>;
    int nsides() const { return (int) super::size(); }

    cxx_mpz n;

    /* skewness from poly file, if given, otherwise 0. It is important
     * that the default is 0, because cado_poly_set_skewness_if_undefined
     * works under this assumption.
     */
    double skew = 0;

    static void configure_switches(cxx_param_list &) {}
    static void configure_aliases(cxx_param_list & pl) {
        param_list_configure_alias(pl, "skew", "S");
    }
    static void declare_usage(cxx_param_list & pl) {
        param_list_decl_usage(pl, "poly", "polynomial file");
        param_list_decl_usage(pl, "skew", "skewness");
    }

    cxx_cado_poly() = default;

    struct plist {};

    cxx_cado_poly(plist const &, cxx_param_list & pl);

    cxx_cado_poly(cxx_param_list & pl) {
        const char *tmp;
        if ((tmp = param_list_lookup_string(pl, "poly")) == NULL) {
            fprintf(stderr, "Error: -poly is missing\n");
            param_list_print_usage(pl, NULL, stderr);
            exit(EXIT_FAILURE);
        }
        if (!read(tmp)) {
            ::fprintf(stderr, "Error reading polynomial file %s\n", tmp);
            exit(EXIT_FAILURE);
        }
        /* -skew (or -S) may override (or set) the skewness given in the
         * polynomial file */
        param_list_parse_double(pl, "skew", &(skew));
        if (skew <= 0.0) {
            fprintf(stderr, "Error, please provide a positive skewness\n");
            exit(EXIT_FAILURE);
        }
    }

    // ~cxx_cado_poly() = default;
    // cxx_cado_poly(cxx_cado_poly const & o) = default;
    // cxx_cado_poly & operator=(cxx_cado_poly const & o) = default;
    // cxx_cado_poly(cxx_cado_poly && o) = default;
    // cxx_cado_poly& operator=(cxx_cado_poly && o) = default;

#if 0
    cxx_mpz_poly const & operator[](int i) const {
        /* This is really a gross hack. Since an mpz_poly is really just
         * the same thing as an mpz_poly in memory, just return the
         * mpz_poly by pretending that it's a cxx_mpz_poly. I confess
         * it's ugly, and I'd of course like to make cxx_mpz_poly the
         * first class citizen instead, but I can't because of
         * polyselect.
         */
        ASSERT_ALWAYS(i < x.nsides());
        mpz_poly_srcptr p = x[i];
        // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
        return *reinterpret_cast<cxx_mpz_poly const *>(p);
    }
#endif

    // This reads a file created by polyselect and fill in the structure
    // accordingly. Return 1 if success, 0 if failure (and diagnostic on
    // stderr)
    int read(const char *filename);
    int read(FILE *);
    static int read_next_poly_from_stream (cxx_cado_poly &, FILE *);

    void fprintf (FILE *, const char * prefix = "") const;


    /* More functions for printing cado_poly are defined in polyselect/
     * as only binaries in polyselect/ use them and some functions (like
     * L2_skewness, ...) only defined in polyselect/ are needed.
     */

    void provision_new_poly() { super::emplace_back(); }
    void reset () {
        n = 0;
        skew = 1;
        super::clear();
    }

    private:
    int check_mapping_internal(cxx_mpz_poly * G, cxx_mpz const & N) const;
    public:
    int check_mapping(cxx_mpz_poly & G, cxx_mpz const & N) const
    {
        return check_mapping_internal(&G, N);
    }
    int check_mapping(cxx_mpz const & N) const
    {
        return check_mapping_internal(nullptr, N);
    }

    // Compute m as the common root of f and g mod N.
    // N is taken as a third argument; it can be a strict factor of the N
    // stored in the polynomial.
    // If this fails, then in most of the cases we have found a factor of N
    // that is given instead of m.
    // The return value tells whether it worked (and then m is the common
    // root) or it failed (and then m is a factor of N).
    int getm(cxx_mpz & m, cxx_mpz const & N) const;

    /* Return the rational side or -1 if two algebraic side */
    int get_ratside () const;

    std::string as_string(std::string const & prefix = {}) const;

    friend std::ostream& operator<<(std::ostream& os, cxx_cado_poly const & cpoly) {
        return os << cpoly.as_string();
    }
};

extern void cado_poly_fprintf_MurphyE (FILE *fp, const char * prefix, int side,
        double MurphyE, double bound_f, double bound_g, double area);

#endif	/* CADO_POLY_HPP */
