#include "cado.h" // IWYU pragma: keep

#include <cstring>
#include <cstdio>
#include <cstdlib>

#include <string>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "macros.h"  // for ASSERT_ALWAYS, ASSERT
#include "params.h"
#include "cado_poly.hpp"
#include "mpz_poly.h"
#include "portability.h"

// This function is no longer exported
#define BUF_MAX 10000

/* this is not the same as the cxx_cado_poly(cxx_param_list &) ctor.
 * Here, we expect the _contents_ of the poly file to be the parameters.
 */
cxx_cado_poly::cxx_cado_poly(cxx_cado_poly::plist const &, cxx_param_list & pl)
{
    /* Parse skew value. Set to 0.0 to ensure that we get an invalid skewness
       in case it is not given */
    skew = 0.0;
    param_list_parse_double(pl, "skew", &(skew));

    super::clear();

    /* Parse polynomials. Either in line format (poly%d=%Zd,%Zd,...) either given
     * by coefficients. */
    for( ;; ) {
        std::string tag = fmt::format("poly{}", nsides());
        if (!param_list_lookup_string(pl, tag.c_str()))
            break;
        /* sure, we could probably use the result of the string lookup.
         * Well, that's not how param_list_parse_mpz_poly works.
         * Note that param_list_parse_mpz_poly in itself supports two
         * distinct formats: comma-separated lists of coefficients, or
         * algebraic expressions.
         */
        super::emplace_back();
        param_list_parse(pl, tag, super::back());
    }

    if (super::empty())
        /* Convention: c or X means side 1 (often algebraic)
         *             Y means side 0 (often rational) */
    {
        super::emplace_back();
        super::emplace_back();
        mpz_t coeff;
        mpz_init(coeff);
        /* reading polynomials coefficient by coefficient */
        for (int i = 0; i <= MAX_DEGREE; i++)
        {
            char tc[4]; snprintf(tc, sizeof(tc), "c%d", i);
            char tX[4]; snprintf(tX, sizeof(tX), "X%d", i);
            char tY[4]; snprintf(tY, sizeof(tY), "Y%d", i);
            if (param_list_parse_mpz(pl, tc, coeff))
                mpz_poly_setcoeff((*this)[1], i, coeff);
            if (param_list_parse_mpz(pl, tX, coeff))
                mpz_poly_setcoeff((*this)[1], i, coeff);
            if (param_list_parse_mpz(pl, tY, coeff))
                mpz_poly_setcoeff((*this)[0], i, coeff);
        }
        mpz_clear(coeff);
    }

    /* Parse value of N. Two keys possible: n or None. Return 0 if not found. */
    if (!param_list_parse_mpz(pl, "n", n) &&
            !param_list_parse_mpz(pl, "N", n) &&
            !param_list_parse_mpz(pl, NULL, n))
        throw std::runtime_error("Error, no value for N in cado_poly_set_plist\n");

    for (int side = 0; side < nsides(); side++)
        if ((*this)[side]->deg < 0)
            throw std::runtime_error(fmt::format("Error, polynomial on side {} has degree < 0 in "
                        "cado_poly_set_plist\n", side));

    for (int side = 0; side < nsides(); side++) {
        auto & g = (*this)[side];
        cxx_mpz c;
        mpz_poly_content(c, g);
        if (mpz_cmp_ui(c, 1) != 0) {
            fmt::print (stderr, "# Warning: the polynomial on side {}"
                    " has non-trivial content (gcd of coeffs = {}).\n"
                    "# This is a very unexpected situation."
                    " For all useful cases, it seems that\n"
                    "# the correct thing to do is to divide the polynomial"
                    " by its content. Doing it now.\n",
                    side, c);
            fmt::print(stderr, "# Old poly: {}\n", g);
            mpz_poly_divide_by_content(g);
            fmt::print(stderr, "# Reduced poly: {}\n", g);
        }
    }

    /* check that N divides the resultant*/
    if (!check_mapping(n))
        throw std::runtime_error("Error, polynomial file is inconsisent (no common mapping to target ring)\n");
}

// returns 0 on failure, 1 on success.
int cxx_cado_poly::read (FILE * f)
{
  cxx_param_list pl;
  param_list_read_stream (pl, f, 0);
  *this = cxx_cado_poly(plist {}, pl);
  return 1;
}

// returns 0 on failure, 1 on success.
int cxx_cado_poly::read_next_poly_from_stream (cxx_cado_poly & res, FILE * f)
{
    int r;
    cxx_param_list pl;
    r = param_list_read_stream (pl, f, 1);
    if (r && !param_list_empty(pl)) {
        res = cxx_cado_poly(plist {}, pl);
        r = 1;
    } else {
        r = 0;
    }
    return r;
}

// returns 0 on failure, 1 on success.
int cxx_cado_poly::read(const char *filename)
{
    FILE *file;
    int r;
    const char * magic = "inline-poly://";
    if (strncmp(filename, magic, strlen(magic)) == 0) {
        /* interpret the reset as param list items, and create the
         * polynomial like this.
         */
        cxx_param_list pl;
        const char * p = filename + strlen(magic);
        for( ; *p ; ) {
            const char * q = strchr(p, '/');
            char * newitem;
            if (q == NULL) {
                newitem = strdup(p);
                p += strlen(p);
            } else {
                newitem = strndup(p, q - p);
                p = q + 1;
            }
            char * newkey;
            char * newvalue;
            q = strchr(newitem, '=');
            if (q == NULL) {
                fprintf(stderr, "wrong value in inline-poly file\n");
                exit (EXIT_FAILURE);
            }
            newitem[q-newitem] = '\0';
            newkey = newitem;
            newvalue = newitem + (q-newitem+1);
            param_list_add_key(pl, newkey, newvalue, PARAMETER_FROM_FILE);
            free(newitem);
        }
        *this = cxx_cado_poly(plist {}, pl);
        r = 1;
    } else {
        file = fopen(filename, "r");
        if (file == NULL)
        {
          ::fprintf(stderr, "Error, in cado_poly_read: could not open %s\n", filename);
          return 0;
        }
        r = read(file);
        fclose(file);
    }
    return r;
}

/* check that there is a mapping from Z[x] to the target ring (Z/N,
 * GF(p), GF(p^k)) that factors through each of the polynomials.
 *
 * If we have only two sides, rational or not, this is clearly something
 * that we want, for NFS to make any sort of sense.
 *
 * The question of how to generalize this check to multiple polynomials
 * is not entirely clear.
 *
 * If there is a rational side (only one, see cado_poly_get_ratside),
 * having a resultant with the rational polynomial that is a multiple of
 * N is actually a transitive relation, so that it suffices to check the
 * condition with the rational side only.
 *
 * If we have multiple nonlinear polynomials, we have to make sure that
 * there is a compatible mapping. This is done by checking the pseudo-gcd
 * of the set of polynomials in cado_poly
 *
 * => The conclusion of all that is that what we really care about is the
 * computation of the gcd, and this gcd is returned in the polynomial G.
 * Note that if G is null, then we don't store it, but we do the check
 * nevertheless.
 *
 * This returns true (1) on success, and false (0) on error, after an
 * error message is printed to stderr.
 */
int cxx_cado_poly::check_mapping_internal(cxx_mpz_poly * G, cxx_mpz const & N) const
{
    if (nsides() == 1)
    {
        if (G)
            *G = (*this)[0];
        return 1;
    }

    /* scratch space for mpz_poly_pseudogcd_mpz, which destroys its
     * input.
     */
    cxx_mpz_poly S[2];

    /* a polynomial for the gcd */
    cxx_mpz_poly G_loc;

    /* factors of N, if any. This should never happen */
    mpz_t factors_of_N, ftmp;
    mpz_init_set_ui(factors_of_N, 1);
    mpz_init_set_ui(ftmp, 0);

    /* A copy of N with all trivial factors removed */
    mpz_t Nsmall;
    mpz_init(Nsmall);

    for(int ret = 0 ; ret == 0; ) {
        mpz_divexact(Nsmall, N, factors_of_N);

        /* compute all pseudo-gcds */
        for(int i = 1 ; i < nsides() ; i++) {
            mpz_poly_set(S[0], (*this)[0]);
            mpz_poly_set(S[1], (*this)[i]);
            ret = mpz_poly_pseudogcd_mpz (S[0], S[1], Nsmall, ftmp);
            if (ret == 0) {
                mpz_mul(factors_of_N, factors_of_N, ftmp);
                break;
            }
            ret = mpz_poly_pseudogcd_mpz (G_loc, S[0], Nsmall, ftmp);
            if (ret == 0) {
                mpz_mul(factors_of_N, factors_of_N, ftmp);
                break;
            }
        }

        /* If a factor was found, start over, otherwise we're done with
         * this silly "restart if small factors pop up" loop (which will
         * presumably never be activated anyway...)
         */
    }

    int found_mapping = G_loc->deg >= 1;

    if (G)
        *G = G_loc;

    if (mpz_cmp_ui(factors_of_N, 1) != 0) {
        /* I don't think that there's any reason to have a situation
         * where N has non-trivial factors. Well, it happens if you feed
         * random numbers to polyselect, though.
         */
        gmp_fprintf(stderr, "Warning: non-trivial factors (%Zd) of N were found while checking the cpoly file. This should not happen.\n", factors_of_N);
        // found_mapping = 0;
    }


    mpz_clear(Nsmall);
    mpz_clear(ftmp);
    mpz_clear(factors_of_N);

    if (!found_mapping) {
        fprintf (stderr, "Error, we found no compatible mapping that factors through the given number fields. Your poly file is most probably wrong.\n");
    }

    return found_mapping;
}


/* If NULL is passed for m, then, just check that N divides the resultant.
 * Assume there at least two polynomials.
 */
int cxx_cado_poly::getm(cxx_mpz & m, cxx_mpz const & N) const
{
    cxx_mpz_poly G;

    ASSERT_ALWAYS(check_mapping(G, N));

    /* It doesn't make sense to ask for m when the mapping goes to a
     * larger structure. Note that this is not necessarily equivalent to
     * having no rational side: Two algebraic sides may well have only a
     * rational mapping in common (case of Joux-Lercier polynomial
     * selection).
     */

    if (mpz_poly_degree(G) != 1) {
        fprintf(stderr, "Error, cannot determine m given that the common mapping is not linear\n");
        exit (EXIT_FAILURE);
    }

    cxx_mpz inv;
    int ret2 = mpz_invert(inv, mpz_poly_coeff_const(G, 1), N);
    // This inversion should always work.
    // If not, it means that N has a small factor (not sure we want
    // to be robust against that...). This should have been reported,
    // at least as a warning, by cado_poly_check_mapping.
    // Or maybe the polynomial selection was really bogus!
    ASSERT_ALWAYS(ret2);
    mpz_mul(inv, inv, mpz_poly_coeff_const(G, 0));
    mpz_neg(inv, inv);
    mpz_mod(m, inv, N);

    return 1;
}

/* Return the rational side or -1 if only algebraic sides */
/* Assume that there is at most 1 rational side (renumber.cpp does the same
 * assumption).
 *
 * Note that several distinct rational sides would not really make sense.
 * A rational side prescribes the mapping from the indeterminate (the X
 * in a-b*X) to the underlying ring (Z/N, GF(p), ...). So if we have
 * distinct rational sides, we would have distinct mappings, and I don't
 * think that we could do something useful in such a situation.
 */
int
cxx_cado_poly::get_ratside () const
{
  int number_of_rational_sides = 0;
  int ratside = -1;
  for(int side = 0; side < nsides() ; side++)
    if((*this)[side]->deg == 1) {
        ratside = side;
        number_of_rational_sides++;
    }
  ASSERT_ALWAYS(number_of_rational_sides <= 1);
  return ratside;
}

void
cxx_cado_poly::fprintf (FILE *fp, const char * prefix) const
{
    auto s = as_string(prefix);
    fputs(s.c_str(), fp);
}

std::string
cxx_cado_poly::as_string (std::string const & prefix) const
{
    std::string s = prefix;
    s += fmt::format ("n: {}\n", n);

  if (nsides() == 2) {
      char * p0, *p1;
      mpz_poly_asprintf_cado_format (&p0, (*this)[0], 'Y', prefix.c_str());
      mpz_poly_asprintf_cado_format (&p1, (*this)[1], 'c', prefix.c_str());
      s += p0;
      s += p1;
      free(p0);
      free(p1);
  } else {
    for (int side = 0; side < nsides(); side++) {
      s += fmt::format("{}poly{}={}\n", prefix, side, (*this)[side]);
      // ::fprintf (fp, "poly%u=", side);
      // mpz_poly_fprintf_coeffs (fp, (*this)[side], ",");
    }
  }

  s += fmt::format("{}skew: {:1.3f}\n", prefix, skew);

  return s;
}

void
cado_poly_fprintf_MurphyE (FILE *fp, const char * prefix, int side,
        double MurphyE, double bound_f, double bound_g, double area)
{
  if (prefix)
    fputs (prefix, fp);
  /* Always print "# " after the prefix and before the MurphyE line. */
  fprintf (fp, "# side %d MurphyE(Bf=%.3e,Bg=%.3e,area=%.3e)=%.3e\n", side, bound_f, bound_g,
               area, MurphyE);
}
