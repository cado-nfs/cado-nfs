#include "cado.h" // IWYU pragma: keep
#include <stdint.h>    // for uint64_t, int64_t
#include <stdlib.h>    // for EXIT_FAILURE, EXIT_SUCCESS
#include <inttypes.h>
#include <stdio.h>
#include <gmp.h>
#include "sm_utils.hpp" // sm_relset_t
#include "mpz_poly.h"   // mpz_poly_srcptr
#include "macros.h"

#define FREQ 2 // when possible one time out of FREQ we try sm_single_rel

/* Return number of errors */
int
test_sm (FILE * datafile)
{
  int err = 0;
  do
  {
    int ret, degF, degN, degD, nb_ab, nbSM;
    unsigned int nb_test_single_rel = 0;
    cxx_mpz_poly F, N, Nc, D, Dc, SM, SMc;
    cxx_mpz tmp, ell;
    int64_t a, e[MAX_LEN_RELSET];
    uint64_t b, len_relset, r[MAX_LEN_RELSET];
    std::vector<pair_and_sides> ab_polys;
    sm_relset_t relset;

    ret = fscanf(datafile, "in %d", &degF);
    if (ret == EOF)
      break;
    else
      ASSERT_ALWAYS (ret == 1);

    for (int i = 0; i <= degF; i++)
    {
      gmp_fscanf (datafile, " %Zi", (mpz_ptr) tmp);
      mpz_poly_setcoeff(F, i, tmp);
    }

    gmp_fscanf (datafile, " %Zi", (mpz_ptr) ell);

    sm_side_info sm_info(F, ell, 0);

    ret = fscanf(datafile, " %d", &nb_ab);
    ASSERT_ALWAYS (ret == 1);

    for (int i = 0; i < nb_ab; i++)
    {
      ret = fscanf (datafile, " %" SCNd64 " %" SCNu64 "", &a, &b);
      ASSERT_ALWAYS (ret == 2);
      pair_and_sides ps(a, b, 0, 1);
      ab_polys.push_back(ps);
    }

    ret = fscanf(datafile, " %" SCNu64 "", &len_relset);
    ASSERT_ALWAYS (ret == 1);
    ASSERT_ALWAYS (len_relset <= MAX_LEN_RELSET);

    for (uint64_t i = 0; i < len_relset; i++)
    {
      ret = fscanf (datafile, " %" SCNu64 " %" SCNd64 "", &r[i], &e[i]);
      ASSERT_ALWAYS (ret == 2);
    }

    ret = fscanf(datafile, "\nout %d", &degN);
    ASSERT_ALWAYS (ret == 1);
    ASSERT_ALWAYS (degN >= -1);
#ifdef __COVERITY__
    __coverity_mark_pointee_as_sanitized(&degN, GENERIC);
#endif

    for (int i = 0; i <= degN; i++)
    {
      gmp_fscanf (datafile, " %Zi",(mpz_ptr)  tmp);
      mpz_poly_setcoeff(N, i, tmp);
    }
    
    ret = fscanf(datafile, " %d", &degD);
    ASSERT_ALWAYS (ret == 1);

    for (int i = 0; i <= degD; i++)
    {
      gmp_fscanf (datafile, " %Zi", (mpz_ptr) tmp);
      mpz_poly_setcoeff(D, i, tmp);
    }
    
    ret = fscanf(datafile, " %d", &nbSM);
    ASSERT_ALWAYS (ret == 1);
    ASSERT_ALWAYS (nbSM <= degF);

    for (int i = 0; i < nbSM; i++)
    {
      gmp_fscanf (datafile, " %Zi", (mpz_ptr) tmp);
      mpz_poly_setcoeff(SM, i, tmp);
    }

    char c = ' ';
    ret = fscanf(datafile, "%c", &c);
    ASSERT_ALWAYS (ret == 1 && c == '\n');
    

    //     /* Real tests begin here */
    //     /* artificially duplicate data, to test both sides */
    //     mpz_poly_ptr FF[2];
    //     FF[0] = &F[0]; FF[1] = &F[0];
    //     cxx_mpz_poly SMc2;
    //     mpz_poly_ptr SSMc[2];
    //     SSMc[0] = &SMc[0]; SSMc[1] = &SMc2[0];
    if (len_relset == 1 && e[0] == 1 && nb_test_single_rel % FREQ == 0)
    {
      nb_test_single_rel++;
      sm_info.compute_piecewise(SMc, ab_polys[r[0]].ab);
    } else {
      mpz_poly_srcptr FF[2] = {F, F};
      sm_relset_init (relset, FF, 2);
      sm_build_one_relset (relset, r, e, len_relset, ab_polys, FF, 2, sm_info.ell2);
      mpz_poly_set (Nc, relset->num[0]);
      mpz_poly_set (Dc, relset->denom[0]);
      mpz_poly_reduce_frac_mod_f_mod_mpz (relset->num[0], relset->denom[0],
              F, sm_info.ell2);
      sm_info.compute_piecewise (SMc, relset->num[0]);
      sm_relset_clear (relset);
    }
    // mpz_poly_clear(SMc2);

    /* In case of error, print all relevant information */
    if (mpz_poly_cmp(SM, SMc) != 0)
    {
      err++;
      fprintf (stderr, "### ERROR: computation of SM is wrong with:\nF = ");
      mpz_poly_fprintf(stderr, F);
      gmp_fprintf(stderr, "ell = %Zi\nell2 = %Zi\n\n",
              (mpz_srcptr) ell, (mpz_srcptr) sm_info.ell2);
      sm_info.print(stderr);
      fprintf (stderr, "# Relation-set is:\n%" PRIu64 "", len_relset);
      for (uint64_t i = 0; i < len_relset; i++)
        fprintf (stderr, " %" PRIu64 ":%" PRId64 "", r[i], e[i]);
      fprintf (stderr, "\n# (a,b) pairs are:\n");
      for (int i = 0; i < nb_ab; i++)
      {
        cxx_mpz tmp;
        mpz_neg(tmp, mpz_poly_coeff_const(ab_polys[i].ab, 1));
        gmp_fprintf (stderr, "%d %" SCNd64 ",%" SCNu64 "\n", i,
                mpz_poly_coeff_const(ab_polys[i].ab, 0),
                (mpz_srcptr) tmp);
      }
      if (mpz_poly_cmp (N, Nc) != 0) {
        fprintf (stderr, "# Expected numerator in fraction corresponding to "
                         "the relation-set:\n");
        mpz_poly_fprintf (stderr, N);
        fprintf (stderr, "# Instead computed numerator is:\n");
        mpz_poly_fprintf (stderr, Nc);
      }
      if (mpz_poly_cmp (D, Dc) != 0)
      {
        fprintf (stderr, "# Expected denominator in fraction corresponding to "
                         "the relation-set:\n");
        mpz_poly_fprintf (stderr, D);
        fprintf (stderr, "# Instead computed denominator is:\n");
        mpz_poly_fprintf (stderr, Dc);
      }
      fprintf (stderr, "# Values of SM should be:\n");
      for (int i = 0; i < nbSM; i++) {
        gmp_fprintf (stderr, "%Zi ", mpz_poly_coeff_const(SM, i));
      }
      fprintf (stderr, "\n# but computed values of SM are:\n");
      for (int i = 0; i < nbSM; i++) {
        gmp_fprintf (stderr, "%Zi ", mpz_poly_coeff_const(SMc, i));
      }
      fprintf (stderr, "\n#######################\n");
    }
  } while (1);

  return err;
}

int main(int argc, char const * argv[])
{
  FILE *datafile = NULL;

  ASSERT_ALWAYS (argc == 2);
  const char * datafilename = argv[1];
  datafile = fopen (datafilename, "r");

  int const err = test_sm (datafile);

  fclose (datafile); 
  if (err)
    fprintf (stderr, "# %d erro%s found\n", err, (err == 1) ? "r" : "rs");
  return (err) ? EXIT_FAILURE : EXIT_SUCCESS;
}
