#include "cado.h"

#include <stdio.h>

#include "portability.h"
#include "utils.h"

#define TEST_MAX_AB 16
#define FREQ 2 // when possible one time out of FREQ we try sm_single_rel

void
mpz_poly_getcoeff_wrapper (mpz_t res, int i, const mpz_poly_t f)
{
  if (i <= f->deg)
    mpz_poly_getcoeff (res, i, f);
  else
    mpz_set_ui (res, 0);
}


/* Return number of errors */
int
test_sm (FILE * datafile)
{
  int err = 0;
  do
  {
    int ret, degF, degN, degD, nb_ab, nbSM;
    unsigned int nb_test_single_rel = 0;
    mpz_poly_t F, N, Nc, D, Dc, SM, SMc;
    mpz_t tmp, ell;
    int64_t a, e[MAX_LEN_RELSET];
    uint64_t b, len_relset, r[MAX_LEN_RELSET];
    mpz_poly_t ab_polys[TEST_MAX_AB];
    sm_relset_t relset;

    ret = fscanf(datafile, "in %d", &degF);
    if (ret == EOF)
      break;
    else
      ASSERT_ALWAYS (ret == 1);

    mpz_poly_init (F, degF);
    mpz_init (tmp);
    mpz_init (ell);

    for (int i = 0; i <= degF; i++)
    {
      gmp_fscanf (datafile, " %Zi", tmp);
      mpz_poly_setcoeff(F, i, tmp);
    }

    gmp_fscanf (datafile, " %Zi", ell);

    sm_side_info sm_info;
    sm_side_info_init(sm_info, F, ell);

    ret = fscanf(datafile, " %d", &nb_ab);
    ASSERT_ALWAYS (ret == 1);

    for (int i = 0; i < nb_ab; i++)
    {
      ret = fscanf (datafile, " %" SCNd64 " %" SCNu64 "", &a, &b);
      ASSERT_ALWAYS (ret == 2);
      mpz_poly_init_set_ab (ab_polys[i], a, b);
    }

    ret = fscanf(datafile, " %" SCNu64 "", &len_relset);
    ASSERT_ALWAYS (ret == 1);

    for (uint64_t i = 0; i < len_relset; i++)
    {
      ret = fscanf (datafile, " %" SCNu64 " %" SCNd64 "", &r[i], &e[i]);
      ASSERT_ALWAYS (ret == 2);
    }

    ret = fscanf(datafile, "\nout %d", &degN);
    ASSERT_ALWAYS (ret == 1);
    mpz_poly_init (N, degN);

    for (int i = 0; i <= degN; i++)
    {
      gmp_fscanf (datafile, " %Zi", tmp);
      mpz_poly_setcoeff(N, i, tmp);
    }
    
    ret = fscanf(datafile, " %d", &degD);
    ASSERT_ALWAYS (ret == 1);
    mpz_poly_init (D, degD);

    for (int i = 0; i <= degD; i++)
    {
      gmp_fscanf (datafile, " %Zi", tmp);
      mpz_poly_setcoeff(D, i, tmp);
    }
    
    ret = fscanf(datafile, " %d", &nbSM);
    ASSERT_ALWAYS (ret == 1);
    mpz_poly_init (SM, nbSM-1);

    for (int i = 0; i < nbSM; i++)
    {
      gmp_fscanf (datafile, " %Zi", tmp);
      mpz_poly_setcoeff(SM, i, tmp);
    }

    char c = ' ';
    ret = fscanf(datafile, "%c", &c);
    ASSERT_ALWAYS (ret == 1 && c == '\n');
    

    mpz_poly_init(SMc, F->deg);

    //     /* Real tests begin here */
    //     barrett_init(invl2, ell2);
    //     mpz_poly_init (SMc, degF);
    mpz_poly_init (Nc, degF);
    mpz_poly_init (Dc, degF);
    //     /* artificially duplicate data, to test both sides */
    //     mpz_poly_ptr FF[2];
    //     FF[0] = &F[0]; FF[1] = &F[0];
    //     mpz_poly_t SMc2;
    //     mpz_poly_init(SMc2, degF);
    //     mpz_poly_ptr SSMc[2];
    //     SSMc[0] = &SMc[0]; SSMc[1] = &SMc2[0];
    if (len_relset == 1 && e[0] == 1 && nb_test_single_rel % FREQ == 0)
    {
      nb_test_single_rel++;
      if (nb_test_single_rel % (FREQ + 1))
          compute_sm_piecewise(SMc, ab_polys[r[0]], sm_info);
      else
          compute_sm_straightforward(SMc, ab_polys[r[0]], sm_info);

    } else {
        mpz_poly_ptr FF[2] = {F, F};
        int dd[2] = {degF, degF};
      sm_relset_init (relset, dd, 2);
      sm_build_one_relset (relset, r, e, len_relset, ab_polys, FF, 2, sm_info->ell2);
      mpz_poly_set (Nc, relset->num[0]);
      mpz_poly_set (Dc, relset->denom[0]);
      mpz_poly_reduce_frac_mod_f_mod_mpz (relset->num[0], relset->denom[0],
              F, sm_info->ell2, sm_info->invl2);
      compute_sm_straightforward (SMc, relset->num[0], sm_info);
      sm_relset_clear (relset, 2);
    }
    // mpz_poly_clear(SMc2);

    /* In case of error, print all relevant information */
    if (mpz_poly_cmp(SM, SMc) != 0)
    {
      err++;
      fprintf (stderr, "### ERROR: computation of SM is wrong with:\nF = ");
      mpz_poly_fprintf(stderr, F);
      gmp_fprintf(stderr, "ell = %Zi\nell2 = %Zi\n\n", ell, sm_info->ell2);
      sm_side_info_print(stderr, sm_info);
      fprintf (stderr, "# Relation-set is:\n%" PRIu64 "", len_relset);
      for (uint64_t i = 0; i < len_relset; i++)
        fprintf (stderr, " %" PRIu64 ":%" PRId64 "", r[i], e[i]);
      fprintf (stderr, "\n# (a,b) pairs are:\n");
      for (int i = 0; i < nb_ab; i++)
      {
        mpz_poly_getcoeff_wrapper (tmp, 0, ab_polys[i]);
        a = mpz_get_si (tmp);
        mpz_poly_getcoeff_wrapper (tmp, 1, ab_polys[i]);
        b = mpz_get_ui (tmp);
        fprintf (stderr, "%d %" SCNd64 ",%" SCNu64 "\n", i, a, b);
      }
      if (mpz_poly_cmp (N, Nc) != 0)
      {
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
      for (int i = 0; i < nbSM; i++)
      {
        mpz_poly_getcoeff_wrapper (tmp, i, SM);
        gmp_fprintf (stderr, "%Zi ", tmp);
      }
      fprintf (stderr, "\n# but computed values of SM are:\n");
      for (int i = 0; i < nbSM; i++)
      {
        mpz_poly_getcoeff_wrapper (tmp, i, SMc);
        gmp_fprintf (stderr, "%Zi ", tmp);
      }
      fprintf (stderr, "\n#######################\n");
    }


    for (int i = 0; i < nb_ab; i++)
      mpz_poly_clear (ab_polys[i]);
    mpz_clear (tmp);
    mpz_clear (ell);
    mpz_poly_clear (F);
    mpz_poly_clear (N);
    mpz_poly_clear (Nc);
    mpz_poly_clear (D);
    mpz_poly_clear (Dc);
    mpz_poly_clear (SMc);
    mpz_poly_clear (SM);
  } while (1);

  return err;
}

int
main (int argc, char **argv)
{
  FILE *datafile = NULL;

  ASSERT_ALWAYS (argc == 2);
  const char * datafilename = argv[1];
  datafile = fopen (datafilename, "r");

  int err = test_sm (datafile);

  fclose (datafile); 
  if (err)
    fprintf (stderr, "# %d erro%s found\n", err, (err == 1) ? "r" : "rs");
  return (err) ? EXIT_FAILURE : EXIT_SUCCESS;
}
