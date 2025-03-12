/**
 * \file cado-nfs/polyselect/polyselect_test.c
 * 
 * \date 21/08/2014
 * \author Aurore Guillevic
 * \email guillevic@lix.polytechnique.fr
 * \brief test gfpkdlpolyselect functions
 *
 * \test ask for polynomials for p of 60, 80, 100... dd for Fp2.
 * (TODO)
 */

#include "cado.h" // IWYU pragma: keep

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <gmp.h>

#include "gfpkdlpolyselect.h"
#include "gfpkdlpolyselect_impl.h"
#include "macros.h"
#include "portability.h"

extern const tPyf_t tab_t_Py_f_deg4_type0_h1[26];
extern const unsigned int table_f4_size;
extern const tPyf_poly_t table_f4;

// maybe: in a separate file.
// gfpkdlpolyselect_utils.c --> with all the auxiliary functions
// gfpkdlpolyselect.c       --> main function (interface with Python ?)
// gfpkdlpolyselect_test.c  --> with also a main test function



static void
usage (){
    fprintf (stderr, "./gfpkdlpolyselect -p xxx -k xxx [otional: -label xxx] \n");
    fprintf (stderr, "    with p a prime and k >= 2 the extension field degree\n");
    fprintf (stderr, "    and optionally, a label for the output filename.\n");
    exit (1);
}

/**
 * \brief test function
 * 
 * \param[in] p: size in decimal digits of p
 * \param[in] k: extension degree
 * \param<optional>[in] label: test label
 * \param[out] 
 *
 * \brief For a given prime p and extension degree k (=2 for the moment)
 *        compute the two polynomials (f,g) and save them in file p<k>dd<dd>[_<label>].poly
 *
 * inspired by main function of dlpolyselect.c
 */
int main(int argc, char const * argv[])
{
  // first, catch k and p
  // then, compute the two polynomials (f,g)
  int i, k = 0;
#define LABEL_SIZE_MAX 48
  char label[LABEL_SIZE_MAX];
  label[0] = '\0';
  mpz_t p;
  mpz_init (p);
  
    /* printf command-line. argv[0] = command name */
    printf ("#");
    for (i = 0; i < argc; i++)
        printf (" %s", argv[i]);
    printf ("\n");
    fflush (stdout);

    /* parsing */
    while (argc >= 3 && argv[1][0] == '-')
    {
        if (argc >= 3 && strcmp (argv[1], "-p") == 0) {
            mpz_set_str (p, argv[2], 10);
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-k") == 0) {
            char * p;
            k = (int) strtol(argv[2], &p, 0);
            ASSERT_ALWAYS(*p == '\0');
            argv += 2;
            argc -= 2;
        }
        else if (argc >= 3 && strcmp (argv[1], "-label") == 0) {
	  label[0] = '_';
	  int rc = strlcpy(label+1, argv[2], LABEL_SIZE_MAX-1);
          ASSERT_ALWAYS(rc < LABEL_SIZE_MAX-1);
          argv += 2;
          argc -= 2;
        }
        else {
            fprintf (stderr, "Invalid option: %s\n", argv[1]);
            usage();
            exit (1);
        }
    }

    if (mpz_cmp_ui (p, 0) <= 0) {
        fprintf (stderr, "Error, missing input prime number (-p option)\n");
        usage ();
    }
    if (k == 0) {
        fprintf (stderr, "Error, missing input or bad value extension degree (-k option)\n");
        usage ();
    }
    if (k < 2) {
        fprintf (stderr, "Error, error extension degree (-k option), should be >= 2\n");
        usage ();
    }

    gfpkdlpolyselect(k, p, NULL, 0, NULL);
    mpz_clear (p);
}
