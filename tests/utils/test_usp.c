#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "macros.h"
#include "mpz_poly.h"
#include "tests_common.h"
#include "usp.h"

#define MAX_DEGREE 100

void
test_usp ()
{
    usp_root_interval R[MAX_DEGREE];

    for (int i = 0; i < MAX_DEGREE; i++)
        usp_root_interval_init(R[i]);

    mpz_poly p;
    mpz_poly_init(p, -1);

    {
        mpz_poly_set_from_expression(p, "x");
        int n = mpz_poly_number_of_real_roots_extra (p, 2, R);
        ASSERT_ALWAYS (n == 1);
        /* check [a/2^ka, b/2^kb] contains root 0 */
        ASSERT_ALWAYS (mpz_cmp_ui (R[0]->a, 0) <= 0);
        ASSERT_ALWAYS (0 <= mpz_cmp_ui (R[0]->b, 0));
    }

    {
        mpz_poly_set_from_expression(p, "2*x");
        int n = mpz_poly_number_of_real_roots_extra (p, 2, R);
        ASSERT_ALWAYS (n == 1);
        double root = usp_root_interval_refine (R[0], p, 1e-9);
        ASSERT_ALWAYS(-0.000001 <= root && root <= 0.000001);
    }

    struct {
        const char * expression;
        int expected_nroots;
    } examples[] = {
        { "x+1", 1},
        { "x-1", 1},
        { "x^2+1", 0},
        { "x^2+2", 0},
        { "x^2-1", 2},
        { "x^2-3*x+2", 2},
        { "x^2-2*x+2", 0},
        { "x^2+1000*x+2", 2},
        { "(x-1)*(x-2)*(x-3)", 3},
        { "(x-1)*(x-2)*(x-3)*(x-4)", 4},
        { "(x-1)*(x-2)*(x-3)*(x-4)*(x-5)", 5},
        { "(x-1)*(x-2)*(x-3)*(x-4)*(x-5)*(x-6)", 6},
        { "(x-1)*(x-2)*(2*x-3)", 3},
        { "(x-1)*(x-2)*(2*x-3)*(4*x-5)", 4},
        { "(x-1)*(2*x-3)*(x-2)*(x-3)*(x-4)", 5},
        { "2*x-3", 1},
        { "2*(x-17)*(5*x-17)", 2},
        { "2000*x^3 - 35000*x^2 + 2*x - 35", 1},
        { "33*x^4-124*x^3+137*x^2-32*x-11", 4},
        { "2*x^4 + 4*x^3 - 4*x^2 + 2*x + 3", 2},
        { "2*x^4 + x^3 - 3*x^2 + 2*x + 3", 2},
        { NULL, 0},
    };

    for(int i = 0 ; examples[i].expression ; i++) {
        mpz_poly_set_from_expression(p, examples[i].expression);
        int n = mpz_poly_number_of_real_roots(p);
        ASSERT_ALWAYS(n == examples[i].expected_nroots);
    }


    /* infinite loop that appeared on 32-bit computers, 2nd April 2014 */
    {
        mpz_poly_set_from_expression(p, "x^0*202156081496031128767910903808-x*1014774808369763543554925264896+x^2*305592973565950609559904059392+x^3*532741928198739321928042938368-x^4*265275048860303928171627020288+x^5*152821965764546794524860481536-x^6*40951570592501524318484692992");
        usp_root_interval_set_ui(R[0], 2, 0, 4, 0);
        double root = usp_root_interval_refine (R[0], p, 1e-9); /* root is near 3.00763029864372 */
        ASSERT_ALWAYS(3.00763 <= root && root <= 3.00764);
    }


    for (int d = 1; d < MAX_DEGREE; d++)
    {
        mpz_poly_set_randomb(p, d, state, 128,
                MPZ_POLY_SIGNED_COEFFICIENTS |
                MPZ_POLY_RRANDOM |
                MPZ_POLY_DEGREE_EXACT);
        int n = mpz_poly_number_of_real_roots(p);
        ASSERT_ALWAYS (0 <= n && n <= d);
    }

    {
        /* check with large coefficients */
        mpz_poly_set_randomb(p, 3, state, 2048,
                MPZ_POLY_UNSIGNED_COEFFICIENTS |
                MPZ_POLY_RRANDOM |
                MPZ_POLY_DEGREE_EXACT);
        int n = mpz_poly_number_of_real_roots(p);
        ASSERT_ALWAYS (0 <= n && n <= 3);
    }


    /* bug #18879 */
    {
        mpz_poly_set_from_expression(p,"x^0*-399370140104138021004049555602050048725505629855192896740786176+x*6152454519896677691538310427796089495118818146087045879381884928-x^2*44005628564777919654831813189277128400998342025364587400350662656+x^3*193693587647452249746878956308340478214944212045235974402799042560-x^4*586131524042923039573328295315283165577996285519208621852726394880+x^5*1289948018788260792668169964588439385725345410925233394477273448448-x^6*2129171214885956678206724080558522465836045133020822627924230275072+x^7*2677624344970803918386735215388265920949327908075617878580737867776-x^8*2578130043460498861049699336930493887253985498977331087052720046080+x^9*1891301125005133082116243679283866801297483539999500792756728496128+x^10*-1040585580532476627932685878679921764622233387744782090307763175424+x^11*416382268948646574906631753282253840199545851564089327142421135360+x^12*-114545848993024762527402482645186791330188904180948233736106803200+x^13*19391576494070251913587093980476025365869714617333619797109243904+x^14*-1524165754529617102335904328302294916528952457696967340985417728");
        int n = mpz_poly_number_of_real_roots_extra (p, 0, R);
        ASSERT_ALWAYS (n == 2);
        /* we only care about not hitting an infinite loop */
        usp_root_interval_refine (R[0], p, 6.103516e-05);
        usp_root_interval_refine (R[1], p, 6.103516e-05);
    }

    {
        mpz_poly_set_from_expression(p,"-1907753428620166438197776483820542501644580146356003274752+x*15452436600339983267007229667323593777349076211560511176704-x^2*54758274613009131316199852690823331935974939002324115259392+x^3*110882878534980127672983781446528375367066401782722337964032-x^4*140332817723613361637526395534759643288947238768788031668224+x^5*113666888827955590415811847475276525451763140364809078308864-x^6*57542498902827796543386171783637664156309699559757732380672+x^7*16645828445978446377629358188283303394024754747896961695744-x^8*2106687741183691798122020284370616107326915629763006562304");
        usp_root_interval_set_ui(R[0], 63, 6, 64, 6);
        usp_root_interval_refine (R[0], p, 0.000244);
    }

    for (int i = 0; i < MAX_DEGREE; i++)
        usp_root_interval_clear(R[i]);
    mpz_poly_clear(p);
}

int main(int argc, char const * argv[])
{
    tests_common_cmdline(&argc, &argv, PARSE_SEED | PARSE_ITER);
    test_usp ();
    tests_common_clear();
    exit (EXIT_SUCCESS);
}
