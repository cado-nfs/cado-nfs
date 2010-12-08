#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h>
#include "auxiliary.h"
#include "murphyE.h"
#include "rootsieve.h"
#include "rho.h"


//print a floating point GMP number using scientific notation
//strangely this function doesn't seem to be implemented in GMP
void printscientific(mpf_t x) {

    mpf_t mantisse;
    mpf_init2(mantisse, 10);
    mpf_set(mantisse, x);
    int exposant = 0;
    while(mpf_cmp_ui(mantisse, 10) > 0) {
        mpf_div_ui(mantisse, mantisse, 10);
        exposant++;
    }
    gmp_printf ("%Ff", mantisse);
    printf("e%d", exposant);

    mpf_clear(mantisse);

}


//returns the square root of N modulo p
//p is here a long integer. we may have to adapt this procedure to the case p is too big to be an unsigned int
//We use Tonelli-Shanks' algorithm
unsigned long int modular_square_root(mpz_t N, unsigned long int p) {

    mpz_t P;
    mpz_init_set_ui(P, p);

    assert(mpz_kronecker(N, P) == 1);

    mpz_t n;
    mpz_init(n);
    mpz_mod(n, N, P);

    mpz_t z;
    mpz_init_set_ui(z, 2);

    mpz_t c;
    mpz_init(c);

    mpz_t R;
    mpz_init(R);

    mpz_t t;
    mpz_init(t);

    mpz_t t0;
    mpz_init(t0);

    mpz_t b;
    mpz_init(b);

    if (p == 2) {
        return mpz_get_ui(n);
    }

    int Q = p-1;
    int S = 0;
    while (Q%2 == 0) {
        S++;
        Q = Q/2;
    }

    int M = S;

    while (mpz_legendre(z, P) != -1) {
        mpz_add_ui(z, z, 1);
    }
    mpz_powm_ui(c, z, Q, P);
    mpz_powm_ui(R, n, (Q+1)/2, P);
    mpz_powm_ui(t, n, Q, P);

    while (1 == 1) {
        if (mpz_cmp_ui(t, 1) == 0) {
            mpz_clear(P);
            mpz_clear(n);
            mpz_clear(z);
            mpz_clear(c);
            mpz_clear(t);
            mpz_clear(t0);
            mpz_clear(b);
            unsigned long int RR = mpz_get_ui(R);
            mpz_clear(R);
            return RR;
        }
        int i = 0;
        mpz_set(t0, t);
        while(mpz_cmp_ui(t, 1) != 0) {
            mpz_pow_ui(t, t, 2);
            mpz_mod(t, t, P);
            i++;
        }
        mpz_set(b, c);
        int j;
        for (j = 0; j < M - i - 1; j++) {
            mpz_pow_ui(b, b, 2);
            mpz_mod(b, b, P);
        }
        mpz_mul(R, R, b);
        mpz_mod(R, R, P);
        mpz_mul(c, b, b);
        mpz_mod(c, c, P);
        mpz_mul(t, t0, c);
        mpz_mod(t, t, P);
        M = i;
    }



}


//assign to no the value of the euclidean norm of a
void norm(mpz_t no, mpz_t a[3]) {
    
    mpz_mul(no, a[0], a[0]);
    mpz_addmul(no, a[1], a[1]);
    mpz_addmul(no, a[2], a[2]);

}


//a single step of the LLL reduction
void reduction(mpz_t k, mpz_t a[3], mpz_t b[3]) {

    mpz_t n;
    mpz_init(n);
    mpz_mul(n, a[0], b[0]);
    mpz_addmul(n, a[1], b[1]);
    mpz_addmul(n, a[2], b[2]);

    mpz_t d;
    mpz_init(d);
    norm(d, b);

    mpz_fdiv_q(k, n, d);


    mpz_t X1[3];
    int i;
    for (i = 0; i < 3; i++) {
        mpz_init_set(X1[i], a[i]);
        mpz_submul(X1[i], k, b[i]);
    }
    mpz_t norm1;
    mpz_init(norm1);
    norm(norm1, X1);

    mpz_t X2[3];
    for (i = 0; i < 3; i++) {
        mpz_init_set(X2[i], X1[i]);
        mpz_sub(X2[i], X2[i], b[i]);
    }
    mpz_t norm2;
    mpz_init(norm2);
    norm(norm2, X2);

    if (mpz_cmp(norm1, norm2) > 0) {
        mpz_add_ui(k, k, 1);
        for (i = 0; i < 3; i++) {
            mpz_set(a[i], X2[i]);
        }
    }
    else {
        for (i = 0; i < 3; i++) {
            mpz_set(a[i], X1[i]);
        }
    }

    mpz_clear(n);
    mpz_clear(d);
    for (i = 0; i < 3; i++) {
        mpz_clear(X1[i]);
        mpz_clear(X2[i]);
    }
    mpz_clear(norm1);
    mpz_clear(norm2);

}


//the complete LLL reduction on two mpz_t vectors
void LLL(mpz_t a[3], mpz_t b[3]) {

    mpz_t k;
    mpz_init_set_ui(k, 1);

    mpz_t l;
    mpz_init_set_ui(l, 1);

    while ( (mpz_cmp_si(k, 0) != 0) || (mpz_cmp_si(l, 0) != 0) ) {

        reduction(k, a, b);
        reduction(l, b, a);
    }

    mpz_clear(k);
    mpz_clear(l);

}


//assigns to a and b the coefficients of the two polynomials found via Montgomery's Two Quadratics Method,
//given n and p
void TQponctuel(mpz_t N, unsigned long int p, mpz_t a[3], mpz_t b[3], int i) {
assert((i == 0)||(i == 1));
  
    mpz_t P;
    mpz_init_set_ui(P, p);

    mpz_t z;
    mpz_init(z);
    mpz_sqrt(z, N);  

    
    //We compute the vectors a, b, c the way Peter Montgomery explained, starting from a root of n modulo p
    //These vectors span a lattice of Z^3

    mpz_t c0;
    mpz_init_set_ui(c0, p);

    mpz_t c1;
    mpz_init_set_ui(c1, modular_square_root(N, p));
    if (i == 1) {mpz_neg(c1, c1);}

    mpz_t X;
    mpz_init(X);
    mpz_add_ui(X, z, p/2);
    mpz_sub(X, X, c1);
    mpz_fdiv_q_ui(X, X, p);
    mpz_mul_ui(X, X, p);
    mpz_add(c1, c1, X);

    mpz_t c2;
    mpz_init(c2);
    mpz_mul(c2, c1, c1);
    mpz_sub(c2, c2, N);
    mpz_fdiv_q_ui(c2, c2, p);

    mpz_t c[3];
    mpz_init_set(c[0], c0);
    mpz_init_set(c[1], c1);
    mpz_init_set(c[2], c2);

    mpz_set(a[0], c1);
    mpz_set_si(a[1], -p);
    mpz_set_si(a[2], 0);

    mpz_invert(X, c1, P);
    mpz_mul(X, X, c2);
    mpz_mod(X, X, P);

    mpz_t b1;
    mpz_init_set(b1, X);
    mpz_neg(b1, b1);

    mpz_mul(X, c1, X);
    mpz_sub(X, X, c2);
    mpz_fdiv_q_ui(X, X, p);

    mpz_t b0;
    mpz_init_set(b0, X);

    mpz_set(b[0], b0);
    mpz_set(b[1], b1);
    mpz_set_ui(b[2], 1);

   //We then perform a LLL reduction on a and b

    LLL(a, b);

    mpz_clear(c0);
    mpz_clear(c1);
    mpz_clear(c2);
    mpz_clear(b0);
    mpz_clear(b1);
    mpz_clear(X);
    mpz_clear(z);
    mpz_clear(P);
    for (i = 0; i < 3; i++) {
        mpz_clear(c[i]);
    }

}


//TTT(x, a[3], b[3]) = 4/1575*(
//                             175*a0^2*b0^2 
//                           + 75*(2*a0^2*b0*b2 + a0^2*b1^2 + 4*a0*a1*b0*b1 + (2*a0*a2 + a1^2)*b0^2)*x
//                           + 63*(a0^2*b2^2 + 4*a1*a2*b0*b1 + a2^2*b0^2 + (2*a0*a2 + a1^2)*b1^2 + 2*(2*a0*a1*b1 + (2*a0*a2 + a1^2)*b0)*b2)*x^2 
//                           + 75*(a2^2*b1^2 + (2*a0*a2 + a1^2)*b2^2 + 2*(2*a1*a2*b1 + a2^2*b0)*b2)*x^3 
//                           + 175*a2^2*b2^2*x^4 
//                             )/x^2
//                   = 4/1575*(Y[0] + Y[1]*x + Y[2]*x^2 + Y[3]*x^3 + Y[4]*x^4)/x^2
void TTT(mpf_t resultat, mpf_t x, mpz_t a[3], mpz_t b[3]) {

    mpz_t a0;
    mpz_init_set(a0, a[0]);
    mpz_t a1;
    mpz_init_set(a1, a[1]);
    mpz_t a2;
    mpz_init_set(a2, a[2]);

    mpz_t b0;
    mpz_init_set(b0, b[0]);
    mpz_t b1;
    mpz_init_set(b1, b[1]);
    mpz_t b2;
    mpz_init_set(b2, b[2]);

    mpz_t Y[5];
    int i = 0;
    for (i = 0; i < 5; i++) {
        mpz_init(Y[i]);
    }

    mpz_t Z;
    mpz_init(Z);


    ////Y[0]

    mpz_set_ui(Y[0], 175);
    mpz_mul(Y[0], Y[0], a0);
    mpz_mul(Y[0], Y[0], a0);
    mpz_mul(Y[0], Y[0], b0);
    mpz_mul(Y[0], Y[0], b0);

    ////Y[1]

    mpz_mul(Z, a0, a0);
    mpz_mul(Z, Z, b0);
    mpz_mul(Z, Z, b2);
    mpz_mul_ui(Z, Z, 2);

    mpz_set(Y[1], Z);

    mpz_mul(Z, a0, a0);
    mpz_mul(Z, Z, b1);
    mpz_mul(Z, Z, b1);

    mpz_add(Y[1], Y[1], Z);

    mpz_set_ui(Z, 4);
    mpz_mul(Z, Z, a0);
    mpz_mul(Z, Z, a1);
    mpz_mul(Z, Z, b0);
    mpz_mul(Z, Z, b1);

    mpz_add(Y[1], Y[1], Z);

    mpz_mul(Z, a0, a2);
    mpz_mul_ui(Z, Z, 2);
    mpz_addmul(Z, a1, a1);
    mpz_t macro;
    mpz_init_set(macro, Z);  //macro = 2*a0*a2 + a1^2
    mpz_mul(Z, Z, b0);
    mpz_mul(Z, Z, b0);

    mpz_add(Y[1], Y[1], Z);

    mpz_mul_ui(Y[1], Y[1], 75);

    ////Y[2]

    mpz_mul(Z, a0, a1);
    mpz_mul(Z, Z, b1);
    mpz_mul_ui(Z, Z, 2);
    mpz_addmul(Z, macro, b0);
    mpz_mul(Z, Z, b2);
    mpz_mul_ui(Z, Z, 2);

    mpz_set(Y[2], Z);

    mpz_mul(Z, b1, b1);
    mpz_mul(Z, Z, macro);

    mpz_add(Y[2], Y[2], Z);

    mpz_mul(Z, a2, b0);
    mpz_mul(Z, Z, Z);

    mpz_add(Y[2], Y[2], Z);

    mpz_mul(Z, a1, a2);
    mpz_mul(Z, Z, b0);
    mpz_mul(Z, Z, b1);
    mpz_mul_ui(Z, Z, 4);

    mpz_add(Y[2], Y[2], Z);

    mpz_mul(Z, a0, b2);
    mpz_mul(Z, Z, Z);

    mpz_add(Y[2], Y[2], Z);

    mpz_mul_ui(Y[2], Y[2], 63);


    ////Y[3]

    mpz_mul(Z, a1, b1);
    mpz_mul_ui(Z, Z, 2);
    mpz_addmul(Z, a2, b0);
    mpz_mul(Z, Z, a2);
    mpz_mul(Z, Z, b2);
    mpz_mul_ui(Z, Z, 2);

    mpz_set(Y[3], Z);

    mpz_mul(Z, b2, b2);
    mpz_mul(Z, Z, macro);

    mpz_add(Y[3], Y[3], Z);

    mpz_mul(Z, a2, b1);
    mpz_mul(Z, Z, Z);

    mpz_add(Y[3], Y[3], Z);

    mpz_mul_ui(Y[3], Y[3], 75);

    ////Y[4]

    mpz_mul(Z, a2,  b2);
    mpz_mul(Z, Z, Z);
    mpz_mul_ui(Z, Z, 175);

    mpz_set(Y[4], Z);


    ////Y[0] + Y[1]*x + Y[2]*x^2 + Y[3]*x^3 + Y[4]*x^4

    mpf_t F;
    mpf_init(F);

    mpf_t RES;
    mpf_init(RES);

    mpf_set_z(F, Y[4]);
    mpf_set(RES, F);

    for(i = 0; i < 4; i++) {
        mpf_mul(RES, RES, x);
        mpf_set_z(F, Y[3-i]);
        mpf_add(RES, RES, F);
    }

    ////TTT(x, a, b)

    mpf_mul_ui(RES, RES, 4);
    mpf_div_ui(RES, RES, 1575);
    mpf_div(RES, RES, x);
    mpf_div(RES, RES, x);

    mpf_set(resultat, RES);


    mpz_clear(a0);
    mpz_clear(a1);
    mpz_clear(a2);
    mpz_clear(b0);
    mpz_clear(b1);
    mpz_clear(b2);
    mpz_clear(Z);
    for(i = 0; i < 5; i++) {
        mpz_clear(Y[i]);
    }

    mpf_clear(F);
    mpf_clear(RES);
    mpz_clear(macro);

}


/// norme := min TTT(x, a, b)  (a and b fixed)
/// s0 := argmin(TTT(x, a, b)) 
void normeL2(mpf_t s0, mpf_t norme, mpz_t a[3], mpz_t b[3]) {
    
    mpf_t x0;
    mpf_t x1;
    mpf_t x2;
    mpf_t x3;

    mpf_init_set_d(x0, 0.000001);
    mpf_init_set_d(x1, 1);
    mpf_init_set_d(x2, 1);
    mpf_init_set_d(x3, 2);

    mpf_t T0;
    mpf_t T1;
    mpf_t T2;
    mpf_t T3;

    mpf_init(T0);
    mpf_init(T1);
    mpf_init(T2);
    mpf_init(T3);

    TTT(T0, x0, a, b);
    TTT(T1, x1, a, b);
    TTT(T2, x2, a, b);
    TTT(T3, x3, a, b);


    while(mpf_cmp(T2, T3) > 0) {
        mpf_set(x2, x3);
        mpf_set(T2, T3);
        mpf_add(x3, x3, x3);
        TTT(T3, x3, a, b);
    }

    mpf_t ecart;
    mpf_init(ecart);
    mpf_sub(ecart, x3, x0);

    mpf_t F;
    mpf_init(F);


    while(mpf_cmp_d(ecart, 0.001) > 0) {
        mpf_div_ui(F, ecart, 3);
        mpf_add(x1, x0, F);
        TTT(T1, x1, a, b);
        mpf_add(x2, x1, F);
        TTT(T2, x2, a, b);

        int j = 0;
        TTT(F, x0, a, b);

        if (mpf_cmp(F, T1) > 0) {
            j = 1;
            mpf_set(F, T1);
        }
        if (mpf_cmp(F, T2) > 0) {
            j = 2;
            mpf_set(F, T2);
        }
        if (mpf_cmp(F, T3) > 0) {
            j = 3;
            mpf_set(F, T3);
        }

        if (j == 0) {
            mpf_set(x3, x1);
        }
        if (j == 1) {
            mpf_set(x3, x2);
        }
        if (j == 2) {
            mpf_set(x0, x1);
        }
        if (j == 3) {
            mpf_set(x0, x2);
        }
        mpf_sub(ecart, x3, x0);
    }

    mpf_sqrt(s0, x0);
    mpf_sqrt(s0, s0);

    TTT(norme, x0, a, b);

    mpf_clear(x0);
    mpf_clear(x1);
    mpf_clear(x2);
    mpf_clear(x3);
    mpf_clear(T0);
    mpf_clear(T1);
    mpf_clear(T2);
    mpf_clear(T3);
    mpf_clear(ecart);
    mpf_clear(F);

}

//x := exp(z)
void expo(mpf_t x, mpf_t z) {

    mpf_t x0;
    mpf_init_set_ui(x0, 0);

    mpf_t zz;
    mpf_init_set_ui(zz, 1);

    int i;
    for (i = 1; i < 100; i++ ) {
        mpf_add(x0, x0, zz);
        mpf_mul(zz, zz, z);
        mpf_div_ui(zz, zz, i);
    }

    mpf_set(x, x0);

    mpf_clear(x0);
    mpf_clear(zz);

}


// note := exp(alpha(a,2000))*exp(alpha(b,2000))*norme(a,b)     (cf. normeL2)
void notation(mpf_t note, mpz_t a[3], mpz_t b[3], mpf_t norme) {

    mpf_t alpha_a;
    double aa = get_alpha(a, 2, 2000);
    mpf_init_set_d(alpha_a, aa);

    mpf_t alpha_b;
    double bb = get_alpha(b, 2, 2000);
    mpf_init_set_d(alpha_b, bb);

    mpf_t X;
    mpf_init(X);
    mpf_sqrt(X, norme);

    mpf_t exp_a;
    mpf_init(exp_a);

    mpf_t exp_b;
    mpf_init(exp_b);

    expo(exp_a, alpha_a);
    expo(exp_b, alpha_b);

    mpf_mul(X, X, exp_a);
    mpf_mul(X, X, exp_b);

    mpf_set(note, X);

    mpf_clear(alpha_a);
    mpf_clear(alpha_b);
    mpf_clear(X);
    mpf_clear(exp_a);
    mpf_clear(exp_b);

}



//Given a big number N and an integer k,
//tests k pairs of polynomials and print the pair having the best rating,
//along with its skewness, norm, alpha(f), alpha(g), and the prime numer which gave it
void
twoquadratics (mpz_t N, unsigned long int k)
{
    unsigned int compteur = 0;
    unsigned long int p = 2;
    
    mpz_t P;
    mpz_init_set_ui(P, p);

    mpz_t a[3];
    mpz_t b[3];
    mpz_t a0[3];
    mpz_t b0[3];

    mpf_t no1;
    mpf_init(no1);

    mpf_t no;
    mpf_init(no);

    mpf_t s1;
    mpf_init(s1);

    mpf_t s0;
    mpf_init(s0);

    mpf_t alpha_a;
    mpf_init(alpha_a);

    mpf_t alpha_b;
    mpf_init(alpha_b);

    mpf_t note0;
    mpf_init(note0);

    mpf_t note1;
    mpf_init(note1);

    mpz_t P0;
    mpz_init(P0);


    int i;
    for (i = 0; i < 3; i++ ) {
        mpz_init(a[i]);
        mpz_init(b[i]);
        mpz_init(a0[i]);
        mpz_init(b0[i]);
    }

    while(compteur < k) {
        while(mpz_kronecker(N, P) != 1) {
            mpz_nextprime(P, P);
        }
        p = mpz_get_ui(P);
        TQponctuel(N, p, a, b, 0);
        //gmp_printf ("\n\na = %Zd + %Zd * X + %Zd * X^2 \n", a[0], a[1], a[2]);
        //gmp_printf ("b = %Zd + %Zd * X + %Zd * X^2 \n\n", b[0], b[1], b[2]);
        normeL2(s1, no1, a, b);
        notation(note1, a, b, no1);
        if(compteur == 0) {
            mpz_set(a0[0], a[0]);
            mpz_set(a0[1], a[1]);
            mpz_set(a0[2], a[2]);
            mpz_set(b0[0], b[0]);
            mpz_set(b0[1], b[1]);
            mpz_set(b0[2], b[2]);
            mpf_set(s0, s1);
            mpf_set(no, no1);
            mpf_set(note0, note1);
        }
        //printf ("note = ");
        //printscientific(note1);
        //printf("\n");
        if (mpf_cmp(note0, note1) > 0) {
            mpz_set(a0[0], a[0]);
            mpz_set(a0[1], a[1]);
            mpz_set(a0[2], a[2]);
            mpz_set(b0[0], b[0]);
            mpz_set(b0[1], b[1]);
            mpz_set(b0[2], b[2]);
            mpf_set(s0, s1);
            mpf_set(note0, note1);
            mpz_set(P0, P);
        }
        TQponctuel(N, p, a, b, 1);
        //gmp_printf ("\na = %Zd + %Zd * X + %Zd * X^2 \n", a[0], a[1], a[2]);
        //gmp_printf ("b = %Zd + %Zd * X + %Zd * X^2 \n\n", b[0], b[1], b[2]);
        normeL2(s1, no1, a, b);
        notation(note1, a, b, no);
        if (mpf_cmp(note0, note1) > 0) {
            mpz_set(a0[0], a[0]);
            mpz_set(a0[1], a[1]);
            mpz_set(a0[2], a[2]);
            mpz_set(b0[0], b[0]);
            mpz_set(b0[1], b[1]);
            mpz_set(b0[2], b[2]);
            mpf_set(s0, s1);
            mpf_set(note0, note1);
            mpz_set(P0, P);
        }
        //printf ("note = ");
        //printscientific(note1);
        //printf("\n");
        mpz_nextprime(P, P);
        compteur += 2;
        //printf("\n%d\n", compteur);
    }

    printf("\n\nbest polynomials found:\n");
    gmp_printf ("a = %Zd + %Zd * X + %Zd * X^2 \n", a0[0], a0[1], a0[2]);
    gmp_printf ("b = %Zd + %Zd * X + %Zd * X^2 \n\n", b0[0], b0[1], b0[2]);

    printf ("note = ");
    printscientific(note0);
    gmp_printf ("\ns0 = %Ff\n", s0);
    double aaa = get_alpha(a0, 2, 2000);
    printf ("alpha(a) = %f\n", aaa);
    double bbb = get_alpha(b0, 2, 2000);
    printf ("alpha(b) = %f\n", bbb);
    normeL2(s0, no, a0, b0);
    printf ("norme = ");
    printscientific(no);
    printf("\n");
    gmp_printf ("P = %Zd\n", P0);


    mpf_clear(no1);
    mpf_clear(no);
    mpf_clear(s1);
    mpf_clear(s0);
    mpf_clear(alpha_a);
    mpf_clear(alpha_b);
    mpf_clear(note0);
    mpf_clear(note1);
    mpz_clear(P0);
    mpz_clear(P);

    for (i = 0; i < 3; i++ ) {
        mpz_clear(a[i]);
        mpz_clear(b[i]);
        mpz_clear(a0[i]);
        mpz_clear(b0[i]);
    }

}

void
usage ()
{
  fprintf (stderr, "Usage: twoquadratics [-k nnn] [N]\n");
  exit (EXIT_FAILURE);
}

int
main (int argc, char *argv[])
{
    mpz_t N;
    int k = 100;

    while (argc >= 2 && argv[1][0] == '-')
      {
        if (argc >= 3 && strcmp (argv[1], "-k") == 0)
          {
            k = atoi (argv[2]);
            argc -= 2;
            argv += 2;
          }
        else
          {
            fprintf (stderr, "Unknown option: %s\n", argv[1]);
            usage ();
          }
      }

    if (argc == 1) /* test on the c59 */
      mpz_init_set_str (N, "71641520761751435455133616475667090434063332228247871795429", 10);
    else
      mpz_init_set_str (N, argv[1], 10);

    twoquadratics (N, k);

    mpz_clear (N);

    return 0;
}

