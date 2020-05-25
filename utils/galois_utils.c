#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <string.h>
#include "mod_ul.h"     // modulusul_t
#include "galois_utils.h"
#include "macros.h" // ASSERT_ALWAYS

void automorphism_init(int *ord, int mat[4], const char *galois_autom){
    int A, B, C, D; // x -> (A*x+B)/(C*x+D)

    if(strcmp(galois_autom, "autom2.1") == 0 
       || strcmp(galois_autom, "autom2.1g") == 0
       || strcmp(galois_autom, "1/y") == 0){
	*ord = 2; A = 0; B = 1; C = 1; D = 0;
    }
    else if(strcmp(galois_autom, "autom2.2") == 0
	    || strcmp(galois_autom, "autom2.2g") == 0
	    || strcmp(galois_autom, "_y") == 0){
	*ord = 2; A = -1; B = 0; C = 0; D = 1;
    }
    else if(strcmp(galois_autom, "autom3.1") == 0
	    || strcmp(galois_autom, "autom3.1g") == 0){
	// 1-1/x = (x-1)/x
	*ord = 3; A = 1; B = -1; C = 1; D = 0;
    }
    else if(strcmp(galois_autom, "autom3.2") == 0){
	// -1-1/x = (-x-1)/x
	*ord = 3; A = -1; B = -1; C = 1; D = 0;
    }
    else if(strcmp(galois_autom, "autom4.1") == 0){
	// -(x+1)/(x-1)
	*ord = 4; A = -1; B = -1; C = 1; D = -1;
    }
    else if(strcmp(galois_autom, "autom6.1") == 0){
	// -(2*x+1)/(x-1)
	*ord = 6; A = -2; B = -1; C = 1; D = -1;
    }
    else{
	fprintf(stderr, "Unknown automorphism: %s\n", galois_autom);
	ASSERT_ALWAYS(0);
    }
    mat[0] = A; mat[1] = B; mat[2] = C; mat[3] = D;
}

/* OUTPUT: sigma(r).
   Since we can have r == qq, we need unsigned long's.
*/
unsigned long automorphism_apply(residueul_t mat[4], unsigned long r,
				 const modulusul_t mm, const unsigned long qq)
{
    residueul_t xx;
    unsigned long sigma_r;

    modul_init(xx, mm);
    if(r == qq){
	// sigma(oo) = A/C
	// TODO: test C w.r.t. +/- 1
	modul_inv(xx, mat[2], mm);
	modul_mul(xx, xx, mat[0], mm);
	sigma_r = modul_get_ul(xx, mm);
    }
    else{
	residueul_t rr;

	modul_init(rr, mm);
	modul_set_ul(rr, r, mm);
	// denominator: C*r+D
	modul_mul(xx, mat[2], rr, mm);
	modul_add(xx, xx, mat[3], mm);
	if(modul_is0(xx, mm))
	    // 1/0 = oo
	    sigma_r = qq;
	else{
	    residueul_t yy;

	    modul_init(yy, mm);
	    modul_inv(yy, xx, mm);
	    // numerator: A*r+B
	    modul_mul(xx, mat[0], rr, mm);
	    modul_add(xx, xx, mat[1], mm);
	    modul_mul(yy, yy, xx, mm);
	    sigma_r = modul_get_ul(yy, mm);
	    modul_clear(yy, mm);
	}
	modul_clear(rr, mm);
    }
    modul_clear(xx, mm);
    return sigma_r;
}
