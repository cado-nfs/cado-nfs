#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "utils.h"
//#include "poly.h"

#error "dead code"

#define DEBUG 0

void polymodF_from_ab(polymodF_t tmp, long a, unsigned long b) {
  tmp->v = 0;
  tmp->p->deg = 1;
  mpz_set_ui(tmp->p->coeff[1], b);
  mpz_neg(tmp->p->coeff[1], tmp->p->coeff[1]);
  mpz_set_si(tmp->p->coeff[0], a);
}

int main(int argc, char **argv)
{
    char *ratside = argv[2], *algside = argv[4], *magma = argv[5];
    FILE *abfile, *algfile;
    cado_poly pol;
    long a;
    unsigned long b;
    poly_t F;
    polymodF_t prd, tmp;
    int i, ret, nab;

    ret = read_polynomial(pol, argv[1]);
    ASSERT (ret == 1);

    abfile = fopen(argv[3], "r");
    ASSERT (abfile != NULL);

    algfile = fopen(algside, "w");
    ASSERT (algfile != NULL);

    // Init F to be the algebraic polynomial
    poly_alloc(F, pol->degree);
    for (i = pol->degree; i >= 0; --i)
	poly_setcoeff(F, i, pol->f[i]);
    
    // Init prd to 1.
    poly_alloc(prd->p, pol->degree);
    mpz_set_ui(prd->p->coeff[0], 1);
    prd->p->deg = 0;
    prd->v = 0;
    
    // Allocate tmp
    poly_alloc(tmp->p, 1);
    nab = 0;
    while(fscanf(abfile, "%ld %lu", &a, &b) != EOF){
	if(!(nab % 100000))
	    fprintf(stderr, "# Reading ab pair #%d at %2.2lf\n",nab,seconds());
	if((a == 0) && (b == 0))
	    break;
	polymodF_from_ab(tmp, a, b);
	polymodF_mul(prd, prd, tmp, F);
#if DEBUG >= 1
	fprintf(stderr, "// a = %ld, b = %lu => %d\n", a, b, prd->v);
	fprintf(stderr, "prodab:=1/fd^%d*(", prd->v);
	for (i = prd->p->deg; i >= 0; i--)
	    gmp_fprintf(stderr, "+x^%d*(%Zd)\n", i, prd->p->coeff[i]);
	fprintf(stderr, ");\n");
#endif
	nab++;
    }

    fprintf(stderr, "#(a, b) = %d, v = %d\n", nab, prd->v);
    fprintf(stderr, "sizeinbit of cst term = %d\n", (int)mpz_sizeinbase(prd->p->coeff[0], 2));
    fprintf(stderr, "sizeinbit of leading term = %d\n", (int)mpz_sizeinbase(prd->p->coeff[pol->degree-1], 2));
    
    for (i = 0; i < pol->degree; ++i) 
	gmp_fprintf(algfile, "%Zx\n", prd->p->coeff[i]);
    fprintf(algfile, "%d\n", prd->v);
    fclose(algfile);
    
    {
	FILE * formagma;
	mpz_t m;
	
	formagma = fopen(magma, "w");
	if (formagma == NULL) {
	    fprintf(stderr, "cannot create %s\n", magma);
	    return 1;
	}
	fprintf(formagma, "NFSRATFILE := \"%s\";\n", ratside);
	fprintf(formagma, "NFSALGFILE := \"%s\";\n", algside);
	fprintf(formagma, "f := ");
	for (i = 0; i <= pol->degree; ++i)
	    gmp_fprintf(formagma, "%+Zd*x^%d", pol->f[i], i);
	fprintf(formagma, ";\nm := ");
	mpz_init(m);
	mpz_neg(m, pol->g[0]);
	gmp_fprintf(formagma, "%Zd;\nN := %Zd;\n", m, pol->n);
	mpz_clear(m);
	fclose(formagma);
    }
    
    poly_free(F);
    poly_free(prd->p);
    poly_free(tmp->p);
    
    return 0;
}
