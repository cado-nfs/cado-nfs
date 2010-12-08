/* 
 * Program: free relations
 * Author : F. Morain
 * Purpose: creating free relations in a suitable format
 * 
 */

#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gmp.h>
#include "mod_ul.c"
#include "utils.h"
#include "sieve/fb.h"

#include "hashpair.h"

// The best reference for free relations is Huizing (Exp. Math. 1996,
// Section 4). Handling them correctly in the factorization is also detailed.
// If f(X)=cK*X^dK+... and g(X)=cL*X^dL+... are the two polynomials used,
// then p is a free prime iff f(X) and g(X) split completely mod p into
// dK (resp. dL) *distinct* linear factors. In other words, we have to check
// that p does not divide the discriminants of f and g, and p doesn't divide
// cK*cL.
//
// TODO: make this definition go into the actual code, which is not the
// case currently.

// When p1 == p2, sort r's in increasing order
int
compare_ul2(const void *v1, const void *v2)
{
    unsigned long w1 = *((unsigned long*) v1);
    unsigned long w2 = *((unsigned long*) v2);

    if(w1 > w2)
	return 1;
    if(w1 < w2)
	return -1;
    w1 = *(1 + (unsigned long*) v1);
    w2 = *(1 + (unsigned long*) v2);
    return (w1 >= w2 ? 1 : -1);
}

// find free relations by inspection.
int
findFreeRelations(hashtable_t *H, cado_poly pol, int nprimes)
{
    unsigned long *tmp = (unsigned long *)malloc((2+(nprimes<<1)) * sizeof(unsigned long));
    unsigned int i, j, k, ntmp = 0;
    int pdeg, nfree;

    for(i = 0; i < H->hashmod; i++)
	if(H->hashcount[i] > 0){
          tmp[ntmp++] = (unsigned long) GET_HASH_P(H,i);
          tmp[ntmp++] = GET_HASH_R(H,i);
	}
    qsort(tmp, (ntmp>>1), 2 * sizeof(unsigned long), compare_ul2);
    // add a sentinel
    tmp[ntmp] = 0;
    tmp[ntmp+1] = 0;
    // now, inspect
    for(i = 0; ; i += 2)
	if(((long)tmp[i+1]) >= 0)
	    break;
    j = i+2;
    pdeg = 1;
    nfree = 0;
    while(1){
	if(tmp[j] == tmp[i]){
	    // same prime
	    if(((long)tmp[j+1]) >= 0)
		// another algebraic factor
		pdeg++;
	}
	else{
	    // new prime
	    if(pdeg == pol[0].degree){
		// a free relation
		printf("%lu,0:%lx:", tmp[i], tmp[i]);
		for(k = i; k < i+(pdeg<<1); k += 2){
		    printf("%lx", tmp[k+1]);
		    if(k < i+(pdeg<<1)-2)
			printf(",");
		}
		printf("\n");
		nfree++;
	    }
	    i = j;
	    if(tmp[j] == 0)
		break;
	    pdeg = 1;
	}
	j += 2;
    }
    free(tmp);
    return nfree;
}

/* On largeFreeRelations:
 *
commit 7cb6c0e177ed7ba86c1768942bfb42f2e575522e
Author: morain <morain@de4796e2-0b1e-0410-91b4-86bc23549317>
Date:   Tue Feb 12 12:41:46 2008 +0000

    * Cleaned purge.c
    * added a new option to freerel to look for (small and large) free relations
    by inspection of all the relations found. Quite experimental at the time
    being. Do not use it.
    
    
    git-svn-id: svn+ssh://scm.gforge.inria.fr/svn/cado/trunk@801 de4796e2-0b1e-0
 */

void
largeFreeRelations (cado_poly pol, char **fic, int verbose)
{
    hashtable_t H;
    int Hsizea, nprimes_alg = 0, nfree = 0;
    int need64 = (pol->lpbr > 32) || (pol->lpba > 32);

    ASSERT(fic != NULL);
    /* The number of algebraic large primes is about 1/2*L/log(L)
       where L is the algebraic large prime bound, i.e., L=2^lpba.
       However since we store separately primes p and the corresponding root r
       of f mod p, the number of (p,r) pairs is about L/log(L). */
    Hsizea = (1 << pol[0].lpba) / ((int)((double) pol[0].lpba * log(2.0)));
    hashInit (&H, Hsizea, verbose, need64);
    if (verbose)
      fprintf (stderr, "Scanning relations\n");

    relation_stream rs;
    relation_stream_init(rs);
    // rs->sort_primes = 1;
    // rs->reduce_mod2 = 1;
    for( ; *fic ; fic++) {
        relation_stream_openfile(rs, *fic);
        if (verbose)
            fprintf (stderr, "Adding file %s\n", *fic);
        for( ; relation_stream_get(rs, NULL) >= 0 ; ) {
            if (rs->rel.b == 0) {
                fprintf(stderr, "Ignoring already found free relation...\n");
                continue;
            }
            relation_compress_rat_primes(&rs->rel);
            relation_compress_alg_primes(&rs->rel);
            reduce_exponents_mod2(&rs->rel);
            computeroots(&rs->rel);
            for(int j = 0; j < rs->rel.nb_ap; j++){
                int h = hashInsert(&H, rs->rel.ap[j].p, rs->rel.ap[j].r);
                nprimes_alg += (H.hashcount[h] == 1); // new prime
            }
        }
        relation_stream_closefile(rs);
        if (verbose)
            hashCheck (&H);
    }
    relation_stream_clear(rs);

    if (verbose)
      fprintf (stderr, "nprimes_alg = %d\n", nprimes_alg);
    nfree = findFreeRelations(&H, pol, nprimes_alg);
    if (verbose)
      fprintf (stderr, "Found %d usable free relations\n", nfree);
    hashFree(&H);
}

int
countFreeRelations(int *deg, char *roots)
{
    FILE *ifile = fopen(roots, "r");
    char str[1024], str0[128], str1[128], str2[128], *t;
    int nfree = 0, nroots;

    *deg = -1;
    // look for the degree
    while(1){
	if(!fgets(str, 1024, ifile))
	    break;
	sscanf(str, "%s %s %s", str0, str1, str2);
	if((str0[0] == '#') && !strcmp(str1, "DEGREE:")){
	    *deg = atoi(str2);
	    break;
	}
    }
    if(*deg == -1){
	fprintf(stderr, "No degree found, sorry\n");
	exit(1);
    }
    fprintf(stderr, "DEGREE = %d\n", *deg);
    while(1){
	if(!fgets(str, 1024, ifile))
            break;
	// format is "11: 2,6,10"
	sscanf(str, "%s %s", str0, str1);
	nroots = 1;
	for(t = str1; *t != '\0'; t++)
	    if(*t == ',')
		nroots++;
	if(nroots == *deg)
	    nfree++;
    }
    fclose(ifile);
    return nfree;
}

void
addFreeRelations(char *roots, int deg)
{
    factorbase_degn_t *fb = fb_read(roots, 1.0, 0), *fbptr;
    fbprime_t p;
    int i;

    for(fbptr = fb; fbptr->p != FB_END; fbptr = fb_next(fbptr)){
	p = fbptr->p;
	if(fbptr->nr_roots == deg){
	    // free relation, add it!
	    // (p, 0) -> p-0*m
	    printf("%d,0:%lx:", p, (long)p);
	    for(i = 0; i < fbptr->nr_roots; i++){
		printf("%x", fbptr->roots[i]);
		if(i < fbptr->nr_roots-1)
		    printf(",");
	    }
	    printf("\n");
	}
    }
}

static void
smallFreeRelations(char *fbfilename)
{
    int deg;
    int nfree = countFreeRelations (&deg, fbfilename);
    fprintf (stderr, "# Free relations: %d\n", nfree);
    
    fprintf (stderr, "Handling free relations...\n");
    addFreeRelations (fbfilename, deg);
}

static void
usage (char *argv0)
{
  fprintf (stderr, "Usage: %s -fb xxx.roots\n", argv0);
  fprintf (stderr, "or     %s [-v] -poly xxx.poly xxx.rels1 xxx.rels2 ... xxx.relsk\n", argv0);
  exit (1);
}

int
main(int argc, char *argv[])
{
    char *fbfilename = NULL, *polyfilename = NULL, **fic;
    char *argv0 = argv[0];
    cado_poly cpoly;
    int nfic = 0;
    int verbose = 0;

    while (argc > 1 && argv[1][0] == '-')
      {
        if (argc > 2 && strcmp (argv[1], "-fb") == 0)
          {
            fbfilename = argv[2];
            argc -= 2;
            argv += 2;
          }
        else if (argc > 1 && strcmp (argv[1], "-v") == 0)
          {
            verbose ++;
            argc --;
            argv ++;
          }
        else if (argc > 2 && strcmp (argv[1], "-poly") == 0)
          {
            polyfilename = argv[2];
            argc -= 2;
            argv += 2;
          }
        else
          usage (argv0);
      }

    fic = argv+1;
    nfic = argc-1;

    if (polyfilename == NULL || ((fbfilename == NULL) && (nfic == 0)))
      usage (argv0);

    cado_poly_init(cpoly);
    if (!cado_poly_read (cpoly, polyfilename))
      {
        fprintf (stderr, "Error reading polynomial file\n");
        exit (EXIT_FAILURE);
      }

    /* check that n divides Res(f,g) [might be useful to factor n...] */
    check_polynomials (cpoly);

#if 0
    if (mpz_cmp_ui (cpoly->g[1], 1) != 0)
      {
        fprintf (stderr, "Error, non-monic linear polynomial not yet treated");
        fprintf (stderr, " (more theory needed)\n");
        exit (1);
      }
#endif

    if (nfic == 0)
	smallFreeRelations(fbfilename);
    else
      largeFreeRelations(cpoly, fic, verbose);

    cado_poly_clear (cpoly);

    return 0;
}
