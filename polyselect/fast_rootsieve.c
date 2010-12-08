#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "gmp_aux.h"
#include "utils.h"
#include "auxiliary.h"
#include "mod_ul_default.h"
#include "mod_ul.h"
#include "fast_rootsieve.h"
#include "macros.h"
#define COMP 0
#define VERBOSE 0
#define TRACE 0
#define DEBUG 0
#define DUMP 0
#define PLB 2
#define PUB 13
#define V_LBOUND 0
#define V_UBOUND 289
#define U_LBOUND 0
#define U_UBOUND 837
#define MAXPOW 20
#define CUTOFF 1.0e-5
#define REALLOC_INCREMENT 300
#define EXPECTED_NB_OF_LATTICES 200

/* Rootsieve : we compute the root property for every polynomial of 
   the form f(x)+(u*x+v)*g(x) for u,v in a given (rectangular) range,
   and for given polynomials f and g.

   The code is a C port of Emmanuel Thomé's code in Sage.
*/

// We could eliminate the following by using mod_ul residues directly
unsigned long modular_inverse( unsigned long r, unsigned long N );
unsigned long modular_product( unsigned long a, unsigned long b, unsigned long N ); 

// Rootsieve - Specific
long rotation_inner(rootsieve_dictionary rdict, unsigned long p,poly_t ff,poly_t gg,poly_t fdg_gdf,long u0,long v0,long l0,int ld,int m,long twist_v1,double scale,poly_t dphi);
long rootsieve_handle_p(rootsieve_dictionary rdict, unsigned long p);
long light (double * alpha, long u0, long v0, long us, long vs, long skew, long um, long vm, double contribution);
void extend_ppow_list(rootsieve_dictionary rdict,int newlength, unsigned long p); 
long add_lattice(rootsieve_dictionary rdict, long u0, long v0, long us, long vs, long skew, long U0, long U1, long V0, long V1, double contribution, unsigned long p);
int empty(long u0, long v0, long us, long vs, long skew, long U0, long U1, long V0, long V1);

// Polynomials
void poly_set_identity(poly_t f); 
void compose_psi(poly_t f, const poly_t g,unsigned long l,mpz_t * ppow,int maxpow);
void compose_reduce(poly_t f, const poly_t g,unsigned long l,mpz_t * ppow,int maxpow, mpz_t m);
void print_dictionary(rootsieve_dictionary rdict);

/*
int main() {
  // This main function just runs an example.
  poly_t f,g;

  poly_alloc(f,5);
  poly_alloc(g,1);

  poly_setcoeff_str(f,5,"1008593880",10);
  poly_setcoeff_str(f,4,"115824918113473",10);
  poly_setcoeff_str(f,3,"-78935188484415054956",10);
  poly_setcoeff_str(f,2,"-57788662352740339351859805969",10);
  poly_setcoeff_str(f,1,"-2655040283375779850204251805625017",10);
  poly_setcoeff_str(f,0,"-30504936769336291824389723279850895305",10);

  poly_setcoeff_str(g,1,"1628876881135933",10);
  poly_setcoeff_str(g,0,"161423410516092787320054348810",10);  

  print_lattices(rootsieve(f->coeff,5,g->coeff[0],g->coeff[1],U_LBOUND,U_UBOUND,V_LBOUND,V_UBOUND,PUB,VERBOSE));     
  
  poly_free(f);
  poly_free(g);
  }*/


lattice_list rootsieve(mpz_t * f_coeffs, int degf, mpz_t g_coeff1, mpz_t g_coeff0, long U0, long U1, long V0, long V1, unsigned long prime_bound) {  
  
  poly_t f,g;

  poly_alloc(f,degf);
  poly_set(f,f_coeffs,degf);

  poly_alloc(g,1);
  poly_setcoeff(g,0,g_coeff0);
  poly_setcoeff(g,1,g_coeff1);
  mpz_neg(g->coeff[0],g->coeff[0]);

  if(VERBOSE) {
      printf("Received polynomials f ...\n");
      poly_print(f);
      printf("and g ... \n");
      poly_print(g);
  }
  
  // Counters  
  unsigned long p;  
  
  // We prepare the dictionary
  rootsieve_dictionary rdict;

  // We prepare the lattice list
  rdict->lattice_alloc = EXPECTED_NB_OF_LATTICES;
  rdict->lattices = malloc(rdict->lattice_alloc*sizeof(rootsieve_lattice));
  rdict->L = 0;  
  
  // The rectangle
  rdict->U0 = U0;
  rdict->U1 = U1;
  rdict->V0 = V0;
  rdict->V1 = V1;

  if(VERBOSE) {
    printf("Rotation takes place in the rectangle of lower left corner (%ld,%ld) and upper right corner (%ld,%ld) ...\n",U0,V0,U1,V1);
  }
      
  // Contributions smaller than cutoff are neglected (ie, we stop the recursion)
  rdict->cutoff = CUTOFF;  

  // This gives the exceptional values for u (once divided by g(l)^2)
  poly_t dg,df,aux;  
  poly_alloc(dg,g->deg-1);
  poly_alloc(df,f->deg-1);  
  poly_alloc(aux,-1);
  poly_alloc(rdict->fdg_gdf,-1);
  poly_derivative(dg,g);
  poly_derivative(df,f);

  if(DEBUG) {
    printf("Computed derivative of g ...\n");
    poly_print(dg);
    printf("Computed derivative of f ...\n");
    poly_print(df);    
  }

  poly_mul(rdict->fdg_gdf,f,dg);
  poly_mul(aux,g,df);  
  poly_sub(rdict->fdg_gdf,rdict->fdg_gdf,aux);

  if(DEBUG) {
      printf("Computed f*dg - g*df ...\n");
      poly_print(rdict->fdg_gdf);
  }

  poly_free(aux);
  poly_free(dg);
  poly_free(df);
  
  poly_alloc(rdict->f,f->deg);
  poly_copy(rdict->f,f);
  poly_alloc(rdict->g,g->deg);
  poly_copy(rdict->g,g);

  if(VERBOSE) {
    printf("Added f and g to the dictionary...\n");
    printf("Now entering prime number loop...\n");
  }
    
  for ( p = PLB ; p <= prime_bound ; p += 1 + (1 & p)) {
    if (isprime(p)) {

      if(VERBOSE)
	printf("Finding lattices for prime: %lu\n",p);      

      rootsieve_handle_p(rdict,p);        
    }
  }
  
  poly_free(rdict->f);
  poly_free(rdict->g);
  poly_free(rdict->fdg_gdf);    

  // We pack and return the list of lattices
  lattice_list ll;
  ll.lattices = rdict->lattices;
  ll.length = rdict->L;
  return ll;

}


long rootsieve_handle_p(rootsieve_dictionary rdict,unsigned long p) { // Ready
  

  //"""Fast root sieve modulo p and its powers. The rdict argument must
  // have been setup beforehand"""

  // Counters
  int i;
  long hits;

  long u0,v0,l0;
  double scale;
  int ld,m;
  poly_t ff,gg,fdg_gdf;
  
  if(DEBUG) {    
    printf("Creating polynomials ff and gg ...\n");
  }
  poly_alloc(ff,rdict->f->deg);
  poly_copy(ff,rdict->f);
  
  if(DEBUG) {    
    printf("Polynomial ff ...\n");
    poly_print(ff);
  }
  poly_alloc(gg,rdict->g->deg);
  poly_copy(gg,rdict->g);

  if(DEBUG) {    
    printf("Polynomial gg ...\n");
    poly_print(gg);
  }  

  poly_alloc(fdg_gdf,rdict->fdg_gdf->deg);  
  scale = log(p)*(double)p/(double)(p+1);

  if(DEBUG) {
    printf("The scale of contribution is %f\n",scale);
  }
  
  // The rectanglewise upper bound for the powers of prime p.  
  rdict->m_max = (int)floor(log(MAX(labs(rdict->U1-rdict->U0),labs(rdict->V1-rdict->V0)))/log(p));

  if(DEBUG) {
    printf("The rectangle upper bound is %d ...\n",rdict->m_max);
  }

  l0 = 0;
  ld = m = 0;  
  // In theory, one does not need to initialize these values,
  // we do it now to control warnings.
  u0 = 0; 
  v0 = 0; 
  
  
  // We must be able to compose polynomials with psi(x) = l+p*x, that's the reason of taking the max of the degrees.
  //rdict->maxpow = MAX(MAX(rdict->m_max +1,rdict->g->deg),MAX(rdict->f->deg,rdict->fdg_gdf->deg)); 
  rdict->maxpow = MAXPOW;  

  if(TRACE) {
    printf("\nrdict->maxpow=%d,rdict->m_max=%d\n",rdict->maxpow,rdict->m_max);
  }

  ASSERT_ALWAYS(rdict->maxpow >= rdict->m_max);

  if(DEBUG) {    
    printf("Maximum power of %lu is %d ...\n",p,rdict->maxpow);
  }

  rdict->ppow = malloc((rdict->maxpow+1)*sizeof(mpz_t));

  mpz_init_set_ui(rdict->ppow[0],1);
  for( i=1 ; i<=rdict->maxpow; i++) {
    mpz_init(rdict->ppow[i]);
    mpz_mul_ui(rdict->ppow[i],rdict->ppow[i-1],p);
  }

  if(DEBUG) {    
    printf("Initialized powers of %lu ...\n",p);
  }
  
  poly_alloc(rdict->gmodp,rdict->g->deg);
  poly_reduce_mod_mpz(rdict->gmodp,rdict->g,rdict->ppow[1]);
  
  if(DEBUG) {    
    printf("Reduced polynomial g modulo %lu, the result is\n",p);
     poly_print(rdict->gmodp);
  }

  if(DEBUG) {    
     printf("Composing polynomials ...\n");
  }
  // You must get the universe right.
  //compose_reduce(fdg_gdf,rdict->fdg_gdf,l0,rdict->ppow,rdict->maxpow,rdict->ppow[rdict->maxpow]);
  poly_reduce_mod_mpz(fdg_gdf,rdict->fdg_gdf,rdict->ppow[rdict->maxpow]);

  if(DEBUG) {        
     printf("Composition/reduction of fdg_gdf, result is\n");
     poly_print(fdg_gdf);
     printf("and the original polynomial fdg_gdf is\n");
     poly_print(rdict->fdg_gdf);
  }

  //compose_reduce(ff,ff,l0,rdict->ppow,rdict->maxpow,rdict->ppow[rdict->maxpow]);
  poly_reduce_mod_mpz(ff,ff,rdict->ppow[rdict->maxpow]);

  if(DEBUG) {    
     printf("Composition/reduction of ff, result is\n");
     poly_print(ff);
  }

  //compose_reduce(gg,gg,l0,rdict->ppow,rdict->maxpow,rdict->ppow[rdict->maxpow]);  
  poly_reduce_mod_mpz(gg,gg,rdict->ppow[rdict->maxpow]);
  rdict->l0 = l0;  

  if(DEBUG) {    
     printf("Composition/reduction of gg, result is\n");
     poly_print(gg);
  }

  poly_t id;
  poly_alloc(id,1);
  poly_set_identity(id);  
  
  if(DEBUG) {    
     printf("Entering rootsieving recursion ...\n");     
  }

  hits = rotation_inner(rdict,p,ff,gg,fdg_gdf,u0,v0,l0,ld,m,0,scale,id); 

  if(DEBUG) {    
    printf("Finished rootsieving modulo %lu and its powers...\n",p);    
  }

  poly_free(ff);
  poly_free(gg);
  poly_free(fdg_gdf);
  poly_free(id);

  for( i=0 ; i<=rdict->maxpow; ++i)
    mpz_clear(rdict->ppow[i]);  
  free(rdict->ppow); 

  poly_free(rdict->gmodp);

  return hits;  
}


long rotation_inner(rootsieve_dictionary rdict, unsigned long p,poly_t ff,poly_t gg,poly_t fdg_gdf, long u0, long v0, long l0,int ld,int m,long twist_v1,double scale,poly_t dphi) {

  if(DUMP) {
    printf("\n\n\nInner call:\n");
    printf("p=%ld,u0=%ld,v0=%ld,l0=%ld,ld=%d,m=%d,twist_v1=%ld,scale=%f\n",p,u0,v0,l0,ld,m,twist_v1,scale);    
    printf("ff=\n");
    poly_print(ff);
    printf("gg=\n");
    poly_print(gg);
    printf("fdg_gdf=\n");
    poly_print(fdg_gdf);
    printf("dphi=\n");
    poly_print(dphi);    
    print_dictionary(rdict);
  }

  /*if(COMP) {
    printf("RI:m=%d,ld=%d\n",m,ld);
    }*/

  if(TRACE) {
    //printf("\n%lu %d %d %2.8f %ld %ld %ld %ld ",p,m,ld,scale,u0,v0,l0,twist_v1);
    printf("\n--------------------------------------------------------\n\nm=%d, ld=%d\n ",m,ld);
  }
  
  // Some necessary assertions
  ASSERT_ALWAYS(m>=ld);

  if(DEBUG) {
    printf("(Rotation Inner Start) m=%d lambda=%d\n",m,ld);
  }
  
  // Use your precomputations.
  // Update your precomputations

  //if len(hist) > 6:
  //      nhist=copy(hist)
  //      nhist.append([dphi,u0,v0,dphi(0),ld,m])
  //      print "Found tree of height %d" % len(nhist)
  //      for hi in nhist:
  //          print hi
  //

  if (scale < rdict->cutoff) {    
    printf("Cutoff aborted\n");    
    return 0;
  }

  // Hit counter
  long thits,hits;
  thits = 0;
  hits = 0;

  // Abbreviations
  long igl;
  long pm,pm1,pmax;
    
  // As igl is a residue mod p, the following value is bogus 
  // and it just means "undefined". Added to avoid a warning.
  igl = -1; 

  pm    = mpz_get_si(rdict->ppow[m]);
  pm1   = mpz_get_si(rdict->ppow[m+1]);
  pmax  = mpz_get_si(rdict->ppow[rdict->maxpow]);  
 
  if(DEBUG) {
    printf("Charged powers of %lu : %ld, %ld, %ld\n",p,pm,pm1,pmax);
  }
  
  double scale1,scale2,scale3;
  int look_many_roots;
  poly_t minus_f_over_g_modp;
  
  // gg is never scaled down, and it is always evaluated at points where
  // the value is not zero mod p. So the valuation of g mod p never
  // vanishes.
  
  // This is always useful.  

  scale1 = scale/(p-1);
  scale2 = scale/p;
  scale3 = scale1-scale2;  
  
  if(DEBUG) {
    printf("Scales 1,2,3 are %f, %f, %f.\n",scale1,scale2,scale3);
  }
  
  poly_alloc(minus_f_over_g_modp, ff->deg);
  poly_reduce_mod_mpz(minus_f_over_g_modp, ff,rdict->ppow[1]);  

  if(DEBUG) {
    printf("Polynomial -f/g mod %lu :\n",p);
    poly_print(minus_f_over_g_modp);
  }

  if (m > 0) {
    if(DEBUG) {
      printf("Entering positive m section...\n");
    }
    l0 = rdict->l0;    
    igl = rdict->igl;
    if(DEBUG || TRACE) {
      printf("\nl0=%ld, 1/g(l) mod p=%ld\n",l0,igl);
      
    }    
    ASSERT_ALWAYS(0<igl && igl<(long)p);
    mpz_t aux;
    mpz_init_set_si(aux,-igl);
    poly_mul_mpz(minus_f_over_g_modp,minus_f_over_g_modp,aux);
    mpz_clear(aux);
    poly_reduce_mod_mpz(minus_f_over_g_modp,minus_f_over_g_modp,rdict->ppow[1]);    
    if(DEBUG||TRACE) {
      printf("\n-f/g mod %lu= ",p);
      poly_print(minus_f_over_g_modp);
      printf("\n");
    }
  }
  

  // Deciding whether we're in the general case or not has some
  // subtleties, unfortunately.
  if(DEBUG) {
    printf("Computing the \"Possibly multiple root\" condition.  \n");
    printf("The polynomial -f/g mod %lu is%s constant\n",p,(poly_is_constant(minus_f_over_g_modp))?"":" not");
    if(m==0) {
      printf("The polynomial g mod %lu\n",p);
      poly_print(rdict->gmodp);
      printf("is%s constant\n",(poly_is_constant(rdict->gmodp))?"":" not");
    }
  }
  
  // The following must hold, otherwise g is not irreducible
  //ASSERT_ALWAYS(rdict->gmodp->deg >=0);
  if(rdict->gmodp->deg == -1) {

    if(DEBUG) {
      printf("Error : g can be factored by %lu, hence it is not irreducible. Halting.\n",p);
    }

    ASSERT_ALWAYS(thits==0);
    return thits;
  }

  look_many_roots = (!poly_is_constant(minus_f_over_g_modp)) ||  (m == 0 &&  ! poly_is_constant(rdict->gmodp));

  if(TRACE) {
    printf("\nlook_many_roots=%d\n ",look_many_roots);
  }

  if(DEBUG) {
    printf("We%s expect multiple roots.\n",(look_many_roots)?"":" don't");
  }
  //nhist = copy(hist);
    
  //if (!look_many_roots)
  //  nhist.append([dphi,u0,v0,l0+p*dphi(0),ld,m]);
  
  long dv0,u1,v1;
  long du,du0,new_twist_v1;
  poly_t nfdg_gdf,nff,ngg,ndphi,twist_f,tmp_poly;
  mpz_t rhs;
  
  poly_alloc(nfdg_gdf,fdg_gdf->deg);
  poly_alloc(nff,ff->deg);
  poly_alloc(ngg,gg->deg);
  poly_alloc(ndphi,1);
  poly_alloc(twist_f,0);
  poly_alloc(tmp_poly,0);
  mpz_init(rhs);

  if (look_many_roots) {
    // Then it's the typical case, as encountered for instance at the
    // beginning of the root sieve.
    // Each possible l value must be tried, and will give rise to
    // potentially intersecting lattices. 
    
    if(DEBUG) {
      printf("Entering the multiple root section...\n");
    }
    
    long l;

    
    for ( l=0 ; l<(long)p ; l++) {
      
      //nhist=copy(hist);
      //nhist.append([dphi,u0,v0,l0+p*dphi(l),ld,m]);
      long minus_f_over_gmodp_l,gvl;
      mpz_t mp_tmp, mpl;
      mpz_init(mp_tmp);
      mpz_init_set_si(mpl,l);
      
      poly_eval_mod_mpz(mp_tmp, minus_f_over_g_modp,mpl,rdict->ppow[1]);
      minus_f_over_gmodp_l = mpz_get_si(mp_tmp);     

      if(DEBUG) {
	printf("Evaluated -f(%ld)/g(%ld) mod %lu to ",l,l,p);
	mpz_out_str(stdout,10,mp_tmp);
	printf("\n");
      }
  
      if (m==0) {

	// evaluation 
	//gvl = rdict->gmodp(l);
	
	if(DEBUG) {
	  printf("Entered m==0 section ...\n");
	}
	
	poly_eval_mod_mpz(mp_tmp,rdict->gmodp,mpl,rdict->ppow[1]);
	gvl = mpz_get_si(mp_tmp);

	if(DEBUG) {
	  printf("Evaluated g(%ld) mod %lu to %ld ...\n",l,p,gvl);
	}

	if (gvl == 0) {
	  if(DEBUG) {
	    printf("Singular g(%ld) mod %lu aborted...\n",l,p);
	  }
	  mpz_clear(mp_tmp);
	  mpz_clear(mpl);	  
	  continue;
	}
		
	igl = (long)modular_inverse((unsigned long)gvl,p);
	
	rdict->l0 = l;
	rdict->igl = igl;
	twist_v1 = -rdict->l0;

	if(DEBUG) {
	  printf(" 1/g(%ld) mod %lu = %ld, twist_v1 = %ld ...\n",l,p,igl,twist_v1);	  
	}
	
	ASSERT_ALWAYS(0<igl && igl<(long)p);

	// The next line divides f(l) by g(l) mod p
	minus_f_over_gmodp_l = (long)modular_product((unsigned long)((minus_f_over_gmodp_l==0)?(minus_f_over_gmodp_l):((long)p-minus_f_over_gmodp_l)),(unsigned long)igl,p);
	
	if(DEBUG) {
	  printf("-f(%ld)/g(%ld) mod %lu = %ld ...\n",l,l,p,minus_f_over_gmodp_l);
	}

	ASSERT_ALWAYS(0<=minus_f_over_gmodp_l && minus_f_over_gmodp_l<(long)p);
      }

      mpz_clear(mp_tmp);
      
      dv0 = minus_f_over_gmodp_l;
      u1  = u0;
      v1  = v0 + pm*dv0;

      if(DEBUG) {
	printf("The lattice : dv0 = %ld, u1 = %ld, v1 = %ld...\n",dv0,u1,v1);
      }

      hits = add_lattice(rdict,u1,v1, pm, pm1, twist_v1, rdict->U0,rdict->U1,rdict->V0,rdict->V1,-scale1,p);      

      if(COMP) {
	printf("1:m=%d ld=%d hits>0=%d u0=%ld v0=%ld us=%ld vs=%ld skew=%ld contrib=%f\n",m,ld,(hits>0)?1:0,u1,v1,pm,pm1,twist_v1,-scale1);
	//printf("1:%d \n",(hits>0)?1:0);
      }

      if (hits == 0) {
	if(DEBUG) {
	  printf("No hits, continuing to next iteration ...\n");
	}		
	mpz_clear(mpl);	
	continue;	
      }

      // print "main, level %d : %d hits" % (m, hits)
      thits += hits;
      // We're optimistic here ; we're counting on the fact that we
      // expect the derivative not to cancel.
            
      // However, there is the possibility that we reach the
      // cancellation point. This only once per (p,p) square,
      // meaning that u,v are constrained.

      // Try to guess the right value for du.
      // Rp = Integers(ppow[m-ld+1])
      // RpP=Rp['w']
      // rhs = Z(Z((RpP(fdg_gdf) - u0*RpP(gg)^2)(l))/ppow[m-ld])
      // Slightly faster when not doing reductions while evaluating.
      // Most probably due to the cost of computing remainders all
      // the way to the final value, which ends up being more
      // expensive than it should.

      // TODO: There has to be a way to remove valuations very early
      // on here. In fact, my guess is that it could very probably
      // be done every time m-ld increases...
      
      
      //   This can be simplified.
      //      rhs = Z(K((fdg_gdf(l)-u0*gg(l)^2)/ppow[m-ld]));
     
      mpz_t aux;      
      mpz_init(aux);
      poly_eval(rhs,fdg_gdf,mpl);
      poly_eval(aux,gg,mpl);
      mpz_mul(aux,aux,aux);
      mpz_mul_si(aux,aux,u0);      
      mpz_sub(rhs,rhs,aux);

      ASSERT_ALWAYS(m>=ld);

      mpz_fdiv_q(rhs,rhs,rdict->ppow[m-ld]);
      mpz_mod_ui(rhs,rhs,p);
      mpz_clear(aux);
      mpz_clear(mpl);
      
      
      if ((ld > 0) && (mpz_sgn(rhs)!= 0)) {

	if(TRACE>1) {
	  gmp_printf("\nrhs=%Zd\n",rhs);
	}

	if(DEBUG) {
	  printf("Root is simple, aborting the anticontribution section...\n");
	}
	continue;	
      }

      if(TRACE>1) {
	gmp_printf("\nrhs=%Zd\n",rhs);
      }
	
      long nspots;
      nspots = -1; // This bogus initial value is used for testing that assignment is actually done.

      compose_reduce(nfdg_gdf,fdg_gdf,l,rdict->ppow,rdict->maxpow,rdict->ppow[rdict->maxpow]);      
      compose_reduce(nff,ff,l,rdict->ppow,rdict->maxpow,rdict->ppow[rdict->maxpow]);
      compose_reduce(ngg,gg,l,rdict->ppow,rdict->maxpow,rdict->ppow[rdict->maxpow]);  

      if(DEBUG || TRACE>2) {
	  printf("\nComposed reduced polynomials\n nfdg_gdf=");
	  poly_print(nfdg_gdf);
	  printf("\n\nnff=");
	  poly_print(nff);
	  printf("\n\n ngg=");
	  poly_print(ngg);
	  printf("\n");
      }

      
      if (m == 0) {
	// XXX This is fairly bizarre. Results are correct with
	// this, but I don't see at the moment why we don't have
	// psi instead of x here.

	if(DEBUG) {
	  printf("Fairly bizarre (m==0)...\n");
	}
	
	poly_set_identity(ndphi);

	if(DEBUG) {
	  printf("Polynomial ndphi is ...\n");
	  poly_print(ndphi);
	}		
      }
      else {        
	compose_psi(ndphi,dphi,l,rdict->ppow,rdict->maxpow);
	if(DEBUG) {
	  printf("Polynomial ndphi is ...\n");
	  poly_print(ndphi);
	}	
      }
      
      poly_mul_ui(tmp_poly,ngg,(unsigned long)dv0);       
      poly_add(tmp_poly,nff,tmp_poly);      
      poly_div_ui(nff,tmp_poly,p);            
      if(DEBUG) {
	  printf("Polynomial nff in the lattice is ...\n");
	  poly_print(nff);
      }

      poly_mul(twist_f,ndphi,ngg);				

      if(DEBUG) {
	  printf("Polynomial twist_f is ...\n");
	  poly_print(twist_f);
      }

      if (ld > 0) {
	if(DEBUG) {
	  printf("Entered lambda>0 spots section...\n");
	}
	// Then all u's will do. Note that this case happens
	// only rather deep inside the recursion. So indeed
	// we're going to re-hit the places which have been
	// hit before, but that's not too worrying.

	// This test was redundant, left here for reference.
	ASSERT_ALWAYS(mpz_sgn(rhs) == 0);
	nspots=p;		
      }
      else {	
	if(DEBUG) {
	  printf("Entered lambda=0 spots section...\n");
	  printf("Value of g(l) mod p is %ld ...\n",igl);	
	}
      	
	ASSERT_ALWAYS(igl>0 && igl<(long)p);

	du0 = (igl*igl*mpz_get_si(rhs)) % p;
	u1 += du0 * pm;
	v1 += du0 * twist_v1;
	
	if(DEBUG) {
	  printf("du0=%ld, u1=%ld, v1=%ld\n",du0,u1,v1);
	}

	poly_mul_ui(tmp_poly,twist_f,(unsigned long)du0);
	if(DEBUG) {
	  printf("du0*twistf=");
	  poly_print(tmp_poly);
	  printf("\n");
	}
	poly_add(nff,nff,tmp_poly);
	
	if(DEBUG) {
	  printf("Polynomial nff twisted when ld=0 is ...\n");
	  poly_print(nff);
	}	
	
	if(DEBUG) {
	  printf("Constructed new f(l+p*x) polynomial :\n");
	  poly_print(nff);
	}
	nspots=1;
      }	          

      if(DEBUG) {
	printf("Number of spots to visit : %ld\n",nspots);	
      }
      
      ASSERT_ALWAYS(nspots>0);      
      new_twist_v1 = p * twist_v1;

      if(DEBUG) {
	printf("The new twist is new_twist_v1=%ld \n",new_twist_v1);
      }

      if(DEBUG) {
	printf("Entering anticontribution loop...\n");
      }
      
      
      for ( du=0 ; du<nspots; du++) {

	if(DEBUG) {
	  printf("Distributing contribution\n");
	}
	hits = add_lattice(rdict, u1,v1, pm1, pm1, 0, rdict->U0, rdict->U1, rdict->V0, rdict->V1, scale3,p);
	
	if(COMP) {
	  printf("2:m=%d ld=%d hits>0=%d u0=%ld v0=%ld us=%ld vs=%ld skew=%ld contrib=%f\n",m,ld,(hits>0)?1:0,u1,v1,pm1,pm1,0L,scale3);
	  //printf("2:%d \n",(hits>0)?1:0);
	}

	if(DEBUG) {
	  printf("Had %ld hits...\n",hits);
	}

	// print "main/excep, level %d : %d hits" % (m, hits)
	thits += hits;
	// scale3 is scale/(p-1)-scale/p ; hence it brings
	// back the cells to the contribution -scale/p, which is the
	// contribution from a _single_ zero at this location, not
	// liftable because of ramification. Later recursive calls
	// will investigate the possibility that despite the multiple
	// root mod p, we still get roots at higher precision.

	

	if (hits >0) {	      
	  
	  if (1) //m<rdict->m_max && m<rdict->maxpow) 
	    thits += rotation_inner(rdict,p,nff,ngg,nfdg_gdf,u1,v1,l0,ld+1,m+1,new_twist_v1,scale2,ndphi);	      	    
	  else {  	    
	    if(DEBUG) {
	      if(m == rdict->m_max)
		printf("Rectangle aborted ");
	      
	      if(m == rdict->maxpow)
		printf("Power aborted");
	      
	      printf("\n");
	    }	    	    
	  }	    
	}
	else {
	  if(DEBUG) {
	    printf("No recursion because no hits...\n");
	  }	  
	}

	poly_add(nff,nff,twist_f);		
	u1 += pm;
	v1 += twist_v1;
	
	if(DEBUG) {
	  printf("Updated lattice...\n");
	  printf("u1=%ld, v1=%ld, and f(l+p*x) is\n",u1,v1);
	  poly_print(nff);
	}
      }
	
      
    
      if(DEBUG) {
	printf("Completed multiple root section: exiting...\n");
      }    
    }
  }
  else {
    // f mod p is a constant polynomial.  This means that our
    // expression has an extra root only for specific u,v pairs.
    // Basically, we have a linear term in u and v which must be
    // cancelled.

    if(DEBUG) {
      printf("Entering simple root section...\n");
      poly_t auxpoly;
      poly_alloc(auxpoly,0);
      poly_reduce_mod_mpz(auxpoly,ff,rdict->ppow[1]);
      ASSERT_ALWAYS(poly_is_constant(auxpoly));
      //printf("f mod %lu is supposed to be constant, and the test gives %s \n",p,(poly_is_constant(auxpoly))?"constant":"not constant");
      poly_free(auxpoly);
    }
 
    mpz_t mp_tmp;
    mpz_init(mp_tmp);

    if (m==0) {
      
      // This case is rare in practice, but it appeared in the example in the
      // example given in the file kleinjung.sage, for p=37. 
      //
      // We need to initialize 1/g(l) mod p, knowing that g is constant.     
      // IMPORTANT : we know that this inverse igl will be valid for all values
      // of l, that's why we can compute it here, instead of in the main loop.

      long gvl,igl;             
	
      if(DEBUG) {
	printf("Entered m==0 section ...\n");
      }

      poly_getcoeff(mp_tmp,0,rdict->gmodp);      
      gvl = mpz_get_si(mp_tmp);
      ASSERT_ALWAYS(0<=gvl && gvl<(long)p);

      if(DEBUG) {
	printf("Constant polynomial g mod %lu evaluates to %ld ...\n",p,gvl);
      }
      
      ASSERT_ALWAYS(gvl!=0); // Otherwise g is reducible
		
      igl = (long)modular_inverse((unsigned long)gvl,p);
	
      rdict->igl = igl;

      if(DEBUG) {
	printf(" 1/g mod %lu = %ld, ...\n",p,igl);	  
      }
	
      ASSERT_ALWAYS(0<igl && igl<(long)p);
      
    }
        
    poly_getcoeff(mp_tmp,0,minus_f_over_g_modp);

    if(DEBUG) {
      printf("The constant function -f/g mod %lu equals ",p);
      mpz_out_str(stdout,10,mp_tmp);
      printf("\n");
    }

    dv0 = mpz_get_si(mp_tmp);
    mpz_clear(mp_tmp);
    u1  = u0;
    v1  = v0 + pm * dv0;
    
    if(DEBUG) {
      printf("New lattice : u1=%ld, v1=%ld\n Adding contribution...",u1,v1);
    }
    
    hits=add_lattice(rdict, u1, v1, pm, pm1, twist_v1, rdict->U0, rdict->U1, rdict->V0, rdict->V1, -scale,p);
    
    if(COMP) {      
      printf("3:m=%d ld=%d hits>0=%d u0=%ld v0=%ld us=%ld vs=%ld skew=%ld contrib=%f\n",m,ld,(hits>0)?1:0,u1,v1,pm,pm1,twist_v1,-scale);
      //printf("3:%d \n",(hits>0)?1:0);
    }

    if(DEBUG) {
      printf("Counted %ld hits in the rectangle.\n",hits);
    }

    if (hits == 0) {
      
      if(DEBUG) {
	printf("No hits, returning total hits and terminating (Rotation Inner)...\n");
      }
                 
      poly_free(minus_f_over_g_modp);
      poly_free(nfdg_gdf);
      poly_free(nff);
      poly_free(ngg);
      poly_free(ndphi);
      poly_free(twist_f);
      poly_free(tmp_poly);
      mpz_clear(rhs);

      return thits;

    }
    // print "secondary, level %d : %d hits" % (m, hits)

    thits += hits;    
            
    poly_mul_ui(tmp_poly,gg,(unsigned long)dv0);
    poly_add(tmp_poly,ff,tmp_poly);
    poly_div_ui(nff,tmp_poly,p);    
    poly_mul(twist_f,dphi,gg);				
    new_twist_v1 = p * twist_v1;
    
    if(DEBUG) {
      printf("New twisted polynomial :\n");
      poly_print(twist_f);
      printf("New lattice twist : new_twist_v1=%ld\n",new_twist_v1);
    }
    
    if (1) { //m<rdict->m_max && m<rdict->maxpow) {
      for ( du=0 ; du<(long)p ; du++) {      

	if(DEBUG) {
	  printf("Recalling (Rotation Inner)...\n");
	}
	hits = rotation_inner(rdict,p,nff,gg,fdg_gdf,u1,v1,l0,ld,m+1,new_twist_v1,scale,dphi);
	
	// print "secondary/excep, level %d : %d hits" % (m, hits)
	if(DEBUG) {
	  printf("Rotation Inner returned %ld hits.\n",hits);
	}
	thits += hits;
	u1 += pm;
	v1 += twist_v1;      
	poly_add(nff,nff,twist_f);
	
	if(DEBUG) {
	  printf("New lattice : u1=%ld, v1=%ld, and the twisted polynomial is\n",u1,v1);
	  poly_print(nff);
	}
      }
    }
    else {
	
      if(DEBUG) {	  

	if(m == rdict->m_max)
	  printf("Rectangle aborted ");
	  
	if(m == rdict->maxpow)
	  printf("Power aborted");
	  
	printf("\n");
      }	
    }
     
    if(DEBUG) {
      printf("Completed simple root section: exiting...\n");
    } 
    
  }
  poly_free(minus_f_over_g_modp);
  poly_free(nfdg_gdf);
  poly_free(nff);
  poly_free(ngg);
  poly_free(ndphi);
  poly_free(twist_f);
  poly_free(tmp_poly);
  mpz_clear(rhs);

  return thits;
}

/* Returns 1 when the lattice generated by {(us,skew),(0,vs)} has some point in [U0,U1]x[V0,V1], and it adds the lattice to the 
   list rdict->lattice. Otherwise it just returns 0.
*/

long add_lattice(rootsieve_dictionary rdict,long u0, long v0, long us, long vs, long skew, long U0, long U1, long V0, long V1, double contribution, unsigned long p) {
  long u,v;
  long q;
  long hits;

  hits=0;

  // We ensure that the lattice basis is positively oriented.
  if(vs<0)
    vs = -vs;
  if(us<0) {
    us = -us;
    skew = -skew;
  }
  while(skew<0)
    skew += vs;  
  
  ASSERT_ALWAYS(vs>=0 && skew >= 0 && us>=0);  

  q = (long)floor(((double)(u0-U0))/((double)us));
  u = u0-q*us;
  v = v0-q*skew;

  q = (long)floor(((double)(v-V0))/((double)vs));
  v = v - q*vs;

  ASSERT_ALWAYS(v>=V0);
  ASSERT_ALWAYS(u>=U0);

  // If the lattice has some point in the rectangle, we add the lattice to the list.
	// Here 1 hit == 1 lattice
  if (!empty(u,v,us,vs,skew,U0,U1,V0,V1)) {
    if(DEBUG>9)
      printf("Initial point u=%ld,v=%ld, belongs to the rectangle. \n",u,v);
    rootsieve_lattice rl;
    rl.u0 = u;
    rl.v0 = v;
    rl.us = us;
    rl.vs = vs;
    rl.skew = skew;
    rl.contribution = contribution;
    rl.p = p;
    if(rdict->lattice_alloc == rdict->L) {
      rdict->lattice_alloc += REALLOC_INCREMENT;
      rdict->lattices = realloc(rdict->lattices,rdict->lattice_alloc*sizeof(rootsieve_lattice));
    }
    rdict->lattices[rdict->L] = rl;
    rdict->L++;
    return 1;
  }
  else {
    if(DEBUG>9)
      printf("Initial point u=%ld,v=%ld, does NOT belong to the rectangle. \n",u,v);
    return 0; 
  }
     
}

int empty(long u0, long v0, long us, long vs, long skew, long U0, long U1, long V0, long V1) {

  ASSERT_ALWAYS(u0>=U0 && v0>=V0);
  long u,v;
  long q;

  if (u0<=U1) {
    if (v0<=V1)
      return 0;
    else {
      if (skew==0)
	return 1;      
      u = u0;
      v = v0-vs;
      while(u<U1) {
	u+=us;
	v+=skew;
	q = (long)floor(((double)(v-V0))/(double)vs);
	v = v-q*vs;
	if (v<=V1 && u<=U1){
	  //fprintf(stderr,"u0=%d v0=%d u=%d v=%d i=%d\n",u0,v0,u,v,i);
	  return 0;
	}	
      }
      return 1;
    }
  }
  else 
    return 1;
  
}

unsigned long modular_product( unsigned long a, unsigned long b, unsigned long N ) {

  residue_t res,opa,opb;
  unsigned long result;
  modulus_t modu;
	
  modul_initmod_ul(modu,N);
  modul_init(res,modu);
  modul_init(opa,modu);
  modul_init(opb,modu);
  modul_set_ul_reduced(opa,a,modu);
  modul_set_ul_reduced(opb,b,modu);	
  modul_mul(res,opa,opb,modu);
  result = modul_get_ul(res,modu);
  modul_clear(opa,modu);
  modul_clear(opb,modu);
  modul_clear(res,modu);
  modul_clearmod(modu);

  return result;
}


unsigned long modular_inverse( unsigned long r,  unsigned long N ) {

  //ASSERT_ALWAYS(gcd(r,N)==1);
  residue_t res,inv;
  modulus_t modu;

  unsigned long result;

  modul_initmod_ul(modu,N);
  modul_init(res,modu);
  modul_init(inv,modu);
  modul_set_ul_reduced(res,r,modu);
  modul_inv(inv,res,modu);
  result = modul_get_ul(inv,modu);
  modul_clear(inv,modu);
  modul_clear(res,modu);
  modul_clearmod(modu);

  return result;
}

// Polynomials
// Ensures that f(X)=X
void poly_set_identity(poly_t f) {
  mpz_set_ui(f->coeff[0],0);
  mpz_set_ui(f->coeff[1],1);
  cleandeg(f,1);
}

// psi(x) = l + p*x
// f <- g(psi) 
// f must be allocated, f and g can be the same. 
// TODO : Ensure that no constraints on the degree of f are necessary.
void compose_psi(poly_t f, const poly_t g, unsigned long l, mpz_t * ppow, int maxpow) {
  
  // Suspicious case
  if(g->deg <=0) {
    if(DEBUG) {
      printf("(Compose psi) Warning : received constant polynomial of degree %d\n",g->deg);
    }
    poly_copy(f,g);
    return;
  }

  ASSERT_ALWAYS(g->deg <=maxpow); // Otherwise we will get a segfault.
  ASSERT_ALWAYS(g->alloc >= g->deg +1); // Idem  

  unsigned int d;  
  d = (unsigned int)(g->deg);

  unsigned long i,k;
  mpz_t tmpsum,tmpterm,tmplpow;
  poly_t aux;
  poly_alloc(aux,d);


  mpz_init(tmpsum);
  mpz_init(tmpterm);
  mpz_init(tmplpow);    

  //printf("Degree of polynomial to be reduced is %d ...\n",d);
  //printf("(Compose Psi) Entering compose loop ...\n");

  for( k=0 ; k<=d ; k++) {

    mpz_set_ui(tmplpow,1);
    mpz_set_ui(tmpsum,0);
       
    for( i=k ; i<=d ; i++) {
      // printf("Computing term %lu of the sum...\n",i);
      // Compute the thing
      mpz_bin_uiui(tmpterm,i,k);
      //printf("Before coeff...\n");
      mpz_mul(tmpterm,tmpterm,g->coeff[i]);
      //printf("Before ppow...\n");      
      mpz_mul(tmpterm,tmpterm,ppow[k]);
      //mpz_mul(tmpterm,tmpterm,rdict->ppow[k]);
      //printf("After both\n");
      mpz_mul(tmpterm,tmpterm,tmplpow);      
      // Update pows
      mpz_mul_ui(tmplpow,tmplpow,l);         
      // Add to the sum
      mpz_add(tmpsum,tmpsum,tmpterm);      
    }
    //printf("Computed the sum...\n");
    //printf("(Compose Psi) Coeff %lu of aux is \n",k);    
    //mpz_out_str(stdout,10,tmpsum);
    //printf("\n");
    poly_setcoeff(aux,k,tmpsum);
  }

  //printf("(Compose Psi) The aux polynomial, before copying into f :\n");
  //poly_print(aux);
  poly_copy(f,aux);

  //printf("(Compose Psi) The output's f polynomial :\n");
  //poly_print(f);

  mpz_clear(tmpsum);
  mpz_clear(tmpterm);
  mpz_clear(tmplpow);  
  poly_free(aux);

}

/* We assume (rdict->delta) > 0. No problem when f and g are the same. */
void compose_reduce(poly_t f, const poly_t g,unsigned long l,mpz_t * ppow, int maxpow, mpz_t m) {
  poly_t aux;
  poly_alloc(aux,g->deg);
  compose_psi(aux,g,l,ppow,maxpow);
  poly_reduce_mod_mpz(f,aux,m);
  //printf("(Compose reduce) The original and reduced polynomials are\n");
  //poly_print(g);
  //printf("(Compose reduce) and\n");
  //poly_print(f);
  poly_free(aux);
}

/* Prints a   */
void print_lattices(lattice_list l) {
  long i;
  printf("print_lattices:\n");
  for(i=0; i<l.length; i++) 
    printf("u0=%ld v0=%ld us=%ld vs=%ld skew=%ld contrib=%f p=%lu\n",l.lattices[i].u0,l.lattices[i].v0,l.lattices[i].us,l.lattices[i].vs,l.lattices[i].skew,l.lattices[i].contribution,l.lattices[i].p);  
}

void print_dictionary(rootsieve_dictionary rdict) {

  printf("\n\nDictionary :\n");
  printf("m_max=%d,maxpow=%d,cutoff=%f,l0=%ld,igl=%ld\n\n",rdict->m_max,rdict->maxpow,rdict->cutoff,rdict->l0,rdict->igl);
  printf("f:\n");
  poly_print(rdict->f);
  printf("g:\n");
  poly_print(rdict->g);
  printf("fdg_gdf:\n");
  poly_print(rdict->fdg_gdf);
  printf("gmodp:\n");
  poly_print(rdict->gmodp);  
  
}

