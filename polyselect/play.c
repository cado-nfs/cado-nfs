#include "cado.h"
#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdint.h>
#include <stdlib.h>
#include <time.h>

#define UNLIKELY(x)     __builtin_expect(x, 0)

clock_t delta;

int
naive_search (double **f, int d, int l, double eps)
{
  int *mu, i;
  double *s, fr;
  int found = 0;

  mu = (int*) malloc (l * sizeof (int));
  s = (double*) malloc ((l + 1) * sizeof (double));
  /* s[i] contains f0 + f[0][mu[0]] + ... + f[i-1][mu[i-1]] */

  /* initializes mu[] to (0,...,0) */
  for (i = 0; i < l; i++)
    mu[i] = 0;
  s[0] = 0;
  for (i = 1; i <= l; i++)
    s[i] = s[i-1] + f[i-1][mu[i-1]];

  delta = -clock();

  while (1)
    {
      /* check current sum, the following is a trick to avoid a call to
         the round() function, which is slow */
      fr = s[l] + 6755399441055744.0;
      fr = fr - 6755399441055744.0; /* fr = round(s[l]) */
      fr = fabs (fr - s[l]);
      if (UNLIKELY(fr <= eps)) /* Prob ~ 4e-7 on RSA155 with l=7, degree 5,
                                  M=5e24, pb=256, even less for smaller M */
	{
            found++;
        }
      
      /* go to next combination */
      for (i = l - 1; i >= 0 && mu[i] == d - 1; i--);
      /* now either i < 0 and we are done, or mu[i] < d-1 */
      if (UNLIKELY(i < 0))
	break;
      mu[i] ++;
      s[i+1] = s[i] + f[i][mu[i]];
      while (++i < l)
	{
	  mu[i] = 0;
	  s[i+1] = s[i] + f[i][mu[i]];
	}
    }

  delta += clock();

  free (mu);
  free (s);

  return found;
}

int search2(double **f, int d, int l, double eps)
{
    uint64_t ** g;
    double two_64 = pow(2,64);
    uint64_t lim = 2 * eps * two_64;
    int *mu, i, j;
    uint64_t *s;
    int found = 0;

    mu = (int*) malloc (l * sizeof (int));
    s = (uint64_t*) malloc ((l + 1) * sizeof (uint64_t));
    g = malloc(l * sizeof(uint64_t*));
    for(i = 0 ; i < l ; i++) {
        g[i] = malloc(d * sizeof(uint64_t));
        for(j = 0 ; j < d ; j++) {
            g[i][j] = f[i][j] * two_64;
        }
    }
    /* Offset by epsilon. */
    for(j = 0 ; j < d ; j++) {
        g[0][j] = (f[0][j]+eps) * two_64;
    }


    /* initializes mu[] to (0,...,0) */
    for (i = 0; i < l; i++)
        mu[i] = 0;
    s[0] = 0;
    for (i = 1; i <= l; i++)
        s[i] = s[i-1] + g[i-1][mu[i-1]];

    delta = -clock();

    while (1)
    {
        /* check current sum */
        /* Prob ~ 4e-7 on RSA155 with l=7, degree 5, M=5e24, pb=256,
         * even less for smaller M */
        found += s[l] <= lim;

        /* go to next combination */
        for (i = l - 1; i >= 0 && mu[i] == d - 1; i--);
        /* now either i < 0 and we are done, or mu[i] < d-1 */
        if (UNLIKELY(i < 0))
            break;
        mu[i] ++;
        s[i+1] = s[i] + g[i][mu[i]];
        while (++i < l)
        {
            mu[i] = 0;
            s[i+1] = s[i] + g[i][mu[i]];
        }
    }

    delta += clock();



    free(mu);
    free(s);
    for(i = 0 ; i < l ; i++) {
        free(g[i]);
    }
    free(g);

    return found;
}

typedef int (*sortfunc_t) (const void *a, const void * b);

int cmp(uint64_t * a, uint64_t * b)
{
    if (*a < *b) return -1;
    if (*b < *a) return 1;
    return 0;
}

void save_all_sums(uint64_t * dst, uint64_t ** g, int d, int l)
{
    int *mu, i;
    uint64_t *s;
    mu = (int*) malloc (l * sizeof (int));
    s = (uint64_t*) malloc ((l + 1) * sizeof (uint64_t));
    uint64_t * d0 = dst;

    for (i = 0; i < l; i++)
        mu[i] = 0;
    s[0] = 0;
    for (i = 1; i <= l; i++)
        s[i] = s[i-1] + g[i-1][mu[i-1]];

    while (1)
    {
        *dst++ = s[l];

        /* go to next combination */
        for (i = l - 1; i >= 0 && mu[i] == d - 1; i--);
        /* now either i < 0 and we are done, or mu[i] < d-1 */
        if (UNLIKELY(i < 0))
            break;
        mu[i] ++;
        s[i+1] = s[i] + g[i][mu[i]];
        while (++i < l)
        {
            mu[i] = 0;
            s[i+1] = s[i] + g[i][mu[i]];
        }
    }

    qsort(d0, dst-d0, sizeof(uint64_t), (sortfunc_t) cmp);

    free(s);
    free(mu);
}
int search3(double **f, int d, int l, double eps)
{
    uint64_t ** g;
    double two_64 = pow(2,64);
    uint64_t lim = 2 * eps * two_64;
    int found = 0;
    int i,j;
    unsigned int dll, dlh;

    g = malloc(l * sizeof(uint64_t*));
    for(i = 0 ; i < l ; i++) {
        g[i] = malloc(d * sizeof(uint64_t));
        for(j = 0 ; j < d ; j++) {
            g[i][j] = f[i][j] * two_64;
        }
    }
    /* Offset by epsilon. */
    for(j = 0 ; j < d ; j++) {
        g[0][j] = (f[0][j]+eps) * two_64;
    }

    int lh = l - l/2;
    int ll = l/2;
    dlh = (long) pow(d, lh);
    uint64_t * all_h = malloc(dlh * sizeof(uint64_t));

    dll = (long) pow(d, ll);

    /* This ``extra'' value is some replicated data before and after the
     * all_l array, so that we can avoid having to wrap around in the
     * tight pointer loop. We want to have all_h[0] matched with values
     * which, besides the positions all_l[0..x[, are also found in
     * all_l[dll..dll+x[ for some x. Afterwards, the pointers within
     * all_l merely have to go downwards while the pointer within all_h
     * goes upward.
     */
    int extra = 2 * eps * 2.5 * dll;
    uint64_t * all_l = malloc((dll + 2 * extra) * sizeof(uint64_t));


    delta = -clock();

    save_all_sums(all_h, g, d, lh);
    save_all_sums(all_l + extra, g + lh, d, ll);

    /* wrap around, so that we can simplify the
     * inner loop */
    for(i = 0 ; i < extra ; i++) {
        all_l[extra + dll+i]=all_l[extra + i];
        all_l[i]=all_l[dll + i];
    }

    uint64_t * lv;
    uint64_t * rv;

    rv = all_l + dll + 2 * extra - 1;

    if (!(*all_h + *rv >= lim)) {
        fprintf(stderr, "Bad heuristic !\n");
        exit(1);
    }

    /* This integer stores the difference between the two pointers within
     * all_l
     */
    int pd = 0;

    for(lv = all_h ; dlh-- ; lv++) {
        /* arrange so that *rv is the furthermost value < 0 */
        for( ; ((int64_t) (*lv + *rv)) >= 0 ; --rv, pd++);
        /* arrange so that rv[pd] is the furthermost value >=0 and < lim,
         * but no further than rv */
        for( ; pd && *lv + rv[pd] >= lim ; pd--);
        /* The difference between the two pointers is exactly the number
         * of solutions.
         */
        found += pd;
    }

    delta += clock();

    for(i = 0 ; i < l ; i++) {
        free(g[i]);
    }
    free(g);

    return found;
}



double** mkrand(int d, int l)
{
    double ** f;
    int i,j;
    f = malloc(l * sizeof(double*));
    for(i = 0 ; i < l ; i++) {
        f[i] = malloc(d * sizeof(double));
        for(j = 0 ; j < d ; j++) {
            f[i][j] = drand48();
        }
    }
    return f;
}

void clearrand(double ** f, int d __attribute__((unused)), int l)
{
    int i;
    for(i = 0 ; i < l ; i++) {
        free(f[i]);
    }
    free(f);
}

int main(int argc, char * argv[])
{
    int d;
    int l;
    double eps = 1.0e-9;
    long seed = getpid();

    if (argc < 3) {
        fprintf(stderr, "Usage: ./play <d> <l> [<epsilon> [<seed>]]\n");
        exit(1);
    }

    d = atoi(argv[1]);
    l = atoi(argv[2]);
    if (argc >= 4) {
        eps = atof(argv[3]);
        if (argc >= 5) {
            seed = atoi(argv[4]);
        }
    }

    srand48(seed);
    double ** f = mkrand(d, l);
    int r;
    // r = naive_search(f, d, l, eps);
    printf("%d -- %.2f\n", r, (double) (delta)/CLOCKS_PER_SEC);

    // r = search2(f, d, l, eps);
    printf("%d -- %.2f\n", r, (double) (delta)/CLOCKS_PER_SEC);
    
    r = search3(f, d, l, eps);
    printf("%d -- %.2f\n", r, (double) (delta)/CLOCKS_PER_SEC);

    clearrand(f,d,l);
    return 0;
}

