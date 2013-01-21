/*****************************************************************
*                Functions for the factor base                  *
*****************************************************************/

#ifndef FB_H
#define FB_H

#include <stdio.h>
#include <stddef.h>
#include <string.h>
#include <gmp.h>

/* Data types */

typedef unsigned int fbprime_t; /* 32 bits should be enough for everyone */
#define FBPRIME_FORMAT "%u"
#define FBPRIME_MAX 4294967295U
#define FBPRIME_BITS 32
typedef fbprime_t fbroot_t;
#define FBROOT_FORMAT "%u"
typedef unsigned long largeprime_t; /* On IA32 they'll only get 32 bit 
                                       large primes */
#define LARGEPRIME_FORMAT "%lu"

/* Factor base entry with (possibly) several roots */
typedef struct {
  fbprime_t p;            /* A prime or a prime power */
  unsigned char nr_roots; /* how many roots there are for this prime */
  unsigned char plog;     /* logarithm (to some suitable base) of this prime */
  unsigned char size;     /* The length of the struct in bytes */
  unsigned char dummy[1]; /* For dword aligning the roots. In fact uneeded, C99
                             guarantees proper alignment of roots[]. It's only a
                             precaution against -fpack-struct or other
                             ABI-breaking behaviours */
  unsigned long invp;     /* -1/p (mod 2^wordsize) for REDC: although we need
			     only a 32-bit inverse in say redc_32, we need a
			     full-limb inverse on 64-bit machines for trial
			     division */
  /* Note that invp has a stronger alignment constraint than p, thus must
   * not appear before the tiny fields plog and nr_roots which can easily
   * fit inbetween the two */
  fbroot_t roots[0];      /* the actual length of this array is determined
                             by nr_roots */
} factorbase_degn_t;

#define FB_END ((fbprime_t) 1)

void		fb_fprint_entry (FILE *, const factorbase_degn_t *);
void            fb_fprint (FILE *, const factorbase_degn_t *);
void            fb_sortprimes (fbprime_t *, const unsigned int);
unsigned char	fb_log (double, double, double);
factorbase_degn_t * 	fb_make_linear (const mpz_t *, const fbprime_t, 
					const fbprime_t, const double, 
					const int, const int, FILE *);
factorbase_degn_t *	fb_read (const char *, const double, const int, const fbprime_t, fbprime_t);
fbprime_t	*fb_extract_bycost (const factorbase_degn_t *, 
                                    const fbprime_t, const fbprime_t costlim);
size_t          fb_size (const factorbase_degn_t *);                   
size_t          fb_nroots_total (const factorbase_degn_t *fb);

/* Some inlined functions which need to be fast */
  
/* Hack to get around C's automatic multiplying constants to be added to 
   pointers by the pointer's base data type size */

__attribute__ ((unused))
static inline factorbase_degn_t *
fb_skip (const factorbase_degn_t *fb, const size_t s)
{
  return (factorbase_degn_t *)((char *)fb + s);
}
  
__attribute__ ((unused))
static inline factorbase_degn_t *
fb_next (const factorbase_degn_t *fb)
{
  return (factorbase_degn_t *)((char *)fb + fb->size);
}

__attribute__ ((unused))
static size_t
fb_entrysize (const factorbase_degn_t *fb)
{
  return (sizeof (factorbase_degn_t) + fb->nr_roots * sizeof (fbroot_t));
}

__attribute__ ((unused))
static size_t
fb_entrysize_uc (const unsigned char n)
{
  return (sizeof (factorbase_degn_t) + n * sizeof (fbroot_t));
}

/* Most often, we're in fact considering const iterators, but we refuse
 * to bother with having two interfaces */

struct fb_iterator_s {
    factorbase_degn_t * fb;
    int i;
};

typedef struct fb_iterator_s fb_iterator[1];
typedef struct fb_iterator_s *fb_iterator_ptr;
typedef const struct fb_iterator_s *fb_iterator_srcptr;

static inline void fb_iterator_init_set_fb(fb_iterator_ptr t, factorbase_degn_t * fb)
{
    memset(t, 0, sizeof(*t));
    t->fb = fb;
}

static inline void fb_iterator_init_set(fb_iterator_ptr t, fb_iterator_srcptr u)
{
    memcpy(t, u, sizeof(*u));
}

static inline void fb_iterator_clear(fb_iterator_ptr t)
{
    memset(t, 0, sizeof(*t));
}

static inline void fb_iterator_next(fb_iterator_ptr t)
{
    if (++(t->i) < t->fb->nr_roots)
        return;
    t->fb = fb_next(t->fb);
    t->i = 0;
}

static inline void fb_iterator_add_n(fb_iterator_ptr t, int n)
{
    if (t->i) {
        n += t->i;
        t->i = 0;
    }
    for( ; t->fb->nr_roots <= n ; ) {
        n -= t->fb->nr_roots;
        t->fb = fb_next(t->fb);
    }
    t->i = n;
}

static inline fbprime_t fb_iterator_get_r(fb_iterator_srcptr t)
{
    return t->fb->roots[t->i];
}

static inline int fb_iterator_lessthan_fb(fb_iterator_srcptr t, const factorbase_degn_t * fb)
{
    return (char*)(t->fb) < (char*)fb;
}

static inline int fb_iterator_lessthan(fb_iterator_srcptr t, fb_iterator_srcptr u)
{
    int r = (char*)(u->fb) - (char*)(t->fb);
    if (r > 0) return 1;
    if (r < 0) return 0;
    return t->i < u->i;
}

/* Computes t - u. Assumes for the primary branch that u <= t. If t < u,
 * symmetrises by swapping arguments */
static inline int fb_iterator_diff(fb_iterator_srcptr t, fb_iterator_srcptr u)
{
    fb_iterator q;
    if (fb_iterator_lessthan(t, u)) {
        return -fb_iterator_diff(u, t);
    }
    fb_iterator_init_set(q, u);
    int n = -q->i;
    q->i = 0;
    for( ; fb_iterator_lessthan_fb(q, t->fb) ; ) {
        n += q->fb->nr_roots;
        q->fb = fb_next(q->fb);
    }
    n += t->i;
    fb_iterator_clear(q);
    return n;
}
static inline int fb_iterator_diff_fb(fb_iterator_srcptr t, factorbase_degn_t * u)
{
    fb_iterator qu;
    fb_iterator_init_set_fb(qu, u);
    int n = fb_iterator_diff(t, qu);
    fb_iterator_clear(qu);
    return n;
}

static inline int fb_diff(factorbase_degn_t * t, factorbase_degn_t * u)
{
    fb_iterator qu;
    fb_iterator_init_set_fb(qu, u);
    fb_iterator qt;
    fb_iterator_init_set_fb(qt, t);
    int n = fb_iterator_diff(qt, qu);
    fb_iterator_clear(qt);
    fb_iterator_clear(qu);
    return n;
}
static inline ptrdiff_t fb_iterator_diff_bytes(fb_iterator_srcptr t, fb_iterator_srcptr u)
{
    fb_iterator q;
    if (fb_iterator_lessthan(t, u)) {
        return -fb_iterator_diff_bytes(u, t);
    }
    fb_iterator_init_set(q, u);
    ptrdiff_t n = -q->i * sizeof(fbprime_t);
    q->i = 0;
    for( ; fb_iterator_lessthan_fb(q, t->fb) ; ) {
        n += q->fb->nr_roots * sizeof(fbprime_t) + sizeof(factorbase_degn_t);
        q->fb = fb_next(q->fb);
    }
    n += t->i;
    fb_iterator_clear(q);
    return n;
}
static inline ptrdiff_t fb_iterator_diff_bytes_fb(fb_iterator_srcptr t, factorbase_degn_t * u)
{
    fb_iterator qu;
    fb_iterator_init_set_fb(qu, u);
    ptrdiff_t n = fb_iterator_diff_bytes(t, qu);
    fb_iterator_clear(qu);
    return n;
}

static inline ptrdiff_t fb_diff_bytes(factorbase_degn_t * t, factorbase_degn_t * u)
{
    fb_iterator qu;
    fb_iterator_init_set_fb(qu, u);
    fb_iterator qt;
    fb_iterator_init_set_fb(qt, t);
    ptrdiff_t n = fb_iterator_diff_bytes(qt, qu);
    fb_iterator_clear(qt);
    fb_iterator_clear(qu);
    return n;
}
static inline int fb_iterator_over(fb_iterator_srcptr t)
{
    return t->fb->p == FB_END;
}

#endif
