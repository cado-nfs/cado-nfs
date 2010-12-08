#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <gmp.h>

#include "random_generation.h"

#ifdef  USE_GMP_RANDOM
static gmp_randstate_t random_state;
mp_limb_t myrand()
{
    return gmp_urandomb_ui(random_state, GMP_LIMB_BITS);
}
void myseed(unsigned long int x)
{
    gmp_randseed_ui(random_state, x);
}
#else
mp_limb_t myrand()
{
    return random();
}
void myseed(unsigned long int x)
{
    srand(x);
}
#endif

void myrand_area(void * p, size_t s)
{
    for(size_t i = 0 ; i < s ; i++) {
        ((char*)p)[i]=myrand();
    }
#if 0
    /* FIXME This does not work with seeding ! */
    FILE * f = fopen("/dev/urandom", "r");
    BUG_ON(f == NULL);
    size_t r = fread(p, 1, s, f);
    BUG_ON(r != s);
    fclose(f);
#endif
}

void setup_seeding(unsigned int s)
{
#ifdef  USE_GMP_RANDOM
    gmp_randinit_mt(random_state);
#endif
    if (s == 0) {
        myseed(time(NULL));
    } else {
        myseed(s);
    }
}
