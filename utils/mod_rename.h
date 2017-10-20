/* Renames mod_ function names so that the function for the desired
   arithmetic type/width is used */

#undef mod_intinit
#undef mod_intclear
#undef mod_intset
#undef mod_intset_ul
#undef mod_intset_uls
#undef mod_intget_ul
#undef mod_intget_uls
#undef mod_intget_double
#undef mod_intequal
#undef mod_intequal_ul
#undef mod_intcmp
#undef mod_intcmp_ul
#undef mod_intcmp_uint64
#undef mod_intfits_ul
#undef mod_intadd
#undef mod_intsub
#undef mod_intbits
#undef mod_intshr
#undef mod_intshl
#undef mod_intdivexact
#undef mod_intmod
#undef mod_init
#undef mod_init_noset0
#undef mod_clear
#undef mod_set
#undef mod_set_ul
#undef mod_set_ul_reduced
#undef mod_set_int
#undef mod_set_int_reduced
#undef mod_swap
#undef mod_initmod_ul
#undef mod_initmod_int
#undef mod_getmod_ul
#undef mod_getmod_int
#undef mod_clearmod
#undef mod_get_ul
#undef mod_get_int
#undef mod_equal
#undef mod_is0
#undef mod_is1
#undef mod_add
#undef mod_add1
#undef mod_add_ul
#undef mod_sub
#undef mod_sub_ul
#undef mod_neg
#undef mod_mul
#undef mod_sqr
#undef mod_div2
#undef mod_div3
#undef mod_div5
#undef mod_div7
#undef mod_div11
#undef mod_div13
#undef mod_pow_ul
#undef mod_2pow_ul
#undef mod_pow_mp
#undef mod_2pow_mp
#undef mod_V_ul
#undef mod_V_mp
#undef mod_sprp
#undef mod_sprp2
#undef mod_isprime
#undef mod_gcd
#undef mod_inv
#undef mod_batchinv
#undef mod_jacobi
#undef mod_set0
#undef mod_set1
#undef mod_next
#undef mod_finished

#define mod_intinit          MOD_RENAME(intinit)
#define mod_intclear         MOD_RENAME(intclear)
#define mod_intset           MOD_RENAME(intset)
#define mod_intset_ul        MOD_RENAME(intset_ul)
#define mod_intset_uls       MOD_RENAME(intset_uls)
#define mod_intget_ul        MOD_RENAME(intget_ul)
#define mod_intget_uls       MOD_RENAME(intget_uls)
#define mod_intget_double    MOD_RENAME(intget_double)
#define mod_intequal         MOD_RENAME(intequal)
#define mod_intequal_ul      MOD_RENAME(intequal_ul)
#define mod_intcmp           MOD_RENAME(intcmp)
#define mod_intcmp_ul        MOD_RENAME(intcmp_ul)
#define mod_intcmp_uint64    MOD_RENAME(intcmp_uint64)
#define mod_intfits_ul       MOD_RENAME(intfits_ul)
#define mod_intadd           MOD_RENAME(intadd)
#define mod_intsub           MOD_RENAME(intsub)
#define mod_intbits          MOD_RENAME(intbits)
#define mod_intshr           MOD_RENAME(intshr)
#define mod_intshl           MOD_RENAME(intshl)
#define mod_intdivexact      MOD_RENAME(intdivexact)
#define mod_intmod           MOD_RENAME(intmod)
#define mod_init             MOD_RENAME(init)
#define mod_init_noset0      MOD_RENAME(init_noset0)
#define mod_clear            MOD_RENAME(clear)
#define mod_set              MOD_RENAME(set)
#define mod_set_ul           MOD_RENAME(set_ul)
#define mod_set_ul_reduced   MOD_RENAME(set_ul_reduced)
#define mod_set_int          MOD_RENAME(set_int)
#define mod_set_int_reduced  MOD_RENAME(set_int_reduced)
#define mod_swap             MOD_RENAME(swap)
#define mod_initmod_ul       MOD_RENAME(initmod_ul)
#define mod_initmod_int      MOD_RENAME(initmod_int)
#define mod_getmod_ul        MOD_RENAME(getmod_ul)
#define mod_getmod_int       MOD_RENAME(getmod_int)
#define mod_clearmod         MOD_RENAME(clearmod)
#define mod_get_ul           MOD_RENAME(get_ul)
#define mod_get_int          MOD_RENAME(get_int)
#define mod_equal            MOD_RENAME(equal)
#define mod_is0              MOD_RENAME(is0)
#define mod_is1              MOD_RENAME(is1)
#define mod_add              MOD_RENAME(add)
#define mod_add1             MOD_RENAME(add1)
#define mod_add_ul           MOD_RENAME(add_ul)
#define mod_sub              MOD_RENAME(sub)
#define mod_sub_ul           MOD_RENAME(sub_ul)
#define mod_neg              MOD_RENAME(neg)
#define mod_mul              MOD_RENAME(mul)
#define mod_sqr              MOD_RENAME(sqr)
#define mod_div2             MOD_RENAME(div2)
#define mod_div3             MOD_RENAME(div3)
#define mod_div5             MOD_RENAME(div5)
#define mod_div7             MOD_RENAME(div7)
#define mod_div11            MOD_RENAME(div11)
#define mod_div13            MOD_RENAME(div13)
#define mod_pow_ul           MOD_RENAME(pow_ul)
#define mod_2pow_ul          MOD_RENAME(2pow_ul)
#define mod_pow_mp           MOD_RENAME(pow_mp)
#define mod_2pow_mp          MOD_RENAME(2pow_mp)
#define mod_V_ul             MOD_RENAME(V_ul)
#define mod_V_mp             MOD_RENAME(V_mp)
#define mod_sprp             MOD_RENAME(sprp)
#define mod_sprp2            MOD_RENAME(sprp2)
#define mod_isprime          MOD_RENAME(isprime)
#define mod_gcd              MOD_RENAME(gcd)
#define mod_inv              MOD_RENAME(inv)
#define mod_batchinv         MOD_RENAME(batchinv)
#define mod_jacobi           MOD_RENAME(jacobi)
#define mod_set0             MOD_RENAME(set0)
#define mod_set1             MOD_RENAME(set1)
#define mod_next             MOD_RENAME(next)
#define mod_finished         MOD_RENAME(finished)

/* A function that is not used anywhere. The purpose is solely generating 
   compilation errors if any of the renamed functions, which constitute  
   kind of a definition of the API, are not implemented. */
static inline void * 
mod_test_if_functions_exist()
{
  void *p;
  p = (void*) &mod_intinit;
  p = (void*) &mod_intclear;
  p = (void*) &mod_intset;
  p = (void*) &mod_intset_ul;
  p = (void*) &mod_intset_uls;
  p = (void*) &mod_intget_ul;
  p = (void*) &mod_intget_uls;
  p = (void*) &mod_intequal;
  p = (void*) &mod_intequal_ul;
  p = (void*) &mod_intcmp;
  p = (void*) &mod_intcmp_ul;
  p = (void*) &mod_intcmp_uint64;
  p = (void*) &mod_intfits_ul;
  p = (void*) &mod_intadd;
  p = (void*) &mod_intsub;
  p = (void*) &mod_intbits;
  p = (void*) &mod_intshr;
  p = (void*) &mod_intshl;
  p = (void*) &mod_intdivexact;
  p = (void*) &mod_intmod;
  p = (void*) &mod_init;
  p = (void*) &mod_init_noset0;
  p = (void*) &mod_clear;
  p = (void*) &mod_set;
  p = (void*) &mod_set_ul;
  p = (void*) &mod_set_ul_reduced;
  p = (void*) &mod_set_int;
  p = (void*) &mod_set_int_reduced;
  p = (void*) &mod_swap;
  /* This is implemented only for 1-word arithmetic.
     FIXME: Since these functions are not generally available, they probably 
     should not be renamed, either.
  p = &mod_initmod_ul; 
  p = &mod_getmod_ul;
  */
  p = (void*) &mod_initmod_int;
  p = (void*) &mod_getmod_int;
  p = (void*) &mod_clearmod;
  p = (void*) &mod_get_ul;
  p = (void*) &mod_get_int;
  p = (void*) &mod_equal;
  p = (void*) &mod_is0;
  p = (void*) &mod_is1;
  p = (void*) &mod_add;
  p = (void*) &mod_add1;
  p = (void*) &mod_add_ul;
  p = (void*) &mod_sub;
  p = (void*) &mod_sub_ul;
  p = (void*) &mod_neg;
  p = (void*) &mod_mul;
  p = (void*) &mod_sqr;
  p = (void*) &mod_div2;
  p = (void*) &mod_div3;
  p = (void*) &mod_div5;
  p = (void*) &mod_div7;
  p = (void*) &mod_div11;
  p = (void*) &mod_div13;
  p = (void*) &mod_pow_ul;
  p = (void*) &mod_2pow_ul;
  p = (void*) &mod_pow_mp;
  p = (void*) &mod_2pow_mp;
  p = (void*) &mod_V_ul;
  p = (void*) &mod_V_mp;
  p = (void*) &mod_sprp;
  p = (void*) &mod_sprp2;
  p = (void*) &mod_isprime;
  p = (void*) &mod_gcd;
  p = (void*) &mod_inv;
  p = (void*) &mod_batchinv;
  p = (void*) &mod_jacobi;
  p = (void*) &mod_set0;
  p = (void*) &mod_set1;
  p = (void*) &mod_next;
  p = (void*) &mod_finished;
  return p;
}
