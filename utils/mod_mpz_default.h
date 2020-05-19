#ifndef MOD_MPZ_DEFAULT_H_
#define MOD_MPZ_DEFAULT_H_

#include "mod_mpz.h"

// pragma no prototypes
/*
   Here are typedef's that rename all functions to mod_* instead of 
   modmpz_*, which one day might become an automatic renaming scheme so 
   that different modulus sizes can be used simply by #including different 
   mod*.h files, but without changing anything else in the source code. 
 */

#ifdef  mod_init
#warning "mod_mpz_default.h included after another mod*_default.h ; bindings overwritten"
#endif

#undef MOD_RENAME
#define MOD_RENAME(x) MODMPZ_RENAME(x)
#include "mod_rename.h"

#undef residue_t
#undef modulus_t
#undef modint_t
#undef MOD_SIZE
#undef MOD_MINBITS
#undef MOD_MAXBITS
#undef MOD_APPEND_TYPE
#undef PRIMODu
#undef PRIMODx
#undef MOD_PRINT_INT
#undef MOD_PRINT_MODULUS
#define residue_t            residuempz_t
#define modulus_t            modulusmpz_t
#define modint_t             modintmpz_t
#define MOD_MAXBITS          MODMPZ_MAXBITS
#define MOD_APPEND_TYPE(x)   x##_mpz
#define PRIMODu "Zd"
#define PRIMODx "Zx"
#define MOD_PRINT_INT(x) x
#define MOD_PRINT_MODULUS(m) m

#endif	/* MOD_MPZ_DEFAULT_H_ */
