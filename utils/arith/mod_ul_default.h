#ifndef CADO_MOD_UL_DEFAULT_H
#define CADO_MOD_UL_DEFAULT_H

#include "mod_ul.h"

// pragma no prototypes
// scan-headers: skip

/*
   Here are typedef's that rename all functions to mod_* instead of 
   modul_*, which one day might become an automatic renaming scheme so 
   that different modulus sizes can be used simply by #including different 
   mod*.h files, but without changing anything else in the source code. 
 */

#ifdef  mod_init
#warning "mod_ul_default.h included after another mod*_default.h ; bindings overwritten"
#endif

#undef MOD_RENAME
#define MOD_RENAME(x) MODUL_RENAME(x)
#include "mod_rename.h"

#undef residue_t
#undef modulus_t
#undef modint_t
#undef MOD_SIZE
#undef MOD_MINBITS
#undef MOD_MAXBITS
#define residue_t            residueul_t
#define modulus_t            modulusul_t
#define modint_t             modintul_t
#define MOD_SIZE             MODUL_SIZE
#define MOD_MAXBITS          MODUL_MAXBITS

#endif	/* CADO_MOD_UL_DEFAULT_H */
