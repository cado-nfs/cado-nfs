#ifndef MODREDC_15UL_DEFAULT_H_
#define MODREDC_15UL_DEFAULT_H_
#include "modredc_15ul.h"

// pragma no prototypes
/*
   Here are typedef's that rename all functions to mod_* instead of 
   modredc15ul_*, which one day might become an automatic renaming scheme so 
   that different modulus sizes can be used simply by #including different 
   mod*.h files, but without changing anything else in the source code. 
 */

#ifdef  mod_init
#warning "modredc_15ul_default.h included after another mod*_default.h ; bindings overwritten"
#endif

#undef MOD_RENAME
#define MOD_RENAME(x) MODREDC15UL_RENAME(x)
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
#define residue_t            residueredc15ul_t
#define modulus_t            modulusredc15ul_t
#define modint_t             modintredc15ul_t
#define MOD_SIZE             MODREDC15UL_SIZE
#define MOD_MINBITS          MODREDC15UL_MINBITS
#define MOD_MAXBITS          MODREDC15UL_MAXBITS
#define MOD_APPEND_TYPE(x)   x##_15ul
#define PRIMODu "Nd"
#define PRIMODx "Nx"
#define MOD_PRINT_INT(x) x, MOD_SIZE
#define MOD_PRINT_MODULUS(x) x[0].m, MOD_SIZE

/* This function is used in mod_2ul_common.h with different 
   implementations for 15ul and 2ul2. Maybe there's a prettier way */
#undef mod_divn
#define mod_divn modredc15ul_divn

#endif	/* MODREDC_15UL_DEFAULT_H_ */
