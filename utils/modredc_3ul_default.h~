#include "modredc_3ul.h"

/*
   Here are typedef's that rename all functions to mod_* instead of 
   modredc3ul_*, which one day might become an automatic renaming scheme so 
   that different modulus sizes can be used simply by #including different 
   mod*.h files, but without changing anything else in the source code. 
 */

#ifdef  mod_init
#warning "modredc_3ul_default.h included after another mod*_default.h ; bindings overwritten"
#endif

#undef MOD_RENAME
#define MOD_RENAME(x) MODREDC3UL_RENAME(x)
#include "mod_rename.h"

#undef residue_t
#undef modulus_t
#undef modint_t
#undef MOD_SIZE
#undef MOD_MINBITS
#undef MOD_MAXBITS
#define residue_t            residueredc3ul_t
#define modulus_t            modulusredc3ul_t
#define modint_t             modintredc3ul_t
#define MOD_SIZE             MODREDC3UL_SIZE
#define MOD_MINBITS          MODREDC3UL_MINBITS
#define MOD_MAXBITS          MODREDC3UL_MAXBITS

