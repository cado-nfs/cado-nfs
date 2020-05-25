#include <stdio.h>
#ifndef EC_ARITH_COMMON_H_
#define EC_ARITH_COMMON_H_

#ifndef mod_init
  #error "One of the mod*_default.h headers must be included before this file"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* Types of coordinates */
typedef enum {
  SHORT_WEIERSTRASS_aff,
  SHORT_WEIERSTRASS_proj,
  MONTGOMERY_xz,
  TWISTED_EDWARDS_proj,
  TWISTED_EDWARDS_ext,
} ec_point_coord_type_t;


/* A point on an elliptic curve:
 *   The significant coordinates are:
 *      - x,y     for SHORT_WEIERSTRASS_aff
 *      - x,y,z   for SHORT_WEIERSTRASS_proj
 *      - x,  z   for MONTGOMERY_xz
 *      - x,y,z   for TWISTED_EDWARDS_proj
 *      - x,y,z,t for TWISTED_EDWARDS_ext
 */
struct ec_point_s
{
  residue_t x,y,z,t;
};

typedef struct ec_point_s ec_point_t[1];
typedef struct ec_point_s *ec_point_ptr;
typedef const struct ec_point_s *ec_point_srcptr;

#define ec_point_init MOD_APPEND_TYPE(ec_point_init)
static inline void
ec_point_init (ec_point_t P, const modulus_t m)
{
  mod_init (P->x, m);
  mod_init (P->y, m);
  mod_init (P->z, m);
  mod_init (P->t, m);
}

#define ec_point_init_noset0 MOD_APPEND_TYPE(ec_point_init_noset0)
static inline void
ec_point_init_noset0 (ec_point_t P, const modulus_t m)
{
  mod_init_noset0 (P->x, m);
  mod_init_noset0 (P->y, m);
  mod_init_noset0 (P->z, m);
  mod_init_noset0 (P->t, m);
}

#define ec_point_clear MOD_APPEND_TYPE(ec_point_clear)
static inline void
ec_point_clear (ec_point_t P, const modulus_t m)
{
  mod_clear (P->x, m);
  mod_clear (P->y, m);
  mod_clear (P->z, m);
  mod_clear (P->t, m);
}

#define ec_point_set MOD_APPEND_TYPE(ec_point_set)
static inline void
ec_point_set (ec_point_t Q, const ec_point_t P, const modulus_t m,
              const ec_point_coord_type_t coord)
{
  mod_set (Q->x, P->x, m);
  if (coord != MONTGOMERY_xz)
    mod_set (Q->y, P->y, m);
  if (coord != SHORT_WEIERSTRASS_aff)
    mod_set (Q->z, P->z, m);
  if (coord == TWISTED_EDWARDS_ext)
    mod_set (Q->t, P->t, m);
}

#define ec_point_swap MOD_APPEND_TYPE(ec_point_swap)
static inline void
ec_point_swap (ec_point_t Q, ec_point_t P, const modulus_t m,
               const ec_point_coord_type_t coord)
{
  mod_swap (Q->x, P->x, m);
  if (coord != MONTGOMERY_xz)
    mod_swap (Q->y, P->y, m);
  if (coord != SHORT_WEIERSTRASS_aff)
    mod_swap (Q->z, P->z, m);
  if (coord == TWISTED_EDWARDS_ext)
    mod_swap (Q->t, P->t, m);
}

#define ec_point_fprintf MOD_APPEND_TYPE(ec_point_fprintf)
static inline void
ec_point_fprintf (FILE *out, const ec_point_t P,
                  const ec_point_coord_type_t coord, const modulus_t m)
{
  modint_t x, y, z, t;

  mod_intinit (x);
  mod_intinit (y);
  mod_intinit (z);
  mod_intinit (t);
  
  mod_get_int (x, P->x, m);
  if (coord != MONTGOMERY_xz)
    mod_get_int (y, P->y, m);
  if (coord != SHORT_WEIERSTRASS_aff)
    mod_get_int (z, P->z, m);
  if (coord == TWISTED_EDWARDS_ext)
    mod_get_int (t, P->t, m);

  switch (coord)
  {
    case SHORT_WEIERSTRASS_aff: /* (x, y) */
      mod_fprintf (out, "(0x%" PRIMODx ", 0x%" PRIMODx ")", MOD_PRINT_INT(x),
          MOD_PRINT_INT(y));
      break;
    case TWISTED_EDWARDS_proj: /* (X : Y : Z) */
    case SHORT_WEIERSTRASS_proj: /* (X : Y : Z) */
      mod_fprintf (out, "(0x%" PRIMODx " : 0x%" PRIMODx " : 0x%" PRIMODx ")",
          MOD_PRINT_INT(x), MOD_PRINT_INT(y), MOD_PRINT_INT(z));
      break;
    case MONTGOMERY_xz: /* (X :: Z) */
      mod_fprintf (out, "(0x%" PRIMODx " :: 0x%" PRIMODx ")", MOD_PRINT_INT(x),
          MOD_PRINT_INT(z));
      break;
    case TWISTED_EDWARDS_ext: /* (X : Y : Z : T) */
      mod_fprintf (out, "(0x%" PRIMODx " : 0x%" PRIMODx " : 0x%" PRIMODx " : "
          "0x%" PRIMODx ")", MOD_PRINT_INT(x), MOD_PRINT_INT(y),
          MOD_PRINT_INT(z), MOD_PRINT_INT(t));
      break;
    default:
      fprintf (stderr, "%s: unknown coordinates system\n", __func__);
      abort ();
  }

  mod_intclear (x);
  mod_intclear (y);
  mod_intclear (z);
  mod_intclear (t);
}

#ifdef __cplusplus
}
#endif

#endif /* EC_ARITH_COMMON_H_ */
