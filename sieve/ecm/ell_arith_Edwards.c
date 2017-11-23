#include "ell_arith.h"

/* #define SAFE_EDW_TO_MONTG */


/* (x : y : z) to (x : y : t : z)
   cost: 3m+1s by computing (xz, yz, xy, z^2) */
inline void
edwards_proj_to_ext (ell_point_t Q, const ell_point_t P, const modulus_t m)
{
  mod_mul (Q->x, P->x, P->z, m);
  mod_mul (Q->y, P->y, P->z, m);
  mod_mul (Q->t, P->x, P->y, m);
  mod_sqr (Q->z, P->z, m);
}

inline void
edwards_neg (ell_point_t P, const modulus_t m)
{
  mod_neg (P->x, P->x, m);
  mod_neg (P->t, P->t, m);
}

/* - edwards_add (R:output_flag, P:edwards_ext, Q:edwards_ext, output_flag) */
/*     R <- P+Q */
/*     output_flag can be edwards_proj, edwards_ext or montgomery */
void
edwards_add (ell_point_t R, const ell_point_t P, const ell_point_t Q,
	  const modulus m, const ell_point_coord_type output_type)
{
  /* The "add-2008-hwcd-4" addition formulas */
  /* Cost: 8M + 8add + 2*2. */
  /* Cost: 8M + 6add dependent upon the first point. */
  /* Source: 2008 Hisil–Wong–Carter–Dawson */ 
  /* http://eprint.iacr.org/2008/522, Section 3.2. */

  residue_t t0, t1, A, B, C, D, E, F, G, H;
  
#if COUNT_ELLE_OPS
  edwards_add_count++;
#endif
  
  /* ASSERT (output_flag != AFF); */

  mod_init_noset0 (t0, m);
  mod_init_noset0 (t1, m);
  mod_init_noset0 (A, m);
  mod_init_noset0 (B, m);
  mod_init_noset0 (C, m);
  mod_init_noset0 (D, m);
  mod_init_noset0 (E, m);
  mod_init_noset0 (F, m);
  mod_init_noset0 (G, m);
  mod_init_noset0 (H, m);

  mod_sub (t0, P->y, P->x, m);     // t0 := (Y1-X1)
  mod_add (t1, Q->y, Q->x, m);     // t1 := (Y2+X2)
  mod_mul (A, t0, t1, m);          // A := (Y1-X1)*(Y2+X2) 
  mod_add (t0, P->y, P->x, m);     // t0 := (X1+Y1)
  mod_add (t1, Q->y, Q->x, m);     // t1 := (Y2-X2)
  mod_mul (B, t0, t1, m);          // B := (Y1+X1)*(Y2-X2)
  mod_mul (C, P->z, Q->t, m);      // C := Z1*T2
  mod_add (C, C, C, m);            // C := 2*Z1*T2
  mod_mul (D, P->t, Q->z, m);      // D := T1*Z2
  mod_add (D, D, D, m);            // D := 2*T1*Z2
  mod_add (E, D, C, m);            // E := D+C
  mod_sub (F, B, A, m);            // F := B-A
  mod_add (G, B, A, m);            // G := B+A
  mod_sub (H, D, C, m);            // H := D-C

  if (output_type != MONTG)
    {
      mod_mul (R->x, E, F, m);     // X3 := E*F
      mod_mul (R->y, G, H, m);     // Y3 := G*H
      mod_mul (R->z, F, G, m);     // Z3 := F*G
      if (output_type == EDW_ext)
	mod_mul (R->t, E, H, m);   // T3 := E*H
    }
  else
    {
#ifdef SAFE_EDW_TO_MONTG
      mod_add (R->x, F, H, m);
      mod_mul (R->x, R->x, E, m);
      mod_sub (R->z, F, H, m);
      mod_mul (R->z, R->z, E, m);
#else
      /* CAUTION! */
      /* This may produce unstable results */
      /* But seems to "work" for our purpose */
      /* TODO: COMMENTS */
      mod_add (R->x, F, H, m);
      mod_sub (R->z, F, H, m);
#endif
    }
  
  mod_clear (t0, m);
  mod_clear (t1, m);
  mod_clear (A, m);
  mod_clear (B, m);
  mod_clear (C, m);
  mod_clear (D, m);
  mod_clear (E, m);
  mod_clear (F, m);
  mod_clear (G, m);
  mod_clear (H, m);  
}


/* - edwards_sub (R:output_flag, P:edwards_ext, Q:edwards_ext, output_flag) */
/*     R <- P-Q */
/*     output_flag can be edwards_proj, edwards_ext or montgomery */
void
edwards_sub (ell_point_t R, const ell_point_t P, const ell_point_t Q,
	  const modulus m, const ell_point_coord_type output_type)
{
  ell_point_t QQ;
  ell_point_init (QQ, m);

  edwards_neg (QQ, Q, m);
  edwards_add (R, P, QQ, m);

  ell_point_clear (QQ, m);
}


/* - edwards_dbl (R:output_flag, P:edwards_proj, output_flag) */
/*     R <- 2*P */
/*     output_flag can be edwards_proj, edwards_ext */
void
edwards_dbl (ell_point_t R, const ell_point_t P,
	  const modulus_t m, const ell_point_coord_type output_type)
{
  /* The "dbl-2008-hwcd" doubling formulas */
  /* Cost: 4M + 4S + 1*a + 6add + 1*2. */
  /* Source: 2008 Hisil–Wong–Carter–Dawson */
  /* http://eprint.iacr.org/2008/522, Section 3.3. */
    
  residue_t A, B, C, D, E, F, G, H;

#if COUNT_ELLE_OPS
  edwards_double_count++;
#endif

  /* ASSERT (output_flag != AFF); */
  /* ASSERT (output_flag != MONTG); */
  
  mod_init_noset0 (A, m);
  mod_init_noset0 (B, m);
  mod_init_noset0 (C, m);
  mod_init_noset0 (D, m);
  mod_init_noset0 (E, m);
  mod_init_noset0 (F, m);
  mod_init_noset0 (G, m);
  mod_init_noset0 (H, m);

  mod_sqr (A, P->x, m);                // A := X1^2
  mod_sqr (B, P->y, m);                // B := Y1^2
  mod_sqr (C, P->z, m);                // C := Z1^2
  mod_add (C, C, C, m);                // C := 2*Z1^2
  mod_mul (D, A, a, m);                // D := a*A
  mod_add (E, P->x, P->y, m);          // E := (X1 + Y1)
  mod_sqr (E, E, m);                   // E := (X1 + Y1)^2
  mod_sub (E, E, A, m);                // E := (X1 + Y1)^2-A
  mod_sub (E, E, B, m);                // E := (X1 + Y1)^2-A-B
  mod_add (G, D, B, m);                // G := D+B
  mod_sub (F, G, C, m);                // F := G-C
  mod_sub (H, D, B, m);                // H := D-B
  mod_mul (R->x, E, F, m);             // X3 := E*F
  mod_mul (R->y, G, H, m);             // Y3 := G*H
  mod_mul (R->z, F, G, m);             // Z3 := F*G
  if (output_flag == EDW_ext)
    mod_mul (R->t, E, H, m);           // T3 := E*H
  
  mod_clear (A, m);
  mod_clear (B, m);
  mod_clear (C, m);
  mod_clear (D, m);
  mod_clear (E, m);
  mod_clear (F, m);
  mod_clear (G, m);
  mod_clear (H, m);  
}


/* - edwards_tpl (R:output_flag, P:edwards_proj, output_flag) */
/*     R <- 3*P */
/*     output_flag can be edwards_proj, edwards_ext */

/* The "tpl-2015-c" tripling formulas */
/* Cost: 11M + 3S + 1*a + 7add + 2*2. */
/* Source: 2015 Chuengsatiansup. */
/* https://hyperelliptic.org/EFD/g1p/auto-twisted-extended-1.html#tripling-tpl-2015-c */
void
edwards_tpl (ell_point_t R, const ell_point_t P,
	  const modulus_t m, const ell_point_coord_type output_type)
{

  residue_t YY, aXX, Ap, B, xB, yB, AA, F, G, xE, yH, zF, zG;

#if COUNT_ELLE_OPS
  edwards_triple_count++;
#endif

  /* ASSERT (output_flag != AFF); */
  /* ASSERT (output_flag != MONTG); */
  
  mod_init_noset0 (YY, m);
  mod_init_noset0 (aXX, m);
  mod_init_noset0 (Ap, m);
  mod_init_noset0 (B, m);
  mod_init_noset0 (xB, m);
  mod_init_noset0 (yB, m);
  mod_init_noset0 (AA, m);
  mod_init_noset0 (F, m);
  mod_init_noset0 (G, m);
  mod_init_noset0 (xE, m);
  mod_init_noset0 (yH, m);
  mod_init_noset0 (zF, m);
  mod_init_noset0 (zG, m);

  mod_sqr (YY, P->y, m);                // YY := Y1^2
  mod_sqr (aXX, P->x, m);               // aXX := X1^2
  mod_neg (aXX, aXX, m);                // aXX := -X1^2
  mod_add (Ap, YY, aXX, m);             // Ap := YY+aXX 
  mod_sqr (B, P->z, m);                 // B := Z1^2
  mod_add (B, B, B, m);                 // B := 2*Z1^2
  mod_sub (B, B, Ap, m);                // B := 2*Z1^2-Ap
  mod_add (B, B, B, m);                 // B := 2*(2*Z1^2-Ap)
  mod_mul (xB, aXX, B, m);              // xB := aXX*B
  mod_mul (yB, YY, B, m);               // yB := YY*B
  mod_sub (AA, YY, aXX, m);             // AA := YY-aXX
  mod_mul (AA, Ap, AA, m);              // AA := Ap*(YY-aXX)
  mod_sub (F, AA, yB, m);               // F := AA-yB
  mod_add (G, AA, xB, m);               // G := AA+xB

  mod_add (xE, yB, AA, m);              // xE := yB+AA
  mod_mul (xE, P->x, xE, m);            // xE := X1*(yB+AA)
  mod_sub (yH, xB, AA, m);              // yH := xB-AA
  mod_mul (yH, P->y, yH, m);            // yH := Y1*(xB-AA)
  mod_mul (zF, P->z, F, m);             // zF := Z1*F
  
  switch (output_type)
    {
    case EDW_proj :
      mod_mul (R->x, xE, F, m);         // X3 := xE*F
      mod_mul (R->y, yH, G, m);         // Y3 := yH*G
      mod_mul (R->z, zF, G, m);         // Z3 := zF*G
      break;
    case EDW_ext :
      mod_mul (zG, P->z, G, m);         // zG := Z1*G
      mod_mul (R->x, xE, zF, m);        // X3 := xE*zF
      mod_mul (R->y, yH, zG, m);        // Y3 := yH*zG
      mod_mul (R->z, zF, zG, m);        // Z3 := zF*zG
      mod_mul (R->t, xE, yH, m);        // T3 := xE*yH
      break;
    }
  
  mod_clear (YY, m);
  mod_clear (aXX, m);
  mod_clear (Ap, m);
  mod_clear (B, m);
  mod_clear (xB, m);
  mod_clear (yB, m);
  mod_clear (AA, m);
  mod_clear (F, m);
  mod_clear (G, m);
  mod_clear (xE, m);
  mod_clear (yH, m);
  mod_clear (zF, m);
  mod_clear (zG, m);
}