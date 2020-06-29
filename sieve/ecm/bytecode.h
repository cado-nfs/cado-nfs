#ifndef BYTECODE_H_
#define BYTECODE_H_

#include <stdint.h>
#include <gmp.h>

typedef uint8_t bytecode_elt;
typedef uint8_t * bytecode;
typedef const uint8_t * bytecode_const;

#ifdef __cplusplus
extern "C" {
#endif

/* No-op (because bytecode_elt is an uint8_t for now). */
static inline uint8_t
bytecode_elt_to_uint8 (bytecode_elt b)
{
  return b;
}

static inline bytecode_elt
bytecode_elt_build_4_4 (uint8_t b1, uint8_t b0)
{
  return (b1 << 4) | (b0 & 0x0f);
}

static inline void
bytecode_elt_split_4_4 (uint8_t *b1, uint8_t *b0, bytecode_elt b)
{
  if (b0)
    *b0 = b & 0x0f;
  if (b1)
    *b1 = b >> 4;
}

static inline void
bytecode_elt_split_2_1_1_4 (uint8_t *b3, uint8_t *b2, uint8_t *b1, uint8_t *b0,
                            bytecode_elt b)
{
  if (b0)
    *b0 = b & 0x0f;
  b >>= 4;
  if (b1)
    *b1 = b & 0x01;
  b >>= 1;
  if (b2)
    *b2 = b & 0x01;
  b >>= 1;
  if (b3)
    *b3 = b & 0x03;
}


/******************************************************************************/
/************************************ cost ************************************/
/******************************************************************************/

/* See README.bytecode for more details on the names of the fields */
struct mishmash_cost_s
{
  double DBL;
  double DBLa;
  double TPL;
  double TPLa;
  double ADD;
  double ADDa;
  double ADDd;
  double dDBL;
  double dADD;
};

/******************************************************************************/
/***************************** double base chains *****************************/
/******************************************************************************/
#define DBCHAIN_OP_DBLADD 0x01
#define DBCHAIN_OP_TPLADD 0x02
#define DBCHAIN_OP_TPLDBLADD 0x03

typedef struct mishmash_cost_s dbchain_cost_t;

/******************************************************************************/
/******************************* precomputation *******************************/
/******************************************************************************/
#define PRECOMP_FINAL 0xff
#define PRECOMP_OP_ADD 0x00
#define PRECOMP_OP_DBL 0x01
#define PRECOMP_OP_TPL 0x02

typedef struct mishmash_cost_s precomp_cost_t;

/******************************************************************************/
/************************************ PRAC ************************************/
/******************************************************************************/
#define PRAC_SUBBLOCK_INIT 0x69 /* 'i' */
#define PRAC_SUBBLOCK_FINAL 0x66 /* 'f' */
#define PRAC_BLOCK_FINAL 0x46 /* 'F' */
#define PRAC_SWAP 0x73 /* 's' */

struct prac_cost_s
{
  double DBL; /* cost of a doubling */
  double dADD; /* cost of a differential addition */
};

typedef struct prac_cost_s prac_cost_t;

void bytecode_prac_encode (bytecode *, unsigned int, unsigned int, unsigned int,
                           const prac_cost_t *, int, int);
int bytecode_prac_check (bytecode_const, mpz_srcptr, int);
void bytecode_prac_cache_free ();


/******************************************************************************/
/********************************** MISHMASH **********************************/
/******************************************************************************/
#define MISHMASH_INIT 0x0
#define MISHMASH_DBCHAIN_BLOCK 0x1
#define MISHMASH_PRECOMP_BLOCK 0x2
#define MISHMASH_PRAC_BLOCK 0x8
#define MISHMASH_FINAL 0xff

typedef struct mishmash_cost_s mishmash_cost_t;

void bytecode_mishmash_encode (bytecode *, unsigned int, unsigned int,
                               unsigned int, const mishmash_cost_t *, int, int);
int bytecode_mishmash_check (bytecode_const, mpz_srcptr, int);

#ifdef __cplusplus
}
#endif



#endif /* BYTECODE_H_ */
