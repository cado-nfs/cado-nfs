#ifndef EC_ARITH_COST_H_
#define EC_ARITH_COST_H_

// pragma no prototypes

/* For "a=-1" twisted Edwards curves */
#define EDWARDS_DBL 7.            /* doubling projective -> projective */
#define EDWARDS_DBLext 8.         /* doubling projective -> extended */
#define EDWARDS_TPL 12.           /* tripling projective -> projective */
#define EDWARDS_TPLext 14.        /* tripling projective -> extended */
#define EDWARDS_ADD 7.            /* addition extended,extended -> projective */
#define EDWARDS_ADDext 8.         /* addition extended,extended -> extended */
#define EDWARDS_ADDmontgomery 4.  /* addition extended,extended -> Montgomery */

/* For Montgomery curves */
#define MONTGOMERY_DBL 5.        /* doubling */
#define MONTGOMERY_dADD 6.       /* differential addition */

#endif /* EC_ARITH_COST_H_ */
