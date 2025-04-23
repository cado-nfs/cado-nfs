#ifndef CADO_ROPT_TREE_H
#define CADO_ROPT_TREE_H

#include <stdint.h>     // int16_t
#include <gmp.h>


/**
 * Struct for the nodes in the lift.
 */
typedef struct node_t {
  // char e;
  float val;
  unsigned int nr;
  unsigned int u;
  unsigned int v;
  unsigned int alloc;
  unsigned int *r;
  char *roottype;
  struct node_t *firstchild;
  struct node_t *nextsibling;
  struct node_t *parent;
} node; // sizeof = 64


/**
 * Priority queue to record sublattices (w, u, v)'s alpha's.
 */
typedef struct alpha_pq_t {
  int *w;
  mpz_t *u;
  mpz_t *v;
  mpz_t *modulus;
  float *alpha;
  int len;
  int used;
} alpha_pq;


/**
 * Priority queue to record alpha values.
 */
typedef struct sievescore_pq_t {
  long *i;
  long *j;
  int16_t *alpha;
  int len;
  int used;
} sievescore_pq;


/**
 * Priority queue to record E.
 */
typedef struct MurphyE_pq_t {
  int *w;
  mpz_t *u;
  mpz_t *v;
  mpz_t *modulus;
  double *E;
  int len;
  int used;
} MurphyE_pq;


/* --- declarations --- */

#ifdef __cplusplus
extern "C" {
#endif


/* tree, used in ropt_stage1.c */
void new_tree ( node **root );

node* new_node ( void );

void insert_node ( node *parent,
                   node **currnode,
                   unsigned int u,
                   unsigned int v,
                   unsigned int r,
                   char curr_e,
                   unsigned int p,
                   unsigned int pe,
                   char k );

void free_node ( node **ptr );

/* alpha_pq, used in ropt_stage1.c */
void new_alpha_pq ( alpha_pq **ppqueue,
                    unsigned long len );

void insert_alpha_pq ( alpha_pq *pqueue, 
                       int w,
                       mpz_srcptr u,
                       mpz_srcptr v,
                       mpz_srcptr modulus,
                       double alpha );

void extract_alpha_pq ( alpha_pq *pqueue,
                        int *w,
                        mpz_ptr u,
                        mpz_ptr v,
                        mpz_ptr modulus,
                        double *alpha );

void reset_alpha_pq ( alpha_pq *pqueue );

void remove_rep_alpha ( alpha_pq *pqueue );

void free_alpha_pq ( alpha_pq **ppqueue );

/* sievescore_pq, used in ropt_stage2.c */
void new_sievescore_pq ( sievescore_pq **ppqueue,
                         unsigned long len );

void reset_sievescore_pq ( sievescore_pq *pqueue );

void insert_sievescore_pq ( sievescore_pq *pqueue,
                            long i,
                            long j,
                            int16_t alpha );

void free_sievescore_pq ( sievescore_pq **ppqueue );

/* MurphyE_pq, used in ropt_stage2.c */
void new_MurphyE_pq ( MurphyE_pq **ppqueue,
                      unsigned long len );

void insert_MurphyE_pq ( MurphyE_pq *pqueue,
                         int w,
                         mpz_srcptr u,
                         mpz_srcptr v,
                         mpz_srcptr modulus,
                         double E );

void extract_MurphyE_pq ( MurphyE_pq *pqueue,
                          int *w,
                          mpz_ptr u,
                          mpz_ptr v,
                          mpz_ptr modulus,
                          double *E );

void free_MurphyE_pq ( MurphyE_pq **ppqueue );

void reset_MurphyE_pq ( MurphyE_pq *pqueue );

void remove_rep_MurphyE ( MurphyE_pq *pqueue );


#ifdef __cplusplus
}
#endif


#endif /* CADO_ROPT_TREE_H */
