/**
 * @file ropt_tree.c
 * Liftings and priority queue structs.
 */


#include "cado.h" // IWYU pragma: keep
#include <float.h>      // FLT_MAX DBL_MAX
#include <stdio.h>      // fprintf stderr
#include <stdlib.h>     // free malloc
#include <stdint.h>     /* AIX wants it first (it's a bug) */
#include <gmp.h>
#include "ropt_tree.h"


#if 0 /* Leave them for debug */
/**
 * Print the info for the node
 */
static inline void
print_node ( node *pnode )
{
  fprintf(stderr, "(%u,%u):%u:%d:(%.2f)\n",
          pnode->u, pnode->v, pnode->nr, pnode->e, pnode->val);
  /* int i; */
  /* for (i = 0; i < pnode->nr; i++) */
  /* printf ("pnode->r[%d]: %lu\n", i, pnode->r[i]); */
}

/**
 * Print a tree, non-recursive. Two styles.
 * - If level == 0, print the whole tree;
 * - If level == i, print nodes at height i.
 */
static inline void
print_tree ( node *root,
             int level )
{
  if (level < 0) {
    fprintf(stderr,"Error: level >= 0. \n");
    exit(1);
  }

  int i, currlevel = 0;
  node *ptr = root;
  printf ("\n");

  if (!ptr)  /* if empty */
    return;

  /* find leftmost */
  while (ptr->firstchild) {
    ptr = ptr->firstchild;
    ++ currlevel;
  }
  while (ptr) {
    if (level != 0) {
      if (currlevel == level) {
        for (i = 0; i < currlevel; ++i)
          printf("-   ");
        print_node(ptr);
      }
    }
    else {
      for (i = 0; i < currlevel; ++i)
        printf("-   ");
      print_node(ptr);
    }

    if (ptr->nextsibling) {
      ptr = ptr->nextsibling;
      //++ level; // aligned layout, uncommont for sloped layout
      while (ptr->firstchild) {
        ptr = ptr->firstchild;
        ++ currlevel;
      }
    }
    else {
      ptr = ptr->parent;  /* move up */
      -- currlevel;
    }
  }
  printf ("\n");
  return;
}
#endif


/**
 * Malloc or realloc (double) the memory space of pnode->r[]
 * and pnode->roottype[].
*/
void
alloc_r_node ( node *pnode )
{
  //assert (pnode->alloc <= pnode->nr);

  if (pnode->nr == 0) {
    pnode->r = (unsigned int *)
      malloc ( sizeof (unsigned int) );
    pnode->roottype = (char *)
      malloc ( sizeof (char) );
    pnode->alloc = 1;
  }
  else {
    pnode->r = (unsigned int *)
      realloc ( pnode->r, 2 * pnode->nr * sizeof (unsigned int) );
    pnode->roottype = (char *)
      realloc ( pnode->roottype, 2 * pnode->nr * sizeof (char) );
    pnode->alloc = (pnode->alloc * 2);
  }

  if (pnode->r == NULL) {
    fprintf(stderr, "Error, cannot allocate memory in %s\n", __func__);
    exit (1);
  }
}


/**
 * Some indices, note the queue is shifted by 1.
 */
static inline int
pq_parent ( int i )
{
  return (i>>1);
}


/**
 * Initialise an empty tree.
 */
void
new_tree ( node **root )
{
  *root = NULL;
}


/**
 * Create new empty node.
 */
node *
new_node ( void )
{
  node *pnode;
  pnode = (node *) malloc(sizeof(node));
  if (pnode == NULL) {
    fprintf(stderr,"Error: malloc failed in new_node()\n");
    exit(1);
  }
  pnode->firstchild = pnode->nextsibling = pnode->parent = NULL;
  pnode->r = NULL;
  pnode->roottype = NULL;
  pnode->u = pnode->v = pnode->nr = 0; // pnode->e = 0;
  pnode->val = 0.0;
  return pnode;
}


/**
 *  Insert child node to the current node, where *parent is
 *  the ptr to the parent of the current node.
 *
 *  - If (u, v) exits in any childrent of the parent;
 *  -- Check if r exists;
 *  --- If not, add r and/or "is_multiple=k".
 *  - If (u, v) doesnot exit;
 *  -- Add a node with (u, v, r, curr_e, val) and/or "is_multiple=k".
 */
#define DEBUG_INSERT_NODE 0
void
insert_node ( node *parent,
              node **currnode,
              unsigned int u,
              unsigned int v,
              unsigned int r,
              char curr_e,
              unsigned int p,
              unsigned int pe,
              char k )
{
  unsigned int i;
  node *lastnode = NULL;
  node *nextnode = NULL;
  nextnode = parent->firstchild;
  lastnode = nextnode;

  /* add r to some existing node (u, v) */
  while ((nextnode != NULL)) {

    if ((nextnode->u == u) && (nextnode->v == v)) {

      /* if (u,v) pair exists, see whether
         the root r already exists in the node. */
      for (i = 0; i < nextnode->nr; i ++) {
        if ( nextnode->r[i] == r )
          return;
      }

      /* otherwise, insert this root */
      if (nextnode->alloc <= nextnode->nr)
        alloc_r_node (nextnode);
      nextnode->roottype[nextnode->nr] = k;
      nextnode->r[nextnode->nr] = r;
      nextnode->nr += 1;

      if (k == 2)
        nextnode->val += 1.0 / (double) pe;
      else {
        if (curr_e == 1)
          nextnode->val += 1.0 / ((double)p - 1.0);
      }

#if DEBUG_INSERT_NODE
      fprintf ( stderr, "add to (u: %u, v: %u), nr: %u, "
                "r: %u, e: %d, val: %f\n",
                u, v, nextnode->nr, nextnode->r[nextnode->nr-1],
                curr_e, nextnode->val );
#endif
      break;
    }
    lastnode = nextnode;
    nextnode = nextnode->nextsibling;
  }

  /* add new (u, v) node */
  if (nextnode == NULL) {

    *currnode = new_node();

    /* has left sibling? */
    if (lastnode != NULL)
      lastnode->nextsibling = (*currnode);
    else
      parent->firstchild = (*currnode);

    (*currnode)->parent = parent;
    (*currnode)->u = u;
    (*currnode)->v = v;

    alloc_r_node (*currnode);
    (*currnode)->r[0] = r;
    (*currnode)->roottype[0] = k;
    (*currnode)->nr += 1;
    // (*currnode)->e = curr_e;
    /* such (u, v) is new, hence we inheritate the val from its parent. */
    if (k == 2)
      (*currnode)->val += 1.0 / (double) pe + (*currnode)->parent->val;
    else {
      if (curr_e == 1)
        (*currnode)->val += 1.0 / ((double)p - 1.0);
      else
        (*currnode)->val = (*currnode)->parent->val;
    }

#if DEBUG_INSERT_NODE
    fprintf ( stderr, "creating (u: %u, v: %u), nr: %u, "
              "r: %u, e: %d, val: %f\n", u, v,
              (*currnode)->nr, (*currnode)->r[0], curr_e,
              (*currnode)->val );
#endif
  }
  // print_tree(parent, 0);
}


/**
 * Free current node ptr.
 */
void
free_node ( node **ptr )
{
  if (*ptr) {
    if ((*ptr)->r)
      free ((*ptr)->r);
    if ((*ptr)->roottype)
      free ((*ptr)->roottype);
    (*ptr)->parent = (*ptr)->firstchild =
      (*ptr)->nextsibling = NULL;
    free (*ptr);
  }
}

/**
 * Create priority queue for sublattices with best alpha.
 */
void
new_alpha_pq ( alpha_pq **ppqueue,
               unsigned long len )
{
  if ( len < 2 ) {
    fprintf(stderr,"Error: len < 2 in new_alpha_pq()\n");
    exit(1);
  }

  (*ppqueue) = (alpha_pq *) malloc (sizeof (alpha_pq));
  if ( (*ppqueue) == NULL) {
    fprintf(stderr,"Error: malloc failed in new_alpha_pq()\n");
    exit(1);
  }

  (*ppqueue)->len = len;
  (*ppqueue)->w = (int *) malloc (len* sizeof (int));
  (*ppqueue)->alpha = (float *) malloc (len* sizeof (float));
  (*ppqueue)->u = (mpz_t *) malloc (len* sizeof (mpz_t));
  (*ppqueue)->v = (mpz_t *) malloc (len* sizeof (mpz_t));
  (*ppqueue)->modulus = (mpz_t *) malloc (len* sizeof (mpz_t));

  if ( (*ppqueue)->w == NULL ||
       (*ppqueue)->alpha == NULL ||
       (*ppqueue)->modulus == NULL ||
       (*ppqueue)->u == NULL ||
       (*ppqueue)->v == NULL ) {
    fprintf(stderr,"Error: malloc failed in new_sub_alpha_pq()\n");
    exit(1);
  }

  int i;
  for (i = 0; i < (*ppqueue)->len; i++)
  {
    mpz_init ( (*ppqueue)->u[i] );
    mpz_init ( (*ppqueue)->v[i] );
    mpz_init ( (*ppqueue)->modulus[i] );
  }

  mpz_set_ui ( (*ppqueue)->u[0], 0 );
  mpz_set_ui ( (*ppqueue)->v[0], 0 );
  mpz_set_ui ( (*ppqueue)->modulus[0], 0 );
  (*ppqueue)->alpha[0] = FLT_MAX;
  (*ppqueue)->w[0] = 0;
  (*ppqueue)->used = 1;
}


/**
 * Free
 */
void
free_alpha_pq ( alpha_pq **ppqueue )
{
  int i;
  for (i = 0; i < (*ppqueue)->len; i++)
  {
    mpz_clear ( (*ppqueue)->u[i] );
    mpz_clear ( (*ppqueue)->v[i] );
    mpz_clear ( (*ppqueue)->modulus[i] );
  }

  free ( (*ppqueue)->u );
  free ( (*ppqueue)->v );
  free ( (*ppqueue)->modulus );
  free ( (*ppqueue)->w );
  free ( (*ppqueue)->alpha );
  free ( *ppqueue );
}


/**
 * Sift-up to add, if the queue is not full.
 */
static inline void
insert_alpha_pq_up ( alpha_pq *pqueue,
                     int w,
                     mpz_srcptr u,
                     mpz_srcptr v,
                     mpz_srcptr modulus,
                     double alpha )
{
  int k;

  for ( k = pqueue->used;
        alpha > pqueue->alpha[pq_parent(k)];
        k /= 2 )
  {
    mpz_set ( pqueue->u[k], pqueue->u[pq_parent(k)] );
    mpz_set ( pqueue->v[k], pqueue->v[pq_parent(k)] );
    mpz_set ( pqueue->modulus[k], pqueue->modulus[pq_parent(k)] );
    pqueue->w[k] = pqueue->w[pq_parent(k)];
    pqueue->alpha[k] = pqueue->alpha[pq_parent(k)];
  }
  mpz_set (pqueue->u[k], u);
  mpz_set (pqueue->v[k], v);
  mpz_set (pqueue->modulus[k], modulus);
  pqueue->w[k] = w;
  pqueue->alpha[k] = alpha;
  pqueue->used ++;
}


/**
 * Sift-down, if the heap is full.
 */
static inline void
insert_alpha_pq_down ( alpha_pq *pqueue,
                       int w,
                       mpz_srcptr u,
                       mpz_srcptr v,
                       mpz_srcptr modulus,
                       double alpha )
{
  int k, l;

  for (k = 1; k*2 < pqueue->used; k = l) {

    l = (k << 1);

    /* right > left ? */
    if ( (l+1) < pqueue->used &&
         (pqueue->alpha[l+1] > pqueue->alpha[l]) )
      l ++;

    /* switch larger child with parent */
    if ( pqueue->alpha[l] > alpha ) {
      mpz_set (pqueue->u[k], pqueue->u[l]);
      mpz_set (pqueue->v[k], pqueue->v[l]);
      mpz_set (pqueue->modulus[k], pqueue->modulus[l]);
      pqueue->w[k] = pqueue->w[l];
      pqueue->alpha[k] = pqueue->alpha[l];
    }
    else
      break;
  }

  mpz_set (pqueue->u[k], u);
  mpz_set (pqueue->v[k], v);
  mpz_set (pqueue->modulus[k], modulus);
  pqueue->w[k] = w;
  pqueue->alpha[k] = alpha;
}


/**
 * Extract the max of the priority queue.
 */
void
extract_alpha_pq ( alpha_pq *pqueue,
                   int *w,
                   mpz_ptr u,
                   mpz_ptr v,
                   mpz_ptr modulus,
                   double *alpha )
{
  /* don't extract u[0] since it is just a placeholder. */
  pqueue->used --;
  mpz_set (u, pqueue->u[1]);
  mpz_set (v, pqueue->v[1]);
  mpz_set (modulus, pqueue->modulus[1]);
  *alpha = pqueue->alpha[1];
  *w = pqueue->w[1];

  insert_alpha_pq_down ( pqueue,
                         pqueue->w[pqueue->used],
                         pqueue->u[pqueue->used],
                         pqueue->v[pqueue->used],
                         pqueue->modulus[pqueue->used],
                         pqueue->alpha[pqueue->used] );
}


/**
 * Insert to the priority queue.
 */
void
insert_alpha_pq ( alpha_pq *pqueue,
                  int w,
                  mpz_srcptr u,
                  mpz_srcptr v,
                  mpz_srcptr modulus,
                  double alpha )
{

  /* queue is full,  */
  if (pqueue->len == pqueue->used) {
    if ( pqueue->alpha[1] > alpha) {
      insert_alpha_pq_down (pqueue, w, u, v, modulus, alpha);
    }
  }

  /* queue is not full, sift-up */
  else if (pqueue->len > pqueue->used) {
    insert_alpha_pq_up (pqueue, w, u, v, modulus, alpha);
  }
  else {
    fprintf(stderr,"Error: error (pqueue->len < pqueue->used) "
            "in insert_alpha_pq()\n");
    exit(1);
  }
}


/**
 * Reset priority queue for alpha.
 */
void
reset_alpha_pq ( alpha_pq *pqueue )
{

  if ( pqueue == NULL) {
    fprintf(stderr,"Error: malloc failed in reset_alpha_pq()\n");
    exit(1);
  }

  if ( pqueue->w == NULL || pqueue->u == NULL ||
       pqueue->v == NULL || pqueue->modulus == NULL ||
       pqueue->alpha == NULL ) {
    fprintf(stderr,"Error: malloc failed in reset_alpha_pq()\n");
    exit(1);
  }

  mpz_set_ui ( pqueue->u[0], 0 );
  mpz_set_ui ( pqueue->v[0], 0 );
  mpz_set_ui ( pqueue->modulus[0], 0 );
  pqueue->w[0] = 0;
  pqueue->alpha[0] = FLT_MAX;
  pqueue->used = 1;
}


/**
 * Remove duplicates in alpha queue.
 */
void
remove_rep_alpha ( alpha_pq *pqueue )
{
  int flag;
  alpha_pq *new;
  new_alpha_pq (&new, pqueue->len);
  for (int i=1; i<pqueue->used; i++) {
    flag = 0;
    for (int k=i+1; k<pqueue->used; k++) {
      if ( (pqueue->w[i] == pqueue->w[k]) &&
           mpz_cmp(pqueue->u[i],pqueue->u[k])==0 &&
           mpz_cmp(pqueue->v[i],pqueue->v[k])==0 &&
           mpz_cmp(pqueue->modulus[i], pqueue->modulus[k])==0) {
        flag = 1;
        break;
      }
    }
    if (flag==0) {
      insert_alpha_pq (new, pqueue->w[i], pqueue->u[i], pqueue->v[i],
                       pqueue->modulus[i], pqueue->alpha[i]);
    }
  }
  reset_alpha_pq (pqueue);
  pqueue->used = new->used;
  for (int i=1; i<pqueue->used; i++) {
    pqueue->w[i] = new->w[i];
    mpz_set(pqueue->u[i],new->u[i]);
    mpz_set(pqueue->v[i],new->v[i]);
    mpz_set(pqueue->modulus[i],new->modulus[i]);
    pqueue->alpha[i] = new->alpha[i];
  }
  free_alpha_pq (&new);
}


/**
 * Create priority queue for best root scores (in MAT[] for a sublattice).
 */
void
new_sievescore_pq ( sievescore_pq **ppqueue,
                    unsigned long len )
{
  if ( len < 2 ) {
    fprintf(stderr,"Error: len < 2 in new_sievescore_pq()\n");
    exit(1);
  }

  (*ppqueue) = (sievescore_pq *) malloc (sizeof (sievescore_pq));
  if ( (*ppqueue) == NULL) {
    fprintf(stderr,"Error: malloc failed in new_sievescore_pq()\n");
    exit(1);
  }

  (*ppqueue)->len = len;

  (*ppqueue)->i = (long *) malloc (len* sizeof (long));
  (*ppqueue)->j = (long *) malloc (len* sizeof (long));
  (*ppqueue)->alpha = (int16_t *) malloc (len* sizeof (int16_t));

  if ( (*ppqueue)->i == NULL ||
       (*ppqueue)->j == NULL ||
       (*ppqueue)->alpha == NULL ) {
    fprintf(stderr,"Error: malloc failed in new_sievescore_pq()\n");
    exit(1);
  }

  /* i[0] and j[0] are null elements */
  (*ppqueue)->i[0] = 0;
  (*ppqueue)->j[0] = 0;
  (*ppqueue)->alpha[0] = INT16_MAX;
  (*ppqueue)->used = 1;
}


/**
 * Reset priority queue for best root scores (in MAT[] for a sublattice).
 */
void
reset_sievescore_pq ( sievescore_pq *pqueue )
{

  if ( pqueue == NULL) {
    fprintf(stderr,"Error: malloc failed in new_sievescore_pq()\n");
    exit(1);
  }

  if ( pqueue->i == NULL ||
       pqueue->j == NULL ||
       pqueue->alpha == NULL ) {
    fprintf(stderr,"Error: malloc failed in new_sievescore_pq()\n");
    exit(1);
  }

  /* i[0] and j[0] are null elements */
  pqueue->i[0] = 0;
  pqueue->j[0] = 0;
  pqueue->alpha[0] = INT16_MAX;
  pqueue->used = 1;
}


/**
 * Free
 */
void
free_sievescore_pq ( sievescore_pq **ppqueue )
{
  free ( (*ppqueue)->i );
  free ( (*ppqueue)->j );
  free ( (*ppqueue)->alpha );
  free ( *ppqueue );
}


/**
 * Sift-up to add, if the queue is not full.
 */
static inline void
insert_sievescore_pq_up ( sievescore_pq *pqueue,
                          long i,
                          long j,
                          int16_t alpha )
{
  int k;

  for ( k = pqueue->used;
        alpha > pqueue->alpha[pq_parent(k)];
        k /= 2 )
  {
    pqueue->i[k] = pqueue->i[pq_parent(k)];
    pqueue->j[k] = pqueue->j[pq_parent(k)];
    pqueue->alpha[k] = pqueue->alpha[pq_parent(k)];
  }

  pqueue->i[k] = i;
  pqueue->j[k] = j;
  pqueue->alpha[k] = alpha;

  pqueue->used ++;
}


/**
 * Sift-down, if the heap is full.
 */
static inline void
insert_sievescore_pq_down ( sievescore_pq *pqueue,
                            long i,
                            long j,
                            int16_t alpha )
{
  int k, l;

  for (k = 1; k*2 < pqueue->used; k = l) {

    l = (k << 1);

    /* right > left ? */
    if ( (l+1) < pqueue->used &&
         (pqueue->alpha[l+1] > pqueue->alpha[l]) )
      l ++;

    /* switch larger child with parent */
    if ( pqueue->alpha[l] > alpha ) {
      pqueue->i[k] = pqueue->i[l];
      pqueue->j[k] = pqueue->j[l];
      pqueue->alpha[k] = pqueue->alpha[l];
    }
    else
      break;
  }

  pqueue->i[k] = i;
  pqueue->j[k] = j;
  pqueue->alpha[k] = alpha;
}


/**
 * Insert to the priority queue.
 */
void
insert_sievescore_pq ( sievescore_pq *pqueue,
                       long i,
                       long j,
                       int16_t alpha )
{
  /* fprintf (stderr, "# Debug: inserting (%ld, %ld), "
     "alpha: %d, used: %d, len: %d\n", */

  /* queue is full,  */
  if (pqueue->len == pqueue->used) {
    if ( alpha < pqueue->alpha[1] ) {

      insert_sievescore_pq_down (pqueue, i, j, alpha);
    }
  }
  /* queue is not full, sift-up */
  else if (pqueue->len > pqueue->used) {
    insert_sievescore_pq_up (pqueue, i, j, alpha);
  }
  else {
    fprintf(stderr,"Error: error (pqueue->len < pqueue->used) "
            "in insert_sublattice_pq()\n");
    exit(1);
  }
}


/**
 * Create priority queue for MurphyE
 */
void
new_MurphyE_pq ( MurphyE_pq **ppqueue,
                 unsigned long len )
{
  if ( len < 2 ) {
    fprintf(stderr,"Error: len < 2 in new_MurphyE_pq()\n");
    exit(1);
  }

  (*ppqueue) = (MurphyE_pq *) malloc (sizeof (MurphyE_pq));
  if ( (*ppqueue) == NULL) {
    fprintf(stderr,"Error: malloc failed in new_MurphyE_pq()\n");
    exit(1);
  }

  (*ppqueue)->len = len;
  (*ppqueue)->w = (int *) malloc (len* sizeof (int));
  (*ppqueue)->u = (mpz_t *) malloc (len* sizeof (mpz_t));
  (*ppqueue)->v = (mpz_t *) malloc (len* sizeof (mpz_t));
  (*ppqueue)->modulus = (mpz_t *) malloc (len* sizeof (mpz_t));
  (*ppqueue)->E = (double *) malloc (len* sizeof (double));

  if ( (*ppqueue)->u == NULL || (*ppqueue)->v == NULL ||
       (*ppqueue)->w == NULL || (*ppqueue)->E == NULL ||
       (*ppqueue)->modulus == NULL ) {
    fprintf(stderr,"Error: malloc failed in new_MurphyE_pq()\n");
    exit(1);
  }

  int i;
  for (i = 0; i < (*ppqueue)->len; i++)
  {
    mpz_init ( (*ppqueue)->u[i] );
    mpz_init ( (*ppqueue)->v[i] );
    mpz_init ( (*ppqueue)->modulus[i] );
  }

  mpz_set_ui ( (*ppqueue)->u[0], 0 );
  mpz_set_ui ( (*ppqueue)->v[0], 0 );
  mpz_set_ui ( (*ppqueue)->modulus[0], 0 );
  (*ppqueue)->w[0] = 0;
  (*ppqueue)->E[0] = -DBL_MAX;
  (*ppqueue)->used = 1;
}


/**
 * Free
 */
void
free_MurphyE_pq ( MurphyE_pq **ppqueue )
{

  int i;
  for (i = 0; i < (*ppqueue)->len; i++) {
    mpz_clear ( (*ppqueue)->u[i] );
    mpz_clear ( (*ppqueue)->v[i] );
    mpz_clear ( (*ppqueue)->modulus[i] );
  }

  free ( (*ppqueue)->w );
  free ( (*ppqueue)->u );
  free ( (*ppqueue)->modulus );
  free ( (*ppqueue)->v );
  free ( (*ppqueue)->E );
  free ( *ppqueue );
}


/**
 * Sift-up to add, if the queue is not full.
 */
static inline void
insert_MurphyE_pq_up ( MurphyE_pq *pqueue,
                       int w,
                       mpz_srcptr u,
                       mpz_srcptr v,
                       mpz_srcptr modulus,
                       double E )
{
  int k;

  for ( k = pqueue->used;
        E < pqueue->E[pq_parent(k)];
        k /= 2 )
  {
    mpz_set ( pqueue->u[k], pqueue->u[pq_parent(k)] );
    mpz_set ( pqueue->v[k], pqueue->v[pq_parent(k)] );
    mpz_set ( pqueue->modulus[k], pqueue->modulus[pq_parent(k)] );
    pqueue->w[k] = pqueue->w[pq_parent(k)];
    pqueue->E[k] = pqueue->E[pq_parent(k)];
  }

  mpz_set (pqueue->u[k], u);
  mpz_set (pqueue->v[k], v);
  mpz_set (pqueue->modulus[k], modulus);
  pqueue->w[k] = w;
  pqueue->E[k] = E;

  pqueue->used ++;
}


/**
 * Sift-down, if the heap is full.
 */
static inline void
insert_MurphyE_pq_down ( MurphyE_pq *pqueue,
                         int w,
                         mpz_srcptr u,
                         mpz_srcptr v,
                         mpz_srcptr modulus,
                         double E )
{
  int k, l;

  for (k = 1; k*2 < pqueue->used; k = l) {

    l = (k << 1);

    /* right < left ? */
    if ( (l+1) < pqueue->used &&
         (pqueue->E[l+1] < pqueue->E[l]) )
      l ++;

    /* switch smaller child with parent */
    if ( pqueue->E[l] < E ) {

      mpz_set (pqueue->u[k], pqueue->u[l]);
      mpz_set (pqueue->v[k], pqueue->v[l]);
      mpz_set (pqueue->modulus[k], pqueue->modulus[l]);
      pqueue->w[k] = pqueue->w[l];
      pqueue->E[k] = pqueue->E[l];
    }
    else
      break;
  }

  mpz_set (pqueue->u[k], u);
  mpz_set (pqueue->v[k], v);
  mpz_set (pqueue->modulus[k], modulus);
  pqueue->w[k] = w;
  pqueue->E[k] = E;
}



/**
 * Extract the max of the priority queue.
 */
void
extract_MurphyE_pq ( MurphyE_pq *pqueue,
                     int *w,
                     mpz_ptr u,
                     mpz_ptr v,
                     mpz_ptr modulus,
                     double *E )
{
  /* don't extract u[0] since it is just a placeholder. */
  pqueue->used --;
  mpz_set (u, pqueue->u[1]);
  mpz_set (v, pqueue->v[1]);
  mpz_set (modulus, pqueue->modulus[1]);
  *w = pqueue->w[1];
  *E = pqueue->E[1];

  insert_MurphyE_pq_down ( pqueue,
                           pqueue->w[pqueue->used],
                           pqueue->u[pqueue->used],
                           pqueue->v[pqueue->used],
                           pqueue->modulus[pqueue->used],
                           pqueue->E[pqueue->used] );
}


/**
 * Insert to the priority queue.
 */
void
insert_MurphyE_pq ( MurphyE_pq *pqueue,
                    int w,
                    mpz_srcptr u,
                    mpz_srcptr v,
                    mpz_srcptr modulus,
                    double E )
{

  /* queue is full,  */
  if (pqueue->len == pqueue->used) {
    if ( E > pqueue->E[1] ) {
      insert_MurphyE_pq_down (pqueue, w, u, v, modulus, E);
    }
  }
  /* queue is not full, sift-up */
  else if (pqueue->len > pqueue->used) {
    insert_MurphyE_pq_up (pqueue, w, u, v, modulus, E);
  }
  else {
    fprintf(stderr,"Error: error (pqueue->len < pqueue->used) "
            "in insert_MurphyE_pq()\n");
    exit(1);
  }
}


/**
 * Reset priority queue for alpha.
 */
void
reset_MurphyE_pq ( MurphyE_pq *pqueue )
{

  if ( pqueue == NULL) {
    fprintf(stderr,"Error: malloc failed in reset_MurphyE_pq()\n");
    exit(1);
  }

  if ( pqueue->w == NULL || pqueue->u == NULL ||
       pqueue->v == NULL || pqueue->modulus == NULL ||
       pqueue->E == NULL ) {
    fprintf(stderr,"Error: malloc failed in reset_MurphyE_pq()\n");
    exit(1);
  }

  mpz_set_ui ( pqueue->u[0], 0 );
  mpz_set_ui ( pqueue->v[0], 0 );
  mpz_set_ui ( pqueue->modulus[0], 0 );
  pqueue->w[0] = 0;
  pqueue->E[0] = -DBL_MAX;
  pqueue->used = 1;
}


/**
 * Remove duplicates in murphyE queue.
 */
void
remove_rep_MurphyE ( MurphyE_pq *pqueue )
{
  int flag;
  MurphyE_pq *new;
  new_MurphyE_pq (&new, pqueue->len);
  for (int i=1; i<pqueue->used; i++) {
    flag = 0;
    for (int k=i+1; k<pqueue->used; k++) {
      if ( (pqueue->w[i] == pqueue->w[k]) &&
           mpz_cmp(pqueue->u[i],pqueue->u[k])==0 &&
           mpz_cmp(pqueue->v[i],pqueue->v[k])==0 &&
           mpz_cmp(pqueue->modulus[i], pqueue->modulus[k])==0) {
        flag = 1;
        break;
      }
    }
    if (flag==0) {
      insert_MurphyE_pq (new, pqueue->w[i], pqueue->u[i], pqueue->v[i],
                         pqueue->modulus[i], pqueue->E[i]);
    }
  }
  reset_MurphyE_pq (pqueue);
  pqueue->used = new->used;
  for (int i=1; i<pqueue->used; i++) {
    pqueue->w[i] = new->w[i];
    mpz_set(pqueue->u[i],new->u[i]);
    mpz_set(pqueue->v[i],new->v[i]);
    mpz_set(pqueue->modulus[i],new->modulus[i]);
    pqueue->E[i] = new->E[i];
  }
  free_MurphyE_pq (&new);
}
