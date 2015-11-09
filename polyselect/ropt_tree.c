/**
 * @file ropt_tree.c
 * Liftings and priority queue structs.
 */


#include "cado.h"
#include <stdint.h>     /* AIX wants it first (it's a bug) */
#include "ropt_tree.h"
#include "portability.h"


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
    fprintf (stderr, "Error, cannot reallocate memory in alloc_r_node().\n");
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
  pnode->u = pnode->v = pnode->nr = pnode->e = 0;
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
    (*currnode)->e = curr_e;
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
 * Create priority queue for best sublattices of length len.
 */
void
new_sublattice_pq ( sublattice_pq **ppqueue,
                    unsigned long len )
{
  if ( len < 2 ) {
    fprintf(stderr,"Error: len < 2 in new_sublattice_pq()\n");
    exit(1);
  }

  (*ppqueue) = (sublattice_pq *) malloc (sizeof (sublattice_pq));
  if ( (*ppqueue) == NULL) {
    fprintf(stderr,"Error: malloc failed in new_sublattice_pq()\n");
    exit(1);
  }

  (*ppqueue)->len = len;

  (*ppqueue)->u = (mpz_t *) malloc (len* sizeof (mpz_t));
  (*ppqueue)->v = (mpz_t *) malloc (len* sizeof (mpz_t));
  (*ppqueue)->val = (float *) malloc (len* sizeof (float));
  (*ppqueue)->modulus = (mpz_t *) malloc (len* sizeof (mpz_t));
  if ( (*ppqueue)->u == NULL || (*ppqueue)->v == NULL ||
       (*ppqueue)->v == NULL || (*ppqueue)->val == NULL ) {
    fprintf(stderr,"Error: malloc failed in new_sublattice_pq()\n");
    exit(1);
  }

  int i;
  for (i = 0; i < (*ppqueue)->len; i++) {
    mpz_init ( (*ppqueue)->u[i] );
    mpz_init ( (*ppqueue)->v[i] );
    mpz_init ( (*ppqueue)->modulus[i] );
  }

  mpz_set_ui ( (*ppqueue)->u[0], 0 );
  mpz_set_ui ( (*ppqueue)->v[0], 0 );
  mpz_set_ui ( (*ppqueue)->modulus[0], 0 );
  (*ppqueue)->val[0] = -FLT_MAX;
  (*ppqueue)->used = 1; // u[0] and v[0] are null elements
}


/**
 * Create priority queue for best sublattices of length len.
 */
void
free_sublattice_pq ( sublattice_pq **ppqueue )
{
  int i;
  for (i = 0; i < (*ppqueue)->len; i++)
  {
    mpz_clear ( (*ppqueue)->u[i] );
    mpz_clear ( (*ppqueue)->v[i] );
    mpz_clear ( (*ppqueue)->modulus[i] );
  }

  free ( (*ppqueue)->val );
  free ( (*ppqueue)->u );
  free ( (*ppqueue)->v );
  free ( (*ppqueue)->modulus );
  free ( *ppqueue );
}


/**
 * Sift-up to add, if the queue is not full.
 */
static inline void
insert_sublattice_pq_up ( sublattice_pq *pqueue,
                          mpz_t u,
                          mpz_t v,
                          mpz_t mod,
                          float val )
{
  int i;

  for ( i = pqueue->used;
        val < pqueue->val[pq_parent(i)];
        i /= 2 ) {

    mpz_set ( pqueue->u[i], pqueue->u[pq_parent(i)] );
    mpz_set ( pqueue->v[i], pqueue->v[pq_parent(i)] );
    mpz_set ( pqueue->modulus[i], pqueue->modulus[pq_parent(i)] );
    pqueue->val[i] = pqueue->val[pq_parent(i)];
  }

  mpz_set (pqueue->u[i], u);
  mpz_set (pqueue->v[i], v);
  mpz_set (pqueue->modulus[i], mod);
  pqueue->val[i] = val;

  pqueue->used ++;
}


/**
 * Sift-down, if the heap is full.
 */
static inline void
insert_sublattice_pq_down ( sublattice_pq *pqueue,
                            mpz_t u,
                            mpz_t v,
                            mpz_t mod,
                            float val )
{
  int i, l;

  for (i = 1; i*2 < pqueue->used; i = l) {

    l = (i << 1);

    /* right > left ? */
    if ( (l+1) < pqueue->used &&
         pqueue->val[l+1] <= pqueue->val[l] )
      l ++;

    /* switch larger child with parent */
    if ( (pqueue->val[l] < val) ||
         (pqueue->val[l] == val && mpz_cmp (pqueue->u[l], u) > 0) ) {
      mpz_set (pqueue->u[i], pqueue->u[l]);
      mpz_set (pqueue->v[i], pqueue->v[l]);
      mpz_set (pqueue->modulus[i], pqueue->modulus[l]);
      pqueue->val[i] = pqueue->val[l];
    }
    else
      break;
  }

  mpz_set (pqueue->u[i], u);
  mpz_set (pqueue->v[i], v);
  mpz_set (pqueue->modulus[i], mod);
  pqueue->val[i] = val;
}


/**
 * Insert to the priority queue.
 */
void
insert_sublattice_pq ( sublattice_pq *pqueue,
                       mpz_t u,
                       mpz_t v,
                       mpz_t mod,
                       float val )
{
  /*
    gmp_fprintf (stderr, "# Debug: inserting (%Zd, %Zd), val: %.2f, "
    "used: %d, len: %d\n", u, v, val,
    pqueue->used, pqueue->len);
  */
  /* queue is full,  */
  if (pqueue->len == pqueue->used) {
    if ( val >= pqueue->val[1] ) {
      insert_sublattice_pq_down (pqueue, u, v, mod, val);
    }
  }

  /* queue is not full, sift-up */
  else if (pqueue->len > pqueue->used) {
    insert_sublattice_pq_up (pqueue, u, v, mod, val);
  }
  else {
    fprintf(stderr,"Error: error (pqueue->len < pqueue->used) "
            "in insert_sublattice_pq()\n");
    exit(1);
  }
}


/**
 * Create priority queue for sublattices over a single p^e.
 */
void
new_single_sublattice_pq ( single_sublattice_pq **ppqueue,
                           unsigned long len )
{
  if ( len < 1 ) {
    fprintf(stderr,"Error: len < 1 in new_single_sublattice_pq()\n");
    exit(1);
  }

  (*ppqueue) = (single_sublattice_pq *) malloc (
    sizeof(single_sublattice_pq) );
  if ( (*ppqueue) == NULL) {
    fprintf(stderr,"Error: malloc failed in new_single_sublattice_pq()\n");
    exit(1);
  }

  (*ppqueue)->len = len;
  (*ppqueue)->u = (unsigned int *) malloc (len * sizeof (unsigned int));
  (*ppqueue)->v = (unsigned int *) malloc (len * sizeof (unsigned int));
  (*ppqueue)->e = (char *) malloc (len * sizeof (char));
  (*ppqueue)->val = (float *) malloc (len * sizeof (float));

  if ( (*ppqueue)->u == NULL ||
       (*ppqueue)->v == NULL ||
       (*ppqueue)->e == NULL ||
       (*ppqueue)->val == NULL ) {
    fprintf(stderr,"Error: malloc failed in new_single_sublattice_pq()\n");
    exit(1);
  }

  /* u[0] and v[0] are null elements */
  (*ppqueue)->u[0] = 0;
  (*ppqueue)->v[0] = 0;
  (*ppqueue)->e[0] = 0;
  (*ppqueue)->val[0] = -FLT_MAX;
  (*ppqueue)->used = 1;
}


/**
 * Free
 */
void
free_single_sublattice_pq ( single_sublattice_pq **ppqueue )
{
  free ( (*ppqueue)->u );
  free ( (*ppqueue)->v );
  free ( (*ppqueue)->e );
  free ( (*ppqueue)->val );
  free ( *ppqueue );
}


/**
 * Sift-up to add, if the queue is not full.
 */
static inline void
insert_single_sublattice_pq_up ( single_sublattice_pq *pqueue,
                                 unsigned int u,
                                 unsigned int v,
                                 float val,
                                 char e )
{
  int k;

  for ( k = pqueue->used;
        val < pqueue->val[pq_parent(k)];
        k /= 2 )
  {
    pqueue->u[k] = pqueue->u[pq_parent(k)];
    pqueue->v[k] = pqueue->v[pq_parent(k)];
    pqueue->e[k] = pqueue->e[pq_parent(k)];
    pqueue->val[k] = pqueue->val[pq_parent(k)];
  }

  pqueue->u[k] = u;
  pqueue->v[k] = v;
  pqueue->e[k] = e;
  pqueue->val[k] = val;

  pqueue->used ++;
}


/**
 * Sift-down, if the heap is full.
 */
static inline void
insert_single_sublattice_pq_down ( single_sublattice_pq *pqueue,
                                   unsigned int u,
                                   unsigned int v,
                                   float val,
                                   char e )
{
  int k, l;

  for (k = 1; k*2 < pqueue->used; k = l) {

    l = (k << 1);

    /* right < left ? */
    if ( (l+1) < pqueue->used &&
         (pqueue->val[l+1] < pqueue->val[l]) )
      l ++;

    /* switch smaller child with parent */
    if ( pqueue->val[l] < val ) {

      pqueue->u[k] = pqueue->u[l];
      pqueue->v[k] = pqueue->v[l];
      pqueue->e[k] = pqueue->e[l];
      pqueue->val[k] = pqueue->val[l];
    }
    else
      break;
  }
  pqueue->u[k] = u;
  pqueue->v[k] = v;
  pqueue->e[k] = e;
  pqueue->val[k] = val;
}


/**
 * Insert to the priority queue.
 */
void
insert_single_sublattice_pq ( single_sublattice_pq *pqueue,
                              unsigned int u,
                              unsigned int v,
                              float val,
                              char e )
{

  /* queue is full,  */
  if (pqueue->len == pqueue->used) {
    if ( val > pqueue->val[1] ) {
      insert_single_sublattice_pq_down (pqueue, u, v, val, e);
    }
  }
  /* queue is not full, sift-up */
  else if (pqueue->len > pqueue->used) {
    insert_single_sublattice_pq_up (pqueue, u, v, val, e);
  }
  else {
    fprintf (stderr, "Error: error (pqueue->len < pqueue->used) "
             "in insert_single_sublattice_pq()\n");
    exit(1);
  }
}


/**
 * Extract the max of the priority queue.
 */
void
extract_single_sublattice_pq ( single_sublattice_pq *pqueue,
                               unsigned int *u,
                               unsigned int *v,
                               float *val,
                               char *e )
{
  /* don't extract u[0] since it is just a placeholder. */
  pqueue->used --;
  (*u) = pqueue->u[1];
  (*v) = pqueue->v[1];
  (*val) = pqueue->val[1];
  (*e) = pqueue->e[1];

  insert_single_sublattice_pq_down ( pqueue,
                                     pqueue->u[pqueue->used],
                                     pqueue->v[pqueue->used],
                                     pqueue->val[pqueue->used],
                                     pqueue->e[pqueue->used] );
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
                     mpz_t u,
                     mpz_t v,
                     mpz_t modulus,
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
                       mpz_t u,
                       mpz_t v,
                       mpz_t modulus,
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
                   mpz_t u,
                   mpz_t v,
                   mpz_t modulus,
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
                  mpz_t u,
                  mpz_t v,
                  mpz_t modulus,
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
                       mpz_t u,
                       mpz_t v,
                       mpz_t modulus,
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
                         mpz_t u,
                         mpz_t v,
                         mpz_t modulus,
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
                     mpz_t u,
                     mpz_t v,
                     mpz_t modulus,
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
                    mpz_t u,
                    mpz_t v,
                    mpz_t modulus,
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

