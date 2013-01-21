#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define ASSERT(x)       assert(x)
#define croak__(x,y) do {                                               \
        fprintf(stderr,"%s in %s at %s:%d -- %s\n",                     \
                (x),__func__,__FILE__,__LINE__,(y));                    \
    } while (0)
#define ASSERT_ALWAYS(x)                                                \
    do {                                                                \
        if (!(x)) {                                                     \
            croak__("code BUG() : condition " #x " failed",             \
                    "Abort");                                           \
            abort();                                                    \
        }                                                               \
    } while (0)

#ifndef MAX
#define MAX(h,i) ((h) > (i) ? (h) : (i))
#endif

#include "types.h"
#include "methods.h"


// The main structure that contains the strategy tree: this is just a dynamical
// array of nodes.
// We use a wrapping structure to simplify the memory management.

typedef struct {
    int n;
    int alloc;
    node_t *nodes;
} tree_struct;

typedef       tree_struct  tree_t[1];
typedef       tree_struct *tree_ptr;
typedef const tree_struct *tree_srcptr;

void init_tree(tree_ptr tree) {
    tree->alloc = 100;
    tree->nodes = (node_t *)malloc(tree->alloc*sizeof(node_t));
    ASSERT_ALWAYS(tree->nodes != NULL);
    tree->n = 0;
}

// Add a node to the tree. Return its index in the array of nodes.
int insert_node(tree_ptr tree, node_srcptr node) {
    if (tree->n == tree->alloc) {
        tree->alloc += tree->alloc/2;
        tree->nodes = (node_t *)realloc(tree->nodes, tree->alloc*sizeof(node_t));
        ASSERT_ALWAYS(tree->nodes != NULL);
    }
    memcpy(&tree->nodes[tree->n][0], node, sizeof(node_struct));
    tree->n++;
    return tree->n-1;
}

// Merge similar nodes according to some distance, but without creating
// loops.
void compress_tree(tree_ptr tree) {

}

// Return a set of contiguous cofac ranges.
// These are the more likely ranges to occur for the state+method of the
// given node. If everything is too unlikely, then return an empty list.
// nsizes is the set of ranges.
// sizes is an array of nsizes+1 integers, which give the bound of the
// ranges:
//   [ sizes[0]..sizes[1] [
//   [ sizes[1]..sizes[2] [
//   ...
//   [ sizes[nsizes-1]..sizes[nsizes] [
void most_likely_cofac_sizes(int *sizes, int *nsizes, node_srcptr node) {

}


// Create the subtree starting from the given nodes. At each
// step, choose a naive strategy.
void create_sub_tree_naive(tree_ptr tree, int node_idx) {
    node_ptr curr_node = tree->nodes[node_idx];
    // Choose the method
    cofac_method_srcptr method = get_method_naive(curr_node->state);
    
    // Create "leaf" childs
    curr_node->childs[0].type=CHILD_ABORT;
    curr_node->childs[0].node=NULL;
    curr_node->childs[1].type=CHILD_NOT_SMOOTH;
    curr_node->childs[1].node=NULL;
    curr_node->childs[2].type=CHILD_SMOOTH;
    curr_node->childs[2].node=NULL;
    
    // Create sub_tree childs corresponding to successes.
    node_t newnode;
    update_state(newnode->state, curr_node->state, method);
    int nsizes;
    int * sizes;
    most_likely_cofac_sizes(sizes, &nsizes, curr_node);
    curr_node->n_childs = nsizes+3;
    for (int i = 0; i < nsizes; ++i) {
        newnode->cofac_size[0] = sizes[i];
        newnode->cofac_size[1] = sizes[i+1];
        newnode->prior = get_prior(newnode->cofac_size);
        int new_idx = insert_node(tree, newnode);
        curr_node->childs[3+i].type = CHILD_CONTINUE;
        curr_node->childs[3+i].node = &tree->nodes[new_idx][0];
    }

    // Mark the probabilities as uncomputed.
    for (int i = 0; i < curr_node->n_childs; ++i)
        curr_node->childs[i].proba = -1.0; 
}



