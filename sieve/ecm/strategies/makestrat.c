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

#include "methods.h"


// We estimate that branches below a node that found a factor of size 45
// bits are enough, so that it is not necessary to optimize them.
#define MAX_FOUND 45

prior_ptr get_prior(int lo, int hi, int fbb);

// Main structure: the tree of strategies.

struct my_node {
    // data inside the node
    int cofac_range[2];
    unsigned int constraints; // tells whether we already played p-1...
    cofac_method_srcptr method;
    prior_srcptr prior;
    float acc_failure[60]; // that's the accumulated failure probabilities
                           // *before* running the method of this node.
    // childs
    struct my_node *failure_child;
    struct my_node **success_child; // end of array marked by NULL
};

typedef struct my_node strategy_node_struct;

typedef strategy_node_struct  strategy_node_t[1];
typedef strategy_node_struct *strategy_node_ptr;
typedef const strategy_node_struct *strategy_node_srcptr;

typedef strategy_node_t strategy_tree_t;


void strategy_tree_init(strategy_tree_t strategy_tree,
        int cofac_low, int cofac_hi, int fbb) {
    strategy_tree->constraints = 0;
    strategy_tree->cofac_range[0] = cofac_low;
    strategy_tree->cofac_range[1] = cofac_hi;
    strategy_tree->prior = get_prior(cofac_low, cofac_hi, fbb);
    strategy_tree->method = NULL;
    for (int i = 0; i < 60; ++i)
        strategy_tree->acc_failure[i] = 1;
    strategy_tree->failure_child = NULL;
    strategy_tree->success_child = NULL;
}

// forward declaration
strategy_node_ptr strategy_subtree_create(int cofac_range[2],
        prior_srcptr prior, float acc_failure[60], unsigned constraints);

// The creation at the root of the tree is syntactically different, but
// is essentially the same as for subtrees.
void strategy_tree_create(strategy_tree_t st_tree) {
    unsigned constraints = 0;
    strategy_node_ptr self = strategy_subtree_create(st_tree->cofac_range,
            st_tree->prior, st_tree->acc_failure, constraints);
    assert (self != NULL);
    st_tree->method = self->method;
    st_tree->failure_child = self->failure_child;
    st_tree->success_child = self->success_child;
    free(self);
}
    

// FIXME: do this properly, we need something more clever, here.
int is_valid_range(int cofac_range[2], prior_srcptr prior,
        float acc_fail[60])
{
    float prob;
    for (int i = 0; i < 40; ++i) {
        prob = prior->prob[i]*acc_fail[i]/(1-prior->prob[i]*(1-acc_fail[i]));
        if (prob > 0.025)
            return 1;
    }
    return 0;
}

strategy_node_ptr strategy_subtree_create(int cofac_range[2],
        prior_srcptr prior, float acc_failure[60], unsigned constraints)
{
//    printf("Entering subtree_create: %d-%d\n", cofac_range[0], cofac_range[1]);
    int fbb = prior->fbb;
    // easy case where there is no node.
    if (cofac_range[1] < 2*fbb)
        return NULL;
    // validate the node, otherwise return NULL
    // It takes into account the probability of getting a smooth number,
    // and early aborts if this is too small.
    if (!is_valid_range(cofac_range, prior, acc_failure))
        return NULL;
   
    strategy_node_ptr new_node = (strategy_node_ptr) malloc(
            sizeof(strategy_node_struct));
    ASSERT_ALWAYS(new_node != NULL);
    new_node->constraints = constraints;
    new_node->cofac_range[0] = cofac_range[0];
    new_node->cofac_range[1] = cofac_range[1];
    new_node->prior = prior;
    for (int i = 0; i < 60; ++i)
        new_node->acc_failure[i] = acc_failure[i];
    new_node->method = get_method_naive(new_node->cofac_range, new_node->prior,
            new_node->acc_failure, &constraints);
    
    float new_fail[60];
    for (int i = 0; i < 60; ++i)
        new_fail[i] = new_node->acc_failure[i]*(1-new_node->method->success[i]);

    // Create failure subtree
    new_node->failure_child = strategy_subtree_create(new_node->cofac_range,
            new_node->prior, new_fail, constraints);

    // Create success subtrees
    // For the moment, propagate the range size of the cofactor.
    int range_size = new_node->cofac_range[1] - new_node->cofac_range[0];
    // If cofactor size is less than 2*fbb, we are done.
    int cofac_limit = MAX(2*fbb, cofac_range[0]-MAX_FOUND);
    assert (range_size < fbb);
    int max_new_cofac = new_node->cofac_range[1] - fbb;
    int hi = new_node->cofac_range[0];
    while (hi >= max_new_cofac + range_size)
        hi -= range_size;
    int new_range[2] = {hi-range_size, hi};
    new_node->success_child = NULL;
    int nb_succ_child = 0;
    while (new_range[1] > cofac_limit) {
        prior_srcptr new_prior = get_prior(new_range[0], new_range[1], fbb);
        strategy_node_ptr child = strategy_subtree_create(new_range, new_prior,
                new_fail, constraints);
        if (child != NULL) {
            nb_succ_child++;
            new_node->success_child = (strategy_node_struct **) realloc(
                    new_node->success_child,
                    (1+nb_succ_child)*sizeof(strategy_node_struct *));
            ASSERT_ALWAYS(new_node->success_child != NULL);
            new_node->success_child[nb_succ_child-1] = child;
        }
        new_range[0] = new_range[0]-range_size;
        new_range[1] = new_range[1]-range_size;
    }
    if (nb_succ_child)
        new_node->success_child[nb_succ_child] = NULL;

    return new_node;
}

static int node_id=0;

void strategy_node_print_dot(strategy_node_srcptr node, const char * id)
{
    printf("%s [label = \"%s(%d,%d),(%d,%d,%d)\"];\n", id, id,
            node->cofac_range[0], node->cofac_range[1],
            node->method->type, node->method->B1, node->method->B2);
    
    char str2[1000];
    if (node->failure_child != NULL) {
        sprintf(str2, "%d", node_id++);
        printf("%s -> %s;\n", id, str2);
        strategy_node_print_dot(node->failure_child, str2);
    }
    if (node->success_child != NULL) {
        int i = 0;
        while (node->success_child[i] != NULL) {
            sprintf(str2, "%d", node_id++);
            printf("%s -> %s;\n", id, str2);
            strategy_node_print_dot(node->success_child[i], str2);
            ++i;
        }
    }
}

void strategy_tree_print_dot(strategy_node_srcptr st)
{
    printf("digraph G {\n");
    strategy_node_print_dot(st, "Init");
    printf("}\n");
}


int main(int argc, char **argv)
{
    strategy_tree_t st;
    strategy_tree_init(st, 150, 160, 16);

    strategy_tree_create(st);

    strategy_tree_print_dot(st);

}

