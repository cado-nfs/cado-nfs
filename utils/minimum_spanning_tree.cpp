#include "cado.h"
#include <vector>
#include <queue>
#include <cstdint>
#include <limits>
#include <climits>

#include "minimum_spanning_tree.hpp"
#include "macros.h"

using namespace std;

typedef int weight_t;
typedef int vertex_t;

/* naive implementation of Prim's algorithm.  */
pair<weight_t, vector<pair<vertex_t, vertex_t>>>
minimum_spanning_tree_PrimNaive (matrix<weight_t> weights)
{
    weight_t w = 0;
    vertex_t m = weights.dimension(0);
    vector<pair<vertex_t, vertex_t>> edges;

    /* S is the set of vertices in the current tree */
    vector<vertex_t> S(1,0);

    /* T is the set of remaining vertices */
    vector<vertex_t> T(m-1);
    for (vertex_t i = 0; i < m-1; i++)
        T[i] = i+1; /* T = {1, 2, ..., m-1} */

    for( ; !T.empty() ; T.pop_back()) {
        /* find the edge with minimal weight from S to T */
        weight_t wmin = std::numeric_limits<weight_t>::max();
        vertex_t imin = -1;
        vertex_t jmin = -1;
        for (vertex_t i = 0; i < (vertex_t) S.size(); i++)
            for (vertex_t j = 0; j < (vertex_t) T.size(); j++)
                if (weights(S[i],T[j]) < wmin) {
                    imin = i;
                    jmin = j;
                    wmin = weights(S[i],T[j]);
                }
        if (wmin == std::numeric_limits<weight_t>::max()) {
            /* graph is disconnected. */
            edges.clear();
            return make_pair(wmin, edges);
        }
        ASSERT(imin != -1 && jmin != -1);
        w += wmin;
        S.push_back(T[jmin]);
        edges.push_back(make_pair(S[imin], T[jmin]));
        T[jmin] = T.back();
    }
    return make_pair(w, edges);
}

struct queue_item {
    weight_t weight;
    vertex_t start;
    vertex_t end;
    queue_item(weight_t weight, vertex_t start, vertex_t end) : weight(weight), start(start), end(end) {}
    bool operator<(queue_item const& o) const { return weight < o.weight; }
    bool operator>(queue_item const& o) const { return weight > o.weight; }
};

pair<weight_t, vector<pair<vertex_t, vertex_t>>>
minimum_spanning_tree_WithPrim(matrix<weight_t> weights)
{
    weight_t w = 0;
    vertex_t m = weights.dimension(0);
    std::priority_queue<queue_item, std::vector<queue_item>, std::greater<queue_item>> Q;

    vector<vertex_t> V(m);
    vector<vertex_t> index_in_V(m); /* index in V */

    for(vertex_t i = 0; i < m; i++)
        index_in_V[i] = V[i] = i;

    /* we maintain:
     *  - index_in_V[V[i]]==i for all i < V.size().
     *  - V[index_in_V[i]]==i for i unattached.
     * Attached vertices have
     * index_in_V[i]==attached.
     */

    vector<pair<vertex_t, vertex_t>> edges;

    const vertex_t attached = std::numeric_limits<vertex_t>::max();
    const weight_t no_edge = std::numeric_limits<weight_t>::max();

#define ATTACH_VERTEX(u) do {						\
    V[index_in_V[u]] = V.back(); 					\
    index_in_V[V.back()] = index_in_V[u]; 				\
    index_in_V[u] = attached;						\
    V.pop_back();							\
} while (0)

#define ENQUEUE_EDGES(u) do {						\
    /* push to Q the edges from u to the unattached vertices */	\
    for(auto const & v : V) {						\
        if (weights(u, v) != no_edge)					\
        Q.push(queue_item(weights(u, v), u, v));		        \
    }									\
} while (0)

    ATTACH_VERTEX(0);
    ENQUEUE_EDGES(0);

    /* V becomes empty sooner than Q */
    for( ; !V.empty() ; ) {
        if (Q.empty()) {
            /* graph is disconnected. */
            edges.clear();
            return make_pair(no_edge, edges);
        }
        auto edge = Q.top();
        Q.pop();
        /* edge.start is attached by construction. when edge.end is
         * attached too, we have a cycle, so we ditch it. */
        if (index_in_V[edge.end] == attached) continue;
        w += edge.weight;
        edges.push_back(make_pair(edge.start, edge.end));
        ATTACH_VERTEX(edge.end);
        ENQUEUE_EDGES(edge.end);
    }

#undef ATTACH_VERTEX
#undef ENQUEUE_EDGES

    return make_pair(w, edges);
}

pair<weight_t, vector<pair<vertex_t, vertex_t>>>
minimum_spanning_tree(matrix<weight_t> weights)
{
    vertex_t m = weights.dimension(0);

    if (m <= 25)
        return minimum_spanning_tree_PrimNaive (weights);
    else
        return minimum_spanning_tree_WithPrim (weights);
}
