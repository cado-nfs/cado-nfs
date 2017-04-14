#include "cado.h"
#include <vector>
#include <queue>
#include <cstdint>
#include <limits>
#include <climits>

#include "minimum_spanning_tree.hpp"
#include "macros.h"

struct queue_item {
    uint16_t weight;
    uint8_t start;
    uint8_t end;
    bool operator<(queue_item const& o) const { return weight < o.weight; }
};

using namespace std;

/* naive implementation of Prim's algorithm.  */
pair<int, vector<pair<int, int>>>
minimalSpanningTreePrimNaive (vector<int> weights, size_t m)
{
    ASSERT_ALWAYS(weights.size() == m * m);
    int w = 0;
    vector<pair<int, int>> edges;

    /* S is the set of vertices in the current tree */
    vector<int> S(1,0);

    /* T is the set of remaining vertices */
    vector<int> T(m-1);
    for (size_t i = 0; i < m-1; i++)
        T[i] = i+1; /* T = {1, 2, ..., m-1} */

    for( ; !T.empty() ; T.pop_back()) {
        /* find the edge with minimal weight from S to T */
        int wmin = INT_MAX;
        int imin = -1;
        int jmin = -1;
        for (size_t i = 0; i < S.size(); i++)
            for (size_t j = 0; j < T.size(); j++)
                if (weights[S[i]*m + T[j]] < wmin) {
                    imin = i;
                    jmin = j;
                    wmin = weights[S[i]*m + T[j]];
                }
        ASSERT(imin != -1 && jmin != -1);
        w += wmin;
        S.push_back(T[jmin]);
        edges.push_back(make_pair(S[imin], T[jmin]));
        T[jmin] = T.back();
    }
    return make_pair(w, edges);
}

#if 0
/* Return the weight of the minimal spanning tree.
   For each (s,t) which is part of the tree, we have
   start[i] = s and end[i] = t for 0 <= i < m-1. */
static pair<int, vector<pair<int, int>>>
minimalSpanningTreeWithPrim (vector<int> weights, size_t m)
{
    /* the queue contains at most m-1 edges when nU=1 (those connected to the
       root node), then at most 2*(m-2) when nU=2 (we remove one node and add
       m-2), ... More generally when nU=k it is at most k*(m-k) + (k-1)*(k-2)/2.
       The maximum is for k=m-1 or m-2 when it equals m^2/2-3/2*m+2 < m^2/2. */
    index_t Q[MAX_QUEUE_SIZE];
    int u, s, t, i, nU, nV, w;
    int V[MERGE_LEVEL_MAX]; /* edges remaining to be dealt with */
    int index[MERGE_LEVEL_MAX];

    // nodes of T
    for(i = 0; i < m; i++){
        V[i] = i;
        index[i] = i; /* V[index[i]] = i if i is in V, -1 otherwise */
    }
    u = 0;
    index[u] = -1;
    index[V[m-1]] = u;
    V[u] = V[m-1];
    nU = 1;   /* number of edges already in the MST */
    nV = m-1; /* number of edges remaining to be dealt with */
    Q[0] = 0; /* heap is empty initially */
    addAllEdges (Q, u, V, nV, A);
#if DEBUG >= 1
    printQueue(Q);
#endif
    w = 0;
    while (!isQueueEmpty(Q)) /* while queue is non empty */
    {
        popQueue (&s, &t, Q); // pop queue
#if DEBUG >= 1
        fprintf(stderr, "Popping a = (%d, %d) of weight %d\n", s, t, A[s][t]);
#endif
        if (index[t] != -1){
            // t does not close a cycle
            // T[u] <- T[u] union (s, t)
#if DEBUG >= 1
            fprintf(stderr, "new edge: (%d, %d) of weight %d\n",s,t,A[s][t]);
#endif
            w += A[s][t];
            start[nU - 1] = s;
            end[nU - 1] = t;
            ASSERT(V[index[t]] == t);
            index[V[nV-1]] = index[t];
            V[index[t]] = V[nV-1];
            index[t] = -1;
            nV--;
            nU++;
            if (nV == 0)
                break;
            addAllEdges (Q, t, V, nV, A);
#if DEBUG >= 1
            printQueue(Q);
#endif
        }
    }
    return w;
}

    int
minimalSpanningTree(int *start, int *end,
        int m, int A[MERGE_LEVEL_MAX][MERGE_LEVEL_MAX])
{
    int ret;

    if (m <= 25)
        ret = minimalSpanningTreePrimNaive (start, end, A, m);
    else
        ret = minimalSpanningTreeWithPrim (start, end, A, m);
    return ret;
}

#endif
