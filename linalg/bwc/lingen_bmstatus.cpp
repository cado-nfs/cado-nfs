#include "cado.h" // IWYU pragma: keep

#include <climits>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <map>
#include <tuple>

#include "macros.h"
#include "lingen_bmstatus.hpp"
#include "lingen_call_companion.hpp"
#include "lingen_hints.hpp"
#include "select_mpi.h"

/* Attention: reloading a checkpoint invalidates this reference !! */
template<bool is_binary>
lingen_call_companion & bmstatus<is_binary>::companion(int depth, size_t L)/*{{{*/
{
    lingen_hints::key_type const K { depth, L };

    if (hints.find(K) != hints.end())
        return hints[K];

    int rank;
    MPI_Comm_rank(com[0], &rank);

    if (!rank)
        fprintf(stderr, "# No tuned configuration for"
                " depth=%d input_length=%zu\n",
                depth, L);

    struct same_depth {
        int d;
        bool operator()(lingen_hints::value_type const & kv) const {
            return kv.first.depth == d;
        }
    };

    auto it = std::find_if(hints.begin(), hints.end(), same_depth { depth } );
    if (it == hints.end()) {
        if (!rank)
            fprintf(stderr, "# No tuned configuration for"
                    " depth=%d !!!\n", depth);
        exit(EXIT_FAILURE);
    }
    if (!rank)
        fprintf(stderr, "# Using nearby configuration for"
                " depth=%d input_length=%zu\n",
                it->first.depth, it->first.L);

    hints[K]=it->second;

    ASSERT_ALWAYS(it->second.complete);

    return it->second;
}/*}}}*/


template<bool is_binary>
void bmstatus<is_binary>::display_deltas() const /*{{{*/
{
    unsigned int const m = d.m;
    unsigned int const n = d.n;

    int rank;
    MPI_Comm_rank(com[0], &rank);

    if (!rank) {
        /*
        printf("Final, t=%u: delta =", t);
        for(unsigned int j = 0; j < m + n; j++) {
            printf(" %u", delta[j]);
            if (lucky[j] < 0) {
                printf("(*)");
            }
        }
        printf("\n");
        */
        printf("Final, t=%4u: delta =", t);
        unsigned int last = UINT_MAX;
        unsigned int nrep = 0;
        int overflow = INT_MAX;
        for(unsigned int i = 0 ; i < m + n ; i++) {
            unsigned int const d = delta[i];
            if (d == last && (lucky[i] < 0) == overflow) {
                nrep++;
                continue;
            }
            // Flush the pending repeats
            if (last != UINT_MAX) {
                printf(" %u", last);
                if (overflow)
                    printf(" (*)");
                if (nrep > 1)
                    printf(" [%u]", nrep);
            }
            last = d;
            overflow = lucky[i] < 0;
            nrep = 1;
        }
        ASSERT_ALWAYS(last != UINT_MAX);
        printf(" %u", last);
        if (overflow)
            printf(" (*)");
        if (nrep > 1)
            printf(" [%u]", nrep);
        printf("\n");
    }
}/*}}}*/

template<bool is_binary>
std::tuple<unsigned int, unsigned int> bmstatus<is_binary>::get_minmax_delta_on_solutions() const /*{{{*/
{
    unsigned int maxdelta = 0;
    unsigned int mindelta = UINT_MAX;
    for(unsigned int j = 0 ; j < d.m + d.n ; j++) {
        if (lucky[j] <= 0) continue;
        maxdelta = std::max(maxdelta, delta[j]);
        mindelta = std::min(maxdelta, delta[j]);
    }
    return std::make_tuple(mindelta, maxdelta);
}/*}}}*/
template<bool is_binary>
unsigned int bmstatus<is_binary>::get_max_delta_on_solutions() const/*{{{*/
{
    unsigned int mindelta, maxdelta;
    std::tie(mindelta, maxdelta) = get_minmax_delta_on_solutions();
    return maxdelta;
}/*}}}*/

#ifdef LINGEN_BINARY
template struct bmstatus<true>;
#else
template struct bmstatus<false>;
#endif
