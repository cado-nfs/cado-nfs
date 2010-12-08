#include "cado.h"
#include <stdio.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

#include "bwc_config.h"
#include "debug.h"
#include "params.h"
#include "async.h"
#include "filenames.h"

#include "rolling.h"
#include "bw-common.h"

typedef int (*sortfunc_t) (const void *, const void *);

int uint_cmp(const unsigned int * a, const unsigned int * b)
{
    return (*b < *a) - (*a < *b);
}

void keep_rolling_checkpoints(balancing_ptr bal, const char * stem, unsigned int v)
{
    if (bw->keep_rolling_checkpoints == 0)
        return;

    DIR * d;
    d = opendir(".");
    struct dirent * de;

    char * spat;
    int rc = asprintf(&spat, COMMON_VECTOR_ITERATE_PARSE_PATTERN, stem, bal->h->checksum);
    ASSERT_ALWAYS(rc >= 0);

    unsigned int * vs = NULL;
    size_t svs = 0;
    size_t avs = 0;
    
    avs = 32;
    vs = realloc(vs, avs * sizeof(unsigned int));

    for( ; (de = readdir(d)) != NULL ; ) {
        unsigned int k;
        if (sscanf(de->d_name, spat, &k) != 1)
            continue;
        if (v && k > v)
            continue;
        if (svs >= avs) {
            avs += avs / 4;
            vs = realloc(vs, avs * sizeof(unsigned int));
        }
        vs[svs++] = k;
    }
    closedir(d);
    ASSERT_ALWAYS(svs);
    qsort(vs, svs, sizeof(unsigned int), (sortfunc_t) &uint_cmp);
    if (svs <= (size_t) bw->keep_rolling_checkpoints)
        return;
    for(size_t i = 0 ; i < svs - bw->keep_rolling_checkpoints ; i++) {
        unsigned int k = vs[i];
        if (bw->checkpoint_precious && (k % bw->checkpoint_precious == 0))
            continue;
        char * v;
        rc = asprintf(&v, spat, k);
        ASSERT_ALWAYS(rc >= 0);
        printf("Discarding old checkpoint %s\n", v);
        rc = unlink(v);
        if (rc < 0) {
            perror(v);
        }
        free(v);
    }
    free(spat);
}

