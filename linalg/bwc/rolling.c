#include "cado.h" // IWYU pragma: keep
#include <errno.h>      // for errno, ENOENT
#include <stdio.h>      // for printf, asprintf, size_t, NULL, perror, sscanf
#include <stdlib.h>     // for free, realloc, qsort
#include <string.h>     // for strerror
#include <time.h>       // for time

#include <unistd.h>     // for unlink
#include <sys/stat.h>   // for stat, st_mtime
#include <sys/types.h>  // for time_t
#include <dirent.h>     // for closedir, opendir, readdir, DIR, dirent

#include "rolling.h"
#include "bw-common.h"  // for bw
#include "macros.h"     // for ASSERT_ALWAYS
#include "portability.h" // asprintf // IWYU pragma: keep


typedef int (*sortfunc_t) (const void *, const void *);

int uint_cmp(const unsigned int * a, const unsigned int * b)
{
    return (*b < *a) - (*a < *b);
}

void keep_rolling_checkpoints(const char * stem, unsigned int v)
{
    if (bw->keep_rolling_checkpoints == 0)
        return;

    DIR * d;
    d = opendir(".");
    struct dirent * de;

    char * spat;
    int rc = asprintf(&spat, "%s.%%u", stem);
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
    if (svs <= (size_t) bw->keep_rolling_checkpoints) {
        free(vs);
        free(spat);
        return;
    }
    qsort(vs, svs, sizeof(unsigned int), (sortfunc_t) &uint_cmp);
    for(size_t i = 0 ; i < svs - bw->keep_rolling_checkpoints ; i++) {
        unsigned int k = vs[i];
        if (bw->checkpoint_precious && (k % bw->checkpoint_precious == 0))
            continue;
        if (k == 0)
            continue;
        char * v;
        rc = asprintf(&v, spat, k);
        ASSERT_ALWAYS(rc >= 0);
        struct stat sbuf[1];
        rc = stat(v, sbuf);
        if (rc < 0) {
            if (errno == ENOENT) {
                printf("Old checkpoint %s is gone already\n", v);
            } else {
                printf("Old checkpoint %s: %s\n", v, strerror(errno));
            }
        } else {
            ASSERT_ALWAYS(rc == 0);
            time_t now = time(NULL);
            int age = now - sbuf->st_mtime;
            if (age < bw->keep_checkpoints_younger_than) {
                printf("Not discarding old checkpoint %s, too recent (%d s < %d)\n", v, age, bw->keep_checkpoints_younger_than);
            } else {
                printf("Discarding old checkpoint %s\n", v);
                rc = unlink(v);
                if (rc < 0) perror(v);
            }
        }
        free(v);
    }
    free(spat);
    free(vs);
}

