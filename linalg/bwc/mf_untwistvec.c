#include "cado.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "balancing.h"
#include "utils.h"

int main(int argc, char * argv[])
{
    balancing bal;
    param_list pl;
    param_list_init(pl);
    argv++, argc--;
    // specifying nullspace-left or nullspace-right has no impact here,
    // as we have only one permutation to choose from (we're still
    // forcing the conjugated permnutations case). However, future
    // evolutions might require this.
    int nullspace_left = 0;
    param_list_configure_knob(pl, "--nullspace-left", &nullspace_left);
    unsigned int wild = 0;
    const char * bname = NULL;
    const char * vname = NULL;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        if (argv[0][0] != '-' && wild == 0) {
            bname = argv[0];
            wild++;
            argv++, argc--;
            continue;
        }
        if (argv[0][0] != '-' && wild == 1) {
            vname = argv[0];
            wild++;
            argv++, argc--;
            continue;
        }
    }
    balancing_read(bal, bname);
    FILE * v;
    v = fopen(vname, "r");
    if (!v) {
        perror(vname);
        exit(1);
    }
    struct stat sbuf[1];
    if (stat(vname, sbuf) < 0) {
        perror(vname);
        exit(1);
    }
    unsigned long expected = nullspace_left ? bal->trows : bal->tcols;
    assert(sbuf->st_size % expected == 0);
    size_t chunk = sbuf->st_size / expected;
    char * area = malloc(sbuf->st_size);
    memset(area, 0, sbuf->st_size);
    char * tmp = malloc(chunk);
    memset(tmp, 0, chunk);
    // in reality, it's the same permutations as long as only conjugated
    // permutations are implemented.
    uint32_t * perm = nullspace_left ? bal->rowperm : bal->colperm;
    for(unsigned long k = 0 ; k < expected ; k++) {
        if (fread(tmp, 1, chunk, v) < chunk) {
            fprintf(stderr, "short read\n");
            exit(1);
        }
        memcpy(area + perm[k] * chunk, tmp, chunk);
    }
    fclose(v);
    free(tmp);
    {
        const char * tmp;
        if ((tmp = param_list_lookup_string(pl, "out")) != NULL) {
            FILE * o;
            o= fopen(tmp, "w");
            fwrite(area, chunk, expected, o);
            fclose(o);
        } else {
            fwrite(area, chunk, expected, stdout);
        }
    }
    free(area);
    return 0;
}

