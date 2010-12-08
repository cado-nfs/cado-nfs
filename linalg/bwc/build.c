#include "cado.h"

#include <stdio.h>
#include "bwc_config.h"
#include "matmul.h"
#include "abase.h"
#include "macros.h"

void usage()
{
    fprintf(stderr, "Usage: ./build <file> [<impl>] [right|left]\n");
    exit(1);
}
int main(int argc, char * argv[])
{
    abobj_t xx MAYBE_UNUSED;
    abobj_init(xx);

    matmul_t mm;

    const char * impl = "basic";

    setvbuf(stdout,NULL,_IONBF,0);
    setvbuf(stderr,NULL,_IONBF,0);

    if (argc != 2 && argc != 3 && argc != 4) {
        usage();
    }

    if (argc >= 3)
        impl = argv[2];

    int d = 0;
    int nd = 2;

    if (argc == 4) {
        if (strcmp(argv[3],"left") == 0)  {
            d = 0;
        } else if(strcmp(argv[3],"right") == 0)  {
            d = 1;
        } else {
            usage();
        }
        nd = 1;
    }

    for( ; nd-- ; d ^= 1) {
        if (d == 1) {
            fprintf(stderr, "Saving cache for matrix-times-vector\n");
        } else {
            fprintf(stderr, "Saving cache for vector-times-matrix\n");
        }
        mm = matmul_init(xx, 0, 0, argv[1], argv[2], NULL, d);
        matmul_build_cache(mm, NULL);
        matmul_save_cache(mm);
        matmul_clear(mm);
    }
}
