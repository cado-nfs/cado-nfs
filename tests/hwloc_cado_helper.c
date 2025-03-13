#include "cado.h" // IWYU pragma: keep

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <hwloc.h>

void usage()
{
    fprintf(stderr,
            "Usage: hwloc-get-xml [-o [output filename] | --ncores ]\n");
    exit(EXIT_FAILURE);
}

int main(int argc, char const * argv[])
{
    if (argc < 2)
        usage();

    hwloc_topology_t top;
    hwloc_topology_init(&top);
    hwloc_topology_load(top);
    int const depth = hwloc_topology_get_depth(top);

    argv++, argc--;

    for( ; argc ; ) {
        if (argc >= 2 && strcmp(argv[0], "-o") == 0) {
            hwloc_topology_export_xml(top, argv[1], 0);
            argv++, argc--;
            argv++, argc--;
        } else if (argc >= 1 && strcmp(argv[0], "--ncores") == 0) {
            unsigned int const npu = hwloc_get_nbobjs_by_depth(top, depth-1);
            printf("%u\n", npu);
            argv++, argc--;
        } else {
            usage();
        }
    }

    hwloc_topology_destroy(top);
    return EXIT_SUCCESS;
}
