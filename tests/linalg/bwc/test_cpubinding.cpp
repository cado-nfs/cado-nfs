#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cctype>
#include <cstring>
#include <cstdlib>

#include <hwloc.h>

#include "cpubinding.hpp"
#include "macros.h"
#include "params.h"
#include "misc.h"

/* Example quick use:
 * $build_tree/tests/linalg/bwc/test_cpubinding  -s "NUMANode:4 Socket:1 Core:12 PU:2" thr=8x8   cpubinding="NUMANode=>2x2 Core*12=>4x4"
 */

static void usage() {
    fprintf(stderr, "cpubinding example program\n"
            "Options:\n"
            "--input-topology-file <filename>     take <filename> as an hwloc hardware description\n"
            "--input-topology-string <string>       take <string> as an hwloc synthetic hardware description\n"
            "thr=<int>x<int>   give results for this target mapping\n"
            "cpubinding=<string>   use this mapping (or file)\n"
           );
    exit(1);
}

static int verbose = 0;

static int do_cpubinding_tests(const char * cpubinding_conf)
{
    FILE * f = fopen(cpubinding_conf, "r");
    if (f == nullptr) {
        fprintf(stderr, "%s: fopen failed\n", cpubinding_conf);
        exit(1);
    }

    char line[1024];

    // int nb_ok = 0;
    int nb_nok = 0;
    int idx = 0;
    for( ; fgets(line, sizeof(line), f) ; ) {
        int pos = 0;
        int want = 0;
        const char * find_magic2 =  "# EXPECT_FIND2";
        const char * find_magic =  "# EXPECT_FIND";
        const char * fail_magic =  "# EXPECT_FAIL";
        if (strncmp(line, find_magic2, strlen(find_magic2)) == 0) {
#if HWLOC_API_VERSION >= 0x020000
            pos += strlen(find_magic2);
            want = 1;
#else
            continue;
#endif
        } else if (strncmp(line, find_magic, strlen(find_magic)) == 0) {
            pos += strlen(find_magic);
            want = 1;
        } else if (strncmp(line, fail_magic, strlen(fail_magic)) == 0) {
            pos += strlen(fail_magic);
            want = 0;
        } else {
            continue;
        }

        int pos2;
        unsigned int t[2];
        int rc = sscanf(line + pos, "%u %u %n", &t[0], &t[1], &pos2);
        for(int n = strlen(line); n && isspace(line[n-1]); line[--n]='\0');
        ASSERT_ALWAYS(rc >= 2);

        if (verbose) printf("doing subtest %d: %s\n", idx, line);
        idx++;

        cxx_param_list pl2;
        param_list_add_key(pl2, "cpubinding", cpubinding_conf, PARAMETER_FROM_CMDLINE);
        param_list_add_key(pl2, "input-topology-string", line + pos + pos2, PARAMETER_FROM_CMDLINE);

        char * msg;
        void * cc = cpubinding_get_info(&msg, pl2, t[0], t[1]);
        /* don't make the tests too verbose */
        if (msg) {
            if (verbose) puts(msg);
            free(msg);
        }
        int ok = want == (cc != nullptr);
        if (verbose) printf("result: %s\n", ok_NOK(ok));
        // nb_ok += ok;
        nb_nok += !ok;

        cpubinding_do_pinning(cc, 0, 0);
        cpubinding_free_info(cc, t[0], t[1]);
    }
    fclose(f);
    return !nb_nok;
}

int main(int argc, char const * argv[])
{
    const char * cpubinding_conf = nullptr;
    cxx_param_list pl;

    argv++,argc--;
    param_list_configure_alias(pl, "input-topology-file", "-i");
    param_list_configure_alias(pl, "input-topology-string", "-s");
    param_list_configure_alias(pl, "cpubinding", "-c");
    param_list_configure_switch(pl, "-v", &verbose);
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Do perhaps some other things on the argument that haven't
         * been eaten at all. Like check whether it is a valid file to
         * source in order to get more options. See
         * param_list_read_stream and param_list_read_file for that. */
        fprintf (stderr, "Unknown option: %s\n", argv[0]);
        usage();
    }

    int seen_i = param_list_lookup_string(pl, "input-topology-file") != nullptr;
    int seen_s = param_list_lookup_string(pl, "input-topology-string") != nullptr;

    if (seen_i && seen_s) {
        fprintf(stderr, "Cannot have both -i and -s\n");
        exit(1);
    }

    cpubinding_conf = param_list_lookup_string(pl, "cpubinding");

    unsigned int thr[2] = {1,1};
    int parsed_thr = param_list_parse_uint_and_uint(pl, "thr", thr, "x");

    if (param_list_warn_unused(pl)) {
        usage();
    }

    char * msg;

    int rc = 0;
    if (parsed_thr) {
        /* This mode is here because this binary also serves as a quick
         * debug program, just to see how a given mapping string is
         * interpreted.
         */
        void * cc = cpubinding_get_info(&msg, pl, thr[0], thr[1]);
        if (msg) {
            puts(msg);
            free(msg);
        }
        cpubinding_do_pinning(cc, 0, 0);
        cpubinding_free_info(cc, thr[0], thr[1]);
    } else if (cpubinding_conf) {
        rc = do_cpubinding_tests(cpubinding_conf) ? EXIT_SUCCESS : EXIT_FAILURE;
    } else {
        fprintf(stderr, "don't know what to do !\n");
        exit(1);
    }

    return rc;
}
