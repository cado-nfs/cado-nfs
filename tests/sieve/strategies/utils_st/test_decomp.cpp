#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include "fmt/base.h"

#include "decomp.hpp"
#include "macros.h"
#include "tab_decomp.hpp"

int main()
{
    decomp el1 {1000, { 1U, 2U, 3U} };
    tabular_decomp t {el1};
    if (!(el1 == t[0])) {
        fprintf(stderr, "error with the test(1)!!!\n");
        return EXIT_FAILURE;
    }
    t.push_back(el1);
    // test realloc()
    const decomp el2 { 10000, { 11U, 10U, 9U, 8U, 7U }};
    t.push_back(el2);
    if (t[2] != el2) {
        fprintf(stderr, "error with the test(2)!!!\n");
        return EXIT_FAILURE;
    }
    if (t[2] == el1) {
        fprintf(stderr, "error with the test(3)!!!\n");
        return EXIT_FAILURE;
    }
    // set and get
    el1 = el2;
    if (el1 != el2) {
        fprintf(stderr, "error with the test(4)!!!\n");
        return EXIT_FAILURE;
    }
    // fprint fscan

    // coverity complains about insecure temp files. For tests, I don't
    // think it's a problem, really.
    // coverity[secure_temp]
    FILE * file = tmpfile();
    DIE_ERRNO_DIAG(file == nullptr, "tmpfile(%s)", "");
    fmt::print(file, "{}\n", t);

    /* rewind, and read again. This depends on some parsing code that I
     * haven't implemented yet. TODO
     */
#if 0
    fseek(file, 0, SEEK_SET);
    tabular_decomp_t * t2 = tabular_decomp_fscan(file);
    if (t2 == NULL) {
        fprintf(stderr, "read error on temp file\n");
        exit(EXIT_FAILURE);
    }
    fclose(file);

    if (tabular_decomp_are_equal(t2, t) != 1) {
        fprintf(stderr, "error with the test(5)!!!\n");
        return EXIT_FAILURE;
    }
    decomp_set_nb_elem(t2->tab[1], 21);
    if (tabular_decomp_are_equal(t2, t) != 0) {
        fprintf(stderr, "error with the test(6)!!!\n");
        return EXIT_FAILURE;
    }
    // free
#endif

    return EXIT_SUCCESS;
}
