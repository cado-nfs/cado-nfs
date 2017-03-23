#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include "renumber.h"


int main(int argc, char * argv[])
{
    const char * tablefile = NULL;
    const char * polyfile = NULL;
    argc--,argv++;
    for( ; argc ; argc--,argv++) {
        if (strcmp(argv[0], "--table") == 0) {
            argc--,argv++;
            tablefile = argv[0];
        } else if (strcmp(argv[0], "--poly") == 0) {
            argc--,argv++;
            polyfile = argv[0];
        } else {
            fprintf(stderr, "Unexpected argument %s\n", argv[0]);
            exit(EXIT_FAILURE);
        }
    }
    if (!polyfile || !tablefile) {
        fprintf(stderr, "Need both a polynomial and a renumbering table\n");
        exit(EXIT_FAILURE);
    }

    cado_poly cpoly;
    cado_poly_init(cpoly);
    cado_poly_read(cpoly, polyfile);

    renumber_t tab;
    renumber_init_for_reading(tab);
    renumber_read_table(tab, tablefile);

    renumber_iterator it;
    renumber_iterator_init(it, tab, cpoly);

    for( ; !renumber_iterator_done(it) ; renumber_iterator_next(it)) {
        if (it->side >= 0) {
            if (it->p == 0) {
                printf("%" PRid ": added column\n", it->i);
            } else if (cpoly->pols[it->side]->deg == 1){
                /* don't print r for rational side */
                printf("%" PRid ": %d %" PRpr "\n", it->i, it->side, it->p);
            } else {
                printf("%" PRid ": %d %" PRpr " %" PRpr "\n", it->i, it->side, it->p, it->r);
            }
        } else {
            /* This is for bad ideals only */
            printf("%" PRid ": %d %" PRpr " %" PRpr " [bad ideal]\n", it->i, -1-it->side, it->p, it->r);
        }
    }


    renumber_iterator_clear(it);
    renumber_clear(tab);
    cado_poly_clear(cpoly);

    return 0;
}
