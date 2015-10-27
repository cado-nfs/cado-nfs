#include "cado.h"

#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <stdio.h>
#include <gmp.h>
#include "bwc_config.h"
#include "cado_config.h"
#include "bw-common.h"
#include "params.h"
#include "filenames.h"
#include "utils.h"
#include "select_mpi.h"
#include "timing.h"
#include "portability.h"



struct bw_params bw[1];

const char * dirtext[] = { "left", "right" };

typedef int (*sortfunc_t) (const void*, const void*);

static int intcmp(const int * a, const int * b)
{
    return (*a>*b)-(*b>*a);
}

void bw_common_decl_usage(param_list pl)/*{{{*/
{
    /* We declare here the doc parameters which are parsed *in this file* !  */

    /* {{{ Parameters related to the problem we are solving */
    param_list_decl_usage(pl, "prime", "prime defining the field over which we work");
    param_list_decl_usage(pl, "nullspace", "whether we solve xM=0 (nullspace=left), or Mx=0 (nullspace=right). Default is left for p=2, right for p>2.");
    /* }}} */

    /* {{{ general BW parameters */
    param_list_decl_usage(pl, "mn", "set the block Wiedemann parameters m and n to the same given value");
    param_list_decl_usage(pl, "m", "set the block Wiedemann parameter m to this value");
    param_list_decl_usage(pl, "n", "set the block Wiedemann parameter n to this value");
    /* }}} */

    /* {{{ Parameters which are related to the interaction with the OS */
    param_list_decl_usage(pl, "wdir", "working directory, created if it does not exist. All file accesses are relative to this directory.");
    param_list_decl_usage(pl, "cfg", "path to a file containing more parameters to be interpreted.");
    param_list_decl_usage(pl, "seed", "set seed for all pseudo-random number generation");
    param_list_decl_usage(pl, "v", "More verbose output");
    /* }}} */

    /* {{{ Parameters related to checkpoints for krylov/mksol */
    param_list_decl_usage(pl, "interval", "frequency of the checkpoints within kkrylov/mksol");
    param_list_decl_usage(pl, "start", "start krylov or mksol at this checkpoint");
    param_list_decl_usage(pl, "end", "end krylov or mksol at this checkpoint");
    param_list_decl_usage(pl, "skip_online_checks", "skip consistency checks after each iteration. Use at your own risk.");
    param_list_decl_usage(pl, "keep_rolling_checkpoints", "keep only this number of checkpoints, and remove the others");
    param_list_decl_usage(pl, "keep_checkpoints_younger_than", "assuming keep_rolling_checkpoints is on, keep some checkpoints nevertheless if recent enough");
    param_list_decl_usage(pl, "checkpoint_precious", "assuming keep_rolling_checkpoints is on, never delete checkpoints of index which is a multiple of that number");
    param_list_decl_usage(pl, "yes_i_insist", "do what I say, even it seems stupid or dangerous");
    param_list_decl_usage(pl, "check_stops", "for secure, create check files for all these values of the interval. This enables checking more checkpoint files against eachother, once offline checking is functional.");
    /* }}} */

    /* {{{ Parameters related to multi-sequences */
    param_list_decl_usage(pl, "ys", "indicates which sequence(s) should be worked on. Syntax is <n0>..<n1> for working with vectors of indices i with n0<=i<n1, with of course 0<=n0<n1<=n.");
    param_list_decl_usage(pl, "nsolvecs", "for mksol and gather, produce this many independent solutions. Default is n");
    /* }}} */

    verbose_decl_usage(pl);
}
/*}}}*/

#if 0
const char * bw_common_usage_string()
{
    static char t[]=
        "Common options:\n"
        "\twdir=<path>\tchdir to <path> beforehand\n"
        "\tcfg=<file>\timport many settings from <file>\n"
        "\tm=<int>\tset bw->m blocking factor\n"
        "\tn=<int>\tset bw->n blocking factor\n"
        "\tmn=<int>\tset both bw->m and bw->n (exclusive with the two above)\n"
        "\tmpi=<int>x<int>\tset number of mpi jobs. Must agree with mpiexec\n"
        "\tthr=<int>x<int>\tset number of threads.\n"
        "\tcheckpoints=<bool>\tsave checkpoints.\n"
        "\tmatrix=<filename>\tset matrix\n"
        "\tinterval=<int>\tset checking bw->interval\n"
        "\tseed=<int>\tseed value for picking random numbers\n"
        "\tys=<int>..<int>\tcoordinate range for krylov/mksol task\n"
        ;
    return t;
}
#endif

void bw_common_parse_cmdline(struct bw_params * bw, param_list pl, int * p_argc, char *** p_argv)/*{{{*/
{
    bw->original_argc = *p_argc;
    bw->original_argv = *p_argv;
    bw->wct_base = wct_seconds();

    (*p_argv)++, (*p_argc)--;
    param_list_configure_switch(pl, "-v", &bw->verbose);
    for( ; (*p_argc) ; ) {
        if (param_list_update_cmdline(pl, p_argc, p_argv)) { continue; }
        if (strcmp((*p_argv)[0],"--") == 0) {
            (*p_argv)++, (*p_argc)--;
            break;
        }
        fprintf(stderr, "Unhandled parameter %s\n", (*p_argv)[0]);
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    if (bw->can_print) {
        param_list_print_command_line(stderr, pl);
        param_list_print_command_line(stdout, pl);
    }
}
/*}}}*/

void bw_common_interpret_parameters(struct bw_params * bw, param_list pl)/*{{{*/
{
    verbose_interpret_parameters(pl);

    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "wdir")) != NULL) {
        /* We now do mkdir -p beforehand on all jobs. Note that at the
         * point where the current function is being called, we're not
         * multithreaded yet */
        mkdir_with_parents(tmp, 1);
        if (chdir(tmp) < 0) {
            fprintf(stderr, "chdir(%s): %s\n", tmp, strerror(errno));
            exit(1);
        }
    }

    const char * cfg;

    if ((cfg = param_list_lookup_string(pl, "cfg"))) {
        param_list_read_file(pl, cfg);
    }


    param_list_parse_int(pl, "seed", &bw->seed);
    param_list_parse_int(pl, "interval", &bw->interval);
    param_list_parse_int_and_int(pl, "ys", bw->ys, "..");
    param_list_parse_int(pl, "start", &bw->start);
    param_list_parse_int(pl, "end", &bw->end);
    param_list_parse_int(pl, "skip_online_checks", &bw->skip_online_checks);
    param_list_parse_int(pl, "keep_rolling_checkpoints", &bw->keep_rolling_checkpoints);
    param_list_parse_int(pl, "keep_checkpoints_younger_than", &bw->keep_checkpoints_younger_than);
    param_list_parse_int(pl, "checkpoint_precious", &bw->checkpoint_precious);

    int yes_i_insist = 0;
    param_list_parse_int(pl, "yes_i_insist", &yes_i_insist);

    if (bw->skip_online_checks && bw->keep_rolling_checkpoints) {
        fprintf(stderr, "The combination of skip_online_checks and keep_rolling_checkpointsis a dangerous match.");
        if (!yes_i_insist) {
            printf("\n");
            exit(1);
        }
        fprintf(stderr, " Proceeding anyway\n");
    }

    mpz_init_set_ui(bw->p, 2);
    param_list_parse_mpz(pl, "prime", bw->p);
    int nullspace_forced = 0;

    if ((tmp = param_list_lookup_string(pl, "nullspace")) != NULL) {
        char * tmp_l = strdup(tmp);
        for(unsigned int i = 0 ; i < strlen(tmp_l) ; i++) {
            char c = tmp[i];
            char cl = tolower(c);
            tmp_l[i] = cl;
            nullspace_forced |= c != cl;
        }
        if (strcmp(tmp_l, dirtext[0]) == 0) {
            bw->dir = 0;
        } else if (strcmp(tmp_l, dirtext[1]) == 0) {
            bw->dir = 1;
        } else {
            fprintf(stderr, "Parameter nullspace may only be %s|%s\n",
                    dirtext[0], dirtext[1]);
            exit(1);
        }
        free(tmp_l);
    } else {
        /* Default is right nullspace for p>2, and left for p==2 */
        bw->dir = mpz_cmp_ui(bw->p, 2) != 0;
        param_list_add_key(pl, "nullspace", dirtext[bw->dir], PARAMETER_FROM_FILE);
    }

    if ((mpz_cmp_ui(bw->p, 2) == 0) != (bw->dir == 0)) {
        if (mpz_cmp_ui(bw->p, 2) == 0) {
            fprintf(stderr, "p==2 seems appropriate for factoring. Yet, the nullspace parameter has been passed as nullspace=right. This looks odd.\n");
            if (nullspace_forced) {
                fprintf(stderr, "Proceeding anyway (uppercase nullspace argument)\n");
            } else {
                fprintf(stderr, "Aborting. Pass nullspace=LEFT if this is really intended.\n");
                exit(1);
            }
        } else {
            fprintf(stderr, "p>2 seems appropriate for discrete logarithm. Yet, the nullspace parameter has been passed as nullspace=left. This looks odd.\n");
            if (nullspace_forced) {
                fprintf(stderr, "Proceeding anyway (uppercase nullspace argument)\n");
            } else {
                fprintf(stderr, "Aborting. Pass nullspace=RIGHT if this is really intended.\n");
                exit(1);
            }
        }
    }

    
    int okm=0, okn=0;
    int mn;
    if (param_list_parse_int(pl, "mn", &mn)) {
        bw->m=mn;
        bw->n=mn;
        okm++;
        okn++;
    }
    okm += param_list_parse_int(pl, "m", &bw->m);
    okn += param_list_parse_int(pl, "n", &bw->n);
    if (!okm || !okn) {
        fprintf(stderr, "parameter m and/or n is missing\n");
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }
    bw->nsolvecs = bw->n;
    param_list_parse_int(pl, "nsolvecs", &bw->nsolvecs);

    bw->number_of_check_stops = param_list_parse_int_list(pl, "check_stops", bw->check_stops, MAX_NUMBER_OF_CHECK_STOPS, ",");
    int interval_already_in_check_stops = 0;
    for(int i = 0 ; i < bw->number_of_check_stops ; i++) {
        if (bw->check_stops[i] == bw->interval) {
            interval_already_in_check_stops = 1;
            break;
        }
    }
    if (!interval_already_in_check_stops) {
        ASSERT_ALWAYS(bw->number_of_check_stops < MAX_NUMBER_OF_CHECK_STOPS - 1);
        bw->check_stops[bw->number_of_check_stops++] = bw->interval;
    }
    qsort(bw->check_stops, bw->number_of_check_stops, sizeof(int), (sortfunc_t) &intcmp);

    if (bw->verbose && bw->can_print)
        param_list_display (pl, stderr);
}
/*}}}*/

static int bw_common_init_defaults(struct bw_params * bw)/*{{{*/
{
    /*** defaults ***/
    memset(bw, 0, sizeof(*bw));
    bw->interval = 1000;
    bw->can_print = 1;
    bw->ys[0] = bw->ys[1] = -1;
    bw->dir = 1;

    return 0;
}
/*}}}*/

int bw_common_init_new(struct bw_params * bw, int * p_argc, char *** p_argv)/*{{{*/
{
    /* First do MPI_Init */
#if defined(MPI_LIBRARY_MT_CAPABLE)
    int req = MPI_THREAD_MULTIPLE;
    int prov;
    MPI_Init_thread(p_argc, p_argv, req, &prov);
    if (req != prov) {
        fprintf(stderr, "Cannot init mpi with MPI_THREAD_MULTIPLE ;"
                " got %d != req %d\n",
                prov, req);
        exit(1);
    }
#if 0   /* was: elif OMPI_VERSION_ATLEAST(1,8,2) */
    /* This is just a try, right. In practice, we do rely on the
     * SERIALIZED model, so let's at least do as we cared about telling
     * the MPI implementation about it.
     */
    int req = MPI_THREAD_SERIALIZED;
    int prov;
    MPI_Init_thread(p_argc, p_argv, req, &prov);
    if (req != prov) {
        fprintf(stderr, "Cannot init mpi with MPI_THREAD_SERIALIZED ;"
                " got %d != req %d\n",
                prov, req);
        exit(1);
    }
#endif  /* #if 0 */
#else
    MPI_Init(p_argc, p_argv);
#endif
    int rank;
    int size;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    bw_common_init_defaults(bw);

    bw->can_print = rank == 0 || getenv("CAN_PRINT");

    setvbuf(stdout,NULL,_IONBF,0);
    setvbuf(stderr,NULL,_IONBF,0);

    return 0;
}
/*}}}*/
int bw_common_clear_new(struct bw_params * bw)/*{{{*/
{
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double wct = wct_seconds() - bw->wct_base;
    double cpu = seconds();
    char * ptr = strrchr(bw->original_argv[0], '/');
    if (ptr) {
        ptr++;
    } else {
        ptr = bw->original_argv[0];
    }
    MPI_Allreduce(MPI_IN_PLACE, &cpu, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    if (bw->can_print) {
        /* valgrind has a tendency to complain about this code depending
         * on unitialized data in the variable "cpu". This is most
         * probably due to MPI_Allreduce, and there's not much we can do,
         * unfortunately.
         */
        printf("Timings for %s: wct=%.2f cpu=%.2f (aggregated over all threads and %d MPI jobs)\n",
                ptr,
                wct, cpu,
                size);
    }
    mpz_clear(bw->p);
    MPI_Finalize();
    return 0;
}/*}}}*/



int get_rhs_file_header_stream(FILE * f, uint32_t * p_nrows, unsigned int * p_nrhs, mpz_ptr p_p)
{
    int rc;
    if (p_nrows) {
        rc = fscanf(f, "%" SCNu32, p_nrows);
        ASSERT_ALWAYS(rc == 1);
    } else {
        rc = fscanf(f, "%*" SCNu32);
        ASSERT_ALWAYS(rc == 0);
    }
    if (p_nrhs) {
        rc = fscanf(f, "%d", p_nrhs);
        ASSERT_ALWAYS(rc == 1);
    } else {
        rc = fscanf(f, "%*d");
        ASSERT_ALWAYS(rc == 0);
    }
    if (p_p) {
        rc = gmp_fscanf(f, "%Zd", p_p);
        ASSERT_ALWAYS(rc == 1);
    } else {
        rc = gmp_fscanf(f, "%*Zd");
        ASSERT_ALWAYS(rc == 0);
    }
    return 1;
}

int get_rhs_file_header(const char * filename, uint32_t * p_nrows, unsigned int * p_nrhs, mpz_ptr p_p)
{
    FILE * f = NULL;
    f = fopen(filename, "r");
    ASSERT_ALWAYS(f != NULL);
    int rc = get_rhs_file_header_stream(f, p_nrows, p_nrhs, p_p);
    fclose(f);
    return rc;
}

