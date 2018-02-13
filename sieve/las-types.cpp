
#include "cado.h"
#include <stdio.h>
#include <stdarg.h>
#include <gmp.h>
#include "las-types.hpp"
#include "las-config.h"
#include "las-norms.hpp"
/* sieve_info stuff */

/* This function creates a new sieve_info structure, taking advantage of
 * structures which might already exist in las.sievers
 *  - for sieving, if factor base parameters are similar, we're going to
 *    share_factor_bases()
 *  - for cofactoring, if large prime bounds and mfbs are similar, we're
 *    going to reuse the strategies.
 * 
 * The siever_config structure to be passed to this function is not
 * permitted to lack anything.
 *
 * This function differs from las_info::get_sieve_info_from_config(),
 * since the latters also registers the returned object within
 * las.sievers (while the function here only *reads* this list).
 */
sieve_info::sieve_info(siever_config const & sc, cado_poly_srcptr cpoly, std::list<sieve_info> & sievers, cxx_param_list & pl) /*{{{*/
    : cpoly(cpoly), conf(sc)
{
    I = 1UL << sc.logI_adjusted;

    std::list<sieve_info>::iterator psi;
    
    /*** Sieving ***/

    psi = find_if(sievers.begin(), sievers.end(), sc.same_fb_parameters());

    if (psi != sievers.end()) {
        sieve_info & other(*psi);
        verbose_output_print(0, 1, "# copy factor base data from previous siever\n");
        share_factor_bases(other);
    } else {
        verbose_output_print(0, 1, "# bucket_region = %zu\n",
                BUCKET_REGION);
        init_factor_bases(pl);
        for (int side = 0; side < 2; side++) {
            print_fb_statistics(side);
        }
    }

    // Now that fb have been initialized, we can set the toplevel.
    toplevel = -1;
    for(int side = 0 ; side < 2 ; side++) {
        if (!sides[side].fb) continue;
        int level = sides[side].fb->get_toplevel();
        if (level > toplevel) toplevel = level;
    }

#if 0
    /* Initialize the number of buckets */
    /* (it's now done in sieve_info::update, which is more timely) */
    /* set the maximal value of the number of buckets. This directly
     * depends on A */
    uint32_t XX[FB_MAX_PARTS] = { 0, NB_BUCKETS_2, NB_BUCKETS_3, 0};
    uint64_t BRS[FB_MAX_PARTS] = BUCKET_REGIONS;
    uint64_t A = UINT64_C(1) << conf.logA;
    XX[toplevel] = iceildiv(A, BRS[toplevel]);
    for (int i = toplevel+1; i < FB_MAX_PARTS; ++i)
        XX[i] = 0;
    if (toplevel > 1 && XX[toplevel] == 1) {
        XX[toplevel-1] = iceildiv(A, BRS[toplevel-1]);
        ASSERT_ALWAYS(XX[toplevel-1] != 1);
    }
    for (int i = 0; i < FB_MAX_PARTS; ++i) {
        nb_buckets[i] = XX[i];
    }
#endif

    for(int side = 0 ; side < 2 ; side++) {
	init_trialdiv(side); /* Init refactoring stuff */
        if (!sides[side].fb) continue;
        init_fb_smallsieved(side);
        verbose_output_print(0, 2, "# small side-%d factor base", side);
        size_t nr_roots;
        sides[side].fb->get_part(0)->count_entries(NULL, &nr_roots, NULL);
        verbose_output_print(0, 2, " (total %zu)\n", nr_roots);
    }

    /*** Cofactoring ***/

    psi = find_if(sievers.begin(), sievers.end(), sc.same_cofactoring());

    if (psi != sievers.end()) {
        sieve_info & other(*psi);
        verbose_output_print(0, 1, "# copy cofactoring strategies from previous siever\n");
        strategies = other.strategies;
    } else {
        init_strategies(pl);
    }
}
/*}}}*/
/* This ctor is a simplified version of the former. Of course I'd be
 * happy to use delegated constructors here.
 */
sieve_info::sieve_info(siever_config const & sc, cado_poly_srcptr cpoly, cxx_param_list & pl) /*{{{*/
    : cpoly(cpoly), conf(sc)
{
    I = 1UL << sc.logI_adjusted;
    std::list<sieve_info>::iterator psi;
    verbose_output_print(0, 1, "# bucket_region = %zu\n", BUCKET_REGION);
    init_factor_bases(pl);
    for (int side = 0; side < 2; side++) {
        print_fb_statistics(side);
    }
    // Now that fb have been initialized, we can set the toplevel.
    toplevel = -1;
    for(int side = 0 ; side < 2 ; side++) {
        if (!sides[side].fb) continue;
        int level = sides[side].fb->get_toplevel();
        if (level > toplevel) toplevel = level;
    }

    for(int side = 0 ; side < 2 ; side++) {
        init_trialdiv(side); /* Init refactoring stuff */
        init_fb_smallsieved(side);
        verbose_output_print(0, 2, "# small side-%d factor base", side);
        size_t nr_roots;
        sides[side].fb->get_part(0)->count_entries(NULL, &nr_roots, NULL);
        verbose_output_print(0, 2, " (total %zu)\n", nr_roots);
    }
    init_strategies(pl);
}
/*}}}*/

void sieve_info::update_norm_data ()/*{{{*/
{
    for(int side = 0 ; side < 2 ; side++) {
        sides[side].lognorms = std::make_shared<lognorm_smart>(conf, cpoly, side, qbasis, J);
    }
}

/*}}}*/
void sieve_info::update (size_t nr_workspaces)/*{{{*/
{
    uint64_t A = UINT64_C(1) << conf.logA;

    /* update number of buckets at toplevel */
    size_t (&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;
    nb_buckets[toplevel] = iceildiv(A, BRS[toplevel]);

    // maybe there is only 1 bucket at toplevel and less than 256 at
    // toplevel-1, due to a tiny J.
    if (toplevel > 1) {
        if (nb_buckets[toplevel] == 1) {
            nb_buckets[toplevel-1] = iceildiv(A, BRS[toplevel - 1]);
            // we forbid skipping two levels.
            ASSERT_ALWAYS(nb_buckets[toplevel-1] != 1);
        } else {
            nb_buckets[toplevel-1] = BRS[toplevel]/BRS[toplevel-1];
        }
    }

    /* Update the slices of the factor base according to new log base */
    for(int side = 0 ; side < 2 ; side++) {
        /* The safety factor controls by how much a single slice should fill a
           bucket array at most, i.e., with .5, a single slice should never fill
           a bucket array more than half-way. */
        const double safety_factor = .5;
        sieve_info::side_info & sis(sides[side]);
        double max_weight[FB_MAX_PARTS];
        for (int i_part = 0; i_part < FB_MAX_PARTS; i_part++) {
            max_weight[i_part] = sis.max_bucket_fill_ratio[i_part] / nr_workspaces
                * safety_factor;
        }
        if (!sis.fb) continue;
        sis.fb->make_slices(sis.lognorms->scale, max_weight);
    }
}/*}}}*/


static void las_verbose_enter(cxx_param_list & pl, FILE * output, int verbose)
{
    verbose_interpret_parameters(pl);
    verbose_output_init(NR_CHANNELS);
    verbose_output_add(0, output, verbose + 1);
    verbose_output_add(1, stderr, 1);
    /* Channel 2 is for statistics. We always print them to las' normal output */
    verbose_output_add(2, output, 1);
    if (param_list_parse_switch(pl, "-stats-stderr")) {
        /* If we should also print stats to stderr, add stderr to channel 2 */
        verbose_output_add(2, stderr, 1);
    }
#ifdef TRACE_K
    const char *trace_file_name = param_list_lookup_string(pl, "traceout");
    FILE *trace_file = stderr;
    if (trace_file_name != NULL) {
        trace_file = fopen(trace_file_name, "w");
        DIE_ERRNO_DIAG(trace_file == NULL, "fopen", trace_file_name);
    }
    verbose_output_add(TRACE_CHANNEL, trace_file, 1);
#endif
}

static void las_verbose_leave()
{
    verbose_output_clear();
}

las_augmented_output_channel::las_augmented_output_channel(cxx_param_list & pl)
{
    output = stdout;
    outputname = param_list_lookup_string(pl, "out");
    if (outputname) {
	if (!(output = fopen_maybe_compressed(outputname, "w"))) {
	    fprintf(stderr, "Could not open %s for writing\n", outputname);
	    exit(EXIT_FAILURE);
	}
    }
    verbose = param_list_parse_switch(pl, "-v");
    setvbuf(output, NULL, _IOLBF, 0);      /* mingw has no setlinebuf */
    las_verbose_enter(pl, output, verbose);
}
las_augmented_output_channel::~las_augmented_output_channel()
{
    if (outputname)
        fclose_maybe_compressed(output, outputname);
    las_verbose_leave();
}


/* las_info stuff */

las_info::las_info(cxx_param_list & pl)
    : las_augmented_output_channel(pl),
        config_pool(pl)
#ifdef  DLP_DESCENT
      , dlog_base(pl)
#endif
      /*{{{*/
{
    /* We strive to initialize things in the exact order they're written
     * in the struct */
    // ----- general operational flags {{{
    nb_threads = 1;		/* default value */
    param_list_parse_int(pl, "t", &nb_threads);
    if (nb_threads <= 0) {
	fprintf(stderr,
		"Error, please provide a positive number of threads\n");
	exit(EXIT_FAILURE);
    }

    galois = param_list_lookup_string(pl, "galois");
    suppress_duplicates = param_list_parse_switch(pl, "-dup");

    param_list_print_command_line(output, pl);


    las_display_config_flags();
    /*  Parse polynomial */
    cado_poly_init(cpoly);
    const char *tmp;
    if ((tmp = param_list_lookup_string(pl, "poly")) == NULL) {
        fprintf(stderr, "Error: -poly is missing\n");
        param_list_print_usage(pl, NULL, stderr);
	cado_poly_clear(cpoly);
        exit(EXIT_FAILURE);
    }
    if (!cado_poly_read(cpoly, tmp)) {
	fprintf(stderr, "Error reading polynomial file %s\n", tmp);
	cado_poly_clear(cpoly);
	exit(EXIT_FAILURE);
    }
    // sc.skewness = cpoly->skew;
    /* -skew (or -S) may override (or set) the skewness given in the
     * polynomial file */
    param_list_parse_double(pl, "skew", &(cpoly->skew));
    if (cpoly->skew <= 0.0) {
	fprintf(stderr, "Error, please provide a positive skewness\n");
	cado_poly_clear(cpoly);
	exit(EXIT_FAILURE);
    }
    gmp_randinit_default(rstate);
    unsigned long seed = 0;
    if (param_list_parse_ulong(pl, "seed", &seed))
        gmp_randseed_ui(rstate, seed);

    if (const char * tmp = param_list_lookup_string(pl, "bkmult")) {
        bk_multiplier = bkmult_specifier(tmp);
    }

    // }}}

    // ----- stuff roughly related to the descent {{{
#ifdef  DLP_DESCENT
    descent_helper = NULL;
#endif
    // }}}

    // ----- todo list and such {{{
    nq_pushed = 0;
    nq_max = UINT_MAX;
    random_sampling = 0;
    if (param_list_parse_uint(pl, "random-sample", &nq_max)) {
        random_sampling = 1;
    } else if (param_list_parse_uint(pl, "nq", &nq_max)) {
        if (param_list_lookup_string(pl, "q1") || param_list_lookup_string(pl, "rho")) {
            fprintf(stderr, "Error: argument -nq is incompatible with -q1 or -rho\n");
            exit(EXIT_FAILURE);
        }
    }

    /* Init and parse info regarding work to be done by the siever */
    /* Actual parsing of the command-line fragments is done within
     * las_todo_feed, but this is an admittedly contrived way to work */
    const char * filename = param_list_lookup_string(pl, "todo");
    if (filename) {
        todo_list_fd = fopen(filename, "r");
        if (todo_list_fd == NULL) {
            fprintf(stderr, "%s: %s\n", filename, strerror(errno));
            /* There's no point in proceeding, since it would really change
             * the behaviour of the program to do so */
            cado_poly_clear(cpoly);
            exit(EXIT_FAILURE);
        }
    } else {
        todo_list_fd = NULL;
    }

    /* composite special-q ? */
    allow_composite_q = param_list_parse_switch(pl, "-allow-compsq");
    if (allow_composite_q) {
        if (galois) {
            fprintf(stderr, "-galois and -allow-compsq are incompatible options at the moment");
            exit(EXIT_FAILURE);
        }
        if (!param_list_parse_uint64(pl, "qfac-min", &qfac_min)) {
            qfac_min = 1024;
        }
        if (!param_list_parse_uint64(pl, "qfac-max", &qfac_max)) {
            qfac_max = UINT64_MAX;
        }
    }

    // }}}

    // ----- batch mode {{{
    batch = param_list_parse_switch(pl, "-batch");
    batch_print_survivors = param_list_parse_switch(pl, "-batch-print-survivors");
    cofac_list_init (L);
    // }}} 

    dump_filename = param_list_lookup_string(pl, "dumpfile");

    init_cof_stats(pl);
}/*}}}*/

las_info::~las_info()/*{{{*/
{
    // ----- general operational flags {{{
    gmp_randclear(rstate);
    cado_poly_clear(cpoly);
    // }}}

    // ----- todo list and such {{{
    if (todo_list_fd) {
        fclose(todo_list_fd);
        todo_list_fd = NULL;
    }
    // }}}
 
    // ----- batch mode: very little
    cofac_list_clear (L);

    clear_cof_stats();
}/*}}}*/

// {{{ las_info::{init,clear,print}_cof_stats
void las_info::init_cof_stats(param_list_ptr pl)
{
    const char *statsfilename = param_list_lookup_string (pl, "stats-cofact");
    if (statsfilename != NULL) { /* a file was given */
        if (config_pool.default_config_ptr == NULL) {
            fprintf(stderr, "Error: option stats-cofact works only "
                    "with a default config\n");
            exit(EXIT_FAILURE);
#ifdef DLP_DESCENT
        } else if (param_list_lookup_string(pl, "descent-hint-table")) {
            verbose_output_print(0, 1, "# Warning: option stats-cofact "
                    "only applies to the default siever config\n");
#endif
        }
        siever_config const & sc(*config_pool.default_config_ptr);

        cof_stats_file = fopen (statsfilename, "w");
        if (cof_stats_file == NULL)
        {
            fprintf (stderr, "Error, cannot create file %s\n",
                    statsfilename);
            exit (EXIT_FAILURE);
        }
        //allocate cof_call and cof_succ
        int mfb0 = sc.sides[0].mfb;
        int mfb1 = sc.sides[1].mfb;
        cof_call = (uint32_t**) malloc ((mfb0+1) * sizeof(uint32_t*));
        cof_succ = (uint32_t**) malloc ((mfb0+1) * sizeof(uint32_t*));
        for (int i = 0; i <= mfb0; i++) {
            cof_call[i] = (uint32_t*) malloc ((mfb1+1) * sizeof(uint32_t));
            cof_succ[i] = (uint32_t*) malloc ((mfb1+1) * sizeof(uint32_t));
            for (int j = 0; j <= mfb1; j++)
                cof_call[i][j] = cof_succ[i][j] = 0;
        }
    } else {
        cof_stats_file = NULL;
        cof_call = NULL;
        cof_succ = NULL;
    }
}
void las_info::print_cof_stats()
{
    if (!cof_stats_file) return;
    ASSERT_ALWAYS(config_pool.default_config_ptr);
    siever_config const & sc0(*config_pool.default_config_ptr);
    int mfb0 = sc0.sides[0].mfb;
    int mfb1 = sc0.sides[1].mfb;
    for (int i = 0; i <= mfb0; i++) {
        for (int j = 0; j <= mfb1; j++)
            fprintf (cof_stats_file, "%u %u %u %u\n", i, j,
                    cof_call[i][j],
                    cof_succ[i][j]);
    }
}

void las_info::clear_cof_stats()
{
    if (!cof_stats_file) return;
    ASSERT_ALWAYS(config_pool.default_config_ptr);
    siever_config const & sc0(*config_pool.default_config_ptr);
    for (int i = 0; i <= sc0.sides[0].mfb; i++) {
        free (cof_call[i]);
        free (cof_succ[i]);
    }
    free (cof_call);
    free (cof_succ);
    fclose (cof_stats_file);
    cof_stats_file = NULL;
}
//}}}

sieve_info & sieve_info::get_sieve_info_from_config(siever_config const & sc, cado_poly_srcptr cpoly, std::list<sieve_info> & registry, cxx_param_list & pl)/*{{{*/
{
    std::list<sieve_info>::iterator psi;
    psi = find_if(registry.begin(), registry.end(), sc.same_config());
    if (psi != registry.end()) {
        sc.display();
        return *psi;
    }
    registry.push_back(sieve_info(sc, cpoly, registry, pl));
    sieve_info & si(registry.back());
    verbose_output_print(0, 1, "# Creating new sieve configuration for q~2^%d on side %d (logI=%d)\n",
            sc.bitsize, sc.side, si.conf.logI_adjusted);
    sc.display();
    return registry.back();
}/*}}}*/
