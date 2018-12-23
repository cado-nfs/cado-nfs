#include "cado.h"
#include <stdio.h>
#include <stdarg.h>
#include <gmp.h>
#include "las-info.hpp"
#include "las-config.h"
#include "las-norms.hpp"
#include "misc.h"
#include "memusage.h"

/* las_info stuff */

void las_info::configure_aliases(cxx_param_list & pl)
{
    cxx_cado_poly::configure_aliases(pl);
}

void las_info::configure_switches(cxx_param_list & pl)
{
    cxx_cado_poly::configure_switches(pl);
    param_list_configure_switch(pl, "-allow-compsq", NULL);
    param_list_configure_switch(pl, "-dup", NULL);
    param_list_configure_switch(pl, "-batch", NULL);
}

void las_info::declare_usage(cxx_param_list & pl)
{
    cxx_cado_poly::declare_usage(pl);
    siever_config_pool::declare_usage(pl);
    sieve_shared_data::declare_usage(pl);
#ifdef  DLP_DESCENT
    las_dlog_base::declare_usage(pl);
#endif
    cofactorization_statistics::declare_usage(pl);


    param_list_decl_usage(pl, "seed", "Use this seed for random state seeding (currently used only by --random-sample)");

    param_list_decl_usage(pl, "galois", "depending on the specified galois automorphism, sieve only part of the q's");

    /* Note: also declared by las_todo_list ! */
    param_list_decl_usage(pl, "allow-compsq", "allows composite special-q");
    param_list_decl_usage(pl, "qfac-min", "factors of q must be at least that");
    param_list_decl_usage(pl, "qfac-max", "factors of q must be at most that");




    param_list_decl_usage(pl, "dup", "suppress duplicate relations");
    param_list_decl_usage(pl, "dup-qmin", "lower limit of global q-range for 2-sided duplicate removal");
    param_list_decl_usage(pl, "dup-qmax", "upper limit of global q-range for 2-sided duplicate removal");
    param_list_decl_usage(pl, "adjust-strategy", "strategy used to adapt the sieving range to the q-lattice basis (0 = logI constant, J so that boundary is capped; 1 = logI constant, (a,b) plane norm capped; 2 = logI dynamic, skewed basis; 3 = combine 2 and then 0) ; default=0");


    param_list_decl_usage(pl, "batch", "use batch cofactorization");
    param_list_decl_usage(pl, "batch0", "side-0 batch file");
    param_list_decl_usage(pl, "batch1", "side-1 batch file");
    param_list_decl_usage(pl, "batchmfb0", "cofactor bound on side 0 to be considered after batch cofactorization. After primes below 2^batchlpb0 have been extracted, cofactors below this bound will go through ecm. Defaults to lpb0.");
    param_list_decl_usage(pl, "batchmfb1", "cofactor bound on side 1 to be considered after batch cofactorization. After primes below 2^batchlpb1 have been extracted, cofactors below this bound will go through ecm. Defaults to lpb1.");
    param_list_decl_usage(pl, "batchlpb0", "large prime bound on side 0 to be considered by batch cofactorization. Primes between lim0 and 2^batchlpb0 will be extracted by product trees. Defaults to lpb0.");
    param_list_decl_usage(pl, "batchlpb1", "large prime bound on side 1 to be considered by batch cofactorization. Primes between lim1 and 2^batchlpb1 will be extracted by product trees. Defaults to lpb1.");
    param_list_decl_usage(pl, "batch-print-survivors", "just print survivors to the specified file (or pipe) for an external cofactorization");


    param_list_decl_usage(pl, "dumpfile", "Dump entire sieve region to file for debugging.");
}


void las_info::prepare_sieve_shared_data(cxx_param_list & pl)
{
#ifdef HAVE_HWLOC
    shared_structure_cache.clear();
    set_loose_binding();
    for(int k = 0 ; k < number_of_subjobs_total() ; k+= number_of_subjobs_per_memory_binding_zone()) {
        set_subjob_mem_binding(k);
        /* gcc pre 4.7 does not have map::emplace ; while it's
         * admittedly not a c++11-complete compiler anyway, we still have
         * hope to use such an ancient thing on a few machines (old *bsd,
         * for example).
         *
         * See:
         * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=44436#c30
         */
#if 0
        shared_structure_cache.emplace(
                current_memory_binding(),
                sieve_shared_data(cpoly, pl));
#else
        std::pair<cxx_hwloc_nodeset, sieve_shared_data> v(
                current_memory_binding(),
                sieve_shared_data(cpoly, pl));
        shared_structure_cache.insert(std::move(v));
        /* for this one, the default ctor is good enough */
        las_memory_accessor_cache[current_memory_binding()];
#endif
    }
    set_loose_binding();
#else
    shared_structure_private = sieve_shared_data(cpoly, pl);
    /* the default-constructed memory accessor is fine */
#endif
}

void las_info::load_factor_base(cxx_param_list & pl)
{
#ifdef HAVE_HWLOC
    set_loose_binding();
    for(int k = 0 ; k < number_of_subjobs_total() ; k+= number_of_subjobs_per_memory_binding_zone()) {
        set_subjob_mem_binding(k);
        /* right, so at this point we would probably need to 
         * compute the factor base just once, and copy it in ram, isn't
         * it ?
         */
        shared_structure_cache.at(current_memory_binding()).load_factor_base(pl, number_of_threads_loose());
    }
    set_loose_binding();
#else
    shared_structure_private.load_factor_base(pl, number_of_threads_loose());
#endif
}

las_info::las_info(cxx_param_list & pl)
    : cpoly(pl),
      config_pool(pl),
#ifdef HAVE_HWLOC
      shared_structure_cache(),
#else
      shared_structure_private(cpoly, pl),
#endif
#ifdef  DLP_DESCENT
      dlog_base(pl),
#endif
      cofac_stats(pl)
      /*{{{*/
{
    /* We strive to initialize things in the exact order they're written
     * in the struct */
    // ----- general operational flags {{{
    gmp_randinit_default(rstate);
    unsigned long seed = 0;
    if (param_list_parse_ulong(pl, "seed", &seed))
        gmp_randseed_ui(rstate, seed);

    galois = param_list_lookup_string(pl, "galois");
    suppress_duplicates = param_list_parse_switch(pl, "-dup");

    if (const char * tmp = param_list_lookup_string(pl, "bkmult")) {
        bk_multiplier = bkmult_specifier(tmp);
    }

    param_list_parse_int(pl, "adjust-strategy", &adjust_strategy);

    // }}}


    /* composite special-q ? Note: this block is present both in
     * las-todo-list.cpp and las-info.cpp */
    if ((allow_composite_q = param_list_parse_switch(pl, "-allow-compsq"))) {
        /* defaults are set in the class description */
        param_list_parse_uint64(pl, "qfac-min", &qfac_min);
        param_list_parse_uint64(pl, "qfac-max", &qfac_max);
    }

    // ----- stuff roughly related to the descent {{{
#ifdef  DLP_DESCENT
    descent_helper = NULL;
#endif
    // }}}

    /* {{{ duplicate suppression */
    dupqmin = {{ 0, 0 }};
    dupqmax = {{ ULONG_MAX, ULONG_MAX}};
    if (!param_list_parse_ulong_and_ulong(pl, "dup-qmin", &dupqmin[0], ",") && suppress_duplicates) {
        fprintf(stderr, "Error: -dup-qmin is mandatory with -dup\n");
        exit(EXIT_FAILURE);
    }
    param_list_parse_ulong_and_ulong(pl, "dup-qmax", &dupqmax[0], ",");
    /* Change 0 (not initialized) into LONG_MAX */
    for (auto & x : dupqmin) if (x == 0) x = ULONG_MAX;

    /* }}} */

    // ----- batch mode {{{
    batch = param_list_parse_switch(pl, "-batch");

    if (batch) {
        ASSERT_ALWAYS(config_pool.default_config_ptr);
        siever_config const & sc0(*config_pool.default_config_ptr);
	batchlpb[0] = sc0.sides[0].lpb;
	batchlpb[1] = sc0.sides[1].lpb;
	batchmfb[0] = sc0.sides[0].lpb;
	batchmfb[1] = sc0.sides[1].lpb;
        param_list_parse_int(pl, "batchlpb0", &(batchlpb[0]));
        param_list_parse_int(pl, "batchlpb1", &(batchlpb[1]));
        param_list_parse_int(pl, "batchmfb0", &(batchmfb[0]));
        param_list_parse_int(pl, "batchmfb1", &(batchmfb[1]));
	batch_file[0] = param_list_lookup_string (pl, "batch0");
	batch_file[1] = param_list_lookup_string (pl, "batch1");

        for(int side = 0 ; side < 2 ; side++) {
            // the product of primes up to B takes \log2(B)-\log\log 2 /
            // \log 2 bits. The added constant is 0.5287.
            if (batchlpb[side] + 0.5287 >= 31 + log2(GMP_LIMB_BITS)) {
                fprintf(stderr, "Gnu MP cannot deal with primes product that large (max 37 bits, asked for batchlpb%d=%d)\n", side, batchlpb[side]);
                abort();
            } else if (batchlpb[side] + 0.5287 >= 34) {
                fprintf(stderr, "Gnu MP's mpz_inp_raw and mpz_out_raw functions are limited to integers of at most 34 bits (asked for batchlpb%d=%d)\n",side,batchlpb[side]);
                abort();
            }
        }
    }


    const char * bps = param_list_lookup_string(pl, "batch-print-survivors");
    if (bps) {
        batch_print_survivors = fopen(bps, "w");
        ASSERT_ALWAYS(batch_print_survivors != NULL);
    } else {
        batch_print_survivors = NULL;
    }
    // }}} 

    dump_filename = param_list_lookup_string(pl, "dumpfile");
}/*}}}*/


las_info::~las_info()/*{{{*/
{

    // ----- general operational flags {{{
    gmp_randclear(rstate);
    // }}}

    // ----- batch mode: very little
    if (batch_print_survivors) fclose(batch_print_survivors);
}/*}}}*/
