#include "cado.h"

#include "las-sieve-info.hpp"

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
sieve_info::sieve_info(siever_config const & sc, cxx_cado_poly const & cpoly, std::list<sieve_info> & sievers, cxx_param_list & pl, bool try_fbc) /*{{{*/
    : cpoly_ptr(&cpoly), conf(sc)
{
    I = 1UL << sc.logI;

    std::list<sieve_info>::iterator psi;

    /*** factor base ***/
    
    /* This is screwed.
     *
     * We wish to address the case where depending on the special-q,
     * there may be different sieving bounds (lim). And then, this causes
     * us to read new factor bases.
     *
     * This is screwed because morally, there should be one single factor
     * base. It should be truncated if needed, when creating a slicing.
     *
     * How exactly this should be done is not totally clear. Maybe along
     * the following lines.
     *
     * The factor base can be initialized essentially with the
     * polynomials, lim, and powlim.
     *
     * The siever_config that comes now *may* be used to restrict the
     * factor base: it defines a new fb_factorbase::key_type, which then
     * is used to create a slicing.
     */

    psi = find_if(sievers.begin(), sievers.end(), sc.same_fb_parameters());

    if (psi != sievers.end()) {
        sieve_info & other(*psi);
        verbose_output_print(0, 1, "# copy factor base data from previous siever\n");
        sides[0].fb = other.sides[0].fb;
        sides[1].fb = other.sides[1].fb;
    } else {
        const char * fbc_filename = param_list_lookup_string(pl, "fbc");
        for(int side = 0 ; side < 2 ; side++) {
            unsigned long lim = sc.sides[side].lim;
            unsigned long powlim = sc.sides[side].powlim;
            sides[side].fb=std::make_shared<fb_factorbase>(cpoly, side, lim, powlim, pl, try_fbc ? fbc_filename : NULL);
        }
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

void sieve_info::update_norm_data ()/*{{{*/
{
    for(int side = 0 ; side < 2 ; side++) {
        sides[side].lognorms = std::make_shared<lognorm_smart>(conf, *cpoly_ptr, side, qbasis, conf.logI, J);
    }
}

/*}}}*/

/* This is one of the most terrible misnomers in cado-nfs... */
void sieve_info::update (unsigned int nb_threads)
{
    uint64_t A = UINT64_C(1) << conf.logA;

    for(int side = 0 ; side < 2 ; side++) {
        sieve_info::side_info & sis(sides[side]);
        if (!sis.fb) continue;

        fb_factorbase::key_type K = conf.instantiate_thresholds(side);

        K.scale = sis.lognorms->scale,
        /* The size of the slices must be in accordance to our
         * multithread setting. Hopefully, that one does not change...
         * We want to divide in small enough pieces, so that the amout of
         * work sums up more or less evenly across the threads.
         *
         * This "small enough" criterion used to be computed around this
         * place, and has now moved to inside the make_slices call. It is
         * inherently tied to the count of the entries in each part.
         */
        K.nb_threads = nb_threads;

        sis.fbK = K;
        sis.fbs = &(*sis.fb)[K];
    }

    for(int side = 0 ; side < 2 ; side++) {
	init_trialdiv(side); /* Init refactoring stuff */
        sieve_info::side_info & sis(sides[side]);
        if (!sis.fb) continue;

        {
            /* TODO: In this block, we're doing just a copy. We should
             * kill that, and use the field from the slicing directly. */
            const fb_factorbase::slicing * fbs = sis.fbs;
            /* put the resieved primes first. */
            auto s = std::make_shared<std::vector<fb_entry_general>>();
            std::vector<fb_entry_general> const & RS(fbs->small_sieve_entries.resieved);
            s->insert(s->end(), RS.begin(), RS.end());
            sis.resieve_start_offset = 0;
            sis.resieve_end_offset = s->size();
            std::vector<fb_entry_general> const & R(fbs->small_sieve_entries.rest);
            s->insert(s->end(), R.begin(), R.end());
            sis.fb_smallsieved = s;
        }

        /*
        verbose_output_print(0, 2, "# small side-%d factor base", side);
        size_t nr_roots;
        sides[side].fb->get_part(0)->count_entries(NULL, &nr_roots, NULL);
        verbose_output_print(0, 2, " (total %zu)\n", nr_roots);
        */
    }

    // Now that fb have been initialized, we can set the toplevel.
    toplevel = -1;
    for(int side = 0 ; side < 2 ; side++) {
        if (!sides[side].fb) continue;
        int level = sides[side].fbs->get_toplevel();
        if (level > toplevel) toplevel = level;
    }

    /* update number of buckets at toplevel */
    size_t (&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;

    for(int i = 0 ; i < FB_MAX_PARTS ; ++i) nb_buckets[i] = 0;

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
}/*}}}*/
