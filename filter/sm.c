/* Shirokauer maps 
   
Input:

* A list of the (npurged) a,b pairs. This is obtained from the
  purgedfile.
* A matrix of (small_nrows) rows and (npurged) cols, which indicates
  the contents of each relation-set. This is obtained from the
  indexfile.
* The sub-group order (ell) such that ell | p-1
  Note: All computations are done mod ell^2.
* (eps): the exponent used in the computation of the Shirokauer maps.
  Note: eps = ppcm(eps_i), where eps_i = ell^(deg(f_i)) - 1 and f = f_1 ... f_k mod ell
  
Output

* A matrix of (small_nrows) rows and (nmaps)=deg(f) cols (mpz_t).  For each
  relation (rel) the (nmaps) Shirokauer maps are computed as the second
  least-significant digit of the ell-adic representation of the polynomial 
  equal to (rel^eps - 1) / ell.

  In case of two algebraic sides, SM's are computed for sides 0..1 in that
  order.

*/

#include "cado.h"

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pthread.h>
#include <errno.h>

#include "macros.h"
#include "utils_with_io.h"
#include "filter_config.h"

stats_data_t stats; /* struct for printing progress */

void *
thread_sm (void * context_data, earlyparsed_relation_ptr rel)
{
  mpz_poly_t * abpolys = (mpz_poly_t *) context_data;
  mpz_poly_init_set_ab(abpolys[rel->num], rel->a, rel->b);

  return NULL;
}

sm_relset_ptr build_rel_sets(const char * purgedname, const char * indexname,
                             uint64_t * small_nrows, mpz_poly_ptr *F,
			     int nb_polys,
                             const mpz_t ell2)
{
  uint64_t nrows, ncols, small_ncols, len_relset;
  uint64_t r[MAX_LEN_RELSET];
  int64_t e[MAX_LEN_RELSET];
  int ret;

  /* array of (a,b) pairs from (purgedname) file */
  mpz_poly_t *pairs;

  purgedfile_read_firstline (purgedname, &nrows, &ncols);
  pairs = (mpz_poly_t *) malloc (nrows * sizeof(mpz_poly_t));
  ASSERT_ALWAYS (pairs != NULL);
  /* For each rel, read the a,b-pair and init the corresponding poly pairs[] */
  fprintf(stdout, "\n# Reading %" PRIu64 " (a,b) pairs\n", nrows);
  fflush(stdout);
  char *fic[2] = {(char *) purgedname, NULL};
  filter_rels (fic, (filter_rels_callback_t) thread_sm, pairs,
          EARLYPARSE_NEED_AB_HEXA, NULL, NULL);


  /* Array of (small_nrows) relation-sets built from array (pairs) and
     (indexname) file  */
  sm_relset_ptr rels;
  FILE * ix = fopen_maybe_compressed(indexname, "r");

  /* small_ncols isn't used here: we don't care */
  ret = fscanf(ix, "%" SCNu64 " %" SCNu64 "", small_nrows, &small_ncols);
  ASSERT(ret == 2);

  rels = (sm_relset_ptr) malloc (*small_nrows * sizeof(sm_relset_t));
  ASSERT_ALWAYS (rels != NULL);

  fprintf(stdout, "\n# Building %" PRIu64 " relation-sets\n", *small_nrows);
  fflush(stdout);
  uint64_t i;
  stats_init (stats, stdout, &i, nbits(*small_nrows)-5, "Computed",
              "relation-sets", "", "relsets");
  for(i = 0 ; i < *small_nrows ; i++)
  {
    ret = fscanf(ix, "%" SCNu64 "", &len_relset);
    ASSERT_ALWAYS(ret == 1 && len_relset < MAX_LEN_RELSET);

    for (uint64_t k = 0 ; k < len_relset ; k++)
    {
      ret = fscanf(ix, " %" SCNx64 ":%" SCNd64 "", &r[k], &e[k]); 
      ASSERT_ALWAYS(ret == 2);
    }
    
    int dF[NB_POLYS_MAX];
    for (int s = 0; s < nb_polys; ++s) {
        if (F[s] == NULL)
	    dF[s] = 0;
	else
	    dF[s] = F[s]->deg;
    }
    sm_relset_init (&rels[i], dF, nb_polys);
    sm_build_one_relset (&rels[i], r, e, len_relset, pairs, F, nb_polys, ell2);

    if (stats_test_progress(stats))
      stats_print_progress (stats, i, 0, 0, 0);
  }
  stats_print_progress (stats, *small_nrows, 0, 0, 1);
  fclose_maybe_compressed(ix, indexname);

  for (uint64_t i = 0; i < nrows; i++)
    mpz_poly_clear (pairs[i]);
  free (pairs);
  
  return rels;
}

struct thread_info {
  int offset;
  int nb, nb_polys;
  sm_relset_ptr rels;
  sm_side_info * sm_info;
  /* where we are supposed to write our result */
  mpz_poly_t ** dst;
};

void * thread_start(void *arg) {
    struct thread_info *ti = (struct thread_info *) arg;
    sm_relset_ptr rels = ti->rels;
    int offset = ti->offset;

    for (int i = 0; i < ti->nb; i++) {
        for(int side = 0 ; side < ti->nb_polys ; side++) {
            if (ti->sm_info[side]->nsm == 0)
                continue;

            mpz_poly_reduce_frac_mod_f_mod_mpz(
                    rels[offset+i].num[side],
                    rels[offset+i].denom[side],
                    ti->sm_info[side]->f0,
                    ti->sm_info[side]->ell2,
                    ti->sm_info[side]->invl2
                    );

            compute_sm_piecewise(ti->dst[i][side],
                    rels[offset+i].num[side],
                    ti->sm_info[side]);
        }
    }
    return NULL;
}

uint64_t print_thread_result(FILE * out, struct thread_info * ti)
{
    uint64_t out_cpt = 0;
    for (int k = 0; k < ti->nb; ++k, ++out_cpt) {
        for(int side = 0, c = 0 ; side < ti->nb_polys ; side++) {
            if (ti->sm_info[side]->nsm == 0)
                continue;
            if (c++) fprintf(out, " ");
            print_sm (out,
                    ti->dst[k][side],
                    ti->sm_info[side]->nsm,
                    ti->sm_info[side]->f->deg);
        }
        fprintf(out, "\n");
    }
    return out_cpt;
}

#define SM_BATCH_SIZE 512

void mt_sm (int nt, const char * outname, sm_relset_ptr rels, 
	    uint64_t nb_relsets, mpz_srcptr ell, sm_side_info * sm_info,
	    int nb_polys)
{
  // We'll use a rotating buffer of thread id.
  pthread_t *threads;
  threads = (pthread_t *) malloc(nt*sizeof(pthread_t));
  int active_threads = 0;  // number of running threads
  int threads_head = 0;    // next thread to wait / restart.
  
  // Prepare the main loop
  uint64_t i = 0; // counter of relation-sets.
  uint64_t out_cpt = 0; // counter of already printed relation-sets;
  FILE * out = outname ? fopen(outname, "w") : stdout;
  DIE_ERRNO_DIAG(out==NULL, "fopen", outname);
  int nsm_total=0;
  for (int side = 0; side < nb_polys; side++) {
      nsm_total += sm_info[side]->nsm;
  }
  /*
  gmp_fprintf(out, "%" PRIu64 " %d %Zd\n", nb_relsets, nsm_total, ell);
  */
  /* mingw's gmp chokes on the %I64u format string which is used by
   * windows as a real value for PRIu64...
   */
  fprintf(out, "%" PRIu64 " %d", nb_relsets, nsm_total);
  gmp_fprintf(out, " %Zd\n", ell);

  
  // Arguments for threads
  struct thread_info *tis;
  tis = (struct thread_info*) malloc(nt*sizeof(struct thread_info));
  for (int i = 0; i < nt; ++i) {
    tis[i].rels = rels;
    tis[i].dst = (mpz_poly_t **) malloc(SM_BATCH_SIZE*sizeof(mpz_poly_t*));
    tis[i].nb_polys = nb_polys;
    for (int j = 0; j < SM_BATCH_SIZE; ++j) {
        tis[i].dst[j] = (mpz_poly_t *) malloc(nb_polys*sizeof(mpz_poly_t));
        memset(tis[i].dst[j], 0, nb_polys*sizeof(mpz_poly_t));
        for(int side = 0 ; side < nb_polys ; side++) {
            if (sm_info[side]->nsm != 0)
                mpz_poly_init(tis[i].dst[j][side],
                        sm_info[side]->f->deg);
        }
    }
    tis[i].sm_info = sm_info;
    // offset and nb must be adjusted.
  }

  // Main loop
  stats_init (stats, stdout, &out_cpt, nbits(nb_relsets)-5, "Computed", "SMs", "", "SMs");
  while ((i < nb_relsets) || (active_threads > 0)) {
    // Start / restart as many threads as allowed
    if ((active_threads < nt) && (i < nb_relsets)) { 
      tis[threads_head].offset = i;
      tis[threads_head].nb = MIN(SM_BATCH_SIZE, nb_relsets-i);
      pthread_create(&threads[threads_head], NULL, 
          &thread_start, (void *)(&tis[threads_head]));
      i += SM_BATCH_SIZE;
      active_threads++;
      threads_head++; 
      if (threads_head == nt) 
        threads_head = 0;
      continue;
    }
    // Wait for the next thread to finish in order to print result.
    pthread_join(threads[threads_head], NULL);
    active_threads--;

    out_cpt += print_thread_result(out, &(tis[threads_head]));

    // report
    if (stats_test_progress(stats))
      stats_print_progress (stats, out_cpt, 0, 0, 0);

    // If we are at the end, no job will be restarted, but head still
    // must be incremented.
    if (i >= nb_relsets) { 
      threads_head++;
      if (threads_head == nt) 
        threads_head = 0;
    }
  }
  stats_print_progress (stats, nb_relsets, 0, 0, 1);

  if (outname) fclose(out);
  for (int i = 0; i < nt; ++i) {
      for (int j = 0; j < SM_BATCH_SIZE; ++j) {
          for(int side = 0 ; side < nb_polys ; side++) {
              if (sm_info[side]->nsm != 0)
                  mpz_poly_clear(tis[i].dst[j][side]);
          }
          free(tis[i].dst[j]);
      }
      free(tis[i].dst);
  }
  free(tis);
  free(threads);
}

#if 0
/* kept just because it avoids the hairy logic of the one above. But
 * this is obsolete */
void sm (const char * outname, sm_relset_ptr rels, uint64_t nb_relsets,
        mpz_srcptr ell, sm_side_info * sm_info)
{
  // Prepare the main loop
  FILE * out = outname ? fopen(outname, "w") : stdout;
  DIE_ERRNO_DIAG(out==NULL, "fopen", outname);
  int nsm_total=0;
  for (int side = 0; side < 2; side++) {
      nsm_total += sm_info[side]->nsm;
  }
  /*
  gmp_fprintf(out, "%" PRIu64 " %d %Zd\n", nb_relsets, nsm_total, ell);
  */
  /* mingw's gmp chokes on the %I64u format string which is used by
   * windows as a real value for PRIu64...
   */
  fprintf(out, "%" PRIu64 " %d", nb_relsets, nsm_total);
  gmp_fprintf(out, " %Zd\n", ell);

  // Main loop
  uint64_t i;
  stats_init (stats, stdout, &i, nbits(nb_relsets)-5, "Computed", "SMs", "", "SMs");
  mpz_poly_t SM;
  mpz_poly_init(SM, -1);

  for (i = 0; i < nb_relsets; i++) {
    for (int side = 0, c = 0; side < 2; side++) {
      if (sm_info[side]->nsm == 0) continue;
      mpz_poly_reduce_frac_mod_f_mod_mpz (
              rels[i].num[side],
              rels[i].denom[side],
              sm_info[side]->f0,
              sm_info[side]->ell2,
              sm_info[side]->invl2
              );
      compute_sm_straightforward (SM,
              rels[i].num[side],
              sm_info[side]->f0,
              sm_info[side]->ell,
              sm_info[side]->exponent,
              sm_info[side]->ell2,
              sm_info[side]->invl2);
      if (c++) fprintf(out, " ");
      print_sm (out, SM, sm_info[side]->nsm, sm_info[side]->f->deg);
    }
    fprintf(out, "\n");
    // report
    if (stats_test_progress(stats))
      stats_print_progress (stats, i, 0, 0, 0);
  }
  stats_print_progress (stats, nb_relsets, 0, 0, 1);

  if (outname) fclose(out);
}
#endif

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "poly", "(required) poly file");
  param_list_decl_usage(pl, "purged", "(required) purged file");
  param_list_decl_usage(pl, "index", "(required) index file");
  param_list_decl_usage(pl, "out", "output file (stdout if not given)");
  param_list_decl_usage(pl, "ell", "(required) group order");
  param_list_decl_usage(pl, "nsm", "number of SM on side 0,1,... (default is "
                                   "computed by the program)");
  param_list_decl_usage(pl, "t", "number of threads (default 1)");
  verbose_decl_usage(pl);
}

static void usage (const char *argv, const char * missing, param_list pl)
{
  if (missing) {
    fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
        missing);
  }
  param_list_print_usage(pl, argv, stderr);
  exit (EXIT_FAILURE);
}

/* -------------------------------------------------------------------------- */

int main (int argc, char **argv)
{
  char *argv0 = argv[0];

  const char *polyfile = NULL;
  const char *purgedfile = NULL;
  const char *indexfile = NULL;
  const char *outfile = NULL;

  param_list pl;
  cado_poly pol;
  mpz_poly_ptr F[NB_POLYS_MAX];

  sm_relset_ptr rels = NULL;
  uint64_t nb_relsets;
  mpz_t ell, ell2;
  int nsm_arg[NB_POLYS_MAX];
  int mt = 1;
  double t0;

  /* negative value means that the value that will be used is the value
   * computed later by sm_side_info_init */
  for (int side = 0; side < NB_POLYS_MAX; side++)
    nsm_arg[side] = -1;

  /* read params */
  param_list_init(pl);
  declare_usage(pl);

  if (argc == 1)
    usage (argv[0], NULL, pl);

  argc--,argv++;
  for ( ; argc ; ) {
    if (param_list_update_cmdline (pl, &argc, &argv)) { continue; }
    fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
    usage (argv0, NULL, pl);
  }
  /* print command-line arguments */
  verbose_interpret_parameters(pl);
  param_list_print_command_line (stdout, pl);

  /* Read poly filename from command line */
  if ((polyfile = param_list_lookup_string(pl, "poly")) == NULL) {
      fprintf(stderr, "Error: parameter -poly is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  /* Read purged filename from command line */
  if ((purgedfile = param_list_lookup_string(pl, "purged")) == NULL) {
      fprintf(stderr, "Error: parameter -purged is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  /* Read index filename from command line */
  if ((indexfile = param_list_lookup_string(pl, "index")) == NULL) {
      fprintf(stderr, "Error: parameter -index is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  /* Read outfile filename from command line ; defaults to stdout. */
  outfile = param_list_lookup_string(pl, "out");

  /* Read ell from command line (assuming radix 10) */
  mpz_init (ell);
  if (!param_list_parse_mpz(pl, "ell", ell)) {
      fprintf(stderr, "Error: parameter -ell is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_int(pl, "t", &mt);
  if (mt < 1) {
    fprintf(stderr, "Error: parameter mt must be at least 1\n");
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }

  /* Init polynomial */
  cado_poly_init (pol);
  if (!cado_poly_read (pol, polyfile))
  {
    fprintf (stderr, "Error reading polynomial file\n");
    exit (EXIT_FAILURE);
  }

  /* Read number of sm to be printed from command line */
  param_list_parse_int_list (pl, "nsm", nsm_arg, pol->nb_polys, ",");

  for(int side = 0; side < pol->nb_polys; side++)
  {
    F[side] = pol->pols[side];
    if (nsm_arg[side] > F[side]->deg)
    {
      fprintf(stderr, "Error: nsm%d=%d can not exceed the degree=%d\n",
                      side, nsm_arg[side], F[side]->deg);
      exit (EXIT_FAILURE);
    }
  }

  if (param_list_warn_unused(pl))
    usage (argv0, NULL, pl);

  /* Print ell and ell^2 */
  mpz_init(ell2);
  mpz_mul(ell2, ell, ell);
  gmp_fprintf(stdout, "# Sub-group order:\nell = %Zi\n# Computation is done "
                      "modulo ell2 = ell^2:\nell2 = %Zi\n", ell, ell2);

  sm_side_info sm_info[NB_POLYS_MAX];

  for(int side = 0 ; side < pol->nb_polys ; side++) {
      sm_side_info_init(sm_info[side], F[side], ell);
  }

  for (int side = 0; side < pol->nb_polys; side++) {
      fprintf(stdout, "\n# Polynomial on side %d:\n# F[%d] = ", side, side);
      mpz_poly_fprintf(stdout, F[side]);
      printf("# SM info on side %d:\n", side);
      sm_side_info_print(stdout, sm_info[side]);
      if (nsm_arg[side] >= 0)
        sm_info[side]->nsm = nsm_arg[side]; /* command line wins */
      printf("# Will compute %d SMs on side %d\n", sm_info[side]->nsm, side);

      /* do some consistency checks */
      if (sm_info[side]->unit_rank != sm_info[side]->nsm)
      {
        fprintf(stderr, "# On side %d, unit rank is %d, computing %d SMs ; "
                        "weird.\n", side, sm_info[side]->unit_rank,
                        sm_info[side]->nsm);
        /* for the 0 case, we haven't computed anything: prevent the
         * user from asking SM data anyway */
        ASSERT_ALWAYS(sm_info[side]->unit_rank != 0);
      }
  }
  fflush(stdout);

  t0 = seconds();

  // If nsm is 0 on one side, then set F[side] to NULL to desactivate the
  // corresponding computations.
  // TODO: this will go.
  for (int side = 0; side < pol->nb_polys; ++side) {
      if (sm_info[side]->nsm == 0)
          F[side] = NULL;
  }
  rels = build_rel_sets(purgedfile, indexfile, &nb_relsets, F, pol->nb_polys, ell2);

  /* adjust the number of threads based on the number of relations */
  double ntm = ceil((nb_relsets + 0.0)/SM_BATCH_SIZE);
  if (mt > ntm)
    mt = (int) ntm;

  fprintf(stdout, "\n# Computing Shirokauer maps for %" PRIu64 " relation-sets "
                  "using %d thread(s)\n", nb_relsets, mt);
  fflush(stdout);

  mt_sm(mt, outfile, rels, nb_relsets, ell, sm_info, pol->nb_polys);
  // sm(outfile, rels, nb_relsets, ell, sm_info);

  fprintf(stdout, "\n# sm completed in %2.2lf seconds\n", seconds() - t0);
  fflush(stdout);

  for (uint64_t i = 0; i < nb_relsets; i++)
      sm_relset_clear (&rels[i], pol->nb_polys);
  free(rels);

  for (int side = 0 ; side < pol->nb_polys ; side++)
    sm_side_info_clear(sm_info[side]);

  mpz_clear(ell);
  mpz_clear(ell2);
  cado_poly_clear(pol);
  param_list_clear(pl);

  return 0;
}
