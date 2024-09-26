/* Schirokauer maps 
   
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

#include "cado.h" // IWYU pragma: keep

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>    // for PRIu64, SCNu64, SCNd64, SCNx64
#include <stdint.h>      // for uint64_t, int64_t
#include <gmp.h>
#include "cado_poly.h"   // for NB_POLYS_MAX, cado_poly_clear, cado_poly_init
#include "filter_io.h"  // earlyparsed_relation_ptr
#include "gmp_aux.h"    // nbits
#include "gzip.h"       // fopen_maybe_compressed
#include "macros.h"
#include "mpz_poly.h"   // mpz_poly_ptr
#include "omp_proxy.h" // IWYU pragma: keep
#include "params.h"
#include "purgedfile.h" // purgedfile_read_firstline
#include "select_mpi.h"
#include "sm_utils.h"   // sm_relset_ptr
#include "stats.h"      // stats_data_t
#include "timing.h"      // for seconds
#include "verbose.h"    // verbose_output_print

stats_data_t stats; /* struct for printing progress */

void *
thread_sm (void * context_data, earlyparsed_relation_ptr rel)
{
  pair_and_sides * ps = (pair_and_sides *) context_data;
  mpz_poly_set_ab(ps[rel->num]->ab, rel->a, rel->b);
  ps[rel->num]->active_sides[0] = rel->active_sides[0];
  ps[rel->num]->active_sides[1] = rel->active_sides[1];

  return NULL;
}

sm_relset_ptr build_rel_sets(const char * purgedname, const char * indexname,
                             uint64_t * small_nrows, const mpz_poly_srcptr * F,
			     int nb_polys,
                             const mpz_t ell2)
{
  uint64_t nrows, ncols, len_relset;
  uint64_t r[MAX_LEN_RELSET];
  int64_t e[MAX_LEN_RELSET];
  int ret;

  /* array of (a,b) pairs from (purgedname) file */

  purgedfile_read_firstline (purgedname, &nrows, &ncols);

  pair_and_sides * pairs = (pair_and_sides *) malloc (nrows * sizeof(pair_and_sides));
  for (uint64_t i = 0; i < nrows; i++)
    mpz_poly_init (pairs[i]->ab, -1);

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
  ASSERT_ALWAYS (ix != NULL);

  ret = fscanf(ix, "%" SCNu64 "\n", small_nrows);
  ASSERT(ret == 1);

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
    
    sm_relset_init (&rels[i], F, nb_polys);
    sm_build_one_relset (&rels[i], r, e, len_relset, pairs, F, nb_polys, ell2);

    if (stats_test_progress(stats))
      stats_print_progress (stats, i, 0, 0, 0);
  }
  stats_print_progress (stats, *small_nrows, 0, 0, 1);
  fclose_maybe_compressed(ix, indexname);

  for (uint64_t i = 0; i < nrows; i++)
    mpz_poly_clear (pairs[i]->ab);
  free (pairs);
  
  return rels;
}


void print_all_sm(FILE *out, sm_side_info *sm_info, int nb_polys,
    mpz_poly *sm) {
  for(int side = 0, c = 0 ; side < nb_polys ; side++) {
    if (sm_info[side]->nsm == 0)
      continue;
    if (c++) fprintf(out, " ");
    print_sm (out, sm_info[side], sm[side]);
  }
  fprintf(out, "\n");
}



// Basic MPI communications

void MPI_Send_mpz(mpz_ptr z, int dst) {
  mp_size_t nlimbs = mpz_size(z);
  MPI_Send(&nlimbs, 1, CADO_MPI_MP_SIZE_T, dst, 0, MPI_COMM_WORLD);
  MPI_Send(&z->_mp_size, 1, CADO_MPI_MPZ_INTERNAL_SIZE_T, dst, 0, MPI_COMM_WORLD);
  MPI_Send(z->_mp_d, nlimbs, CADO_MPI_MP_LIMB_T, dst, 0, MPI_COMM_WORLD);
}

void MPI_Recv_mpz(mpz_ptr z, int src) {
  mp_size_t nlimbs;
  MPI_Recv(&nlimbs, 1, CADO_MPI_MP_SIZE_T, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  _mpz_realloc(z, nlimbs);
  MPI_Recv(&z->_mp_size, 1, CADO_MPI_MPZ_INTERNAL_SIZE_T, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(z->_mp_d, nlimbs, CADO_MPI_MP_LIMB_T, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}


void MPI_Send_mpz_poly(mpz_poly_ptr poly, int dst) {
  MPI_Send(&poly->deg, 1, MPI_INT, dst, 0, MPI_COMM_WORLD);
  for (int i = 0; i <= poly->deg; ++i)
    MPI_Send_mpz(poly->coeff[i], dst);
}

void MPI_Recv_mpz_poly(mpz_poly_ptr poly, int src) {
  MPI_Recv(&poly->deg, 1, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  ASSERT_ALWAYS(poly->deg + 1 >= 0);
  /* FIXME -- what about already allocated coefficients ? */
  ASSERT_ALWAYS(poly->alloc == 0);
  if (poly->alloc < ((unsigned int) poly->deg+1)) {
    poly->alloc = poly->deg+1;
    poly->coeff = (mpz_t *) realloc(poly->coeff, poly->alloc*sizeof(mpz_t));
    ASSERT_ALWAYS(poly->coeff != NULL);
  }
  /* It pretty much seems that if poly->alloc is non zero on entry, then
   * we have a leak here.
   */
  for (int i = 0; i <= poly->deg; ++i)
    MPI_Recv_mpz(poly->coeff[i], src);
}

void MPI_Send_relset(sm_relset_ptr relset, int dst, int nb_polys) {
  ASSERT_ALWAYS(relset->nb_polys == nb_polys);
  for (int i = 0; i < nb_polys; ++i) {
    MPI_Send_mpz_poly(relset->num[i], dst);
    MPI_Send_mpz_poly(relset->denom[i], dst);
  }
}

void MPI_Recv_relset(sm_relset_ptr relset, int src, int nb_polys) {
  relset->nb_polys = nb_polys;
  for (int i = 0; i < nb_polys; ++i) {
    MPI_Recv_mpz_poly(relset->num[i], src);
    MPI_Recv_mpz_poly(relset->denom[i], src);
  }
}

void MPI_Send_res(mpz_poly * res, int dst, sm_side_info *sm_info,
    int nb_polys) {
  for(int side = 0 ; side < nb_polys ; side++) {
    if (sm_info[side]->nsm == 0)
      continue;
    MPI_Send_mpz_poly(res[side], dst);
  }
}

void MPI_Recv_res(mpz_poly * res, int src, sm_side_info * sm_info,
    int nb_polys) {
  for(int side = 0 ; side < nb_polys ; side++) {
    if (sm_info[side]->nsm == 0)
      continue;
    MPI_Recv_mpz_poly(res[side], src);
  }
}


// Pthread part: on each node, we use shared memory instead of mpi

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "poly", "(required) poly file");
  param_list_decl_usage(pl, "purged", "(required) purged file");
  param_list_decl_usage(pl, "index", "(required) index file");
  param_list_decl_usage(pl, "out", "output file (stdout if not given)");
  param_list_decl_usage(pl, "ell", "(required) group order");
  param_list_decl_usage(pl, "sm-mode", "SM mode (see sm-portability.h)");
  param_list_decl_usage(pl, "nsms", "number of SM on side 0,1,... (default is "
                                   "computed by the program)");
  param_list_decl_usage(pl, "t", "number of threads on each mpi job (default 1)");
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
  MPI_Init(&argc, &argv);
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int idoio = (rank == 0); // Am I the job allowed to do I/O ?
  double t0 = seconds();

  char *argv0 = argv[0];

  const char *polyfile = NULL;
  const char *purgedfile = NULL;
  const char *indexfile = NULL;
  const char *outfile = NULL;

  param_list pl;
  cado_poly cpoly;

  sm_relset_ptr rels = NULL;
  uint64_t nb_relsets;
  mpz_t ell, ell2;

  /* read params */
  param_list_init(pl);
  declare_usage(pl);

  if (argc == 1)
    usage (argv[0], NULL, pl);

  argc--,argv++;
  for ( ; argc ; ) {
    if (param_list_update_cmdline (pl, &argc, &argv)) { continue; }
    if (idoio) {
      fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
      usage (argv0, NULL, pl);
    }
  }
  /* print command-line arguments */
  verbose_interpret_parameters(pl);
  param_list_print_command_line (stdout, pl);

  /* Read poly filename from command line */
  if ((polyfile = param_list_lookup_string(pl, "poly")) == NULL) {
    if (idoio) {
      fprintf(stderr, "Error: parameter -poly is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
    }
    exit(EXIT_FAILURE);
  }

  /* Read purged filename from command line */
  if ((purgedfile = param_list_lookup_string(pl, "purged")) == NULL) {
    if (idoio) {
      fprintf(stderr, "Error: parameter -purged is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
    }
    exit(EXIT_FAILURE);
  }

  /* Read index filename from command line */
  if ((indexfile = param_list_lookup_string(pl, "index")) == NULL) {
    if (idoio) {
      fprintf(stderr, "Error: parameter -index is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
    }
    exit(EXIT_FAILURE);
  }

  /* Read outfile filename from command line ; defaults to stdout. */
  outfile = param_list_lookup_string(pl, "out");

  /* Read ell from command line (assuming radix 10) */
  mpz_init (ell);
  if (!param_list_parse_mpz(pl, "ell", ell)) {
    if (idoio) {
      fprintf(stderr, "Error: parameter -ell is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
    }
    exit(EXIT_FAILURE);
  }

  /* Init polynomial */
  cado_poly_init (cpoly);
  if (!cado_poly_read (cpoly, polyfile))
  {
    if (idoio) {
      fprintf (stderr, "Error reading polynomial file\n");
    }
    exit (EXIT_FAILURE);
  }

  int * nsm_arg = malloc(cpoly->nb_polys * sizeof(int));

  /* negative value means that the value that will be used is the value
   * computed later by sm_side_info_init */
  for (int side = 0; side < cpoly->nb_polys; side++)
    nsm_arg[side] = -1;
  /* Read number of sm to be printed from command line */
  param_list_parse_int_args_per_side(pl, "nsm", nsm_arg, cpoly->nb_polys,
          ARGS_PER_SIDE_DEFAULT_AS_IS);

  mpz_poly_srcptr * F = malloc(cpoly->nb_polys * sizeof(mpz_poly_srcptr));

  for(int side = 0; side < cpoly->nb_polys; side++)
  {
    F[side] = cpoly->pols[side];
    if (nsm_arg[side] > F[side]->deg) {
      if (idoio) {
        fprintf(stderr, "Error: nsm%d=%d can not exceed the degree=%d\n",
                      side, nsm_arg[side], F[side]->deg);
      }
      exit (EXIT_FAILURE);
    }
  }

  const char * sm_mode_string = param_list_lookup_string(pl, "sm-mode");

  if (param_list_warn_unused(pl)) {
    if (idoio) {
      usage (argv0, NULL, pl);
    } else {
      exit (EXIT_FAILURE);
    }
  }

  /* Print ell and ell^2 */
  mpz_init(ell2);
  mpz_mul(ell2, ell, ell);
  if (idoio)
    gmp_fprintf(stdout, "# Sub-group order:\nell = %Zi\n# Computation is done "
                      "modulo ell2 = ell^2:\nell2 = %Zi\n", ell, ell2);

  sm_side_info * sm_info = malloc(cpoly->nb_polys * sizeof(sm_side_info));

  for(int side = 0 ; side < cpoly->nb_polys ; side++) {
      sm_side_info_init(sm_info[side], F[side], ell);
      sm_side_info_set_mode(sm_info[side], sm_mode_string);
  }

  for (int side = 0; side < cpoly->nb_polys; side++) {
    if (idoio) {
      fprintf(stdout, "\n# Polynomial on side %d:\n# F[%d] = ", side, side);
      mpz_poly_fprintf(stdout, F[side]);
      printf("# SM info on side %d:\n", side);
      sm_side_info_print(stdout, sm_info[side]);
    }
    if (nsm_arg[side] >= 0) {
      if (idoio && nsm_arg[side] < sm_info[side]->nsm)
        printf("# Warning: as per the command line parameters, we will only compute %d SMs instead of %d on side %d. Be sure to know what you are doing!\n", nsm_arg[side], sm_info[side]->nsm, side);
      sm_info[side]->nsm = nsm_arg[side]; /* command line wins */
    }
    if (idoio)
      printf("# Will compute %d SMs on side %d\n", sm_info[side]->nsm, side);

    /* do some consistency checks */
    if (sm_info[side]->unit_rank != sm_info[side]->nsm)
    {
      if (idoio)
        fprintf(stderr, "# On side %d, unit rank is %d, computing %d SMs ; "
          "weird.\n", side, sm_info[side]->unit_rank,
          sm_info[side]->nsm);
      /* for the 0 case, we haven't computed anything: prevent the
       * user from asking SM data anyway */
      ASSERT_ALWAYS(sm_info[side]->unit_rank != 0);
    }
  }
  fflush(stdout);

  // If nsm (or nsm_arg) is 0 on one side, then set F[side] to NULL to
  // desactivate the corresponding computations.
  // (note that it does not seem to be a good idea to forcibly deactivate
  // the SMs on a given side, I think).
  for (int side = 0; side < cpoly->nb_polys; ++side) {
      if (sm_info[side]->nsm == 0) {
          F[side] = NULL;
      }
  }

  ///////////////////////
  // Only process 0 constructs the relation sets.
#ifdef HAVE_OPENMP
      unsigned int thmax = omp_get_max_threads();
#else
      const unsigned int thmax = 1;
#endif
  if (rank == 0) {
    rels = build_rel_sets(purgedfile, indexfile, &nb_relsets, F, cpoly->nb_polys, ell2);
    fprintf(stdout, "\n# Computing Shirokauer maps for %" PRIu64
        " relation-sets, using %d threads and %d jobs.\n", nb_relsets, thmax, size);
    fflush(stdout);
  }
  MPI_Bcast(&nb_relsets, 1, CADO_MPI_UINT64_T, 0, MPI_COMM_WORLD);

  ///////////////////////
  // Send a share of the rel sets to each process (round Robin)
  uint64_t nb_parts = (nb_relsets - 1) / size + 1; // ceiling
  sm_relset_ptr part_rels = (sm_relset_ptr)malloc(nb_parts*sizeof(sm_relset_t));
  ASSERT_ALWAYS(part_rels != NULL);
  for (uint64_t i = 0; i < nb_parts; ++i) {
    sm_relset_init(&part_rels[i], F, cpoly->nb_polys);
  }
  if (rank == 0) {
    for (uint64_t i = 0; i < nb_parts; ++i) {
      sm_relset_copy(&part_rels[i], &rels[i*size]);
      for (int j = 1; j < size; ++j) {
        if (i*size+j < nb_relsets)
          MPI_Send_relset(&rels[i*size+j], j, cpoly->nb_polys);
      }
    }
  } else {
    for (uint64_t i = 0; i < nb_parts; ++i) {
      if (i*size+rank < nb_relsets)
        MPI_Recv_relset(&part_rels[i], 0, cpoly->nb_polys);
    }
  }

  // Can now free the original rels on process 0
  if (rank == 0) {
    for (uint64_t i = 0; i < nb_relsets; i++)
      sm_relset_clear (&rels[i]);
    free(rels);
  }

  ///////////////////////
  // Process the relsets.
  t0 = seconds();

  mpz_poly **dst = (mpz_poly **) malloc(nb_parts*sizeof(mpz_poly*));
  for (uint64_t j = 0; j < nb_parts; ++j) {
    dst[j] = (mpz_poly *) malloc(cpoly->nb_polys*sizeof(mpz_poly));
    memset(dst[j], 0, cpoly->nb_polys*sizeof(mpz_poly));
    for(int side = 0 ; side < cpoly->nb_polys ; side++) {
      if (sm_info[side]->nsm != 0)
        mpz_poly_init(dst[j][side], sm_info[side]->f->deg);
    }
  }

  /* updated only by thread 0 */
  uint64_t count_processed_sm = 0;
  stats_init(stats, stdout, &count_processed_sm, nbits(nb_relsets)-5,
          "Computed", "SMs", "", "SMs");

#ifdef HAVE_OPENMP
#pragma omp parallel
#endif
  {
#ifdef HAVE_OPENMP
      unsigned int thid = omp_get_thread_num();
      unsigned int thnb = omp_get_num_threads();
#else
      const unsigned int thid = 0;
      const unsigned int thnb = 1;
#endif
#ifdef HAVE_OPENMP
#pragma omp for schedule(static)
#endif
      for(uint64_t i = 0 ; i < nb_relsets ; i++) {
          for(int side = 0 ; side < cpoly->nb_polys ; side++) {
              if (sm_info[side]->nsm == 0)
                  continue;
              mpz_poly_reduce_frac_mod_f_mod_mpz(
                      part_rels[i].num[side],
                      part_rels[i].denom[side],
                      sm_info[side]->f0,
                      sm_info[side]->ell2);
              compute_sm_piecewise(dst[i][side],
                      part_rels[i].num[side],
                      sm_info[side]);
          }
          if (thid == 0) {
              /* Static schedule, all threads can reasonably be expected
               * to progress at the same speed */
              count_processed_sm+=thnb;
              if (stats_test_progress(stats))
                  stats_print_progress (stats, count_processed_sm, 0, 0, 0);
          }
      }
  }
  stats_print_progress (stats, nb_parts, 0, 0, 1);

  fprintf(stderr, "Job %d: processed all relsets in %f s\n",
      rank, seconds()-t0);

  // Send back results and print
  if (rank != 0) { // sender
    for (uint64_t i = 0; i < nb_parts; ++i) {
      if (i*size+rank < nb_relsets)
        MPI_Send_res(dst[i], 0, sm_info, cpoly->nb_polys);
    }
  } else { // rank 0 receives and prints. (round Robin again)
    FILE *out = outfile ? fopen(outfile, "w") : stdout;
    ASSERT_ALWAYS(out != NULL);
    int nsm_total=0;
    for (int side = 0; side < cpoly->nb_polys; side++) {
      nsm_total += sm_info[side]->nsm;
    }
    fprintf(out, "%" PRIu64 " %d", nb_relsets, nsm_total);
    gmp_fprintf(out, " %Zd\n", ell);
    mpz_poly *res;
    res = (mpz_poly *) malloc(cpoly->nb_polys*sizeof(mpz_poly));
    for(int side = 0 ; side < cpoly->nb_polys ; side++) {
      mpz_poly_init(res[side], sm_info[side]->f->deg);
    }
    for (uint64_t i = 0; i < nb_parts; ++i) {
      print_all_sm(out, sm_info, cpoly->nb_polys, dst[i]);
      for (int j = 1; j < size; ++j) {
        if (i*size+j < nb_relsets) {
          MPI_Recv_res(res, j, sm_info, cpoly->nb_polys);
          print_all_sm(out, sm_info, cpoly->nb_polys, res);
        }
      }
    }
    for(int side = 0 ; side < cpoly->nb_polys ; side++)
      mpz_poly_clear(res[side]);
    free(res);
    if (outfile) fclose(out);
  }

  // It's time to free...
  for (uint64_t j = 0; j < nb_parts; ++j) {
    for(int side = 0 ; side < cpoly->nb_polys ; side++) {
      mpz_poly_clear(dst[j][side]);
    }
    free(dst[j]);
  }
  free(dst);
  for (uint64_t i = 0; i < nb_parts; ++i) {
    sm_relset_clear(&part_rels[i]);
  }
  free(part_rels);
  free(F);
  free(nsm_arg);

  for (int side = 0 ; side < cpoly->nb_polys ; side++)
    sm_side_info_clear(sm_info[side]);
  free(sm_info);

  mpz_clear(ell);
  mpz_clear(ell2);
  cado_poly_clear(cpoly);
  param_list_clear(pl);

  MPI_Finalize();

  return 0;
}
