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
#include "select_mpi.h"

stats_data_t stats; /* struct for printing progress */

void *
thread_sm (void * context_data, earlyparsed_relation_ptr rel)
{
  mpz_poly * abpolys = (mpz_poly *) context_data;
  mpz_poly_init(abpolys[rel->num], EARLYPARSED_RELATION_MAX_AB - 1);
  mpz_poly_setcoeffs_int64(abpolys[rel->num], rel->ab, EARLYPARSED_RELATION_MAX_AB-1);

  return NULL;
}

sm_relset_ptr build_rel_sets(const char * purgedname, const char * indexname,
                             uint64_t * small_nrows, mpz_poly_ptr *F,
			     int nb_polys,
                             const mpz_t ell2)
{
  uint64_t nrows, ncols, len_relset;
  uint64_t r[MAX_LEN_RELSET];
  int64_t e[MAX_LEN_RELSET];
  int ret;

  /* array of (a,b) pairs from (purgedname) file */
  mpz_poly *pairs;

  purgedfile_read_firstline (purgedname, &nrows, &ncols);
  pairs = (mpz_poly *) malloc (nrows * sizeof(mpz_poly));
  ASSERT_ALWAYS (pairs != NULL);
  /* For each rel, read the a,b-pair and init the corresponding poly pairs[] */
  fprintf(stdout, "\n# Reading %" PRIu64 " (a,b) pairs\n", nrows);
  fflush(stdout);
  char *fic[2] = {(char *) purgedname, NULL};
  filter_rels (fic, (filter_rels_callback_t) thread_sm, pairs,
          EARLYPARSE_NEED_AB, NULL, NULL);


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


void print_all_sm(FILE *out, sm_side_info *sm_info, int nb_polys,
    mpz_poly *sm) {
  for(int side = 0, c = 0 ; side < nb_polys ; side++) {
    if (sm_info[side]->nsm == 0)
      continue;
    if (c++) fprintf(out, " ");
    print_sm (out,
        sm[side],
        sm_info[side]->nsm,
        sm_info[side]->f->deg);
  }
  fprintf(out, "\n");
}



// Basic MPI communications

#define MPI_MY_MP_SIZE_T        MPI_LONG
#define MPI_MY_MP_LIMB_T        MPI_UNSIGNED_LONG
#define MPI_MY_GMP_INTERNAL_SIZE_FIELD_T MPI_INT
#define MPI_MY_UINT64_T         MPI_UNSIGNED_LONG

void MPI_Send_mpz(mpz_ptr z, int dst) {
  mp_size_t nlimbs = mpz_size(z);
  MPI_Send(&nlimbs, 1, MPI_MY_MP_SIZE_T, dst, 0, MPI_COMM_WORLD);
  MPI_Send(&z->_mp_size, 1, MPI_MY_GMP_INTERNAL_SIZE_FIELD_T, dst, 0, MPI_COMM_WORLD);
  MPI_Send(z->_mp_d, nlimbs, MPI_MY_MP_LIMB_T, dst, 0, MPI_COMM_WORLD);
}

void MPI_Recv_mpz(mpz_ptr z, int src) {
  mp_size_t nlimbs;
  MPI_Recv(&nlimbs, 1, MPI_MY_MP_SIZE_T, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  _mpz_realloc(z, nlimbs);
  MPI_Recv(&z->_mp_size, 1, MPI_MY_GMP_INTERNAL_SIZE_FIELD_T, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(z->_mp_d, nlimbs, MPI_MY_MP_LIMB_T, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}


void MPI_Send_mpz_poly(mpz_poly_ptr poly, int dst) {
  MPI_Send(&poly->deg, 1, MPI_INT, dst, 0, MPI_COMM_WORLD);
  for (int i = 0; i <= poly->deg; ++i)
    MPI_Send_mpz(poly->coeff[i], dst);
}

void MPI_Recv_mpz_poly(mpz_poly_ptr poly, int src) {
  MPI_Recv(&poly->deg, 1, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  if (poly->alloc < poly->deg+1) {
    poly->alloc = poly->deg+1;
    poly->coeff = (mpz_t *) realloc(poly->coeff, poly->alloc*sizeof(mpz_t));
    ASSERT_ALWAYS(poly->coeff != NULL);
  }
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

struct th_info_s {
  sm_side_info * sm_info;
  sm_relset_ptr rels;
  mpz_poly ** dst;
  uint64_t tot_relset;
  int nb_polys;
};
typedef struct th_info_s th_info_t;

uint64_t next_relset = 0; // All relsets below this are done, or being processed
uint64_t count_processed_sm = 0; // Number of already computed relsets.
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER; 

void * thread_process(void *th_arg) {
  th_info_t * args = (th_info_t *)th_arg;
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  fprintf(stderr, "Starting thread on mpi job %d\n", rank);

  const uint64_t BLOCK = 10;

  while(1) {
    // get a new block of relsets to compute
    pthread_mutex_lock(&mutex);
    if (next_relset >= args->tot_relset) {
      pthread_mutex_unlock(&mutex);
      fprintf(stderr, "Finishing thread on mpi job %d\n", rank);
      return NULL;
    }
    uint64_t first = next_relset;
    uint64_t last = MIN(first + BLOCK, args->tot_relset);
    next_relset = last;
    pthread_mutex_unlock(&mutex);

    // Process the relsets
    for (uint64_t i = first; i < last; ++i) {
      for(int side = 0 ; side < args->nb_polys ; side++) {
        if (args->sm_info[side]->nsm == 0)
          continue;
        mpz_poly_reduce_frac_mod_f_mod_mpz(
            args->rels[i].num[side],
            args->rels[i].denom[side],
            args->sm_info[side]->f0,
            args->sm_info[side]->ell2,
            args->sm_info[side]->invl2
            );
        compute_sm_piecewise(args->dst[i][side],
            args->rels[i].num[side],
            args->sm_info[side]);
      }
    }
    pthread_mutex_lock(&mutex);
    count_processed_sm += BLOCK;
    pthread_mutex_unlock(&mutex);
  }
}


static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "poly", "(required) poly file");
  param_list_decl_usage(pl, "purged", "(required) purged file");
  param_list_decl_usage(pl, "index", "(required) index file");
  param_list_decl_usage(pl, "out", "output file (stdout if not given)");
  param_list_decl_usage(pl, "ell", "(required) group order");
  param_list_decl_usage(pl, "nsm", "number of SM on side 0,1,... (default is "
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
  int mt = 1;
  double t0 = seconds();

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

  param_list_parse_int(pl, "t", &mt);
  if (mt < 1) {
    if (idoio) {
      fprintf(stderr, "Error: parameter -mt must be at least 1\n");
      param_list_print_usage(pl, argv0, stderr);
    }
    exit(EXIT_FAILURE);
  }

  /* Init polynomial */
  cado_poly_init (pol);
  if (!cado_poly_read (pol, polyfile))
  {
    if (idoio) {
      fprintf (stderr, "Error reading polynomial file\n");
    }
    exit (EXIT_FAILURE);
  }

  /* Read number of sm to be printed from command line */
  param_list_parse_int_list (pl, "nsm", nsm_arg, pol->nb_polys, ",");

  for(int side = 0; side < pol->nb_polys; side++)
  {
    F[side] = pol->pols[side];
    if (nsm_arg[side] > F[side]->deg)
    {
      if (idoio) {
        fprintf(stderr, "Error: nsm%d=%d can not exceed the degree=%d\n",
                      side, nsm_arg[side], F[side]->deg);
      }
      exit (EXIT_FAILURE);
    }
  }

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

  sm_side_info sm_info[NB_POLYS_MAX];

  for(int side = 0 ; side < pol->nb_polys ; side++) {
      sm_side_info_init(sm_info[side], F[side], ell);
  }

  for (int side = 0; side < pol->nb_polys; side++) {
    if (idoio) {
      fprintf(stdout, "\n# Polynomial on side %d:\n# F[%d] = ", side, side);
      mpz_poly_fprintf(stdout, F[side]);
      printf("# SM info on side %d:\n", side);
      sm_side_info_print(stdout, sm_info[side]);
    }
    if (nsm_arg[side] >= 0)
      sm_info[side]->nsm = nsm_arg[side]; /* command line wins */
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

  // If nsm is 0 on one side, then set F[side] to NULL to desactivate the
  // corresponding computations.
  // TODO: this will go.
  int dF[NB_POLYS_MAX];
  for (int side = 0; side < pol->nb_polys; ++side) {
    if (sm_info[side]->nsm == 0) {
      F[side] = NULL;
      dF[side] = 0;
    } else {
      dF[side] = F[side]->deg;
    }
  }

  ///////////////////////
  // Only process 0 constructs the relation sets.
  if (rank == 0) {
    rels = build_rel_sets(purgedfile, indexfile, &nb_relsets, F, pol->nb_polys, ell2);
    fprintf(stdout, "\n# Computing Shirokauer maps for %" PRIu64
        " relation-sets.\n", nb_relsets);
    fflush(stdout);
  }
  MPI_Bcast(&nb_relsets, 1, MPI_MY_UINT64_T, 0, MPI_COMM_WORLD);

  ///////////////////////
  // Send a share of the rel sets to each process (round Robin)
  uint64_t nb_parts = (nb_relsets - 1) / size + 1; // ceiling
  sm_relset_ptr part_rels = (sm_relset_ptr)malloc(nb_parts*sizeof(sm_relset_t));
  ASSERT_ALWAYS(part_rels != NULL);
  for (uint64_t i = 0; i < nb_parts; ++i) {
    sm_relset_init(&part_rels[i], dF, pol->nb_polys);
  }
  if (rank == 0) {
    for (uint64_t i = 0; i < nb_parts; ++i) {
      sm_relset_copy(&part_rels[i], &rels[i*size]);
      for (int j = 1; j < size; ++j) {
        if (i*size+j < nb_relsets)
          MPI_Send_relset(&rels[i*size+j], j, pol->nb_polys);
      }
    }
  } else {
    for (uint64_t i = 0; i < nb_parts; ++i) {
      if (i*size+rank < nb_relsets)
        MPI_Recv_relset(&part_rels[i], 0, pol->nb_polys);
    }
  }

  // Can now free the original rels on process 0
  if (rank == 0) {
    for (uint64_t i = 0; i < nb_relsets; i++)
      sm_relset_clear (&rels[i], pol->nb_polys);
    free(rels);
  }

  ///////////////////////
  // Process the relsets.
  t0 = seconds();

  mpz_poly **dst = (mpz_poly **) malloc(nb_parts*sizeof(mpz_poly*));
  for (uint64_t j = 0; j < nb_parts; ++j) {
    dst[j] = (mpz_poly *) malloc(pol->nb_polys*sizeof(mpz_poly));
    memset(dst[j], 0, pol->nb_polys*sizeof(mpz_poly));
    for(int side = 0 ; side < pol->nb_polys ; side++) {
      if (sm_info[side]->nsm != 0)
        mpz_poly_init(dst[j][side], sm_info[side]->f->deg);
    }
  }

  // parameters for the threads (same for all).
  th_info_t th_args;
  th_args.sm_info = sm_info;
  th_args.rels = part_rels;
  th_args.dst = dst;
  th_args.nb_polys = pol->nb_polys;
  th_args.tot_relset = nb_parts;
  if ((nb_parts-1)*size + rank > nb_relsets)
    th_args.tot_relset = nb_parts - 1;
  // start the threads
  next_relset = 0;
  pthread_t *th_id;
  th_id = (pthread_t *) malloc(mt*sizeof(pthread_t));
  ASSERT_ALWAYS(th_id != NULL);
  for (int i = 0; i < mt; ++i) {
    int ret;
    ret = pthread_create(&th_id[i], NULL, &thread_process, (void *)&th_args);
    ASSERT_ALWAYS(ret == 0);
  }

  stats_init(stats, stdout, &count_processed_sm, nbits(nb_relsets)-5,
      "Computed", "SMs", "", "SMs");
  while (count_processed_sm < nb_parts) {
    if (stats_test_progress(stats))
      stats_print_progress (stats, count_processed_sm, 0, 0, 0);
    usleep(1);
  }
  stats_print_progress (stats, nb_parts, 0, 0, 1);

  // join the threads.
  for (int i = 0; i < mt; ++i)
    pthread_join(th_id[i], NULL);
  free(th_id);
  fprintf(stderr, "Job %d: processed all relsets in %f s\n",
      rank, seconds()-t0);

  // Send back results and print
  if (rank != 0) { // sender
    for (uint64_t i = 0; i < nb_parts; ++i) {
      if (i*size+rank < nb_relsets)
        MPI_Send_res(dst[i], 0, sm_info, pol->nb_polys);
    }
  } else { // rank 0 receives and prints. (round Robin again)
    FILE *out = outfile ? fopen(outfile, "w") : stdout;
    ASSERT_ALWAYS(out != NULL);
    int nsm_total=0;
    for (int side = 0; side < pol->nb_polys; side++) {
      nsm_total += sm_info[side]->nsm;
    }
    fprintf(out, "%" PRIu64 " %d", nb_relsets, nsm_total);
    gmp_fprintf(out, " %Zd\n", ell);
    mpz_poly *res;
    res = (mpz_poly *) malloc(pol->nb_polys*sizeof(mpz_poly));
    for(int side = 0 ; side < pol->nb_polys ; side++) {
      mpz_poly_init(res[side], sm_info[side]->f->deg);
    }
    for (uint64_t i = 0; i < nb_parts; ++i) {
      print_all_sm(out, sm_info, pol->nb_polys, dst[i]);
      for (int j = 1; j < size; ++j) {
        if (i*size+j < nb_relsets) {
          MPI_Recv_res(res, j, sm_info, pol->nb_polys);
          print_all_sm(out, sm_info, pol->nb_polys, res);
        }
      }
    }
    for(int side = 0 ; side < pol->nb_polys ; side++)
      mpz_poly_clear(res[side]);
    free(res);
    if (outfile) fclose(out);
  }

  // It's time to free...
  for (uint64_t j = 0; j < nb_parts; ++j) {
    for(int side = 0 ; side < pol->nb_polys ; side++) {
      mpz_poly_clear(dst[j][side]);
    }
    free(dst[j]);
  }
  free(dst);

  for (int side = 0 ; side < pol->nb_polys ; side++)
    sm_side_info_clear(sm_info[side]);

  mpz_clear(ell);
  mpz_clear(ell2);
  cado_poly_clear(pol);
  param_list_clear(pl);

  MPI_Finalize();

  return 0;
}
