#include "cado.h"

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <pthread.h>

#include "filter_common.h"

typedef struct
{
  int64_t *a;
} read_data_t;

void *
thread_read (void * context_data, earlyparsed_relation_ptr rel)
{
  read_data_t *data = (read_data_t *) context_data;
  data->a[rel->num] = rel->a;

  return NULL;
}

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "purged", "(required) purged file");
  param_list_decl_usage(pl, "index", "(required) index file");
  param_list_decl_usage(pl, "out", "output file");
  param_list_decl_usage(pl, "gorder", "(required) group order");
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

  const char *purgedfile = NULL;
  const char *indexfile = NULL;
  const char *outfile = NULL;
  const char *group_order = NULL;

  param_list pl;
  mpz_t ell;
  double t0;

  /* read params */
  param_list_init(pl);
  declare_usage(pl);

  if (argc == 1)
    usage (argv[0], NULL, pl);

  argc--,argv++;
  for ( ; argc ; )
  {
    if (param_list_update_cmdline (pl, &argc, &argv)) { continue; }
    fprintf (stderr, "Unhandled parameter %s\n", argv[0]);
    usage (argv0, NULL, pl);
  }


  if ((purgedfile = param_list_lookup_string(pl, "purged")) == NULL)
  {
    fprintf(stderr, "Error: parameter -purged is mandatory\n");
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }

  if ((indexfile = param_list_lookup_string(pl, "index")) == NULL)
  {
    fprintf(stderr, "Error: parameter -index is mandatory\n");
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }

  if ((group_order = param_list_lookup_string(pl, "gorder")) == NULL)
  {
    fprintf(stderr, "Error: parameter -gorder is mandatory\n");
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }

  if ((outfile = param_list_lookup_string(pl, "out")) == NULL)
  {
    fprintf(stderr, "Error: parameter -out is mandatory\n");
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
  }

  if (param_list_warn_unused(pl))
    usage (argv0, NULL, pl);
  param_list_print_command_line (stdout, pl);

  /* read ell from command line (assuming radix 10) */
  mpz_init_set_str(ell, group_order, 10);

  fprintf(stderr, "\nSub-group order:\n\tell = ");
  mpz_out_str(stderr, 10, ell);
  fprintf(stderr, "\n");

  t0 = seconds();

  FILE * out = fopen_maybe_compressed(outfile, "w");

  int64_t *a;
  uint64_t nrows, ncols;
  purgedfile_read_firstline (purgedfile, &nrows, &ncols);
  a = (int64_t *) malloc(nrows * sizeof(int64_t));
  fprintf(stderr, "%" PRIu64 "\n", nrows);

  /* For each rel, read the a,b-pair */
  read_data_t data = {.a = a};
  char *fic[2] = {(char *) purgedfile, NULL};
  filter_rels (fic, (filter_rels_callback_t) thread_read, &data,
          EARLYPARSE_NEED_AB_HEXA, NULL, NULL);

  FILE * ix = fopen_maybe_compressed(indexfile, "r");
  /* small_ncols isn't used here: we don't care */
  int ret = fscanf(ix, "%" SCNu64 " %" SCNu64 "\n", &nrows, &ncols);
  ASSERT(ret == 2);
  mpz_t tmp;
  mpz_init (tmp);

  gmp_fprintf(out, "%" PRIu64 "\n", nrows);
  for(uint64_t i = 0 ; i < nrows ; i++)
  {
    uint64_t nc;
    mpz_set_ui(tmp, 0);
    ret = fscanf(ix, "%" SCNu64 "", &nc);
    ASSERT_ALWAYS(ret == 1);
    for(uint64_t k = 0 ; k < nc ; k++)
    {
      unsigned int ridx;
      int32_t e;
      ret = fscanf(ix, "%x:%" SCNd32 "", &ridx, &e); 
      ASSERT_ALWAYS(ret == 2);
      ASSERT_ALWAYS(e != 0);

      int64_t t = e * a[ridx];
      if (t > 0)
        mpz_add_ui (tmp, tmp, t);
      else
        mpz_sub_ui (tmp, tmp, -t);
      mpz_mod (tmp, tmp, ell);
    }
    gmp_fprintf(out, "%Zd\n", tmp);
  }

  mpz_clear(tmp);
  fclose_maybe_compressed(ix, indexfile);
  free (a);

  fprintf(stderr, "\nCompleted in %2.2lf seconds\n", seconds() - t0);

  mpz_clear(ell);
  fclose_maybe_compressed(out, outfile);
  param_list_clear(pl);

  return 0;
}
