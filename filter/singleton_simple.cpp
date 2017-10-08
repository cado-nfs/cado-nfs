#include "cado.h"

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <string.h>

#include "macros.h"
#include "utils_with_io.h"
#include "filter_config.h"


// Some globals... laziness.
FILE *out;
weight_t *table;
uint64_t count, count_ideals = 0;
uint64_t col_min_index = 0;

/* -------------------------------------------------------------------------- */

void *
update_table(void * dummy MAYBE_UNUSED, earlyparsed_relation_ptr rel)
{
  for (weight_t i = 0; i < rel->nb; i++) {
    index_t ind = rel->primes[i].h;
    if (ind < col_min_index)
      continue;
    count_ideals += table[ind] == 0;
    table[ind]++;
    if (table[ind] == 0) {   // saturate
      table[ind] = UMAX(weight_t);
    }
  }
  return NULL;
}

void *
print_survivors(void * dummy MAYBE_UNUSED, earlyparsed_relation_ptr rel)
{
  for (weight_t i = 0; i < rel->nb; i++) {
    index_t ind = rel->primes[i].h;
    if (ind < col_min_index)
      continue;
    if (table[ind] == 1) {
      return NULL; // There is a singleton in this relation; don't print.
    }
  }
  fputs(rel->line, out);
  count++;
  return NULL;
}

/* -------------------------------------------------------------------------- */

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "col-min-index", "lower index on the considered "
          "columns (default 0)");
  param_list_decl_usage(pl, "col-max-index", "upper bound on the number of "
          "columns (must be at least the number\n"
          "                   of prime ideals in renumber table)");
  param_list_decl_usage(pl, "out", "output file");
  verbose_decl_usage(pl);
}

static void
usage (param_list pl, const char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}


/* -------------------------------------------------------------------------- */

int main (int argc, char **argv)
{
  const char *argv0 = argv[0];
  uint64_t col_max_index = 0;
  const char *outfile = NULL;

  param_list pl;

  /* read params */
  param_list_init(pl);
  declare_usage(pl);
  argv++,argc--;

  if (argc == 0)
    usage(pl, argv0);

  for ( ; argc ; ) {
    if (param_list_update_cmdline (pl, &argc, &argv)) { continue; }
    /* Since we accept file names freeform, we decide to never
     * abort on unrecognized options */
    break;
  }

  outfile = param_list_lookup_string(pl, "out");
  if (outfile == NULL) {
    usage(pl, argv0);
  }

  param_list_parse_uint64(pl, "col-min-index", &col_min_index);

  param_list_parse_uint64(pl, "col-max-index", &col_max_index);
  if (col_max_index == 0) {
    usage(pl, argv0);
  }

  table = (weight_t *)malloc(col_max_index*sizeof(weight_t));
  ASSERT_ALWAYS(table != NULL);

  filter_rels(argv,
      (filter_rels_callback_t) &update_table,
      NULL, EARLYPARSE_NEED_INDEX, NULL, NULL);

  out = fopen(outfile, "w");
  ASSERT_ALWAYS(out != NULL);
  count = 0;

  filter_rels(argv,
      (filter_rels_callback_t) &print_survivors,
      NULL, EARLYPARSE_NEED_LINE | EARLYPARSE_NEED_INDEX, NULL, NULL);

  printf("# %" PRIu64 " ideals read\n", count_ideals);
  printf("# %" PRIu64 " relations written\n", count);

  fclose(out);
  free(table);
  param_list_clear(pl);

  return EXIT_SUCCESS;
}
