#include "cado.h" // IWYU pragma: keep
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cinttypes>   // for PRIu64
#include <cstdint>     // for uint64_t
#include "filter_io.hpp"  // earlyparsed_relation_ptr
#include "macros.h"
#include "misc.h"       // for UMAX
#include "params.h"     // param_list_parse_*
#include "typedefs.h"   // for weight_t, prime_t, index_t

#include "verbose.hpp"


// Some globals... laziness.
static FILE *out;
static weight_t *table;
static uint64_t count, count_ideals = 0;
static uint64_t col_min_index = 0;

/* -------------------------------------------------------------------------- */

static void *
update_table(void * dummy MAYBE_UNUSED, earlyparsed_relation_ptr rel)
{
  for (weight_t i = 0; i < rel->nb; i++) {
    index_t const ind = rel->primes[i].h;
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

static void *
print_survivors(void * dummy MAYBE_UNUSED, earlyparsed_relation_ptr rel)
{
  for (weight_t i = 0; i < rel->nb; i++) {
    index_t const ind = rel->primes[i].h;
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

static void declare_usage(cxx_param_list & pl)
{
  param_list_decl_usage(pl, "col-min-index", "lower index on the considered "
          "columns (default 0)");
  param_list_decl_usage(pl, "col-max-index", "upper bound on the number of "
          "columns (must be at least the number\n"
          "                   of prime ideals in renumber table)");
  param_list_decl_usage(pl, "out", "output file");
  verbose_decl_usage(pl);
}

/* -------------------------------------------------------------------------- */

// coverity[root_function]
int main (int argc, char const **argv)
{
  uint64_t col_max_index = 0;
  const char *outfile = NULL;

  cxx_param_list pl;

  declare_usage(pl);


  param_list_process_command_line(pl, &argc, &argv, true);

  outfile = param_list_lookup_string(pl, "out");
  if (outfile == NULL)
      pl.fail("missing argument -out");

  param_list_parse_uint64(pl, "col-min-index", &col_min_index);

  param_list_parse_uint64(pl, "col-max-index", &col_max_index);
  if (col_max_index == 0)
      pl.fail("missing argument col-max-index");

  table = (weight_t *)malloc(col_max_index*sizeof(weight_t));
  ASSERT_ALWAYS(table != NULL);
  memset(table, 0, col_max_index*sizeof(weight_t));

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

  return EXIT_SUCCESS;
}
