#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mod_ul.c"
#include "portability.h"
#include "utils.h"

#ifdef FOR_FFS
#include "utils-ffs.h"
#endif

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "poly", "input polynomial file");
  param_list_decl_usage(pl, "renumber", "input file for renumbering table");
  param_list_decl_usage(pl, "check", "(switch) check the renumbering table");
  verbose_decl_usage(pl);
}

static void
usage (param_list pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}

int
main (int argc, char *argv[])
{
    int check = 0;
    char *argv0 = argv[0];
    cado_poly cpoly;
    renumber_t tab;

    param_list pl;
    param_list_init(pl);
    declare_usage(pl);
    param_list_configure_switch(pl, "check", &check);

    argv++, argc--;
    if (argc == 0)
      usage (pl, argv0);

    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv))
            continue;
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage (pl, argv0);
    }
    /* print command-line arguments */
    verbose_interpret_parameters(pl);
    param_list_print_command_line (stdout, pl);
    fflush(stdout);

    const char *polyfilename = param_list_lookup_string(pl, "poly");
    const char *renumberfilename = param_list_lookup_string(pl, "renumber");

    if (polyfilename == NULL)
    {
      fprintf (stderr, "Error, missing -poly command line argument\n");
      usage (pl, argv0);
    }
    if (renumberfilename == NULL)
    {
      fprintf (stderr, "Error, missing -renumber command line argument\n");
      usage (pl, argv0);
    }

    cado_poly_init(cpoly);
#ifndef FOR_FFS
    if (!cado_poly_read (cpoly, polyfilename))
#else
    if (!ffs_poly_read (cpoly, polyfilename))
#endif
    {
      fprintf (stderr, "Error reading polynomial file\n");
      exit (EXIT_FAILURE);
    }

  renumber_init_for_reading (tab);
  renumber_read_table (tab, renumberfilename);

  for (index_t i = 0; i < tab->size; i++)
  {
    if (tab->table[i] == RENUMBER_SPECIAL_VALUE)
    {
      int side = renumber_get_side_from_index (tab, i, cpoly);
      if (renumber_is_additional_column (tab, i))
        printf ("i=%" PRid " tab[i]=#   added column for side %d\n", i, side);
      else
        printf ("i=%" PRid " tab[i]=#   above a bad ideals on side %d\n",
                i, side);
    }
    else
    {
      p_r_values_t p, r;
      int side;

      renumber_get_p_r_from_index (tab, &p, &r, &side, i, cpoly);
      printf ("i=%" PRid " tab[i]=%" PRpr " p=%" PRpr "", i, tab->table[i], p);
      if (tab->rat == -1)
      {
        if (r == p)
          printf (" r=%" PRpr " side %d proj\n", r, side);
        else
          printf (" r=%" PRpr " side %d\n", r, side);
      }
      else
      {
        if (side == tab->rat)
          printf (" rat side\n");
        else if (r == p)
          printf (" r=%" PRpr " alg side proj\n", r);
        else
          printf (" r=%" PRpr " alg side\n", r);
      }
    }
  }

  if (tab->bad_ideals.n != 0) {
    printf ("Bad ideals:\n");
    for (int i = 0; i < tab->bad_ideals.n; ++i) {
      p_r_values_t p = tab->bad_ideals.p[i];
      p_r_values_t r = tab->bad_ideals.r[i];
      int nb = tab->bad_ideals.nb[i];
      printf("p=%" PRpr " r=%" PRpr " nb=%d\n", p, r, nb);
    }
  }

  /* Check for all index i if it can retrieved p and r from it and if it can
   * retrieved the same index from this p and r
   */
  uint64_t nerrors = 0;
  if (check)
  {

    index_t i, j;
    p_r_values_t p, r;
    int side;
    for (i = 0; i < tab->size; i++)
    {
      if (tab->table[i] != RENUMBER_SPECIAL_VALUE)
      {
        renumber_get_p_r_from_index(tab, &p, &r, &side, i, cpoly);
        j = renumber_get_index_from_p_r (tab, p, r, side);
        if (i == j)
          fprintf (stderr, "## %" PRid ": Ok\n", i);
        else
        {
          fprintf (stderr, "#### %" PRid ": Error, p = %" PRpr " r=%" PRpr " "
                           "give index j = %" PRid "\n", i, p, r, j);
          nerrors++;
        }
      }
    }
    fprintf (stderr, "Number of errors: %" PRIu64 "\n", nerrors);
  }

  renumber_clear (tab);
  cado_poly_clear (cpoly);
  param_list_clear(pl);
  return (nerrors == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
