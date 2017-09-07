#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for strcmp() */
#include <math.h> /* for sqrt and floor and log and ceil */
#include <pthread.h>
#include "portability.h"
#include "utils.h"
#include <vector>
#include <algorithm>
#include "cxx_mpz.hpp"

/*
 * The goal of this binary is to produce relations that try to be good
 * approximations of what las would produce, with respect to filtering.
 * The idea is to feed the purge/merge steps of the filtering with these
 * relations to get an idea of the final matrix that will come out of the
 * sieving step for a given set of parameters.
 *
 * The shell script cado-nfs/misc/estimate_matsize.sh is an attempt to
 * run all the required steps for such a simulation: sampling with las,
 * fake relation generation with this binary, and filter.
 *
 * The binary takes as input:
 *   - poly file
 *   - lpb on each side
 *   - a range of special-q for a given side
 *   - a sample of relations (output of las) for this range
 *   - the renumber table.
 * The renumber table is required, because this binary will produce
 * relations as if they were coming out of dup2 (hence renumbered).
 *
 * The sample of relations must also be de-duplicated, therefore the -dup
 * option of las must be used for the sampling.
 *
 */

using namespace std;

// A global mutex for I/O
pthread_mutex_t io_mutex = PTHREAD_MUTEX_INITIALIZER;

// We need a re-entrant random. Let's take GMP.
static inline uint64_t long_random(gmp_randstate_t buf) {
#if ULONG_BITS == 64
    return gmp_urandomb_ui(buf, 64);
#elif ULONG_BITS == 32
    cxx_mpz z;
    mpz_urandomb(z, buf, 64);
    return mpz_get_uint64(z);
#endif
}

gmp_randstate_t global_rstate_non_mt;
long myrandom_non_mt() {
    return long_random(global_rstate_non_mt);
}

// Structure that contains indices (in the sense of renumber.[ch]) for
// one side, with the possibility to get a randomized ideal on the same
// side, and more or less of the same size.
// Usage:
//  - init()
//  - do many append() (increasing order of indices)
//  - finalize()
//  - can call random_index()
//  - can call iterator_from_index()

struct indexrange {
    vector <index_t>ind;

    void init() { }
    void finalize() { }

    void append(index_t z) {
        ind.push_back(z);
    }

    index_t random_index(index_t z, gmp_randstate_t buf) {
        // Find position of z
        // Exact version:
        //  vector<index_t>::iterator it;
        //  it = lower_bound (ind.begin(), ind.end(), z);
        //  uint64_t position = it - ind.begin();
        // Fast, approximate version:
        uint64_t position = z >> 1; // half of indices on each side
        uint64_t range = uint64_t((double(position))*0.2);
        uint64_t low = MAX((int64_t)(position - range), 0);
        uint64_t high = MIN(position + range, ind.size());
        if (high == low)
            high++;
        return ind[low + uint64_t(long_random(buf)%(high-low))];
    }

    vector <index_t>::iterator iterator_from_index(index_t z) {
        vector<index_t>::iterator it;
        it = lower_bound(ind.begin(), ind.end(), z);
        return it;
    }
};

// Fill in the indexrange data structure from the renumber table,
// gathering indices in two arrays, one for each side.
void prepare_indexrange(indexrange *Ind, renumber_t ren_tab, 
        cado_poly cpoly) {
    Ind[0].init();
    Ind[1].init();
    for (index_t i = 0; i < ren_tab->size; i++) {
        int side = renumber_get_side_from_index(ren_tab, i, cpoly);
        Ind[side].append(i);
    }
    Ind[0].finalize();
    Ind[1].finalize();
}

#define MAXFACTORS 50  // on each side. Should be enough?

struct fake_rel {
    index_t ind[2][MAXFACTORS];
    int nb_ind[2];
};

void read_rel(fake_rel& rel, uint64_t q, int sqside, const char *str,
        renumber_t ren_tab) {
    rel.nb_ind[0] = 0;
    rel.nb_ind[1] = 0;
    int side = 0;
    // read a and b
    int64_t a;
    uint64_t b;
    const char *pstr = str;
    {
        char *endpstr = NULL;
        a = strtoll(pstr, &endpstr, 10);
        ASSERT_ALWAYS (endpstr != pstr);
        pstr = endpstr;
        ASSERT_ALWAYS (pstr[0]==',');
        pstr++;
        b = strtoull(pstr, &endpstr, 10);
        ASSERT_ALWAYS (endpstr != pstr);
        pstr = endpstr;
        // skip ':'
        ASSERT_ALWAYS (pstr[0]==':');
        pstr++;
    }
    while (pstr[0] != '\0' && pstr[0] != '\n') {
        if (pstr[0] == ':') {
            // change side
            side++;
            pstr++;
            continue;
        }
        if (pstr[0] == ',') {
            pstr++;
            continue;
        }
        char *endpstr = NULL;
        uint64_t p = strtoull(pstr, &endpstr, 16);
        ASSERT_ALWAYS (endpstr != pstr);
        if (side != sqside || q != p) {
            p_r_values_t r = relation_compute_r(a, b, p);
            index_t index;
            int nb;
            if (renumber_is_bad (&nb, &index, ren_tab, p, r, side)) {
                // bad ideal: just pick a random ideal above this prime
                index += (myrandom_non_mt() % nb);
            } else {
                index = renumber_get_index_from_p_r(ren_tab, p, r, side);
            }
            rel.ind[side][rel.nb_ind[side]] = index;
            rel.nb_ind[side]++;
            ASSERT_ALWAYS (rel.nb_ind[side] <= MAXFACTORS);
        }
        pstr = endpstr;
    }
}

void read_sample_file(vector<unsigned int> &nrels, vector<fake_rel> &rels,
        int sqside, const char *filename, renumber_t ren_tab)
{
    FILE * file;
    file = fopen(filename, "r");
    ASSERT_ALWAYS (file != NULL);
    uint64_t q = 0;
    unsigned int nr = 0;

    char line[1024];

    do {
        char * ret = fgets(line, 1024, file);
        if (ret == NULL) {
            nrels.push_back(nr);
            break;
        }
        if (line[0] == '#') {
            char *ptr;
            ptr = strstr(line, "# Sieving side-");
            if (ptr != NULL) {
                // starting a new special-q
                ASSERT_ALWAYS (sqside == (ptr[15] - '0'));
                ptr += 19; // skip "# Sieving side-0 q="
                q = strtoull(ptr, NULL, 10);
                nrels.push_back(nr);
                nr = 0;
            }
        } else {
            fake_rel rel;
            read_rel(rel, q, sqside, line, ren_tab);
            rels.push_back(rel);
            nr++;
        }
    } while(1);

    fclose(file);
}

static
int index_cmp(const void *p1, const void *p2)
{
    index_t pp1 = ((index_t *)(p1))[0];
    index_t pp2 = ((index_t *)(p2))[0];
    if (pp1 == pp2)
        return 0;
    if (pp1 < pp2)
        return -1;
    return 1;
}

void reduce_mod_2(index_t *frel, int *nf) {
    int i = 1;
    index_t curr = frel[0];
    int nb = 1;
    int j = 0;
    while (i < *nf) {
        if (frel[i] != curr) {
            if (nb & 1) {
                frel[j++] = curr;
            }
            curr = frel[i];
            nb = 1;
        } else {
            nb++;
        }
        i++;
    }
    if (nb & 1) {
        frel[j++] = curr;
    }
    *nf = j;
}

void print_fake_rel_manyq(
        vector<index_t>::iterator list_q, uint64_t nq,
        vector<fake_rel> *rels, vector<unsigned int> *nrels,
        indexrange *Ind, int dl,
        gmp_randstate_t buf)
{
    index_t frel[MAXFACTORS];
#define MAX_STR 256  // string containing a printable relation.
    char str[MAX_STR];
    char *pstr;
    int len = MAX_STR;
    for (uint64_t ii = 0; ii < nq; ++ii) {
        index_t indq = list_q[ii];
        int nr = int((*nrels)[long_random(buf)%nrels->size()]);
        for (; nr > 0; --nr) {
            pstr = str;
            len = MAX_STR;
            // pick fake a,b
            int nc = snprintf(pstr, len, "%ld,%lu:",
                    (long)long_random(buf), (unsigned long)long_random(buf));
            pstr += nc;
            len -= nc;
            // pick a random relation as a model and prepare a string to
            // be printed.
            int i = long_random(buf) % rels->size();
            frel[0] = indq;
            int nf = 1;
            for (int side = 0; side < 2; ++side) {
                int np = (*rels)[i].nb_ind[side];
                for (int j = 0; j < np; ++j) {
                    index_t ind = Ind[side].random_index((*rels)[i].ind[side][j], buf);
                    frel[nf] = ind;
                    nf++;
                    ASSERT(nf < MAXFACTORS);
                }
            }
            qsort(frel, nf, sizeof(index_t), index_cmp);
            if (!dl) {
                reduce_mod_2(frel, &nf); // update nf
            }
            for (int i = 0; i < nf; ++i) {
                nc = snprintf(pstr, len, "%" PRid, frel[i]);
                pstr += nc;
                len -= nc;
                if (i != nf-1) {
                    snprintf(pstr, len, ",");
                    pstr++; len--;
                }

            }
            snprintf(pstr, len, "\n");
            // Get the mutex and print
            pthread_mutex_lock(&io_mutex);
            printf("%s", str);
            pthread_mutex_unlock(&io_mutex);
        }
    }
#undef MAX_STR
}

struct th_args {
    vector<index_t>::iterator list_q;
    uint64_t nq;
    vector<fake_rel> *rels;
    vector<unsigned int> *nrels;
    indexrange *Ind;
    int dl;
    gmp_randstate_t rstate;
};


void * do_thread(void * rgs) {
    struct th_args * args = (struct th_args *) rgs;
    print_fake_rel_manyq(args->list_q, args->nq, args->rels, args->nrels,
            args->Ind, args->dl, args->rstate);
    return NULL;
}

void advance_prime_in_fb(int *mult, uint64_t *q, uint64_t *roots,
        cado_poly cpoly, int sqside, prime_info pdata)
{
    int nr;
    unsigned long newp;
    do {
        newp = getprime_mt (pdata);
        ASSERT_ALWAYS (newp > *q);
        nr = mpz_poly_roots_uint64(roots, cpoly->pols[sqside], newp);
    } while (nr == 0);

    *mult = nr;
    *q = newp;
}

static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "poly", "polynomial file");
    param_list_decl_usage(pl, "lpb0", "factor base bound on side 0");
    param_list_decl_usage(pl, "lpb1", "factor base bound on side 1");
    param_list_decl_usage(pl, "q0", "lower bound of the qrange");
    param_list_decl_usage(pl, "q1", "upper bound of the qrange");
    param_list_decl_usage(pl, "sqside", "side of the special-q");
    param_list_decl_usage(pl, "sample", "file where to find a sample of relations");
    param_list_decl_usage(pl, "renumber", "renumber table");
    param_list_decl_usage(pl, "dl", "(switch) dl mode");
    param_list_decl_usage(pl, "t", "number of threads to use");
    verbose_decl_usage(pl);
}

int
main (int argc, char *argv[])
{
  param_list pl;
  cado_poly cpoly;
  int sqside = -1;
  char *argv0 = argv[0];
  int lpb[2] = {0, 0};
  uint64_t q0 = 0;
  uint64_t q1 = 0;
  int dl = 0;
  int mt = 1;

  param_list_init(pl);
  declare_usage(pl);
  param_list_configure_switch(pl, "-dl", &dl);
  
  cado_poly_init(cpoly);

  argv++, argc--;
  for( ; argc ; ) {
      if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }

      /* Could also be a file */
      FILE *f;
      if ((f = fopen(argv[0], "r")) != NULL) {
          param_list_read_stream(pl, f, 0);
          fclose(f);
          argv++,argc--;
          continue;
      }

      fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
      param_list_print_usage(pl, argv0, stderr);
      exit (EXIT_FAILURE);
  }
  verbose_interpret_parameters(pl);
  param_list_print_command_line(stdout, pl);

  const char * filename;
  if ((filename = param_list_lookup_string(pl, "poly")) == NULL) {
      fprintf(stderr, "Error: parameter -poly is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_int(pl, "lpb0", &lpb[0]);
  if (lpb[0] == 0) {
      fprintf(stderr, "Error: parameter -lpb0 is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }
  param_list_parse_int(pl, "lpb1", &lpb[1]);
  if (lpb[1] == 0) {
      fprintf(stderr, "Error: parameter -lpb1 is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_uint64(pl, "q0", &q0);
  if (q0 == 0) {
      fprintf(stderr, "Error: parameter -q0 is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_uint64(pl, "q1", &q1);
  if (q1 == 0) {
      fprintf(stderr, "Error: parameter -q1 is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_int(pl, "t", &mt);

  if (!cado_poly_read(cpoly, filename))
    {
      fprintf (stderr, "Error reading polynomial file %s\n", filename);
      exit (EXIT_FAILURE);
    }

  param_list_parse_int(pl, "sqside", &sqside);
  if (sqside == -1 || sqside > 2) {
      fprintf(stderr, "Error: sqside must be 0 or 1\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  const char * renumberfile;
  if ((renumberfile = param_list_lookup_string(pl, "renumber")) == NULL) {
      fprintf(stderr, "Error: parameter -renumber is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }
  renumber_t ren_table;
  renumber_read_table(ren_table, renumberfile);

  // read sample file
  const char * sample;
  if ((sample = param_list_lookup_string(pl, "sample")) == NULL) {
      fprintf(stderr, "Error: parameter -sample is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  gmp_randinit_default(global_rstate_non_mt);
  gmp_randseed_ui(global_rstate_non_mt, 0);


  vector<fake_rel> rels;
  vector<unsigned int> nrels;
  read_sample_file(nrels, rels, sqside, sample, ren_table);

  param_list_warn_unused(pl);

  // Two index ranges, one for each side
  indexrange Ind[2];
  prepare_indexrange(Ind, ren_table, cpoly);

  // Precompute ideals in [q0-q1]
  prime_info pdata;
  prime_info_init(pdata);
  // fast forward until we reach q0
  uint64_t q = 2;
  while (q < q0) {
      q = getprime_mt(pdata);
  }
  uint64_t roots[MAX_DEGREE];
  int mult = mpz_poly_roots_uint64(roots, cpoly->pols[sqside], q); 
  while (mult == 0) {
      advance_prime_in_fb(&mult, &q, roots, cpoly, sqside, pdata);
  }
  index_t indq = renumber_get_index_from_p_r(ren_table, q, roots[0], sqside);
  vector<index_t>::iterator first_indq = Ind[sqside].iterator_from_index(indq);
  // Same for q1:
  while (q < q1) {
      q = getprime_mt(pdata);
  }
  mult = mpz_poly_roots_uint64(roots, cpoly->pols[sqside], q); 
  while (mult == 0) {
      advance_prime_in_fb(&mult, &q, roots, cpoly, sqside, pdata);
  }
  indq = renumber_get_index_from_p_r(ren_table, q, roots[0], sqside);
  vector<index_t>::iterator last_indq = Ind[sqside].iterator_from_index(indq);

  // go multi-thread
  uint64_t block = (last_indq - first_indq) / mt;
  pthread_t * thid = (pthread_t *)malloc(mt*sizeof(pthread_t));
  struct th_args * args = (struct th_args *)malloc(mt*sizeof(struct th_args));
  for (int i = 0; i < mt; ++i) {
      args[i].list_q = first_indq + block*i;
      if (i < mt-1) {
          args[i].nq = block;
      } else {
          args[i].nq = last_indq - args[i].list_q;
      }
      args[i].rels = &rels;
      args[i].nrels = &nrels;
      args[i].Ind = &Ind[0];
      args[i].dl = dl;
      gmp_randinit_default(args[i].rstate);
      gmp_randseed_ui(args[i].rstate, 171717+i);

      pthread_create(&thid[i], NULL, do_thread, (void *)(&args[i]));
  }
  for (int i = 0; i < mt; ++i) {
      pthread_join(thid[i], NULL);
      gmp_randclear(args[i].rstate);
  }

  gmp_randclear(global_rstate_non_mt);
  free(thid);
  free(args);
  prime_info_clear(pdata);
  renumber_clear(ren_table);
  cado_poly_clear(cpoly);
  param_list_clear(pl);

  return 0;
}
