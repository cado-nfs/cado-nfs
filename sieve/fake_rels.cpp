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
 *     NOTE: las should be run with the -v option
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

// A global variable for the number of relations printed
unsigned long rels_printed = 0;

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
//  - do many append() (increasing order of indices) and append_prime()
//  - finalize()
//  - can call random_index()
//  - can call iterator_from_index()
//  - can call pos_from_p() and p_from_pos() if append_prime() has been done

struct indexrange {
    vector <index_t>ind;
    vector <p_r_values_t>prime;

    void init() { }
    void finalize() { }

    void append(index_t z) {
        ind.push_back(z);
    }
    void append_prime(p_r_values_t p) {
        prime.push_back(p);
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

    p_r_values_t p_from_pos(uint64_t position) {
        return prime[position];
    }

    // Compute one position corresponding to the given p. In the case where
    // there are several positions, the choice is arbitrary.
    // If p is not in the table, compute a position with a close-enough
    // prime.
    uint64_t pos_from_p(p_r_values_t p) {
        uint64_t low, high;
        p_r_values_t pp;
        low = 0; high = prime.size()-1;
        pp = prime[low];
        if (pp == p)
            return low;
        pp = prime[high];
        if (pp == p)
            return high;
        do {
            uint64_t middle = (low+high)>>1;
            pp = prime[middle];
            if (p == pp)
                return middle;
            if (p < pp)
                high = middle;
            else
                low = middle;
        } while (high > low+1);
        return low;
    }
 
    vector <index_t>::iterator iterator_from_index(index_t z) {
        vector<index_t>::iterator it;
        it = lower_bound(ind.begin(), ind.end(), z);
        return it;
    }
};

// Fill in the indexrange data structure from the renumber table,
// gathering indices in two arrays, one for each side.
// In case of composite special-qs, also fill-in the list of the
// corresponding primes on the sqside.
void prepare_indexrange(indexrange *Ind, renumber_t ren_tab, 
        cado_poly cpoly, int sqside, int compsq) {
    Ind[0].init();
    Ind[1].init();
    for (index_t i = 0; i < ren_tab->size; i++) {
        if (renumber_is_additional_column(ren_tab, i)) {
            continue;
        }
        int side = renumber_get_side_from_index(ren_tab, i, cpoly);
        Ind[side].append(i);
        if (compsq && (side == sqside)) {
            p_r_values_t p, r;
            if (ren_tab->table[i] == RENUMBER_SPECIAL_VALUE) {
                renumber_badideal_get_p_r_below(ren_tab, &p, &r, &side, i);
            } else {
                renumber_get_p_r_from_index(ren_tab, &p, &r, &side, i, cpoly);
            }
            Ind[sqside].append_prime(p);
        }
    }
    Ind[0].finalize();
    Ind[1].finalize();
}

#define MAXFACTORS 50  // on each side. Should be enough?

struct fake_rel {
    index_t ind[2][MAXFACTORS];
    int nb_ind[2];
    /* only used while reading, so that we get a stable sort */
    int64_t a;
    uint64_t b;
};

bool operator<(fake_rel const& x, fake_rel const& y) {
    if (x.a < y.a) return true;
    if (x.a > y.a) return false;
    return x.b < y.b;
}

int p_coprimeto_q(uint64_t p, uint64_t q, vector<uint64_t> facq) {
    if (facq.size() == 0) {
        return q != p;
    } else {
        for (auto f : facq) {
            if (f == p)
                return 0;
        }
        return 1;
    }
}


void read_rel(fake_rel& rel, uint64_t q, vector<uint64_t> facq,
        int sqside, const char *str, renumber_t ren_tab) {
    rel.nb_ind[0] = 0;
    rel.nb_ind[1] = 0;
    int side = 0;
    // read a and b
    int64_t a;
    uint64_t b;
    const char *pstr = str;
    {
        char *endpstr = NULL;
        rel.a = a = strtoll(pstr, &endpstr, 10);
        ASSERT_ALWAYS (endpstr != pstr);
        pstr = endpstr;
        ASSERT_ALWAYS (pstr[0]==',');
        pstr++;
        rel.b = b = strtoull(pstr, &endpstr, 10);
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
        if (side != sqside || p_coprimeto_q(p, q, facq)) {
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
        int sqside, const char *filename, renumber_t ren_tab, int compsq)
{
    FILE * file;
    file = fopen_maybe_compressed(filename, "r");
    ASSERT_ALWAYS (file != NULL);
    uint64_t q = 0;
    vector<uint64_t> facq;
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
                uint64_t nq = strtoull(ptr, &ptr, 10);
                if (compsq) {
                    // read also the factorization of q
                    facq.clear();
                    ASSERT_ALWAYS(ptr[0] == '=');
                    do {
                        ptr++;
                        uint64_t fac = strtoull(ptr, &ptr, 10);
                        facq.push_back(fac);
                    } while (ptr[0] == '*');
                }
                if (nq != q) {
                    nrels.push_back(nr);
                    nr = 0;
                    q = nq;
                }
            }
        } else {
            fake_rel rel;
            read_rel(rel, q, facq, sqside, line, ren_tab);
            rels.push_back(rel);
            nr++;
        }
    } while(1);

    size_t s = 0;
    for(auto const & n : nrels) {
        ASSERT_ALWAYS((s + n) <= rels.size());
        std::sort(rels.begin() + s, rels.begin() + s + n);
        s += n;
    }

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

void shrink_indices(index_t *frel, int nf, int shrink_factor) {
    // Indices below this threshold are not shrinked
    // FIXME: I am not sure we should keep the heavy weight columns
    // un-shrinked. The answer might be different in DL and in facto...
    const index_t noshrink_threshold = 0;
    for (int i = 0; i < nf; ++i) {
//        if (frel[i] >= noshrink_threshold) {
        if (1) {
            frel[i] = noshrink_threshold +
                (frel[i] - noshrink_threshold) / shrink_factor;
        }
    }
}

void print_fake_rel_manyq(
        vector<index_t>::iterator list_q, int nfacq, uint64_t nq,
        vector<fake_rel> *rels, vector<unsigned int> *nrels,
        indexrange *Ind, int dl, int shrink_factor,
        gmp_randstate_t buf)
{
    index_t frel[2*MAXFACTORS];
#define BUF_SIZE 100 // buffer containing several relations
#define MAX_STR 256  // string containing a printable relation.
    char str[BUF_SIZE*MAX_STR];
    char *pstr;
    int len = MAX_STR;
    int size = 0; // number of relations in buffer
    unsigned long nrels_thread = 0;
    pstr = str; // initialize buffer
    for (uint64_t ii = 0; ii < nfacq*nq; ii += nfacq) {
        int nr;
        if (shrink_factor == 1) {
            nr = int((*nrels)[long_random(buf)%nrels->size()]);
        } else {
            double nr_dble = double((*nrels)[long_random(buf)%nrels->size()])
                / double(shrink_factor);
            // Do probabilistic rounding, in case nr_dble is small (maybe < 1)
            double trunc_part = trunc(nr_dble);
            double frac_part = nr_dble - trunc_part;
            double rnd = double(long_random(buf)) / double(UINT64_MAX);
            nr = int(trunc_part) + int(rnd < frac_part);
        }
        //        fprintf(stdout, "%u, %u\n", list_q[ii], list_q[ii+1]);
	nrels_thread += nr; /* we will output nr fake relations */
        for (; nr > 0; --nr) {
            len = MAX_STR;
            // pick fake a,b
            int nc = snprintf(pstr, len, "%lx,%lx:",
                    (unsigned long)long_random(buf),
                    (unsigned long)long_random(buf));
            pstr += nc;
            len -= nc;
            // pick a random relation as a model and prepare a string to
            // be printed.
            int i = long_random(buf) % rels->size();
            int nf = 0;
            for (int j = 0; j < nfacq; j++) {
                frel[nf] = list_q[ii + j];
                nf++;
            }
            for (int side = 0; side < 2; ++side) {
                int np = (*rels)[i].nb_ind[side];
                for (int j = 0; j < np; ++j) {
                    index_t ind = Ind[side].random_index((*rels)[i].ind[side][j], buf);
                    ASSERT(nf < 2*MAXFACTORS);
                    frel[nf] = ind;
                    nf++;
                }
            }
            qsort(frel, nf, sizeof(index_t), index_cmp);
            if (!dl) {
                reduce_mod_2(frel, &nf); // update nf
            }
            if (shrink_factor > 1) {
                shrink_indices(frel, nf, shrink_factor);
                if (!dl) {
                    reduce_mod_2(frel, &nf);
                }
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
	    pstr++;
	    size++;
	    if (size == BUF_SIZE)
	      {
		// Get the mutex and print
		pthread_mutex_lock(&io_mutex);
		printf("%s", str);
		pthread_mutex_unlock(&io_mutex);
		pstr = str;
		size = 0;
	      }
        }
    }
    pthread_mutex_lock(&io_mutex);
    if (size > 0) // print leftover relations if any
      printf("%s", str);
    rels_printed += nrels_thread;
    pthread_mutex_unlock(&io_mutex);
#undef MAX_STR
#undef BUF_SIZE
}

struct th_args {
    vector<index_t>::iterator list_q_prime;
    uint64_t nq;
    vector<index_t>::iterator list_q_comp2;
    uint64_t nq2;
    vector<index_t>::iterator list_q_comp3;
    uint64_t nq3;
    vector<fake_rel> *rels;
    vector<unsigned int> *nrels;
    indexrange *Ind;
    int dl;
    int shrink_factor;
    gmp_randstate_t rstate;
};


void * do_thread(void * rgs) {
    struct th_args * args = (struct th_args *) rgs;
    if (args->nq > 0)
        print_fake_rel_manyq(args->list_q_prime, 1, args->nq, args->rels,
                args->nrels, args->Ind, args->dl, args->shrink_factor,
                args->rstate);
    
    if (args->nq2 > 0)
        print_fake_rel_manyq(args->list_q_comp2, 2, args->nq2, args->rels,
                args->nrels, args->Ind, args->dl, args->shrink_factor,
                args->rstate);

    if (args->nq3 > 0)
        print_fake_rel_manyq(args->list_q_comp3, 3, args->nq3, args->rels,
                args->nrels, args->Ind, args->dl, args->shrink_factor,
                args->rstate);

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

// All composite special-q in [q0, q1] with 2 prime factors in the range
// [qfac_min, qfac_max] in Ind.
vector<index_t> all_comp_sq_2(uint64_t q0, uint64_t q1, uint64_t qfac_min,
        uint64_t qfac_max, indexrange &Ind)
{
    vector<index_t> list;

    uint64_t l1min = MAX(q0/qfac_max, qfac_min);
    uint64_t pos_l1min = Ind.pos_from_p(l1min);
    uint64_t l1max = MIN(qfac_max, round(sqrt(q1)));
    uint64_t pos_l1max = Ind.pos_from_p(l1max);

    for (uint64_t pos1 = pos_l1min; pos1 < pos_l1max; ++pos1) {
        uint64_t l1 = Ind.p_from_pos(pos1);
        uint64_t l2min = MAX(l1, q0/l1);
        uint64_t pos_l2min = Ind.pos_from_p(l2min);
        uint64_t l2max = MIN(qfac_max, q1/l1);
        uint64_t pos_l2max = Ind.pos_from_p(l2max);
        for (uint64_t pos2 = pos_l2min; pos2 < pos_l2max; ++pos2) {
            uint64_t l2 = Ind.p_from_pos(pos2);
            uint64_t L = l1*l2;
            if (L >= q0 && L <= q1) {
                list.push_back(Ind.ind[pos1]);
                list.push_back(Ind.ind[pos2]);
//                fprintf(stderr, "%lu*%lu : (%u, %u)\n",
//                        l1, l2, Ind.ind[pos1], Ind.ind[pos2]);
            } else {
//                fprintf(stderr, "Wooops !!!\n");
            }
        }
    }
    fprintf(stderr, "Got %zu 2-composite sq\n", (size_t)(list.size()>>1));
    return list;
}

// All composite special-q in [q0, q1] with 3 prime factors in the range
// [qfac_min, qfac_max] in Ind.
vector<index_t> all_comp_sq_3(uint64_t q0, uint64_t q1, uint64_t qfac_min,
        uint64_t qfac_max, indexrange &Ind)
{
    vector<index_t> list;

    uint64_t l1min = MAX(q0/(qfac_max*qfac_max), qfac_min);
    uint64_t pos_l1min = Ind.pos_from_p(l1min);
    uint64_t l1max = MIN(qfac_max, round(pow(q1, 0.333333333333)));
    uint64_t pos_l1max = Ind.pos_from_p(l1max);

    for (uint64_t pos1 = pos_l1min; pos1 < pos_l1max; ++pos1) {
        uint64_t l1 = Ind.p_from_pos(pos1);
        uint64_t l2min = MAX(l1, q0/(l1*qfac_max));
        uint64_t pos_l2min = Ind.pos_from_p(l2min);
        uint64_t l2max = MIN(qfac_max, round(sqrt(q1/l1)));
        uint64_t pos_l2max = Ind.pos_from_p(l2max);
        for (uint64_t pos2 = pos_l2min; pos2 < pos_l2max; ++pos2) {
            uint64_t l2 = Ind.p_from_pos(pos2);
            uint64_t l3min = MAX(l2, q0/(l1*l2));
            uint64_t pos_l3min = Ind.pos_from_p(l3min);
            uint64_t l3max = MIN(qfac_max, q1/(l1*l2));
            uint64_t pos_l3max = Ind.pos_from_p(l3max);
            for (uint64_t pos3 = pos_l3min; pos3 < pos_l3max; ++pos3) {
                uint64_t l3 = Ind.p_from_pos(pos3);
                uint64_t L = l1*l2*l3;
                if (L >= q0 && L <= q1) {
                    list.push_back(Ind.ind[pos1]);
                    list.push_back(Ind.ind[pos2]);
                    list.push_back(Ind.ind[pos3]);
                } else {
              //      fprintf(stderr, "Wooops !!!\n");
                }
            }
        }
    }
    fprintf(stderr, "Got %zu 3-composite sq\n", (size_t)(list.size()/3));
    return list;
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
    param_list_decl_usage(pl, "shrink-factor", "simulate with a matrix that number (integer >= 1) times smaller");
    param_list_decl_usage(pl, "dl", "dl mode");
    param_list_decl_usage(pl, "allow-compsq", "(switch) allows composite sq");
    param_list_decl_usage(pl, "qfac-min", "factors of q must be at least that");
    param_list_decl_usage(pl, "qfac-max", "factors of q must be at most that");
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
  int compsq = 0;
  uint64_t qfac_min = 1024;
  uint64_t qfac_max = UINT64_MAX;
  int shrink_factor = 1; // by default, no shrink

  param_list_init(pl);
  declare_usage(pl);
  param_list_configure_switch(pl, "-dl", &dl);
  param_list_configure_switch(pl, "-allow-compsq", &compsq);
  
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
  
  param_list_parse_int(pl, "shrink-factor", &shrink_factor);
  if (shrink_factor < 1) {
      fprintf(stderr, "Error: shrink factor must be an integer >= 1\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_int(pl, "t", &mt);
  param_list_parse_uint64(pl, "qfac-min", &qfac_min);
  param_list_parse_uint64(pl, "qfac-max", &qfac_max);

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
  printf ("# Start reading renumber table\n");
  fflush (stdout);
  renumber_read_table(ren_table, renumberfile);
  printf ("# Done reading renumber table\n");
  fflush (stdout);

  for (int side = 0; side < 2; ++side) {
      if (ren_table->lpb[side] != (unsigned long)lpb[side]) {
          fprintf(stderr, "Error: on side %d, lpb on the command-line is different from the one in the renumber file\n", side);
          exit(EXIT_FAILURE);
      }
  }

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
  printf ("# Start reading sample file\n");
  fflush (stdout);
  read_sample_file(nrels, rels, sqside, sample, ren_table, compsq);
  printf ("# Done reading sample file\n");
  fflush (stdout);

  param_list_warn_unused(pl);

  // Two index ranges, one for each side
  indexrange Ind[2];
  printf ("# Start preparing index ranges\n");
  fflush (stdout);
  prepare_indexrange(Ind, ren_table, cpoly, sqside, compsq);
  printf ("# Done preparing index ranges\n");
  fflush (stdout);


  /****** Prime special-q ******/
  vector<index_t>::iterator first_indq, last_indq;
  if (!compsq) {
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
      first_indq = Ind[sqside].iterator_from_index(indq);
      // Same for q1:
      while (q < q1) {
          q = getprime_mt(pdata);
      }
      mult = mpz_poly_roots_uint64(roots, cpoly->pols[sqside], q); 
      while (mult == 0) {
          advance_prime_in_fb(&mult, &q, roots, cpoly, sqside, pdata);
      }
      indq = renumber_get_index_from_p_r(ren_table, q, roots[0], sqside);
      last_indq = Ind[sqside].iterator_from_index(indq);
      prime_info_clear(pdata);
  } else {
      // TODO: we might want to implement this, one day.
      ASSERT_ALWAYS(q0 > qfac_max);
  }

  /****** Composite special-q ******/
  vector<index_t>::iterator first_indq2, last_indq2;
  vector<index_t>::iterator first_indq3, last_indq3;
  vector<index_t> comp2, comp3;
  
  if (compsq) {
      comp2 = all_comp_sq_2(q0, q1, qfac_min, qfac_max, Ind[sqside]);
      first_indq2 = comp2.begin();
      last_indq2 = comp2.end();
      comp3 = all_comp_sq_3(q0, q1, qfac_min, qfac_max, Ind[sqside]);
      first_indq3 = comp3.begin();
      last_indq3 = comp3.end();
  }

  /****** go multi-thread ******/
  uint64_t block = (last_indq - first_indq) / mt;
  uint64_t block2 = (last_indq2 - first_indq2) / (2*mt);
  uint64_t block3 = (last_indq3 - first_indq3) / (3*mt);

  double t0 = seconds ();
  double wct_t0 = wct_seconds ();

  pthread_t * thid = (pthread_t *)malloc(mt*sizeof(pthread_t));
  struct th_args * args = (struct th_args *)malloc(mt*sizeof(struct th_args));
  for (int i = 0; i < mt; ++i) {
      args[i].nq = 0;
      args[i].nq2 = 0;
      args[i].nq3 = 0;
      if (compsq) {
          args[i].list_q_comp2 = first_indq2 + 2*block2*i;
          args[i].list_q_comp3 = first_indq3 + 3*block3*i;
      } else {
          args[i].list_q_prime = first_indq + block*i;
      }
      if (i < mt-1) {
          if (compsq) {
              args[i].nq2 = block2;
              args[i].nq3 = block3;
          } else {
              args[i].nq = block;
          }
      } else {
          if (compsq) {
              args[i].nq2 = (last_indq2 - args[i].list_q_comp2)/2;
              args[i].nq3 = (last_indq3 - args[i].list_q_comp3)/3;
          } else {
              args[i].nq = last_indq - args[i].list_q_prime;
          }
      }
      args[i].rels = &rels;
      args[i].nrels = &nrels;
      args[i].Ind = &Ind[0];
      args[i].dl = dl;
      args[i].shrink_factor = shrink_factor;
      gmp_randinit_default(args[i].rstate);
      gmp_randseed_ui(args[i].rstate, 171717+i);

      pthread_create(&thid[i], NULL, do_thread, (void *)(&args[i]));
  }
  for (int i = 0; i < mt; ++i) {
      pthread_join(thid[i], NULL);
      gmp_randclear(args[i].rstate);
  }

  /* print statistics */
  t0 = seconds () - t0;
  printf ("# Output %lu relations in %.2fs cpu (%.0f rels/s)\n",
	  rels_printed, t0, (double) rels_printed / t0);
  wct_t0 = wct_seconds () - wct_t0;
  printf ("# Output %lu relations in %.2fs wct (%.0f rels/s)\n",
	  rels_printed, wct_t0, (double) rels_printed / wct_t0);

  gmp_randclear(global_rstate_non_mt);
  free(thid);
  free(args);
  renumber_clear(ren_table);
  cado_poly_clear(cpoly);
  param_list_clear(pl);

  return 0;
}
