#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <cstring> /* for strcmp() */
#include <cmath> /* for sqrt and floor and log and ceil */
#include <cstdint>

#include <algorithm>
#include <iosfwd>
#include <iostream>
#include <iterator>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <pthread.h>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "cado_poly.h"
#include "cxx_mpz.hpp"
#include "gzip.h"
#include "indexed_relation.hpp"
#include "special-q.hpp"
#include "macros.h"
#include "misc.h"
#include "params.h"
#include "relation.hpp"
#include "renumber.hpp"
#include "timing.h"
#include "typedefs.h"
#include "verbose.h"
#include "utils_cxx.hpp"

static int verbose = 0; /* verbosity level */

/*
 * The goal of this binary is to produce relations that try to be good
 * approximations of what las would produce, with respect to filtering.
 * The idea is to feed the purge/merge steps of the filtering with these
 * relations to get an idea of the final matrix that will come out of the
 * sieving step for a given set of parameters.
 *
 * The shell script cado-nfs/scripts/estimate_matsize/estimate_matsize.sh
 * is an attempt to run all the required steps for such a simulation:
 * sampling with las, fake relation generation with this binary, and
 * filter.
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

// A global mutex for I/O
static std::mutex io_mutex;

// A global variable for the number of relations printed
static unsigned long rels_printed = 0;

// Structure that contains indices (in the sense of renumber.[ch]) for
// one side, with the possibility to get a randomized ideal on the same
// side, and more or less of the same size.
// Usage:
//  - init()
//  - do many append() (increasing order of indices) and append_prime()
//  - finalize()
//  - can call random_index()
//  - can call iterator_from_index()    (XXX unused code !)
//  - can call pos_from_p() and p_from_pos() if append_prime() has been done

struct indexrange {
    std::vector <index_t>ind;
    std::vector <p_r_values_t>prime;

    void init() { }
    void finalize() { }

    void append(index_t z) {
        ind.push_back(z);
    }
    void append_prime(p_r_values_t p) {
        prime.push_back(p);
    }

    index_t random_index(index_t z, gmp_randstate_t buf) const {
        // Find position of z
        // Exact version:
        //  vector<index_t>::iterator it;
        //  it = lower_bound (ind.begin(), ind.end(), z);
        //  uint64_t position = it - ind.begin();
        // Fast, approximate version:
        uint64_t const position = z >> 1; // half of indices on each side
        auto const range = uint64_t((double(position))*0.2);
        uint64_t const low = position > range ? position - range : 0;
        uint64_t high = MIN(position + range, ind.size());
        if (high == low)
            high++;
        return ind[low + (u64_random(buf)%(high-low))];
    }

    p_r_values_t p_from_pos(uint64_t position) const {
        return prime[position];
    }

    size_t pos_from_p(p_r_values_t p) const {
        return std::ranges::lower_bound(prime, p) - prime.begin();
    }

    std::vector<index_t> all_composites(
            uint64_t q0,
            uint64_t q1,
            uint64_t qfac_min,
            uint64_t qfac_max,
            int nfactors) const;
    std::vector<std::vector<index_t>> all_composites(
            uint64_t q0,
            uint64_t q1,
            uint64_t qfac_min,
            uint64_t qfac_max) const;
};

// Fill in the indexrange data structure from the renumber table,
// gathering indices in two arrays, one for each side.
// In case of composite special-qs, also fill-in the list of the
// corresponding primes on the sqside.
static std::vector<indexrange> prepare_indexrange(renumber_t const & ren_tab, 
        int sqside, int compsq) {
    std::vector<indexrange> Ind(ren_tab.get_nb_polys());
    for(int side = 0 ; side < ren_tab.get_nb_polys() ; side++)
        Ind[side].init();
    index_t i = 0;
    for (auto it = ren_tab.begin() ; it != ren_tab.end() ; ++it, ++i) {
        if (ren_tab.is_additional_column(i))
            continue;
        renumber_t::p_r_side const x = *it;
        Ind[x.side].append(i);
        if (compsq && (x.side == sqside)) {
            Ind[sqside].append_prime(x.p);
        }
    }
    for(int side = 0 ; side < ren_tab.get_nb_polys() ; side++)
        Ind[side].finalize();

    return Ind;
}

static void remove_special_q(relation & rel, special_q const & Q)
{
    auto & V = rel.sides[Q.side];
    auto nn = V.begin();
    for(auto const & pr : V) {
        if (!Q.is_coprime_to(mpz_get_ui(pr.p)))
            continue;
        *nn++ = pr;
    }
    V.erase(nn, V.end());
}


struct model_relation : public indexed_relation_byside {
    template<typename... Args>
        model_relation(Args&&... args)
        : indexed_relation_byside { std::forward<Args>(args)... }
    {}
    model_relation perturb(std::vector<indexrange> const & Ind, gmp_randstate_t buf) const
    {
        relation_ab const ab(i64_random(buf), u64_random(buf));
        model_relation rel(ab);

        rel.set_active_sides(get_active_sides());
        for(size_t side = 0 ; side < sides.size() ; side++) {
            for(auto i : sides[side]) {
                rel[side].push_back(Ind[side].random_index(i, buf));
            }
        }
        return rel;
    }
};

// static std::map<special_q, std::vector<model_relation> >
static std::pair<std::vector<size_t>, std::vector<model_relation>>
read_sample_file(int sqside, const char *filename, renumber_t & ren_tab)
{
    ifstream_maybe_compressed in(filename);

    ASSERT_ALWAYS (in);

    std::set<special_q> current;
    std::map<special_q, std::vector<model_relation> > sample;

    int nbegin = 0, nend = 0, maxdepth = 0;

    for(std::string line ; std::getline(in, line) ; ) {
        if (line.rfind("# Now sieving side-", 0) != std::string::npos) {
            std::istringstream is(line.c_str() + 14);
            special_q Q;
            is >> Q;
            if (!is)
                throw std::runtime_error(fmt::format("parse error at line: {}", line));
            ASSERT_ALWAYS(sqside == Q.side);
            sample[Q];  // auto-vivify
            if (current.insert(Q).second) {
                /* When we start over, there's a "second" beginning. */
                nbegin++;
                maxdepth = std::max(maxdepth, nbegin - nend);
            }
        } else if (line.rfind("# Time for side-", 0) != std::string::npos) {
            nend++;
            for(auto & x : line) if (x == ':') x = ' ';
            std::istringstream is(line.c_str() + 11);
            special_q Q;
            is >> Q;
            if (!is)
                throw std::runtime_error(fmt::format("parse error at line: {}", line));
            current.erase(Q);
        } else if (line[0] == '#') {
            continue;
        } else {
            relation rel;
            std::istringstream(line) >> rel;
            for(auto const & Q : current) {
                if (current.size() > 1) {
                    /* If current.size() > 1 it's more difficult, as have to
                     * look up which Q is the good one. There may even be
                     * more than one, but we consider this unlikely here,
                     * since we don't expect more than a few "active"
                     * special-Qs to choose from
                     */

                    /* So do we have a-br = 0 mod p, with the usual
                     * conventions ? Those are actually quite annoying
                     * conventions in general. r represents a point in
                     * P1(Z/p), and we want to check if a:b matches r.
                     *
                     * In practice,
                     *  - Q is always square-free
                     *  - unless allow_compsq is set, Q is prime
                     *  - projective roots above special-q's are not
                     *  used (and this is true even for "partially
                     *  projective" roots over divisors of Q).
                     *
                     * So that it's actually easy, we just have to check
                     * that a-br = 0 mod p
                     */
                    cxx_mpz z;
                    mpz_mul(z, rel.bz, Q.r);
                    mpz_sub(z, rel.az, z);
                    mpz_mod(z, z, Q.p);
                    if (mpz_sgn(z) != 0)
                        continue;
                }
                remove_special_q(rel, Q);
                sample[Q].emplace_back(rel, ren_tab);
                break;
            }
        }
    }

    for(auto & S : sample)
        std::ranges::sort(S.second);

    if (nbegin == 0) {
        fmt::print(stderr, "# The sample file {} was apparently"
                " created without -v, but -v is mandatory"
                " for fake_rels\n", filename);
        exit(EXIT_FAILURE);
    }

    size_t nq = 0;
    size_t nr = 0;

    // return sample;
    std::pair<std::vector<size_t>, std::vector<model_relation>> ret;
    for(auto const & x : sample) {
        nq++;
        nr+=x.second.size();
        ret.first.push_back(x.second.size());
        std::ranges::copy(x.second, std::back_inserter(ret.second));
    }

    fmt::print("# {}: {} special-q's, {} relations, max {} concurrent special-q's\b",
            filename, nq, nr, maxdepth);
    return ret;
}

static unsigned long print_fake_rel_manyq(
        std::ostream& os,
        std::vector<index_t>::const_iterator qbegin,
        std::vector<index_t>::const_iterator qend,
        int nq,
        // std::vector<std::pair<special_q, std::vector<model_relation> > > const & sample,
        std::pair<std::vector<size_t>, std::vector<model_relation>> const & sample,
        std::vector<indexrange> const & Ind,
        int dl, double shrink_factor,
        gmp_randstate_t buf)
{
    /* There's a question of whether we print at each special-q, or only
     * once at the end. The former has the advantage of avoiding user
     * boredom.
     */
    unsigned long nrels_thread = 0;

    auto R = [buf](unsigned long n) { return u64_random(buf) % n; };

    ASSERT_ALWAYS((qend - qbegin) % nq == 0);

    for (auto it = qbegin ; it != qend ; ) {
        std::ostringstream oss;
        /* This is the part that will go in _all_ relations */
        std::vector<index_t> qpart;
        for(int n = nq ; n-- ; )
            qpart.push_back(*it++);

        auto const & model_nrels = sample.first[R(sample.first.size())];
        // special_q const & model_q = model.first;
        // auto const & model_nrels = model.second.size();

        int nr;
        if (shrink_factor == 1) {
            nr = model_nrels;
        } else {
            double const nr_dble = double_ratio(model_nrels, shrink_factor);
            // Do probabilistic rounding, in case nr_dble is small (maybe < 1)
            double const trunc_part = trunc(nr_dble);
            double const frac_part = nr_dble - trunc_part;
            double const rnd = double_ratio(u64_random(buf), UINT64_MAX);
            nr = int(trunc_part) + int(rnd < frac_part);
        }

        for( ; nr-- ; ) {
            // auto const & model_rel = model.second[R(model_nrels)];
            auto const & model_rel = sample.second[R(sample.second.size())];

            model_relation rel = model_rel.perturb(Ind, buf);

            /* Note that we always add the indices that correspond to q
             * to the **END** of the relation, irrespective of which side
             * q is on ! This is because:
             *   - it doesn't make any difference * down the line,
             *   - and keeping track of the proper side for q (or,
             *   conceivably, for the sides of all divisors of q !) would
             *   be a bit annoying here.
             */
            std::ranges::copy(qpart, std::back_inserter(rel.sides.back()));

            if (shrink_factor > 1)
                rel.shrink(shrink_factor);
            rel.sort();
            rel.compress(dl);

            oss << rel << "\n";
            nrels_thread++;
            // rels_printed++;
        }
        std::lock_guard<std::mutex> const dummy(io_mutex);
        os << oss.str();
    }
    return nrels_thread;
}

std::vector<index_t> indexrange::all_composites(uint64_t q0, uint64_t q1,
        uint64_t qfac_min,
        uint64_t qfac_max,
        int n) const
{
    std::vector<index_t> ret;

    if (n == 0)
        return ret;

    auto jt = std::back_inserter(ret);

    /* primes divisors of n-factor composite special-qs can never be
     * smaller than this:
     */
    uint64_t l1min = MAX(q0/pow(qfac_max, n-1), qfac_min);
    /* XXX prime at this position might be below l1min ! */
    uint64_t pos_l1min = pos_from_p(l1min);
    for( ; p_from_pos(pos_l1min) < l1min ; pos_l1min++);
    l1min = p_from_pos(pos_l1min);

    /* and never bigger than this: */
    uint64_t const l1max = MIN(qfac_max, round(pow(q1, 1/(double) n)));
    uint64_t const pos_l1max = pos_from_p(l1max);

    for (uint64_t pos1 = pos_l1min; pos1 < pos_l1max; ++pos1) {
        /* look for cases where _this_ prime is the smallest one */
        uint64_t const l1 = p_from_pos(pos1);
        /* q0 <= l1 * x < q1
         * implies q0/l1 <= x < q1/l1
         * the left part is easy, but for the right part, we rewrite as:
         * l1 * x <= q1-1
         * x <= (q1-1)/l1 ---> x <= floor((q1-1)/l1)  --> x < ceil(q1/l1)
         */
        if (n == 1) {
            /* no point in recursing to determine a list of zero-length
             * continuations.
             */
            *jt++ = pos1;
        } else {
            auto tail = all_composites(iceildiv(q0,l1), iceildiv(q1,l1), l1, qfac_max, n-1);
            for(auto it = tail.begin() ; it != tail.end() ; ) {
                *jt++ = pos1;
                for(int j=n-1 ; j-- ; )
                    *jt++ = *it++;
            }
        }
    }

    return ret;
}

std::vector<std::vector<index_t>> indexrange::all_composites(uint64_t q0, uint64_t q1,
        uint64_t qfac_min,
        uint64_t qfac_max) const
{
    std::vector<std::vector<index_t>> list;
    list.emplace_back();
    list.emplace_back();
    for(int n = 2 ; ; n++) {
        list.emplace_back(all_composites(q0, q1, qfac_min, qfac_max, n));
        if (list.back().empty()) {
            list.pop_back();
            break;
        }
        fmt::print("# Got {} {}-composite sq\n",
                (size_t)(list.back().size()/n), n);
    }
    return list;
}

static void worker(int tnum, int nt,
        std::vector<indexrange> const & Ind,
        // std::vector<std::pair<special_q, std::vector<model_relation>>> const & sample,
        std::pair<std::vector<size_t>, std::vector<model_relation>> const & sample,
        std::vector<std::vector<index_t>> const & qs,
        double shrink_factor, int dl, unsigned long seed)
{
    gmp_randstate_t buf;
    gmp_randinit_default(buf);
    gmp_randseed_ui(buf, seed + tnum);
    unsigned long ret = 0;
    for(size_t n = 1 ; n < qs.size() ; n++) {
        ASSERT_ALWAYS(qs[n].size() % n == 0);
        auto it0 = qs[n].begin() + n * ((tnum * qs[n].size() / n) / nt);
        auto it1 = qs[n].begin() + n * (((tnum + 1) * qs[n].size() / n) / nt);
        ret += print_fake_rel_manyq(
                std::cout,
                it0, it1,
                n,
                sample,
                Ind,
                dl, shrink_factor,
                buf);
    }
    gmp_randclear(buf);
    std::lock_guard<std::mutex> const dummy(io_mutex);
    rels_printed += ret;
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
    param_list_decl_usage(pl, "seed", "random seed");
    param_list_decl_usage(pl, "t", "number of threads to use");
    param_list_decl_usage(pl, "v", "verbose mode");
    verbose_decl_usage(pl);
}

// coverity[root_function]
int main(int argc, char const * argv[])
{
  cxx_param_list pl;
  cxx_cado_poly cpoly;
  int sqside = -1;
  const char *argv0 = argv[0];
  int lpb[2] = {0, 0};
  uint64_t q0 = 0;
  uint64_t q1 = 0;
  int dl = 0;
  int mt = 1;
  int compsq = 0;
  uint64_t qfac_min = 1024;
  uint64_t qfac_max = UINT64_MAX;
  double shrink_factor = 1; // by default, no shrink

  declare_usage(pl);
  param_list_configure_switch (pl, "-v", &verbose);
  param_list_configure_switch(pl, "-dl", &dl);
  param_list_configure_switch(pl, "-allow-compsq", &compsq);

  argv++, argc--;
  for( ; argc ; ) {
      if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }

      /* Could also be a file */
      FILE *f;
      if ((f = fopen(argv[0], "r")) != nullptr) {
          param_list_read_stream(pl, f, 0);
          fclose(f);
          argv++,argc--;
          continue;
      }

      fmt::print(stderr, "Unhandled parameter {}\n", argv[0]);
      param_list_print_usage(pl, argv0, stderr);
      exit (EXIT_FAILURE);
  }
  verbose_interpret_parameters(pl);
  param_list_print_command_line(stdout, pl);

  const char * filename;
  if ((filename = param_list_lookup_string(pl, "poly")) == nullptr) {
      fmt::print(stderr, "Error: parameter -poly is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_int(pl, "lpb0", &lpb[0]);
  if (lpb[0] == 0) {
      fmt::print(stderr, "Error: parameter -lpb0 is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }
  param_list_parse_int(pl, "lpb1", &lpb[1]);
  if (lpb[1] == 0) {
      fmt::print(stderr, "Error: parameter -lpb1 is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_uint64(pl, "q0", &q0);
  if (q0 == 0) {
      fmt::print(stderr, "Error: parameter -q0 is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_uint64(pl, "q1", &q1);
  if (q1 == 0) {
      fmt::print(stderr, "Error: parameter -q1 is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }
  
  param_list_parse_double(pl, "shrink-factor", &shrink_factor);
  if (shrink_factor < 1) {
      fmt::print(stderr, "Error: shrink factor must be an integer >= 1\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  param_list_parse_int(pl, "t", &mt);
  param_list_parse_uint64(pl, "qfac-min", &qfac_min);
  param_list_parse_uint64(pl, "qfac-max", &qfac_max);

  unsigned long seed = 171717;
  param_list_parse_ulong(pl, "seed", &seed);

  if (!cado_poly_read(cpoly, filename))
    {
      fmt::print (stderr, "Error reading polynomial file {}\n", filename);
      exit (EXIT_FAILURE);
    }

  param_list_parse_int(pl, "sqside", &sqside);
  if (sqside == -1 || sqside > 2) {
      fmt::print(stderr, "Error: sqside must be 0 or 1\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  const char * renumberfile;
  if ((renumberfile = param_list_lookup_string(pl, "renumber")) == nullptr) {
      fmt::print(stderr, "Error: parameter -renumber is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }
  fmt::print ("# Start reading renumber table\n");
  fflush (stdout);
  renumber_t ren_table(cpoly);
  ren_table.read_from_file(renumberfile, dl);
  fmt::print ("# Done reading renumber table\n");
  fflush (stdout);

  for (int side = 0; side < ren_table.get_nb_polys(); ++side) {
      if (ren_table.get_lpb(side) != (unsigned long)lpb[side]) {
          fmt::print(stderr, "Error: on side {}, lpb on the command-line is different from the one in the renumber file\n", side);
          exit(EXIT_FAILURE);
      }
  }

  // read sample file
  const char * samplefile;
  if ((samplefile = param_list_lookup_string(pl, "sample")) == nullptr) {
      fmt::print(stderr, "Error: parameter -sample is mandatory\n");
      param_list_print_usage(pl, argv0, stderr);
      exit(EXIT_FAILURE);
  }

  fmt::print ("# Start reading sample file\n");
  fflush (stdout);

  std::pair<std::vector<size_t>, std::vector<model_relation>> const sample = read_sample_file(sqside, samplefile, ren_table);

  /*
  std::vector<std::pair<special_q, std::vector<model_relation>>>
      sample;
  for(auto& x : read_sample_file(sqside, samplefile, ren_table))
      sample.emplace_back(std::move(x));
  if (verbose) {
      for(auto const & x : sample) {
          std::cout << "# " << x.first << ": " << x.second.size() << " relations\n";
      }
  }
  */

  fmt::print ("# Done reading sample file\n");
  fflush (stdout);

  param_list_warn_unused(pl);

  // Two index ranges, one for each side
  fmt::print ("# Start preparing index ranges\n");
  fflush (stdout);
  std::vector<indexrange> Ind = prepare_indexrange(ren_table, sqside, compsq);
  fmt::print ("# Done preparing index ranges\n");
  fflush (stdout);

  std::vector<std::vector<index_t>> qs;

  /* the indexranges contain the prime-to-index conversion functionality
   * only if compsq holds. In the non-composite case, we fill qs[1]
   * differently, because we don't want to 
   */
  if (!compsq) {
      qs.emplace_back(); /* 0 -- placeholder */
      qs.emplace_back(); /* 1 */

      /* jump-start in the renumber table, and enumerate the large primes
       * from there.
       */
      index_t i = ren_table.index_from_p(q0, sqside);
      renumber_t::const_iterator it(ren_table, i);
      for( ; it != ren_table.end() && (*it).p < q1 ; ++it, ++i) {
          if ((*it).side == sqside)
              qs[1].push_back(i);
      }
  } else {
      /* In theory, this is overkill. We don't _really_ need the prime[]
       * cache in the indexrange type. Except that the renumber
       * "traditional" format has expensive *_from_index lookups, and
       * that would get in the way if we want to do away with
       * indexrange::prime[]. Also, we don't have bidirectional iterators
       * on the renumber table, so it's not really practical to use only
       * the renumber table.
       */
      qs = Ind[sqside].all_composites(q0, q1, qfac_min, qfac_max);
  }


  std::vector<std::thread> threads;

  double t0 = seconds ();
  double wct_t0 = wct_seconds ();

  threads.reserve(mt);
  for(int i = 0 ; i < mt ; i++)
      threads.emplace_back(worker, i, mt, 
              Ind, sample, qs,
              shrink_factor, dl, seed);
  for(auto & t : threads)
      t.join();


  t0 = seconds () - t0;
  wct_t0 = wct_seconds () - wct_t0;

  /* print statistics */

  fmt::print ("# Output {} relations in {:.2f}s cpu ({:.0f} rels/s)\n",
	  rels_printed, t0, double_ratio(rels_printed, t0));
  fmt::print ("# Output {} relations in {:.2f}s wct ({:.0f} rels/s)\n",
	  rels_printed, wct_t0, double_ratio(rels_printed, wct_t0));

  return 0;
}
