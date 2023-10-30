#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <ext/alloc_traits.h>
// IWYU pragma: no_include <memory>
#include <cerrno>              // for ENOENT, errno
#include <climits>             // for UINT_MAX
#include <cstdlib>             // for exit, EXIT_FAILURE
#include <cstring>             // for strncmp, strlen, strrchr
#include <cstdio>               // for size_t, sscanf, fclose, fopen, fread

#include <algorithm>            // for sort, min
#include <map>                  // for map<>::mapped_type, _Rb_tree_iterator
#include <memory>
#include <sstream>              // for ostringstream, operator<<, basic_ostream
#include <stdexcept>            // for runtime_error
#include <string>               // for string, basic_string, operator<<, cha...
#include <utility>              // for pair, make_pair
#include <vector>               // for vector

#include <sys/stat.h>           // for stat
#include <gmp.h>                // for mpz_cmp_ui
#include "bw-common.h"          // for bw, bw_common_clear, bw_common_decl_u...

#include "fmt/core.h"           // for check_format_string
#include "fmt/format.h"         // for basic_buffer::append, basic_parse_con...
#include "macros.h"             // for ASSERT_ALWAYS, MAYBE_UNUSED
#include "misc.h"               // ok_NOKNOK
#include "arith-generic.hpp"
#include "arith-cross.hpp"
#include "params.h"             // for param_list_clear, param_list_init
#include "portability.h" // asprintf // IWYU pragma: keep
#include "select_mpi.h"         // for MPI_Abort, MPI_Comm_rank, MPI_COMM_WORLD

using namespace std;

// arguments we need
// m
// n
// prime
//
// wdir only as a way to shorten file names.
// nullspace ??

const char * my_basename(const char * x)
{
    const char * p = strrchr(x, '/');
    if (p) {
        p++;
    } else {
        p = x;
    }
    return p;
}

struct Cfile : public string {
    unsigned int j0, j1;
    unsigned int stretch;
    Cfile(const char * x) : string(x) {
        int rc = sscanf(my_basename(x), "Cv%u-%u.%u", &j0, &j1, &stretch);
        if (rc != 3) throw std::runtime_error("want Cv%u-%u.%u");
    }
    inline bool operator<(Cfile const& o) {
        return stretch < o.stretch;
    }
};

struct Dfile : public string {
    unsigned int j0, j1;
    unsigned int stretch;
    Dfile(const char * x) : string(x) {
        int rc = sscanf(my_basename(x), "Cd%u-%u.%u", &j0, &j1, &stretch);
        if (rc != 3) throw std::runtime_error("want Cd%u-%u.%u");
    }
    inline bool operator<(Dfile const& o) {
        return stretch < o.stretch;
    }
};

struct Rfile : public string {
    unsigned int nchecks;
    Rfile(const char * x) : string(x) {
        unsigned int nc0, nc1;
        int rc = sscanf(my_basename(x), "Cr0-%u.0-%u", &nc0, &nc1);
        if (rc != 2) throw std::runtime_error("want Cr0-%u.0-%u");
        nchecks = nc0;
        if (nc0 != nc1) throw std::runtime_error("want Cr0-NCHECKS.0-NCHECKS");
    }
    inline bool operator<(Rfile const& o) {
        return nchecks < o.nchecks;
    }
};

struct Tfile : public string {
    unsigned int nchecks;
    int m;
    Tfile(const char * x) : string(x) {
        int rc = sscanf(my_basename(x), "Ct0-%u.0-%d", &nchecks, &m);
        if (rc != 2) throw std::runtime_error("want Ct0-%u.0-%d");
    }
    inline bool operator<(Tfile const& o) {
        if (nchecks != o.nchecks) return nchecks < o.nchecks;
        return m < o.m;
    }
};

struct Vfile : public string {
    unsigned int j0, j1, n;
    int checks = 0;
    std::pair<unsigned int, unsigned int> seq_id() const {
        return make_pair(j0, j1);
    }
    Vfile(const char * x) : string(x) {
        int rc = sscanf(my_basename(x), "V%u-%u.%u", &j0, &j1, &n);
        if (rc != 3) throw std::runtime_error("want " "V%u-%u.%u");
    }
    inline bool operator<(Vfile const& o) {
        if (j0 != o.j0) return j0 < o.j0;
        if (j1 != o.j1) return j1 < o.j1;
        if (n != o.n) return n < o.n;
        return false;
    }
};

struct Afile : public string {
    unsigned int j0, j1, n0, n1;
    int checks = 0;
    std::pair<unsigned int, unsigned int> seq_id() const {
        return make_pair(j0, j1);
    }
    std::pair<unsigned int, unsigned int> iter_range() const {
        return make_pair(n0, n1);
    }
    Afile(const char * x) : string(x) {
        int rc = sscanf(my_basename(x), "A%u-%u.%u-%u", &j0, &j1, &n0, &n1);
        if (rc != 4) throw std::runtime_error("want " "A%u-%u.%u-%u");
    }
    inline bool operator<(Afile const& o) {
        if (j0 != o.j0) return j0 < o.j0;
        if (j1 != o.j1) return j1 < o.j1;
        if (n0 != o.n0) return n0 < o.n0;
        if (n1 != o.n1) return n1 < o.n1;
        return false;
    }
};

#if 0
struct Ffile : public string {
    unsigned int s0, s1, j0, j1;
    int checks = 0;
    std::pair<unsigned int, unsigned int> seq_id() const {
        return make_pair(j0, j1);
    }
    std::pair<unsigned int, unsigned int> sol_id() const {
        return make_pair(s0, s1);
    }
    Ffile(const char * x) : string(x) {
        int rc = sscanf(my_basename(x), F_FILE_SLICE_PATTERN2, &s0, &s1, &j0, &j1);
        if (rc != 4) throw std::runtime_error("want " F_FILE_SLICE_PATTERN2);
    }
    inline bool operator<(Ffile const& o) {
        if (s0 != o.s0) return s0 < o.s0;
        if (s1 != o.s1) return s1 < o.s1;
        if (j0 != o.j0) return j0 < o.j0;
        if (j1 != o.j1) return j1 < o.j1;
        return false;
    }
};

struct Sfile : public string {
    int s0, s1, j0, j1, n0, n1;
    int checks = 0;
    std::pair<unsigned int, unsigned int> seq_id() const {
        return make_pair(j0, j1);
    }
    std::pair<unsigned int, unsigned int> sol_id() const {
        return make_pair(s0, s1);
    }
    std::pair<unsigned int, unsigned int> iter_range() const {
        return make_pair(n0, n1);
    }
    Sfile(const char * x) : string(x) {
        int rc = sscanf(my_basename(x), S_FILE_BASE_PATTERN ".%u-%u", &s0, &s1, &j0, &j1, &n0, &n1);
        if (rc != 6) throw std::runtime_error("want " S_FILE_BASE_PATTERN ".%u-%u");
    }
    inline bool operator<(Sfile const& o) {
        if (s0 != o.s0) return s0 < o.s0;
        if (s1 != o.s1) return s1 < o.s1;
        if (j0 != o.j0) return j0 < o.j0;
        if (j1 != o.j1) return j1 < o.j1;
        if (n0 != o.n0) return n0 < o.n0;
        if (n1 != o.n1) return n1 < o.n1;
        return false;
    }
};
#endif

void vec_alloc(arith_generic * A, arith_generic::elt *& z, size_t vsize)
{
    z = A->alloc(vsize, ALIGNMENT_ON_ALL_BWC_VECTORS);
    A->vec_set_zero(z, vsize);
}

void vec_free(arith_generic * A, arith_generic::elt *& z, size_t vsize MAYBE_UNUSED)
{
    A->free(z);
}

int vec_read(arith_generic * A, void * z, string const & v, size_t vsize, const char * prefix = NULL)
{
    fmt::print(FMT_STRING("{} {} ..."), prefix, v);
    FILE * f;
    if ((f = fopen(v.c_str(), "rb")) != NULL) {
        int rc = fread(z, A->elt_stride(), vsize, f);
        fclose(f);
        if (rc >= 0 && (size_t) rc == vsize) {
            fmt::print(FMT_STRING("{}"), " done\n");
            return rc;
        }
    }
    fmt::print(FMT_STRING("{}"), " failed\n");
    return -1;
}

size_t vec_items(arith_generic * A, string const & v)
{
    struct stat sbuf[1];
    int rc = stat(v.c_str(), sbuf);
    if (rc < 0 && errno == ENOENT)
        return 0;
    ASSERT_ALWAYS(rc == 0);
    return sbuf->st_size / A->elt_stride();
}

template<typename T>
size_t common_size(arith_generic * Ac, std::vector<T> const & Cfiles, const char * name)
{
    size_t vsize = 0;
    std::string vsize_first;

    for(auto & C : Cfiles) {
        size_t items = vec_items(Ac, C);
        if (items == 0) {
            fmt::print(FMT_STRING("{} has disappeared\n"), C);
            continue;
        }
        if (vsize == 0) {
            vsize = items;
            vsize_first = C;
        } else if (vsize != items) {
            fmt::print(stderr,
                    FMT_STRING("File sizes disagree for {} ({} items) and {} ({} items)\n"),
                    vsize_first, vsize, C, items);
            exit(EXIT_FAILURE);
        }
    }
    if (vsize) fmt::print(FMT_STRING("{} files have {} coordinates\n"), name, vsize);
    return vsize;
}

typedef std::map<pair<unsigned int, unsigned int>, vector<Vfile> > vseq_t;
void check_V_files(arith_generic * Ac, vseq_t & Vsequences, std::vector<Cfile> & Cfiles, int & nfailed)/*{{{*/
{
    if (Cfiles.empty()) return;

    int nchecks = Ac->simd_groupsize();
    size_t vsize = common_size(Ac, Cfiles, "Cv");

    for(unsigned int i0 = 0 ; i0 + 1 < Cfiles.size() ; i0++) {
        Cfile& C_i0(Cfiles[i0]);
        arith_generic::elt * Cv_i0;
        vec_alloc(Ac, Cv_i0, vsize);

        int has_read_Cv_i0 = 0;

        for(unsigned int i1 = i0 + 1 ; i1 < Cfiles.size() ; i1++) {
            Cfile& C_i1(Cfiles[i1]);
            const char * c = C_i1.c_str();

            fmt::print(FMT_STRING("Doing checks for distance {} using {} and {}\n"),
                    C_i1.stretch - C_i0.stretch,
                    C_i0, C_i1
                    );

            arith_generic::elt * Cv_i1;
            vec_alloc(Ac, Cv_i1, vsize);

            int has_read_Cv_i1 = 0;

            /* {{{ check all V files together */
            for(vseq_t::iterator it = Vsequences.begin(); it != Vsequences.end(); it++)
            {
                vector<Vfile> & Vs(it->second);

                fmt::print(FMT_STRING(" checks on V files for sequence {}-{}\n"),
                        it->first.first, it->first.second);

                std::unique_ptr<arith_generic> Av(arith_generic::instance(bw->p,it->first.second - it->first.first));

                std::unique_ptr<arith_cross_generic> AcxAv(arith_cross_generic::instance(Ac, Av.get()));


                /* {{{ Check that all V files here have the proper size */
                for(unsigned int i = 0 ; i < Vs.size() ; i++) {
                    size_t items = vec_items(Av.get(), Vs[i]);
                    if (items != vsize) {
                        fmt::print(stderr, FMT_STRING("{} has {} coordinates, different from expected {}\n"), Vs[i], items, vsize);
                        exit(EXIT_FAILURE);
                    }
                }
                /* }}} */

                arith_generic::elt * Vv;
                vec_alloc(Av.get(), Vv, vsize);
                unsigned int Vv_iter = UINT_MAX;

                arith_generic::elt * dotprod_scratch[2];
                vec_alloc(Av.get(), dotprod_scratch[0], nchecks);
                vec_alloc(Av.get(), dotprod_scratch[1], nchecks);

                unsigned int j = 0;
                for(unsigned int i = 0 ; i < Vs.size() ; i++) {
                    for(j = i + 1 ; j < Vs.size() ; j++) {
                        if (Vs[j].n + C_i0.stretch == Vs[i].n + C_i1.stretch)
                            break;
                    }
                    if (j == Vs.size()) continue;

                    if (!has_read_Cv_i0++)
                        if (vec_read(Ac, Cv_i0, C_i0, vsize, " ") < 0)
                            continue;

                    if (!has_read_Cv_i1++)
                        if (vec_read(Ac, Cv_i1, C_i1, vsize, " ") < 0)
                            continue;

                    Vs[j].checks++;
                    const char * vi = Vs[i].c_str();
                    const char * vj = Vs[j].c_str();
                    if (strncmp(vi, c, my_basename(c) - c) == 0) {
                        vi += my_basename(c) - c;
                    }
                    if (strncmp(vj, c, my_basename(c) - c) == 0) {
                        vj += my_basename(c) - c;
                    }
                    fmt::print(FMT_STRING("  check {} against {}\n"), vi, vj);
                    if (Vs[i].n != Vv_iter) {
                        if (vec_read(Ac, Vv, Vs[i].c_str(), vsize, "   ") < 0)
                            continue;
                        Vv_iter = Vs[i].n;
                    }

                    Av->vec_set_zero(dotprod_scratch[0], nchecks);

                    /* compute the dot product */
                    AcxAv->add_dotprod(
                            dotprod_scratch[0],
                            Cv_i1,
                            Vv,
                            vsize);

                    if (vec_read(Ac, Vv, Vs[j].c_str(), vsize, "   ") < 0)
                        continue;

                    Vv_iter = Vs[j].n;

                    Av->vec_set_zero(dotprod_scratch[1], nchecks);
                    AcxAv->add_dotprod(
                            dotprod_scratch[1],
                            Cv_i0,
                            Vv,
                            vsize);

                    int cmp = Av->vec_cmp(dotprod_scratch[0], dotprod_scratch[1], nchecks);

                    std::string diag = fmt::format(
                            FMT_STRING("  check {} against {} -> {}\n"),
                            vi, vj, ok_NOKNOK(cmp == 0));

                    fmt::print(FMT_STRING("{}"), diag);

                    if (cmp != 0) {
                        nfailed++;
                        fmt::print(stderr, FMT_STRING("{}"), diag);
                    }
                }
                vec_free(Av.get(), dotprod_scratch[0], nchecks);
                vec_free(Av.get(), dotprod_scratch[1], nchecks);
                vec_free(Ac, Vv, vsize);

            }
            /* }}} */

            vec_free(Ac, Cv_i1, vsize);
        }
        vec_free(Ac, Cv_i0, vsize);
    }
}/*}}}*/

void check_A_files(arith_generic * Ac, std::vector<Vfile> const & Vfiles, std::vector<Afile> const & Afiles, std::vector<Dfile> const & Dfiles, Rfile & R, Tfile & T, int & nfailed)
{
    if (Dfiles.empty())
        return;
    arith_generic::elt * Dv = NULL;
    size_t vsize = common_size(Ac, Dfiles, "Cd");
    int nchecks = Ac->simd_groupsize();

    size_t rsize = vec_items(Ac, R) / nchecks;
    fmt::print(FMT_STRING("Cr file has {} coordinates\n"), rsize);

    ASSERT_ALWAYS(vec_items(Ac, T) ==  (size_t) bw->m);

    vec_alloc(Ac, Dv, vsize);

    int rc;

    arith_generic::elt * Tdata = NULL;
    Tdata = Ac->alloc(bw->m, ALIGNMENT_ON_ALL_BWC_VECTORS);
    FILE * Tfile = fopen(T.c_str(), "rb");
    rc = fread(Tdata, Ac->vec_elt_stride(bw->m), 1, Tfile);
    ASSERT_ALWAYS(rc == 1);
    fclose(Tfile);

    arith_generic::elt * Rdata = NULL;
    Rdata = Ac->alloc(Ac->vec_elt_stride(nchecks) * rsize, ALIGNMENT_ON_ALL_BWC_VECTORS);
    FILE * Rfile = fopen(R.c_str(), "rb");
    rc = fread(Rdata, Ac->vec_elt_stride(nchecks), rsize, Rfile);
    ASSERT_ALWAYS(rc == (int) rsize);
    fclose(Rfile);

    for(auto & D : Dfiles) {
        if (D.stretch == 0)
            continue;
        if (D.stretch > rsize) {
            fmt::print(stderr, FMT_STRING("Cannot do checks using {}, too few items in R file\n"), R);
            continue;
        }
        fmt::print(FMT_STRING("Doing A file checks for distance {} using {}\n"
                "  (as well as {} and {})\n"),
                D.stretch, D, T, R);
        int has_read_D = 0;
        /* first scan potential base files V, and the restrict to cases
         * where we have all the required A files to do a check...
         */
        for(auto & V0 : Vfiles) {
            /* Look for "A%u-%u.%u-%u" % (V0.j0, V0.j1, V0.n,
             * V0.n+D.stretch), or any collection of files that would
             * concatenate to exactly this file */
            std::ostringstream a_list;
            size_t n_reach = V0.n;
            /* A_files are sorted */
            for(auto const & A : Afiles) {
                if (A.j1 <= V0.j0) continue;
                if (A.j0 > V0.j0) break;
                if (A.n0 > n_reach) break;
                if (A.n1 <= n_reach) continue;
                if (A.n0 <= n_reach) {
                    ASSERT_ALWAYS(n_reach == V0.n || n_reach == A.n0);
                    /* Since we (still) collate A files, it's important
                     * to be able to do this check even if V starts in
                     * the middle of the A file -- this accounts for the
                     * n_reach == V0.n possibility. Otherwise, we must
                     * go from one A file to the next, hence the
                     * requirement that n_reach == A.n0 .*/
                    n_reach = A.n1;
                    a_list << " " << A;
                }
                if (n_reach >= V0.n + D.stretch)
                    break;
            }

            if (n_reach < V0.n + D.stretch)
                continue;

            fmt::print(FMT_STRING("  check {} against {} entries of{}\n"),
                    V0, D.stretch, a_list.str());

            std::unique_ptr<arith_generic> Av(arith_generic::instance(bw->p, V0.j1 - V0.j0));
            std::unique_ptr<arith_cross_generic> AcxAv(arith_cross_generic::instance(Ac, Av.get()));

            arith_generic::elt * dotprod_scratch[3];
            vec_alloc(Av.get(), dotprod_scratch[0], nchecks);
            vec_alloc(Av.get(), dotprod_scratch[1], nchecks);
            vec_alloc(Av.get(), dotprod_scratch[2], nchecks);

            /* read data from the A files. We redo the detection loop
             * that we had above */
            n_reach = V0.n;
            for(auto const & A : Afiles) {
                if (A.j1 <= V0.j0) continue;
                if (A.j0 > V0.j0) break;
                if (A.n0 > n_reach) break;
                if (A.n1 <= n_reach) continue;
                if (A.n0 <= n_reach) {
                    ASSERT_ALWAYS(n_reach == V0.n || n_reach == A.n0);
                    fmt::print(FMT_STRING("   read {} small {}*{} matrices from {}\n"),
                            std::min(A.n1, V0.n + D.stretch) - n_reach,
                            bw->m, bw->n, A);
                    FILE * a = fopen(A.c_str(), "rb");
                    for(unsigned int p = n_reach ; p < A.n1 && p < V0.n + D.stretch ; p++) {
                        Av->vec_set_zero(dotprod_scratch[1], nchecks);
                        for(int c = 0 ; c < bw->m ; c += nchecks) {
                            int rc;
                            for(int r = 0 ; r < nchecks ; r++) {
                                size_t simd = Av->simd_groupsize();
                                size_t rowsize = (A.j1 - A.j0) / simd;
                                size_t matsize = bw->m * rowsize;
                                size_t nmats = p - A.n0;

                                rc = fseek(a,
                                            nmats * Av->vec_elt_stride(matsize) +
                                            (c + r) * Av->vec_elt_stride(rowsize) +
                                            Av->vec_elt_stride((V0.j0 - A.j0) / simd),
                                        SEEK_SET);
                                ASSERT_ALWAYS(rc == 0);
                                rc = fread(Av->vec_subvec(dotprod_scratch[2], r), Av->elt_stride(), 1, a);
                                ASSERT_ALWAYS(rc == 1);
                            }
                            AcxAv->add_dotprod(
                                    dotprod_scratch[1],
                                    Ac->vec_subvec(Tdata, c),
                                    dotprod_scratch[2],
                                    nchecks);
                        }
                        AcxAv->add_dotprod(
                                dotprod_scratch[0],
                                Ac->vec_subvec(Rdata, (p - V0.n) * nchecks),
                                dotprod_scratch[1],
                                nchecks);
                    }
                    fclose(a);
                    n_reach = A.n1;
                }
                if (n_reach >= V0.n + D.stretch)
                    break;
            }

            int can_check = 1;

            if (!has_read_D++)
                if (vec_read(Ac, Dv, D.c_str(), vsize, "   ") < 0)
                    can_check = 0;

            if (can_check) {
                arith_generic::elt * Vv;

                vec_alloc(Av.get(), Vv, vsize);
                if (vec_read(Av.get(), Vv, V0.c_str(), vsize, "   ") < 0)
                    can_check = 0;

                if (can_check) {
                    Av->vec_set_zero(dotprod_scratch[1], nchecks);
                    AcxAv->add_dotprod(
                            dotprod_scratch[1],
                            Dv,
                            Vv,
                            vsize);
                }

                vec_free(Av.get(), Vv, vsize);
            }

            if (can_check) {
                int cmp = Av->vec_cmp(dotprod_scratch[0], dotprod_scratch[1], nchecks);

                std::string diag = fmt::format(
                        FMT_STRING("  check {} against {} entries of{} -> {}\n"),
                        V0, D.stretch, a_list.str(),
                        ok_NOKNOK(cmp == 0));
                fmt::print(FMT_STRING("{}"), diag);

                if (cmp != 0) {
                    nfailed++;
                    fmt::print(stderr, FMT_STRING("{}"), diag);
                }
            }

            if (!can_check)
                fmt::print(FMT_STRING("{}"), "  (check aborted because of missing files)\n");

            vec_free(Av.get(), dotprod_scratch[2], nchecks);
            vec_free(Av.get(), dotprod_scratch[1], nchecks);
            vec_free(Av.get(), dotprod_scratch[0], nchecks);
        }
    }
    Ac->free(Rdata);
    Ac->free(Tdata);
    vec_free(Ac, Dv, vsize);
}

/* This is *not* a parallel program, so we depart significantly from the
 * way programs such as krylov or mksol are written.
 *
 */
void * check_prog(param_list pl MAYBE_UNUSED, int argc, char * argv[])
{
    int withcoeffs = mpz_cmp_ui(bw->p, 2) > 0;
    int nchecks = withcoeffs ? NCHECKS_CHECK_VECTOR_GFp : NCHECKS_CHECK_VECTOR_GF2;
    std::unique_ptr<arith_generic> Ac(arith_generic::instance(bw->p, nchecks));

    vector<Cfile> Cfiles;
    vector<Dfile> Dfiles;
    vector<Rfile> Rfiles;
    vector<Tfile> Tfiles;
    vector<Vfile> Vfiles;
    vector<Afile> Afiles;
    // vector<Ffile> Ffiles;
    // vector<Sfile> Sfiles;

    for(int i = 0 ; i < argc ; i++) {
        try {
            const char * p = my_basename(argv[i]);
#ifdef __COVERITY__
            __coverity_mark_pointee_as_sanitized__(p, PATH);
#endif
            if (*p && p[strlen(p)-1] == '~')
                continue;
            switch(*p) {
                case 'C': 
                    switch(p[1]) {
                        case 'v': Cfiles.push_back(argv[i]); break;
                        case 'd': Dfiles.push_back(argv[i]); break;
                        case 'r': Rfiles.push_back(argv[i]); break;
                        case 't': Tfiles.push_back(argv[i]); break;
                    }
                    break;
                case 'A': Afiles.push_back(argv[i]); break;
                case 'V': Vfiles.push_back(argv[i]); break;
#if 0
                case 'F':
                          if (strncmp(p, "F.sols", 6) == 0)
                              Ffiles.push_back(argv[i]);
                          /* We purposefully skip the F files.
                          */
                          break;
                case 'S': Sfiles.push_back(argv[i]); break;
#endif
                          /* since we changed the format, and since we
                           * were not doing checks of any sort yet
                           * anyway, just discard them right away */
                case 'F': case 'S': break;
                default:
                        fmt::print(stderr, FMT_STRING("File name not recognized: {}\n"), argv[i]);
                        exit(EXIT_FAILURE);
            }
        } catch (std::runtime_error& e) {
            fmt::print(stderr, FMT_STRING("Parse error on {}: {}\n"), argv[i], e.what());
            exit(EXIT_FAILURE);
        }
    }

    /* {{{ some consistency checks */
    for(auto & C : Cfiles) {
        ASSERT_ALWAYS(C.j0 == 0);
        ASSERT_ALWAYS(C.j1 == (unsigned int) nchecks);
    }
    for(auto & D : Dfiles) {
        ASSERT_ALWAYS(D.j0 == 0);
        ASSERT_ALWAYS(D.j1 == (unsigned int) nchecks);
    }
    ASSERT_ALWAYS(Rfiles.size() <= 1);
    ASSERT_ALWAYS(Tfiles.size() <= 1);
    for(auto & R : Rfiles) {
        ASSERT_ALWAYS(R.nchecks == (unsigned int) nchecks);
    }
    for(auto & T : Tfiles) {
        ASSERT_ALWAYS(T.nchecks == (unsigned int) nchecks);
        /* T files for different m's could maybe coexist, at least in
         * theory. However since both C and D depend on T, this does not
         * seem very viable */
        ASSERT_ALWAYS(T.m == bw->m);
    }
    /* }}} */

    /* {{{ sort */
    std::sort(Cfiles.begin(), Cfiles.end());
    std::sort(Dfiles.begin(), Dfiles.end());
    std::sort(Rfiles.begin(), Rfiles.end());
    std::sort(Afiles.begin(), Afiles.end());
    std::sort(Vfiles.begin(), Vfiles.end());
    // std::sort(Ffiles.begin(), Ffiles.end());
    // std::sort(Sfiles.begin(), Sfiles.end());
    /* }}} */

    /* {{{ split V files in sequences -- what for ? */
    /* How many (j0-j1) ranges do we have for V files */
    vseq_t Vsequences;
    for(unsigned int i = 0 ; i < Vfiles.size() ; i++) {
        vseq_t::iterator it = Vsequences.find(Vfiles[i].seq_id());
        if (it == Vsequences.end()) {
            it = Vsequences.insert(make_pair(Vfiles[i].seq_id(), vseq_t::mapped_type())).first;
        }
        it->second.push_back(Vfiles[i]);
    }
    for(vseq_t::iterator it = Vsequences.begin(); it != Vsequences.end(); it++) {
        std::sort(it->second.begin(), it->second.end());
    }
    /* }}} */

    int nfailed = 0;

    check_V_files(Ac.get(), Vsequences, Cfiles, nfailed);

    /* Check A files using V, D, T, and R */
    if ((Tfiles.empty() || Rfiles.empty()) && !Dfiles.empty()) {
        fmt::print(stderr, FMT_STRING("{}"), "It makes no sense to provide Cd files and no Cr and Ct file\n");
        exit(EXIT_FAILURE);
    } else if (!Tfiles.empty() && !Rfiles.empty() && !Dfiles.empty()) {
        check_A_files(Ac.get(), Vfiles, Afiles, Dfiles, Rfiles.front(), Tfiles.front(), nfailed);
    }

    if (nfailed) {
        std::string diag = fmt::format(FMT_STRING("{} checks FAILED !!!!!!!!!!!!!!!!!\n"), nfailed);
        fmt::print(FMT_STRING("{}"), diag);
        fmt::print(stderr, FMT_STRING("{}"), diag);
        exit(EXIT_FAILURE);
    }


#if 0
    uint32_t * gxvecs = NULL;
    unsigned int nx = 0;
    if (!fake) {
        load_x(&gxvecs, bw->m, &nx, pi);
    } else {
        set_x_fake(&gxvecs, bw->m, &nx, pi);
    }

    if (!bw->skip_online_checks) {
        /* We do the dot product by working on the local vector chunks.
         * Therefore, we must really understand the check vector as
         * playing a role in the very same direction of the y vector!
         */
        mmt_vec_init(mmt, Ac, Ac_pi,
                check_vector, bw->dir, THREAD_SHARED_VECTOR, mmt->n[bw->dir]);
        if (tcan_print) { printf("Loading check vector..."); fflush(stdout); }
        mmt_vec_load(check_vector, CHECK_FILE_BASE, bw->interval,  mmt->n0[bw->dir]);
        if (tcan_print) { printf("done\n"); }
    }

    if (!bw->skip_online_checks) {
        ahead = Ac->alloc(nchecks, ALIGNMENT_ON_ALL_BWC_VECTORS);
    }

    /* We'll store all xy matrices locally before doing reductions. Given
     * the small footprint of these matrices, it's rather innocuous.
     */
    void * xymats;

    if (tcan_print) {
        printf("Each thread allocates %zd kb for the Ac matrices\n",
                Ac->vec_elt_stride(Ac, bw->m*bw->interval) >> 10);
    }
    if (!bw->skip_online_checks) {
        Ac->vec_set_zero(Ac, ahead, nchecks);
        AxAc->dotprod(Ac, Ac, ahead,
                mmt_my_own_subvec(check_vector),
                mmt_my_own_subvec(ymy[0]),
                mmt_my_own_size_in_items(ymy[0]));
    }

    for(int i = 0 ; i < bw->interval ; i++) {
        /* Compute the product by x */
        x_dotprod(Ac->vec_subvec(Ac, xymats, i * bw->m),
                gxvecs, bw->m, nx, ymy[0], 1);

        matmul_top_mul(mmt, ymy, timing);

        timing_check(pi, timing, s+i+1, tcan_print);
    }
    serialize(pi->m);

    /* See remark above. */
    pi_interleaving_flip(pi);
    pi_interleaving_flip(pi);

    if (!bw->skip_online_checks) {
        /* Last dot product. This must cancel ! */
        x_dotprod(ahead, gxvecs, nchecks, nx, ymy[0], -1);

        pi_allreduce(NULL, ahead, nchecks, mmt->pitype, BWC_PI_SUM, pi->m);
        if (!Ac->vec_is_zero(Ac, ahead, nchecks)) {
            printf("Failed check at iteration %d\n", s + bw->interval);
            exit(1);
        }
    }

    mmt_vec_untwist(mmt, ymy[0]);

    /* Now (and only now) collect the xy matrices */
    pi_allreduce(NULL, xymats,
            bw->m * bw->interval,
            mmt->pitype, BWC_PI_SUM, pi->m);

    if (pi->m->trank == 0 && pi->m->jrank == 0) {
        char * tmp;
        int rc;
        rc = asprintf(&tmp, A_FILE_PATTERN, ys[0], ys[1], s, s+bw->interval);
        FILE * f = fopen(tmp, "wb");
        rc = fwrite(xymats, Ac->vec_elt_stride(Ac, 1), bw->m*bw->interval, f);
        if (rc != bw->m*bw->interval) {
            fprintf(stderr, "Ayee -- short write\n");
            // make sure our input data won't be deleted -- this
            // chunk will have to be redone later, maybe the disk
            // failure is temporary (?)
        }
        fclose(f);
        free(tmp);
    }

    mmt_vec_save(ymy[0], v_name, s + bw->interval, unpadded);

    if (pi->m->trank == 0 && pi->m->jrank == 0)
        keep_rolling_checkpoints(v_name, s + bw->interval);

    serialize(pi->m);

    // reached s + bw->interval. Count our time on cpu, and compute the sum.
    timing_disp_collective_oneline(pi, timing, s + bw->interval, tcan_print, "check");

    free(gxvecs);
    free(v_name);

    for(int i = 0 ; i < mmt->nmatrices + nmats_odd ; i++) {
        mmt_vec_clear(mmt, ymy[i]);
    }
    free(ymy);
#endif

    return NULL;
}


int main(int argc, char * argv[])
{
    param_list pl;

    bw_common_init(bw, &argc, &argv);
    param_list_init(pl);

    param_list_usage_header(pl,
            "Usage: %s [options] -- [list of file names]\n"
            "\n"
            "File names are checked together (all relevant checks tying files to eachother within the provided argument list are done). The program tells which checks have been done\n"
            "Options are as follows. Note that not all are relevant to this program specifically:\n",
            argv[0]);

    bw_common_decl_usage(pl);
    /* declare local parameters and switches: none here (so far). */

    bw_common_parse_cmdline(bw, pl, &argc, &argv);

    bw_common_interpret_parameters(bw, pl);

    if (param_list_warn_unused(pl)) {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (!rank) param_list_print_usage(pl, bw->original_argv[0], stderr);
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    check_prog(pl, argc, argv);

    param_list_clear(pl);
    bw_common_clear(bw);

    return 0;
}

