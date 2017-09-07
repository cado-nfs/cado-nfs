#include "cado.h"
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <vector>
#include <algorithm>
#include <stdexcept>
#include <map>
#include <string>
#include <sstream>

#include "bwc_config.h"
#include "parallelizing_info.h"
#include "matmul_top.h"
#include "select_mpi.h"
#include "params.h"
#include "xvectors.h"
#include "portability.h"
#include "misc.h"
#include "bw-common.h"
#include "async.h"
#include "xdotprod.h"
#include "rolling.h"
#include "mpfq/mpfq.h"
#include "mpfq/mpfq_vbase.h"
#include "cheating_vec_init.h"


using namespace std;

// arguments we need
// m
// n
// prime
//
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
    unsigned int stretch;
    Cfile(const char * x) : string(x) {
        int rc = sscanf(my_basename(x), "C.%u", &stretch);
        if (rc != 1) throw std::runtime_error("want C.%u");
    }
    inline bool operator<(Cfile const& o) {
        return stretch < o.stretch;
    }
};

struct Vfile : public string {
    unsigned int j0, j1, n;
    int checks;
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
    int checks;
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
    int checks;
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
    int checks;
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

void vec_alloc(mpfq_vbase_ptr A, void *& z, size_t vsize)
{
    cheating_vec_init(A, &z, vsize);
    A->vec_set_zero(A, z, vsize);
}

void vec_free(mpfq_vbase_ptr A, void *& z, size_t vsize)
{
    cheating_vec_clear(A, &z, vsize);
}

void vec_read(mpfq_vbase_ptr A, void * z, string const & v, size_t vsize, const char * prefix = NULL)
{
    if (prefix) printf("%sload %s\n", prefix, v.c_str());
    FILE * f = fopen(v.c_str(), "rb");
    ASSERT_ALWAYS(f != NULL);
    int rc = fread(z, A->vec_elt_stride(A, 1), vsize, f);
    ASSERT_ALWAYS(rc >= 0);
    ASSERT_ALWAYS((size_t) rc == vsize);
    fclose(f);
}

size_t vec_items(mpfq_vbase_ptr A, string const & v)
{
    struct stat sbuf[1];
    int rc = stat(v.c_str(), sbuf);
    ASSERT_ALWAYS(rc == 0);
    return sbuf->st_size / A->vec_elt_stride(A, 1);
}

/* This is *not* a parallel program, so we depart significantly from the
 * way programs such as krylov or mksol are written.
 *
 */
void * check_prog(param_list pl MAYBE_UNUSED, int argc, char * argv[])
{
    int withcoeffs = mpz_cmp_ui(bw->p, 2) > 0;
    int nchecks = withcoeffs ? NCHECKS_CHECK_VECTOR_GFp : NCHECKS_CHECK_VECTOR_GF2;
    mpfq_vbase Ac;
    mpfq_vbase_oo_field_init_byfeatures(Ac, 
            MPFQ_PRIME_MPZ, bw->p,
            MPFQ_GROUPSIZE, nchecks,
            MPFQ_DONE);

    vector<Cfile> Cfiles;
    vector<Vfile> Vfiles;
    vector<Afile> Afiles;
    // vector<Ffile> Ffiles;
    // vector<Sfile> Sfiles;

    for(int i = 0 ; i < argc ; i++) {
        try {
            const char * p = my_basename(argv[i]);
            if (*p && p[strlen(p)-1] == '~')
                continue;
            switch(*p) {
                case 'C': Cfiles.push_back(argv[i]); break;
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
                  fprintf(stderr, "File name not recognized: %s\n", argv[i]);
                  exit(EXIT_FAILURE);
            }
        } catch (std::runtime_error& e) {
            fprintf(stderr, "Parse error on %s: %s\n", argv[i], e.what());
            exit(EXIT_FAILURE);
        }
    }

    std::sort(Cfiles.begin(), Cfiles.end());
    std::sort(Afiles.begin(), Afiles.end());
    std::sort(Vfiles.begin(), Vfiles.end());
    // std::sort(Ffiles.begin(), Ffiles.end());
    // std::sort(Sfiles.begin(), Sfiles.end());

    /* How many (j0-j1) ranges do we have for V files */
    typedef std::map<pair<unsigned int, unsigned int>, vector<Vfile> > vseq_t;
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

    /* Check that C files have consistent size */
    size_t vsize = 0;

    for(unsigned int cc = 0 ; cc < Cfiles.size() ; cc++) {
        Cfile& C(Cfiles[cc]);

        size_t items = vec_items(Ac, C);

        if (vsize == 0) {
            vsize = items;
        } else if (vsize != items) {
            fprintf(stderr, "File sizes disagree for %s (%zu items) and %s (%zu items)\n",
                    Cfiles[0].c_str(), vsize,
                    C.c_str(), items);
            exit(EXIT_FAILURE);
        }
    }
    printf("C files have %zu coordinates\n", vsize);

    int nok=0;

    for(unsigned int i0 = 0 ; i0 < Cfiles.size() - 1 ; i0++) {
        Cfile& C_i0(Cfiles[i0]);
        void * Cv_i0;
        vec_alloc(Ac, Cv_i0, vsize);

        int has_read_Cv_i0 = 0;

        for(unsigned int i1 = i0 + 1 ; i1 < Cfiles.size() ; i1++) {
            Cfile& C_i1(Cfiles[i1]);
            const char * c = C_i1.c_str();

            printf("Doing checks for distance %d using %s and %s\n",
                    C_i1.stretch - C_i0.stretch,
                    C_i0.c_str(), C_i1.c_str()
                    );

            void * Cv_i1;
            vec_alloc(Ac, Cv_i1, vsize);

            int has_read_Cv_i1 = 0;

            /* {{{ check all V files together */
            for(vseq_t::iterator it = Vsequences.begin(); it != Vsequences.end(); it++)
            {
                vector<Vfile>& Vs(it->second);

                printf(" checks on V files for sequence %u-%u\n",
                        it->first.first, it->first.second);

                mpfq_vbase Av;
                mpfq_vbase_oo_field_init_byfeatures(Av, 
                        MPFQ_PRIME_MPZ, bw->p,
                        MPFQ_GROUPSIZE, it->first.second - it->first.first,
                        MPFQ_DONE);


                mpfq_vbase_tmpl AvxAc;
                mpfq_vbase_oo_init_templates(AvxAc, Av, Ac);

                /* {{{ Check that all V files here have the proper size */
                for(unsigned int i = 0 ; i < Vs.size() ; i++) {
                    size_t items = vec_items(Av, Vs[i]);
                    if (items != vsize) {
                        fprintf(stderr, "%s has %zu coordinates, different from expected %zu\n", Vs[i].c_str(), items, vsize);
                        exit(EXIT_FAILURE);
                    }
                }
                /* }}} */

                void * Vv;
                vec_alloc(Av, Vv, vsize);
                unsigned int Vv_iter = UINT_MAX;

                void * dotprod_scratch[2];
                vec_alloc(Av, dotprod_scratch[0], nchecks);
                vec_alloc(Av, dotprod_scratch[1], nchecks);

                unsigned int j = 0;
                for(unsigned int i = 0 ; i < Vs.size() ; i++) {
                    for(j = i + 1 ; j < Vs.size() ; j++) {
                        if (Vs[j].n + C_i0.stretch == Vs[i].n + C_i1.stretch)
                            break;
                    }
                    if (j == Vs.size()) continue;

                    if (!has_read_Cv_i0++)
                        vec_read(Ac, Cv_i0, C_i0, vsize, " ");

                    if (!has_read_Cv_i1++)
                        vec_read(Ac, Cv_i1, C_i1, vsize, " ");

                    Vs[j].checks++;
                    const char * vi = Vs[i].c_str();
                    const char * vj = Vs[j].c_str();
                    if (strncmp(vi, c, my_basename(c) - c) == 0) {
                        vi += my_basename(c) - c;
                    }
                    if (strncmp(vj, c, my_basename(c) - c) == 0) {
                        vj += my_basename(c) - c;
                    }
                    printf("  check %s against %s\n", vi, vj);
                    if (Vs[i].n != Vv_iter) {
                        vec_read(Ac, Vv, Vs[i].c_str(), vsize, "   ");
                        Vv_iter = Vs[i].n;
                    }

                    Av->vec_set_zero(Av, dotprod_scratch[0], nchecks);

                    /* compute the dot product */
                    AvxAc->dotprod(Av, Ac, 
                            dotprod_scratch[0],
                            Cv_i1, Vv, vsize);

                    vec_read(Ac, Vv, Vs[j].c_str(), vsize, "   ");
                    Vv_iter = Vs[j].n;

                    AvxAc->dotprod(Av, Ac, 
                            dotprod_scratch[1],
                            Cv_i0, Vv, vsize);

                    int cmp = Av->vec_cmp(Av, dotprod_scratch[0], dotprod_scratch[1], nchecks);

                    printf("  check %s against %s -> %s\n",
                            vi, vj, cmp == 0 ? "ok" : "NOK NOK NOK NOK NOK");

                    if (cmp != 0) {
                        nok++;
                        fprintf(stderr, " check %s against %s -> %s\n",
                                vi, vj, cmp == 0 ? "ok" : "NOK NOK NOK NOK NOK");
                    }
                }
                vec_free(Av, dotprod_scratch[0], nchecks);
                vec_free(Av, dotprod_scratch[1], nchecks);
                vec_free(Ac, Vv, vsize);

                Av->oo_field_clear(Av);
            }
            /* }}} */

            /* {{{ check A files */

            /* }}} */

            vec_free(Ac, Cv_i1, vsize);
        }
        vec_free(Ac, Cv_i0, vsize);
    }

    if (nok) {
        printf("%d checks FAILED !!!!!!!!!!!!!!!!!\n", nok);
        fprintf(stderr, "%d checks FAILED !!!!!!!!!!!!!!!!!\n", nok);
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
        cheating_vec_init(Ac, &ahead, nchecks);
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
    }

    free(gxvecs);
    free(v_name);

    for(int i = 0 ; i < mmt->nmatrices + nmats_odd ; i++) {
        mmt_vec_clear(mmt, ymy[i]);
    }
    free(ymy);
#endif

    Ac->oo_field_clear(Ac);

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
        param_list_print_usage(pl, bw->original_argv[0], stderr);
        exit(EXIT_FAILURE);
    }

    check_prog(pl, argc, argv);

    param_list_clear(pl);
    bw_common_clear(bw);

    return 0;
}

