/* This standalone program checks matrix products produced by lingen
   (as stored in the cp/ subdirectory).

  Example: check -prime ... -dim 13 -k 3 pi.3.2.127 pi.3.127.251 pi.2.2.251

  Optional arguments: [-seed xxx] [-v]
*/

#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <climits>


#include <iostream>
#include <fstream>
#include <vector>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include <unistd.h>
#include <gmp.h>
#include "fmt/format.h"

#include "gmp_aux.h"
#include "cxx_mpz.hpp"
#include "gmp-hacks.h"
#include "lingen_bmstatus.hpp"
#include "lingen_checkpoints.hpp"
#include "lingen_hints.hpp"
#include "macros.h"
#include "omp_proxy.h"
#include "params.h"
#include "select_mpi.h"
#include "subdivision.hpp"

/* define WARNING to get warning for non-zero padding coefficients */
// #define WARNING

static struct
{
    unsigned int m, n;
} bw_parameters;

static cxx_mpz prime;          /* prime modulus */
static unsigned long lingen_p; /* number of limbs per coefficient */
static int mpi_k = 1;          /* matrix is cut in k x k submatrices */
static cxx_gmp_randstate state;
static int verbose = 0;
static unsigned long seed;
static unsigned long global_batch = 128;
static unsigned int restrict_E = UINT_MAX;
static lingen_hints hints;

struct matrix
{
    unsigned long nrows; /* matrix dimension */
    unsigned long ncols; /* matrix dimension */
    std::vector<cxx_mpz> coeff;

    matrix(unsigned long nrows,
           unsigned int ncols)
      : nrows(nrows)
      , ncols(ncols)
      , coeff(nrows * ncols)

    {}
    void zero()
    {
        for (auto& x : coeff)
            mpz_set_ui(x, 0);
    }
};

class matrix_reader
{
    std::string stem;
    unsigned long nrows; /* matrix dimension */
    unsigned long ncols; /* matrix dimension */
    unsigned long k;     /* matrix is cut into k x k submatrices */
    unsigned long deg;
    bool reverse = false;
    std::vector<std::ifstream> files;
    std::vector<bool> warned_padding;

    std::string get_filename_ij(unsigned long k, unsigned long i, unsigned long j) const {
        int nij = k * i + j;
        std::string filename;
        if (k > 1) {
            filename = stem + fmt::format(".{}.data", nij);
        } else {
            filename = stem + ".single.data";
        }
        return filename;
    }

    bool has_single_file_data() const {
        std::string const filename = get_filename_ij(1, 0, 0);
        return access(filename.c_str(), R_OK) == 0;
    }
    unsigned int has_mpi_data() const {
        /* return 0 if no mpi data is found. If consistent mpi data is
         * found, return the corresponding splitting */
        std::string filename;
        unsigned int c = 0;
        for (unsigned int k = 1 ; k < 64 ; k++) {
            for (unsigned int s = 0 ; c < k * k ; c++, s++) {
                filename = stem + fmt::format(".{}.data", c);
                fmt::print(stderr, "test {}\n", filename.c_str());
                if (access(filename.c_str(), R_OK) != 0) {
                    if (s == 0) {
                        if (k-1 == 1) {
                            throw std::runtime_error(fmt::format("weird: we have 1-node mpi data for checkpoint {}, which in theory we shouldn't produce\n", filename));
                            /* anyway it's going to fail with the present
                             * code, because the meaning of k==1 is
                             * ambiguous */
                        }
                        return k-1;
                    } else {
                        throw std::runtime_error("got non-square mpi data for checkpoint " + filename);
                    }
                }
            }
        }
        return 0;
    }

public:
    matrix_reader(unsigned long nrows,
            unsigned int ncols,
            unsigned long deg,
            std::string const& stem,
            bool reverse)
        : stem(stem)
          , nrows(nrows)
          , ncols(ncols)
          , deg(deg)
          , reverse(reverse)
          , warned_padding(nrows * ncols, false)
    {
        if ((k = has_mpi_data()) != 0) {
            /* do nothing, it's all fine! */
        } else if (has_single_file_data()) {
            /* use the single file data instead. */
            k = 1;
        } else {
            throw std::runtime_error("neither single-file nor mpi data found for " + stem);
        }

        for (unsigned long i = 0; i < k; i++) {
            for (unsigned long j = 0; j < k; j++) {
                std::string const filename = get_filename_ij(k, i, j);
                files.emplace_back(filename, std::ios_base::in);
                if (!files.back().good())
                    throw std::runtime_error("cannot open " + filename);
            }
        }
    }

    private:
    /* read 1 coefficient (degree idx, or deg-idx in reverse case),
     * from block (block_i, block_j), multiply it by x, and accumulate to
     * M at the right place.
     */
    void read_n_accumulate(matrix& M,
            unsigned long block_i,
            unsigned long block_j,
            cxx_mpz & x,
            cxx_mpz const & cx,
            unsigned int idx,
            unsigned int n
            )
    {
        int const nij = k * block_i + block_j;
        ASSERT_ALWAYS(nrows == M.nrows);
        ASSERT_ALWAYS(ncols == M.ncols);
        unsigned int i0, i1;
        unsigned int j0, j1;
        subdivision const Srows(M.nrows, k);
        subdivision const Scols(M.ncols, k);
        std::tie(i0, i1) = Srows.nth_block(block_i);
        std::tie(j0, j1) = Scols.nth_block(block_j);
        unsigned int const nix = Srows.block_size_upper_bound();
        unsigned int const njx = Scols.block_size_upper_bound();
        cxx_mpz tmp;
        mpz_realloc2(tmp, lingen_p * mp_bits_per_limb);
        std::ifstream& F(files[nij]);
        for(unsigned int s = 0 ; s < n && (idx + s <= deg) ; s++) {
            unsigned long cidx = idx + s;
            if (reverse)
                cidx = deg - cidx;
            if (reverse) {
                F.seekg(cidx * nix * njx * lingen_p * sizeof(mp_limb_t));
                ASSERT_ALWAYS(F.good());
            }
            for (unsigned int di = 0; di < nix; di++) {
                for (unsigned int dj = 0; dj < njx; dj++) {
                    ASSERT_ALWAYS((unsigned long) ALLOC(tmp) >= lingen_p);
                    SIZ(tmp) = lingen_p;
                    size_t const sz = lingen_p * sizeof(mp_limb_t);
                    {
                    F.read((char*)PTR(tmp), sz);
                    bool const good_read = F.good();
                    if (!good_read) {
                        fmt::print(stderr, "read error on file ({},{}) [{}] when reading coefficient of degree {}, local position ({},{}). File offset of coefficient is {}\n",
                                block_i,
                                block_j,
                                get_filename_ij(k, block_i, block_j),
                                cidx,
                                di,
                                dj,
                                cidx * nix * njx * lingen_p * sizeof(mp_limb_t)
                                );
                    }
                    ASSERT_ALWAYS(good_read);
                    }
                    MPN_NORMALIZE(PTR(tmp), SIZ(tmp));
                    if (i0 + di < i1 && j0 + dj < j1) {
                        cxx_mpz& dst(M.coeff[(i0 + di) * ncols + j0 + dj]);
                        mpz_mul(tmp, tmp, x);
                        mpz_add(tmp, tmp, dst);
                        mpz_mod(dst, tmp, prime);
                    } else {
#ifdef WARNING
                        if (SIZ(tmp) != 0 &&
                                !warned_padding[(i0 + di) * ncols + j0 + dj]) {
                            fmt::print(
                                    stderr,
                                    "Warning, padding coefficient {},{} is not zero\n",
                                    i0 + di,
                                    j0 + dj);
                            warned_padding[(i0 + di) * ncols + j0 + dj] = true;
                        }
#endif
                    }
                }
            }
            ASSERT(
                F.tellg() == (ssize_t) ((cidx+1) * nix * njx * lingen_p * sizeof(mp_limb_t)));
            mpz_mul(x, x, cx);
            mpz_mod(x, x, prime);
        }
    }

    public:
    /* read 1 coefficient (degree idx, or deg-idx in reverse case),
     * multiply it by x, and accumulate to M at the right place.
     */
    void read1_accumulate(matrix& M, cxx_mpz & x, cxx_mpz const & cx, unsigned int idx)
    {
        read_n_accumulate(M, x, cx, idx, 1);
    }
    void read_n_accumulate(matrix& M, cxx_mpz & x, cxx_mpz const & cx, unsigned int idx, unsigned int n)
    {
#ifdef HAVE_OPENMP
#pragma omp parallel 
#endif
        {
            cxx_mpz xpriv;
            bool adv = false;
#ifdef HAVE_OPENMP
#pragma omp for collapse(2)
#endif
            for (unsigned long i = 0; i < k; i++)
                for (unsigned long j = 0; j < k; j++) {
                    xpriv = x;
                    read_n_accumulate(M, i, j, xpriv, cx, idx, n);
                    adv = true;
                }
#ifdef HAVE_OPENMP
#pragma omp barrier
#endif
            if (adv) {
                /* Can't do conditional single, apparently */
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
                {
                    /* only one thread stores the final value. */
                    x = std::move(xpriv);
                }
            }
        }
    }
};

static void mpz_urandomm_nz(mpz_ptr a, cxx_gmp_randstate & state, mpz_srcptr prime)
{
    ASSERT_ALWAYS(mpz_cmp_ui(prime, 1) > 0);
    do {
        mpz_urandomm(a, state, prime);
    } while (mpz_cmp_ui(a, 0) == 0);
}

/* return a vector of n random numbers mod p */
static void
fill_random(std::vector<cxx_mpz>& u)
{
    for (auto& a : u)
        mpz_urandomm_nz(a, state, prime);
}

template<bool is_binary>
typename lingen_checkpoint<is_binary>::header_info read_cp_aux(std::string const& prefix)
{
    auto filename = prefix + ".aux";
    typename lingen_checkpoint<is_binary>::header_info h;
    if (!(std::ifstream(filename) >> h))
        throw typename lingen_checkpoint<is_binary>::invalid_aux_file(fmt::format(
                    "Reading header from {} failed", filename));
    return h;
}

/* read a matrix of dimension n, divided into kxk submatrices.
 * return the evaluation of the matrix polynomial at x.
 * */
template<bool is_binary>
static matrix
read_matrix(std::string const & s,
            unsigned long nrows,
            unsigned long ncols,
            cxx_mpz const& x)
{
    auto const cp = read_cp_aux<is_binary>(s);
    unsigned long const deg = cp.ncoeffs - 1;
    matrix M(nrows, ncols);
    matrix_reader R(nrows, ncols, deg, s, false);
    cxx_mpz x_power_k;
    mpz_set_ui(x_power_k, 1);
    const unsigned long batch = global_batch;
    for (unsigned long k = 0; k <= deg; k += batch) {
        /* invariant: x_power_k = x^k mod prime */
        R.read_n_accumulate(M, x_power_k, x, k, batch);
    }
    return M;
}

/* w <- v*M evaluated at x and modulo p */
static void
mul_left(std::vector<cxx_mpz>& w,
         std::vector<cxx_mpz> const& v,
         matrix const& M)
{
    ASSERT_ALWAYS(v.size() == M.nrows);
    ASSERT_ALWAYS(w.size() == M.ncols);
    cxx_mpz tmp;
    for (unsigned long j = 0; j < M.ncols; j++) {
        mpz_set_ui(w[j], 0);
        for (unsigned long i = 0; i < M.nrows; i++) {
            /* w[j] += v[i]*M[i,j] */
            mpz_mul(tmp, v[i], M.coeff[i * M.ncols + j]);
            mpz_add(w[j], w[j], tmp);
        }
        mpz_mod(w[j], w[j], prime);
    }
}

static std::vector<cxx_mpz> operator*(std::vector<cxx_mpz> const& v, matrix const& M)
{
    std::vector<cxx_mpz> res(M.ncols);
    mul_left(res, v, M);
    return res;
}

/* w <- M*v evaluated */
static void
mul_right(std::vector<cxx_mpz>& w,
          matrix const& M,
          std::vector<cxx_mpz> const& v)
{
    ASSERT_ALWAYS(v.size() == M.ncols);
    ASSERT_ALWAYS(w.size() == M.nrows);
    unsigned long const n = M.ncols;
    cxx_mpz tmp;
    for (unsigned long i = 0; i < M.nrows; i++) {
        mpz_set_ui(w[i], 0);
        for (unsigned long j = 0; j < n; j++) {
            /* w[i] += M[i,j]*v[j] */
            mpz_mul(tmp, M.coeff[i * M.ncols + j], v[j]);
            mpz_add(w[i], w[i], tmp);
        }
        mpz_mod(w[i], w[i], prime);
    }
}

static std::vector<cxx_mpz> operator*(matrix const& M, std::vector<cxx_mpz> const& v)
{
    std::vector<cxx_mpz> res(M.nrows);
    mul_right(res, M, v);
    return res;
}

static void
add_scalar_product(cxx_mpz& res,
                   std::vector<cxx_mpz> const& u,
                   std::vector<cxx_mpz> const& v)
{
    ASSERT_ALWAYS(u.size() == v.size());
    for (unsigned long i = 0; i < u.size(); i++)
        mpz_addmul(res, u[i], v[i]);
    mpz_mod(res, res, prime);
}

static cxx_mpz
scalar_product(std::vector<cxx_mpz> const& u, std::vector<cxx_mpz> const& v)
{
    ASSERT_ALWAYS(u.size() == v.size());
    cxx_mpz res;
    mpz_set_ui(res, 0);
    add_scalar_product(res, u, v);
    return res;
}

static void
declare_usage(cxx_param_list& pl)
{
    param_list_usage_header(
      pl,
      "Usage: lingen_verify_checkpoints [options] -- [list of file names]\n"
      "\n"
      "The list of file names can have one of the following formats:\n"
      " - pi0 pi1 pi2 : check that pi0*pi1==pi2\n"
      " - E pi : check that E*pi=O(x^length(E))\n"
      "Options are as follows.\n");
    param_list_decl_usage(pl, "prime", "characteristic of the base field");
    param_list_decl_usage(pl, "mpi", "mpi geometry");
    param_list_decl_usage(pl, "seed", "random seed");
    param_list_decl_usage(pl, "m", "block Wiedemann parameter m");
    param_list_decl_usage(pl, "n", "block Wiedemann parameter n");
    param_list_decl_usage(pl, "v", "More verbose output");
    param_list_decl_usage(pl, "restrict_E", "(E*pi check only) check only mod this power of X");
    param_list_decl_usage(pl, "sanity-check", "Do sanity check on the checkpoint auxiliary file");
    param_list_decl_usage(pl, "tuning_schedule_filename",
                "load tuning schedule from this file (sanity checks only)");
    param_list_decl_usage(pl, "cpdir",
                "load checkpoints from this directory");
}

template<bool is_binary>
static int
do_check_pi(std::string const & pi_left_filename,
            std::string const & pi_right_filename,
            std::string const & pi_filename)
{
    int ret;
    unsigned long const nrows = bw_parameters.m + bw_parameters.n;
    unsigned long const ncols = bw_parameters.m + bw_parameters.n;

    std::vector<cxx_mpz> u(nrows);
    std::vector<cxx_mpz> v(nrows);

    fill_random(u);
    fill_random(v);

    cxx_mpz x;
    mpz_urandomm_nz(x, state, prime); /* random variable */
    // fmt::print("Using seed {}\n", seed);
    if (verbose)
        fmt::print("x={}\n", x);

    auto const cp = read_cp_aux<is_binary>(pi_filename);

    std::string const check_name =
      fmt::format("check (seed={}, depth {}, t={}, pi_left*pi_right=pi)",
                   seed,
                   cp.level,
                   cp.t);

    std::vector<cxx_mpz> u_times_piab(ncols);
    std::vector<cxx_mpz> pibc_times_v(nrows);
    std::vector<cxx_mpz> piac_times_v(nrows);

    /* we used to have omp sections here. It feels wrong.  */
    {
        {
            matrix const Mab = read_matrix<is_binary>(pi_left_filename, nrows, ncols, x);
            mul_left(u_times_piab, u, Mab);
        }

        {
            matrix const Mbc = read_matrix<is_binary>(pi_right_filename, nrows, ncols, x);
            mul_right(pibc_times_v, Mbc, v);
        }

        {
            matrix const Mac = read_matrix<is_binary>(pi_filename, nrows, ncols, x);
            mul_right(piac_times_v, Mac, v);
        }
    }

    cxx_mpz res_left = scalar_product(u_times_piab, pibc_times_v);
    cxx_mpz res_right = scalar_product(u, piac_times_v);

    if (mpz_cmp(res_left, res_right) != 0) {
        fmt::print(stderr, "FAILED {}\n", check_name);
        fmt::print("res_left  = {}\n", res_left);
        fmt::print("res_right = {}\n", res_right);
        ret = 0;
    } else {
        fmt::print("ok {}\n", check_name);
        ret = 1;
    }

    return ret;
}

/* check that E*pi = O(x^length(E)) at a given level. */
template<bool is_binary>
static int
do_check_E_short(std::string const& E_filename, std::string const& pi_filename)
{
    int ret;
    unsigned int const m = bw_parameters.m;
    unsigned int const n = bw_parameters.n;

    std::vector<cxx_mpz> u(m);
    fill_random(u);
    std::vector<cxx_mpz> v(m + n);
    fill_random(v);

    cxx_mpz x;
    mpz_urandomm_nz(x, state, prime); /* random variable */
    if (verbose)
        fmt::print("x={}\n", x);

    /* Note that all the useful info is in the aux file for pi, really.
     * The one for E is stored at t0, and is not really useful.
     */
    auto const cp = read_cp_aux<is_binary>(pi_filename);

    unsigned long const deg_pi = cp.ncoeffs - 1;
    unsigned long const t = cp.t;
    unsigned long const t0 = cp.t0;
    unsigned long const t1 = cp.t1;
    unsigned long deg_E = t - t0 - 1;
    deg_E = MIN(deg_E, (unsigned long) restrict_E);

    std::string check_name = fmt::format(
            "check (seed={}, depth {}, t={}, E*pi=O(X^{}))",
            seed, cp.level, cp.t0, deg_E);

    if (t < t1)
        check_name += " [truncated cp at end]";


    matrix pi(m + n, m + n);
    matrix_reader Rpi(m + n, m + n, deg_pi, pi_filename, true);
    matrix E(m, m + n);
    matrix_reader RE(m, m + n, deg_E, E_filename, false);

    /* We'll compute the evaluation at x of the short product of E*pi,
     * capped to degree deg_E
     *
     * Note that pi has degree smaller than E.
     *
     * We use the following expression of the result.
     *
     * res = (E_0 + ... + E_{deg_E} * x^deg_E) * pi_0
     *     + (E_0 + ... + E_{deg_E-1} * x^{deg_E-1}) * pi_1 * x
     *     + (E_0 + ... + E_{deg_E-2} * x^{deg_E-2}) * pi_2 * x^2
     *     + ...
     *     + (E_0 + ... + E_{deg_E-deg_pi} * x^{deg_E-deg_pi}) * pi_deg_pi *
     * x^deg_pi and the evaluation will be performed from the last line to the
     * first. We get a similarly significant result by multiplying
     * everything by x^(-deg_pi), which gives the following formula:
     *
     * res = (E_0 + ... + E_{deg_E-deg_pi} * x^{deg_E-deg_pi}) * pi_deg_pi
     *     + ...
     *     + (E_0 + ... + E_{deg_E-2} * x^{deg_E-2}) * pi_2 * xinv^{deg_pi-2}
     *     + (E_0 + ... + E_{deg_E-1} * x^{deg_E-1}) * pi_1 * xinv^{deg_pi-1}
     *     + (E_0 + ... + E_{deg_E} * x^deg_E) * pi_0 * xinv^{deg_pi}
     */

    cxx_mpz xinv;
    mpz_invert(xinv, x, prime);

    cxx_mpz res;

    cxx_mpz x_inc;
    mpz_set_ui(x_inc, 1);

    const unsigned long batch = global_batch;
    for (unsigned long k = 0; deg_pi + k < deg_E; k += batch) {
        /* invariant: x_power_k = x^k mod prime */
        RE.read_n_accumulate(E, x_inc, x, k, MIN(batch, deg_E - deg_pi- k));
    }
    cxx_mpz x_dec;
    mpz_set_ui(x_dec, 1);
    /* we don't really know how to batch this one, do we ? */
    for (unsigned long k = deg_pi - MIN(deg_pi, deg_E); k <= deg_pi; k++) {
        RE.read1_accumulate(E, x_inc, x, deg_E - deg_pi + k);
        pi.zero();
        Rpi.read1_accumulate(pi, x_dec, xinv, k);
        std::vector<cxx_mpz> const u_E = u * E;
        std::vector<cxx_mpz> const pi_v = pi * v;
        add_scalar_product(res, u_E, pi_v);
    }

    if (mpz_cmp_ui(res, 0) != 0) {
        fmt::print(stderr, "FAILED {}\n", check_name);
        fmt::print("res  = {}\n", res);
        ret = 0;
    } else {
        fmt::print("ok {}\n", check_name);
        ret = 1;
    }

    return ret;
}

template<bool is_binary>
static int sanity_check(std::string const & filename)
{
    auto const cp = read_cp_aux<is_binary>(filename);
    bmstatus<is_binary> bm(bw_parameters.m,bw_parameters.n, prime);
    bm.set_t0(cp.t0);
    bm.hints= hints;
    lingen_checkpoint<is_binary> lcp(bm, cp.t0, cp.t1, mpi_k > 1, filename);
    size_t Xsize;
    try {
        if (!lcp.load_aux_file(Xsize)) {
            fmt::print(stderr, "{} is missing\n", lcp.auxfile);
            return false;
        }
        int const sdata_ok = access(lcp.sdatafile.c_str(), R_OK) == 0;
        int const gdata_ok = lcp.rank || access(lcp.gdatafile.c_str(), R_OK) == 0;
        fmt::print(stderr, "scattered datafile {}: {}\n", lcp.sdatafile, sdata_ok ? "ok" : "not found");
        fmt::print(stderr, "gathered datafile {}: {}\n", lcp.gdatafile, gdata_ok ? "ok" : "not found");
        return sdata_ok || gdata_ok;
    } catch (typename lingen_checkpoint<is_binary>::invalid_aux_file const & inv) {
        fmt::print(stderr, "Invalid checkpoint aux file {} [{}]\n",
                lcp.auxfile, inv.what());
        return false;
    }
}

template<bool is_binary>
static int all_tests(cxx_param_list & pl,
        std::vector<std::pair<const char **, int> > const & todo,
        std::string const & cpdir)
{
    const char * tmp;
    int ret = 0;
    /* it doesn't seem to make sense to compare
     *      MP(E_{level,t0}, pi_{level+1,t0,t}) and E_{level+1,t}
     * ? This is covered by the test of the short product
     *      E_{level,t0} * pi_{level,t0}
     * since we can validate the fact that pi_{level,t0} is the
     * product of the two pi matrices at the level below, and that
     * the short product E_{level+1,t0}*pi_{level+1,t0} is also zero.
     */
    if ((tmp = param_list_lookup_string(pl, "sanity-check")) != NULL) {
        ret = sanity_check<is_binary>(tmp);
    } else {
        for(auto const & x : todo) {
            const char ** argv = x.first;
            int argc = x.second;
            if (argc == 3) {
                ret = do_check_pi<is_binary>(
                        cpdir + argv[0],
                        cpdir + argv[1],
                        cpdir + argv[2]);
            } else if (argc == 2) {
                ret = do_check_E_short<is_binary>(
                        cpdir + argv[0],
                        cpdir + argv[1]);
            }
            if (!ret) break;
        }
    }
    return ret;
}


// coverity[root_function]
int main(int argc, char const * argv[])
{
    /* We're not really mpi, but we link code that _is_ mpi */
    MPI_Init(&argc, (char ***) &argv);
    cxx_param_list pl;

    seed = getpid();

    declare_usage(pl);

    const char* argv0 = argv[0];

    std::vector<std::pair<const char **, int> > todo;
    std::string cpdir;

    param_list_configure_switch(pl, "-v", &verbose);
    for (argc--, argv++; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }
        if (strcmp(argv[0], "--") == 0) {
            argc--, argv++;
            if (argc == 2 || argc == 3) {
                todo.emplace_back(argv, argc);
                break;
            } else if (argc >= 3 && strcmp(argv[2], "--") == 0) {
                todo.emplace_back(argv, 2);
                argc--, argv++;
                argc--, argv++;
                continue;
            } else if (argc >= 4 && strcmp(argv[3], "--") == 0) {
                todo.emplace_back(argv, 3);
                argc--, argv++;
                argc--, argv++;
                argc--, argv++;
                continue;
            } else {
                fmt::print(stderr, "bad argument list\n");
                exit(EXIT_FAILURE);
            }
        }
        fmt::print(stderr, "Unhandled parameter {}\n", argv[0]);
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    param_list_parse(pl, "cpdir", cpdir);
    if (!cpdir.empty() && cpdir.back() != '/')
            cpdir += '/';

    if (!param_list_parse_mpz(pl, "prime", (mpz_ptr)prime)) {
        fprintf(stderr, "Missing parameter: prime\n");
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }
    int mpi_dims[2] = { 0, 0 };
    if (param_list_parse_int_and_int(pl, "mpi", mpi_dims, "x")) {
        ASSERT_ALWAYS(mpi_dims[0] == mpi_dims[0]);
        mpi_k = mpi_dims[0];
    }
    param_list_parse_ulong(pl, "seed", &seed);
    if (!param_list_parse_uint(pl, "m", &bw_parameters.m)) {
        fmt::print(stderr, "Missing parameter: m\n");
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }
    if (!param_list_parse_uint(pl, "n", &bw_parameters.n)) {
        fmt::print(stderr, "Missing parameter: n\n");
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }
    param_list_parse_uint(pl, "restrict_E", &restrict_E);
    param_list_lookup_string(pl, "sanity-check");

    const char * tmp;

    if ((tmp = param_list_lookup_string(pl, "tuning_schedule_filename")) != NULL) {
        std::ifstream is(tmp);
        if (is && is >> hints) {
            /* This one _always_ goes to stdout */
            std::cout << fmt::format("# Read tuning schedule from {}\n", tmp);
        } else {
            std::cerr << fmt::format("# Failed to read tuning schedule from {}\n", tmp);
            hints = lingen_hints();
        }
    }

    if (param_list_warn_unused(pl)) {
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    gmp_randseed_ui(state, seed);

    lingen_p = mpz_size(prime);

#ifdef LINGEN_BINARY
    constexpr bool is_binary = true;
#else
    constexpr bool is_binary = false;
#endif

    int ret = all_tests<is_binary>(pl, todo, cpdir);

    MPI_Finalize();

    return ret ? EXIT_SUCCESS : EXIT_FAILURE;
}
/* vim: set sw=4 sta et: */
