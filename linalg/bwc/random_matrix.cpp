#include "cado.h" // IWYU pragma: keep

#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cfloat>
#include <cstring>
#include <cstdint>
#include <cmath>
#include <ctime>

#include <algorithm>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "cxx_mpz.hpp"
#include "gmp_aux.h"
#include "macros.h"
#include "misc.h"
#include "parallelizing_info.hpp"
#include "portability.h"
#include "params.h"
#include "random_distributions.hpp"
#include "matrix_u32.hpp"
#include "utils_cxx.hpp"

/* the files below are not useful for the standalone program */
#ifndef WANT_MAIN
#include "random_matrix.hpp"
#include "balancing.hpp"
#include "verbose.h"
#endif

#ifdef WANT_MAIN
#include <cinttypes>
#endif

/* The random generation works as follows.
 *
 * We consider that for each coefficient of the matrix, the probability
 * of being non-zero is independent from the others, and given by a
 * probability distribution function which is scale/(i+offset)^alpha.
 *
 * In order to sample this, we approximate as follows.
 * For each row, the expectation of the row weight follows a binomial
 * distribution. We approximate it as a poisson distribution. Then, given
 * the count for the number of non-zero coefficients in the row, and
 * knowing the cumulative distribution function from above, we do inverse
 * transform sampling to find the location of the non-zero coefficients.
 *
 */

struct random_matrix_process_data; /* {{{ */

struct rhs_writer {// {{{
    int n = 0;
    cxx_mpz p;
    std::unique_ptr<FILE, delete_FILE> f;
    operator bool() const { return n; }
    rhs_writer() = default;
    ~rhs_writer() = default;
    rhs_writer(rhs_writer const &) = delete;
    rhs_writer& operator=(rhs_writer const &) = delete;
    rhs_writer(rhs_writer &&) = default;
    rhs_writer& operator=(rhs_writer &&) = default;
    rhs_writer(random_matrix_process_data const &, cxx_param_list &);

    static void declare_usage(cxx_param_list & pl) {
        param_list_decl_usage(pl, "rhs", "rhs output (comma-separated <nrhs>,<prime>,<filename>[,<nullspace_direction>])");
    }
};
// }}}
/* {{{ generic_params_process_loop
 *
 * This reads the full parameter list -- not only the param_list
 * structure --, and fills r with all argument which has been found
 * relevant. This can primarily be seen as a function dedicated to the
 * standalone program, even though the random_matrix= mechanism uses it
 * too as a back-end.
 *
 * This function does *NOT* check that all arguments in pl have been
 * consumed.
 */

static std::vector<int> generic_params_process_loop(cxx_param_list & pl,
        int argc, char const ** argv)
{
    std::vector<int> wild; // nrows ncols coeffs_per_row
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        if (argv[0][0] != '-' && wild.size() < 3) {
            char * tmp;
            wild.push_back((int) strtoul(argv[0], &tmp, 0));
            if (*tmp != '\0') {
                fprintf(stderr, "Parse error for parameter %s\n", argv[0]);
                exit(1);
            }
            argv++, argc--;
            continue;
        }
        throw parameter_error(fmt::format("Unhandled {}", argv[0]));
    }
    return wild;
}

template<typename iterator>
static std::vector<int> generic_params_process_loop(cxx_param_list & pl,
        iterator begin, iterator end)
{
    std::vector<std::string> tmp;
    for(iterator it = begin ; it != end ; ++it)
        tmp.emplace_back(*it);
    std::vector<const char *> argv;
    argv.reserve(tmp.size());
    for(auto const & s : tmp)
        argv.emplace_back(s.c_str());
    return generic_params_process_loop(pl, (int) argv.size(), argv.data());
}
/* }}} */


/* {{{ random_matrix_process_data */
/* This data type gathers the internal state of the random generation */
struct random_matrix_process_data {
    unsigned long nrows = 0;
    unsigned long ncols = 0;
    int density = 0;
    unsigned long seed = 0;
    int maxcoeff = 0;
    int ascii = 0;
    std::unique_ptr<FILE, delete_FILE> owned_out;
    FILE * out = nullptr;

    rhs_writer rhs;

    std::unique_ptr<FILE, delete_FILE> cw, rw;

    static void configure_aliases(cxx_param_list & pl) {
        param_list_configure_alias(pl, "output", "o");
        param_list_configure_alias(pl, "density", "d");
        param_list_configure_alias(pl, "seed", "s");
    }

    static void configure_switches(cxx_param_list & pl) {
        param_list_configure_switch(pl, "binary", nullptr);
        param_list_configure_switch(pl, "freq", nullptr);
    }

    static void process_arguments(cxx_param_list & pl, int argc, char * argv[]);

    static void declare_usage(cxx_param_list & pl) {
        param_list_decl_usage(pl, "density", "desired density per row");
        param_list_decl_usage(pl, "seed", "seed");
        param_list_decl_usage(pl, "c", "add coefficients");
        param_list_decl_usage(pl, "output", "output file name");
        param_list_decl_usage(pl, "binary", "output in binary");
        param_list_decl_usage(pl, "kleft", "ensure at least a left kernel of dimension d");
        param_list_decl_usage(pl, "kright", "ditto for right kernel");
        param_list_decl_usage(pl, "freq", "output row and column weight matrices");
        rhs_writer::declare_usage(pl);
    }

    random_matrix_process_data(
        cxx_param_list & pl,
        std::vector<int> const & wild);

    private:
    struct ctor_helper {
        mutable cxx_param_list pl;
        std::vector<int> wild;
        explicit ctor_helper(std::string const & description)
        {
            auto tokens = split(description, ",");
            wild = generic_params_process_loop(pl, tokens.begin(), tokens.end());
        }
    };

    explicit random_matrix_process_data(ctor_helper const & h)
        : random_matrix_process_data(h.pl, h.wild)
    {}

    public:
    /*
     * This is primarily used for the random_matrix= hack. The standalone
     * program does not follow this route. Here we check that all parts of
     * the provided string are understood as legitimate arguments to
     * random_matrix=
     */
    explicit random_matrix_process_data(const char * str)
        : random_matrix_process_data(ctor_helper(str))
    {
        if (!owned_out)
            out = nullptr;
        /* the default is out == stdout, EXCEPT when we init from a
         * string, where out == NULL is preferred */
    }
};
// }}}

// {{{ rhs_writer::rhs_writer
rhs_writer::rhs_writer(random_matrix_process_data const & R, cxx_param_list & pl)
{
    const char * description = param_list_lookup_string(pl, "rhs");

    if (!description)
        return;

    auto tokens = split(description, ",");
    if (tokens.size() != 3 && tokens.size() != 4)
        throw std::runtime_error("rhs arg invalid");
    std::istringstream(tokens[0]) >> n;
    std::istringstream(tokens[1]) >> p;

    if (n == 0)
        throw std::runtime_error("--rhs argument requires setting more than 0 vectors !");

    const char * rhsname = tokens[2].c_str();
    f.reset(fopen(rhsname, "w"));

    DIE_ERRNO_DIAG(!bool(f), "fopen(%s)", rhsname);

    size_t rhs_rows = R.nrows;
    if (tokens.size() == 4) {
        auto const & l = tokens[3];
        if (l == "left" || l == "LEFT")
            rhs_rows = R.ncols;
        else if (l == "right" || l == "RIGHT")
            rhs_rows = R.nrows;
        else
            throw std::runtime_error("bad nullspace_direction argument in rhs");
    }

    fmt::print(f.get(), "{} {} {}\n", rhs_rows, n, p);
}


random_matrix_process_data::random_matrix_process_data(
        cxx_param_list & pl,
        std::vector<int> const & wild)
{
    /* {{{ parse r->nrows, r->ncols, density */
    if (wild.empty() || (nrows = wild[0]) == 0)
        throw std::runtime_error("Please specify r->nrows");
    if (wild.size() < 2 || (ncols = wild[1]) == 0)
        ncols = nrows;
    if (param_list_parse<int>(pl, "density", density)) {
        if (wild.size() >= 3)
            throw std::runtime_error("density specified twice");
    } else if (wild.size() < 3 || (density == wild[2]) == 0) {
        density = MIN(100, MAX(ncols / 10, MIN(4, ncols)));
    }
    ASSERT_ALWAYS(ncols > 10 && nrows > 10);
    /* }}} */

    param_list_parse(pl, "seed", seed);
    if (!seed) seed = static_cast<unsigned long>(time(nullptr));
    param_list_parse(pl, "c", maxcoeff);

    bool const binary = param_list_parse_switch(pl, "binary");
    bool const freq = param_list_parse_switch(pl, "freq");

    ascii = !binary;

    /* {{{ try to parse the rhs info */
    rhs = rhs_writer(*this, pl);
    ASSERT_ALWAYS(!rhs || maxcoeff > 0);
    /* }}} */

    out = stdout;

    const char * ofilename = nullptr;

    if ((ofilename = param_list_lookup_string(pl, "output"))) {
        owned_out.reset(fopen(ofilename, binary ? "wb" : "w"));
        DIE_ERRNO_DIAG(!bool(owned_out), "fopen(%s)", ofilename);
        out = owned_out.get();
    } else {
        if (binary)
            throw std::runtime_error("--binary requires --output");
        if (freq)
            throw std::runtime_error("--freq requires --output");
    }

    if (freq) {
        std::unique_ptr<char, free_delete<char>> const cwname(
                derived_filename(ofilename, "cw", binary ? ".bin" : ".txt"));
        cw.reset(fopen(cwname.get(), binary ? "wb" : "w"));
        DIE_ERRNO_DIAG(!bool(cw), "fopen(%s)", cwname.get());

        std::unique_ptr<char, free_delete<char>> const rwname(
                derived_filename(ofilename, "rw", binary ? ".bin" : ".txt"));
        rw.reset(fopen(rwname.get(), binary ? "wb" : "w"));
        DIE_ERRNO_DIAG(!bool(rw), "fopen(%s)", rwname.get());
    }
}


/* }}} */
/* }}} */

/* {{{ random_matrix_ddata type -- characteristics of the distribution */
struct random_matrix_ddata : public matrix_column_distribution {
    double alpha = 0;
    double offset = 32; /* this controls the peakedness for the leftmost
                           columns. It is difficult to make this much
                           smaller than 32 presently. Quite unsafe to
                           change. */
    double scale = 1;   /* event function is scale/(x+offset)^alpha */
    unsigned int maxcoeff = 0;      /* 0 for factorization matrices */
    double coeff_alpha = DBL_MAX; /* computed */
    double coeff_n0 = DBL_MAX;    /* computed */
    double mean = DBL_MAX;        /* computed */
    double sdev = DBL_MAX;        /* computed */
    double spread = 1;
    unsigned long ncols = ULONG_MAX;     /* only for constraint correction */
    unsigned long nrows = ULONG_MAX;     /* informational */
    unsigned long padcols = ULONG_MAX;   /* only for constraint correction */
    unsigned long padrows = ULONG_MAX;   /* informational */
    bool print = false;  /* 1 if we should print */

    uint64_t total_coeffs = UINT64_MAX;   /* informational, after generation */
    double row_avg = DBL_MAX;     /* informational, after generation */
    double row_sdev = DBL_MAX;    /* informational, after generation */

    static random_matrix_ddata default_parameters () {
        random_matrix_ddata F;
        F.alpha = 0.94;
        F.offset = 32;
        return F;
    }
    void adjust(random_matrix_process_data const & r, parallelizing_info_srcptr pi, unsigned long padded_nrows, unsigned long padded_ncols);
    void adjust_force_kernel(random_matrix_process_data const & r, parallelizing_info_srcptr pi, unsigned long padded_nrows, unsigned long padded_ncols, int kernel_left, int kernel_right);
    void info(FILE * out) const;

    /* probability mass function */
    double p(double x) const
    {
        return alpha<=0 ? scale : scale*pow(x*spread+offset,-alpha);
    }

    /* cumulative distribution function */
    double q(double x) const override
    {
        if (alpha <= 0) {
            return x * scale;
        }
        double const beta = 1 - alpha;
        double const u = scale / beta / spread;
        return u * (pow(x*spread + offset, beta) - pow(offset, beta));
    }

    /* reciprocal of the cumulative distribution function */
    double qrev(double y) const override
    {
        if (alpha <= 0) {
            return y / scale;
        }
        double const beta = 1 - alpha;
        double const u = scale / beta / spread;
        double const r = pow(y / u + pow(offset, beta), 1 / beta) - offset;
        return r;
    }

    /* variance for the count of successes */
    double qq(double x) const
    {
        if (alpha < 0) {
            /* don't need it */
            abort();
        }
        double const gamma = 1 - 2 * alpha;
        double const v = scale * scale / gamma / spread;
        return v * (pow(x + offset, gamma) - pow(offset, gamma));
    }

    std::vector<uint32_t> generate_row(cxx_gmp_randstate & rstate, punched_interval::pool_t & pool) const;
    std::vector<uint32_t> generate_row(cxx_gmp_randstate & rstate) const;
    uint32_t generate_row(cxx_gmp_randstate & rstate, uint32_t * ptr, punched_interval::pool_t & pool) const;
    int32_t generate_coefficient(cxx_gmp_randstate & rstate, unsigned long j MAYBE_UNUSED) const;

    /* get random matrices, _AND_ fill the stats */
    matrix_u32 get_byrows(cxx_gmp_randstate & rstate);
    matrix_u32 get_bycolumns(cxx_gmp_randstate & rstate);
    matrix_u32 get_u32(parallelizing_info_ptr pi,
            cxx_param_list & pl,
            unsigned long data_nrows, unsigned long data_ncols,
            unsigned long padded_nrows, unsigned long padded_ncols,
            bool transpose) const;
};

/* }}} */

/* the probability mass function for value i is scale/(i+offset)^alpha.
 * the cumulative distribution function (sum on [0,j[ ) is thus:
 * we'll do inverse transform sampling for computing position of non-zero
 * coefficients.
 *
 * scale * ((j+offset)^beta - offset^beta)/beta
 *
 * with beta = 1-alpha
 *
 * the reverse cumulative function is:
 *
 * (x*beta/scale+offset^beta)^(1/beta)-offset
 *
 * the variance is \sum p(i) - \sum p(i)^2  [proof left to reader]
 *
 * and based on this, the last term is:
 *
 * scale^2*((j+offset)^gamma- offset^gamma)/gamma
 *
 * with gamma = 1-2*alpha
 */

/* {{{ more random_matrix_ddata things */

void random_matrix_ddata::info(FILE * out) const
{
    /* some checking and info */
    double const p0 = p(0);
    double const mean0 = double(nrows) * p0;
    double const sdev0 = sqrt(double(nrows) * p0 * (1-p0));
    double const pn = p(double(ncols-1));
    double const mean_n = double(nrows) * pn;
    double const sdev_n = sqrt(double(nrows) * pn * (1-pn));
    fmt::print(out, "Expected row weight: {:.3f}, sdev {:.3f}\n", mean, sdev);
    fmt::print(out, "Expected weight for first column is {:.3f} (sdev {:.3f}, m/sdev=%.1f)\n",
            mean0, sdev0, mean0 / sdev0);
    fmt::print(out, "Expected weight for last column is {:.3f} (sdev {:.3f}, m/sdev=%.1f)\n",
            mean_n, sdev_n, mean_n / sdev_n);
    fmt::print(out, "Worst-case expectation for last column weight by normal approximation: {:.3f}\n",
            extreme_normal(double(nrows), mean_n, -sdev_n));
}

/* in the mmt structures, because of the balancing work, we promised that
 * matrices of size padded_nrows*padded_ncols would be generated on each
 * node. However, we know that the real matrix has to have one particular
 * shape, which means that on the current node, it might be that we'll
 * have some padding rows and columns to generate.
 *
 * Note that the on-the-fly random_matrix setup omits the balancing
 * permutations, so that all padding rows are on the last blocks.
 */
void random_matrix_ddata::adjust_force_kernel(random_matrix_process_data const & R, parallelizing_info_srcptr pi, unsigned long padded_nrows, unsigned long padded_ncols, int kernel_left, int kernel_right)
{
    print = pi ? pi->m->jrank == 0 && pi->m->trank == 0 : true;
    /* Adapt to the parallelizing_info structure : divide */
    /* note that padding has to still be padding. */
    nrows = (R.nrows - kernel_right) / (pi ? pi->wr[1]->totalsize : 1);
    ncols = (R.ncols - kernel_left) / (pi ? pi->wr[0]->totalsize : 1);

#define ADJUST(pi, items, comm, ker) do {				\
    if (pi) {								\
        unsigned int const rk = (comm)->jrank * (comm)->ncores + (comm)->trank;	\
        if (rk * padded_n ## items >= R.n ## items - (ker)) {		\
            n ## items = 0;						\
        } else if ((rk+1) * padded_n ## items >= R.n ## items - (ker)) {	\
            n ## items = R.n ## items - (ker) - rk * padded_n ## items; \
        } else {							\
            n ## items = padded_n ## items;				\
        }								\
        pad ## items = padded_n ## items - n ## items;		\
    } else {								\
        n ## items =  R.n ## items - (ker);				\
    }									\
} while (0)

    ADJUST(pi, rows, pi->wr[1], kernel_right);
    ADJUST(pi, cols, pi->wr[0], kernel_left);

    // experimental: don't scale. Somehow it seems that I'm doing this scaling
    // twice. I shouldn't. Alas, I see no obvious place where this seems to
    // happen.
    //
    // double density = r.density / (pi ? pi->wr[0]->totalsize : 1);
    double const density = R.density;

    /* sets the scale parameter so that the expected row weight matches
     * our desired density target */
    scale = density / q(double(ncols));
    spread = pi ? pi->wr[0]->totalsize : 1;
    mean = q(double(ncols));
    sdev = sqrt(mean * mean - qq(double(ncols)));
    maxcoeff = R.maxcoeff;
    if (maxcoeff) {
        /* Compute n0, which is used to generate coefficients. It
         * essentially counts, in the heaviest column, the number of
         * coefficients equal to 1 */
        double n0 = double(nrows) / 2.0;
        double old_n0 = INFINITY;
        for(int spin = 0 ; n0 != old_n0 && spin < 100  ; spin++) {
            /* How many rows in total would be fit for that n0 ? */
            old_n0 = n0;
            double const alpha = pow(n0, -1.0 / maxcoeff);
            double const y = 2*(n0-1) / (1-alpha);
            n0 /= y / double(nrows);
        }
        coeff_n0 = n0;
        coeff_alpha = pow(n0, -1.0 / maxcoeff);
    }
    if (p(0) >= 1.0) {
        fprintf(stderr, "Error: this density is not acceptable for the current distribution equation. Please adjust the internal offset parameter to something larger.\nrows");
        exit(1);
    }
}

void random_matrix_ddata::adjust(random_matrix_process_data const & r, parallelizing_info_srcptr pi, unsigned long padded_nrows, unsigned long padded_ncols)
{
    adjust_force_kernel(r, pi, padded_nrows, padded_ncols, 0, 0);
}

/* }}} */


std::vector<uint32_t> random_matrix_ddata::generate_row(cxx_gmp_randstate & rstate, punched_interval::pool_t & pool) const
{
    /* pick a row weight */
    /*
       unsigned long weight = random_normal_constrained(rstate, mean, sdev, 0, ncols);
       */
    uint32_t weight;
    std::vector<uint32_t> ret;
    // NOLINTNEXTLINE(bugprone-narrowing-conversions,cppcoreguidelines-narrowing-conversions)
    for( ; (weight = random_poisson(rstate, mean)) >= ncols ; );

    auto range = punched_interval::alloc(pool, 0, mean);
    for(uint32_t i = 0 ; i < weight ; i++) {
        // punched_interval_print(stdout, range);
        uint32_t k = range->pick(pool, *this, rstate);
        if (k >= ncols)
            k = ncols - 1;
        ret.push_back(k);
    }
    std::ranges::sort(ret);
    punched_interval::recycle(std::move(range), pool);
    return ret;
}

std::vector<uint32_t> random_matrix_ddata::generate_row(cxx_gmp_randstate & rstate) const
{
    punched_interval::pool_t pool;
    auto ret = generate_row(rstate, pool);
    return ret;
}

uint32_t random_matrix_ddata::generate_row(cxx_gmp_randstate & rstate, uint32_t * ptr, punched_interval::pool_t & pool) const
{
    auto const v = generate_row(rstate, pool);
    std::ranges::copy(v, ptr);
    return v.size();
}

int32_t random_matrix_ddata::generate_coefficient(cxx_gmp_randstate & rstate, unsigned long j MAYBE_UNUSED) const
{
    unsigned long x = gmp_urandomm_ui(rstate, nrows);
    long neg;
    if ((neg = (x >= nrows/2))) { x -= nrows/2; }

    int c = 1;
    if (j < 100) {
        double const alpha = coeff_alpha;
        auto xd = double(x);
        for(double b = coeff_n0; xd >= b; xd -= b, b *= alpha, c++) ;
    } else {
        c += double(x) < log(coeff_n0);
    }

    if (neg) c = -c;

    return c;
}

#ifndef WANT_MAIN
void random_matrix_fill_fake_balancing_header(balancing & bal, parallelizing_info_ptr pi, const char * rtmp)
{
    random_matrix_process_data const r(rtmp);
    bal.nh = pi->wr[1]->totalsize;
    bal.nv = pi->wr[0]->totalsize;
    bal.nrows = r.nrows;
    bal.ncols = r.ncols;
    bal.ncoeffs = 0; /* FIXME ; what should I do ? */
    bal.checksum = 0;
    bal.flags = FLAG_COLPERM;
    if (bal.nrows == bal.ncols)
        bal.flags |= FLAG_REPLICATE;
    bal.pshuf[0] = 1;
    bal.pshuf[1] = 0;
    bal.pshuf_inv[0] = 1;
    bal.pshuf_inv[1] = 0;
    balancing_set_row_col_count(bal);
}

/*{{{ borrowed from balancing_workhorse.c*/
struct progress_info {
    time_t t;   /* last printed time */
    size_t z;   /* last printed data amount */
};

static int should_print_now(struct progress_info * last_printed, size_t z)
{
    if (z >= last_printed->z + (10UL<<20) || time(nullptr) >= last_printed->t + 10) {
        last_printed->z = z;
        last_printed->t = time(nullptr);
        return 1;
    }
    return 0;
}
/*}}}*/

#if 0
/* This is totally dumb. */
uint32_t * matrix_transpose(uint32_t * p, size_t size, unsigned long nrows, unsigned long ncols)
{
    size_t ncoeffs = size - nrows;
    uint32_t * big = malloc(ncoeffs * 2 * sizeof(uint32_t));
    char buf[16];
    size_disp(ncoeffs * 2 * sizeof(uint32_t), buf);
    printf("allocating temp area of size %s\n", buf);
    uint32_t * q = p;
    uint32_t * w = big;
    uint32_t * fence = big + ncoeffs * 2;
    for(unsigned long i = 0 ; i < nrows ; i++) {
        unsigned long weight = *q++;
        for(unsigned long j = 0 ; j < weight ; j++) {
            /* we'll sort by column indices */
            *w++=*q++;
            *w++=i;
        }
    }
    ASSERT_ALWAYS(w == fence);

    free(p);
    time_t t0 = time(NULL);
    printf("now sorting\n");
    qsort(big, ncoeffs, 2 * sizeof(uint32_t), (sortfunc_t) &cmp_2u32);
    printf("sort time %d s\n", (int) (time(NULL)-t0));
    p = malloc((ncoeffs + ncols) * sizeof(uint32_t));
    q = p;
    w = big;
    for(unsigned long j = 0 ; j < ncols ; j++) {
        uint32_t * qq = q++;
        *qq = 0;
        for( ; w < fence && *w == j ; w += 2) {
            ++*qq;
            *q++ = w[1];
        }
        if (j == ncols-1) ASSERT_ALWAYS(w == fence);
    }
    free(big);
    return p;
}
#endif

matrix_u32 random_matrix_ddata::get_byrows(cxx_gmp_randstate & rstate)
{
    matrix_u32 ret { matrix_u32::withcoeffs_option { maxcoeff > 0 } };

    /* Now we essentially have a copy of printrows above, except that
     * we're outputting binary, and not to a stream but to memory.  */
    total_coeffs = 0;
    double tot_sq = 0;

    time_t const t0 = time(nullptr);
    struct progress_info last_printed[1];
    last_printed->t = t0;
    last_printed->z = 0;

    /* we'd like to avoid constant malloc()'s and free()'s */
    punched_interval::pool_t pool;
    for(unsigned long i = 0 ; i < nrows ; i++) {
        // long v = 0;
        auto v = generate_row(rstate, pool);
        ret.p.push_back(v.size());
        for(auto j : v) {
            ret.p.push_back(j);
            if (ret.withcoeffs)
                ret.p.push_back(generate_coefficient(rstate, j));
        }
        size_t const c = v.size();
        total_coeffs += c;
        tot_sq += (double) c * (double) c;
        if (print && should_print_now(last_printed, ret.p.size() * sizeof(uint32_t))) {
            time_t const dt = last_printed->t - t0;
            char buf[16];
            char buf2[16];
            fmt::print("{}, {} rows in {} s ; {}/s  \n",
                    size_disp(ret.p.size() * sizeof(uint32_t), buf),
                    i, int(dt),
                    size_disp(dt > 0 ? (size_t) (ret.p.size() * sizeof(uint32_t) / dt) : 0, buf2));
            fflush(stdout);
        }
    }
    for(unsigned long j = 0 ; j < padrows ; j++) {
        ret.p.push_back(0);
    }
    if (print) fmt::print("\n");

    double const e = double(total_coeffs) / double(nrows);
    double const s = tot_sq / double(nrows);
    double const sdev = sqrt(s - e*e);
    row_avg = e;
    row_sdev = sdev;
    if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD) && print) {
        printf ("Actual density per row avg %.2f sdev %.2f\n",
                row_avg, row_sdev);
    }

    return ret;
}

matrix_u32 random_matrix_ddata::get_bycolumns(cxx_gmp_randstate & rstate)
{
    matrix_u32 ret { matrix_u32::withcoeffs_option { maxcoeff > 0 } };

    total_coeffs = 0;
    double tot_sq = 0;

    time_t const t0 = time(nullptr);
    struct progress_info last_printed[1];
    last_printed->t = t0;
    last_printed->z = 0;

    /* this will be used as a temporary buffer for the columns being
     * created, before they get pushed to the main matrix (temp area
     * is without coefficients -- those are generated on the second
     * pass).
     */
    std::vector<uint32_t> ptr;

    random_matrix_ddata G;
    /* use a special ddata, for our specially simple process (which
     * still needs the pick-and-punch thing */
    G.alpha = 0;
    G.ncols = nrows; /* yes */

    /* Then in fact it's easier, as we can avoid inverse transform
     * sampling for the computation of the coefficients */
    // int heavy = 1;
    punched_interval::pool_t pool;
    for(unsigned long j = 0 ; j < ncols ; j++) {
        double const p = this->p(double(j));
        G.scale = p;
        ptr.clear();
        unsigned long weight;
        if (p > 0.1) {
            for(unsigned long i = 0 ; i < nrows ; i++) {
                if (random_uniform(rstate) < p)
                    ptr.push_back(i);
            }
        } else {
            weight = (unsigned long) random_binomial(rstate, nrows, p);
            for(unsigned long i = 0 ; i < weight ; i++)
                ptr.push_back(gmp_urandomm_ui(rstate, nrows));
            std::ranges::sort(ptr);
            auto nt = ptr.begin();
            for(auto it = ptr.begin(), jt = it ; it != ptr.end(); it = jt) {
                for(++jt ; jt != ptr.end() && *it == *jt ; ++jt) ;
                *nt++ = *it;
            }
            ptr.erase(nt, ptr.end());
        }
        weight = ptr.size();
        ret.p.push_back(weight);
        for(auto i : ptr) {
            ret.p.push_back(i);
            if (ret.withcoeffs) {
                ret.p.push_back(generate_coefficient(rstate, j));
                // if (r->rhs->n) v += co * (1+ptr[j]);
            }
        }
        total_coeffs += weight;
        tot_sq += (double) weight * (double) weight;
        if (print && should_print_now(last_printed, ret.p.size() * sizeof(uint32_t))) {
            time_t const dt = last_printed->t - t0;
            fmt::print("{}, {} cols in {} s ; {}/s (last weight: {}) \n",
                    size_disp(ret.p.size() * sizeof(uint32_t)),
                    j,
                    int(dt),
                    size_disp(dt > 0 ? (size_t) (ret.p.size() * sizeof(uint32_t) / dt) : 0), weight);
        }
    }
    for(unsigned long j = 0 ; j < padcols ; j++) {
        ret.p.push_back(0);
    }

    double const e = double(total_coeffs) / double(nrows);
    double const s = tot_sq / double(nrows);
    double const sdev = sqrt(s - e*e);
    row_avg = e;
    row_sdev = sdev;
    if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD) && print) {
        printf ("Actual density per row avg %.2f sdev %.2f\n",
                row_avg, row_sdev);
    }

    return ret;
}


matrix_u32 random_matrix_get_u32(parallelizing_info_ptr pi, cxx_param_list & pl, unsigned long data_nrows, unsigned long data_ncols, unsigned long padded_nrows, unsigned long padded_ncols, bool withcoeffs, bool transpose)
{
    const char * rtmp = param_list_lookup_string(pl, "random_matrix");
    ASSERT_ALWAYS(rtmp);
    random_matrix_process_data const r(rtmp);

    /* If not, then the user forgot to add c=<something> in the
     * random_matrix option on the command line.
     */
    ASSERT_ALWAYS(withcoeffs == (r.maxcoeff > 0));

    /* This is not supported here -- mostly because we haven't been
     * extremely careful. */
    ASSERT_ALWAYS(!r.rhs.n);

    random_matrix_ddata F;
    F.adjust(r, pi, data_nrows, data_ncols);

    if (F.print) {
        printf("Each of the %u jobs on %u nodes creates a matrix with %lu rows %lu cols, and %.2f coefficients per row on average. Seed for rank 0 is %lu.\n",
                pi->m->totalsize, pi->m->njobs,
                F.nrows, F.ncols,
                (double) r.density / pi->wr[0]->totalsize, r.seed);
    }

    cxx_gmp_randstate rstate;
    gmp_randseed_ui(rstate, r.seed + pi->m->jrank * pi->m->ncores + pi->m->trank);

    if (transpose) {
        auto mat = F.get_bycolumns(rstate);
        for(unsigned int i = data_ncols ; i < padded_ncols ; i++)
            mat.p.push_back(0);
        return mat;
    } else {
        auto mat = F.get_byrows(rstate);
        for(unsigned int i = data_nrows ; i < padded_nrows ; i++)
            mat.p.push_back(0);
        return mat;
    }
}

#endif

#ifdef  WANT_MAIN

static int avoid_zero_columns = 0;

/* FIXME: this is unholy -- we're actually filling the stats fields in F
 * as a byproduct of reading the matrix.
 */
static void random_matrix_process_print(random_matrix_process_data & r, random_matrix_ddata & F)
{
    int const ascii = r.ascii;
    FILE * out = r.out;
    ASSERT_ALWAYS(out);
    cxx_gmp_randstate rstate;
    gmp_randseed_ui(rstate, r.seed);
    std::unique_ptr<uint32_t[]> const colweights { new uint32_t[r.ncols] };
    memset(colweights.get(), 0, r.ncols * sizeof(uint32_t));
    uint32_t next_priority_col = 0;

    int const has_coeffs = r.maxcoeff > 0;


#define WU32(out, pre, x, post) do {					\
        if (ascii) {							\
            fprintf(out, pre "%" PRIu32 post, (x));			\
        } else {							\
            fwrite(&(x), sizeof(uint32_t), 1, out);			\
        }								\
    } while (0)

#define WS32(out, pre, x, post) do {                                    \
        if (ascii) {							\
            fprintf(out, pre "%" PRId32 post, (x));			\
        } else {							\
            fwrite(&(x), sizeof(uint32_t), 1, out);			\
        }								\
    } while (0)

#define WZa(out, pre, x, post, p) do {					\
            gmp_fprintf(out, pre "%Zd" post, x);			\
    } while (0)


    if (ascii)
        fprintf(out, "%lu %lu\n", r.nrows, r.ncols);
    std::vector<uint32_t> ptr;
    ptr.reserve(r.ncols);
    uint64_t total_coeffs = 0;
    double tot_sq = 0;
    punched_interval::pool_t pool;
    for(unsigned long i = 0 ; i < r.nrows ; i++) {
        long v = 0;
        if (i >= F.nrows)
            ptr.clear();
        else
            ptr = F.generate_row(rstate, pool);
        if (avoid_zero_columns && i >= (unsigned long) (0.9 * double(r.ncols))) {
            for( ; next_priority_col < r.ncols ; next_priority_col++)
                if (!colweights[next_priority_col]) break;
            if (next_priority_col < r.ncols) {
                // don't print anything, because stdout might be our data
                // output...
                // printf("injecting col %" PRIu32 " for row %lu\n", next_priority_col, i);
                ptr.push_back(next_priority_col);
                std::ranges::sort(ptr);
            }
        }
        uint32_t c = ptr.size();
        WU32(out, "", c, "");
        if (r.rw) {
            WU32(r.rw.get(), "", c, "\n");
        }
        for(uint32_t j = 0 ; j < c ; j++) {
            WU32(out, " ", ptr[j], "");
            colweights[ptr[j]]++;
            if (has_coeffs) {
                int32_t co = F.generate_coefficient(rstate, ptr[j]);
                WS32(out, ":", co, "");
                if (r.rhs.n) v += (long) co * (long) (1+ptr[j]);
            }
        }
        if (ascii) { fprintf(out, "\n"); }
        if (r.rhs.n) {
            cxx_mpz x, s;
            mpz_set_si(s, v);

            for(int j = 0 ; j < r.rhs.n - 1 ; j++) {
                mpz_urandomm(x, rstate, r.rhs.p);
                mpz_addmul_ui(s, x, r.ncols + j + 1);
                WZa(r.rhs.f.get(), "", (mpz_srcptr) x, " ", (mpz_srcptr) r.rhs.p);
            }
            mpz_set_si(x, -int(r.ncols + r.rhs.n));
            mpz_invert(x, x, r.rhs.p);
            mpz_mul(s, s, x);
            mpz_mod(s, s, r.rhs.p);
            WZa(r.rhs.f.get(), "", (mpz_srcptr) s, "\n", r.rhs.p);
        }

        total_coeffs += c;
        tot_sq += (double) c * (double) c;
    }
    /* FIXME -- what the hell ? r.nrows is the full length anyway...
    for(unsigned long i = 0 ; i < kernel_right ; i++) {
        if (ascii) {
            fprintf(out, "0\n");
        } else {
            const uint32_t c = 0;
            fwrite(&c, sizeof(uint32_t), 1, out);
        }
    }
    */
    if (r.cw) {
        if (ascii) {
            for(unsigned long j = 0 ; j < r.ncols ; j++) {
                // NOLINTNEXTLINE(bugprone-redundant-branch-condition)
                WU32(r.cw.get(), "", colweights[j], "\n");
            }
        } else {
            fwrite(colweights.get(), sizeof(uint32_t), r.ncols, r.cw.get());
        }
    }
    F.total_coeffs = total_coeffs;
    double const e = (double) total_coeffs / double(r.nrows);
    double const s = (double) tot_sq / double(r.nrows);
    double const sdev = sqrt(s - e*e);
    F.row_avg = e;
    F.row_sdev = sdev;
}

static void usage()
{
    fprintf(stderr, "Usage: ./random_matrix <nrows> [<ncols>] [<density>] [options]\n"
            "Options:\n"
            "\t-d <density> : desired density per row\n"
            "\t-s <seed> : seed\n"
            "\t-c <maxc> : add coefficients\n"
            "\t-v : turn verbosity on\n"
            "\t-Z : avoid zero columns\n"
            "\t-o <matfile> : output file name\n"
            "\t--binary : output in binary\n"
            "\t--kleft <d>: ensure at least a left kernel of dimension d\n"
            "\t--kright <d>: ditto for right kernel\n"
            "\t--rhs <nrhs>,<prime>,<filename>,<nullspace_direction>: rhs output\n"
           );
    exit(1);
}

int main(int argc, char const * argv[])
{
    cxx_param_list pl;
    int verbose = 0;
    int kernel_left = 0;
    int kernel_right = 0;

    random_matrix_process_data::declare_usage(pl);
    random_matrix_process_data::configure_switches(pl);
    random_matrix_process_data::configure_aliases(pl);

    param_list_decl_usage(pl, "v", "turn verbosity on");
    param_list_decl_usage(pl, "Z", "avoid zero columns");
    param_list_configure_switch(pl, "v", nullptr);
    param_list_configure_switch(pl, "Z", nullptr);

    argv++, argc--;
    
    auto wild = generic_params_process_loop(pl, argc, argv);

    verbose = param_list_parse_switch(pl, "v");
    avoid_zero_columns = param_list_parse_switch(pl, "Z");

    random_matrix_process_data r(pl, wild);

    /* {{{ parse kernel size. default is to make the matrix invertible */
    param_list_parse(pl, "kleft", kernel_left);
    param_list_parse(pl, "kright", kernel_right);

    ASSERT_ALWAYS(r.nrows <= INT_MAX);
    ASSERT_ALWAYS(r.ncols <= INT_MAX);

    /* we've been given dimensions for the target matrix, together with a
     * constraint on the kernel size (to be understood as "at least that
     * large").
     *
     * given the row/col unbalance, we may or may not have to generate
     * fewer rows/cols.
     */
    if (r.ncols > r.nrows) {
        if (kernel_right >= int(r.ncols - r.nrows)) {
            kernel_right -= int(r.ncols - r.nrows);
        } else {
            kernel_right = 0;
        }
    }
    if (r.nrows > r.ncols) {
        if (kernel_left >= int(r.nrows - r.ncols)) {
            kernel_left -= int(r.nrows - r.ncols);
        } else {
            kernel_left = 0;
        }
    }
    if (kernel_right > int(r.nrows) / 4) {
        fprintf(stderr, "Warning, right kernel is large."
                " Could trigger misbehaviours\n");
    }
    if (kernel_left > int(r.ncols) / 4) {
        fprintf(stderr, "Warning, left kernel is large."
                " Could trigger misbehaviours\n");
    }
    if (kernel_left >= int(r.ncols)) {
        kernel_left = int(r.ncols - 1);
    }
    if (kernel_right >= int(r.nrows)) {
        kernel_right = int(r.nrows - 1);
    }

    if (r.rhs.n) {
        ASSERT_ALWAYS(kernel_left == 0);
        ASSERT_ALWAYS(kernel_right == 0);
        // kernel_right = 0;
    }
    /* }}} */


    if (param_list_warn_unused(pl)) usage();

    random_matrix_ddata F;
    F.adjust_force_kernel(r, nullptr,
            r.nrows, r.ncols,
            kernel_left, kernel_right);
    if (verbose) F.info(stderr);

    random_matrix_process_print(r, F);

    if (verbose)
        printf ("Actual density per row avg %.2f sdev %.2f\n",
                F.row_avg, F.row_sdev);

    return 0;
}
#endif
