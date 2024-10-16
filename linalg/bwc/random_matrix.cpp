#include "cado.h" // IWYU pragma: keep
                  //
#include <cinttypes>   // PRIu32 // IWYU pragma: keep
#include <climits> // ULONG_MAX
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>              // for uint32_t, int32_t, uint64_t
#include <cmath>
#include <ctime>

#include <gmp.h>

#include "macros.h"
#include "misc.h" // derived_filename mkdir_with_parents next_power_of_2
#include "parallelizing_info.hpp"
#include "portability.h" // strdup // IWYU pragma: keep
#include "params.h"              // for param_list_configure_switch, param_l...
#include "verbose.h"  // verbose_enabled // IWYU pragma: keep
#include "random_distributions.hpp"

/* the files below are not useful for the standalone program */
#ifndef WANT_MAIN
#include "random_matrix.hpp"
#include "balancing.hpp"
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
    std::unique_ptr<FILE> f;
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

std::vector<int> generic_params_process_loop(cxx_param_list & pl,
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
std::vector<int> generic_params_process_loop(cxx_param_list & pl,
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
    std::unique_ptr<FILE> owned_out;
    FILE * out = nullptr;

    rhs_writer rhs;

    std::unique_ptr<FILE> cw, rw;

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
        auto l = tokens[3];
        if (l == "left" || l == "LEFT")
            rhs_rows = R.ncols;
        else if (l == "right" || l == "RIGHT")
            rhs_rows = R.nrows;
        else
            throw std::runtime_error("bad nullspace_direction argument in rhs");

        fmt::print(f.get(), "{} {} {}\n", rhs_rows, n, p);
    }
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
    if (!seed) seed = time(nullptr);
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
        std::unique_ptr<char[]> const cwname(
                derived_filename(ofilename, "cw", binary ? ".bin" : ".txt"));
        cw.reset(fopen(cwname.get(), binary ? "wb" : "w"));
        DIE_ERRNO_DIAG(!bool(cw), "fopen(%s)", cwname.get());

        std::unique_ptr<char[]> const rwname(
                derived_filename(ofilename, "rw", binary ? ".bin" : ".txt"));
        rw.reset(fopen(rwname.get(), binary ? "wb" : "w"));
        DIE_ERRNO_DIAG(!bool(rw), "fopen(%s)", rwname.get());
    }
}


/* }}} */
/* }}} */

/* {{{ random_matrix_ddata type -- characteristics of the distribution */
struct random_matrix_ddata_s {
    double alpha;
    double offset;         /* this controls the peakedness for the leftmost
                           columns. It is difficult to make this much
                           smaller than 32 presently. Quite unsafe to
                           change. */
    double scale;       /* event function is scale/(x+offset)^alpha */
    unsigned int maxcoeff;      /* 0 for factorization matrices */
    double coeff_alpha; /* computed */
    double coeff_n0;    /* computed */
    double mean;        /* computed */
    double sdev;        /* computed */
    double spread;
    unsigned long ncols;        /* only for constraint correction */
    unsigned long nrows;        /* informational */
    unsigned long padcols;        /* only for constraint correction */
    unsigned long padrows;        /* informational */
    int print;  /* 1 if we should print */

    uint64_t total_coeffs;   /* informational, after generation */
    double row_avg;     /* informational, after generation */
    double row_sdev;    /* informational, after generation */
};
typedef struct random_matrix_ddata_s random_matrix_ddata[1];
typedef struct random_matrix_ddata_s * random_matrix_ddata_ptr;
typedef const struct random_matrix_ddata_s * random_matrix_ddata_srcptr;

void random_matrix_ddata_init(random_matrix_ddata_ptr d);
void random_matrix_ddata_set_default(random_matrix_ddata_ptr d);
void random_matrix_ddata_clear(random_matrix_ddata_ptr d);
void random_matrix_ddata_adjust(random_matrix_ddata_ptr f, random_matrix_process_data const & r, parallelizing_info_srcptr pi, unsigned long padded_nrows, unsigned long padded_ncols);
void random_matrix_ddata_adjust_force_kernel(random_matrix_ddata_ptr f, random_matrix_process_data const & r, parallelizing_info_srcptr pi, unsigned long padded_nrows, unsigned long padded_ncols, int kernel_left, int kernel_right);
void random_matrix_ddata_info(FILE * out, random_matrix_ddata_srcptr f);
void random_matrix_ddata_init(random_matrix_ddata_ptr F);
void random_matrix_ddata_clear(random_matrix_ddata_ptr F);
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

/* {{{ direct queries on the random_matrix_ddata type */

/* probability mass function */
double dist_p(random_matrix_ddata_srcptr f, double x)
{
    return f->alpha<=0 ? f->scale : f->scale*pow(x*f->spread+f->offset,-f->alpha);
}

/* cumulative distribution function */
double dist_q(random_matrix_ddata_srcptr f, double x)
{
    if (f->alpha <= 0) {
        return x * f->scale;
    }
    double const beta = 1 - f->alpha;
    double const u = f->scale / beta / f->spread;
    return u * (pow(x*f->spread + f->offset, beta) - pow(f->offset, beta));
}

/* reciprocal of the cumulative distribution function */
double dist_qrev(random_matrix_ddata_ptr f, double y)
{
    if (f->alpha <= 0) {
        return y / f->scale;
    }
    double const beta = 1 - f->alpha;
    double const u = f->scale / beta / f->spread;
    double const r = pow(y / u + pow(f->offset, beta), 1 / beta) - f->offset;
    return r;
}

/* variance for the count of successes */
double dist_qq(random_matrix_ddata_srcptr f, double x)
{
    if (f->alpha < 0) {
        /* don't need it */
        abort();
    }
    double const gamma = 1 - 2 * f->alpha;
    double const v = f->scale * f->scale / gamma / f->spread;
    return v * (pow(x + f->offset, gamma) - pow(f->offset, gamma));
}
/* }}} */

/* {{{ more random_matrix_ddata things */
void random_matrix_ddata_init(random_matrix_ddata_ptr F)
{
    memset(F, 0, sizeof(*F));
    F->scale = 1;
    F->spread = 1;
}
void random_matrix_ddata_clear(random_matrix_ddata_ptr F MAYBE_UNUSED)
{
}

void random_matrix_ddata_set_default(random_matrix_ddata_ptr F)
{
    F->alpha = 0.94;
    F->offset = 32;
    F->scale = 1;
    F->spread = 1;
}

void random_matrix_ddata_info(FILE * out, random_matrix_ddata_srcptr f)
{
    unsigned long const nrows = f->nrows;
    unsigned long const ncols = f->ncols;
    /* some checking and info */
    double const p0 = dist_p(f, 0);
    double const mean0 = nrows * p0;
    double const sdev0 = sqrt(nrows * p0 * (1-p0));
    double const pn = dist_p(f, ncols-1);
    double const mean_n = nrows * pn;
    double const sdev_n = sqrt(nrows * pn * (1-pn));
    fprintf(out, "Expected row weight: %.3f, sdev %.3f\n", f->mean, f->sdev);
    fprintf(out, "Expected weight for first column is %.3f (sdev %.3f, m/sdev=%.1f)\n",
            mean0, sdev0, mean0 / sdev0);
    fprintf(out, "Expected weight for last column is %.3f (sdev %.3f, m/sdev=%.1f)\n",
            mean_n, sdev_n, mean_n / sdev_n);
    fprintf(out, "Worst-case expectation for last column weight by normal approximation: %.3f\n",
            extreme_normal(nrows, mean_n, -sdev_n));
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
void random_matrix_ddata_adjust_force_kernel(random_matrix_ddata_ptr f, random_matrix_process_data const & R, parallelizing_info_srcptr pi, unsigned long padded_nrows, unsigned long padded_ncols, int kernel_left, int kernel_right)
{
    f->print = pi ? pi->m->jrank == 0 && pi->m->trank == 0 : 1;
    /* Adapt to the parallelizing_info structure : divide */
    /* note that padding has to still be padding. */
    f->nrows = (R.nrows - kernel_right) / (pi ? pi->wr[1]->totalsize : 1);
    f->ncols = (R.ncols - kernel_left) / (pi ? pi->wr[0]->totalsize : 1);

#define ADJUST(pi, items, comm, ker) do {				\
    if (pi) {								\
        unsigned int rk = (comm)->jrank * (comm)->ncores + (comm)->trank;	\
        if (rk * padded_n ## items >= R.n ## items - (ker)) {		\
            f->n ## items = 0;						\
        } else if ((rk+1) * padded_n ## items >= R.n ## items - (ker)) {	\
            f->n ## items = R.n ## items - (ker) - rk * padded_n ## items; \
        } else {							\
            f->n ## items = padded_n ## items;				\
        }								\
        f->pad ## items = padded_n ## items - f->n ## items;		\
    } else {								\
        f->n ## items =  R.n ## items - (ker);				\
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
    f->scale = density / dist_q(f, f->ncols);
    f->spread = pi ? pi->wr[0]->totalsize : 1;
    f->mean = dist_q(f, f->ncols);
    f->sdev = sqrt(f->mean * f->mean - dist_qq(f, f->ncols));
    f->maxcoeff = R.maxcoeff;
    if (f->maxcoeff) {
        /* Compute n0, which is used to generate coefficients. It
         * essentially counts, in the heaviest column, the number of
         * coefficients equal to 1 */
        double n0 = f->nrows / 2.0;
        double old_n0 = INFINITY;
        for(int spin = 0 ; n0 != old_n0 && spin < 100  ; spin++) {
            /* How many rows in total would be fit for that n0 ? */
            old_n0 = n0;
            double const alpha = pow(n0, -1.0 / f->maxcoeff);
            double const y = 2*(n0-1) / (1-alpha);
            n0 /= y / f->nrows;
        }
        f->coeff_n0 = n0;
        f->coeff_alpha = pow(n0, -1.0 / f->maxcoeff);
    }
    if (dist_p(f, 0) >= 1.0) {
        fprintf(stderr, "Error: this density is not acceptable for the current distribution equation. Please adjust the internal offset parameter to something larger.\nrows");
        exit(1);
    }
}

void random_matrix_ddata_adjust(random_matrix_ddata_ptr f, random_matrix_process_data const & r, parallelizing_info_srcptr pi, unsigned long padded_nrows, unsigned long padded_ncols)
{
    random_matrix_ddata_adjust_force_kernel(f, r, pi, padded_nrows, padded_ncols, 0, 0);
}

/* }}} */


typedef int (*sortfunc_t)(const void *, const void *);

int cmp_u32(uint32_t const * a, uint32_t const * b)
{
    return (*a > *b) - (*b > *a);
}

uint32_t generate_row(gmp_randstate_t rstate, random_matrix_ddata_ptr f, uint32_t * ptr, punched_interval_ptr range, punched_interval_ptr * pool)
{
    /* pick a row weight */
    /*
       unsigned long weight = random_normal_constrained(rstate, f->mean, f->sdev, 0, f->ncols);
       */
    uint32_t weight;
    for( ; (weight = random_poisson(rstate, f->mean)) >= f->ncols ; );
    // punched_interval_ptr range = punched_interval_alloc(0, f->mean);
    punched_interval_set_full(range, 0, f->mean);
    for(uint32_t i = 0 ; i < weight ; i++) {
        // punched_interval_print(stdout, range);
        uint32_t k = punched_interval_pick(pool, range,
                (double (*)(const void *, double)) dist_q,
                (double (*)(const void *, double)) dist_qrev,
                (const void *) f,
                rstate);
        if (k >= f->ncols)
            k = f->ncols - 1;
        ptr[i] = k;
    }
    qsort(ptr, weight, sizeof(uint32_t), (sortfunc_t) &cmp_u32);
    // punched_interval_free(range);
    return weight;
}

int32_t generate_coefficient(gmp_randstate_t rstate, random_matrix_ddata_ptr F, unsigned long j MAYBE_UNUSED)
{
    unsigned long x = gmp_urandomm_ui(rstate, F->nrows);
    long neg;
    if ((neg = x >= F->nrows/2)) { x -= F->nrows/2; }

    int c = 1;
    if (j < 100) {
        double const alpha = F->coeff_alpha;
        for(double b = F->coeff_n0; x >= b; x -= b, b *= alpha, c++) ;
    } else {
        c += x < log(F->coeff_n0);
    }

    if (neg) c = -c;

    return c;
}

#ifndef WANT_MAIN
void random_matrix_fill_fake_balancing_header(balancing & bal, parallelizing_info_ptr pi, const char * rtmp)
{
    random_matrix_process_data r(rtmp);
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
    if (z >= last_printed->z + (10UL<<20) || time(NULL) >= last_printed->t + 10) {
        last_printed->z = z;
        last_printed->t = time(NULL);
        return 1;
    }
    return 0;
}
/*}}}*/

int cmp_2u32(uint32_t * a, uint32_t * b)
{
    int const r = (*a > *b) - (*b > *a);
    if (r) return r;
    a++;
    b++;
    return (*a > *b) - (*b > *a);
}

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

void random_matrix_get_u32_byrows(gmp_randstate_t rstate, random_matrix_ddata_ptr F, matrix_u32_ptr arg)
{
    int const has_coeffs = F->maxcoeff > 0;

    if (arg->withcoeffs != has_coeffs) {
        fprintf(stderr, "Fatal error, parameters for random matrix generation disagree with the base field.\n");
        if (!has_coeffs)
            fprintf(stderr, "Add the \"c=\" information to the random_matrix argument\n");
        else
            fprintf(stderr, "Remove the \"c=\" information from the random_matrix argument\n");
        exit(EXIT_FAILURE);
    }


    /* Now we essentially have a copy of printrows above, except that
     * we're outputting binary, and not to a stream but to memory.  */
    uint64_t total_coeffs = 0;
    double tot_sq = 0;

    size_t alloc = 0;
    ASSERT_ALWAYS(arg->p == NULL);
    ASSERT_ALWAYS(arg->size == 0);

#define PUSH_P(x) do {    						\
        if (arg->size >= alloc) {					\
            alloc = arg->size + 64 + alloc / 4;			        \
            arg->p = (uint32_t *) realloc(arg->p, alloc * sizeof(uint32_t));	        \
            memset(arg->p + arg->size, 0xFF, (alloc - arg->size) * sizeof(uint32_t)); \
        }								\
        arg->p[arg->size++] = (x);					\
    } while (0)

    time_t const t0 = time(NULL);
    struct progress_info last_printed[1];
    last_printed->t = t0;
    last_printed->z = 0;

        uint32_t * ptr = (uint32_t *) malloc(F->ncols * sizeof(uint32_t));
        /* we'd like to avoid constant malloc()'s and free()'s */
        punched_interval_ptr pool = NULL;
        punched_interval_ptr range = punched_interval_alloc(&pool, 0, 1);
        for(unsigned long i = 0 ; i < F->nrows ; i++) {
            // long v = 0;
            uint32_t const c = generate_row(rstate, F, ptr, range, & pool);
            PUSH_P(c);
            for(unsigned long j = 0 ; j < c ; j++) {
                PUSH_P(ptr[j]);
                if (has_coeffs) {
                    int32_t const co = generate_coefficient(rstate, F, ptr[j]);
                    PUSH_P(co);
                    // if (r->rhs->n) v += co * (1+ptr[j]);
                }
            }
            total_coeffs += c;
            tot_sq += (double) c * (double) c;
            if (F->print && should_print_now(last_printed, arg->size * sizeof(uint32_t))) {
                double const dt = last_printed->t - t0;
                char buf[16];
                char buf2[16];
                printf("%s, %lu rows in %d s ; %s/s  \n",
                        size_disp(arg->size * sizeof(uint32_t), buf), i, (int) dt,
                        size_disp(dt > 0 ? (size_t) (arg->size * sizeof(uint32_t) / dt) : 0, buf2));
                fflush(stdout);
            }
        }
        for(unsigned long j = 0 ; j < F->padrows ; j++) {
            PUSH_P(0);
        }
        if (F->print) printf("\n");
        punched_interval_free(range, &pool);
        punched_interval_free_pool(&pool);
        free(ptr);
#undef PUSH_P

    F->total_coeffs = total_coeffs;
    double const e = (double) total_coeffs / F->nrows;
    double const s = (double) tot_sq / F->nrows;
    double const sdev = sqrt(s - e*e);
    F->row_avg = e;
    F->row_sdev = sdev;
    if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD) && F->print) {
        printf ("Actual density per row avg %.2f sdev %.2f\n",
                F->row_avg, F->row_sdev);
    }
}

void random_matrix_get_u32_bycolumns(gmp_randstate_t rstate, random_matrix_ddata_ptr F, matrix_u32_ptr arg)
{
    uint64_t total_coeffs = 0;
    double tot_sq = 0;

    int const has_coeffs = F->maxcoeff > 0;

    if (arg->withcoeffs != has_coeffs) {
        fprintf(stderr, "Fatal error, parameters for random matrix generation disagree with the base field.\n");
        if (!has_coeffs)
            fprintf(stderr, "Add the \"c=\" information to the random_matrix argument\n");
        else
            fprintf(stderr, "Remove the \"c=\" information from the random_matrix argument\n");
        exit(EXIT_FAILURE);
    }

    size_t alloc = 0;
    ASSERT_ALWAYS(arg->p == NULL);
    ASSERT_ALWAYS(arg->size == 0);

#define PUSH_P(x) do {    						\
        if (arg->size >= alloc) {					\
            alloc = arg->size + 64 + alloc / 4;			        \
            arg->p = (uint32_t *) realloc(arg->p, alloc * sizeof(uint32_t));  \
            memset(arg->p + arg->size, 0xFF, (alloc - arg->size) * sizeof(uint32_t)); \
        }								\
        arg->p[arg->size++] = (x);					\
    } while (0)

    time_t const t0 = time(NULL);
    struct progress_info last_printed[1];
    last_printed->t = t0;
    last_printed->z = 0;

        size_t const size0 = 0;
        /* this will be used as a temporary buffer for the columns being
         * created, before they get pushed to the main matrix (temp area
         * is without coefficients -- those are generated on the second
         * pass).
         */
        uint32_t * ptr = (uint32_t *) malloc(F->nrows * sizeof(uint32_t));
        random_matrix_ddata G;
        random_matrix_ddata_init(G);
        /* use a special ddata, for our specially simple process (which
         * still needs the pick-and-punch thing */
        G->alpha=0;
        G->ncols = F->nrows; /* yes */
        /* Then in fact it's easier, as we can avoid inverse transform
         * sampling for the computation of the coefficients */
        // int heavy = 1;
        punched_interval_ptr pool = NULL;
        for(unsigned long j = 0 ; j < F->ncols ; j++) {
            double const p = dist_p(F, j);
            G->scale = p;
            unsigned long weight;
            if (p > 0.1) {
                weight = 0;
                for(unsigned long i = 0 ; i < F->nrows ; i++) {
                    if (random_uniform(rstate) < p)
                        ptr[weight++]=i;
                }
            } else {
                weight = random_binomial(rstate, F->nrows, p);
                for(unsigned long i = 0 ; i < weight ; i++) {
                    ptr[i] = gmp_urandomm_ui(rstate, F->nrows);
                }
                qsort(ptr, weight, sizeof(uint32_t), (sortfunc_t) &cmp_u32);
                unsigned long nw = 0;
                for(unsigned long i = 0, j ; i < weight ; i=j) {
                    for(j = i + 1; j < weight && ptr[i] == ptr[j] ; j++) ;
                    ptr[nw++] = ptr[i];
                }
                weight = nw;
#if 0
            double wmean = nrows * p; 
            // double wsdev = sqrt(nrows * p * (1-p));
            unsigned long weight = random_binomial(rstate, nrows, p);
            } else if (heavy && weight < sqrt(0.1 * 2 * nrows)) {

            /* pick uniformly a subset of exactly [weight] row
             * indices, within [0..nrows[.  */
                if (pi->m->jrank == 0 && pi->m->trank == 0) {
                    printf("from now on, replacing pick_and_punch by accept-reject\n");
                }
                heavy = 0;
                size0 = arg->size;
                t0 = time(NULL);
            }
            if (heavy) {
                punched_interval_ptr range = punched_interval_alloc(&pool, 0, 1);
                punched_interval_set_full(range, 0, wmean);
                for(unsigned long i = 0 ; i < weight ; i++) {
                    // punched_interval_print(stdout, range);
                    double x = random_uniform(rstate) * (range->b1 - range->holes);
                    unsigned long k = pick_and_punch(G, &pool, range, x);
                    ptr[i] = k;
                }
                punched_interval_free(range, &pool);
                punched_interval_pre_free_pool(&pool, 2 * weight,
                        pi->m->jrank == 0 && pi->m->trank == 0);
                qsort(ptr, weight, sizeof(uint32_t), (sortfunc_t) &cmp_u32);
            } else {
                for(int ok = 0 ; !ok ; ) {
                    for(unsigned long i = 0 ; i < weight ; i++) {
                        ptr[i] = gmp_urandomm_ui(rstate, nrows);
                    }
                    qsort(ptr, weight, sizeof(uint32_t), (sortfunc_t) &cmp_u32);
                    ok=1;
                    for(unsigned long i = 1 ; i < weight ; i++) {
                        if (ptr[i] == ptr[i-1]) {
                            ok=0;
                            break;
                        }
                    }
                }
#endif
            }
            PUSH_P(weight);
            for(unsigned long i = 0 ; i < weight ; i++) {
                PUSH_P(ptr[i]);
                if (has_coeffs) {
                    int32_t const co = generate_coefficient(rstate, F, j);
                    PUSH_P(co);
                    // if (r->rhs->n) v += co * (1+ptr[j]);
                }
            }
            total_coeffs += weight;
            tot_sq += (double) weight * (double) weight;
            if (F->print && should_print_now(last_printed, arg->size * sizeof(uint32_t))) {
                double const dt = last_printed->t - t0;
                char buf[16];
                char buf2[16];
                printf("%s, %lu cols in %d s ; %s/s (last weight: %lu) \n",
                        size_disp(arg->size * sizeof(uint32_t), buf), j, (int) dt,
                        size_disp(dt > 0 ? (size_t) ((arg->size-size0) * sizeof(uint32_t) / dt) : 0, buf2), weight);
                fflush(stdout);
            }
        }
        for(unsigned long j = 0 ; j < F->padcols ; j++) {
            PUSH_P(0);
        }
        punched_interval_free_pool(&pool);
        random_matrix_ddata_clear(G);
        free(ptr);
#undef PUSH_P
    F->total_coeffs = total_coeffs;
    double const e = (double) total_coeffs / F->nrows;
    double const s = (double) tot_sq / F->nrows;
    double const sdev = sqrt(s - e*e);
    F->row_avg = e;
    F->row_sdev = sdev;
    if (verbose_enabled(CADO_VERBOSE_PRINT_BWC_CACHE_BUILD) && F->print) {
        printf ("Actual density per row avg %.2f sdev %.2f\n",
                F->row_avg, F->row_sdev);
    }
}


void random_matrix_get_u32(parallelizing_info_ptr pi, param_list pl, matrix_u32_ptr arg, unsigned long data_nrows, unsigned long data_ncols, unsigned long padded_nrows, unsigned long padded_ncols)
{
    const char * rtmp = param_list_lookup_string(pl, "random_matrix");
    ASSERT_ALWAYS(rtmp);
    random_matrix_process_data const r(rtmp);

    /* This is not supported here -- mostly because we haven't been
     * extremely careful. */
    ASSERT_ALWAYS(!r.rhs.n);

    random_matrix_ddata F;
    random_matrix_ddata_init(F);
    random_matrix_ddata_set_default(F);
    random_matrix_ddata_adjust(F, r, pi, data_nrows, data_ncols);

    if (F->print) {
        printf("Each of the %u jobs on %u nodes creates a matrix with %lu rows %lu cols, and %.2f coefficients per row on average. Seed for rank 0 is %lu.\n",
                pi->m->totalsize, pi->m->njobs,
                F->nrows, F->ncols, (double) r.density / pi->wr[0]->totalsize, r.seed);
    }

    cxx_gmp_randstate rstate;
    gmp_randseed_ui(rstate, r.seed + pi->m->jrank * pi->m->ncores + pi->m->trank);

#define PUSH_P(x) do {    						\
        if (arg->size >= alloc) {					\
            alloc = arg->size + 64 + alloc / 4;			        \
            arg->p = (uint32_t *) realloc(arg->p, alloc * sizeof(uint32_t));	        \
            memset(arg->p + arg->size, 0xFF, (alloc - arg->size) * sizeof(uint32_t)); \
        }								\
        arg->p[arg->size++] = (x);					\
    } while (0)
    /* This is ugly, we should store alloc within arg. But this whole
     * embarrassment of a type is meant to go away someday anyway and I
     * have a branch that kills it, so let's touch only the minimum
     */
    size_t alloc = arg->size;
    if (arg->transpose) {
        random_matrix_get_u32_bycolumns(rstate, F, arg);
        for(unsigned int i = data_ncols ; i < padded_ncols ; i++) {
            PUSH_P(0);
        }
    } else {
        random_matrix_get_u32_byrows(rstate, F, arg);
        for(unsigned int i = data_nrows ; i < padded_nrows ; i++) {
            PUSH_P(0);
        }
    }
#undef PUSH_P

    random_matrix_ddata_clear(F);
}

#endif

#ifdef  WANT_MAIN

int avoid_zero_columns = 0;

void random_matrix_process_print(random_matrix_process_data & r, random_matrix_ddata_ptr F)
{
    int const ascii = r.ascii;
    FILE * out = r.out;
    ASSERT_ALWAYS(out);
    cxx_gmp_randstate rstate;
    gmp_randseed_ui(rstate, r.seed);
    std::unique_ptr<uint32_t[]> colweights { new uint32_t[r.ncols] };
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
    std::unique_ptr<uint32_t[]> ptr { new uint32_t[r.ncols] };
    uint64_t total_coeffs = 0;
    double tot_sq = 0;
    punched_interval_ptr pool = NULL;
    punched_interval_ptr range = punched_interval_alloc(&pool, 0, 1);
    for(unsigned long i = 0 ; i < r.nrows ; i++) {
        long v = 0;
        uint32_t c;
        if (i >= F->nrows)
            c = 0;
        else
            c = generate_row(rstate, F, ptr.get(), range, &pool);
        if (avoid_zero_columns && i >= 0.9 * r.ncols) {
            for( ; next_priority_col < r.ncols ; next_priority_col++)
                if (!colweights[next_priority_col]) break;
            if (next_priority_col < r.ncols) {
                // don't print anything, because stdout might be our data
                // output...
                // printf("injecting col %" PRIu32 " for row %lu\n", next_priority_col, i);
                ptr[c++] = next_priority_col;
                qsort(ptr.get(), c, sizeof(uint32_t), (sortfunc_t) &cmp_u32);
            }
        }
        WU32(out, "", c, "");
        if (r.rw) {
            WU32(r.rw.get(), "", c, "\n");
        }
        for(uint32_t j = 0 ; j < c ; j++) {
            WU32(out, " ", ptr[j], "");
            colweights[ptr[j]]++;
            if (has_coeffs) {
                int32_t co = generate_coefficient(rstate, F, ptr[j]);
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
            mpz_set_si(x, -(r.ncols + r.rhs.n));
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
                WU32(r.cw.get(), "", colweights[j], "\n");
            }
        } else {
            fwrite(colweights.get(), sizeof(uint32_t), r.ncols, r.cw.get());
        }
    }
    punched_interval_free(range, &pool);
    punched_interval_free_pool(&pool);
    F->total_coeffs = total_coeffs;
    double const e = (double) total_coeffs / r.nrows;
    double const s = (double) tot_sq / r.nrows;
    double const sdev = sqrt(s - e*e);
    F->row_avg = e;
    F->row_sdev = sdev;
}

void usage()
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
    unsigned long kernel_left = 0;
    unsigned long kernel_right = 0;

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
    param_list_parse_ulong(pl, "kleft", &kernel_left);
    param_list_parse_ulong(pl, "kright", &kernel_right);

    /* we've been given dimensions for the target matrix, together with a
     * constraint on the kernel size (to be understood as "at least that
     * large").
     *
     * given the row/col unbalance, we may or may not have to generate
     * fewer rows/cols.
     */
    if (r.ncols > r.nrows) {
        if (kernel_right >= r.ncols - r.nrows) {
            kernel_right -= r.ncols - r.nrows;
        } else {
            kernel_right = 0;
        }
    }
    if (r.nrows > r.ncols) {
        if (kernel_left >= r.nrows - r.ncols) {
            kernel_left -= r.nrows - r.ncols;
        } else {
            kernel_left = 0;
        }
    }
    if (kernel_right > r.nrows / 4) {
        fprintf(stderr, "Warning, right kernel is large."
                " Could trigger misbehaviours\n");
    }
    if (kernel_left > r.ncols / 4) {
        fprintf(stderr, "Warning, left kernel is large."
                " Could trigger misbehaviours\n");
    }
    if (kernel_left >= r.ncols) {
        kernel_left = r.ncols - 1;
    }
    if (kernel_right >= r.nrows) {
        kernel_right = r.nrows - 1;
    }

    if (r.rhs.n) {
        ASSERT_ALWAYS(kernel_left == 0);
        ASSERT_ALWAYS(kernel_right == 0);
        // kernel_right = 0;
    }
    /* }}} */


    if (param_list_warn_unused(pl)) usage();

    random_matrix_ddata F;
    random_matrix_ddata_init(F);
    random_matrix_ddata_set_default(F);
    random_matrix_ddata_adjust_force_kernel(F, r, NULL, r.nrows, r.ncols, kernel_left, kernel_right);
    if (verbose) random_matrix_ddata_info(stderr, F);

    random_matrix_process_print(r, F);

    if (verbose)
        printf ("Actual density per row avg %.2f sdev %.2f\n",
                F->row_avg, F->row_sdev);
    random_matrix_ddata_clear(F);

    return 0;
}
#endif
