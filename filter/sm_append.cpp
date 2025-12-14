/* Schirokauer maps

   This is largely copied from sm_simple.cpp ; however the purpose of
   sm_simple was to be kept simple and stupid, so let's keep it like
   this.

   This program is to be used as an mpi accelerator for reconstructlog-dl

   Given a relation file where each line begins with an a,b pair IN
   HEXADECIMAL FORMAT, output the same line, appending the SM values in
   the end.

   SM computation is offloaded to the (many) MPI jobs.


   alternatively, this program can run multithreaded (and not MPI) with -t.

   The code is not nice.
*/

#include "cado.h" // IWYU pragma: keep

#include <cinttypes>
#include <cstdarg>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <condition_variable>
#include <functional>
#include <iterator>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "fmt/format.h"
#include <gmp.h>

#include "cado_poly.h"
#include "gzip.h"
#include "macros.h"
#include "mpz_poly.h"
#include "params.h"
#include "runtime_numeric_cast.hpp"
#include "select_mpi.h"
#include "sm_utils.hpp"
#include "timing.h"
#include "verbose.h"

// {{{ debug print interface
static unsigned int const debug = 0;

#define CSI_RED "\033[00;31m"
#define CSI_GREEN "\033[00;32m"
#define CSI_YELLOW "\033[00;33m"
#define CSI_BLUE "\033[00;34m"
#define CSI_PINK "\033[00;35m"
#define CSI_BOLDRED "\033[01;31m"
#define CSI_BOLDGREEN "\033[01;32m"
#define CSI_BOLDYELLOW "\033[01;33m"
#define CSI_BOLDBLUE "\033[01;34m"
#define CSI_BOLDPINK "\033[01;35m"
#define CSI_RESET "\033[m"

static void sm_append_log(std::string const & s)
{
    if (!debug)
        return;
    fmt::print("{:.3f} {}\n", wct_seconds(), s);
}

// }}}
struct print_wrap { // {{{
    std::string context;
    std::string what;
    std::string csi;
    double t0;
    explicit print_wrap(std::string const & context, std::string const & what,
                        std::string const & csi = CSI_RESET)
        : context(context)
        , what(what)
        , csi(csi)
        , t0(wct_seconds())
    {
        sm_append_log(
            fmt::format("{}{} start {}" CSI_RESET, csi, context, what));
    }
    print_wrap(print_wrap const &) = delete;
    print_wrap(print_wrap &&) = delete;
    print_wrap & operator=(print_wrap const &) = delete;
    print_wrap & operator=(print_wrap &&) = delete;
    ~print_wrap()
    {
        sm_append_log(fmt::format("{}{} done {}" CSI_RESET " [taken {:.1f}]",
                                  csi, context, what, wct_seconds() - t0));
    }
};
// }}}

struct ab_pair {   // {{{
    int64_t a = 0; /* only a is allowed to be negative */
    uint64_t b = 0;
    explicit ab_pair(char const * p);
    ab_pair() = default;
};

ab_pair::ab_pair(char const * p)
{
    char const * p0 = p;
    int64_t sign = 1;
    if (*p == '-') {
        sign = -1;
        p++;
    }
    if (sscanf(p, "%" SCNx64 ",%" SCNx64 ":", &a, &b) < 2) {
        fprintf(stderr, "Parse error at line: %s\n", p0);
        exit(EXIT_FAILURE);
    }
    a *= sign;
}

// }}}

struct sm_capable { // {{{
    std::vector<sm_side_info> const & sm_info;
    int nsm_total = 0;
    size_t limbs_per_ell = 0;
    explicit sm_capable(std::vector<sm_side_info> const & sm_info)
        : sm_info(sm_info)
    {
        for (auto const & S: sm_info) {
            nsm_total += S.nsm;
            if (S.nsm)
                limbs_per_ell = mpz_size(S.ell);
        }
    }

    size_t compute_datasize(std::vector<ab_pair> const & batch) const
    {
        return batch.size() * nsm_total * limbs_per_ell;
    }
    std::unique_ptr<mp_limb_t[]>
    compute(std::vector<ab_pair> const & batch) const
    {
        std::unique_ptr<mp_limb_t[]> res(
            new mp_limb_t[compute_datasize(batch)]);
        std::fill_n(res.get(), compute_datasize(batch), 0);

        {
            cxx_mpz_poly smpol, pol;
            mp_limb_t * dst = res.get();
            for (auto const & B: batch) {
                mpz_poly_setcoeff_int64(pol, 0, B.a);
                mpz_poly_setcoeff_int64(pol, 1, -(int64_t)B.b);
                int smidx = 0;
                for (auto const & S: sm_info) {
                    S.compute_piecewise(smpol, pol);
                    for (int k = 0; k < S.nsm; k++, smidx++) {
                        if (k <= smpol->deg) {
                            for (int j = 0; j < int(limbs_per_ell); j++) {
                                dst[j] = mpz_getlimbn(
                                    mpz_poly_coeff_const(smpol, k), j);
                            }
                        }
                        dst += limbs_per_ell;
                    }
                }
            }
        }
        return res;
    }
}; // }}}

static unsigned int batch_size = 128;

struct task_globals : public sm_capable { // {{{
    FILE * in;
    FILE * out;
    size_t nrels_in = 0;
    size_t nrels_out = 0;
    int npeers = 0;
    task_globals(FILE * in, FILE * out,
                 std::vector<sm_side_info> const & sm_info)
        : sm_capable(sm_info)
        , in(in)
        , out(out)
    {
    }
};
// }}}

struct peer_link_base { // {{{
    int id;
    int nb_tasks_posted = 0;
    int nb_tasks_done = 0;
    explicit peer_link_base(int id)
        : id(id)
    {
    }
    std::string name() const { return fmt::format("peer {}", id); }
};
// }}}

struct downlink_base : public peer_link_base { // {{{
    std::vector<ab_pair> batch;
    std::vector<std::string> rels;
    std::unique_ptr<mp_limb_t[]> returns;
    explicit downlink_base(int id)
        : peer_link_base(id)
    {
    }

    void consume(task_globals & tg);
    int produce(task_globals & tg);
};
// }}}

struct thread_link_common : public downlink_base { // {{{
    std::mutex m;
    std::condition_variable work_to_do;
    std::condition_variable work_done;

    explicit thread_link_common(int id)
        : downlink_base(id)
    {
    }
};
// }}}

/* {{{ uplinks and peer control loops */
/* This is the worker-facing struct, thread flavor and then mpi flavor */
struct uplink_thread { // {{{
    thread_link_common & C;
    int id() const { return C.id; }
    std::vector<ab_pair> const & get_input(int)
    {
        std::unique_lock<std::mutex> lk(C.m);
        for (; C.nb_tasks_posted == C.nb_tasks_done;)
            C.work_to_do.wait(lk);
        return C.batch;
    }
    void put_output(int, std::unique_ptr<mp_limb_t[]> && returns, size_t)
    {
        C.returns = std::move(returns);
        std::unique_lock<std::mutex> const lk(C.m);
        C.nb_tasks_done++;
        C.work_done.notify_one();
    }
    explicit uplink_thread(thread_link_common & C)
        : C(C)
    {
    }
};
// }}}
struct uplink_mpi { // {{{
    int rank = -1;
    ;
    int id() const { return rank; }

    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    std::vector<ab_pair> get_input(int turn) const
    {
        unsigned long bsize;
        MPI_Recv(&bsize, 1, MPI_UNSIGNED_LONG, 0, turn, MPI_COMM_WORLD,
                 MPI_STATUS_IGNORE);
        if (bsize == 0)
            return {};
        std::vector<ab_pair> batch(bsize, ab_pair());
        MPI_Recv((char *)batch.data(),
                 runtime_numeric_cast<int>(bsize * sizeof(ab_pair)), MPI_BYTE,
                 0, turn, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        return batch;
    }

    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    void put_output(int turn, std::unique_ptr<mp_limb_t[]> && returns,
                    size_t sz)
    {
        MPI_Send(returns.get(),
                 runtime_numeric_cast<int>(sz * sizeof(mp_limb_t)), MPI_BYTE, 0,
                 turn, MPI_COMM_WORLD);
    }
    uplink_mpi()
    {
        /* We could possibly construct with a peer_link_base, but
         * then we'd have to keep track of the data fields, do we really
         * want that?
         */
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }
};
// }}}

/* The worker can be implemented generically */
template <typename uplink_base> // {{{
struct uplink
    : public uplink_base
    , public sm_capable {
    template <typename... Args>
    explicit uplink(std::vector<sm_side_info> const & sm_info, Args &&... args)
        : uplink_base(std::forward<Args>(args)...)
        , sm_capable(sm_info)
    {
    }

    using uplink_base::id;

    void loop()
    {
        if (id() == 0)
            return;

        for (int turn = 0;; turn++) {
            auto const me = fmt::format("turn {} peer {}", turn, id());

            auto batch = uplink_base::get_input(turn);
            if (batch.empty()) {
                sm_append_log(fmt::format("{} receive finish", me));
                break;
            }

            double t0 = wct_seconds();
            std::unique_ptr<mp_limb_t[]> returns;
            {
                print_wrap const pw(
                    me, fmt::format("batch of size {}", batch.size()),
                    CSI_BLUE);
                returns = compute(batch);
            }

            fprintf(stderr,
                    "# peer processes batch of %zu in %.1f [%.1f SMs/s]\n",
                    batch.size(), wct_seconds() - t0,
                    double(batch.size()) / (wct_seconds() - t0));

            {
                print_wrap const pw(me, "send return");
                uplink_base::put_output(turn, std::move(returns),
                                        compute_datasize(batch));
            }
        }
        sm_append_log(fmt::format(CSI_RED "peer {} leaving" CSI_RESET, id()));
    }
}; // }}}

static void uplink_loop_thread(std::vector<sm_side_info> const & sm_info,
                               thread_link_common & C)
{
    uplink<uplink_thread> U(sm_info, C);
    U.loop();
}

static void uplink_loop_mpi(std::vector<sm_side_info> const & sm_info)
{
    uplink<uplink_mpi> U(sm_info);
    U.loop();
}
/* }}} */

/* {{{ downlink specifics */
struct downlink_thread : public thread_link_common { // {{{
    // NOLINTNEXTLINE(readability-convert-member-functions-to-static)
    void put_input(int)
    {
        std::unique_lock<std::mutex> const lk(m);
        nb_tasks_posted++;
        work_to_do.notify_one();
    }

    void get_output(int, size_t)
    {
        std::unique_lock<std::mutex> lk(m);
        for (; nb_tasks_posted > nb_tasks_done;)
            work_done.wait(lk);
        batch.clear();
    }

    std::thread worker;

    downlink_thread(downlink_thread const &) = delete;
    downlink_thread(downlink_thread &&) noexcept = delete;
    downlink_thread & operator=(downlink_thread const &) = delete;
    downlink_thread & operator=(downlink_thread &&) noexcept = delete;
    downlink_thread() = delete;
    ~downlink_thread() { worker.join(); }

    downlink_thread(int id, std::vector<sm_side_info> const & sm_info)
        : thread_link_common(id)
        , worker(uplink_loop_thread, std::cref(sm_info),
                 std::ref((thread_link_common &)*this))
    {
    }

    explicit downlink_thread(
        std::pair<int, std::vector<sm_side_info> const &> p)
        : downlink_thread(p.first, p.second)
    {
    }
};
// }}}
struct downlink_mpi : public downlink_base { // {{{
    void put_input(int turn) const
    {
        unsigned long const bsize = batch.size();
        MPI_Send(&bsize, 1, MPI_UNSIGNED_LONG, id, turn, MPI_COMM_WORLD);
        if (batch.empty())
            return;
        MPI_Send((char *)batch.data(),
                 runtime_numeric_cast<int>(batch.size() * sizeof(ab_pair)),
                 MPI_BYTE, id, turn, MPI_COMM_WORLD);
    }

    void get_output(int turn, size_t sz)
    {
        returns.reset(new mp_limb_t[sz]);

        MPI_Recv(returns.get(),
                 runtime_numeric_cast<int>(sz * sizeof(mp_limb_t)), MPI_BYTE,
                 id, turn, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    explicit downlink_mpi(int id, std::vector<sm_side_info> const &)
        : downlink_base(id)
    {
    }
    explicit downlink_mpi(std::pair<int, std::vector<sm_side_info> const &> p)
        : downlink_mpi(p.first, p.second)
    {
    }
};
// }}}
/* }}} */

/* {{{ downlink_base generics */
void downlink_base::consume(task_globals & tg) // {{{
{
    mp_limb_t const * src = returns.get();
    for (auto const & r: rels) {
        fputs(r.c_str(), tg.out);
        bool comma = false;
        for (int j = 0; j < tg.nsm_total; j++) {
            gmp_fprintf(tg.out, "%c%Nd", comma ? ',' : ':', src,
                        tg.limbs_per_ell);
            src += tg.limbs_per_ell;
            comma = true;
        }
        fputc('\n', tg.out);
        tg.nrels_out++;
    }
    rels.clear();
    batch.clear();
}
// }}}
int downlink_base::produce(task_globals & tg) // {{{
{
    int eof = 0;
    char buf[1024];

    ASSERT_ALWAYS(batch.empty());
    ASSERT_ALWAYS(rels.empty());

    while (!eof && batch.size() < batch_size && fgets(buf, 1024, tg.in)) {
        int n = strlen(buf);
        if (!n) {
            fprintf(stderr, "Got 0-sized buffer in fgets, shouldn't happen. "
                            "Assuming EOF.\n");
            eof = true;
            break;
        }
        buf[n - 1] = '\0';

        if (buf[0] == '#') {
            fputs(buf, tg.out);
            fputc('\n', tg.out);
            continue;
        }

        tg.nrels_in++;
        rels.emplace_back(buf);

        batch.emplace_back(buf);
    }
    if (!eof && batch.size() < batch_size) {
        eof = true;
        if (ferror(stdin)) {
            fprintf(stderr, "Error on stdin\n");
        }
    }
    return eof;
}
// }}}
/* }}} */

/* {{{ mrepl2 -- helper struct
 * This is an iterator that can be dereferenced exactly n times, and
 * returns the same reference every time, together with a counter.
 * We use it in order to create a vector will all instances created with
 * the same ctor args
 */
template <typename T> struct repl2 {
    using iterator_category = std::forward_iterator_tag;
    using value_type = std::pair<int, T const &>;
    using difference_type = ptrdiff_t;
    using pointer = value_type *;
    using reference = value_type &;
    T & r;
    int i = 0;
    int n;
    repl2(T & r, int n)
        : r(r)
        , n(n)
    {
    }
    value_type operator*() { return {i, r}; }
    repl2<T> & operator++()
    {
        i++;
        return *this;
    }
    bool operator==(repl2<T> const & o) const { return (n - i) == (o.n - o.i); }
    bool operator!=(repl2<T> const & o) const { return !operator==(o); }
};
template <typename T> static repl2<T> mrepl2(T & tg, int n = 0)
{
    return {tg, n};
}

/* }}} */

/* {{{ main control loop */
template <typename T>
static void sm_append_master(FILE * in, FILE * out,
                             std::vector<sm_side_info> const & sm_info,
                             int size)
{
    /* need to know how many mp_limb_t's we'll get back from each batch */
    task_globals tg(in, out, sm_info);

    std::vector<T> peers(mrepl2(sm_info, size), mrepl2(sm_info));

    int eof = 0;
    /* eof = 1 on first time. eof = 2 when all receives are done */

    fprintf(stderr, "# running master with %d slaves, batch size %u\n",
            size - 1, batch_size);

    fprintf(stderr, "# if running under MPI, make sure you use \"--bind-to "
                    "core\" or equivalent\n");

    double const t0 = wct_seconds();
    int turn;
    for (turn = 0; eof <= 2; turn++, eof += !!eof) {
        auto const me = fmt::format("turn {}", turn);
        print_wrap const tw(me, "", CSI_BOLDGREEN);
        for (auto & peer: peers) {
            if (peer.id == 0)
                continue;

            if (eof && (turn && peer.batch.empty()))
                /* Our last send was a 0-send, so we have nothing to do.
                 * If turn == 0, we haven't yed had a chance to send a
                 * batch of size 0, so we still need to do something,
                 * though.
                 */
                continue;

            if (turn) {
                print_wrap const pw(
                    fmt::format("turn {}", turn - 1),
                    fmt::format("receive from {}", peer.name()));
                peer.get_output(turn - 1, tg.compute_datasize(peer.batch));
                peer.consume(tg);
            }

            if (eof) {
                print_wrap const pw(
                    me, fmt::format("send finish to {}", peer.name()));
                peer.put_input(turn);
            } else if (!eof) {
                print_wrap const pw(me, fmt::format("send to {}", peer.name()));
                eof = peer.produce(tg);
                peer.put_input(turn);
            }
        }
        if (turn && !(turn & (turn + 1))) {
            /* print only when turn is a power of two */
            fprintf(stderr,
                    "# printed %zu rels in %.1f s"
                    " (%.1f / batch, %.1f rels/s)\n",
                    tg.nrels_out, wct_seconds() - t0,
                    (wct_seconds() - t0) / turn,
                    double(tg.nrels_out) / (wct_seconds() - t0));
        }
    }
    fprintf(stderr,
            "# final: printed %zu rels in %.1f s"
            " (%.1f / batch, %.1f rels/s)\n",
            tg.nrels_out, wct_seconds() - t0, (wct_seconds() - t0) / turn,
            double(tg.nrels_out) / (wct_seconds() - t0));
}
/* }}} */

/* {{{ fallback code */
static void sm_append_sync(FILE * in, FILE * out,
                           std::vector<sm_side_info> const & sm_info)
{
    char buf[1024];
    cxx_mpz_poly smpol;
    int maxdeg = 0;
    for (auto const & S: sm_info)
        maxdeg = std::max(maxdeg, S.f->deg);

    while (fgets(buf, 1024, in)) {
        int n = strlen(buf);
        if (!n)
            break;
        buf[n - 1] = '\0';

        if (buf[0] == '#') {
            fputs(buf, out);
            fputc('\n', out);
            continue;
        }

        char * p = buf;
        int64_t a; /* only a is allowed to be negative */
        uint64_t b;
        int64_t sign = 1;
        if (*p == '-') {
            sign = -1;
            p++;
        }
        if (sscanf(p, "%" SCNx64 ",%" SCNx64 ":", &a, &b) < 2) {
            fprintf(stderr, "Parse error at line: %s\n", buf);
            exit(EXIT_FAILURE);
        }

        cxx_mpz_poly pol;

        mpz_poly_set_ab(pol, a * sign, b);

        fputs(buf, out);
        fputc(':', out);
        bool has_sep = true;
        for (auto const & S: sm_info) {
            S.compute_piecewise(smpol, pol);
            if (!has_sep)
                fputc(',', out);
            print_sm2(out, S, smpol, ",");
            has_sep = S.nsm == 0;
        }
        fputc('\n', out);
    }
}
/* }}} */

/* {{{ multiplexer */
static void sm_append(FILE * in, FILE * out,
                      std::vector<sm_side_info> const & sm_info, int nthr)
{
    int rank;
    int size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size > 1) {
        if (nthr > 1)
            throw std::runtime_error("This program won't do MPI **AND** "
                                     "threads at the same time!\n");
        if (rank == 0) {
            sm_append_master<downlink_mpi>(in, out, sm_info, size * nthr);
        } else {
            uplink_loop_mpi(sm_info);
        }
    } else if (nthr > 1) {
        /* Because we mimick what goes on in the MPI setting, peer 0 will
         * do nothing, actually */
        sm_append_master<downlink_thread>(in, out, sm_info, size * nthr);
    } else {
        sm_append_sync(in, out, sm_info);
    }
}
/* }}} */

static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "poly", "(required) poly file");
    param_list_decl_usage(pl, "ell", "(required) group order");
    param_list_decl_usage(pl, "nsm", "number of SMs to use per side");
    param_list_decl_usage(pl, "sm-mode", "SM mode (see sm-portability.h)");
    param_list_decl_usage(pl, "in", "data input (defaults to stdin)");
    param_list_decl_usage(pl, "out", "data output (defaults to stdout)");
    param_list_decl_usage(pl, "b", "batch size for loop");
    param_list_decl_usage(pl, "t", "number of threads (exclusive with MPI!)");
    verbose_decl_usage(pl);
}

static void usage(char const * argv, char const * missing, param_list pl)
{
    if (missing) {
        fprintf(stderr, "\nError: missing or invalid parameter \"-%s\"\n",
                missing);
    }
    param_list_print_usage(pl, argv, stderr);
    exit(EXIT_FAILURE);
}

/* -------------------------------------------------------------------------- */

// coverity[root_function]
int main(int argc, char const * argv[])
{
    MPI_Init(&argc, (char ***)&argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    char const * argv0 = argv[0];

    char const * polyfile = nullptr;

    param_list pl;
    cado_poly cpoly;

    mpz_t ell;

    /* read params */
    param_list_init(pl);
    declare_usage(pl);

    if (argc == 1)
        usage(argv[0], nullptr, pl);

    argc--, argv++;
    for (; argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) {
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        usage(argv0, nullptr, pl);
    }

    /* Read poly filename from command line */
    if ((polyfile = param_list_lookup_string(pl, "poly")) == nullptr) {
        fprintf(stderr, "Error: parameter -poly is mandatory\n");
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    /* Read ell from command line (assuming radix 10) */
    mpz_init(ell);
    if (!param_list_parse_mpz(pl, "ell", ell)) {
        fprintf(stderr, "Error: parameter -ell is mandatory\n");
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    /* Init polynomial */
    cado_poly_init(cpoly);
    cado_poly_read(cpoly, polyfile);

    std::vector<mpz_poly_srcptr> F(cpoly->nb_polys, nullptr);

    for (int side = 0; side < cpoly->nb_polys; side++)
        F[side] = cpoly->pols[side];

    std::vector<int> nsm_arg(cpoly->nb_polys, -1);
    param_list_parse_int_args_per_side(pl, "nsm", nsm_arg.data(),
                                       cpoly->nb_polys,
                                       ARGS_PER_SIDE_DEFAULT_AS_IS);

    FILE * in = rank ? nullptr : stdin;
    FILE * out = rank ? nullptr : stdout;
    char const * infilename = param_list_lookup_string(pl, "in");
    char const * outfilename = param_list_lookup_string(pl, "out");

    if (!rank && infilename) {
        in = fopen_maybe_compressed(infilename, "r");
        ASSERT_ALWAYS(in != nullptr);
    }
    if (!rank && outfilename) {
        out = fopen_maybe_compressed(outfilename, "w");
        ASSERT_ALWAYS(out != nullptr);
    }

    param_list_parse_uint(pl, "b", &batch_size);

    verbose_interpret_parameters(pl);

    char const * sm_mode_string = param_list_lookup_string(pl, "sm-mode");

    int nthr = 1;
    param_list_parse(pl, "t", nthr);

    if (param_list_warn_unused(pl))
        usage(argv0, nullptr, pl);

    if (!rank)
        param_list_print_command_line(stdout, pl);

    std::vector<sm_side_info> sm_info;

    for (int side = 0; side < cpoly->nb_polys; side++) {
        sm_info.emplace_back(F[side], ell, 0);
        sm_info[side].set_mode(sm_mode_string);
        if (nsm_arg[side] >= 0)
            sm_info[side].nsm = nsm_arg[side]; /* command line wins */
        if (!rank)
            printf("# Using %d SMs on side %d\n", sm_info[side].nsm, side);
    }

    /*
       if (!rank) {
       for (int side = 0; side < cpoly->nb_polys; side++) {
       printf("\n# Polynomial on side %d:\nF[%d] = ", side, side);
       mpz_poly_fprintf(stdout, F[side]);

       printf("# SM info on side %d:\n", side);
       sm_info[side].print(stdout);

       fflush(stdout);
       }
       }
       */

    sm_append(in, out, sm_info, nthr);

    /* Make sure we print no footer line, because reconstructlog-dl won't
     * grok it */
    if (!rank) {
        fflush(stdout);
    }

    if (!rank && infilename)
        fclose_maybe_compressed(in, infilename);
    if (!rank && out != stdout)
        fclose_maybe_compressed(out, outfilename);

    mpz_clear(ell);
    cado_poly_clear(cpoly);
    param_list_clear(pl);

    MPI_Finalize();

    return 0;
}
