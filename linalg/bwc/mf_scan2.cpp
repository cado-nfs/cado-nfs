#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>

#include <algorithm>
#include <atomic>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#ifdef HAVE_HWLOC
#include <hwloc.h>
#endif

#include "fmt/base.h"

#include "ringbuf.h"
#include "macros.h"          // for ASSERT_ALWAYS, MAX, MIN
#include "params.h"     // param_list
#include "timing.h"     // wct_seconds
#include "misc.h"       // size_disp
#include "fix-endianness.h" // fwrite32_little
#include "utils_cxx.hpp"        // for unique_ptr<FILE, delete_FILE>

static void mf_scan2_decl_usage(cxx_param_list & pl)
{
    param_list_usage_header(pl,
            "This program make one reading pass through a binary matrix, and produces\n"
            "the companion .rw and .cw files.\n"
            "Typical usage:\n"
            "\tmf_scan2 [<matrix file name> | options...]\n"
            );
    param_list_decl_usage(pl, "withcoeffs", "Handle DLP matrix, with coefficients\n");
    param_list_decl_usage(pl, "mfile", "Input matrix name (free form also accepted)");
    param_list_decl_usage(pl, "rwfile", "Name of the row weight file to write (defaults to auto-determine from matrix name)");
    param_list_decl_usage(pl, "cwfile", "Name of the col weight file to write (defaults to auto-determine from matrix name)");
    param_list_decl_usage(pl, "threads", "Number of threads to use (defaults to auto detect\n");
    param_list_decl_usage(pl, "io-memory", "Amount of RAM to use for rolling buffer memory (in GB, floating point allowed). Defaults to min(16M, 1/64-th of RAM (with hwloc))");
    param_list_decl_usage(pl, "thread-private-count", "Number of columns for which a thread-private zone is used");
    param_list_decl_usage(pl, "thread-read-window", "Chunk size for consumer thread reads from rolling buffer");
    param_list_decl_usage(pl, "thread-write-window", "Chunk size for producer thread writes to rolling buffer");
}

static size_t thread_private_count = 1UL << 20;

/* These two are in bytes as far as the default value is concerned, but
 * they're converted to number of uint32_t's when the program runs.
 */
static size_t thread_read_window = 1UL << 13;
static size_t thread_write_window = 1UL << 10;

static int withcoeffs = 0;

class reporter {
    std::atomic<size_t> produced;
    std::atomic<size_t> consumed;
    double t0 = 0;
    double last_report = 0;
    double delay = 1;
    std::mutex m;
    void report(bool force = false) {
        std::lock_guard<std::mutex> const dummy(m);
        double const tt = wct_seconds();
        bool const want = tt >= last_report + delay;
        if (!force && !want) return;
        fmt::print("read {}, parsed {}, in {:.1f} s\n",
                size_disp(produced),
                size_disp(consumed.load()),
                (last_report = tt) - t0);
        if (want)
            delay *= 1.189207;
    }
    public:
    struct consumer_data {
        double last;
        size_t s = 0;
        consumer_data() : last(wct_seconds()) {}
    };
    /* tells wheter this consumer has reason to schedule a new report */
    void consumer_report(consumer_data & D, size_t s, bool force = false) {
        double const tt = wct_seconds();
        if (!force && tt < D.last + 0.9) {
            D.s += s;
            return;
        }
        consumed += D.s + s;
        D.last = tt;
        D.s = 0;
        report();
    }
    void producer_report(size_t s, bool force = false) {
        produced += s;
        report(force);
    }
    void reset() { t0 = last_report = wct_seconds(); }
};


static reporter report;

static inline unsigned int get_segment_index(uint32_t c)
{
    return 64 - cado_clz64((uint64_t) c);
}
static inline uint32_t get_segment_offset(unsigned int t) { return 1UL << (t-1); }
static inline uint32_t get_segment_size(unsigned int t) { return 1UL << (t-1); }

struct segment {
    std::vector<std::atomic<uint32_t>> data;
    static_assert(sizeof(decltype(data)::value_type) == sizeof(uint32_t), "size mismatch");
    static size_t segment_size(unsigned int t) { return get_segment_size(t); }
    explicit segment(unsigned int t)
        : data(segment_size(t))
    {
    }
    segment(segment const&) = delete;
    segment& operator=(segment const&) = delete;
    segment(segment &&) = delete;
    segment& operator=(segment &&) = delete;
    ~segment() = default;
    void incr(uint32_t c) { data[c]++; }
    uint32_t * non_atomic_data() { return (uint32_t *) (data.data()); }
};

/* It might seem somewhat overkill to use std::atomic here. Some of the
 * associated fencing is quite probably overkill on x86. But I'm not too
 * sure.
 *
 * I've added some loose memory_order constraints below, that seem to
 * improve performance. But I'm on thin ice, I'm not sure of what I'm
 * doing.
 *
 * (the reassuring thing is that I _think_ that the worst that can happen
 * is a seg fault, which would be loud enough, and therefore fine).
 *
 * Some pointers:
 *
 * https://bartoszmilewski.com/2008/12/01/c-atomics-and-memory-ordering/
 * https://bartoszmilewski.com/2008/12/23/the-inscrutable-c-memory-model/
 * http://www.cplusplus.com/reference/atomic/memory_order/
 */
static std::atomic<segment *> segments[64];
static std::mutex segment_mutexes[64];

template<bool> struct nz_coeff; // IWYU pragma: keep
template<> struct nz_coeff<false> { struct type { uint32_t j; }; };
template<> struct nz_coeff<true> { struct type { uint32_t j; int32_t v; }; };

template<bool withcoeffs>
struct parser_thread {
    using coeff_t = typename nz_coeff<withcoeffs>::type;
    std::vector<uint32_t> cw;
    uint32_t colmax=0;
    parser_thread() : cw(thread_private_count, 0) {};
    void loop(ringbuf_ptr R) {
        reporter::consumer_data D;
        coeff_t buffer[thread_read_window];
        for(size_t s ; (s = ringbuf_get(R, (char*) buffer, sizeof(buffer))) != 0 ; ) {
            report.consumer_report(D, s);
            auto * v = (coeff_t *) buffer;
            ASSERT_ALWAYS(s % sizeof(coeff_t) == 0);
            size_t const sv = s / sizeof(coeff_t);
            for(size_t i = 0 ; i < sv ; i++) {
                uint32_t const c = v[i].j;
                colmax = MAX(colmax, c+1);
                if (c < thread_private_count) {
                    cw[c]++;
                } else {
                    /* Get the bit size */
                    unsigned int const t = get_segment_index(c);
                    uint32_t const c1 = c-get_segment_offset(t);
                    ASSERT_ALWAYS(t < 64);
                    ASSERT_ALWAYS(c1 < get_segment_size(t));
                    segment * x;
                    {
                        x = segments[t].load(std::memory_order_relaxed);
                        if (!x) {
                            std::lock_guard<std::mutex> const dummy(segment_mutexes[t]);
                            x = segments[t].load(std::memory_order_relaxed);
                            if (!x)
                                segments[t].store(x = new segment(t), std::memory_order_relaxed);
                        }
                    }
                    x->incr(c1);
                }
            }
        }
        report.consumer_report(D, 0, true);
    }
};

template<bool withcoeffs>
static void master_loop(ringbuf_ptr R, FILE * f_in, FILE * f_rw)
{
    constexpr int c = withcoeffs != 0;
    using coeff_t = typename nz_coeff<withcoeffs>::type;
    ASSERT_ALWAYS(thread_write_window % sizeof(coeff_t) == 0);
    size_t const tw = thread_write_window / sizeof(coeff_t);
    coeff_t buf[tw];
    for( ; ; ) {
        uint32_t row_length;
        {
            auto const rc = fread32_little(&row_length, 1, f_in);
            if (rc != 1)
                break;
        }
        if (((int32_t)row_length) < 0) {
            fmt::print(stderr, "Found row with more than 2G entries."
                    " You most probably omitted the --withcoeffs flag\n");
            exit(EXIT_FAILURE);
        }
        {
            auto const rc = fwrite32_little(&row_length, 1, f_rw);
            ASSERT_ALWAYS(rc == 1);
        }
        for( ; row_length ; ) {
            size_t const s = std::min((size_t) row_length, tw);
            auto k = fread32_little((uint32_t *) buf, s * (1 + c), f_in);
            if (k != s * (1+c)) {
                fmt::print(stderr,
                        "Input error while reading {} 32-bit entries"
                        " from fd {} at position {}:"
                        " Got only {} values\n",
                        s * (1+c),
                        fileno(f_in),
                        ftell(f_in),
                        k);
                abort();
            }
            ASSERT_ALWAYS(k == s * (1 + c));
            ringbuf_put(R, (char *) buf, s * sizeof(coeff_t));
            report.producer_report(s * sizeof(coeff_t));
            row_length -= s;
        }
    }
    ringbuf_mark_done(R);
}

static void finish_write_and_clear_segments(uint32_t c, uint32_t colmax, FILE * f_cw)
{
    for( ; c < colmax ; ) {
        unsigned int const t = get_segment_index(c);
        ASSERT_ALWAYS(t < 64);
        std::lock_guard<std::mutex> const dummy(segment_mutexes[t]);
        uint32_t const c1 = c-get_segment_offset(t);
        uint32_t const max1 = MIN(colmax-get_segment_offset(t), get_segment_size(t));
        uint32_t const n1 = max1 - c1;
        segment * x = segments[t];
        if (!x) x = new segment(t);
        size_t const rc = fwrite32_little(x->non_atomic_data() + c1, n1, f_cw);
        ASSERT_ALWAYS(rc == n1);
        c += n1;
        delete x;
    }
}

template<bool withcoeffs>
static void write_column_weights(std::vector<parser_thread<withcoeffs>> & T, FILE * f_cw)
{
    for(size_t i = 1 ; i < T.size() ; i++) {
        for(size_t j = 0 ; j < thread_private_count ; j++)
            T[0].cw[j] += T[i].cw[j];
        T[0].colmax = MAX(T[0].colmax, T[i].colmax);
    }
    uint32_t const colmax = T[0].colmax;
    uint32_t c = 0;
    for( ; c < thread_private_count && c < colmax ; c++) {
        auto const rc = fwrite32_little(&T[0].cw[c], 1, f_cw);
        ASSERT_ALWAYS(rc == 1);
    }
    finish_write_and_clear_segments(c, colmax, f_cw);
}

template<bool withcoeffs>
static void maincode(ringbuf_ptr R, int nb_consumers, FILE * f_in, FILE * f_rw, FILE * f_cw)
{
    std::vector<parser_thread<withcoeffs>> T(nb_consumers);
    std::thread producer(master_loop<withcoeffs>, R, f_in, f_rw);
    std::vector<std::thread> consumers;
    consumers.reserve(T.size());
    for(auto & t : T)
        consumers.emplace_back(&parser_thread<withcoeffs>::loop, std::ref(t), R);

    producer.join();
    report.producer_report(0, true);
    for(auto & t : consumers)
        t.join();
    write_column_weights(T, f_cw);
}

// coverity[root_function]
int main(int argc, char const * argv[])
{
    const char * argv0 = argv[0];

    cxx_param_list pl;
    std::string rwfile;
    std::string cwfile;
    std::string mfile;

    unsigned int wild =  0;

    argv++,argc--;

    mf_scan2_decl_usage(pl);

    param_list_configure_switch(pl, "--withcoeffs", &withcoeffs);

    for(;argc;) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        if (argv[0][0] != '-' && wild == 0) {
            mfile = argv[0];
            wild++;
            argv++,argc--;
            continue;
        }
        fmt::print(stderr, "unknown option {}\n", argv[0]);
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    param_list_parse(pl, "thread-private-count", thread_private_count);
    param_list_parse(pl, "thread-read-window", thread_read_window);
    param_list_parse(pl, "thread-write-window", thread_write_window);

    param_list_parse(pl, "mfile", mfile);
    param_list_parse(pl, "rwfile", rwfile);
    param_list_parse(pl, "cwfile", cwfile);

    if (mfile.empty()) {
        param_list_print_usage(pl, argv0, stderr);
        exit(EXIT_FAILURE);
    }

    if (!ends_with(mfile, ".bin")) {
        fmt::print(stderr, "Warning: matrix file name should end in .bin\n");
    }

    if (rwfile.empty()) {
        rwfile = std::unique_ptr<char, free_delete<char>>(
                derived_filename(mfile.c_str(), "rw", ".bin")).get();
    }

    if (cwfile.empty()) {
        cwfile = std::unique_ptr<char, free_delete<char>>(
                derived_filename(mfile.c_str(), "cw", ".bin")).get();
    }

#ifdef HAVE_HWLOC
    /* Detect hardware */
    hwloc_topology_t topology;
    hwloc_topology_init(&topology);
    hwloc_topology_load(topology);
    int const depth = hwloc_topology_get_depth(topology);
    unsigned int const npu = hwloc_get_nbobjs_by_depth(topology, depth-1);
    hwloc_obj_t root = hwloc_get_root_obj(topology);
#if HWLOC_API_VERSION < 0x020000
    uint64_t ram = root->memory.total_memory;
#else
    uint64_t ram = root->total_memory;
#endif
    for(uint64_t x = ram >> 4; x ; x >>= 1) ram |= x;
    ram = ram + 1;
    hwloc_topology_destroy(topology);
    int threads = int(npu);
    {
        /* our test environment might force this to a small value. It's a
         * priori also valid for plain threads. */
        const char * tmp = getenv("CADO_NFS_MAX_THREADS");
        if (tmp) {
            int const n = (int) strtoul(tmp, nullptr, 0);
            threads = std::min(threads, n);
        }
    }
#else
    uint64_t ram = uint64_t(1) << 30;
    int threads = 2;
#endif

    size_t ringbuf_size = std::min(ram / 64, UINT64_C(1) << 24);

    param_list_parse(pl, "threads", threads);
    {
        double r;
        if (param_list_parse(pl, "io-memory", r)) {
            ringbuf_size = size_t(r * double(1UL << 30));
        }
    }

    ringbuf R;
    ringbuf_init(R, ringbuf_size);

    /* Start with the input */

    std::unique_ptr<FILE, delete_FILE> const f_in(fopen(mfile.c_str(), "rb"));
    if (!f_in) { perror(mfile.c_str()); exit(EXIT_FAILURE); }
    std::unique_ptr<FILE, delete_FILE> const f_rw(fopen(rwfile.c_str(), "wb"));
    if (!f_rw) { perror(rwfile.c_str()); exit(EXIT_FAILURE); }
    std::unique_ptr<FILE, delete_FILE> const f_cw(fopen(cwfile.c_str(), "wb"));
    if (!f_cw) { perror(cwfile.c_str()); exit(EXIT_FAILURE); }

    ASSERT_ALWAYS(threads >= 2);

    ASSERT_ALWAYS(thread_read_window % sizeof(uint32_t) == 0);
    thread_read_window  /= sizeof(uint32_t);
    report.reset();

    int const consumers = threads-1;

    if (!withcoeffs) {
        maincode<false>(R, consumers, f_in.get(), f_rw.get(), f_cw.get());
    } else {
        maincode<true>(R, consumers, f_in.get(), f_rw.get(), f_cw.get());
    }

    ringbuf_clear(R);
}

