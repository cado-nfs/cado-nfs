#ifndef CADO_PARALLELIZING_INFO_HPP
#define CADO_PARALLELIZING_INFO_HPP

/* don't enable this. Clutters output a lot */
#define xxxCONCURRENCY_DEBUG

#include <cstdio>
#include <cstddef>

#include <array>
#include <map>
#include <vector>
#include <utility>
#include <memory>
#include <set>
#ifdef CONCURRENCY_DEBUG
#include <mutex>
#endif

#include <sys/time.h>

#include "barrier.h"
#include "macros.h"
#include "params.hpp"
#include "select_mpi.h"
#include "arith-generic.hpp"
#include "lock_guarded_container.hpp"

/*
 * The main guy here is the parallelizing_info data type. It is commonly
 * declared as pi.
 *
 * Threads and jobs are arranged in a grid-like manner. Several logical
 * communication groups are defined. One for the global grid, as well as
 * one for each column / row.
 *
 * pi->m denotes the global communicator. There is only one such
 * communicator. All jobs, all threads contribute to it. At the mpi
 * level, the related communicator is MPI_COMM_WORLD (well, unless we're
 * working in interleaving mode).
 *
 * Two other communicators are defined:
 *      pi->wr[0] denotes horizontal, a.k.a. ROW groups.
 *      pi->wr[1] denotes vertical, a.k.a. COLUMN groups.
 * [Note: communicators used to be called "wirings", hence the variable
 * name]
 *
 * When a matrix is mapped to a grid process of this kind, say a matrix
 * having been split in nhs * nvs slices, then there are exactly nhs ROW
 * groups, and nvs COLUMN groups. ROW groups consist of nvs (job,thread)
 * items (as many as one finds COL groups), and conversely.
 */

#ifdef CONCURRENCY_DEBUG
extern std::mutex stdio_mutex;
#endif

/*
 * MPI_LIBRARY_MT_CAPABLE: do mpi communications in a furiously
 * concurrent manner.
 *
 * This isn't ready yet. In presence of an MPI library which correctly
 * supports MPI_THREAD_MULTIPLE, chances are that this gets close to
 * working ok, and even maybe improve performance quite a bit.
 *
 * As of openmpi-1.8, this is not reliable enough when the basic
 * transport layer is openib (on infiniband). Status may of course vary
 * for other mpi implementations.
 *
 * As for the ``fake mpi'' implementation we have, well, it's up to the
 * reader to decide. All collectives are nops, so it's trivially capable.
 * OTOH, not setting the flags ensures that the rest of the code compiles
 * ok.
 *
 * --> we never allow this flag currently. Corresponding pieces of code
 *  have been deleted as the program evolved, given that it had zero
 *  testing.
 */

#if defined(MPICH2) && MPICH2_NUMVERSION >= 10100002
/* In fact, even in this case we might consider disabling it. */
#define xxxMPI_LIBRARY_MT_CAPABLE
#elif defined(OPEN_MPI) && OMPI_VERSION_ATLEAST(1,8,2)
#define xxxMPI_LIBRARY_MT_CAPABLE
/*
 * at present I know of no version of openmpi with MPI_THREAD_MULTIPLE
 * working, but to be honest I haven't tried hard. For sure there are
 * some bugs in my code as well anyway, at least that's what enabling it
 * shows.
 */
#else
/* Assume it does not work */
#endif

/* {{{ definition of the parallelizing_info communicator type */

/* {{{ utility structures for communicators. */

/* {{{ This one is stored in thread-shared memory ; shared locks and so on */
struct pthread_things {
    barrier_t bh[1];
    my_pthread_barrier_t b[1];

    pthread_mutex_t m[1];
    std::string desc;
    void * utility_ptr = nullptr;
    // int count;
};
/* }}} */

/* {{{ logging. To be activated in debug mode only */
struct pi_log_entry {
    struct timeval tv[1];
    char what[80];
};

#define PI_LOG_BOOK_ENTRIES     32
struct pi_log_book {
    struct pi_log_entry t[PI_LOG_BOOK_ENTRIES];
    int hsize = 0;  // history size -- only a count, once the things wraps.
    int next = 0;   // next free pointer.
};
/* }}} */
/* }}} */

struct parallelizing_info;

struct pi_datatype;

struct pi_op_s;
using pi_op = pi_op_s;
using pi_op_ptr = pi_op *;

struct pi_comm { /* {{{ */
    /* njobs : number of mpi jobs concerned by this logical group */
    /* ncores : number of threads concerned by this logical group */
    unsigned int njobs = 0;
    unsigned int ncores = 0;

    /* product njobs * ncores */
    unsigned int totalsize = 0;
    unsigned int jrank = 0;
    unsigned int trank = 0;
    MPI_Comm pals = 0;

    struct pthread_things * th = nullptr;
#ifdef  CONCURRENCY_DEBUG
    int th_count = 0;
#endif
    std::unique_ptr<pi_log_book> log_book;

    /* This is the communicator in "the other direction". If this pointer
     * is not NULL, then there are exactly this->xwr->ncores such that
     * this->trank == 0
     */
    pi_comm * xwr = nullptr;

private:
    pi_comm& operator=(pi_comm const & o)
    {
        if (this == &o)
            return *this;
        njobs = o.njobs;
        ncores = o.ncores;
        totalsize = o.totalsize;
        jrank = o.jrank;
        pals = o.pals;
        th = nullptr;
#ifdef  CONCURRENCY_DEBUG
        th_count;
#endif
        return *this;
    }
public:

    /* pointers must be different on all threads */
    void thread_bcast(void * sendbuf, size_t count, pi_datatype * datatype, unsigned int root);
    void thread_allreduce(void * sendbuf, void * recvbuf, size_t count, pi_datatype * datatype, pi_op_ptr op);
    int thread_data_eq(void * buffer, size_t count, pi_datatype * datatype);

    public:
    void bcast(void * sendbuf, size_t count, pi_datatype * datatype, unsigned int jroot, unsigned int troot);
    void allreduce(void * sendbuf, void * recvbuf, size_t count, pi_datatype * datatype, pi_op_ptr op);
    void allgather(void * sendbuf, size_t sendcount, pi_datatype * sendtype, void * recvbuf, size_t recvcount, pi_datatype * recvtype);
    int data_eq(void * buffer, size_t count, pi_datatype * datatype);

    void abort(int err) const;

    int serialize(const char * file, unsigned int line);
    int serialize_threads(const char * file, unsigned int line);

    /* shared_malloc is like malloc, except that the pointer returned will be
     * equal on all threads (proper access will deserve proper locking of
     * course). shared_malloc_set_zero sets to zero too */
    /* As a side-effect, all shared_* functions serialize threads */
    void * shared_malloc(size_t size);
    void * shared_malloc_set_zero(size_t size);
    void shared_free(void * ptr);

    template<typename T> T * shared_new(size_t size);
    template<typename T> T * shared_new();
    template<typename T> void shared_delete(T *, size_t);
    template<typename T> void shared_delete(T *);

    /* stuff related to log entry printing */
    void log_init();
    void log_clear();
    void log_op(const char * fmt, ...);
    void log_print();
    private:
    void init_pthread_things(const char * desc);
    void destroy_pthread_things();
    friend struct parallelizing_info;
};
using pi_comm_ptr = pi_comm *;
using pi_comm_srcptr = const pi_comm *;
/* }}} */

/* {{{ interleaving two pi structures. */
struct pi_interleaving {
    int idx = 0;                    /* 0 or 1 */
    my_pthread_barrier_t * b = nullptr;   /* not a 1-sized array on purpose --
                                   being next to index, it can't ! */
};
using pi_interleaving_ptr = pi_interleaving *;
/* }}} */

/* {{{ This arbitrary associative array is meant to be very global, even
 * common to two interleaved pi structures. Used to pass lightweight info
 * only */

using pi_dictionary = lock_guarded_container<std::map<std::pair<unsigned long, unsigned long>, void *>>;
/* }}} */

/* {{{ global parallelizing_info handle */
#define PI_NAMELEN      32
struct parallelizing_info {
    // row-wise, column-wise.
    pi_comm wr[2];
    // main.
    pi_comm m;
    pi_interleaving_ptr interleaved = nullptr;
    pi_dictionary * dict = nullptr;
    std::string nodename;
    std::string nodeprefix;
    std::string nodenumber;
    /* This pointer is identical on all threads. It is non-null only in
     * case we happen to have sufficiently recent gcc, together with
     * sufficiently recent hwloc */
    void * cpubinding_info = nullptr;
    std::array<int, 2> thr_orig {{0, 0}};/* when only_mpi is 1, this is what the
                                   thr parameter was set to originally.
                                   Otherwise we have {0,0} here. */

    void hello();
    void log_print_all() const;
    void grid_print(char * buf, size_t siz, int print);
    void store_generic(unsigned long k1, unsigned long k2, void * val);
    void * load_generic(unsigned long k1, unsigned long k2);

    void interleaving_flip();
    void interleaving_enter();
    void interleaving_leave();
    static void init_attribute_things();
    static void clear_attribute_things();

    /* we define new datatypes in a way which diverts from the mpi calling
     * interface, because that interface is slightly awkward for our needs */
    pi_datatype * alloc_arith_datatype(arith_generic * abase);
    void free_arith_datatype(pi_datatype *);

    /* private: */ /* we'd like to! */
    parallelizing_info * grid_init();
    void grid_clear(parallelizing_info *);
    void clear_mpilevel();
    void init_mpilevel(cxx_param_list & pl);

    static void declare_usage(cxx_param_list & pl);
    static void lookup_parameters(cxx_param_list & pl);
};
/* }}} */

/* {{{ collective operations and user-defined types */

struct pi_datatype {
    MPI_Datatype datatype;
    /* two attributes we're really happy to use */
    arith_generic * abase;
    size_t item_size;
};

extern pi_datatype * BWC_PI_INT;
extern pi_datatype * BWC_PI_DOUBLE;
extern pi_datatype * BWC_PI_BYTE;
extern pi_datatype * BWC_PI_UNSIGNED;
extern pi_datatype * BWC_PI_UNSIGNED_LONG;
extern pi_datatype * BWC_PI_UNSIGNED_LONG_LONG;
extern pi_datatype * BWC_PI_LONG;
extern pi_datatype * BWC_PI_SIZE_T;

struct pi_op_s {
    MPI_Op stock;  /* typically MPI_SUM */
    MPI_Op custom = MPI_OP_NULL;  /* for arith types, the mpi-level user-defined op */
    using MPI_Op_t = void (*)(void *, void *, int *, MPI_Datatype *);
    void (*f_stock)(arith_generic::elt const *, arith_generic::elt *, int *, MPI_Datatype *) = nullptr;
    void (*f_custom)(arith_generic::elt const *, arith_generic::elt *, size_t, pi_datatype *) = nullptr;
    pi_op_s(MPI_Op s) : stock(s) {}
};
using pi_op = pi_op_s;
using pi_op_ptr = pi_op *;

extern struct pi_op_s BWC_PI_MIN[1];
extern struct pi_op_s BWC_PI_MAX[1];
extern struct pi_op_s BWC_PI_SUM[1];
extern struct pi_op_s BWC_PI_BXOR[1];
extern struct pi_op_s BWC_PI_BAND[1];
extern struct pi_op_s BWC_PI_BOR[1];

/* This _only_ works if the datatype has been registered with
 * parallelizing_info::alloc_arith_datatype
 */
extern arith_generic * pi_arith_datatype_get_abase(MPI_Datatype datatype);

/* }}} */

/* {{{ I/O layer */
struct pi_file_handle {
    std::string name;        /* just for reference. I doubt we'll need them */
    std::string mode;
    FILE * f = nullptr;   /* meaningful only at root */
    parallelizing_info & pi;
    int inner = 0;
    int outer = 0;

    pi_file_handle(parallelizing_info & pi)
        : pi(pi)
    {}

    int open(int inner, std::string const & name, const char * mode);
    ~pi_file_handle() {
        if (f) close();
    }
    int close();
    size_t write(void * buf, size_t size, size_t sizeondisk);
    size_t read(void * buf, size_t size, size_t sizeondisk);
    size_t write_chunk(void * buf, size_t size, size_t sizeondisk, size_t chunksize, size_t spos, size_t epos);
    size_t read_chunk(void * buf, size_t size, size_t sizeondisk, size_t chunksize, size_t spos, size_t epos);
};
/* }}} */


/* pi_go is the main function. It is responsible of creating all the
 * parallelizing_info data structures, set up the different inter-job and
 * inter-thread conciliation toys (communicators, pthread barriers, and
 * so on), and eventually run the provided function.
 *
 * the param_list is checked for parameters mpi and thr, so as to define
 * the mpi and thr splttings.
 *
 * nhc, nvc are the same for threads (cores).
 */
extern void pi_go(
        void *(&&fcn)(parallelizing_info &, cxx_param_list & pl, void * arg),
        cxx_param_list & pl,
        void * arg);

extern void pi_hello(parallelizing_info & pi);

/* the parallelizing_info layer has some collective operations which
 * deliberately have prototypes simlar or identical to their mpi
 * counterparts (we use size_t for the count arguments, though).
 *
 * Note: These functions use the wr->utility_ptr field.
 */

/* almost similar to the mpi-level reduction functions, except that
 * we added const, dropped the int *, and dropped the multiplexing
 * argument with the datatype (perhaps we shouldn't, but so far we have
 * no use for it) */
typedef void (*thread_reducer_t)(const void *, void *, size_t);

/* These two interfaces are experimental only */

namespace parallelizing_info_experimental {
    void broadcast(std::vector<unsigned int>& v, parallelizing_info & pi);
    void allgather(std::vector<unsigned int>& v, pi_comm & wr);
    void broadcast(std::set<unsigned int>& v, parallelizing_info & pi);
    void allgather(std::set<unsigned int>& v, pi_comm & wr);
}

/* prints the given string in a ascii-form matrix. */
extern void grid_print(parallelizing_info & pi, char * buf, size_t siz, int print);

/* These are the calls for interleaving. The 2n threads are divided into
 * two grous. It is guaranteed that at a given point, the two groups of n
 * threads are separated on either size of the pi_interleaving_flip call.
 *
 * The called function must make sure that alternating blocks (delimited
 * by _flip) either do or don't contain mpi calls, IN TURN.
 *
 * _enter and _leave are called from pi_go, so although they are exposed,
 * one does not have to know about them.
 */
extern void pi_interleaving_flip(parallelizing_info &);
extern void pi_interleaving_enter(parallelizing_info &);
extern void pi_interleaving_leave(parallelizing_info &);

extern void pi_store_generic(parallelizing_info &, unsigned long, unsigned long, void *);
extern void * pi_load_generic(parallelizing_info &, unsigned long, unsigned long);


template<typename T>
inline T * pi_comm::shared_new(size_t size)
{
    void * ptr = nullptr;
    if (trank == 0) ptr = new T[size];
    thread_bcast(&ptr, sizeof(void*), BWC_PI_BYTE, 0);
    return static_cast<T *>(ptr);
}

template<typename T>
inline void pi_comm::shared_delete(T * ptr, size_t)
{
    serialize_threads(__FILE__, __LINE__);
    if (trank == 0) delete[] ptr;
}

template<typename T>
inline T * pi_comm::shared_new()
{
    void * ptr = nullptr;
    if (trank == 0) ptr = new T;
    thread_bcast(&ptr, sizeof(void*), BWC_PI_BYTE, 0);
    return static_cast<T *>(ptr);
}

template<typename T>
inline void pi_comm::shared_delete(T * ptr)
{
    serialize_threads(__FILE__, __LINE__);
    if (trank == 0) delete ptr;
}

template<typename T>
struct shared_delete_deleter {
    pi_comm * wr = nullptr;
    shared_delete_deleter() = default;
    ~shared_delete_deleter() = default;
    explicit shared_delete_deleter(pi_comm * wr) : wr(wr) {}
    shared_delete_deleter(shared_delete_deleter const &) = default;
    shared_delete_deleter& operator=(shared_delete_deleter const &) = default;
    shared_delete_deleter(shared_delete_deleter &&) noexcept = default;
    shared_delete_deleter& operator=(shared_delete_deleter &&) noexcept = default;
    void operator()(T * ptr) const { if (wr) wr->shared_delete(ptr); }
};

template<typename T>
struct shared_delete_deleter<T[]> {
    pi_comm * wr = nullptr;
    shared_delete_deleter() = default;
    ~shared_delete_deleter() = default;
    explicit shared_delete_deleter(pi_comm_ptr wr) : wr(wr) {}
    shared_delete_deleter(shared_delete_deleter const &) = default;
    shared_delete_deleter& operator=(shared_delete_deleter const &) = default;
    shared_delete_deleter(shared_delete_deleter &&) noexcept = default;
    shared_delete_deleter& operator=(shared_delete_deleter &&) noexcept = default;
    void operator()(T * ptr) const { if (wr) wr->shared_delete(ptr, 0); }
};

template<typename T>
struct pi_shared_array : public std::unique_ptr<T[], shared_delete_deleter<T[]>> {
    using super = std::unique_ptr<T[], shared_delete_deleter<T[]>>;
    /* the default ctor attaches no communicator, of course */
    pi_shared_array() = default;
    explicit pi_shared_array(pi_comm & wr, unsigned int n)
        : super(wr.shared_new<T>(n), shared_delete_deleter<T[]>(&wr))
    {}
};
template<typename T>
struct pi_shared_object : public std::unique_ptr<T, shared_delete_deleter<T>> {
    using super = std::unique_ptr<T, shared_delete_deleter<T>>;
    /* the default ctor attaches no communicator, of course */
    pi_shared_object() = default;
    explicit pi_shared_object(pi_comm & wr)
        : super(wr.shared_new<T>(), shared_delete_deleter<T>(&wr))
    {}
};

/* This provides a fairly typical construct, used like this:
 * SEVERAL_THREADS_PLAY_MPI_BEGIN(some pi communicator) {
 *      some code which happens for threads in the pi comm in turn
 * }
 * SEVERAL_THREADS_PLAY_MPI_END();
 */
#ifndef MPI_LIBRARY_MT_CAPABLE
#define SEVERAL_THREADS_PLAY_MPI_BEGIN(comm) do {			\
    for(unsigned int t__ = 0 ; t__ < (comm)->ncores ; t__++) {		\
        (comm)->serialize_threads(__FILE__, __LINE__);			\
        if (t__ != (comm)->trank) continue; /* not our turn. */         \
        do
#define SEVERAL_THREADS_PLAY_MPI_END()    while (0); } } while (0)
/* This construct is used similarly. It differs slightly, in that we
 * guarantee that only one thread (in the communicator) will issue mpi
 * calls */
#define SEVERAL_THREADS_PLAY_MPI_BEGIN2(comm, t__) do { 		\
    (comm)->serialize_threads(__FILE__, __LINE__);                      \
    if ((comm)->trank == 0) {                                           \
        for(unsigned int t__ = 0 ; t__ < (comm)->ncores ; t__++) {	\
            do
#define SEVERAL_THREADS_PLAY_MPI_END2(comm)                             \
            while (0);							\
        }								\
    }									\
    (comm)->serialize_threads(__FILE__, __LINE__);                      \
} while (0)
#else
#define SEVERAL_THREADS_PLAY_MPI_BEGIN(comm)     /**/
#define SEVERAL_THREADS_PLAY_MPI_END()           /**/
#define SEVERAL_THREADS_PLAY_MPI_BEGIN2(comm, t) /**/
#define SEVERAL_THREADS_PLAY_MPI_END2(comm)      /**/
#endif

#endif	/* PARALLELIZING_INFO_HPP_ */
