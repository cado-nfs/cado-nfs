#ifndef PARALLELIZING_INFO_HPP_
#define PARALLELIZING_INFO_HPP_

/* don't enable this. Clutters output a lot */
#define xxxCONCURRENCY_DEBUG

#include <cstdio>            // for FILE

#include <map>
#include <vector>
#include <set>
#ifdef CONCURRENCY_DEBUG
#include <mutex>
#endif

#include <sys/types.h>        // for ssize_t
#include <sys/time.h>         // for struct timeval

#include "barrier.h"          // for barrier_t
#include "macros.h"           // for OMPI_VERSION_ATLEAST
#include "params.h"
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
    char * desc;
    void * utility_ptr;
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
    int hsize;  // history size -- only a count, once the things wraps.
    int next;   // next free pointer.
};
/* }}} */
/* }}} */

struct parallelizing_info_s;

struct pi_comm_s;
typedef struct pi_comm_s * pi_comm_ptr;
typedef const struct pi_comm_s * pi_comm_srcptr;

struct pi_comm_s { /* {{{ */
    /* njobs : number of mpi jobs concerned by this logical group */
    /* ncores : number of threads concerned by this logical group */
    unsigned int njobs;
    unsigned int ncores;

    /* product njobs * ncores */
    unsigned int totalsize;
    unsigned int jrank;
    unsigned int trank;
    MPI_Comm pals;

    struct pthread_things * th;
#ifdef  CONCURRENCY_DEBUG
    int th_count;
#endif
    struct pi_log_book * log_book;

    /* This is the communicator in "the other direction". If this pointer
     * is not NULL, then there are exactly this->xwr->ncores such that
     * this->trank == 0
     */
    pi_comm_ptr xwr;
};
typedef struct pi_comm_s pi_comm[1];
/* }}} */
/* }}} */

/* {{{ interleaving two pi structures. */
struct pi_interleaving_s {
    int idx;                    /* 0 or 1 */
    my_pthread_barrier_t * b;   /* not a 1-sized array on purpose --
                                   being next to index, it can't ! */
};
typedef struct pi_interleaving_s pi_interleaving[1];
typedef struct pi_interleaving_s * pi_interleaving_ptr;
/* }}} */

/* {{{ This arbitrary associative array is meant to be very global, even
 * common to two interleaved pi structures. Used to pass lightweight info
 * only */

typedef lock_guarded_container<std::map<std::pair<unsigned long, unsigned long>, void *>> pi_dictionary;
/* }}} */

/* {{{ global parallelizing_info handle */
#define PI_NAMELEN      32
struct parallelizing_info_s {
    // row-wise, column-wise.
    pi_comm wr[2];
    // main.
    pi_comm m;
    pi_interleaving_ptr interleaved;
    pi_dictionary * dict;
    char nodename[PI_NAMELEN];
    char nodeprefix[PI_NAMELEN];
    char nodenumber_s[PI_NAMELEN];
    /* This pointer is identical on all threads. It is non-null only in
     * case we happen to have sufficiently recent gcc, together with
     * sufficiently recent hwloc */
    void * cpubinding_info;
    int thr_orig[2];            /* when only_mpi is 1, this is what the
                                   thr parameter was set to originally.
                                   Otherwise we have {0,0} here. */
};

typedef struct parallelizing_info_s parallelizing_info[1];
typedef struct parallelizing_info_s * parallelizing_info_ptr;
typedef const struct parallelizing_info_s * parallelizing_info_srcptr;
/* }}} */

/* {{{ collective operations and user-defined types */

struct pi_datatype_s {
    MPI_Datatype datatype;
    /* two attributes we're really happy to use */
    arith_generic * abase;
    size_t item_size;
};
typedef struct pi_datatype_s * pi_datatype_ptr;
extern pi_datatype_ptr BWC_PI_INT;
extern pi_datatype_ptr BWC_PI_DOUBLE;
extern pi_datatype_ptr BWC_PI_BYTE;
extern pi_datatype_ptr BWC_PI_UNSIGNED;
extern pi_datatype_ptr BWC_PI_UNSIGNED_LONG;
extern pi_datatype_ptr BWC_PI_UNSIGNED_LONG_LONG;
extern pi_datatype_ptr BWC_PI_LONG;
extern pi_datatype_ptr BWC_PI_SIZE_T;

struct pi_op_s {
    MPI_Op stock;  /* typically MPI_SUM */
    MPI_Op custom = MPI_OP_NULL;  /* for arith types, the mpi-level user-defined op */
    typedef void (*MPI_Op_t)(void *, void *, int *, MPI_Datatype *);
    void (*f_stock)(arith_generic::elt const *, arith_generic::elt *, int *, MPI_Datatype *) = NULL;
    void (*f_custom)(arith_generic::elt const *, arith_generic::elt *, size_t, pi_datatype_ptr) = NULL;
    pi_op_s(MPI_Op s) : stock(s) {}
};
typedef struct pi_op_s * pi_op_ptr;
extern struct pi_op_s BWC_PI_MIN[1];
extern struct pi_op_s BWC_PI_MAX[1];
extern struct pi_op_s BWC_PI_SUM[1];
extern struct pi_op_s BWC_PI_BXOR[1];
extern struct pi_op_s BWC_PI_BAND[1];
extern struct pi_op_s BWC_PI_BOR[1];

/* we define new datatypes in a way which diverts from the mpi calling
 * interface, because that interface is slightly awkward for our needs */

extern pi_datatype_ptr pi_alloc_arith_datatype(parallelizing_info_ptr pi, arith_generic * abase);
extern void pi_free_arith_datatype(parallelizing_info_ptr pi, pi_datatype_ptr ptr);

/* This _only_ works if the datatype has been registered with
 * pi_alloc_arith_datatype
 */
extern arith_generic * pi_arith_datatype_get_abase(MPI_Datatype datatype);



/* }}} */

/* {{{ I/O layer */
struct pi_file_handle_s {
    char * name;        /* just for reference. I doubt we'll need them */
    char * mode;
    FILE * f;   /* meaningful only at root */
    parallelizing_info_ptr pi;
    int inner;
    int outer;
};
typedef struct pi_file_handle_s pi_file_handle[1];
typedef struct pi_file_handle_s * pi_file_handle_ptr;
/* }}} */

extern void parallelizing_info_init();
extern void parallelizing_info_finish();
extern void parallelizing_info_decl_usage(cxx_param_list & pl);
extern void parallelizing_info_lookup_parameters(cxx_param_list & pl);

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
        void *(&&fcn)(parallelizing_info_ptr, cxx_param_list & pl, void * arg),
        cxx_param_list & pl,
        void * arg);

extern void pi_hello(parallelizing_info_ptr pi);

/* I/O functions */
extern int pi_file_open(pi_file_handle_ptr f, parallelizing_info_ptr pi, int inner, const char * name, const char * mode);
extern int pi_file_close(pi_file_handle_ptr f);
/* totalsize is the size which should be on disk. It may be shorter than
 * the sum of the individual sizes, in case of padding */
extern ssize_t pi_file_write(pi_file_handle_ptr f, void * buf, size_t size, size_t totalsize);
extern ssize_t pi_file_read(pi_file_handle_ptr f, void * buf, size_t size, size_t totalsize);
extern ssize_t pi_file_write_chunk(pi_file_handle_ptr f, void * buf, size_t size, size_t totalsize, size_t chunksize, size_t spos, size_t epos);
extern ssize_t pi_file_read_chunk(pi_file_handle_ptr f, void * buf, size_t size, size_t totalsize, size_t chunksize, size_t spos, size_t epos);

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

/* pointers must be different on all threads */
extern void pi_thread_bcast(void * sendbuf,
        size_t count, pi_datatype_ptr datatype,
        unsigned int root,
        pi_comm_ptr wr);
extern void pi_abort(int err, pi_comm_ptr wr);
extern void pi_thread_allreduce(void * sendbuf, void * recvbuf,
        size_t count, pi_datatype_ptr datatype, pi_op_ptr op,
        pi_comm_ptr wr);
extern int pi_thread_data_eq(void * buffer,
        size_t count, pi_datatype_ptr datatype,
        pi_comm_ptr wr);
extern void pi_bcast(void * sendbuf,
        size_t count, pi_datatype_ptr datatype,
        unsigned int jroot, unsigned int troot,
        pi_comm_ptr wr);
extern void pi_allreduce(void * sendbuf, void *recvbuf,
        size_t count, pi_datatype_ptr datatype, pi_op_ptr op,
        pi_comm_ptr wr);
extern void pi_allgather(void * sendbuf,
        size_t sendcount, pi_datatype_ptr sendtype,
        void *recvbuf,
        size_t recvcount, pi_datatype_ptr recvtype,
        pi_comm_ptr wr);
extern int pi_data_eq(void * buffer,
        size_t count, pi_datatype_ptr datatype,
        pi_comm_ptr wr);

/* These two interfaces are experimental only */

namespace parallelizing_info_experimental {
    void broadcast(std::vector<unsigned int>& v, parallelizing_info_ptr pi);
    void allgather(std::vector<unsigned int>& v, pi_comm_ptr wr);
    void broadcast(std::set<unsigned int>& v, parallelizing_info_ptr pi);
    void allgather(std::set<unsigned int>& v, pi_comm_ptr wr);
}

/* prints the given string in a ascii-form matrix. */
extern void grid_print(parallelizing_info_ptr pi, char * buf, size_t siz, int print);

#define serialize(w)   serialize__(w, __FILE__, __LINE__)
extern int serialize__(pi_comm_ptr, const char *, unsigned int);
#define serialize_threads(w)   serialize_threads__(w, __FILE__, __LINE__)
extern int serialize_threads__(pi_comm_ptr, const char *, unsigned int);

/* stuff related to log entry printing */
extern void pi_log_init(pi_comm_ptr);
extern void pi_log_clear(pi_comm_ptr);
extern void pi_log_op(pi_comm_ptr, const char * fmt, ...);
extern void pi_log_print_all(parallelizing_info_ptr);
extern void pi_log_print(pi_comm_ptr);

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
extern void pi_interleaving_flip(parallelizing_info_ptr);
extern void pi_interleaving_enter(parallelizing_info_ptr);
extern void pi_interleaving_leave(parallelizing_info_ptr);

extern void pi_store_generic(parallelizing_info_ptr, unsigned long, unsigned long, void *);
extern void * pi_load_generic(parallelizing_info_ptr, unsigned long, unsigned long);

/* shared_malloc is like malloc, except that the pointer returned will be
 * equal on all threads (proper access will deserve proper locking of
 * course). shared_malloc_set_zero sets to zero too */
/* As a side-effect, all shared_* functions serialize threads */
extern void * shared_malloc(pi_comm_ptr wr, size_t size);
extern void * shared_malloc_set_zero(pi_comm_ptr wr, size_t size);
extern void shared_free(pi_comm_ptr wr, void * ptr);

/* Use in std::unique_ptr<T, shared_free_deleter<T>>
 *     or std::unique_ptr<T[], shared_free_deleter<T>>
 */
template<typename T>
struct shared_free_deleter {
    pi_comm_ptr wr = nullptr;
    shared_free_deleter() = default;
    ~shared_free_deleter() = default;
    // NOLINTNEXTLINE(hicpp-explicit-conversions)
    shared_free_deleter(pi_comm_ptr wr) : wr(wr) {}
    shared_free_deleter(shared_free_deleter const &) = default;
    shared_free_deleter& operator=(shared_free_deleter const &) = default;
    shared_free_deleter(shared_free_deleter &&) noexcept = default;
    shared_free_deleter& operator=(shared_free_deleter &&) noexcept = default;
    void operator()(T * ptr) const { if (wr) shared_free(wr, ptr); }
};

template<typename T>
T * shared_new(pi_comm_ptr wr, size_t size)
{
    void * ptr = nullptr;
    if (wr->trank == 0) ptr = new T[size];
    pi_thread_bcast(&ptr, sizeof(void*), BWC_PI_BYTE, 0, wr);
    return static_cast<T *>(ptr);
}

template<typename T>
void shared_delete(pi_comm_ptr wr, T * ptr, size_t)
{
    serialize_threads(wr);
    if (wr->trank == 0) delete[] ptr;
}

template<typename T>
T * shared_new(pi_comm_ptr wr)
{
    void * ptr = nullptr;
    if (wr->trank == 0) ptr = new T;
    pi_thread_bcast(&ptr, sizeof(void*), BWC_PI_BYTE, 0, wr);
    return static_cast<T *>(ptr);
}

template<typename T>
void shared_delete(pi_comm_ptr wr, T * ptr)
{
    serialize_threads(wr);
    if (wr->trank == 0) delete ptr;
}

template<typename T>
struct shared_delete_deleter {
    pi_comm_ptr wr = nullptr;
    shared_delete_deleter() = default;
    ~shared_delete_deleter() = default;
    // NOLINTNEXTLINE(hicpp-explicit-conversions)
    shared_delete_deleter(pi_comm_ptr wr) : wr(wr) {}
    shared_delete_deleter(shared_delete_deleter const &) = default;
    shared_delete_deleter& operator=(shared_delete_deleter const &) = default;
    shared_delete_deleter(shared_delete_deleter &&) noexcept = default;
    shared_delete_deleter& operator=(shared_delete_deleter &&) noexcept = default;
    void operator()(T * ptr) const { if (wr) shared_delete(wr, ptr); }
};

template<typename T>
struct shared_delete_deleter<T[]> {
    pi_comm_ptr wr = nullptr;
    shared_delete_deleter() = default;
    ~shared_delete_deleter() = default;
    // NOLINTNEXTLINE(hicpp-explicit-conversions)
    shared_delete_deleter(pi_comm_ptr wr) : wr(wr) {}
    shared_delete_deleter(shared_delete_deleter const &) = default;
    shared_delete_deleter& operator=(shared_delete_deleter const &) = default;
    shared_delete_deleter(shared_delete_deleter &&) noexcept = default;
    shared_delete_deleter& operator=(shared_delete_deleter &&) noexcept = default;
    void operator()(T * ptr) const { if (wr) shared_delete(wr, ptr, 0); }
};

template<typename T>
struct pi_shared_array : public std::unique_ptr<T[], shared_delete_deleter<T[]>> {
    typedef std::unique_ptr<T[], shared_delete_deleter<T[]>> super;
    explicit pi_shared_array(pi_comm_ptr wr, unsigned int n)
        : super(shared_new<T>(wr, n), wr)
    {}
};
template<typename T>
struct pi_shared_object : public std::unique_ptr<T, shared_delete_deleter<T>> {
    typedef std::unique_ptr<T, shared_delete_deleter<T>> super;
    explicit pi_shared_object(pi_comm_ptr wr)
        : super(shared_new<T>(wr), wr)
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
    for(unsigned int t__ = 0 ; t__ < comm->ncores ; t__++) {		\
        serialize_threads(comm);					\
        if (t__ != comm->trank) continue; /* not our turn. */           \
        do
#define SEVERAL_THREADS_PLAY_MPI_END()    while (0); } } while (0)
/* This construct is used similarly. It differs slightly, in that we
 * guarantee that only one thread (in the communicator) will issue mpi
 * calls */
#define SEVERAL_THREADS_PLAY_MPI_BEGIN2(comm, t__) do { 		\
    serialize_threads(comm);                                            \
    if (comm->trank == 0) {                                             \
        for(unsigned int t__ = 0 ; t__ < comm->ncores ; t__++) {	\
            do
#define SEVERAL_THREADS_PLAY_MPI_END2(comm)                             \
            while (0);							\
        }								\
    }									\
    serialize_threads(comm);                                            \
} while (0)
#else
#define SEVERAL_THREADS_PLAY_MPI_BEGIN(comm)     /**/
#define SEVERAL_THREADS_PLAY_MPI_END()           /**/
#define SEVERAL_THREADS_PLAY_MPI_BEGIN2(comm, t) /**/
#define SEVERAL_THREADS_PLAY_MPI_END2(comm)      /**/
#endif

#endif	/* PARALLELIZING_INFO_HPP_ */
