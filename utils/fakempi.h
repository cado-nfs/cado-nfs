#ifndef	FAKEMPI_H_
#define	FAKEMPI_H_

#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "cado_mpi_config.h"
#include "macros.h"
#include "portability.h"
#include "misc.h"

typedef int MPI_Status;
typedef int MPI_Datatype;
typedef int MPI_Comm;
typedef int MPI_Op;
typedef void MPI_User_function(void *invec, void *inoutvec,
             int *len, MPI_Datatype *datatype);
typedef int MPI_Errhandler;
typedef int MPI_Request;

// type keys are made as follows:
// lowest 8 bits: sizeof()
// high 8 bits: a uniqueifying tag
#define MPI_DATATYPE_NULL       0x0000
#define MPI_BYTE                0x0101
#define MPI_INT                 (0x0200 | sizeof(int))
#define MPI_DOUBLE              (0x0300 | sizeof(double))
#define MPI_UNSIGNED_LONG       (0x0400 | sizeof(unsigned long))
#define MPI_UNSIGNED_LONG_LONG       (0x0500 | sizeof(unsigned long long))
#define MPI_LONG                (0x0600 | sizeof(long))
#define MPI_LONG_LONG           (0x0700 | sizeof(long long))
/* It seems that MPI_UNSIGNED_INT is in fact unspecified */
#define MPI_UNSIGNED            (0x0800 | sizeof(unsigned int))

#define fakempi_sizeof_type(x) ((x) & 0xff)

#define MPI_Type_size(x, s)     *(s)=fakempi_sizeof_type(x)

#define MPI_COMM_WORLD	0

#define MPI_THREAD_SINGLE       0
#define MPI_THREAD_FUNNELED     1
#define MPI_THREAD_SERIALIZED   2
#define MPI_THREAD_MULTIPLE     3

/* We define different ops, but since they're collected amongst
 * communicators of size 1 anyway, the operation does not matter much...
 * */
#define MPI_BXOR       0
#define MPI_SUM        1
#define MPI_MAX        2
#define MPI_MIN        3
#define MPI_LAND       4
#define MPI_BAND       5
#define MPI_BOR        6

#define MPI_ERRORS_ARE_FATAl        0
#define MPI_ERRORS_RETURN        1

#define MPI_IN_PLACE    0

#define MPI_MAX_OBJECT_NAME     64

#define MPI_STATUS_IGNORE       0
#define MPI_STATUSES_IGNORE       0

#define MPI_ANY_TAG     -1

/* Adding this pragma would yield the benefit of removing all the
 * MAYBE_UNUSED clutter which is inherent to this file. Unfortunately
 * it's not supported by very old gcc's, and neither by some broken stuff
 * which pretends to be gcc (happens on macs, for instance).
 *
#ifdef  __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif
*/

static inline int MPI_Wait(MPI_Request *request MAYBE_UNUSED, MPI_Status *status MAYBE_UNUSED) { return 0; }
static inline int MPI_Waitall(int count MAYBE_UNUSED, MPI_Request *request MAYBE_UNUSED, MPI_Status *statuses MAYBE_UNUSED) { return 0; }
static inline int MPI_Abort(MPI_Comm comm MAYBE_UNUSED, int s) { exit(s); }
static inline int MPI_Comm_rank(int s MAYBE_UNUSED, int  * p) { *p=0; return 0;}
static inline int MPI_Comm_size(int s MAYBE_UNUSED, int  * p) { *p=1; return 0;}
static inline int MPI_Initialized(int  * p) { *p=1; return 0; }
static inline int MPI_Init(int * argc MAYBE_UNUSED, char *** argv MAYBE_UNUSED) { return 0; }
static inline int MPI_Init_thread(int * argc MAYBE_UNUSED, char *** argv MAYBE_UNUSED, int req, int * prov) { if (prov) *prov=req; return 0; }
static inline int MPI_Finalize() {return 0;}
static inline int MPI_Op_create( MPI_User_function *function MAYBE_UNUSED, int commute MAYBE_UNUSED, MPI_Op *op MAYBE_UNUSED ){return 0;}
static inline int MPI_Op_free(MPI_Op *op MAYBE_UNUSED ){return 0;}
static inline int MPI_Send( void *buf MAYBE_UNUSED, int count MAYBE_UNUSED, MPI_Datatype datatype MAYBE_UNUSED, int dest MAYBE_UNUSED,int tag MAYBE_UNUSED, MPI_Comm comm  MAYBE_UNUSED){return 0;}
static inline int MPI_Sendrecv( void *sbuf MAYBE_UNUSED, int scount MAYBE_UNUSED, MPI_Datatype sdatatype MAYBE_UNUSED, int sdest MAYBE_UNUSED,int stag MAYBE_UNUSED,  void *rbuf MAYBE_UNUSED, int rcount MAYBE_UNUSED, MPI_Datatype rdatatype MAYBE_UNUSED, int rdest MAYBE_UNUSED,int rtag MAYBE_UNUSED, MPI_Comm comm  MAYBE_UNUSED, MPI_Status *status MAYBE_UNUSED){return 0;}
static inline int MPI_Isend( void *buf MAYBE_UNUSED, int count MAYBE_UNUSED, MPI_Datatype datatype MAYBE_UNUSED, int dest MAYBE_UNUSED,int tag MAYBE_UNUSED, MPI_Comm comm  MAYBE_UNUSED, MPI_Request * zz MAYBE_UNUSED){return 0;}
static inline int MPI_Recv( void *buf MAYBE_UNUSED, int count MAYBE_UNUSED, MPI_Datatype datatype MAYBE_UNUSED, int source MAYBE_UNUSED,int tag MAYBE_UNUSED, MPI_Comm comm MAYBE_UNUSED, MPI_Status *status MAYBE_UNUSED ){ abort(); return 0;}
static inline int MPI_Irecv( void *buf MAYBE_UNUSED, int count MAYBE_UNUSED, MPI_Datatype datatype MAYBE_UNUSED, int source MAYBE_UNUSED,int tag MAYBE_UNUSED, MPI_Comm comm MAYBE_UNUSED, MPI_Request * zz MAYBE_UNUSED){ abort(); return 0;}
static inline int MPI_Bcast( void *buffer MAYBE_UNUSED, int count MAYBE_UNUSED, MPI_Datatype datatype MAYBE_UNUSED, int root MAYBE_UNUSED, MPI_Comm comm MAYBE_UNUSED){return 0;}
static inline int MPI_Reduce ( void *sendbuf, void *recvbuf, int count,MPI_Datatype datatype, MPI_Op op MAYBE_UNUSED, int root MAYBE_UNUSED, MPI_Comm comm  MAYBE_UNUSED)
{
    if (sendbuf) memcpy(recvbuf, sendbuf, count * fakempi_sizeof_type(datatype));
    return 0;
}

static inline int MPI_Reduce_scatter(void *sendbuf, void *recvbuf, int *recvcounts,
                    MPI_Datatype datatype, MPI_Op op MAYBE_UNUSED, MPI_Comm comm MAYBE_UNUSED)
{
    if (sendbuf) memcpy(recvbuf, sendbuf, recvcounts[0] * fakempi_sizeof_type(datatype));
    return 0;
}
static inline int MPI_Allreduce ( void *sendbuf, void *recvbuf, int count,MPI_Datatype datatype, MPI_Op op MAYBE_UNUSED, MPI_Comm comm  MAYBE_UNUSED)
{
    if (sendbuf) memcpy(recvbuf, sendbuf, count * fakempi_sizeof_type(datatype));
    return 0;
}
static inline int MPI_Comm_split (MPI_Comm x MAYBE_UNUSED, int color MAYBE_UNUSED, int key MAYBE_UNUSED, MPI_Comm * y)
{
    *y=0;
    return 0;
}
static inline int MPI_Comm_set_errhandler (MPI_Comm x MAYBE_UNUSED, MPI_Errhandler e MAYBE_UNUSED)
{
    return 0;
}
static inline int MPI_Comm_free (MPI_Comm * x MAYBE_UNUSED) { return 0; }
static inline int MPI_Comm_dup (MPI_Comm y, MPI_Comm * x) { *x = y; return 0; }
static inline int MPI_Comm_set_name(MPI_Comm comm MAYBE_UNUSED, char *comm_name MAYBE_UNUSED) { return 0;}
static inline int MPI_Comm_get_name(MPI_Comm comm MAYBE_UNUSED, char *comm_name MAYBE_UNUSED, int * rlen) { *comm_name='\0'; *rlen=0; return 0;}
static inline int MPI_Scatterv(void * sendbuf, int * sendcounts, int * displs,  MPI_Datatype st, void * recvbuf, int recvcount, MPI_Datatype rt, int root MAYBE_UNUSED, MPI_Comm x MAYBE_UNUSED) {
    ASSERT_ALWAYS(sendcounts[0] * fakempi_sizeof_type(st) == recvcount * fakempi_sizeof_type(rt));
    memcpy(recvbuf, ((char *)sendbuf) + displs[0] * fakempi_sizeof_type(st), recvcount * fakempi_sizeof_type(rt));
    return 0;
}

static inline int MPI_Scatter(void * sendbuf, int sendcount, MPI_Datatype st, void * recvbuf, int recvcount, MPI_Datatype rt, int root MAYBE_UNUSED, MPI_Comm x MAYBE_UNUSED) {
    ASSERT_ALWAYS(sendcount * fakempi_sizeof_type(st) == recvcount * fakempi_sizeof_type(rt));
    if (recvbuf && sendbuf)
        memcpy(recvbuf, sendbuf, recvcount * fakempi_sizeof_type(rt));
    return 0;
}

static inline int MPI_Barrier (MPI_Comm x MAYBE_UNUSED) { return 0; }

static inline int MPI_Gather(void * sendbuf, int sendcount,  MPI_Datatype st, void * recvbuf, int recvcount MAYBE_UNUSED, MPI_Datatype rt MAYBE_UNUSED, int root MAYBE_UNUSED, MPI_Comm x MAYBE_UNUSED) {
    if (sendbuf == MPI_IN_PLACE) return 0;
    memcpy(((char *)recvbuf), (char*) sendbuf, sendcount * fakempi_sizeof_type(st));
    return 0;
}
static inline int MPI_Gatherv(void * sendbuf, int sendcount,  MPI_Datatype st, void * recvbuf, int * recvcounts, int * displs, MPI_Datatype rt, int root MAYBE_UNUSED, MPI_Comm x MAYBE_UNUSED) {
    ASSERT_ALWAYS(sendcount * fakempi_sizeof_type(st) == recvcounts[0] * fakempi_sizeof_type(rt));
    memcpy(((char *)recvbuf) + displs[0] * fakempi_sizeof_type(rt), sendbuf, sendcount * fakempi_sizeof_type(st));
    return 0;
}
static inline int MPI_Allgather(void * sendbuf MAYBE_UNUSED, int sendcount MAYBE_UNUSED,  MPI_Datatype st MAYBE_UNUSED, void * recvbuf MAYBE_UNUSED, int recvcount MAYBE_UNUSED, MPI_Datatype rt MAYBE_UNUSED, MPI_Comm x MAYBE_UNUSED) {
    ASSERT_ALWAYS(sendbuf == MPI_IN_PLACE || sendcount * fakempi_sizeof_type(st) == recvcount * fakempi_sizeof_type(rt));
    if (sendbuf) memcpy(recvbuf, sendbuf, sendcount * fakempi_sizeof_type(st));
    return 0;
}

static inline int MPI_Iallgather(void *sendbuf, int  sendcount, MPI_Datatype st, void *recvbuf, int recvcount, MPI_Datatype rt, MPI_Comm comm, MPI_Request *request MAYBE_UNUSED)
{
    return MPI_Allgather(sendbuf, sendcount, st, recvbuf, recvcount, rt, comm);
}

static inline int MPI_Allgatherv(void *sendbuf, int sendcount MAYBE_UNUSED,
            MPI_Datatype sendtype MAYBE_UNUSED, void *recvbuf MAYBE_UNUSED, int *recvcount MAYBE_UNUSED,
            int *displs MAYBE_UNUSED, MPI_Datatype recvtype MAYBE_UNUSED, MPI_Comm comm MAYBE_UNUSED)
{
    ASSERT_ALWAYS(sendbuf == MPI_IN_PLACE);
    return 0;
}

static inline int MPI_Alltoall(const void *sendbuf, int sendcount,
            MPI_Datatype sendtype, void *recvbuf, int recvcount,
            MPI_Datatype recvtype, MPI_Comm comm MAYBE_UNUSED)
{
    if (sendbuf == MPI_IN_PLACE) return 0;
    ASSERT_ALWAYS(sendcount * fakempi_sizeof_type(sendtype) == recvcount * fakempi_sizeof_type(recvtype));
    memcpy(recvbuf, sendbuf, recvcount * fakempi_sizeof_type(recvtype));
    return 0;
}

static inline int MPI_Alltoallv(const void *sendbuf, const int sendcounts[],
            const int sdispls[], MPI_Datatype sendtype,
            void *recvbuf, const int recvcounts[],
            const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm MAYBE_UNUSED)
{
    if (sendbuf == MPI_IN_PLACE) return 0;
    ASSERT_ALWAYS(sendcounts[0] * fakempi_sizeof_type(sendtype) == recvcounts[0] * fakempi_sizeof_type(recvtype));
    memcpy(((char*)recvbuf) + rdispls[0] * fakempi_sizeof_type(recvtype),
            ((const char*)sendbuf) + sdispls[0] * fakempi_sizeof_type(sendtype), 
            sendcounts[0] * fakempi_sizeof_type(sendtype));
    return 0;
}

static inline int MPI_Type_contiguous(int count, MPI_Datatype oldtype, MPI_Datatype * newtype)
{
    *newtype = count * oldtype;
    return 0;
}
static inline int MPI_Type_commit(MPI_Datatype * t MAYBE_UNUSED) { return 0; }
static inline int MPI_Type_free(MPI_Datatype * t MAYBE_UNUSED) { return 0; }
static inline int MPI_Type_set_attr(MPI_Datatype type MAYBE_UNUSED, int key MAYBE_UNUSED, void *value MAYBE_UNUSED)
{
    /* XXX We are *NOT* storing the type attribute here. This is because
     * we expect that the user function that we use, and which exploits
     * this, will actually *never* be called in this context of ``fake''
     * mpi: by assumption, all communicator sizes are equal to 1, and
     * thus mpi reduction just amounts to copying data.
     */
    return 0;
}
static inline int MPI_Type_get_attr(MPI_Datatype type MAYBE_UNUSED, int key MAYBE_UNUSED, void *value, int * flag)
{
    /* Same as above */
    // *(void**)value = NULL;
    memset(value, 0, sizeof(void*));
    *flag=1;
    return 0;
}
static inline int MPI_Type_delete_attr(MPI_Datatype type MAYBE_UNUSED, int key MAYBE_UNUSED) { return 0; }

typedef int MPI_Type_copy_attr_function(MPI_Datatype oldtype,
           int type_keyval, void *extra_state, void *attribute_val_in,
           void *attribute_val_out, int *flag);
typedef int MPI_Type_delete_attr_function(MPI_Datatype type, int type_keyval,
            void *attribute_val, void *extra_state);
#define MPI_TYPE_DUP_FN NULL
#define MPI_TYPE_NULL_DELETE_FN NULL
static inline int MPI_Type_create_keyval(
        MPI_Type_copy_attr_function *type_copy_attr_fn MAYBE_UNUSED,
        MPI_Type_delete_attr_function *type_delete_attr_fn MAYBE_UNUSED,
        int *type_keyval,
        void *extra_state MAYBE_UNUSED)
{
    /* same rationale as above: we don't care about providing usable
     * function. Only the prototypes are barely right. The rest will
     * actually never be called.
     */
    *type_keyval = 0;
    return 0;
}
static inline int MPI_Type_free_keyval(int * x MAYBE_UNUSED)
{
    return 0;
}

#define MPI_MAX_ERROR_STRING    64
static inline int MPI_Error_string(int err, char * msg, int * len)
{
    *len = snprintf(msg, MPI_MAX_ERROR_STRING, "%s", strerror(err));
    return *len >= 0 ? 0 : ENOMEM;
}

/* We advertise version 2.1, because there's an impossible
 * chicken-and-egg problem with MPI_Reduce_local. We can't implement it
 * without having access to the function doing the reduction. And the
 * latter can be a custom one, accessed with... MPI_Reduce_local.
 */
#define MPI_VERSION 2
#define MPI_SUBVERSION 1

#if 0
/* This function is 3.0 only. For the above reason, we don't expose it */
#define MPI_MAX_LIBRARY_VERSION_STRING  32
static inline int MPI_Get_library_version(char * ptr, int * len)
{
    size_t res = strlcpy(ptr, "fake mpi 0.0", MPI_MAX_LIBRARY_VERSION_STRING);
    ASSERT_ALWAYS(res < MPI_MAX_LIBRARY_VERSION_STRING);
    *len = strlen(ptr);
    return 0;
}
#endif

static inline int MPI_Get_version(int * ver, int * subver)
{
    *ver = MPI_VERSION;
    *subver = MPI_SUBVERSION;
    return 0;
}

/* See at the beginning of this file.
#ifdef  __GNUC__
#pragma GCC diagnostic pop
#endif
*/

#endif /* FAKEMPI_H_ */
