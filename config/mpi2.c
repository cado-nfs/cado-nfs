#include <mpi.h>
#if !((MPI_VERSION > 2) || (MPI_VERSION == 2 && MPI_SUBVERSION >= 1))
#error "MPI version 2.1 not supported"
#endif
int main(int argc, char** argv)
{
    MPI_Init(&argc,&argv);
    MPI_Finalize();
}
