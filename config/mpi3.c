#include <mpi.h>
#if !((MPI_VERSION > 3) || (MPI_VERSION == 3 && MPI_SUBVERSION >= 0))
#error "MPI version 3.0 not supported"
#endif
int main(int argc, char** argv)
{
    MPI_Init(&argc,&argv);
    MPI_Finalize();
}
