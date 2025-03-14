#include "cado.h" // IWYU pragma: keep
#include <omp.h>
#include <stdio.h>

int main()
{
    printf("%d\n", omp_get_max_threads());
    /*
    omp_set_dynamic(1);
#pragma omp parallel
    {
        if (omp_get_thread_num() == 0)
        printf("%d\n", omp_get_num_threads());
    }
    */
    return 0;
}

