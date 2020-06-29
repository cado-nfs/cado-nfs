#define _POSIX_C_SOURCE 200809L
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

pthread_barrier_t barrier;

void * t(void * p)
{
    pthread_barrier_wait(&barrier);
    return p;
}
int main(int argc, char * argv[])
{
    int n = 2;
    if (argc > 1) {
        n = atoi(argv[1]);
    }
    pthread_barrier_init(&barrier, NULL, n);
    pthread_t * tid = malloc(n * sizeof(pthread_t));
    void ** r = malloc(n * sizeof(void *));
    int rc;
    for(int i = 0 ; i < n ; i++)
        rc = pthread_create(&tid[i], NULL, t, NULL);
    for(int i = 0 ; i < n ; i++)
        pthread_join(tid[i], &r[i]);
    pthread_barrier_destroy(&barrier);
    return EXIT_SUCCESS;
}

