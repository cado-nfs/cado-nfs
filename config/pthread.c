#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

void * t(void * p)
{
    return p;
}
int main(int argc, char * argv[])
{
    int n = 2;
    if (argc > 1) {
        n = atoi(argv[1]);
    }
    pthread_t * tid = malloc(n * sizeof(pthread_t));
    void ** r = malloc(n * sizeof(void *));
    int rc;
    for(int i = 0 ; i < n ; i++)
        rc = pthread_create(&tid[i], NULL, t, NULL);
    for(int i = 0 ; i < n ; i++)
        pthread_join(tid[i], &r[i]);
    return r[0] == r[1] ? EXIT_SUCCESS : EXIT_FAILURE;
}

