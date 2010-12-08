#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <time.h>
#include <assert.h>

pthread_mutex_t mx;
pthread_barrier_t barrier;

int width;
int n;
unsigned long * * data;;

void * program(void * arg)
{
    int t = ((unsigned long **)arg) - data;
    unsigned long * ptr = * (unsigned long **)arg;
    printf("Thread %d says hi\n", t);

    for(int i = 0 ; i < width ; i++) {
        ptr[i] = t;
    }

    clock_t t0;
    double dt;
    int method=0;

    /* Method 1 */

#define HEADER() do {                   				\
        method++;                           				\
        pthread_barrier_wait(&barrier);     				\
        t0 = clock();                       				\
    } while (0)
#define TRAILER() do {							\
        pthread_barrier_wait(&barrier);					\
        dt = (double) (clock() - t0) / CLOCKS_PER_SEC;			\
        if (t == 0) {							\
            printf("Method %d: Took %.2fs\n", method, dt);		\
        }								\
    } while (0)



    HEADER();
    if (t == 0) {
        for(int j = 1 ; j < n ; j++) {
            for(int k = 0 ; k < width ; k++) {
                data[0][k] ^= data[j][k];
            }
        }
    }
    TRAILER();

    HEADER();
    int k0 = t * width / n;
    int k1 = (t+1) * width / n;
    for(int k = k0 ; k < k1 ; k++) {
        for(int j = 1 ; j < n ; j++) {
            data[0][k] ^= data[j][k];
        }
    }
    TRAILER();

    return NULL;
}



int main(int argc, char * argv[])
{
    if (argc != 3) {
        fprintf(stderr, "Usage: ./test-pthread-reduce <nthreads> <datasize>\n");
        exit(1);
    }

    n = atoi(argv[1]);
    width = atoi(argv[2]);
    pthread_barrier_init(&barrier, NULL, n);
    pthread_mutex_init(&mx, NULL);

    pthread_t * threads = malloc(n * sizeof(pthread_t));
    data = malloc(n * sizeof(unsigned long *));

    for(int i = 0 ; i < n ; i++) {
        data[i] = malloc(width * sizeof(unsigned long));
        pthread_create(&(threads[i]), NULL, &program, data + i);
    }

    for(int i = 0 ; i < n ; i++) {
        pthread_join(threads[i], NULL);
        free(data[i]);
    }
    free(data);
    free(threads);
    pthread_mutex_destroy(&mx);
    pthread_barrier_destroy(&barrier);
}
