#include "cado.h" // IWYU pragma: keep
#include <string.h>
#include <pthread.h>
#include <unistd.h>     // nanosleep is in posix 2001
#include <stdio.h>
#include "tests_common.h"
#include "verbose.h"
#include "portability.h" // sleep // IWYU pragma: keep

int verbose;
volatile size_t conflict;

void *print_stuff (void *data)
{
  const char *text = data;
  const size_t len = strlen(text);
  verbose_output_start_batch();
  for (size_t i = 0; i < len; i++) {
    verbose_output_print(1, verbose, "%c", text[i]);
    conflict = i;
    fflush(stdout);
    struct timespec tv = { 0, 5*1000*1000 };
    struct timespec rem;
    nanosleep(&tv, &rem);
  }
  verbose_output_end_batch();
  return NULL;
}


int main(int argc, const char **argv)
{
  tests_common_cmdline(&argc, &argv, PARSE_VERBOSE);
  verbose = tests_common_get_verbose();

  /* Default outputs */
  verbose_output_print(0, verbose, "1 stdout\n");
  verbose_output_print(1, verbose, "2 stderr\n");

  /* Define a few outputs */
  verbose_output_init(3);
  verbose_output_add(0, stderr, 1);
  verbose_output_add(1, stdout, 2);
  verbose_output_add(2, stdout, 3);
  verbose_output_add(2, stderr, 4);
  /* 0=stderr, 1=stdout, 2=stdout+stderr now */

  verbose_output_print(0, verbose, "3 ch0");
  verbose_output_print(1, verbose, "4 ch1");
  verbose_output_print(2, verbose, "5 ch2");

  int nr_lines = 5;
  for (size_t channel = 0; channel < 3; channel++) {
    FILE *out;
    for (size_t i=0; (out=verbose_output_get(channel, verbose, i)) != NULL; i++)
        fprintf(out, "%d ch%zu out%zu\n", ++nr_lines, channel, i);
  }

  pthread_t threads[2];
  pthread_create(&threads[0], NULL, print_stuff, "the quick brown fox\n");
  pthread_create(&threads[1], NULL, print_stuff, "jumps over the lazy dog\n");
  pthread_join(threads[0], NULL);
  pthread_join(threads[1], NULL);

  /* Back to default channels */
  verbose_output_clear();
  verbose_output_print(0, verbose, "%d stdout\n", ++nr_lines);
  verbose_output_print(1, verbose, "%d stderr\n", ++nr_lines);

  tests_common_clear();
}
