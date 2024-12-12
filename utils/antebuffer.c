/* Very efficient antebuffer to preempt the file(s) load.
   Syntax: antebuffer X file [file] | YourCommand
   with buffer size = 2^X.
   The files are written on stdout.
   Best size for the buffer is about 4MB (local disk) to
   16MB (NFS), eventually until 128 MB.
   Limitation: X <= 31; more than 2GB antebuffer has no sense anyway.
*/

/* To avoid the warning: implicit declaration of nanosleep for c99 compliant */
#include "cado.h" // IWYU pragma: keep

// I don't want to look into this code, really.
// NOLINTBEGIN

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdint.h>
#include <pthread.h>
#include <errno.h>
#include <time.h>
#include "macros.h"     // LIKELY UNLIKELY
#include "timing.h"
#include "portability.h" // sleep // IWYU pragma: keep

#ifdef HAVE_MINGW
int _CRT_fmode = _O_BINARY; /* Binary open for stdin/out/err */
#endif

#ifndef HAVE_NANOSLEEP
  int nanosleep(const struct timespec *req, struct timespec *rem) {
    if (rem == NULL) {
      /* Dummy to shut up the warning */
    }
#ifdef HAVE_USLEEP
    unsigned long usec = req->tv_sec * 1000000UL + req->tv_nsec / 1000UL;
    usleep(usec);
#else
    sleep(req->tv_sec);
#endif
    return 0;
  }
#endif

/* These variables are used by pthread and main */
static volatile uintptr_t ab_cptp = 0, ab_cptc = 0;        /* Main counters */
static volatile int ab_end = 0;                            /* 1 if end */
static char *ab_buf;                                       /* Main buffer */
static const struct timespec waiting = { 0, 1<<13 };       /* About 8 to 20 microseconds */
static int ab_in;                                          /* fd for loading */
static size_t ab_size, ab_sizeio;                          /* size for buffer & in */

/* The 3 shared variables must be considered with non atomical access.
   So, dedicated mutex protection for all.
   NB: spinlocks are not mandatory in POSIX and don't exist in macos.
   But it's so uncommon a mutex is blocking here than the mutex version
   has the same speed than the spin locks version, and both are faster
   than the gcc __builtin_fetch_and_add version.
*/
static pthread_mutex_t mutex_ab_cptp = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mutex_ab_cptc = PTHREAD_MUTEX_INITIALIZER;
static pthread_mutex_t mutex_ab_end = PTHREAD_MUTEX_INITIALIZER;

static void ab_cons () {
  uintptr_t cpy_ab_cptc = 0;
  size_t c, t = 0;
  ssize_t w;
  
  pthread_setcanceltype (PTHREAD_CANCEL_DEFERRED, NULL);
  for ( ; ; ) {
    for ( ; ; ) {
      uintptr_t mut_ab_cptp;
      int mut_ab_end;
      pthread_mutex_lock (&mutex_ab_cptp); mut_ab_cptp = ab_cptp; pthread_mutex_unlock (&mutex_ab_cptp);
      c = (mut_ab_cptp - cpy_ab_cptc);
      if (LIKELY (c)) break;
      pthread_mutex_lock (&mutex_ab_end); mut_ab_end = ab_end; pthread_mutex_unlock (&mutex_ab_end);
      if (UNLIKELY (mut_ab_end)) {
	pthread_mutex_lock (&mutex_ab_cptp); mut_ab_cptp = ab_cptp; pthread_mutex_unlock (&mutex_ab_cptp);
	c = (int) (mut_ab_cptp - cpy_ab_cptc);
	if (LIKELY (c)) break;
	pthread_exit(NULL);
      }
      nanosleep (&waiting, NULL);
    }
    if (c > ab_sizeio) c = ab_sizeio;
    if (c > ab_size - t) c = ab_size - t;
    while ((w = write(1, &(ab_buf[t]), c)) <= 0) nanosleep (&waiting, NULL);
    cpy_ab_cptc += w;
    t = (t + w) & (ab_size - 1);
    pthread_mutex_lock (&mutex_ab_cptc); ab_cptc = cpy_ab_cptc; pthread_mutex_unlock (&mutex_ab_cptc);
  }
}

int main(int argc, char const * argv[])
{
  pthread_t ab_tc;
  pthread_attr_t ab_attr;
  uintptr_t cpy_ab_cptp = 0;
  size_t c, t = 0;
  ssize_t r;
  unsigned int p;
  char *real_malloc;

#ifdef HAVE_MINGW
  _fmode = _O_BINARY;     /* Binary open for all others files */
#endif

  if (argc < 3) {
  error:
    fprintf (stderr, "%s syntax:  %s SIZE file [file]  or  %s SIZE -\n"
	     "This command is an ante buffer which loads the files or stdin in a\n"
	     "buffer and write them on stdout. The size of this buffer is 2^SIZE.\n"
	     "16<=SIZE<=31, best for 20 to 24.\n", argv[0], argv[0], argv[0]);
    exit (1);
  }
  if ((sscanf(argv[1], "%zu", &ab_size)) != 1) goto error;
  if (ab_size < 16 || ab_size > 31) goto error;
  ab_size = (uintptr_t) 1 << ab_size;
  ab_sizeio = ab_size >> 4;
  if (!(real_malloc = malloc (ab_size + 0xFFF))) {
    fprintf (stderr, "%s: malloc error: %s\n", argv[0], strerror(errno));
    exit(1);
  }
  ab_buf = (char *) (((uintptr_t) real_malloc + 0xFFF) & ~ (uintptr_t) 0xFFF);
  pthread_attr_init(&ab_attr);
  pthread_attr_setstacksize(&ab_attr, 1<<12);
  pthread_attr_setdetachstate(&ab_attr, PTHREAD_CREATE_JOINABLE);
  if (pthread_create(&ab_tc, &ab_attr, (void *) ab_cons, NULL)) {
    fprintf (stderr, "%s: pthread_create error: %s\n", argv[0], strerror(errno));
    exit(1);
  }
  if (strcmp(argv[2], "-")) {
    p = 3;
    if ((ab_in = open(argv[2], O_RDONLY)) == -1) {
      fprintf (stderr, "%s: open or load error in file %s: %s\n", argv[0], argv[2], strerror(errno));
      exit (1);
    }
  }
  else {
    p = 0;
    ab_in = 0;
  }
  for ( ; ; ) {
    for ( ; ; ) {
      uintptr_t mut_ab_cptc;
      pthread_mutex_lock (&mutex_ab_cptc); mut_ab_cptc = ab_cptc; pthread_mutex_unlock (&mutex_ab_cptc);
      c = ((mut_ab_cptc + ab_size) - cpy_ab_cptp);
      if (LIKELY (c)) break;
      nanosleep (&waiting, NULL);
    }
    if (c > ab_sizeio)   c = ab_sizeio;
    if (c > ab_size - t) c = ab_size - t;
    r = read(ab_in, &(ab_buf[t]), c);
    if (r > 0) {
      cpy_ab_cptp += r;
      t = (t + r) & (ab_size - 1);
      pthread_mutex_lock (&mutex_ab_cptp); ab_cptp = cpy_ab_cptp; pthread_mutex_unlock (&mutex_ab_cptp);
    }
    else 
      if (!r) {
	if (!p || p == (unsigned int) argc) {
	  pthread_mutex_lock (&mutex_ab_end); ab_end = 1; pthread_mutex_unlock (&mutex_ab_end);
	  break;
	}
	close(ab_in);
	if ((ab_in = open(argv[p++], O_RDONLY)) == -1) {
	  fprintf (stderr, "%s: open or load error in file %s: %s\n", argv[0], argv[p-1], strerror(errno));
	  exit (1);
	}
      }
      else {
	fprintf (stderr, "%s: read error in file %s: %s\n", argv[0], argv[p-1], strerror(errno));
	  exit (1);
	}
  }
  if (p) close(ab_in);
  pthread_join(ab_tc, NULL);
  pthread_mutex_destroy (&mutex_ab_end);
  pthread_mutex_destroy (&mutex_ab_cptc);
  pthread_mutex_destroy (&mutex_ab_cptp);
  free (real_malloc);
  double tt[2];
  seconds_user_sys(tt);
  /*
  fprintf(stderr, "antebuffer exits after having spent %.2fs+%.2fs on cpu\n",
          tt[0], tt[1]);
          */
  exit (0);
}
// NOLINTEND
