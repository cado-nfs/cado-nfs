/*
 * Program: free relations
 * Original author : F. Morain
 * Purpose: creating free relations in a suitable format
 * Modified / rewritten by C. Bouvier (and others)
 * Multi-thread code by A. Filbois

This file is part of CADO-NFS.

CADO-NFS is free software; you can redistribute it and/or modify it under the
terms of the GNU Lesser General Public License as published by the Free
Software Foundation; either version 2.1 of the License, or (at your option)
any later version.

CADO-NFS is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
details.

You should have received a copy of the GNU Lesser General Public License
along with CADO-NFS; see the file COPYING.  If not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* Model of this code : One producer -> Many consumers/producers -> one consumer.

   1. One thread produces NB_PRIMERS_PER_BUFFER primes per buffer, in
      NB_BUFFERS_PER_THREAD buffers (named bufs) for each working thread. This
      thread uses the function pthread_primes_producer.

      When all the primes are produced, this thread produces one buffer with
      only one "false" prime = FREEREL_END_MARKER for each consumer. It's the
      end marker.

   2. 'nthreads' threads load a primes buffer by a classical one producer/many
      consumer models (2 pseudo semaphores BY pair of producer/consumer).
      Each thread produces 2 buffers: one buffer containing the local
      renumbering table on the form of a ASCII array and of one buffer
      containing the information to retrieve the free relations. Those threads
      use the function pthread_primes_consumer.

      When a thread loads a primes buffer which begins and contains only
      FREEREL_END_MARKER, it produces a free relations buffer which contains
      only FREEREL_END_MARKER and exits.

   3. With a many producers/one consumer model, the principal programs (the
      function generate_renumber_and_freerels) loads the data produced by the
      'nthreads' pthread_primes_consumer, writes sequentially the renumbering
      table, and print the free relations.

      When a roots buffer contains exactly one FREEREL_END_MARKER, the job is
      done.
*/


#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>   /* for _O_BINARY */
#include <string.h>
#include <math.h>
#include <pthread.h>

#include <gmp.h>
#include "mod_ul.c"
#include "portability.h"
#include "utils.h"
#include "typedefs.h"

/* 2 buffers are sufficient; in fact it works with only one! */
#define NB_BUFFERS_PER_THREAD (1<<1)
#define NB_PRIMES_PER_BUFFER (1<<10)
#define CACHELINESIZE 64
#define INIT_SIZE_BUF_CHAR (1<<15) /* Init size for the char buffer. */
#if INIT_SIZE_BUF_CHAR < RENUMBER_MAX_SIZE_PER_PRIME
  #error "INIT_SIZE_BUF_CHAR must be greater than RENUMBER_MAX_SIZE_PER_PRIME."
#endif
#define FREEREL_END_MARKER ((unsigned long) (-1))

/************************** Unnamed semaphores ********************************/
/* Unnamed Posix semaphores don't exist in MacOS, named Posix semaphores are
 * persistent and need a interrupt handler to clean them, so unnamed semaphores
 * are implemented here with pthread_cond_t and pthread_mutex_t, and of course a
 * value (unsigned int).
 */
struct sema_s
{
  unsigned int value;
  pthread_mutex_t mutex;
  pthread_cond_t conditional;
};

typedef struct sema_s sema_t[1];
typedef struct sema_s * sema_ptr;

/* These functions are sem_init, sem_post, sem_wait and sem_destroy with mutex &
 * pthread conditional variables.
 * Careful:
 *  1. it's not possible to share these pseudo semaphores between processus, but
 *     only between threads.
 *  2. I use pthread_cond_signal here, not pthread_cond_broadcast! So sema_post
 *     wakes only one sleeper!
 */

void sema_init (sema_ptr sema, unsigned int value_init)
{
  sema->value = value_init;
  if (pthread_mutex_init(&(sema->mutex), NULL))
  {
    perror ("Error in pthread_mutex_init\n");
    exit (1);
  }
  if (pthread_cond_init(&(sema->conditional), NULL))
  {
    perror ("Error in pthread_cond_init\n");
    exit (1);
  }
}

void sema_post (sema_ptr sema)
{
  pthread_mutex_lock (&(sema->mutex));
  if (!(sema->value)++)
    pthread_cond_signal (&(sema->conditional)); /* CAREFUL: no broadcast here! */
  pthread_mutex_unlock (&(sema->mutex));
}

void sema_wait (sema_ptr sema)
{
  pthread_mutex_lock (&(sema->mutex));
  /* The loop is need for the "spurious wakeup", but these is only at most 1
   * sleeper by sema... so, not really useful. */
  while (!(sema->value))
    pthread_cond_wait (&(sema->conditional), &(sema->mutex));
  sema->value--;
  pthread_mutex_unlock (&(sema->mutex));
}

void sema_destroy (sema_ptr sema)
{
  pthread_mutex_destroy (&(sema->mutex));
  pthread_cond_destroy (&(sema->conditional));
}

/*************************** char buffer **************************************/
struct char_buffer_s
{
  char *begin, *cur, *end;
};

typedef struct char_buffer_s char_buffer_t[1];
typedef struct char_buffer_s * char_buffer_ptr;
typedef const struct char_buffer_s * char_buffer_srcptr;

static inline void
char_buffer_init (char_buffer_ptr buf, size_t size)
{
  buf->begin = malloc_aligned (size * sizeof (char), CACHELINESIZE);
  ASSERT_ALWAYS (buf->begin != NULL);
  buf->cur = buf->begin;
  buf->end = buf->begin + size;
}

static inline void
char_buffer_reset (char_buffer_ptr buf)
{
  buf->cur = buf->begin;
}

/* Write the char_buffer in one ASCII block */
static inline void
char_buffer_write (FILE *out, char_buffer_ptr buf)
{
  fwrite ((void *) buf->begin, (buf->cur-buf->begin) * sizeof(char), 1, out);
}

/* This function grows the char * if the remaining space is lower than
 * min_needed.
 */
static inline void
char_buffer_resize (char_buffer_ptr buf, size_t min_needed)
{
  if (UNLIKELY(buf->cur + min_needed >= buf->end))
  {
    size_t save_current = buf->cur - buf->begin,
    new_size = (buf->end - buf->begin) << 1;
    buf->begin = (char *) realloc (buf->begin, new_size * sizeof (char));
    ASSERT_ALWAYS (buf->begin != NULL);
    buf->cur = buf->begin + save_current;
    buf->end = buf->begin + new_size;
  }
}

static inline void
char_buffer_clear (char_buffer_ptr buf)
{
  free (buf->begin);
  buf->end = buf->cur = buf->begin;
}

/************************** Thread buffer *************************************/
struct th_buf_s
{
  /* protected by primes_full/primes_empty */
  unsigned long primes[NB_PRIMES_PER_BUFFER]; /* array of prime number */
  uint64_t nprimes_in; /* nb of primes in the array on input */
  /* protected by roots_full/roots_empty */
  char_buffer_t roots_char; /* The local renumbering table is kept as a char*. */
  uint64_t nroots_tot; /* nb of entry in the local renumbering table */
  uint64_t nprimes_out; /* nb of treated primes */
  uint64_t nfreerels;
  unsigned long freerels_p[NB_PRIMES_PER_BUFFER];
  uint64_t freerels_first_index[NB_PRIMES_PER_BUFFER];
};

typedef struct th_buf_s th_buf_t[1];
typedef struct th_buf_s * th_buf_ptr;
typedef const struct th_buf_s * th_buf_srcptr;

static inline void
th_buf_init (th_buf_ptr buf, size_t size)
{
  char_buffer_init (buf->roots_char, size);
}

static inline void
th_buf_reset (th_buf_ptr buf)
{
  buf->nprimes_out = 0;
  buf->nfreerels = 0;
  buf->nroots_tot = 0;
  char_buffer_reset (buf->roots_char);
}

static inline void
th_buf_clear (th_buf_ptr buf)
{
  char_buffer_clear (buf->roots_char);
}

/******************** Data struct for sharing among threads *******************/
/* The main structure for the working pthreads pool */
typedef struct freerel_th_data_s {
  /* 4 pseudo semaphores. roots_* are also used for free_rel buffer. */
  sema_t primes_full, primes_empty, roots_full, roots_empty;
  /* Read only part */
  unsigned int nthreads;
  unsigned int th_id;
  unsigned long *lpb;
  unsigned long pmin, pmax, lpbmax;
  renumber_ptr tab;
  mpz_poly *pols;
  /* Read-write part */
  th_buf_t bufs[NB_BUFFERS_PER_THREAD];
} freerel_th_data_t;

/*************************** Main functions ***********************************/

/* This function produces all the primes from 2 to arg->lpbmax,
   for the workers pthreads pool (see next functions). The primes are in
   buffers; each buffer contains NB_PRIMES_PER_BUFFER primes (optimal: 1024),
   except the last one of course.
   After all buffers have been produced, this function creates again
   nb_pthreads buffers with only the end marker FREEREL_END_MARKER,
   in order to stop each worker.
*/
static void *
pthread_primes_producer (void *arg)
{
  /* This thread get the array of all the thread's data. */
  freerel_th_data_t *data = (freerel_th_data_t *) arg;
  unsigned long lpbmax = data->lpbmax, p = 2;
  unsigned int nthreads = data->nthreads, cur_buf = 0, cur_th = 0;

  prime_info pi;
  prime_info_init (pi);
  while (p <= lpbmax)
  {
    /* Here, we are sure that we can write at least one prime. */
    sema_wait (data[cur_th].primes_empty); /* Need an empty primes buffer */
    unsigned long *primes = data[cur_th].bufs[cur_buf]->primes;
    uint64_t nprimes = 0;

    do /* Main loop that produces primes */
    {
      primes[nprimes++] = p;
      p = getprime_mt (pi); /* get next prime */
    } while (nprimes < NB_PRIMES_PER_BUFFER && p <= lpbmax);

    data[cur_th].bufs[cur_buf]->nprimes_in = nprimes;
    sema_post (data[cur_th].primes_full); /* Drop the primes buffer */

    if (UNLIKELY (++cur_th == nthreads))
    {
      cur_th = 0;
      if (UNLIKELY (++cur_buf == NB_BUFFERS_PER_THREAD))
        cur_buf = 0;
    }
  }
  prime_info_clear (pi);

  /* We have to produce a special end buffer in each pth[].primes */
  unsigned int bak_th = cur_th;
  do
  {
    sema_wait (data[cur_th].primes_empty);
    data[cur_th].bufs[cur_buf]->primes[0] = FREEREL_END_MARKER;
    data[cur_th].bufs[cur_buf]->nprimes_in = 1;
    sema_post (data[cur_th].primes_full);

    if (UNLIKELY (++cur_th == nthreads))
    {
      cur_th = 0;
      if (UNLIKELY (++cur_buf == NB_BUFFERS_PER_THREAD))
        cur_buf = 0;
    }
  } while (cur_th != bak_th);

  return NULL;
}

/* This function takes a primes buffer and generates the corresponding
   local renumbering table and free relations.
   The local renumbering table is directly on the form of a ASCII buffer.
   The free relations must be recomputed by the main program, so a free relation
   contains only its prime and the first index (with respect to the local
   renumbering table).
   FIXME: for now it only works with 2 polynomials.
   This functions stops when the primes buffer contains only one
   FREEREL_END_MARKER, the end marker.
*/
static void *
pthread_primes_consumer (void *arg)
{
  freerel_th_data_t *my_data = arg;
  int nroots[NB_POLYS_MAX];
  unsigned long roots[NB_POLYS_MAX][MAX_DEGREE];

  /* Make local copy of some global (read only) variables and set some pointer
   * to global (read only) data.
   */
  int deg[NB_POLYS_MAX]; /* degree of the polynomials */
  mpz_srcptr lc[NB_POLYS_MAX]; /* pointer to the leading coefficients */
  unsigned int npolys = my_data->tab->nb_polys;
  for (unsigned int side = 0; side < npolys; side++)
  {
    deg[side] = my_data->pols[side]->deg;
    lc[side] = my_data->pols[side]->coeff[deg[side]];
  }

  for (unsigned int cur_buf = 0 ; ; )
  {
    sema_wait (my_data->roots_empty); /* Need an empty roots buffer */
    sema_wait (my_data->primes_full); /* Need a produced primes buffer */

    th_buf_ptr local_buf = my_data->bufs[cur_buf];
    unsigned long *my_p = local_buf->primes;

    /* Have we received the special end buffer in primes buffer ? */
    if (UNLIKELY(local_buf->nprimes_in == 1 && my_p[0] == FREEREL_END_MARKER))
    {
      /* Yes: we create a special end buffer and exit */
      local_buf->freerels_p[0] = FREEREL_END_MARKER;
      local_buf->nfreerels = 1;
      sema_post (my_data->primes_empty); /* Drop my now consumed primes buffer */
      sema_post (my_data->roots_full); /* Give my produced roots buffer */
      return NULL; /* The job is done: exit the endless loop and the pthread. */
    }

    th_buf_reset (local_buf);
    char_buffer_ptr out = my_data->bufs[cur_buf]->roots_char;

    while (local_buf->nprimes_out < local_buf->nprimes_in)
    {
      unsigned long p = *my_p++;  /* p is the current prime */
      /* Compute the roots on each side */
      for (unsigned int side = 0; side < npolys; side++)
      {
        if (UNLIKELY(p > my_data->lpb[side]))
          nroots[side] = 0;
        else if (deg[side] == 1)
          nroots[side] = 1;
        else
        {
          nroots[side] = mpz_poly_roots_ulong (roots[side],
                                               my_data->pols[side], p);
          /* Check for a projective root.
           * The projective root is assigned the value p and is inserted in
           * first place in the roots array (for the next sort)
           */
          if (nroots[side] != deg[side] && mpz_divisible_ui_p (lc[side], p))
          {
            roots[side][nroots[side]++] = roots[side][0];
            roots[side][0] = p;
          }
          if (p <= my_data->tab->bad_ideals.max_p) /* can it be a bad ideal ? */
          {
            for (int i = 0; i < nroots[side]; )
            {
              unsigned long r = roots[side][i];
              if (renumber_is_bad (NULL, NULL, my_data->tab, p, r, side))
              {
                /* bad ideals -> remove this root from the list */
                for (int j = i+1; j < nroots[side]; j++)
                  roots[side][j-1] = roots[side][j];
                nroots[side]--;
              }
              else
                i++;
            }
          }
        }
      }

      /* Write in the temporary local renumbering table */
      out->cur += renumber_write_buffer_p (out->cur, out->end - out->cur, my_data->tab, p, roots, nroots);
      char_buffer_resize (out, RENUMBER_MAX_SIZE_PER_PRIME);

      /* Check if p corresponds to a freerel */
      /* FIXME: What should be done when npolys != 2 ? */
      if (npolys == 2)
      {
        if (UNLIKELY(nroots[0] == deg[0] && nroots[1] == deg[1] &&
                     p >= my_data->pmin && p <= my_data->pmax))
        {
          local_buf->freerels_p[local_buf->nfreerels] = p;
          local_buf->freerels_first_index[local_buf->nfreerels] =
                                              local_buf->nroots_tot;
          local_buf->nfreerels++;
        }
      }

      for (unsigned int side = 0; side < npolys; side++)
        local_buf->nroots_tot += nroots[side];
      local_buf->nprimes_out++;
    } /* Next p in current primes buffer */

    sema_post(my_data->primes_empty); /* Drop my now consumed primes buffer */
    sema_post(my_data->roots_full);   /* Give my produced roots buffer */
    if (UNLIKELY (++(cur_buf) == NB_BUFFERS_PER_THREAD)) /* Go to next buffer */
      cur_buf = 0;
  } /* Endless loop */
}

/* Generate all the free relations and the renumbering table.
 * Return the number of free relations */
static uint64_t
generate_renumber_and_freerels (const char *outfilename,
                                renumber_t renumber_table, cado_poly pol,
                                unsigned long pmin, unsigned long pmax,
                                unsigned int nthreads)
{
  /* open outfile */
  FILE *outfile = NULL;
  outfile = fopen_maybe_compressed (outfilename, "w");
  ASSERT_ALWAYS (outfile != NULL);

  /* precompute sum(pols[k]->deg for k in [1..npolys]) */
  unsigned int sum_degs = 0; /* sum of the degrees of all polynomials */
  for(int k = 0; k < pol->nb_polys; k++)
      sum_degs += pol->pols[k]->deg;

  /* Compute the large primes bounds from their log in base 2. */
  unsigned long lpb[NB_POLYS_MAX] = { 0 };
  for(int k = 0; k < pol->nb_polys; k++)
    lpb[k] = 1UL << renumber_table->lpb[k];
  unsigned long lpbmax = 1UL << renumber_table->max_lpb;

  printf ("Generating renumber table for 2 <= p <= %lu\n", lpbmax);
  printf ("Considering freerels for %lu <= p <= %lu\n", pmin, pmax);
  fflush (stdout);

  /* Init statistics */
  uint64_t nprimes_total = 0; /* Total number of the primes */
  uint64_t nfreerels_total = 0; /* Total number of the free relations */
  stats_data_t stats; /* struct for printing progress */
  /* will print report at 2^10, 2^11, ... 2^23 computed primes and every
   * 2^23 primes after that */
  stats_init (stats, stdout, &nprimes_total, 23, "Treated", "primes", "", "p");

  /* Init threads and threads data */
  freerel_th_data_t *th_data = NULL; /* data of the threads */
  pthread_t primes_producer_th; /* pthread for primes producer */
  pthread_t *primes_consumer_th = NULL;

  th_data = (freerel_th_data_t *) malloc (nthreads * sizeof (freerel_th_data_t));
  ASSERT_ALWAYS (th_data != NULL);
  primes_consumer_th = (pthread_t *) malloc (nthreads * sizeof (pthread_t));
  ASSERT_ALWAYS (primes_consumer_th != NULL);

  for (unsigned int i = 0; i < nthreads; i++)
  {
    sema_init (th_data[i].roots_full, 0);
    sema_init (th_data[i].roots_empty, NB_BUFFERS_PER_THREAD);
    sema_init (th_data[i].primes_full, 0);
    sema_init (th_data[i].primes_empty, NB_BUFFERS_PER_THREAD);
    th_data[i].th_id = i;
    th_data[i].nthreads = nthreads;
    th_data[i].pmin = pmin;
    th_data[i].pmax = pmax;
    th_data[i].lpbmax = lpbmax;
    th_data[i].tab = renumber_table;
    th_data[i].lpb = lpb;
    th_data[i].pols = pol->pols;
    for (unsigned int j = 0; j < NB_BUFFERS_PER_THREAD; j++)
      th_buf_init (th_data[i].bufs[j], INIT_SIZE_BUF_CHAR);
  }

  if (pthread_create (&(primes_producer_th), NULL, &pthread_primes_producer,
                      (void *) th_data))
  {
    perror ("pthread_primes_producer creation failed\n");
    exit (1);
  }
  for (unsigned int i = 0; i < nthreads; i++)
  {
    if (pthread_create (&(primes_consumer_th[i]), NULL,
                        &pthread_primes_consumer,
                        (void *) &(th_data[i])))
    {
      perror ("pthread_primes_consumer pthread creation failed\n");
      exit (1);
    }
  }

  unsigned int cur_th = 0;
  unsigned int cur_buf = 0;
  /* Main loop: load the roots & free rels buffers, and print them */
  for (;;)
  {
    sema_wait (th_data[cur_th].roots_full); /* Need a full roots buffer */
    th_buf_ptr local_buf = th_data[cur_th].bufs[cur_buf];

    /* Is it the special end buffer (with only the end marker ) ? */
    if (local_buf->nfreerels == 1 &&
        local_buf->freerels_p[0] == FREEREL_END_MARKER)
    {
      /* All is done: when a pthread which produces roots & free rels gives this
       * special end buffer, all the next pthreads have finished and give also
       * another special end buffer. So it's not useful to read all; one is
       * sufficient. But it's dirty, so I prefer "read" (consume the semaphores)
       * them all. Morever, some ASSERTs after the loop verify the values of the
       * pseudo semaphores; if all is OK, these values must be equal to their
       * initializations.
       */
      sema_post (th_data[cur_th].roots_empty);
      for (unsigned int th_bak = cur_th;;)
      {
        if (UNLIKELY(++cur_th == nthreads))
          cur_th = 0;
        if (UNLIKELY (cur_th == th_bak))
          break;
        sema_wait (th_data[cur_th].roots_full);
        /* Nothing to do after this, we drop directly this special end buffer. */
        sema_post (th_data[cur_th].roots_empty);
      }
      break; /* End of the main loop */
    }

    /* We write the part of the renumbering table computed by the thread */
    char_buffer_write (renumber_table->file, local_buf->roots_char);

    /* We have to recompute the real index in the renumber table for the free
       relations */
    for (uint64_t k = 0; k < local_buf->nfreerels; k++)
    {
      uint64_t index = local_buf->freerels_first_index[k] + renumber_table->size;
      uint64_t last_index = index + sum_degs;
      fprintf (outfile, "%lx,0:%" PRIx64, local_buf->freerels_p[k], index);
      index++;
      for (; index < last_index; index++)
        fprintf (outfile, ",%" PRIx64, index);
      fputc ('\n', outfile);
    }
    renumber_table->size += local_buf->nroots_tot;
    nfreerels_total += local_buf->nfreerels;
    nprimes_total += local_buf->nprimes_out;

    /* Drop the consumed root & free rels buffer */
    sema_post (th_data[cur_th].roots_empty);

    if (stats_test_progress(stats))
      stats_print_progress (stats, nprimes_total, 0, 0, 0);

    if (UNLIKELY (++cur_th == nthreads))
    {
      cur_th = 0;
      if (UNLIKELY (++cur_buf == NB_BUFFERS_PER_THREAD))
        cur_buf = 0;
    }
  }

  /* All is done! */
  stats_print_progress (stats, nprimes_total, 0, 0, 1);
  if (pthread_join (primes_producer_th, NULL))
  {
    perror ("Error, pthread primes producer stops abnormally\n");
    exit (1);
  }
  for (unsigned int i = 0; i < nthreads; i++)
  {
    if (pthread_join (primes_consumer_th[i], NULL))
    {
      perror ("Error, one of the threads pthread_primes_consumer stops "
              "abnormally\n");
      exit (1);
    }
    ASSERT(th_data[i].roots_full->value == 0);
    ASSERT(th_data[i].roots_empty->value == NB_BUFFERS_PER_THREAD);
    ASSERT(th_data[i].primes_full->value == 0);
    ASSERT(th_data[i].primes_empty->value == NB_BUFFERS_PER_THREAD);
    sema_destroy (th_data[i].roots_full);
    sema_destroy (th_data[i].roots_empty);
    sema_destroy (th_data[i].primes_full);
    sema_destroy (th_data[i].primes_empty);
    for (unsigned int j = 0; j < NB_BUFFERS_PER_THREAD; j++)
      th_buf_clear (th_data[i].bufs[j]);
  }
  free (th_data);
  free (primes_consumer_th);
  fclose_maybe_compressed (outfile, outfilename);
  return nfreerels_total;
}

static void declare_usage(param_list pl)
{
  param_list_decl_usage(pl, "poly", "input polynomial file");
  param_list_decl_usage(pl, "renumber", "output file for renumbering table");
  param_list_decl_usage(pl, "out", "output file for free relations");
  param_list_decl_usage(pl, "lpb0", "large primes bound on side 0");
  param_list_decl_usage(pl, "lpb1", "large primes bound on side 1");
  param_list_decl_usage(pl, "lpbs", "large primes bounds (comma-separated list) "
                                    "(for MNFS)");
  param_list_decl_usage(pl, "pmin", "do not create freerel below this bound");
  param_list_decl_usage(pl, "pmax", "do not create freerel beyond this bound");
  param_list_decl_usage(pl, "badideals", "file describing bad ideals (for DL)");
  param_list_decl_usage(pl, "lcideals", "(switch) Add ideals for the leading "
                                        "coeffs of the polynomials (for DL)");
  param_list_decl_usage(pl, "t", "number of threads");
}

static void
usage (param_list pl, char *argv0)
{
    param_list_print_usage(pl, argv0, stderr);
    exit(EXIT_FAILURE);
}

int
main (int argc, char *argv[])
{
  char *argv0 = argv[0];
  cado_poly cpoly;
  unsigned long pmin = 2, pmax = 0;
  uint64_t nfree;
  renumber_t renumber_table;
  int lcideals = 0;
  int lpb_arg[NB_POLYS_MAX] = { 0 };
  unsigned long lpb[NB_POLYS_MAX] = { 0 };
  unsigned int nthreads = 1;

  param_list pl;
  param_list_init(pl);
  declare_usage(pl);
  param_list_configure_switch(pl, "-lcideals", &lcideals);

#ifdef HAVE_MINGW
  _fmode = _O_BINARY; /* Binary open for all files */
#endif

  argv++, argc--;
  if (argc == 0)
    usage (pl, argv0);

  for( ; argc ; )
  {
    if (param_list_update_cmdline(pl, &argc, &argv))
      continue;
    FILE *f;
    if ((f = fopen(argv[0], "r")) != NULL)
    {
      param_list_read_stream(pl, f, 0);
      fclose(f);
      argv++,argc--;
      continue;
    }
    fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
    usage (pl, argv0);
  }
  /* print command-line arguments */
  verbose_interpret_parameters (pl);
  param_list_print_command_line (stdout, pl);
  fflush (stdout);

  const char * polyfilename = param_list_lookup_string(pl, "poly");
  const char * outfilename = param_list_lookup_string(pl, "out");
  const char * badidealsfilename = param_list_lookup_string(pl, "badideals");
  const char * renumberfilename = param_list_lookup_string(pl, "renumber");

  param_list_parse_ulong (pl, "lpb0", &(lpb[0]));
  param_list_parse_ulong (pl, "lpb1", &(lpb[1]));
  int narg = param_list_parse_int_list (pl, "lpbs", lpb_arg, NB_POLYS_MAX, ",");
  param_list_parse_ulong (pl, "pmin", &pmin);
  param_list_parse_ulong (pl, "pmax", &pmax);
  param_list_parse_uint (pl, "t", &nthreads);

  if (param_list_warn_unused(pl))
    usage (pl, argv0);

  if (polyfilename == NULL)
  {
    fprintf (stderr, "Error, missing -poly command line argument\n");
    usage (pl, argv0);
  }
  if (renumberfilename == NULL)
  {
    fprintf (stderr, "Error, missing -renumber command line argument\n");
    usage (pl, argv0);
  }
  if (outfilename == NULL)
  {
    fprintf (stderr, "Error, missing -out command line argument\n");
    usage (pl, argv0);
  }
  if (nthreads == 0 || nthreads > 512)
  {
    fprintf (stderr, "Error, the number of threads is incorrect, it must be "
                      "between 1 and 512.\n");
  }

  cado_poly_init(cpoly);
  if (!cado_poly_read (cpoly, polyfilename))
  {
    fprintf (stderr, "Error reading polynomial file\n");
    exit (EXIT_FAILURE);
  }

  if (narg == 0) /* lpbs were given as -lpb0 and -lpb1 */
  {
    if (cpoly->nb_polys > 2) /* With more than 2 polys, must use -lpbs. */
    {
      fprintf (stderr, "Error, missing -lpbs command line argument\n");
      usage (pl, argv0);
    }
    if (lpb[0] == 0 || lpb[1] == 0)
    {
      fprintf (stderr, "Error, missing -lpb0 and/or -lpb1 command line "
                       "argument\n");
      usage (pl, argv0);
    }
  }
  else /* lpbs were given as -lpbs x,x,x,... */
  {
    if (narg != cpoly->nb_polys)
    {
      fprintf (stderr, "Error, the number of values given in -lpbs does not "
                       "correspond to the number of polynomials\n");
      usage (pl, argv0);
    }
    if (lpb[0] != 0)
      fprintf (stderr, "Warning, the value given by -lpb0 will be ignored, "
                       "the one given by -lpbs will be used\n");
    if (lpb[1] != 0)
      fprintf (stderr, "Warning, the value given by -lpb1 will be ignored, "
                       "the one given by -lpbs will be used\n");
    for (int i = 0; i < cpoly->nb_polys; i++)
    {
      if (lpb_arg[i] <= 0)
      {
        fprintf (stderr, "Error, -lpbs command line argument cannot contain "
                         "non-positive values\n");
        usage (pl, argv0);
      }
      lpb[i] = lpb_arg[i];
    }
  }

  int ratside = cado_poly_get_ratside (cpoly);
  uint64_t nonmonic = 0;
  for (unsigned int i = cpoly->nb_polys; i > 0 ; i--) 
    nonmonic = (nonmonic << 1) + mpz_poly_is_nonmonic (cpoly->pols[i-1]);

  renumber_init_for_writing (renumber_table, cpoly->nb_polys, ratside, lcideals,
                             nonmonic, lpb);

  /* if pmax is not equal to 0 (i.e., was not given on the command line), 
   * set pmax to the *minimum* of the large prime bounds, since larger primes
   * will never occur on both sides.
   * We generate the renumbering table from the first prime (2) and
   * up to the *maximum* of the large prime bounds.
   */
  unsigned long lpbmin = lpb[0];
  for(int k = 1; k < cpoly->nb_polys; k++)
    lpbmin = MIN(lpbmin, lpb[k]);
  if (!pmax)
    pmax = 1UL << lpbmin;

  renumber_write_open (renumber_table, renumberfilename, badidealsfilename,
                       cpoly);

  nfree = generate_renumber_and_freerels (outfilename, renumber_table, cpoly,
                                          pmin, pmax, nthreads);

  /* /!\ Needed by the Python script. /!\ */
  fprintf (stderr, "# Free relations: %" PRIu64 "\n", nfree);
  fprintf (stderr, "Renumbering struct: nprimes=%" PRIu64 "\n",
                   renumber_table->size);

  /* produce an error when index_t is too small to represent all ideals */
  if (renumber_table->size >> (8 * __SIZEOF_INDEX__))
    {
      fprintf (stderr, "Error, please increase __SIZEOF_INDEX__\n");
      fprintf (stderr, "(see local.sh.example)\n");
      exit (1);
    }

  renumber_write_close (renumber_table, renumberfilename);
  renumber_clear (renumber_table);
  cado_poly_clear (cpoly);
  param_list_clear(pl);

  return 0;
}
