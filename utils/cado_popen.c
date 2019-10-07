#include "cado.h"
#include <stdio.h>
#include <errno.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/time.h>
#include <sys/wait.h>
#include <pthread.h>

#include "macros.h"
#include "cado_popen.h"

/* We need to close file descriptors underlying other popen()-ed calls on
 * the children processes.  Not the underlying streams though, since
 * closing is done in userland, and only once.
 */
static struct {
    pthread_mutex_t m[1];
    int n;
    struct { int fd; pid_t kid; } p[1024];
} popenlist[1] = {{{PTHREAD_MUTEX_INITIALIZER}, 0, {{0,0},}}};

FILE * cado_popen(const char * command, const char * mode)
{
    /* Is this a pipe for reading or for writing ? */
    int imode = -1; /* 0: reading; 1: writing */
    if (strcmp(mode, "r") == 0) {
        imode = 0;
    } else if (strcmp(mode, "rb") == 0) {
        imode = 0;
    } else if (strcmp(mode, "w") == 0) {
        imode = 1;
    } else if (strcmp(mode, "wb") == 0) {
        imode = 1;
    } else {
        fprintf(stderr, "Please fix %s\n", __func__);
    }
    int pipefd[2];
    if (pipe(pipefd) < 0) {
        perror("pipe");
        return NULL;
    }
    /* pipefd[0] is the read end, pipefd[1] is the write end */
    pthread_mutex_lock(popenlist->m);

    popenlist->p[popenlist->n].fd = pipefd[imode];
    pid_t * kid = &(popenlist->p[popenlist->n].kid);
    popenlist->n++;

    pid_t child = fork();
    if (child < 0) { perror("fork"); return NULL; }
    if (child) {
        pthread_mutex_unlock(popenlist->m);
        /* I'm the father. I only want to use pipefd[imode]. */
        close(pipefd[!imode]);
        *kid = child;
        /*
        fprintf(stderr, "%s child %d (%s) through fd %d\n",
                imode ?  "Writing to" : "Reading from",
                child, command, pipefd[imode]);
                */
        return fdopen(pipefd[imode], mode);
    } else {
        /* if father wants to read, we close our standard input
         * (0==imode), and bind our standard output (1==!imode) to the
         * other end. */
        /* We still have the lock here. We don't care much about
         * releasing it, since we're going to do exec() anyway */
        close(imode);
        close(pipefd[imode]);
        dup2(pipefd[!imode], !imode);
        for(int i = 0 ; i < popenlist->n - 1 ; i++) {
            int fd = popenlist->p[i].fd;
            int kid = popenlist->p[i].kid;
            int rc = close(fd);
            if (rc < 0) {
                fprintf(stderr, "Process %d closing fd %d (%s pid %d)"
                        ": %s\n",
                        getpid(), fd,
                        imode ? "->" : "<-", kid, strerror(errno));
            }
        }
        popenlist->n = 0;       /* who cares, we're exec()'ing anyway */
        /*
        fprintf(stderr, "Child process %d (parent %d) executing %s\n",
                getpid(),
                getppid(),
                command);
                */
        execl("/bin/sh", "sh", "-c", command, NULL);
        char const *error_name;
        switch (errno)
          {
          case E2BIG: error_name = "E2BIG"; break;
          case EACCES: error_name = "EACCES"; break;
          case EFAULT: error_name = "EFAULT"; break;
          case EIO: error_name = "EIO"; break;
          case ELOOP: error_name = "ELOOP"; break;
          case ENAMETOOLONG: error_name = "ENAMETOOLONG"; break;
          case ENOENT: error_name = "ENOENT"; break;
          case ENOEXEC: error_name = "ENOEXEC"; break;
          case ENOMEM: error_name = "ENOMEM"; break;
          case ENOTDIR: error_name = "ENOTDIR"; break;
          case ETXTBSY: error_name = "ETXTBSY"; break;
          default: error_name ="unknown"; break;
          }
        fprintf (stderr, "%s\n", error_name);
        fprintf (stderr, "command size is %lu\n", strlen (command)
                 + strlen ("/bin/sh -c "));
        perror("execve() failed");
        exit(1);
    }
    return NULL;
}

#ifdef HAVE_GETRUSAGE
int cado_pclose2(FILE * stream, struct rusage * nr)
#else
int cado_pclose2(FILE * stream, void * nr MAYBE_UNUSED)
#endif
{
    pid_t kid = 0;
    int fd = fileno(stream);
    pthread_mutex_lock(popenlist->m);
    int nn = 0;
    for(int i = 0 ; i < popenlist->n ; i++) {
        popenlist->p[nn] = popenlist->p[i];
        if (popenlist->p[i].fd == fd) {
            if (kid) {
                fprintf(stderr, "Error: two or more child processes (%d and %d) hold a reference to fd %d\n", kid, popenlist->p[i].kid, fd);
                abort();
            }
            ASSERT_ALWAYS(kid == 0);
            kid = popenlist->p[i].kid;
        } else {
            nn++;
        }
    }
    popenlist->n = nn;
    pthread_mutex_unlock(popenlist->m);
    /* close the kid's fd, which will trigger its termination -- at least
     * if it's reading from it. If it's writing to it, at this point we
     * expect the child's source to have been consumed, and thus the
     * write end of the pipe to have been closed already */
    fclose(stream);
    close(fd);
    /* we must wait() for the kid now */
    int status, error;
#ifdef HAVE_GETRUSAGE
    struct rusage r[1];
    error = wait4(kid, &status, 0, r);
#else
    /* if we have no getrusage(), we probably don't have wait4 either */
    error = waitpid(kid, &status, 0);
#endif
    if (error == -1) {
        return -1;
    }
    char long_status[80] = {'\0'};

    if (WIFEXITED(status)) {
        snprintf(long_status, sizeof(long_status),
                "exited with code %d", WEXITSTATUS(status));
    } else if (WIFSIGNALED(status)) {
        snprintf(long_status, sizeof(long_status),
                "terminated by signal %d%s",
                WEXITSTATUS(status), WCOREDUMP(status) ? ", with core" : "");
    } else {
        snprintf(long_status, sizeof(long_status), "[weird status %d]", status);
    }

    /*
    double u = r->ru_utime.tv_sec + (r->ru_utime.tv_usec / 1.0e6);
    double s = r->ru_stime.tv_sec + (r->ru_stime.tv_usec / 1.0e6);
    fprintf(stderr, "Child process %d %s, having spent %.2fs+%.2fs on cpu\n",
            kid, long_status, u, s);
            */
#ifdef HAVE_GETRUSAGE
    if (nr) memcpy(nr, r, sizeof(struct rusage));
#endif
    /* Linux man page: 
       "returns the exit status of the command as returned by wait4(2)"
    */
    return status;
}
