#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <bits/types/struct_rusage.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <sys/time.h> // IWYU pragma: keep
#include <sys/wait.h> // IWYU pragma: keep
#ifdef HAVE_GETRUSAGE
#include <sys/resource.h> // IWYU pragma: keep
#endif
#include <pthread.h>

#include "macros.h"     // MAYBE_UNUSED
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

int cado_fd_popen(const char * command, const char * mode)
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
        abort();
    }
    int pipefd[2];
    if (pipe(pipefd) < 0) {
        perror("pipe");
        return -1;
    }
    /* pipefd[0] is the read end, pipefd[1] is the write end */
    pthread_mutex_lock(popenlist->m);

    popenlist->p[popenlist->n].fd = pipefd[imode];
    pid_t * kid = &(popenlist->p[popenlist->n].kid);
    popenlist->n++;

    pid_t child = fork();
    if (child < 0) {
        pthread_mutex_unlock(popenlist->m);
        perror("fork");
        return -1;
    }
    if (child) {
        /* I'm the father. I only want to use pipefd[imode]. */
        close(pipefd[!imode]);
        *kid = child;
        pthread_mutex_unlock(popenlist->m);
        /*
        fprintf(stderr, "%s child %d (%s) through fd %d\n",
                imode ?  "Writing to" : "Reading from",
                child, command, pipefd[imode]);
                */
        return pipefd[imode];
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
        perror("execl() failed");
        fprintf (stderr, "execl command size is %zu\n",
                strlen (command) + strlen ("sh -c "));
        exit(EXIT_FAILURE);
    }
    return -1;
}

FILE * cado_popen(const char * command, const char * mode)
{
    int fd = cado_fd_popen(command, mode);
    if (fd < 0) return NULL;
    return fdopen(fd, mode);
}

/* remove the fd from our list of popen'ed files, and return the kid id.
 */
static int reap_fd(int fd)
{
    pid_t kid = 0;
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
    return kid;
}

#if defined(HAVE_GETRUSAGE) && defined(HAVE_WAIT4)
int cado_fd_pclose2(int fd, struct rusage * nr)
#else
int cado_fd_pclose2(int fd, void * nr MAYBE_UNUSED)
#endif
{
    int kid = reap_fd(fd);
    /* close the kid's fd, which will trigger its termination -- at least
     * if it's reading from it. If it's writing to it, at this point we
     * expect the child's source to have been consumed, and thus the
     * write end of the pipe to have been closed already */
    close(fd);
    /* we must wait() for the kid now */
    int status, error;
#if defined(HAVE_GETRUSAGE) && defined(HAVE_WAIT4)
    struct rusage r[1];
    error = wait4(kid, &status, 0, r);
#else
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
#if defined(HAVE_GETRUSAGE)
#if defined(HAVE_WAIT4)
    if (nr) memcpy(nr, r, sizeof(struct rusage));
#else
    if (nr) memset(nr, 0, sizeof(struct rusage));
#endif
#endif
    /* Linux man page: 
       "returns the exit status of the command as returned by wait4(2)"
    */
    return status;
}

#if defined(HAVE_GETRUSAGE) && defined(HAVE_WAIT4)
int cado_pclose2(FILE * stream, struct rusage * nr)
#else
int cado_pclose2(FILE * stream, void * nr MAYBE_UNUSED)
#endif
{
    int fd = fileno(stream);
    fclose(stream);
    return cado_fd_pclose2(fd, nr);
}

