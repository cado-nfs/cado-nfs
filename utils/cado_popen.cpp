#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include <mutex>
#include <map>
#include <string>
#include <stdexcept>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#ifdef HAVE_GETRUSAGE
#include <sys/resource.h> // IWYU pragma: keep
#endif
#include <fcntl.h>

#include "fmt/base.h"

#include "macros.h"     // MAYBE_UNUSED
#include "cado_popen.hpp"

/* We go to a certain level of complication here because pclose()ing a
 * file stream implicitly entails wait()ing for the child process
 * *BUT* we would like to get back the child status as well as the rusage
 * struct that are returned by that wait() call (actually wait4() or
 * possibly waitpid() for just the status)
 */

/* The mutex cado_popen_lock has three distinct uses:
 *
 * 1 - protect accesses to the fd_to_pid table. We want this table so
 * that we can wait() the children as soon as we close the file
 * descriptors that are tied to them.
 *
 * 2 - ensure serialization of the pipe()+fcntl() and fork() calls, in
 * order to avoid the situation where a forked child lacks the
 * close-on-exec flags. (pipe2 would alleviate this)
 *
 * 3 - ensure serialization of the close()+access of the fd_to_pid table,
 * otherwise the table can be modified inbetween.
 *
 */

/* we need to prevent static destruction order fiasco here :-( */
static std::mutex& get_popen_lock() {
    // Heap allocation prevents the destruction order fiasco.
    // The pointer itself is static, initialized thread-safely on first use.
    static auto * lock = new std::mutex();
    return *lock;
}

static std::map<int, pid_t>& get_fd_to_pid() {
    // Safe from both initialization and destruction fiascos.
    static auto * table = new std::map<int, pid_t>();
    return *table;
}

int cado_fd_popen(const char * command, const char * mode)
{
    /* Is this a pipe for reading or for writing ? */
    int imode = -1; /* 0: reading; 1: writing */
    if (std::string(mode) == "r") {
        imode = 0;
    } else if (std::string(mode) == "rb") {
        imode = 0;
    } else if (std::string(mode) == "w") {
        imode = 1;
    } else if (std::string(mode) == "wb") {
        imode = 1;
    } else {
        fmt::print(stderr, "Please fix {}\n", __func__);
        abort();
    }
    int pipefd[2];

    pid_t child;

    try {
        const std::scoped_lock dummy(get_popen_lock());
        if (pipe(pipefd) < 0) {
            perror("pipe");
            return -1;
        }

        fcntl(pipefd[0], F_SETFD, FD_CLOEXEC);
        fcntl(pipefd[1], F_SETFD, FD_CLOEXEC);

        /* reserve a slot in the map. Throws if OOM. */
        auto [iter, fresh] = get_fd_to_pid().emplace(pipefd[imode], -1);
        if (!fresh)
            throw std::logic_error("key clash in map");

        /* We _must_ serialize pipe() and fork(), otherwise the forked
         * child might not have the FD_CLOEXEC flag set!! */
        child = fork();
        if (child < 0) {
            get_fd_to_pid().erase(iter);
            perror("fork");
            close(pipefd[0]);
            close(pipefd[1]);
            return -1;
        } else {
            iter->second = child;
        }
    } catch(...) {
        /* either OOM or logic_error. Both deserve to be re-thrown */
        close(pipefd[0]);
        close(pipefd[1]);
        throw;
    }

    /* pipefd[0] is the read end, pipefd[1] is the write end */

    if (child) {
        /* I'm the father. I only want to use pipefd[imode]. */
        close(pipefd[imode^1]);
        /*
        fmt::print(stderr, "{} child {} ({}) through fd {}\n",
                imode ?  "Writing to" : "Reading from",
                child, command, pipefd[imode]);
                */
        return pipefd[imode];
    } else {
        /* if father wants to read, we close our standard input
         * (0==imode), and bind our standard output (1==imode^1) to the
         * other end. */
        close(imode);
        close(pipefd[imode]);
        dup2(pipefd[imode^1], imode^1);
        close(pipefd[imode^1]);
        /*
        fmt::print(stderr, "Child process {} (parent {}) executing {} ; closed [{} {} {}] open [{}]\n",
                getpid(),
                getppid(),
                command,
                imode, pipefd[imode], pipefd[imode^1],
                imode^1
                );
                */
        execl("/bin/sh", "sh", "-c", command, NULL);
        perror("execl() failed");
        exit(EXIT_FAILURE);
    }
    return -1;
}

FILE * cado_popen(const char * command, const char * mode)
{
    const int fd = cado_fd_popen(command, mode);
    if (fd < 0) return nullptr;
    return fdopen(fd, mode);
}

/* assume the fd has has just been closed. remove it from our list of
 * popen'ed files, and reap the kid.
 */
#if defined(HAVE_GETRUSAGE) && defined(HAVE_WAIT4)
static int reap_pid(pid_t kid, struct rusage * nr)
#else
static int reap_pid(pid_t kid, void * nr MAYBE_UNUSED)
#endif
{
    /* we must wait() for the kid now */
    int status, error;
#if defined(HAVE_GETRUSAGE) && defined(HAVE_WAIT4)
    if (nr)
        error = wait4(kid, &status, 0, nr);
    else
        error = waitpid(kid, &status, 0);
#else
    error = waitpid(kid, &status, 0);
#endif
    if (error == -1) {
        return -1;
    }

    /*
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

    double u = r->ru_utime.tv_sec + (r->ru_utime.tv_usec / 1.0e6);
    double s = r->ru_stime.tv_sec + (r->ru_stime.tv_usec / 1.0e6);
    fprintf(stderr, "Child process %d %s, having spent %.2fs+%.2fs on cpu\n",
            kid, long_status, u, s);
            */
    return status;
}

#if defined(HAVE_GETRUSAGE) && defined(HAVE_WAIT4)
int cado_fd_pclose2(int fd, struct rusage * nr)
#else
int cado_fd_pclose2(int fd, void * nr MAYBE_UNUSED)
#endif
{
    pid_t kid;
    {
        const std::scoped_lock dummy(get_popen_lock());
        if (close(fd) < 0)
            perror("fork");
        auto iter = get_fd_to_pid().find(fd);
        ASSERT_ALWAYS(iter != get_fd_to_pid().end());
        kid = iter->second;
        get_fd_to_pid().erase(iter);
    }

    return reap_pid(kid, nr);
}

#if defined(HAVE_GETRUSAGE) && defined(HAVE_WAIT4)
int cado_pclose2(FILE * stream, struct rusage * nr)
#else
int cado_pclose2(FILE * stream, void * nr MAYBE_UNUSED)
#endif
{
    pid_t kid;
    {
        const std::scoped_lock dummy(get_popen_lock());
        const int fd = fileno(stream);
        /* This also closes fd */
        fclose(stream);
        auto iter = get_fd_to_pid().find(fd);
        ASSERT_ALWAYS(iter != get_fd_to_pid().end());
        kid = iter->second;
        get_fd_to_pid().erase(iter);
    }

    return reap_pid(kid, nr);
}

