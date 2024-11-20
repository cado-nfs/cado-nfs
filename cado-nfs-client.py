#!/usr/bin/env python3

# pylint: disable=too-many-lines
# pylint: disable=deprecated-module
# pylint: disable=invalid-name
# pylint: disable=useless-object-inheritance
# pylint: disable=too-many-arguments
# pylint: disable=too-few-public-methods
# pylint: disable=import-error
# pylint: disable=wrong-import-position

# {{{ libs
import sys
import io
import os
import copy
import random
import errno
import stat
import optparse
import time
import subprocess
import hashlib
import logging
import socket
import signal
import re
import json
from io import BytesIO
from string import Template
import ssl
import requests
from urllib.parse import urlparse

if sys.hexversion < 0x03060000:
    sys.exit("Python 3.6 or newer is required to run this program.")

# THIS PART MUST BE EXACTLY IDENTICAL IN cado-nfs.py and cado-nfs-client.py

# Three possible locations for this script
#
# The installed path. Then we should rely on the @-escaped paths, and
# determine the location of all our utility stuff based on that.
#
# The build directory. We rely on the presence of a magic file called
# source-location.txt, and we fetch the python source code from there
#
# The source tree. We call the ./scripts/build_environment.sh script to
# determine where the binaries are being put.

import locale

pathdict = dict()

one_pyfile_example_subpath = "scripts/cadofactor/workunit.py"


def detect_installed_tree(pathdict):
    mydir = os.path.normpath(os.path.dirname(sys.argv[0]))
    md = mydir.split(os.path.sep)
    install_tree = mydir
    bd = "@BINSUFFIX@".split(os.path.sep)
    while md and bd and md[-1] == bd[-1]:
        md.pop()
        bd.pop()
        install_tree = os.path.normpath(os.path.join(install_tree, ".."))
    if bd:
        if os.environ.get("CADO_NFS_DEBUG_PATHDETECT"):
            print("{} does not end in @BINSUFFIX@".format(mydir))
        return False
    example = os.path.join(install_tree,
                           "@LIBSUFFIX@",
                           one_pyfile_example_subpath)
    t = os.path.exists(example)
    if not t:
        if os.environ.get("CADO_NFS_DEBUG_PATHDETECT"):
            print("{} does not exist".format(example))
        return False

    # make all this relocatable, it doesn't cost us much.
    # (note though that the rpaths in the binaries are likely to still
    # contain absolute paths)
    pathdict["pylib"] = os.path.join(install_tree, "@LIBSUFFIX@/scripts")
    pathdict["data"] = os.path.join(install_tree, "@DATASUFFIX@")
    pathdict["lib"] = os.path.join(install_tree, "@LIBSUFFIX@")
    pathdict["bin"] = os.path.join(install_tree, "@BINSUFFIX@")

    if os.environ.get("CADO_NFS_DEBUG_PATHDETECT"):
        print("cado-nfs running in installed tree")
    return True


def detect_build_tree(pathdict):
    # source-location.txt is created by our build system, and can be used
    # *ONLY* when we call this script from the build directory. We don't
    # want this to perspire in any installed file, of course.
    mydir = os.path.normpath(os.path.dirname(sys.argv[0]))
    source_location_subpath = "source-location.txt"
    source_location_file = os.path.join(mydir, source_location_subpath)
    if not os.path.exists(source_location_file):
        if os.environ.get("CADO_NFS_DEBUG_PATHDETECT"):
            print("{} does not exist".format(source_location_file))
        return False

    # ok, we're in the build tree, apparently
    source_tree = open(source_location_file, "r").read().strip()

    pathdict["pylib"] = os.path.join(source_tree, "scripts")
    pathdict["data"] = os.path.join(source_tree, "parameters")
    pathdict["lib"] = mydir
    pathdict["bin"] = mydir

    if os.environ.get("CADO_NFS_DEBUG_PATHDETECT"):
        print("cado-nfs running in build tree")
    return True


def detect_source_tree(pathdict):
    mydir = os.path.normpath(os.path.dirname(sys.argv[0]))
    helper = os.path.join(mydir, "scripts/build_environment.sh")
    if not os.path.exists(helper):
        if os.environ.get("CADO_NFS_DEBUG_PATHDETECT"):
            print("{} does not exist".format(helper))
        return False
    pipe = subprocess.Popen([helper, "--show"], stdout=subprocess.PIPE)
    loc = locale.getlocale()[1]
    if not loc:
        loc = "ascii"
    output = pipe.communicate()[0].decode(loc)
    cado_bin_path = [x.split("=", 2)[1]
                     for x in output.split("\n")
                     if re.match("^build_tree", x)][0]
    cado_bin_path = re.sub("^\"(.*)\"$", "\\1", cado_bin_path)

    pathdict["pylib"] = os.path.join(mydir, "scripts")
    pathdict["data"] = os.path.join(mydir, "parameters")
    pathdict["lib"] = cado_bin_path
    pathdict["bin"] = cado_bin_path

    if os.environ.get("CADO_NFS_DEBUG_PATHDETECT"):
        print("cado-nfs running in source tree")
    return True


if detect_installed_tree(pathdict):
    pass
elif detect_build_tree(pathdict):
    pass
elif detect_source_tree(pathdict):
    pass
else:
    raise RuntimeError("We're unable to determine"
                       " the location of the cado-nfs binaries"
                       " and python files")

sys.path.append(pathdict["pylib"])

# END OF THE PART THAT MUST BE IDENTICAL IN cado-nfs.py and cado-nfs-client.py

from cadofactor.workunit import Workunit    # noqa: E402
from cadofactor import cadologger           # noqa: E402

# }}}


def pid_exists(pid):
    try:
        os.kill(pid, 0)
    except OSError as e:
        return e.errno == errno.EPERM
    else:
        return True


# {{{ locking plumbing.
# File locking functions are specific to Unix/Windows/MacOS platforms.
# The FileLock class is an Interface with static methods.

# we used to have mostly-placeholder mostly-unimplemented code for
# Windows, but Windows is obsolete now, so we don't need to bother.

# could replace "posix" by "xxx" here if os.name is "posix" but you still get
# the error message "IOError: [Errno 37] No locks available"
# See this thread on the mailing list.
# https://sympa.inria.fr/sympa/arc/cado-nfs/2016-05/msg00010.html

# Note however that in the absence of proper file locking, any
# computation that puts some load on the work unit server is bound to
# fail.

if os.name == "posix":
    import fcntl

    class FileLock(object):
        @staticmethod
        def lock(filehandle, exclusive=False, blocking=True):
            """ Lock a file

            If exclusive is True, lock for exclusive (a.k.a "write") access,
            otherwise lock for shared (a.k.a. "read") access.
            If blocking is False, don't block in case of already-locked
            file, but raise IOError with EACCES or EAGAIN (depending on OS).
            """
            mode = fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH
            mode |= 0 if blocking else fcntl.LOCK_NB
            fcntl.flock(filehandle.fileno(), mode)

        @staticmethod
        def unlock(filehandle):
            """ Unlock a file """
            fcntl.flock(filehandle.fileno(), fcntl.LOCK_UN)
else:

    # No file locking. FIXME: What about MacOS?
    class FileLock(object):
        @staticmethod
        def lock(filehandle, exclusive=False, blocking=True):
            """ Do nothing """
            pass

        @staticmethod
        def unlock(filehandle):
            """ Do nothing """
            pass
# }}}


def create_daemon(workdir=None, umask=None, logfile=None):  # {{{
    """Run a sub-process, detach it from the control tty.

    This is a simplified version of the code found there.

    http://code.activestate.com/recipes/278731-creating-a-daemon-the-python-way/

    Changes: workdir is now a parameter, daemon changes CWD only if workdir
    parameter is specified. umask is also a parameter, and the process' umask
    is set only if a value is specified.
    """

    # Default maximum for the number of available file descriptors.
    maxfd_default = 1024

    # The standard I/O file descriptors are redirected to /dev/null by default.
    if hasattr(os, "devnull"):
        redirect_to = os.devnull
    else:
        redirect_to = "/dev/null"

    try:
        # Fork a child process so the parent can exit.  This returns
        # control to the command-line or shell.  It also guarantees that
        # the child will not be a process group leader, since the child
        # receives a new process ID and inherits the parent's process
        # group ID.  This step is required to insure that the next call
        # to os.setsid is successful.
        pid = os.fork()
    except OSError as e:
        raise Exception("%s [%d]" % (e.strerror, e.errno))

    if pid > 0:  # master
        sys.stdout.write("PID: %d\n" % pid)
        sys.stdout.flush()
        sys.exit()

    # To become the session leader of this new session and the process
    # group leader of the new process group, we call os.setsid().  The
    # process is also guaranteed not to have a controlling terminal.
    os.setsid()

    # Since the current working directory may be a mounted filesystem, we
    # avoid the issue of not being able to unmount the filesystem at
    # shutdown time by changing it to the root directory.
    if workdir is not None:
        os.chdir(workdir)

    # We probably don't want the file mode creation mask inherited from
    # the parent, so we give the child complete control over permissions.
    if umask is not None:
        os.umask(umask)

    if logfile is not None:
        # must remove the intermediary handlers that the logging system
        # uses, otherwise we get inconsistent file position and python
        # gets nuts.
        logger = logging.getLogger()
        for handler in list(logger.handlers):  # Remove old handlers
            logger.removeHandler(handler)

    # Iterate through and close all file descriptors.
    fdlist = None

    if fdlist is None:
        try:
            fdlist = [int(c) for c in os.listdir('/proc/self/fd')]
        except FileNotFoundError:
            pass

    if fdlist is None:
        try:
            fdlist = [int(c) for c in os.listdir('/dev/fd')]
        except FileNotFoundError:
            pass

    if fdlist is None:
        import resource     # Resource usage information.
        maxfd = resource.getrlimit(resource.RLIMIT_NOFILE)[1]
        if maxfd == resource.RLIM_INFINITY:
            maxfd = maxfd_default

        fdlist = range(0, maxfd)

    for fd in fdlist:
        try:
            if logfile is None or fd != logfile.fileno():
                os.close(fd)
        except OSError:  # ERROR, fd wasn't open to begin with (ignored)
            pass

    # Redirect the standard I/O file descriptors to the specified file.
    # Since the daemon has no controlling terminal, most daemons redirect
    # stdin, stdout, and stderr to /dev/null.  This is done to prevent
    # side-effects from reads and writes to the standard I/O file
    # descriptors.

    # This call to open is guaranteed to return the lowest file
    # descriptor, which will be 0 (stdin), since it was closed above.
    fd0 = os.open(redirect_to, os.O_RDWR)   # standard input (0)

    fd12 = fd0
    if logfile is not None:
        fd12 = logfile.fileno()

    # Duplicate standard input to standard output and standard error.
    os.dup2(fd12, 1)            # standard output (1)
    os.dup2(fd12, 2)            # standard error (2)

    if logfile is None:
        return

    # Now re-plug the logging system to the same file descriptor as
    # stderr. we have three file descriptors open to the same file, by
    # the way. We might as well decide to do away with one of them (e.g.,
    # logfile.fileno())
    logger.addHandler(logging.StreamHandler(sys.stderr))
    # os.close(logfile.fileno())
# }}}


# {{{ exclusive open/close
class FileLockedException(IOError):
    """ Locking a file for exclusive access failed """
    pass


def open_exclusive(filename):
    """ Open a file and get an exclusive lock on it """
    fileobj = open(filename, "r+")
    try:
        FileLock.lock(fileobj, exclusive=True, blocking=False)
    except IOError as err:
        if err.errno == errno.EACCES or err.errno == errno.EAGAIN:
            fileobj.close()
            raise FileLockedException(errno.EACCES, "File locked", filename)
        raise
    return fileobj


def close_exclusive(fileobj):
    """ Close a file, releasing any held lock on it """
    FileLock.unlock(fileobj)
    fileobj.close()
# }}}


# {{{ run shell command, capture std streams
def run_command(command, stdin=None, print_error=True, **kwargs):
    """ Run command, wait for it to finish, return exit status, stdout
    and stderr

    If print_error is True and the command exits with a non-zero exit code,
    print stdout and stderr to the log. If a KeyboardInterrupt exception
    occurs while waiting for the command to finish, the command is
    terminated.
    """

    command_str = command if isinstance(command, str) else " ".join(command)
    command_list = command if isinstance(command, list) else command.split(" ")

    # changed close_fds from True to False, since otherwise the 'las'
    # clients are not killed when merge starts
    # see https://gforge.inria.fr/tracker/?func=detail&aid=21718
    close_fds = False

    logging.info("Running %s", command_str)

    child = subprocess.Popen(command_list,
                             stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             close_fds=close_fds,
                             **kwargs)

    logging.info("[%s] Subprocess has PID %d", time.asctime(), child.pid)

    # If we receive SIGTERM (the default signal for "kill") while a
    # subprocess is running, we want to be able to terminate the
    # subprocess, too, so that the system is not kept busy with
    # orphaned processes.
    # Python installs by default a signal handler for SIGINT which
    # raises the KeyboardInterrupt exception. This is convenient, as
    # it lets us simply terminate the child in an exception handler.
    # Thus we install the signal handler of SIGINT for SIGTERM as well,
    # so that SIGTERM likewise raises a KeyboardInterrupt exception.

    sigint_handler = signal.getsignal(signal.SIGINT)
    signal.signal(signal.SIGTERM, sigint_handler)

    # Wait for command to finish executing, capturing stdout and stderr
    # in output tuple
    try:
        (stdout, stderr) = child.communicate()
    except KeyboardInterrupt:
        logging.critical("[%s] KeyboardInterrupt received, killing child "
                         "process with PID %d", time.asctime(), child.pid)
        child.terminate()
        (stdout, stderr) = child.communicate()
        logging.error("[%s] Terminated command resulted in exit code %d",
                      time.asctime(), child.returncode)
        raise  # Re-raise KeyboardInterrupt to terminate cado-nfs-client.py

    # Un-install our handler and revert to the default handler
    signal.signal(signal.SIGTERM, signal.SIG_DFL)

    if print_error and child.returncode != 0:
        logging.error("Command resulted in exit code %d", child.returncode)
    return child.returncode, stdout, stderr
# }}}


class ServerGone(Exception):
    def __init__(self):
        Exception.__init__(self)

    def __str__(self):
        return "Server gone"


class NoMoreServers(Exception):
    def __init__(self):
        Exception.__init__(self)

    def __str__(self):
        return "All servers dropped the connection" \
               " (connection reset or refused)"


class Server(object):
    def __init__(self, index, url, cafile=None, certdigest=None):
        print(f"url={url}")
        self.index = index
        self.url = url
        assert (cafile is None) == (certdigest is None)
        self.cafile = cafile
        self.certdigest = certdigest
        self.enable = True
        self.worked_once = False
        (self.scheme, netloc) = urlparse(url)[:2]
        self.ip, self.port = netloc.split(':')
        self.port = int(self.port)
        self.cert_downloaded = False

    def disable(self):
        self.enable = False

    def get_url(self, ep=None):
        return f"{self.url.rstrip('/')}/{ep}" if ep else self.url

    def __str__(self):
        return self.url

    def get_cafile(self):
        return self.cafile

    def get_index(self):
        return self.index

    def get_certificate(self, retry=0):
        """
        Download the certificate if it is missing and check its fingerprint

        If the file 'self.cafile' already exists, the certificate does
        not get downloaded.

        If the certificate existed or could be downloaded and the
        fingerprint matches, returns True.

        If the fingerprint check fails, exits with error.

        If the server refuses connections and retry is False, returns False;

        If retry is a nonzero value, keeps trying indefinitely.
        """
        if not self.enable:
            return False
        if self.cert_downloaded:
            return True
        if not self.cafile:
            return True

        # get the PEM certificate
        certfile_exists = os.path.isfile(self.cafile)
        if certfile_exists:
            logging.info("Using certificate stored in file %s", self.cafile)
            with open(self.cafile, 'r') as certfile:
                cert = certfile.read()
        else:
            logging.info("Downloading certificate from %s", self.url)
            cert = None
            while True:
                """
                In case of connection refused error, if retry is True, retry
                indefinitely waiting retrytime seconds between tries, and if
                retry is False, return None.
                """
                try:
                    cert = ssl.get_server_certificate((self.ip, self.port),
                                                      ca_certs=None)
                    break
                except socket.error as err:
                    if err.errno != errno.ECONNREFUSED:
                        raise
                    if not retry:
                        return False
                wait = float(retry)
                logging.error("Waiting %s seconds before retrying", wait)
                time.sleep(wait)

        # Note: if you want the sha1 just based on the cert file, it's rather
        # easy:
        # openssl x509 -in $wdir/c60.server.cert -outform DER -out - | sha1sum
        bin_cert = ssl.PEM_cert_to_DER_cert(cert)
        h = hashlib.sha1()
        h.update(bin_cert)
        certdigest = h.hexdigest()
        logging.debug("Certificate has SHA1 fingerprint %s", certdigest)
        if not certdigest.lower() == self.certdigest.lower():
            logging.critical("Server certificate's SHA1 fingerprint (%s)"
                             " differs from fingerprint specified"
                             " on command line (%s). Aborting.",
                             certdigest, self.certdigest)
            logging.critical("Possible reason: several factorizations with "
                             "same download directory.")
            sys.exit(1)
        logging.info("Certificate SHA1 hash matches")

        self.cert_downloaded = True
        if not certfile_exists:
            logging.info("Writing certificate to file %s", self.cafile)
            with open(self.cafile, 'w') as certfile:
                certfile.write(cert)
        return True

    def _req(self, req, ep, *args, **kwargs):
        kwargs['verify'] = self.cafile if self.cafile else False
        logging.debug(f"arguments to {req} {ep}: %s", kwargs)
        resp = getattr(requests, req)(self.get_url(ep), *args, **kwargs)
        self.worked_once = True
        return resp

    def get(self, *args, **kwargs):
        return self._req('get', *args, **kwargs)

    def post(self, *args, **kwargs):
        return self._req('post', *args, **kwargs)


class ServerPool(object):  # {{{
    def __init__(self, settings):
        self.ndisabled = 0
        self.has_https = False
        self.current_index = 0
        self.wait = float(settings["DOWNLOADRETRY"])

        self.has_https = any([urlparse(ss)[0] == "https"
                              for ss in settings["SERVER"]])

        self.servers = []

        if not self.has_https:
            if settings["CERTSHA1"] is not None:
                logging.warning("Option --certsha1 makes sense only with"
                                " https URLs, ignoring it.")
            for i, ss in enumerate(settings["SERVER"]):
                self.servers.append(Server(i, ss))
            return

        if settings["CERTSHA1"] is None and not settings["NO_CN_CHECK"]:
            logging.error("no certificate hashes provided,"
                          " pass --nocncheck if this is intentional")
            sys.exit(1)

        if settings["CERTSHA1"] is None:
            logging.warning("https URLs were given"
                            " but no --certsha1 option,"
                            " NO SSL VALIDATION WILL BE PERFORMED.")
            for i, ss in enumerate(settings["SERVER"]):
                self.servers.append(Server(i, ss))
            return

        if len(settings["CERTSHA1"]) != len(settings["SERVER"]):
            logging.critical("Exactly one --certsha1 option"
                             " must be provided per server URL"
                             " (use --certsha1 None for http URLs)")
            sys.exit(1)

        for i, (ss, certsha1) in enumerate(zip(settings["SERVER"],
                                               settings["CERTSHA1"])):
            (scheme, netloc) = urlparse(ss)[0:2]
            cafile = None
            if scheme == "https":
                cafile = os.path.join(settings["DLDIR"],
                                      "server.%s.pem" % certsha1)
            self.servers.append(Server(i, ss, cafile, certsha1))

            # Try downloading the certificate once. If connection is
            # refused, proceed to daemonizing - hopefully server will
            # come up later
            if not self._try_download_certificate(i):
                logging.info("Could not download SSL certificate:"
                             " The connection was refused.")
                logging.info("Assuming the server will come up later.")
                if options.daemon:
                    logging.info("Will keep trying.")
                else:
                    logging.info("Will keep trying after daemonizing.")

    def number_of_active_servers(self):
        return len(self.servers) - self.ndisabled

    def _try_download_certificate(self, server_index):
        return self.servers[server_index].get_certificate(retry=0)

    def get_default_server(self):
        """returns an arbitrary server in the list, really. We have a
        preference towards keeping the server we've been using in the
        recent past.  At any rate, we return a server only if we
        succeeded in downloading the ssl certificate !
        """
        while not self._try_download_certificate(self.current_index):
            self.current_index = (self.current_index + 1) % len(self.servers)
        return self.servers[self.current_index]

    def mark_current_server_alive(self):
        self.servers[self.current_index].worked_once = True

    def current_server_was_successful_once(self):
        return self.servers[self.current_index].worked_once

    def change_server(self):
        """we're not happy with the current server for some reason.
        return a new one
        """
        self.current_index = (self.current_index + 1) % len(self.servers)
        while not self._try_download_certificate(self.current_index):
            self.current_index = (self.current_index + 1) % len(self.servers)
        S = self.servers[self.current_index]
        logging.error("Going to next backup server: %s", S)
        return S

    def disable_server(self, S):
        """multiple errors with this server, disable it permanently.
        Raises an exception if all servers are dead."""
        self.servers[S.get_index()].disable()
        self.ndisabled += 1
        if self.ndisabled == len(self.servers):
            raise NoMoreServers()
        if self.current_index == S.get_index():
            self.change_server()

    def get_current_server_index(self):
        return self.current_index

    def get_current_server(self):
        return self.servers[self.current_index]

    def get_unique_server(self):
        assert len(self.servers) == 1
        return self.servers[0]

# }}}


# {{{ WorkunitProcessor: this object processes once workunit, and owns
# the result files until they get collected by the server.
class WorkunitProcessor(object):
    def __init__(self, workunit, settings):
        self.settings = settings
        self.workunit = workunit

        # self.errorcode gets set if any command exits with code != 0, in
        # which case self.failedcommand is the index of the failed
        # command
        self.errorcode = 0
        self.failedcommand = None

        self.stdio = {"stdout": [], "stderr": []}
        self._answer = None

    def __str__(self):
        return "Processor for Workunit:\n%s" % super(WorkunitProcessor, self)

    def renice(self):
        os.nice(int(self.settings["NICENESS"]))

    @staticmethod
    def is_executable(filename):
        """ Test that the file exists and, if the stat object knows the
        "executable by user" constant, that it is executable
        """
        if not os.path.isfile(filename):
            return False
        if hasattr(stat, "S_IXUSR"):
            return (os.stat(filename).st_mode & stat.S_IXUSR) != 0
        else:
            # perhaps every file is executable?
            return True

    @staticmethod
    def find_binary(filename, searchpath):
        """ Given a search path (array of strings), find the directory which
        contains an executable "filename". If not found, return None.
        """
        # If filename contains any path information (e.g., "./foo"), then
        # try only filename itself, like the shell does
        if os.path.basename(filename) != filename:
            return filename if WorkunitProcessor.is_executable(filename) \
                    else None
        for trydir in searchpath:
            # An empty directory name results in tryname == filename, so it
            # will search in the current working directory, like the shell
            # PATH does
            tryname = os.path.join(trydir, filename)
            if WorkunitProcessor.is_executable(tryname):
                return tryname
        return None

    def apply_overrides(self, command):
        # to override several parameters, use:
        # --override t 1 --override bkthresh1 15000000
        if self.settings["OVERRIDE"] is None:
            return command

        mangled = []
        orig = re.split(' +', command)
        used_overrides = {}
        while orig:
            a = orig.pop(0)
            krepl = None
            for sub in self.settings["OVERRIDE"]:
                if re.match('^-{1,2}' + sub[0] + '$', a):
                    krepl = sub
                    used_overrides[sub[0]] = True
            mangled.append(a)
            if krepl is not None:
                k, repl = krepl
                oldvalue = orig.pop(0)
                logging.info("Overriding argument %s %s"
                             " by %s %s in command line"
                             " (substitution %s %s)",
                             a, oldvalue, a, repl, k, repl)
                mangled.append(repl)
        # apply the overrides even to flags which were *NOT* present in
        # the initial command line.
        for f, v in self.settings["OVERRIDE"]:
            if f in used_overrides:
                continue
            mangled.append('-' + f)
            mangled.append(v)

        return ' '.join(mangled)

    def _locate_binary_file(self, f):
        filename = f['filename']
        if self.settings["BINDIR"]:
            searchpath = self.settings["BINDIR"].split(';')
            suggest = f.get("suggest_path")
            if suggest is not None:
                searchpath += [os.path.join(x, suggest) for x in searchpath]
            binfile = self.find_binary(filename, searchpath)
            if binfile is None:
                raise Exception("Binary file %s not found,"
                                " search path is %s."
                                % (filename, searchpath))
        else:
            binfile = os.path.join(self.settings["DLDIR"], filename)
        logging.info('file %s resolved to %s' % (filename, binfile))
        return binfile

    def run_commands(self):
        if self.result_exists():
            if self.settings["KEEPOLDRESULT"]:
                return True
            self.cleanup()

        files = {}

        # To which directory do workunit files map?
        dirs = {"FILE": self.settings["DLDIR"],
                "RESULT": self.settings["WORKDIR"],
                "STDOUT": self.settings["WORKDIR"],
                "STDERR": self.settings["WORKDIR"],
                "STDIN": self.settings["WORKDIR"],
                }

        for fid, f in self.workunit.get('files', {}).items():
            fc = copy.copy(f)
            for key, d in dirs.items():
                if fid.startswith(key):
                    fc['filename'] = os.path.join(d, f['filename'])
                    break
            if fid.startswith("EXECFILE"):
                fc['filename'] = self._locate_binary_file(f)

            files[fid] = fc

        for (counter, command) in enumerate(self.workunit.get("commands", [])):
            command = command.replace("'", "")  # 21827
            command = Template(command).safe_substitute(
                {i: v['filename'] for i, v in files.items()})

            my_stdin_filename = "STDIN%d" % counter
            my_stdout_filename = "STDOUT%d" % counter
            my_stderr_filename = "STDERR%d" % counter

            # If niceness command line parameter was set, call self.renice()
            # in child process, before executing command
            if int(self.settings["NICENESS"]) > 0:
                renice_func = self.renice
            else:
                renice_func = None

            command = self.apply_overrides(command)

            stdin = None
            if my_stdin_filename in files:
                raise ValueError("I think that this code is dead"
                                 " and does not work")
                # run_command pretty much seems to ignore stdin anyway
                with open(files[my_stdin_filename], "r") as f:
                    stdin = f.read()

            if stdin is not None:
                stdin = stdin.encode()

            rc, stdout, stderr = run_command(command,
                                             stdin=stdin,
                                             preexec_fn=renice_func)

            # steal stdout/stderr, put them to files.
            out = files.get(my_stdout_filename)
            if out is not None:
                if stdout is not None:
                    with open(out['filename'], "wb") as f:
                        f.write(stdout)
                stdout = None
            if stdout:
                self.stdio["stdout"].append(stdout)

            err = files.get(my_stderr_filename)
            if err is not None:
                if stderr is not None:
                    with open(err['filename'], "wb") as f:
                        f.write(stderr)
                stderr = None
            if stderr:
                self.stdio["stderr"].append(stderr)

            if rc != 0:
                self.failedcommand = counter
                self.errorcode = rc
                return False

            logging.debug("Command exited successfully")

        return True

    def result_filepaths(self):
        for fid, f in self.workunit.get('files', {}).items():
            if not fid.startswith("RESULT"):
                continue
            yield os.path.join(self.settings["WORKDIR"], f['filename'])

    def result_exists(self):
        '''
        Check whether all result files already exist.
        returns True of False
        '''
        results = list(self.result_filepaths())
        if not results:
            # If there is no RESULT line in the workunit, always run commands
            return False

        for filepath in results:
            if not os.path.isfile(filepath):
                logging.info("Result file %s does not exist", filepath)
                return False
            logging.info("Result file %s already exists", filepath)
        logging.info("All result files already exist")
        return True

    def cleanup(self):
        '''
        Delete uploaded result files
        '''
        logging.info("Cleaning up for workunit %s", self.workunit.get_id())
        for filepath in self.result_filepaths():
            logging.info("Removing result file %s", filepath)
            try:
                os.remove(filepath)
            except OSError as err:
                # The file won't exist if the program failed too early on.
                logging.error("Could not remove file: %s", err)

    def prepare_answer(self):
        assert self._answer is None
        data = dict(WUid=self.workunit.get_id(),
                    clientid=self.settings["CLIENTID"])
        if self.errorcode:
            data["errorcode"] = self.errorcode
        if self.failedcommand:
            data["failedcommand"] = self.failedcommand

        files = {}
        fileinfo = {}

        for fid, f in self.workunit.get('files', {}).items():
            if not f.get('upload'):
                continue
            filepath = os.path.join(self.settings["WORKDIR"], f['filename'])
            logging.info("Attaching file %s [fid=%s] to upload", filepath, fid)
            basename = os.path.split(filepath)[1]
            files[basename] = open(filepath, 'rb')
            fileinfo[basename] = dict(WUid=self.workunit.get_id(),
                                      key=fid)

        for name, blobs in self.stdio.items():
            for (counter, blob) in enumerate(blobs):
                if blob:
                    logging.info("Attaching %s for command %d to upload",
                                 name, counter)
                    fid = "%s%d" % (name, counter)
                    filename = "%s.%s" % (self.workunit.get_id(), fid)
                    files[filename] = BytesIO(blob)
                    fileinfo[filename] = dict(WUid=self.workunit.get_id(),
                                              key=fid,
                                              command=counter)

        data['fileinfo'] = json.dumps(fileinfo)
        self._answer = dict(files=files, data=data)

    def get_answer(self):
        assert self._answer is not None
        return self._answer
# }}}


class WorkunitParseError(ValueError):
    """ Parsing the workunit failed """
    pass

# class WorkunitClientBrokenConnection(Exception):
#     """ Got "hard" errors several times (ECONNREFUSED or ECONNRESET) """
#     def __init__(self, url, msg):
#         self.text = "Broken connection to %s: %s" % (url, msg)
#     def __str__(self):
#         return self.text
#


class WorkunitClientHalfDownload(Exception):
    """ Timeout """
    def __init__(self, path):
        Exception.__init__(self)
        self.text = "Timed out while downloading %s" % path

    def __str__(self):
        return self.text


class WorkunitClientWrongChecksum(Exception):
    """ Checksum was wrong several times in a row """
    def __init__(self, path, peer, filesum):
        Exception.__init__(self)
        self.text = "Downloaded file %s" \
                    " from server %s has" \
                    " same wrong checksum %s again." % \
                    (path, peer, filesum)

    def __str__(self):
        return self.text


class WorkunitClientToFinish(Exception):
    """ we received a 410 (probably while attempting to download a WU) """
    def __init__(self, explanation):
        Exception.__init__(self)
        self.text = explanation

    def __str__(self):
        return self.text


class PrivateFileAlreadyExists(Exception):
    def __init__(self, oldname, newname):
        Exception.__init__(self)
        self.text = "cannot move %s to %s : destination already exists" \
                    % (oldname, newname)

    def __str__(self):
        return self.text


class WorkunitWrapper(Workunit):
    """ wraps a workunit with info on the originating server, and the
    filename where the WU is stored on the filesystem.
    """
    def __init__(self, filename, peer):
        self.wu_filename = filename
        self.peer = peer
        # may throw FileLockedException, which will be fatal
        try:
            self.wu_file = open_exclusive(self.wu_filename)
            try:
                logging.debug("Parsing workunit from file %s",
                              self.wu_filename)
                # super(Workunit, self).__init__(self.wu_file.read())
                Workunit.__init__(self, **json.load(self.wu_file))
            except Exception:
                close_exclusive(self.wu_file)
                raise
        except Exception:
            os.remove(self.wu_filename)
            raise

#     def __str__(self):
#         # normal __str__ for workunits prints the text in full. In truth,
#         # we don't need it.
#         return self.get_id()
#
    def get_peer(self):
        return self.peer

    def get_filename(self):
        return self.wu_filename

    def cleanup(self):
        close_exclusive(self.wu_file)
        os.remove(self.wu_filename)
        self.wu_file = None

    def move_to_server_private_file(self):
        """ XXX can we do this race-free ? And is it valid to do that
        while we hold a lock on self.wu_file ?
        """
        newname = "%s.%x" % (self.wu_filename, abs(hash(self.peer.get_url())))
        if os.path.isfile(newname):
            raise PrivateFileAlreadyExists(self.wu_filename, newname)
        os.rename(self.wu_filename, newname)
        self.wu_filename = newname

    def is_stale(self):
        d = self.get("deadline")
        if d is not None:
            return time.time() > float(d)
        return False


# {{{ InputDownloader -- persistent class that downloads WUs together
# with their companion files, and provides them when they're ready.
# Half-downloaded WUs are saved in memory, and downloads of companion
# files are retried later on if the peer goes off at the wrong time.
class InputDownloader(object):
    def __init__(self, settings, server_pool):
        self.settings = settings
        self.server_pool = server_pool
        self.wu_filename = os.path.join(self.settings["DLDIR"],
                                        self.settings["WU_FILENAME"])
        self.wu_backlog = []
        self.wu_backlog_alt = []
        self.exit_on_server_gone = self.settings["EXIT_ON_SERVER_GONE"]

    # {{{ download -- this goes through several steps.
    @staticmethod
    def do_checksum(filename):
        """
        Computes the SHA1 checksum for a file.
        """
        h = hashlib.sha1()  # pylint: disable=E1101

        infile = open(filename, "rb")
        FileLock.lock(infile)

        blocksize = 65536
        data = infile.read(blocksize)
        while data:
            h.update(data)
            data = infile.read(blocksize)
        FileLock.unlock(infile)
        infile.close()
        return h.hexdigest()

    def get_file(self, urlpath,
                 dlpath=None,
                 options=None,
                 is_wu=False,
                 executable=False,
                 mandatory_server=None):
        """
        gets a file from the server (of from one of the failover servers,
        for WUs), and wait until we succeed.

        returns the identification of the server that answered if we got
        an answer, or None if we didn't get one. (the latter can happen
        only if we've been told to use one server exclusively, and that
        happens only if we timed out downloading companion files).

        Raises NoMoreServers if we get multiple consecutive connection
        failures on all servers.
        """
        assert is_wu or dlpath is not None
        if dlpath is None:
            filename = urlpath.split("/")[-1]
            dlpath = os.path.join(self.settings["DLDIR"], filename)
        urlpath = urlpath.lstrip("/")

        wait_stable = float(self.settings["DOWNLOADRETRY"])

        wait_min = self.settings["DOWNLOADRETRYMIN"]
        if wait_min is None:
            wait_min = wait_stable / 16
        else:
            wait_min = float(wait_min)

        wait = wait_min

        waiting_since = 0
        # this knowingly mixes http status codes in the 400- 500- range
        # with errno errors. It's ugly.
        last_error = ""
        # record the number of connection failures
        connfailed = 0
        maxconnfailed = int(self.settings["MAX_CONNECTION_FAILURES"])
        silent_wait = self.settings["SILENT_WAIT"]

        current_server = mandatory_server
        if current_server is None:
            current_server = self.server_pool.get_default_server()

        max_loops = self.server_pool.number_of_active_servers()
        cap = is_wu and (self.wu_backlog or self.wu_backlog_alt)
        spin = 0
        if dlpath is None:
            dlpath_tmp = None
        else:
            dlpath_tmp = "%s%d" \
                         % (dlpath, random.randint(0, 2**30) ^ os.getpid())
        while True:
            logging.debug("spin=%d is_wu=%s blog=%d",
                          spin, is_wu,
                          len(self.wu_backlog)+len(self.wu_backlog_alt))
            if cap and spin > max_loops:
                # we've had enough. Out of despair, we'll try our old
                # WUs, but there seems to be veeery little we can do, to
                # be honest. We'll quickly return back here.
                logging.error("Cannot get a fresh WU. Trying our old backlog")
                return None

            logging.info("Downloading %s to %s",
                         current_server.get_url(urlpath),
                         dlpath_tmp)

            ex = None
            try:
                # stream=True / resp.raw.data is to avoid magic
                # decompression of .gz files
                resp = current_server.get(urlpath,
                                          data=options,
                                          stream=True)
                if resp.status_code == 200:
                    open(dlpath_tmp, 'wb').write(resp.raw.data)
                    break
                elif resp.status_code == 410:
                    # We interpret error code 410 as the work unit server
                    # being gone for good. This instructs us to terminate
                    # the workunit client, which we do by letting an
                    # exception pop up a few levels up (eeek)
                    raise WorkunitClientToFinish("Received 410 from server")
                error_str = resp.content.decode()
            except requests.exceptions.RequestException as e:
                ex = e
                error_str = str(e)
                connfailed += 1

            if all([ex,
                    self.exit_on_server_gone,
                    self.server_pool.current_server_was_successful_once()]):
                logging.error(f"Disabling {current_server}"
                              " because of --exit-on-server-gone")
                try:
                    self.server_pool.disable_server(current_server)
                except NoMoreServers:
                    raise ServerGone()
                spin += 1
                current_server = self.server_pool.get_current_server()
                waiting_since = 0
                last_error = ""
                connfailed = 0
                continue

            # otherwise we enter the wait loop

            if any([not silent_wait,
                    waiting_since == 0,
                    error_str != last_error]):
                logging.error("Download failed%s, %s",
                              " with hard error" if ex else "",
                              error_str)
                if waiting_since > 0:
                    logging.error("Waiting %s seconds before retrying"
                                  " (I have been waiting for %s seconds)",
                                  wait, waiting_since)
                else:
                    logging.error("Waiting %s seconds before retrying", wait)
            last_error = error_str

            if connfailed > maxconnfailed:
                logging.error("Connection to %s failed %d times",
                              current_server, connfailed)
                self.server_pool.disable_server(current_server)
                spin += 1
                current_server = self.server_pool.get_current_server()
                waiting_since = 0
                wait = wait_min
                last_error = ""
                connfailed = 0
                continue

            # 4 means that we'll try 5 times.
            if waiting_since >= 4 * wait_stable:
                if mandatory_server is None:
                    current_server = self.server_pool.change_server()
                    spin += 1
                    waiting_since = 0
                    wait = wait_min
                    last_error = ""
                    connfailed = 0
                else:
                    # we failed to download from the mandatory server,
                    # too bad.
                    current_server = self.server_pool.change_server()
                    spin += 1
                    return None
            else:
                time.sleep(wait)
                waiting_since += wait
                wait *= 2
                if wait >= wait_stable:
                    wait = wait_stable

        if waiting_since > 0:
            logging.info("Opened URL %s after %s seconds wait",
                         current_server.get_url(), waiting_since)
        wait = wait_min

        if executable:
            m = dlpath_tmp if dlpath_tmp is not None else dlpath
            mode = os.stat(m).st_mode
            if mode & stat.S_IXUSR == 0:
                logging.info("Setting executable flag for %s", dlpath)
                os.chmod(m, mode | stat.S_IXUSR)

        if dlpath_tmp is not None:
            # We can't atomically rename-unless-dst-does-not-exist-yet.
            os.rename(dlpath_tmp, dlpath)

        return current_server

    def get_missing_file(self, urlpath, filename,
                         checksum=None,
                         options=None,
                         is_wu=False,
                         executable=False,
                         mandatory_server=None
                         ):
        """ Downloads a file if it does not exist already, from one of
        the servers configured for failover.

        Also checks the checksum, if specified; if the file already exists and
        has a wrong checksum, it is deleted and downloaded anew. If the
        downloaded file has the wrong checksum, it is deleted and downloaded
        anew.

        Returns the answering server, or None if the file was found
        locally.

        In case of error, an exception is thrown. It may be:
            WorkunitClientWrongChecksum
            WorkunitClientHalfDownload
            # (no longer) WorkunitClientBrokenConnection
            NoMoreServers (fatal)
        """

        # We either have
        #   filename=the name of a WU file
        #   checksum=None
        #   is_wu=True
        #   mandatory_server=None
        # or
        #   filename=another file name
        #   checksum=something (or None if NOSHA1CHECK)
        #   is_wu=False
        #   mandatory_server=something
        assert filename is not None
        assert (is_wu or self.settings["NOSHA1CHECK"]) == (checksum is None)
        assert is_wu == (mandatory_server is None)

        # print('get_missing_file(%s, %s, %s)' % (urlpath, filename, checksum))
        if os.path.isfile(filename):
            if len(self.server_pool.servers) > 1 and is_wu:
                # can't reuse WUs on disk if multiple servers are
                # specified, unless we save the server id in the WU
                # (which we could consider, actually)
                logging.info("%s already exists,"
                             " removing because of server ambiguity",
                             filename)
                os.remove(filename)
            else:
                if checksum is None:
                    logging.info("%s already exists, not downloading",
                                 filename)
                    if is_wu:
                        return self.server_pool.get_unique_server()
                    return None
                filesum = self.do_checksum(filename)
                if filesum.lower() == checksum.lower():
                    logging.info("%s already exists, not downloading",
                                 filename)
                    if is_wu:
                        return self.server_pool.get_unique_server()
                    return None
                logging.error("Existing file %s has wrong checksum %s, "
                              "workunit specified %s. Deleting file.",
                              filename, filesum, checksum)
                os.remove(filename)

        # If checksum is wrong and does not change during two downloads, exit
        # with failue, as apparently the file on the server and checksum in
        # workunit do not agree
        last_filesum = None
        while True:
            # we were catching HTTPError here previously. Useless now ?
            peer = self.get_file(urlpath, filename, options=options,
                                 is_wu=is_wu,
                                 executable=executable,
                                 mandatory_server=mandatory_server)
            if peer is None:
                if is_wu:
                    return None
                assert mandatory_server is not None
                if os.path.isfile(filename):
                    logging.error("Removing %s since download failed",
                                  filename)
                    os.remove(filename)
                raise WorkunitClientHalfDownload(filename)
            if checksum is None:
                return peer
            filesum = self.do_checksum(filename)
            if filesum.lower() == checksum.lower():
                return peer
            os.remove(filename)
            if last_filesum is not None and filesum == last_filesum:
                u = peer.get_url()
                raise WorkunitClientWrongChecksum(filename, u, filesum)
            logging.error("Downloaded file %s has wrong checksum %s, "
                          "workunit from %s specified %s. Deleting file.",
                          filename, filesum,
                          mandatory_server,
                          checksum)
            last_filesum = filesum
        # never reach here

    def get_files(self, wu):
        server = wu.get_peer()
        for fid, f in wu.get('files', {}).items():
            if not f.get('download'):
                continue
            templ = Template(f['filename'])
            checksum = f.get('checksum')
            urlpath = templ.safe_substitute({"ARCH": self.settings["ARCH"]})
            urlpath = os.path.join('file', urlpath)
            dlname = templ.safe_substitute({"ARCH": ""})
            dlpath = os.path.join(self.settings["DLDIR"], dlname)
            if self.settings["NOSHA1CHECK"]:
                checksum = None
            # If we fail to download the file, we'll deal with it at the
            # level above

            self.get_missing_file(urlpath, dlpath, checksum,
                                  executable=fid.startswith('EXECFILE'),
                                  mandatory_server=server)
            # Try to lock the file once to be sure that download has finished
            # if another cado-nfs-client is doing the downloading
            with open(dlpath) as file_to_lock:
                FileLock.lock(file_to_lock)
                FileLock.unlock(file_to_lock)
        return True

    # }}}

    def _get_fresh_wu(self):
        """
        returns a WorkunitWrapper.  We don't know whether the companion
        files are present at this point. There is still a possibility
        that we timeout while downloading them.

        may throw WorkunitParseError
        """

        while True:
            # Download the WU file if none exists
            url = self.settings["GETWUPATH"]
            # we could maybe add more options, like architecture, qrange
            # size, whatnot.
            options = dict(clientid=self.settings["CLIENTID"])

            # will not throw, or maybe NoMoreServers
            peer = self.get_missing_file(url,
                                         self.wu_filename,
                                         is_wu=True,
                                         options=options)

            if peer is None:
                # could not download a WU...
                return None

            # Get an exclusive lock to avoid two clients working on the same
            # workunit
            try:
                real_peer = peer
                workunit = WorkunitWrapper(self.wu_filename, real_peer)
            except FileLockedException:
                # this one is fatal.
                logging.error("File '%s' is already locked. This may "
                              "indicate that two clients "
                              "with clientid '%s' are "
                              "running. Terminating.",
                              self.wu_filename, self.settings["CLIENTID"])
                raise
            # otherwise the WU constructor itself failed.
            except Exception as err:
                logging.error("Invalid workunit file: %s", err)
                raise WorkunitParseError()

            # Don't do deadline checks on WUs we received from a server.
            if peer is None and workunit.is_stale():
                dline = workunit.get("deadline")
                dline = time.asctime(time.localtime(float(dline)))
                logging.warning("Old workunit %s has passed deadline (%s),"
                                " ignoring",
                                workunit.get_id(), dline)
                workunit.cleanup()
            else:
                break

        logging.debug("Workunit ID is %s (downloaded from %s)",
                      workunit.get_id(), workunit.get_peer())

        return workunit

    def _get_wu(self):
        if self.wu_backlog:
            logging.info("Current backlog of half-downloaded WUs: %s",
                         ", ".join([w.get_id() for w in self.wu_backlog]))
        while self.wu_backlog:
            workunit = self.wu_backlog[0]
            self.wu_backlog = self.wu_backlog[1:]
            if workunit.is_stale():
                dline = workunit.get("deadline")
                dline = time.asctime(time.localtime(float(dline)))
                logging.warning("Old workunit %s has passed deadline (%s),"
                                " ignoring",
                                workunit.get_id(), dline)
                workunit.cleanup()
            else:
                logging.info("Re-attempting previously downloaded workunit %s",
                             workunit.get_id())
                return workunit
        return self._get_fresh_wu()

    def get_wu_full(self):
        """
        returns a workunit object, together with the identification of
        the origin server, and an exclusive file handle on the WU.
        All companion files of the workunit object are
        guaranteed to be there when this function returns.
        """

        self.wu_backlog_alt = []
        while True:
            workunit = self._get_wu()
            if workunit is None:
                self.wu_backlog += self.wu_backlog_alt
                self.wu_backlog_alt = []
                continue
            try:
                self.get_files(workunit)
                break
            except WorkunitClientWrongChecksum as ex:
                # discard the WU
                workunit.cleanup()
                logging.error(ex)
            except WorkunitClientHalfDownload as ex:
                logging.error(ex)
                try:
                    if workunit.get_filename() == self.wu_filename:
                        workunit.move_to_server_private_file()
                    # Important: do not append right now, or we're in for
                    # an infinite loop
                    self.wu_backlog_alt.append(workunit)
                except Exception as err:
                    logging.error("Cannot stow workunit %s: %s ; discarding\n",
                                  workunit.get_id(), err)
                    workunit.cleanup()
        self.wu_backlog += self.wu_backlog_alt
        self.wu_backlog_alt = []
        return workunit
# }}}


# {{{ ResultUploader -- persistent class that handles uploads
class ResultUploader(object):
    def __init__(self, settings, server_pool):
        self.settings = settings
        self.server_pool = server_pool
        self.upload_backlog = [[] for S in self.server_pool.servers]
        self.backlog_size = 0
        self.last_active = 0

    def schedule_upload(self, processor):
        """ takes a WorkunitProcessor object
        """
        idx = processor.workunit.get_peer().get_index()
        self.upload_backlog[idx].append(processor)
        self.backlog_size += 1
        self.last_active = idx

    @staticmethod
    def get_content_charset(conn):
        """ Returns the character set of the server's response.

        Defaults to latin-1 if no charset header was sent.
        The encoding may matter if the path names for the uploaded files
        contain special characters, like accents.
        """
        charset = conn.info().get_content_charset()
        if charset is None:
            charset = "latin-1"
        return charset

    def process_uploads(self):
        if self.backlog_size == 0:
            return

        logging.info("Uploading results: backlog has %d entries",
                     self.backlog_size)
        for u in self.upload_backlog:
            for p in u:
                logging.info("\t%s  -->  %s",
                             p.workunit.get_id(),
                             p.workunit.get_peer())
        mention_each = self.backlog_size > 1

        wait_stable = float(self.settings["DOWNLOADRETRY"])

        wait_min = self.settings["DOWNLOADRETRYMIN"]
        if wait_min is None:
            wait_min = wait_stable / 16
        else:
            wait_min = float(wait_min)

        wait = wait_min

        # Now try to purge our backlog, starting from the server we've
        # just got this WU from.
        # did_progress = False

        for i in range(len(self.upload_backlog)):
            index = (self.last_active + i) % len(self.upload_backlog)
            old_backlog = self.upload_backlog[index]
            new_backlog = []
            for p in old_backlog:
                if p.workunit.is_stale():
                    logging.error("Workunit %s has expired, not uploading",
                                  p.workunit.get_id())
                    p.cleanup()
                    continue
                if mention_each:
                    logging.info("Attempt: %s  -->  %s",
                                 p.workunit.get_id(),
                                 p.workunit.get_peer())
                resp = None
                waiting_since = 0
                while True:
                    e = None
                    try:
                        resp = p.workunit.get_peer().post(
                                self.settings["POSTRESULTPATH"],
                                **p.get_answer())
                        if resp.status_code == 200:
                            break
                        e = f'status code {resp.status_code}, {resp.content}'
                        resp = None
                    except Exception as ex:
                        e = ex

                    logging.error("Upload failed, %s", e)
                    if waiting_since > 0:
                        logging.error("Waiting %s seconds before retrying"
                                      " (I have been waiting for %s seconds)",
                                      wait, waiting_since)
                    else:
                        logging.error("Waiting %s seconds before retrying",
                                      wait)
                    time.sleep(wait)
                    waiting_since += wait
                    wait *= 2
                    if wait >= wait_stable:
                        wait = wait_stable
                    if waiting_since >= 4 * wait_stable:
                        logging.error("Giving up on this upload,"
                                      " will retry later")
                        new_backlog.append(p)
                        wait = wait_min
                        break
                if not resp:
                    # we've appended p to the new backlog at this point.
                    continue
                if p.workunit.is_stale():
                    # maybe it went stale since last time we checked.
                    logging.error("Workunit %s has expired, not uploading",
                                  p.workunit.get_id())
                    p.cleanup()
                    continue
                if waiting_since > 0:
                    logging.info("Opened URL %s after %s seconds wait",
                                 p.workunit.get_peer().get_url(
                                     self.settings["POSTRESULTPATH"]),
                                 waiting_since)
                wait = wait_min
                logging.debug("Server response:\n%s", resp.content)
                # did_progress = True
                logging.info("Upload of %s succeeded.", p.workunit.get_id())
                p.cleanup()
            self.backlog_size -= len(old_backlog) - len(new_backlog)
            self.upload_backlog[index] = new_backlog

        # We could return "did_progress", but alas, returning False will
        # cause the client to exit. We'd rather arrange to have our
        # backlog grow.
# }}}


# {{{ WorkunitClient -- gets one WU, runs it, schedules its upload.
class WorkunitClient(object):
    def __init__(self, settings, server_pool, downloader, uploader):
        self.server_pool = server_pool
        self.downloader = downloader
        self.uploader = uploader
        self.settings = settings
        self.workunit = None
        self.exit_on_server_gone = settings["EXIT_ON_SERVER_GONE"]

    def have_terminate_request(self):
        return self.workunit.get("TERMINATE", None) is not None

    def process(self):
        self.workunit = downloader.get_wu_full()

        if self.have_terminate_request():
            self.workunit.cleanup()
            logging.info("Received TERMINATE, exiting")
            return False

        processor = WorkunitProcessor(self.workunit, self.settings)
        # TODO: a select loop, maybe ? We could handle the pending
        # uploads, maybe.
        processor.run_commands()
        # to check the return value of the above command, change the previous
        # line to ret = processor.run_commands(), and there is an error if
        # ret is false (if not ret).
        # Then we can search for a particular string in stderr as follows:
        # ret = processor.run_commands()
        # if not ret and re.search("xyx", str(processor.stdio["stderr"][0])):
        #    output_something_to_some_log_file
        #    sys.exit(1)
        # this is useful if a given error always happens on a given machine

        processor.prepare_answer()

        # Transfer ownership of "processor" to the schedule_upload
        # control flow.
        self.uploader.schedule_upload(processor)

        # And trigger the upload code now.
        self.uploader.process_uploads()

        self.workunit.cleanup()
        self.workunit = None

        return True
# }}}


# Settings which we require on the command line (no defaults)
REQUIRED_SETTINGS = {
    "SERVER": (None,
               "Base URL for WU server."
               " Can be specified multiple times for failover")
    }


# Optional settings with defaults, overrideable on command line,
# and a help text
OPTIONAL_SETTINGS = {
    "WU_FILENAME":     (None, "Filename under which to store WU files"),
    "CLIENTID":        (None,
                        "Unique ID for this client. If not "
                        "specified, a default of "
                        "<hostname>.<random hex number> is used"),
    "DLDIR":           ('download/', "Directory for downloading files"),
    "WORKDIR":         (None, "Directory for result files"),
    "BINDIR":          (None,
                        "Directory with existing executable "
                        "files to use"),
    "BASEPATH":        (None,
                        "Base directory for"
                        " download and work directories"),
    "GETWUPATH":       ("/workunit",
                        "Path segment of URL for"
                        " requesting WUs from server"),
    "POSTRESULTPATH":  ("/upload",
                        "Path segment of URL for"
                        " reporting results to server"),
    "DEBUG":           ("0", "Debugging verbosity"),
    "ARCH":            ("", "Architecture string for this client"),
    "DOWNLOADRETRY":   ("10", "Time to wait before download retries"),
    "DOWNLOADRETRYMIN": (None,
                         "Min time between download retries"
                         " (exponentially increases to DOWNOADRETRY"),
    "CERTSHA1":        (None,
                        "SHA1 of server SSL certificate."
                        " Specify multiple times for failover servers"),
    "SILENT_WAIT":     (None,
                        "Discard repeated messages about"
                        " client waiting for work (does not affect uploads)"),
    "MAX_CONNECTION_FAILURES":  ("999999",
                                 "Maximum number of successive"
                                 " connection failures to tolerate"),
    "NICENESS":        ("0", "Run subprocesses under this niceness"),
    "LOGLEVEL":        ("INFO", "Verbosity of logging"),
    "LOGFILE":         (None,
                        "File to which to write log output. "
                        "In daemon mode, if no file is specified, a "
                        "default of <workdir>/<clientid>.log is used")
    }


# Merge the two, removing help string
def merge_two_dicts(x, y):
    z = x.copy()
    z.update(y)
    return z


SETTINGS = {a: b for (a, (b, c)) in
            merge_two_dicts(REQUIRED_SETTINGS, OPTIONAL_SETTINGS).items()}

BAD_WU_MAX = 3  # Maximum allowed number of bad WUs


if __name__ == '__main__':

    def parse_cmdline():
        # Create command line parser from the keys in SETTINGS
        parser = optparse.OptionParser()
        for (arg, default) in REQUIRED_SETTINGS.items():
            if arg == "SERVER":
                parser.add_option('--' + arg.lower(), help=default[1],
                                  action='append')
            else:
                parser.add_option('--' + arg.lower(), help=default[1])
        for (arg, default) in OPTIONAL_SETTINGS.items():
            if default[0] is not None:
                parser.add_option('--' + arg.lower(),
                                  default=default[0],
                                  help=f"{default[1]} (default: {default[0]})")
            elif arg == "CERTSHA1":
                parser.add_option('--' + arg.lower(), help=default[1],
                                  action='append')
            else:
                parser.add_option('--' + arg.lower(), help=default[1])
        parser.add_option("-d", "--daemon",
                          action="store_true", dest="daemon",
                          help="Daemonize the client")
        parser.add_option("--ping", type="int", dest="ping",
                          help="Checks health of existing client."
                               " Requires clientid")
        parser.add_option("--keepoldresult",
                          default=False, action="store_true",
                          help="Keep and upload old results"
                               " when client starts")
        parser.add_option("--nosha1check",
                          default=False, action="store_true",
                          help="Skip checking the SHA1 for input files")
        parser.add_option("--single", default=False, action="store_true",
                          help="process only a single WU, then exit")
        parser.add_option("--exit-on-server-gone",
                          default=False, action="store_true",
                          help="when the server is gone"
                               " (but was at least seen present once), exit")
        parser.add_option("--nocncheck", default=False, action="store_true",
                          help="Don't check common name/SAN of certificate. ")
        parser.add_option("--override", nargs=2, action='append',
                          metavar=('REGEXP', 'VALUE'),
                          help="Modify command-line arguments which match "
                               "^-{1,2}REGEXP$ to take the given VALUE."
                               " Note that REGEXP cannot start with a dash")
        parser.add_option("--logdate", default=True, action='store_true',
                          help="Include ISO8601 format date in logging")

        # Parse command line
        (options, args) = parser.parse_args()

        if args:
            sys.stderr.write("Did not understand"
                             " command line arguments %s\n" % " ".join(args))
            raise Exception()
        # Copy values to SETTINGS
        for arg, value in SETTINGS.items():
            if hasattr(options, arg.lower()):
                SETTINGS[arg] = getattr(options, arg.lower())
        for arg, value in REQUIRED_SETTINGS.items():
            if SETTINGS[arg] is None:
                raise Exception("Command line parameter --%s is required"
                                % arg.lower())
        return options

    def makedirs(path, mode=None, exist_ok=False):
        # Python 3.2 os.makedirs() has exist_ok, but older Python do not
        if sys.version_info[0:2] >= (3, 2):
            if mode is None:
                os.makedirs(path, exist_ok=exist_ok)
            else:
                os.makedirs(path, mode=mode, exist_ok=exist_ok)
        else:
            try:
                if mode is None:
                    os.makedirs(path)
                else:
                    os.makedirs(path, mode=mode)
            except OSError as e:
                if e.errno == errno.EEXIST and exist_ok:
                    pass
                else:
                    raise

    options = parse_cmdline()

    if options.ping is not None:
        if SETTINGS["CLIENTID"] is None:
            raise ValueError("--ping requires --clientid")
        if not options.daemon and SETTINGS["LOGFILE"] is None:
            raise ValueError("--ping requires --daemon or --logfile")
    # If no client id is given, we use <hostname>.<randomstr>
    if SETTINGS["CLIENTID"] is None:
        hostname = socket.gethostname()
        random.seed()
        random_str = hex(random.randrange(0, 2**32)).strip('0x')
        SETTINGS["CLIENTID"] = "%s.%s" % (hostname, random_str)

    # If no working directory is given, we use <clientid>.work/
    if SETTINGS["WORKDIR"] is None:
        SETTINGS["WORKDIR"] = SETTINGS["CLIENTID"] + '.work/'
    if SETTINGS["BASEPATH"] is not None:
        SETTINGS["WORKDIR"] = os.path.join(SETTINGS["BASEPATH"],
                                           SETTINGS["WORKDIR"])
        SETTINGS["DLDIR"] = os.path.join(SETTINGS["BASEPATH"],
                                         SETTINGS["DLDIR"])
    # If no WU filename is given, we use "WU." + client id
    if SETTINGS["WU_FILENAME"] is None:
        SETTINGS["WU_FILENAME"] = "WU." + SETTINGS["CLIENTID"]

    SETTINGS["KEEPOLDRESULT"] = options.keepoldresult
    SETTINGS["NOSHA1CHECK"] = options.nosha1check
    SETTINGS["NO_CN_CHECK"] = options.nocncheck

    if options.nocncheck:
        from urllib3.exceptions import InsecureRequestWarning
        from requests.packages.urllib3 import disable_warnings
        disable_warnings(category=InsecureRequestWarning)

    SETTINGS["OVERRIDE"] = options.override
    SETTINGS["SINGLE"] = options.single
    SETTINGS["EXIT_ON_SERVER_GONE"] = options.exit_on_server_gone

    # Create download and working directories if they don't exist
    if not os.path.isdir(SETTINGS["DLDIR"]):
        makedirs(SETTINGS["DLDIR"], exist_ok=True)
    if not os.path.isdir(SETTINGS["WORKDIR"]):
        makedirs(SETTINGS["WORKDIR"], exist_ok=True)

    # print(str(SETTINGS))

    loglevel = getattr(logging, SETTINGS["LOGLEVEL"].upper(), None)
    if not isinstance(loglevel, int):
        raise ValueError('Invalid log level: ' + SETTINGS["LOGLEVEL"])
    logfilename = SETTINGS["LOGFILE"]
    if options.daemon and logfilename is None:
        logfilename = "%s/%s.log" % (SETTINGS["WORKDIR"], SETTINGS["CLIENTID"])
        SETTINGS["LOGFILE"] = logfilename

    if options.ping is not None:
        if pid_exists(options.ping):
            sys.exit(0)
        with open(logfilename, "r") as f:
            size = os.stat(f.fileno()).st_size
            if size >= 8192:
                f.seek(size-8192, io.SEEK_SET)
            lines = f.readlines()
            for ell in lines[-20:]:
                sys.stderr.write("CLIENT ERROR: " + ell)
        sys.exit(1)

    logfile = None if logfilename is None else open(logfilename, "a")
    if options.daemon:
        logging.basicConfig(level=loglevel)
    else:
        logging.getLogger().addHandler(cadologger.ScreenHandler(lvl=loglevel))
    if options.logdate:
        logging.basicConfig(
            format='%(asctime)s - %(levelname)s:%(name)s:%(message)s',
            level=loglevel)
    else:
        logging.basicConfig(level=loglevel)

    if logfile:
        logging.getLogger().addHandler(logging.StreamHandler(logfile))

    logging.info("Starting client %s", SETTINGS["CLIENTID"])
    logging.info("Python version is %d.%d.%d", *sys.version_info[0:3])
    logging.info("Args %s", sys.argv)

    serv_pool = ServerPool(SETTINGS)

    if options.daemon:
        # in fact, logfile can never be None, since we force a logfile no
        # matter what.
        create_daemon(logfile=logfile)

    # main control loop.
    client_ok = True
    bad_wu_counter = 0

    downloader = InputDownloader(SETTINGS, serv_pool)

    uploader = ResultUploader(SETTINGS, serv_pool)

    client = WorkunitClient(SETTINGS, serv_pool, downloader, uploader)

    while client_ok:
        try:
            client_ok = client.process()
        except WorkunitParseError:
            bad_wu_counter += 1
            if bad_wu_counter > BAD_WU_MAX:
                logging.critical("Had %d bad workunit files. Aborting.",
                                 bad_wu_counter)
                break
            continue
        except NoMoreServers as e:
            logging.error(e)
            sys.exit(1)
        except WorkunitClientToFinish as e:
            logging.info("Client finishing: %s. Bye.", e)
            sys.exit(0)
        except ServerGone:
            sys.exit(0)

        if options.single:
            logging.info("Client processed its WU."
                         " Finishing now as implied by --single")
            sys.exit(0)
