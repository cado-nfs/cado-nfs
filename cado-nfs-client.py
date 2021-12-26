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
import random
import errno
import stat
import optparse
import shutil
import time
import subprocess
import hashlib
import logging
import socket
import signal
import re
import base64
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import email.encoders
import email.generator
from string import Template
from io import BytesIO
import ssl
import urllib.request as urllib_request
import urllib.error as urllib_error
from http.client import BadStatusLine
from urllib.parse import urlparse

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

import subprocess
import locale

pathdict=dict()

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
    example = os.path.join(install_tree, "@LIBSUFFIX@", one_pyfile_example_subpath)
    t = os.path.exists(example)
    if not t:
        if os.environ.get("CADO_NFS_DEBUG_PATHDETECT"):
            print("{} does not exist".format(example))
        return False

    # make all this relocatable, it doesn't cost us much.
    # (note though that the rpaths in the binaries are likely to still
    # contain absolute paths)
    pathdict["pylib"] = os.path.join(install_tree, "@LIBSUFFIX@/scripts/cadofactor")
    pathdict["data"]  = os.path.join(install_tree, "@DATASUFFIX@")
    pathdict["lib"]   = os.path.join(install_tree, "@LIBSUFFIX@")
    pathdict["bin"]   = os.path.join(install_tree, "@BINSUFFIX@")

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

    pathdict["pylib"] = os.path.join(source_tree, "scripts/cadofactor")
    pathdict["data"] = os.path.join(source_tree, "parameters")
    pathdict["lib"] = mydir
    pathdict["bin"] = mydir

    if os.environ.get("CADO_NFS_DEBUG_PATHDETECT"):
        print("cado-nfs running in build tree")
    return True

def detect_source_tree(pathdict):
    mydir = os.path.normpath(os.path.dirname(sys.argv[0]))
    t = os.path.exists(os.path.join(mydir, one_pyfile_example_subpath))
    helper = os.path.join(mydir, "scripts/build_environment.sh")
    if not os.path.exists(helper):
        if os.environ.get("CADO_NFS_DEBUG_PATHDETECT"):
            print("{} does not exist".format(helper))
        return False
    pipe = subprocess.Popen([helper, "--show"], stdout=subprocess.PIPE)
    loc = locale.getdefaultlocale()[1]
    if not loc:
        loc="ascii"
    output = pipe.communicate()[0].decode(loc)
    cado_bin_path = [x.split("=",2)[1] for x in output.split("\n") if re.match("^build_tree",x)][0]
    cado_bin_path = re.sub("^\"(.*)\"$", "\\1", cado_bin_path)

    pathdict["pylib"] = os.path.join(mydir, "scripts/cadofactor")
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
    raise RuntimeError("We're unable to determine the location of the cado-nfs binaries and python files")

sys.path.append(pathdict["pylib"])

# END OF THE PART THAT MUST BE EXACTLY IDENTICAL IN cado-nfs.py and cado-nfs-client.py

from workunit import Workunit
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

# Now for 150+ lines of anger.{{{
#
# In Python 3.0, 3.1, 3.2.x < 3.2.4, 3.3.x < 3.3.1, use a fixed BytesGenerator
# which accepts a bytes input. The fact that the BytesGenerator in these Python
# versions doesn't is a bug, see http://bugs.python.org/issue16564
#
# Update: the first bugfix committed in that bugtracker and shipped in Python
# versions 3.2.4, 3.2.5, 3.3.2, ... is still buggy, see
# http://bugs.python.org/issue19003
# and we have to use a different work-around...
#
# Update:
# https://gitlab.inria.fr/cado-nfs/cado-nfs/-/issues/21408
#
#
# Rather than keep a list of faulty versions, we'll try to auto-detect
# them at startup.
#
# For the record, here are the version where bugfix 1 is alright.
# (3,0,0), (3,0,1),
# (3,1,0), (3,1,1), (3,1,2), (3,1,3), (3,1,4), (3,1,5),
# (3,2,0), (3,2,1), (3,2,2), (3,2,3),
# (3,3,0)

candidates_for_BytesGenerator = []

if True:
    candidates_for_BytesGenerator.append(email.generator.BytesGenerator)
    class Version1FixedBytesGenerator(email.generator.BytesGenerator):
        # pylint: disable=W0232
        # pylint: disable=E1101
        # pylint: disable=E1102
        # pylint: disable=E1002
        def _handle_bytes(self, msg):
            payload = msg.get_payload()
            if payload is None:
                return
            if isinstance(payload, bytes):
                # Payload is bytes, output is bytes - just write them
                self._fp.write(payload)
            elif isinstance(payload, str):
                super(Version1FixedBytesGenerator, self)._handle_text(msg)
            else:
                # Payload is neither bytes nor string - this can't be right
                raise TypeError('bytes payload expected: %s' % type(payload))
        _writeBody = _handle_bytes
    candidates_for_BytesGenerator.append(Version1FixedBytesGenerator)

    if tuple(sys.version_info)[0:2] == (3, 2):
        from email.message import _has_surrogates
    else:
        from email.utils import _has_surrogates

    fcre = re.compile(r'^From ', re.MULTILINE)
    class Version2FixedBytesGenerator(email.generator.BytesGenerator):
        # pylint: disable=W0232
        # pylint: disable=E1101
        # pylint: disable=E1102
        # pylint: disable=E1002
        def _handle_application(self, msg):
            # If the string has surrogates the original source was bytes,
            # so just write it back out.

            # Python 3.2 does not have the policy attribute; we use the
            # fixed generator in this case
            cte_is_7bit = getattr(self, "policy.cte_type", None) == '7bit'
            if msg._payload is None:
                return
            if _has_surrogates(msg._payload) and not cte_is_7bit:
                if self._mangle_from_:
                    msg._payload = fcre.sub(">From ", msg._payload)
                # DON'T use _write_lines() here as that mangles data
                self.write(msg._payload)
            else:
                super()._handle_text(msg)
    candidates_for_BytesGenerator.append(Version2FixedBytesGenerator)

    class Version3FixedBytesGenerator(email.generator.BytesGenerator):
        # pylint: disable=W0232
        # pylint: disable=E1101
        # pylint: disable=E1102
        # pylint: disable=E1002
        def _handle_application(self, msg):
            # If the string has surrogates the original source was bytes,
            # so just write it back out.

            # Python 3.2 does not have the policy attribute; we use the
            # fixed generator in this case
            cte_is_7bit = getattr(self, "policy.cte_type", None) == '7bit'
            if msg._payload is None:
                return
            if not cte_is_7bit:
                if self._mangle_from_:
                    msg._payload = fcre.sub(">From ", msg._payload)
                # DON'T use _write_lines() here as that mangles data
                self.write(msg._payload)
            else:
                super()._handle_text(msg)
    candidates_for_BytesGenerator.append(Version3FixedBytesGenerator)


def find_working_bytesgenerator():
    """ Return a working bytesgenerator if we have one.
    Otherwise return None.
    """

    # We use several test strings.
    # 32-byte string, so that it displays nicely in hex dumps.
    test_strings = []
    test_bytes = b'MATCH_ME\x0d\x0a--\x0d#\x0a*\xc0\x0b\xaa\xab0123456789ab'
    test_strings.append(test_bytes + test_bytes)
    test_bytes = b'MATCH_ME\x0d\x0a--\x0b#\x0b*'
    test_strings.append(test_bytes + test_bytes)

    wrong = []
    regexp = re.compile(b'MATCH_ME')
    for byte_generator in candidates_for_BytesGenerator:
        # print("Testing with %s" % str(byte_generator))
        def test_one_string(test_bytes, byte_generator):
            enc = email.encoders.encode_noop
            msg = MIMEApplication(test_bytes, _encoder=enc)
            s = BytesIO()
            g = byte_generator(s)
            g.flatten(msg)
            wireform = s.getvalue()
            msg2 = email.message_from_bytes(wireform)
            postdata = msg2.get_payload(decode=True)
            # At this point test_bytes should be a substring of postdata
            s = re.search(regexp, postdata)
            if not s:
                return False, postdata
            if s.start() + len(test_bytes) > len(postdata):
                return False, postdata
            if postdata[s.start() : s.start() + len(test_bytes)] != test_bytes:
                return False, postdata
            return True, None

        def test_all_strings(byte_generator):
            for test_bytes in test_strings:
                t, v = test_one_string(test_bytes, byte_generator)
                if not t:
                    wrong.append((byte_generator, test_bytes, v))
                    return False
            return True

        if test_all_strings(byte_generator):
            # print("Found working encoder: %s" % str(byte_generator))
            return byte_generator

    logging.error("None of our byte generators work")
    logging.error("See bug #21408")
    logging.error("https://gitlab.inria.fr/cado-nfs/cado-nfs/-/issues/21408")
    for gtp in wrong:
        byte_generator, test_bytes, postdata = gtp
        logging.error("Example of a failing test with %s:", byte_generator)
        logging.error("Original payload")
        info = base64.b64encode(test_bytes)
        info = [info[i:i+70] for i in range(0, len(info), 70)]
        for x in info:
            logging.error(x.decode('ascii'))
        logging.error("Encoded payload")
        info = base64.b64encode(postdata)
        info = [info[i:i+70] for i in range(0, len(info), 70)]
        for x in info:
            logging.error(x.decode('ascii'))
    sys.exit(1)

FixedBytesGenerator = find_working_bytesgenerator()
# }}}

def create_daemon(workdir=None, umask=None, logfile=None):# {{{
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
        # Fork a child process so the parent can exit.  This returns control to
        # the command-line or shell.  It also guarantees that the child will not
        # be a process group leader, since the child receives a new process ID
        # and inherits the parent's process group ID.  This step is required
        # to insure that the next call to os.setsid is successful.
        pid = os.fork()
    except OSError as e:
        raise Exception("%s [%d]" % (e.strerror, e.errno))

    if pid > 0:	# master
        sys.stdout.write("PID: %d\n" % pid)
        sys.stdout.flush()
        sys.exit()

    # To become the session leader of this new session and the process group
    # leader of the new process group, we call os.setsid().  The process is
    # also guaranteed not to have a controlling terminal.
    os.setsid()

    # Since the current working directory may be a mounted filesystem,
    # we avoid the issue of not being able to unmount the filesystem at
    # shutdown time by changing it to the root directory.
    if not workdir is None:
        os.chdir(workdir)

    # We probably don't want the file mode creation mask inherited from
    # the parent, so we give the child complete control over
    # permissions.
    if not umask is None:
        os.umask(umask)

    import resource		# Resource usage information.
    maxfd = resource.getrlimit(resource.RLIMIT_NOFILE)[1]
    if maxfd == resource.RLIM_INFINITY:
        maxfd = maxfd_default

    if logfile is not None:
        # must remove the intermediary handlers that the logging system
        # uses, otherwise we get inconsistent file position and python
        # gets nuts.
        logger = logging.getLogger()
        for handler in list(logger.handlers): #Remove old handlers
            logger.removeHandler(handler)

    # Iterate through and close all file descriptors.
    for fd in range(0, maxfd):
        try:
            if logfile is not None and fd != logfile.fileno():
                os.close(fd)
        except OSError:	# ERROR, fd wasn't open to begin with (ignored)
            pass

    # Redirect the standard I/O file descriptors to the specified file.  Since
    # the daemon has no controlling terminal, most daemons redirect stdin,
    # stdout, and stderr to /dev/null.  This is done to prevent side-effects
    # from reads and writes to the standard I/O file descriptors.

    # This call to open is guaranteed to return the lowest file descriptor,
    # which will be 0 (stdin), since it was closed above.
    fd0 = os.open(redirect_to, os.O_RDWR)	# standard input (0)

    fd12 = fd0
    if logfile is not None:
        fd12 = logfile.fileno()

    # Duplicate standard input to standard output and standard error.
    os.dup2(fd12, 1)			# standard output (1)
    os.dup2(fd12, 2)			# standard error (2)

    if logfile is None:
        return

    # Now re-plug the logging system to the same file descriptor as
    # stderr. we have three file descriptors open to the same file, by
    # the way. We might as well decide to do away with one of them
    # (e.g., logfile.fileno())
    logger.addHandler(logging.StreamHandler(sys.stderr))
    # os.close(logfile.fileno())
# }}}

class WuMIMEMultipart(MIMEMultipart):# {{{
    ''' Defines convenience functions for attaching files and data to a
    MIMEMultipart object
    '''

    def attach_data(self, name, filename, data, filetype=None, command=None):
        ''' Attach the data as a file

        name is the string that is sent to the server as the name of the form
        input field for the upload; for us it is always 'result'.
        filename is the string that is sent to the server as the source file
        name, this is the name as given in the RESULT lines, or some generated
        name for captured stdout/stderr.
        data is the content of the file to send.
        filetype is "RESULT" if the file to upload is specified by a RESULT
        line; "stdout" if it is captured stdout, and "stderr" if it is captured
        stderr.
        command is specified only if the data is captured stdout/stderr, and
        gives the index of the COMMAND line that produced this stdout/stderr.
        '''
        result = MIMEApplication(data, _encoder=email.encoders.encode_noop)
        result.add_header('Content-Disposition', 'form-data',
                          name=name, filename=filename)
        if not filetype is None:
            result.add_header("filetype", filetype)
        if not command is None:
            result.add_header("command", str(command))
        self.attach(result)

    def attach_file(self, name, filename, filepath, filetype=None,
                    command=None):
        ''' Attach the file as a file

        Parameters as in attach_data(), but filepath is the path to the file
        whose data should be sent
        '''
        logging.debug("Adding result file %s to upload", filepath)
        try:
            with open(filepath, "rb") as infile:
                filedata = infile.read()
        except IOError as err:
            logging.error("Could not read file %s: %s", filepath, str(err))
            return
        self.attach_data(name, filename, filedata, filetype, command)

    def attach_key(self, key, value):
        ''' Attach a simple key=value pair '''
        attachment = MIMEText(str(value))
        attachment.add_header('Content-Disposition', 'form-data',
                              name=key)
        self.attach(attachment)

    def flatten(self, debug=0):
        ''' Flatten the mimedata with BytesGenerator and return bytes array '''
        if debug >= 2:
            logging.debug("Headers of mimedata as a dictionary: %s",
                          dict(self.items()))
        bio = BytesIO()
        gen = FixedBytesGenerator(bio, mangle_from_=False)
        gen.flatten(self, unixfrom=False)
        postdata = bio.getvalue() + b"\n"
        if debug >= 2:
            logging.debug("Postdata as a bytes array: %s", postdata)
        return postdata
# }}}

# class SharedFile(object):{{{
#     def __init__(filename, mode=0o777):
#         # Try to create and open the file exclusively
#         self.filename = filename
#         flags = os.O_CREAT | os.O_RDWR | os.O_EXCL
#         try:
#             self.fd = os.open(filename, flags, mode)
#         except OSError as err:
#             if err.errno == errno.EEXIST: # If the file already existed
#                 self.existed = True
#                 self.wait_until_positive_filesize(filename)
#                 self.file = open(filename, "r+b")
#                 FileLock.lock(self.file)
#                 return
#             else:
#                 raise
#         self.existed = False
#         self.file = os.fdopen(fd, "r+b")
#         FileLock.lock(self.file, exclusive=True)
#
#     def close():
#         FileLock.unlock(self.file)
#         self.file.close() # This should also close the fd
#
#     def delete():
#         if self.existed:
#             FileLock.unlock(self.file)
#             FileLock.lock(self.file, exclusive=True)
#         try:
#             os.remove(self.filename)
#         except OSError as err:
#             if err.errno == errno.ENOENT:
#                 pass
#             else:
#                 raise
#         self.close()
#
#     def wait_until_positive_filesize(self, timeout = 60):
#         # There is a possible race condition here. If process A creates
#         # the file, then process B tries and finds that the file exists
#         # and immediately get a shared lock for reading, then process A
#         # can never get an exclusive lock for writing.
#         # To avoid this, we let process B wait until the file has
#         # positive size, which implies that process A must have the
#         # lock already. After 60 seconds, assume the file really has 0
#         # bytes and return
#         slept = 0
#         while slept < timeout and os.path.getsize(self.filename) == 0:
#             logging.warning("Sleeping until %s contains data", self.filename)
#             time.sleep(1)
#             slept += 1
#         if slept == timeout:
#             logging.warning("Slept %d seconds, %s still has no data",
#                             timeout, self.filename)
#         return
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

    logging.info ("[%s] Subprocess has PID %d", time.asctime(), child.pid)

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
        raise # Re-raise KeyboardInterrupt to terminate cado-nfs-client.py

    # Un-install our handler and revert to the default handler
    signal.signal(signal.SIGTERM, signal.SIG_DFL)

    if print_error and child.returncode != 0:
        logging.error("Command resulted in exit code %d", child.returncode)
    return child.returncode, stdout, stderr
# }}}


# {{{ wrap around python urllib
class HTTP_connector(object):
    @staticmethod
    def wait_until_positive_filesize(filename, timeout=60):
        slept = 0
        while slept < timeout and os.path.getsize(filename) == 0:
            logging.warning("Sleeping until %s contains data", filename)
            time.sleep(1)
            slept += 1
        if slept == timeout:
            logging.warning("Slept %d seconds, %s still has no data",
                            timeout, filename)

    @staticmethod
    def _urlopen_maybe_https(request, cafile=None, check_hostname=True):
        """ Treat requests for HTTPS differently depending on whether we are
        on Python 2 or Python 3.
        """
        if isinstance(request, urllib_request.Request):
            if sys.version_info[0:2] < (3, 3):
                # In Python 2, get_type() must be used to get the scheme
                scheme = request.get_type().lower()
            else:
                # The .get_type() method was deprecated in 3.3 and removed
                # in 3.4, now the scheme is stored in the .type attribute
                scheme = request.type.lower()
        else:
            # Assume it's a URL string
            scheme = request.split(":")[0].lower()
        if scheme == "https":
            # Python 3 implements HTTPS certificate checks, we can just
            # let urllib do the work for us

            context = ssl.SSLContext()
            if check_hostname:
                context.verify_mode = ssl.CERT_REQUIRED
            else:
                context.verify_mode = ssl.CERT_NONE
            context.check_hostname = bool(check_hostname)
            context.load_verify_locations(cafile=cafile)

            logging.info ("urllib_request.urlopen")
            return urllib_request.urlopen(request, context=context)

        # If we are not using HTTPS, we can just let urllib do it,
        # and there is no need for a cafile parameter (which Python 2
        # urlopen() does not accept)
        return urllib_request.urlopen(request)

    def _urlopen(self, request, cafile=None):
        """ Wrapper around urllib2.urlopen
        """
        hard_error = False
        check_hostname = not self.no_cn_check
        try:
            logging.info("_urlopen_maybe_https")
            conn = HTTP_connector._urlopen_maybe_https(request, cafile=cafile,
                                                       check_hostname=check_hostname)
            logging.info("_urlopen_maybe_https returns")
            # conn is a file-like object with additional methods:
            # geturl(), info(), getcode()
            return conn, None, None
        except urllib_error.HTTPError as error:
            current_error = error.code
            if error.code == 410:
                # We interpret error code 410 as the work unit server
                # being gone for good. This instructs us to terminate
                # the workunit client, which we do by letting an
                # exception pop up a few levels up (eeek)
                raise WorkunitClientToFinish("Received 410 from server")
            error_str = "URL error: %s" % str(error)
            hard_error = (error.errno == errno.ECONNREFUSED or \
                          error.errno == errno.ECONNRESET)
        except urllib_error.URLError as error:
            error_str = "URL error: %s" % str(error)
            current_error = error.errno
            hard_error = (error.errno == errno.ECONNREFUSED or \
                          error.errno == errno.ECONNRESET)
        except BadStatusLine as error:
            error_str = "Bad Status line: %s" % str(error)
        except socket.error as error:
            error_str = "Connection error: %s" % str(error)
        return None, error_str, hard_error

    def _native_get_file(self, url, dlpath, cafile=None, wait=None):
        # NOTE: the "wait" argument is not used by this method
        request, error_str, hard_error = self._urlopen(url, cafile=cafile)
        if request is None:
            return error_str, hard_error
        # Try to open the file exclusively
        try:
            fd = os.open(dlpath, os.O_CREAT | os.O_WRONLY | os.O_EXCL, 0o600)
        except OSError as err:
            if err.errno == 17: # File exists error
                # There is a possible race condition here. If process A creates
                # the file, then process B tries and finds that the file exists
                # and immediately get a shared lock for reading, then process A
                # can never get an exclusive lock for writing.
                # To avoid this, we let process B wait until the file has
                # positive size, which implies that process A must have the
                # lock already. After 60 seconds, assume the file really has 0
                # bytes and return.
                logging.warning("Looks like another process already created "
                                "file %s", dlpath)
                HTTP_connector.wait_until_positive_filesize(dlpath)
                return None, None
            raise
        logging.info("_native_get_file(%s) -> %s" % (url, dlpath))
        logging.info("fd = %d" % fd)
        outfile = os.fdopen(fd, "wb")
        FileLock.lock(outfile, exclusive=True)
        logging.info("lock acquired")
        shutil.copyfileobj(request, outfile)
        logging.info("copy done")
        try:
            FileLock.unlock(outfile)
            logging.info("lock release")
        except OSError as err:
            logging.info("unlock error: %s" % err)
        outfile.close() # This should also close the fd
        request.close()
        return None, None

    def __init__(self, settings):
        # can be overridden below
        self.get_file = self._native_get_file
        self.no_cn_check = settings["NO_CN_CHECK"]

# }}}



# {{{ ssl certificate stuff
def get_ssl_certificate(server, port=443, retry=False, retrytime=0):
    """ Download the SSL certificate from the server.

    In case of connection refused error, if retry is True, retry
    indefinitely waiting retrytime seconds between tries, and if
    retry is False, return None.
    """
    while True:
        try:
            cert = ssl.get_server_certificate((server, int(port)),
                                              ssl_version=ssl.PROTOCOL_TLSv1_2,
                                              ca_certs=None)
            return cert
        except socket.error as err:
            if err.errno != errno.ECONNREFUSED:
                raise
            if not retry:
                return None
        wait = float(retrytime)
        logging.error("Waiting %s seconds before retrying", wait)
        time.sleep(wait)


def get_missing_certificate(certfilename,
                            netloc,
                            fingerprint,
                            retry=False,
                            retrytime=0):
    """ Download the certificate if it is missing and check its fingerprint

    If the file 'certfilename' already exists, the certificate does not
    get downloaded.
    If the certificate existed or could be downloaded and the fingerprint
    matches, returns True. If the fingerprint check fails, exits with error.
    If the server refuses connections and retry is False, returns False;
    if retry is True, it keeps trying indefinitely.
    """
    certfile_exists = os.path.isfile(certfilename)
    if certfile_exists:
        logging.info("Using certificate stored in file %s", certfilename)
        with open(certfilename, 'r') as certfile:
            cert = certfile.read()
    else:
        logging.info("Downloading certificate from %s", netloc)
        address_port = netloc.split(":")
        cert = get_ssl_certificate(*address_port,
                                   retry=retry,
                                   retrytime=retrytime)
        if cert is None:
            return False
    # Note: if you want the sha1 just based on the cert file, it's rather
    # easy:
    # openssl x509 -in $wdir/c60.server.cert -outform DER -out - | sha1sum
    bin_cert = ssl.PEM_cert_to_DER_cert(cert)
    sha1hash = hashlib.sha1()
    sha1hash.update(bin_cert)
    cert_sha1 = sha1hash.hexdigest()
    logging.debug("Certificate has SHA1 fingerprint %s", cert_sha1)
    if not cert_sha1.lower() == fingerprint.lower():
        logging.critical("Server certificate's SHA1 fingerprint (%s) differs "
                         "from fingerprint specified on command line (%s). "
                         "Aborting.", cert_sha1, fingerprint)
        logging.critical("Possible reason: several factorizations with "
                         "same download directory.")
        sys.exit(1)
    logging.info("Certificate SHA1 hash matches")
    if not certfile_exists:
        logging.info("Writing certificate to file %s", certfilename)
        # FIXME: Set umask first?
        with open(certfilename, 'w') as certfile:
            certfile.write(cert)
    return True
# }}}

class NoMoreServers(Exception):
    def __init__(self):
        Exception.__init__(self)
    def  __str__(self):
        return "All servers dropped the connection (connection reset or refused)"

class ServerPool(object): # {{{
    class Server(object):
        def __init__(self, index, url, cafile, certsha1, needcert):
            self.index = index
            self.url = url
            self.cafile = cafile
            self.certsha1 = certsha1
            self.needcert = needcert
            self.enable = True
        def get_url(self):
            return self.url
        def __str__(self):
            return self.url
        def get_cafile(self):
            return self.cafile
        def get_index(self):
            return self.index
        @staticmethod
        def register(servers, url, cafile=None, certsha1=None, needcert=False):
            servers.append(ServerPool.Server(len(servers), url, cafile, certsha1, needcert))

    def __init__(self, settings):
        self.nservers = len(settings["SERVER"])
        self.ndisabled = 0
        self.has_https = False
        self.current_index = 0
        self.wait = float(settings["DOWNLOADRETRY"])

        for ss in settings["SERVER"]:
            scheme, netloc = urlparse(ss)[0:2]
            self.has_https = scheme == "https"

        # self.servers is a list of tuples. Each tuple contains:
        #
        # url
        # certfilename (or None)
        # sha1 (or None)
        # need_cert (boolean)
        self.servers = []

        if not self.has_https:
            if settings["CERTSHA1"] is not None:
                logging.warning("Option --certsha1 makes sense only with"
                             " https URLs,"
                             " ignoring it.")
            for ss in settings["SERVER"]:
                ServerPool.Server.register(self.servers, ss)
            return

        if settings["CERTSHA1"] is None:
            logging.warning("https URLs were given"
                         " but no --certsha1 option,"
                         " NO SSL VALIDATION WILL BE PERFORMED.")
            for ss in settings["SERVER"]:
                ServerPool.Server.register(self.servers, ss)
            return

        if len(settings["CERTSHA1"]) != len(settings["SERVER"]):
            logging.critical("Exactly one --certsha1 option"
                             " must be provided per server URL"
                             " (use --certsha1 None for http URLs)")
            sys.exit(1)

        for server_index in range(self.nservers):
            ss = settings["SERVER"][server_index]
            certsha1 = settings["CERTSHA1"][server_index]
            (scheme, netloc) = urlparse(ss)[0:2]
            cafile = None
            needcert = True
            if scheme == "https":
                cafile = os.path.join(settings["DLDIR"],
                                      "server.%s.pem" % certsha1)
            else:
                needcert = False
            ServerPool.Server.register(self.servers,
                                       ss,
                                       cafile,
                                       certsha1,
                                       needcert)
            # Try downloading the certificate once. If connection is
            # refused, proceed to daemonizing - hopefully server will
            # come up later
            if not self._try_download_certificate(server_index):
                logging.info("Could not download SSL certificate:"
                             " The connection was refused.")
                logging.info("Assuming the server will come up later.")
                if options.daemon:
                    logging.info("Will keep trying.")
                else:
                    logging.info("Will keep trying after daemonizing.")

    def number_of_active_servers(self):
        return self.nservers - self.ndisabled

    def _try_download_certificate(self, server_index):
        S = self.servers[server_index]
        if not S.enable:
            return False
        if not S.needcert:
            return True
        (scheme, netloc) = urlparse(S.get_url())[0:2]
        if get_missing_certificate(S.cafile, netloc, S.certsha1):
            self.servers[server_index].needcert = False
            return True
        logging.error("Waiting %s seconds before retrying", self.wait)
        time.sleep(self.wait)
        return False

    def get_default_server(self):
        """returns an arbitrary server in the list, really. We have a
        preference towards keeping the server we've been using in the
        recent past.  At any rate, we return a server only if we
        succeeded in downloading the ssl certificate !
        """
        while not self._try_download_certificate(self.current_index):
            self.current_index = (self.current_index + 1) % self.nservers
        return self.servers[self.current_index]

    def change_server(self):
        """we're not happy with the current server for some reason.
        return a new one
        """
        self.current_index = (self.current_index + 1) % self.nservers
        while not self._try_download_certificate(self.current_index):
            self.current_index = (self.current_index + 1) % self.nservers
        S = self.servers[self.current_index]
        logging.error("Going to next backup server: %s", S)
        return S

    def disable_server(self, S):
        """multiple errors with this server, disable it permanently.
        Raises an exception if all servers are dead."""
        self.servers[S.get_index()].enable = False
        self.ndisabled += 1
        if self.ndisabled == self.nservers:
            raise NoMoreServers()
        if self.current_index == S.get_index():
            self.change_server()

    def get_current_server(self):
        return self.servers[self.current_index]

    def get_unique_server(self):
        assert self.nservers == 1
        return self.servers[0]

#    def get_server(self, server_index):
#        url, certfilename, certsha1, needcert = self.servers[server_index]
#        return ss, certfilename
#
#    def get(self, server_index):
#        url, certfilename, certsha1, needcert = self.servers[server_index]
#        assert not needcert
#        return server_index, url, certfilename, certsha1
# }}}

# {{{ WorkunitProcessor: this object processes once workunit, and owns
# the result files until they get collected by the server.
class WorkunitProcessor(object):
    def __init__(self, workunit, settings):
        self.settings = settings
        self.origin = workunit.get_peer()
        self.workunit = workunit
        self.errorcode = 0 # will get set if any command exits with code != 0
        self.failedcommand = None # If any command exits with code != 0, this
                                  # get set to the index of the failed command
        self.stdio = {"stdout": [], "stderr": []}
        self._answer = None

    def  __str__(self):
        return "Processor for Workunit:\n%s" % super(WorkunitProcessor, self)

    def renice(self):
        os.nice(int(self.settings["NICENESS"]))

    @staticmethod
    def is_executable(filename):
        """ Test that the file exists and, if the stat object knows the
        "executable by user" constant, that it is executable
        """
        return os.path.isfile(filename) \
               and not (hasattr(stat, "S_IXUSR") \
                        and (os.stat(filename).st_mode & stat.S_IXUSR) == 0)

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
                    used_overrides[sub[0]]=True
            mangled.append(a)
            if krepl is not None:
                k,repl = krepl
                oldvalue = orig.pop(0)
                logging.info("Overriding argument %s %s"
                             " by %s %s in command line"
                             " (substitution %s %s)",
                             a, oldvalue, a, repl, k, repl)
                mangled.append(repl)
        # apply the overrides even to flags which were *NOT* present in
        # the initial command line.
        for f,v in self.settings["OVERRIDE"]:
            if f in used_overrides:
                continue
            mangled.append('-' + f)
            mangled.append(v)


        return ' '.join(mangled)

    def _locate_binary_file(self, workunit, key, filename):
        if not isinstance(filename, str):
            filename = filename[0] # Drop checksum value
        if self.settings["BINDIR"]:
            searchpath = self.settings["BINDIR"].split(';')
            suggest = workunit.get("SUGGEST_" + key, None)
            if suggest:
                searchpath += [os.path.join(x, suggest) for x in searchpath]
            binfile = self.find_binary(filename, searchpath)
            if binfile is None:
                raise Exception("Binary file %s not found" % filename)
        else:
            binfile = os.path.join(self.settings["DLDIR"], filename)
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

        for key in dirs:
            for (index, filename) in enumerate(self.workunit.get(key, [])):
                if not isinstance(filename, str):
                    filename = filename[0] # Drop checksum value
                # index is 0-based, add 1 to make FILE1, FILE2, etc. 1-based
                files["%s%d" % (key, index + 1)] = \
                        os.path.join(dirs[key], filename)

        key = "EXECFILE"
        for (index, filename) in enumerate(self.workunit.get(key, [])):
            binfile = self._locate_binary_file(self.workunit, key, filename)
            files["%s%d" % (key, index + 1)] = binfile

        for (counter, command) in enumerate(self.workunit.get("COMMAND", [])):
            command = command.replace("'", "") # 21827
            command = Template(command).safe_substitute(files)

            my_stdin_filename = "STDIN%d" % (counter+1)
            my_stdout_filename = "STDOUT%d" % (counter+1)
            my_stderr_filename = "STDERR%d" % (counter+1)

            # If niceness command line parameter was set, call self.renice()
            # in child process, before executing command
            if int(self.settings["NICENESS"]) > 0:
                renice_func = self.renice
            else:
                renice_func = None

            command = self.apply_overrides(command)

            stdin = None
            if my_stdin_filename in files:
                with open(files[my_stdin_filename], "r") as f:
                    stdin=f.read()

            if stdin is not None:
                stdin=stdin.encode()

            rc, stdout, stderr = run_command(command,
                                            stdin=stdin,
                                            preexec_fn=renice_func)

            if stdout is not None:
                stdout=stdout.decode()
            if stderr is not None:
                stderr=stderr.decode()

            # steal stdout/stderr, put them to files.
            if my_stdout_filename in files:
                if stdout is not None:
                    with open(files[my_stdout_filename], "w") as f:
                        f.write(stdout)
                stdout = None
            self.stdio["stdout"].append(stdout)

            if my_stderr_filename in files:
                if stderr is not None:
                    with open(files[my_stderr_filename], "w") as f:
                        f.write(stderr)
                stderr = None
            self.stdio["stderr"].append(stderr)

            if rc != 0:
                self.failedcommand = counter
                self.errorcode = rc
                return False

            logging.debug("Command exited successfully")

        return True

    def result_exists(self):
        ''' Check whether all result files already exist.
            returns True of False
        '''
        # If there is no RESULT line in the workunit, always run commands
        if self.workunit.get("RESULT", None) is None:
            return False
        for filename in self.workunit.get("RESULT", []):
            filepath = os.path.join(self.settings["WORKDIR"], filename)
            if not os.path.isfile(filepath):
                logging.info("Result file %s does not exist", filepath)
                return False
            logging.info("Result file %s already exists", filepath)
        logging.info("All result files already exist")
        return True

    def cleanup(self):
        ''' Delete uploaded result files and files from DELETE lines '''
        logging.info("Cleaning up for workunit %s", self.workunit.get_id())
        for filename in self.workunit.get("RESULT", []):
            filepath = os.path.join(self.settings["WORKDIR"], filename)
            logging.info("Removing result file %s", filepath)
            try:
                os.remove(filepath)
            except OSError as err:
                # The file won't exist if the program failed too early
                # on.
                logging.error("Could not remove file: %s", err)
        for filename in self.workunit.get("DELETE", []):
            filepath = os.path.join(self.settings["WORKDIR"], filename)
            logging.info("Removing file %s", filepath)
            os.remove(filepath)

    def prepare_answer(self):
        assert self._answer is None
        # Make POST data
        mimedata = WuMIMEMultipart()
        # Build a multi-part MIME document containing the WU id and result file
        mimedata.attach_key("WUid", self.workunit.get_id())
        mimedata.attach_key("clientid", self.settings["CLIENTID"])
        if self.errorcode:
            mimedata.attach_key("errorcode", self.errorcode)
        if not self.failedcommand is None:
            mimedata.attach_key("failedcommand", self.failedcommand)
        for filename in self.workunit.get("RESULT", []):
            filepath = os.path.join(self.settings["WORKDIR"], filename)
            logging.info("Attaching file %s to upload", filepath)
            mimedata.attach_file("results", filename, filepath, "RESULT")
        for name in self.stdio:
            for (counter, data) in enumerate(self.stdio[name]):
                if data:
                    logging.info("Attaching %s for command %s to upload",
                                 name, counter)
                    filename = "%s.%s%d" % (self.workunit.get_id(), name,
                                            counter)
                    mimedata.attach_data("results", filename, data, name,
                                         counter)
        postdata = mimedata.flatten(debug=int(self.settings["DEBUG"]))

        url = self.origin.get_url().rstrip("/") + "/" + \
              self.settings["POSTRESULTPATH"].lstrip("/")
        request = urllib_request.Request(url, data=postdata,
                                         headers=dict(mimedata.items()))
        self._answer = (request, self.origin.get_cafile())
        assert request.get_full_url() == url

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
#     def  __str__(self):
#         return self.text
#
class WorkunitClientHalfDownload(Exception):
    """ Timeout """
    def __init__(self, path):
        Exception.__init__(self)
        self.text = "Timed out while downloading %s" % path
    def  __str__(self):
        return self.text

class WorkunitClientWrongChecksum(Exception):
    """ Checksum was wrong several times in a row """
    def __init__(self, path, peer, filesum):
        Exception.__init__(self)
        self.text = "Downloaded file %s" \
                    " from server %s has" \
                    " same wrong checksum %s again." % \
                    (path, peer, filesum)
    def  __str__(self):
        return self.text

class WorkunitClientToFinish(Exception):
    """ we received a 410 (probably while attempting to download a WU) """
    def __init__(self, explanation):
        Exception.__init__(self)
        self.text = explanation
    def  __str__(self):
        return self.text

class PrivateFileAlreadyExists(Exception):
    def __init__(self, oldname, newname):
        Exception.__init__(self)
        self.text = "cannot move %s to %s : destination already exists" \
                    % (oldname, newname)
    def  __str__(self):
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
                Workunit.__init__(self, self.wu_file.read())
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
        d = self.get("DEADLINE")
        if d is None:
            return False
        return time.time() > float(d)

# {{{ InputDownloader -- persistent class that downloads WUs together
# with their companion files, and provides them when they're ready.
# Half-downloaded WUs are saved in memory, and downloads of companion
# files are retried later on if the peer goes off at the wrong time.
class InputDownloader(object):
    def __init__(self, settings, server_pool, connector):
        self.settings = settings
        self.server_pool = server_pool
        self.connector = connector
        self.wu_filename = os.path.join(self.settings["DLDIR"],
                                        self.settings["WU_FILENAME"])
        self.wu_backlog = []
        self.wu_backlog_alt = []

    # {{{ download -- this goes through several steps.
    @staticmethod
    def do_checksum(filename, checksum=None):
        """ Computes the SHA1 checksum for a file. If checksum is None, returns
            the computed checksum. If checksum is not None, return whether the
            computed SHA1 sum and checksum agree """
        blocksize = 65536
        sha1hash = hashlib.sha1() # pylint: disable=E1101
        # Like when downloading, we wait until the file has positive size, to
        # avoid getting the shared lock right after the other process created
        # the file but before it gets the exclusive lock
        HTTP_connector.wait_until_positive_filesize(filename)
        infile = open(filename, "rb")
        FileLock.lock(infile)

        data = infile.read(blocksize)
        while data:
            sha1hash.update(data)
            data = infile.read(blocksize)
        FileLock.unlock(infile)
        infile.close()
        filesum = sha1hash.hexdigest()
        if checksum is None:
            return filesum
        return filesum.lower() == checksum.lower()

    def get_file(self, urlpath,
            dlpath=None,
            options=None,
            is_wu=False,
            executable=False,
            mandatory_server=None):
        """ gets a file from the server (of from one of the failover
        servers, for WUs), and wait until we succeed.

        returns the identification of the server that answered if we got
        an answer, or None if we didn't get one. (the latter can happen
        only if we've been told to use one server exclusively, and that
        happens only if we timed out downloading companion files).

        Raises NoMoreServers if we get multiple
        consecutive connection failures on all servers.
        """
        assert is_wu or dlpath is not None
        if dlpath is None:
            filename = urlpath.split("/")[-1]
            dlpath = os.path.join(self.settings["DLDIR"], filename)
        urlpath = urlpath.lstrip("/")
        if options:
            urlpath = urlpath + "?" + options

        wait = float(self.settings["DOWNLOADRETRY"])
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
            dlpath_tmp = "%s%d" % (dlpath, random.randint(0,2**30)^os.getpid())
        while True:
            logging.info("spin=%d is_wu=%s blog=%d", spin, is_wu,
                    len(self.wu_backlog)+len(self.wu_backlog_alt))
            if cap and spin > max_loops:
                # we've had enough. Out of despair, we'll try our old
                # WUs, but there seems to be veeery little we can do, to
                # be honest. We'll quickly return back here.
                logging.error("Cannot get a fresh WU. Trying our old backlog")
                return None
            url = current_server.get_url().rstrip("/") + "/" + urlpath
            cafile = current_server.get_cafile()
            logging.info("Downloading %s to %s (cafile = %s)",
                         url, dlpath_tmp, cafile)
            error_str, hard_error = self.connector.get_file(url,
                                                            dlpath_tmp,
                                                            cafile=cafile,
                                                            wait=wait)
            if error_str is None:
                break
            # otherwise we enter the wait loop
            if not silent_wait or waiting_since == 0 or error_str != last_error:
                givemsg = True
            if hard_error:
                connfailed += 1
            else:
                connfailed = 0
            if givemsg:
                logging.error("Download failed%s, %s",
                              " with hard error" if hard_error else "",
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
                last_error = ""
                connfailed = 0
                continue

            # 4 means that we'll try 5 times.
            if waiting_since >= 4 * wait:
                if mandatory_server is None:
                    current_server = self.server_pool.change_server()
                    spin += 1
                    waiting_since = 0
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

        if waiting_since > 0:
            logging.info("Opened URL %s after %s seconds wait",
                         url, waiting_since)

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
            if self.server_pool.nservers > 1 and is_wu:
                # can't reuse WUs on disk if multiple servers are
                # specified.
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
        files_to_download = wu.get("FILE", [])
        if not self.settings["BINDIR"]:
            files_to_download += wu.get("EXECFILE", [])
        for (filename, checksum) in files_to_download:
            templ = Template(filename)
            archname = templ.safe_substitute({"ARCH": self.settings["ARCH"]})
            dlname = templ.safe_substitute({"ARCH": ""})
            dlpath = os.path.join(self.settings["DLDIR"], dlname)
            if self.settings["NOSHA1CHECK"]:
                checksum = None
            # If we fail to download the file, we'll deal with it at the
            # level above

            executable = os.name != "nt" and \
                    filename in dict(wu.get("EXECFILE", []))
            self.get_missing_file(archname, dlpath, checksum,
                                  executable=executable,
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
            options = "clientid=" + self.settings["CLIENTID"]
            # we could maybe add more options, like architecture, qrange
            # size, whatnot.

            # will not throw, or maybe NoMoreServers
            peer = self.get_missing_file(url,
                                         self.wu_filename,
                                         options=options,
                                         is_wu=True)

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
                dline = workunit.get("DEADLINE")
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
                dline = workunit.get("DEADLINE")
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
    def __init__(self, settings, server_pool, connector):
        self.settings = settings
        self.server_pool = server_pool
        self.connector = connector
        self.nservers = len(settings["SERVER"])
        self.upload_backlog = [[] for i in range(self.server_pool.nservers)]
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

        wait = float(self.settings["DOWNLOADRETRY"])
        # Now try to purge our backlog, starting from the server we've
        # just got this WU from.
        did_progress = False
        last_error = ""
        for i in range(self.nservers):
            index = (self.last_active + i) % self.server_pool.nservers
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
                request, cafile = p.get_answer()
                url = request.get_full_url()
                conn = None
                waiting_since = 0
                while True:
                    conn, error_str, hard_error = self.connector._urlopen(
                        request, cafile=cafile)
                    if conn:
                        break
                    logging.error("Upload failed, %s", error_str)
                    if waiting_since > 0:
                        logging.error("Waiting %s seconds before retrying"
                                      " (I have been waiting for %s seconds)",
                                      wait, waiting_since)
                    else:
                        logging.error("Waiting %s seconds before retrying",
                                      wait)
                    time.sleep(wait)
                    waiting_since += wait
                    if waiting_since >= 4 * wait:
                        logging.error("Giving up on this upload,"
                                      " will retry later")
                        new_backlog.append(p)
                        break
                if not conn:
                    # we've appended p to the new backlog at this point.
                    continue
                if p.workunit.is_stale():
                    # maybe it went stale since last time we checked.
                    logging.error("Workunit %s has expired, not uploading",
                                  p.workunit.get_id())
                    p.cleanup()
                    conn.close()
                    continue
                if waiting_since > 0:
                    logging.info("Opened URL %s after %s seconds wait",
                                 url, waiting_since)
                response = conn.read()
                encoding = self.get_content_charset(conn)
                response_str = response.decode(encoding=encoding)
                logging.debug("Server response:\n%s", response_str)
                conn.close()
                did_progress = True
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
    "SERVER" : (None,
                "Base URL for WU server."
                " Can be specified multiple times for failover")
    }

# Optional settings with defaults, overrideable on command line,
# and a help text
OPTIONAL_SETTINGS = {
    "WU_FILENAME"    : (None, "Filename under which to store WU files"),
    "CLIENTID"       : (None,
                        "Unique ID for this client. If not "
                        "specified, a default of "
                        "<hostname>.<random hex number> is used"),
    "DLDIR"          : ('download/', "Directory for downloading files"),
    "WORKDIR"        : (None, "Directory for result files"),
    "BINDIR"         : (None,
                        "Directory with existing executable "
                        "files to use"),
    "BASEPATH"       : (None,
                        "Base directory for"
                        " download and work directories"),
    "GETWUPATH"      : ("/cgi-bin/getwu",
                        "Path segment of URL for"
                        " requesting WUs from server"),
    "POSTRESULTPATH" : ("/cgi-bin/upload.py",
                        "Path segment of URL for"
                        " reporting results to server"),
    "DEBUG"          : ("0", "Debugging verbosity"),
    "ARCH"           : ("", "Architecture string for this client"),
    "DOWNLOADRETRY"  : ("10", "Time to wait before download retries"),
    "CERTSHA1"       : (None,
                        "SHA1 of server SSL certificate."
                        " Specify multiple times for failover servers"),
    "SILENT_WAIT"    : (None,
                        "Discard repeated messages about"
                        " client waiting for work (does not affect uploads)"),
    "MAX_CONNECTION_FAILURES" : ("999999",
                                 "Maximum number of successive"
                                 " connection failures to tolerate"),
    "NICENESS"       : ("0", "Run subprocesses under this niceness"),
    "LOGLEVEL"       : ("INFO", "Verbosity of logging"),
    "LOGFILE"        : (None,
                        "File to which to write log output. "
                        "In daemon mode, if no file is specified, a "
                        "default of <workdir>/<clientid>.log is used")
    }
# Merge the two, removing help string
def merge_two_dicts(x, y):
    z = x.copy()
    z.update(y)
    return z

def abort_on_python2():
    if int(sys.version_info[0]) < 3:
        logging.error("You are running cado-nfs-client with Python%d.  "
                "Python2 *used to be* supported, but no longer is.  " % int(sys.version[0]))
        sys.exit(1)

# This syntax is weird, but { a:b for [....] } won't work with python 2.6
# -- which I'm not sure we really strive to support, though.
SETTINGS = dict([(a,b) for (a, (b,c)) in
            merge_two_dicts(REQUIRED_SETTINGS, OPTIONAL_SETTINGS).items()])

BAD_WU_MAX = 3 # Maximum allowed number of bad WUs

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
            if not default[0] is None:
                parser.add_option('--' + arg.lower(),
                                  default=default[0],
                                  help=default[1] + \
                                       " (default: " + default[0] + ")")
            elif arg == "CERTSHA1":
                parser.add_option('--' + arg.lower(), help=default[1],
                                  action='append')
            else:
                parser.add_option('--' + arg.lower(), help=default[1])
        parser.add_option("-d", "--daemon", action="store_true", dest="daemon",
                          help="Daemonize the client")
        parser.add_option("--ping", type="int", dest="ping",
                          help="Checks health of existing client.  Requires clientid")
        parser.add_option("--keepoldresult", default=False, action="store_true",
                          help="Keep and upload old results when client starts")
        parser.add_option("--nosha1check", default=False, action="store_true",
                          help="Skip checking the SHA1 for input files")
        parser.add_option("--single", default=False, action="store_true",
                          help="process only a single WU, then exit")
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

    if options.ping != None:
        if SETTINGS["CLIENTID"] is None:
                raise ValueError("--ping requires --clientid")
        if not options.daemon and SETTINGS["LOGFILE"] is None:
                raise ValueError("--ping requires --daemon or --logfile")
    # If no client id is given, we use <hostname>.<randomstr>
    if SETTINGS["CLIENTID"] is None:
        import random
        hostname = socket.gethostname()
        random.seed()
        random_str = hex(random.randrange(0, 2**32)).strip('0x')
        SETTINGS["CLIENTID"] = "%s.%s" % (hostname, random_str)

    # If no working directory is given, we use <clientid>.work/
    if SETTINGS["WORKDIR"] is None:
        SETTINGS["WORKDIR"] = SETTINGS["CLIENTID"] + '.work/'
    if not SETTINGS["BASEPATH"] is None:
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
    SETTINGS["OVERRIDE"] = options.override
    SETTINGS["SINGLE"] = options.single

    # Create download and working directories if they don't exist
    if not os.path.isdir(SETTINGS["DLDIR"]):
        makedirs(SETTINGS["DLDIR"], exist_ok=True)
    if not os.path.isdir(SETTINGS["WORKDIR"]):
        makedirs(SETTINGS["WORKDIR"], exist_ok=True)

    # print (str(SETTINGS))

    loglevel = getattr(logging, SETTINGS["LOGLEVEL"].upper(), None)
    if not isinstance(loglevel, int):
        raise ValueError('Invalid log level: ' + SETTINGS["LOGLEVEL"])
    logfilename = SETTINGS["LOGFILE"]
    if options.daemon and logfilename is None:
        logfilename = "%s/%s.log" % (SETTINGS["WORKDIR"], SETTINGS["CLIENTID"])
        SETTINGS["LOGFILE"] = logfilename

    if options.ping != None:
        if pid_exists(options.ping):
            sys.exit(0)
        with open(logfilename, "r") as f:
            size = os.stat(f.fileno()).st_size
            if size >= 8192:
                f.seek(size-8192,io.SEEK_SET)
            lines=f.readlines()
            for l in lines[-20:]:
                sys.stderr.write("CLIENT ERROR: " + l)
        sys.exit(1)

    logfile = None if logfilename is None else open(logfilename, "a")
    logging.basicConfig(level=loglevel)
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

    abort_on_python2()

    if FixedBytesGenerator != candidates_for_BytesGenerator[0]:
        logging.info("Using work-around %s for buggy BytesGenerator",
                     FixedBytesGenerator)

    serv_pool = ServerPool(SETTINGS)

    connector = HTTP_connector(SETTINGS)

    if options.daemon:
        # in fact, logfile can never be None, since we force a logfile no
        # matter what.
        create_daemon(logfile = logfile)


    # main control loop.
    client_ok = True
    bad_wu_counter = 0

    downloader = InputDownloader(SETTINGS, serv_pool, connector)

    uploader = ResultUploader(SETTINGS, serv_pool, connector)

    client = WorkunitClient(SETTINGS, serv_pool, downloader, uploader)

    while client_ok:
        try:
            client_ok = client.process()
        except WorkunitParseError:
            bad_wu_counter += 1
            if bad_wu_counter > BAD_WU_MAX:
                logging.critical("Had %d bad workunit files. Aborting.", bad_wu_counter)
                break
            continue
        except NoMoreServers as e:
            logging.error(e)
            sys.exit(1)
        except WorkunitClientToFinish as e:
            logging.info("Client finishing: %s. Bye.", e)
            sys.exit(0)
        if options.single:
            logging.info("Client processed its WU."
                         " Finishing now as implied by --single")
            sys.exit(0)
