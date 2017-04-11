#!/usr/bin/env python

import sys
import os
import errno
import stat
import optparse
import shutil
import time
if sys.version_info[0] == 3:
    # pylint: disable=E0611
    # pylint: disable=F0401
    import urllib.request as urllib_request
    import urllib.error as urllib_error
    from http.client import BadStatusLine
    from urllib.parse import urlparse
elif sys.version_info[0] == 2:
    import urllib2 as urllib_request
    import urllib2 as urllib_error
    from httplib import BadStatusLine
    from urlparse import urlparse
import subprocess
import hashlib
import logging
import socket
import signal
import re
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import email.encoders
import email.generator
from string import Template
from io import BytesIO


import re

cado_python_libs_path="@CMAKE_INSTALL_PREFIX@/@LIBSUFFIX@/scripts/cadofactor"
if not re.search("^/", cado_python_libs_path):
    cado_python_libs_path=os.path.join(os.path.dirname(sys.argv[0]),"scripts","cadofactor")

# note that even though we do have cado-nfs.py and cado-nfs-client.py in
# the build tree, we make *NO PROMISE* as to whether calling these
# scripts works.
sys.path.append(cado_python_libs_path)



from workunit import Workunit
import ssl



# File locking functions are specific to Unix/Windows/MacOS platforms.
# The FileLock class is an Interface with static methods.

if os.name == "nt":
    import msvcrt
    class FileLock(object):
        @staticmethod
        def lock(f, exclusive=False, blocking=True):
            """ Lock a file """
            # How do you get shared locks under Windows? There is a Python
            # for Windows extension available at
            # http://sourceforge.net/projects/pywin32/files/?source=navbar
            # which offers many functions from the Win32 API, including file
            # locking, but it is a binary extension for Python, and you have
            # to choose the exact right package for your Python version and
            # platform... can't be imported just like that.
            # The msvcrt.locking function appears to support only one kind
            # of locking, i.e., an exclusive lock. This mostly defeats the
            # purpose of locking here, as we want to prevent clients from
            # reading files that another client is currently downloading.
            if not exclusive:
                # Don't have shared lock - bail out. FIXME
                return
            else:
                # For now, do nothing, until I figure out file locking under
                # Windows
                return
            # Lock one byte from the start of the file
            mode = msvcrt.LK_LOCK if blocking else msvcrt.LK_NBLCK
            pos = f.tell()
            f.seek(0)
            msvcrt.locking(f.fileno(), mode, 1)
            f.seek(pos)
        @staticmethod
        def unlock(f):
            """ Unlock a file """
            # Unlock one byte from the start of the file
            # For now, do nothing, until I figure out file locking under
            # Windows
            return
            pos = f.tell()
            f.seek(0)
            msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, 1)
            f.seek(pos)
# could replace "posix" by "xxx" here if os.name is "posix" but you still get
# the error message "IOError: [Errno 37] No locks available"
# https://lists.gforge.inria.fr/pipermail/cado-nfs-discuss/2016-May/000634.html
# https://lists.gforge.inria.fr/pipermail/cado-nfs-discuss/2016-May/000636.html
elif os.name == "posix":
    import fcntl
    class FileLock(object):
        @staticmethod
        def lock(f, exclusive=False, blocking=True):
            """ Lock a file
            
            If exclusive is True, lock for exclusive (a.k.a "write") access,
            otherwise lock for shared (a.k.a. "read") access.
            If blocking is False, don't block in case of already-locked
            file, but raise IOError with EACCES or EAGAIN (depending on OS).
            """
            mode = fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH
            mode |= 0 if blocking else fcntl.LOCK_NB
            fcntl.flock(f.fileno(), mode)
        @staticmethod
        def unlock(f):
            """ Unlock a file """
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)
else:
    # No file locking. FIXME: What about MacOS?
    class FileLock(object):
        @staticmethod
        def lock(f, exclusive=False, blocking=True):
            """ Do nothing """
            pass
        @staticmethod
        def unlock(f):
            """ Do nothing """
            pass

# In Python 3.0, 3.1, 3.2.x < 3.2.4, 3.3.x < 3.3.1, use a fixed BytesGenerator
# which accepts a bytes input. The fact that the BytesGenerator in these Python
# versions doesn't is a bug, see http://bugs.python.org/issue16564
# Update: the first bugfix committed in that bugtracker and shipped in Python
# versions 3.2.4, 3.2.5, 3.3.2, ... is still buggy, see
# http://bugs.python.org/issue19003
# and we have to use a different work-around...

# These Python version have bug type #1
BUGGY_MIMEENCODER1 = (
    (3,0,0), (3,0,1),
    (3,1,0), (3,1,1), (3,1,2), (3,1,3), (3,1,4), (3,1,5),
    (3,2,0), (3,2,1), (3,2,2), (3,2,3),
    (3,3,0)
)
BUGGY_MIMEENCODER2 = object()

def test_buggy_mime_encoder():
    """ Return which kind of bug the MIME encoder in this Python version has

    If it is bug type 1 (the original bug), return BUGGY_MIMEENCODER1.
    If it is bug type 2 (the incorrectly fixed one, which mangles newlines),
    return BUGGY_MIMEENCODER2.
    Otherwise return None.
    """

    # The version 2 encoder is ok, afaik
    if tuple(sys.version_info)[0] == 2:
        return None

    # This bug type 1 was replaced by bug type 2 in all Python branches,
    # so it should occur only in a well-known set of versions
    if tuple(sys.version_info)[0:3] in BUGGY_MIMEENCODER1:
        return BUGGY_MIMEENCODER1

    # The "bug 2" type of encoder inserts a '\n' character after the '\x0b'
    bytesdata = b'a\x0bb'
    msg = MIMEApplication(bytesdata, _encoder=email.encoders.encode_noop)
    s = BytesIO()
    g = email.generator.BytesGenerator(s)
    g.flatten(msg)
    wireform = s.getvalue()
    msg2 = email.message_from_bytes(wireform)
    decoded = msg2.get_payload(decode=True)
    # If the decoded data is neither correct, nor mangled in the way the
    # bug 2 does, then it's some bug we don't know about yet - bail out
    assert decoded == bytesdata or decoded == b'a\x0b\nb'
    if decoded != bytesdata:
        return BUGGY_MIMEENCODER2
    else:
        return None

HAVE_WGET = False
HAVE_CURL = False

ENCODERBUG = test_buggy_mime_encoder()
if ENCODERBUG is BUGGY_MIMEENCODER1:
    class FixedBytesGenerator(email.generator.BytesGenerator):
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
                super(FixedBytesGenerator, self)._handle_text(msg)
            else:
                # Payload is neither bytes nor string - this can't be right
                raise TypeError('bytes payload expected: %s' % type(payload))
        _writeBody = _handle_bytes
elif ENCODERBUG is BUGGY_MIMEENCODER2:
    if tuple(sys.version_info)[0:2] == (3, 2):
        from email.message import _has_surrogates
    else:
        from email.utils import _has_surrogates
    import re
    fcre = re.compile(r'^From ', re.MULTILINE)
    class FixedBytesGenerator(email.generator.BytesGenerator):
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
elif sys.version_info[0] == 2:
    # In Python 2.x, use the regular email generator
    # pylint: disable=C0103
    FixedBytesGenerator = email.generator.Generator
else:
    # In other Python versions, use the (hopefully) bug-fixed bytes generator
    # pylint: disable=E1101
    # pylint: disable=C0103
    FixedBytesGenerator = email.generator.BytesGenerator


# Under Python 2, subclass urllib2.HTTPSHandler and httplib.HTTPSConnection
# and check the certificate and server subject identity when opening a HTTPS
# connection
if sys.version_info[0] == 2:
    import httplib
    class MyHTTPSConnection(httplib.HTTPSConnection):
        """ HTTPS connections with certificate subject identity check """
        ca_file = None
        check_hostname = True

        def connect(self):
            """ Open a connection, then wrap the socket with SSL, verify the
            server certificate, and the the server certificate's subject.
            """
            # Python 2.6 does not have the source_address attribute
            if sys.version_info[0:2] == (2, 7):
                sock = socket.create_connection((self.host, self.port),
                                                self.timeout, self.source_address)
            else:
                sock = socket.create_connection((self.host, self.port),
                                                self.timeout)

            if self._tunnel_host:
                self.sock = sock
                self._tunnel()
            self.sock = ssl.wrap_socket(sock, self.key_file, self.cert_file,
                                        ca_certs=self.ca_file,
                                        cert_reqs=ssl.CERT_REQUIRED)

            cert = self.sock.getpeercert()
            host = self.host.split(":")[0]

            if not self.check_hostname:
                return

            had_DNS_name = False
            certhost = ""
            if "subjectAltName" in cert:
                # print(cert["subjectAltName"])
                for (typeid, value) in cert["subjectAltName"]:
                    # Wildcard "*" not implemented, but we don't use it anyhow
                    assert not "*" in value
                    if typeid == "DNS" and value == host:
                        return
                    if typeid == "DNS":
                        had_DNS_name = True
                    if typeid == "IP Address" and value == host:
                        return
            for field in cert['subject']:
                if field[0][0] == 'commonName':
                    certhost = field[0][1]
            # Check common name only if there was no SAN DNS entry
            if not had_DNS_name and certhost == host:
                return

            raise ssl.SSLError("Host name '%s' doesn't match certificate host '%s'"
                               % (host, certhost))	

    class MyHTTPSHandler(urllib_request.HTTPSHandler):
        """ HTTPS handler that uses MyHTTPSConnection for verifying the
        certificate's subject identity
        """
        def https_open(self, req):
            return self.do_open(MyHTTPSConnection, req)

    myOpenerDirector = urllib_request.build_opener(MyHTTPSHandler)
    urllib_request.install_opener(myOpenerDirector)


def create_daemon(workdir=None, umask=None, keepfd=None):
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
    if (hasattr(os, "devnull")):
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
    if (maxfd == resource.RLIM_INFINITY):
        maxfd = maxfd_default
  
    # Iterate through and close all file descriptors.
    for fd in range(0, maxfd):
        try:
            if keepfd is None or not fd in keepfd:
                os.close(fd)
        except OSError:	# ERROR, fd wasn't open to begin with (ignored)
            pass

    # Redirect the standard I/O file descriptors to the specified file.  Since
    # the daemon has no controlling terminal, most daemons redirect stdin,
    # stdout, and stderr to /dev/null.  This is done to prevent side-effects
    # from reads and writes to the standard I/O file descriptors.

    # This call to open is guaranteed to return the lowest file descriptor,
    # which will be 0 (stdin), since it was closed above.
    os.open(redirect_to, os.O_RDWR)	# standard input (0)

    # Duplicate standard input to standard output and standard error.
    os.dup2(0, 1)			# standard output (1)
    os.dup2(0, 2)			# standard error (2)

    return(0)

class WuMIMEMultipart(MIMEMultipart):
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
        logging.debug ("Adding result file %s to upload", filepath)
        try:
            with open(filepath, "rb") as infile:
                filedata = infile.read()
        except IOError as err:
            logging.error ("Could not read file %s: %s", filepath, str(err))
            return
        self.attach_data(name, filename, filedata, filetype, command)

    def attach_key(self, key, value):
        ''' Attach a simple key=value pair '''
        attachment = MIMEText(str(value))
        attachment.add_header('Content-Disposition', 'form-data', 
                              name=key)
        self.attach(attachment)

    def flatten(self, debug = 0):
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


class SharedFile(object):
    def __init__(filename, mode=0o777):
        # Try to create and open the file exclusively
        self.filename = filename
        flags = os.O_CREAT | os.O_RDWR | os.O_EXCL
        try:
            self.fd = os.open(filename, flags, mode)
        except OSError as err:
            if err.errno == errno.EEXIST: # If the file already existed
                self.existed = True
                self.wait_until_positive_filesize(filename)
                self.file = open(filename, "r+b")
                FileLock.lock(self.file)
                return
            else:
                raise
        self.existed = False
        self.file = os.fdopen(fd, "r+b")
        FileLock.lock(self.file, exclusive=True)

    def close():
        FileLock.unlock(self.file)
        self.file.close() # This should also close the fd

    def delete():
        if self.existed:
            FileLock.unlock(self.file)
            FileLock.lock(self.file, exclusive=True)
        try:
            os.remove(self.filename)
        except OSError as err:
            if err.errno == errno.ENOENT:
                pass
            else:
                raise
        self.close()

    def wait_until_positive_filesize(self, timeout = 60):
        # There is a possible race condition here. If process A creates 
        # the file, then process B tries and finds that the file exists
        # and immediately get a shared lock for reading, then process A 
        # can never get an exclusive lock for writing.
        # To avoid this, we let process B wait until the file has 
        # positive size, which implies that process A must have the 
        # lock already. After 60 seconds, assume the file really has 0 
        # bytes and return
        slept = 0
        while slept < timeout and os.path.getsize(self.filename) == 0:
            logging.warning("Sleeping until %s contains data", self.filename)
            time.sleep(1)
            slept += 1
        if slept == timeout:
            logging.warning("Slept %d seconds, %s still has no data", 
                            timeout, self.filename)
        return

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

def run_command(command, print_error=True, **kwargs):
    """ Run command, wait for it to finish, return exit status, stdout
    and stderr

    If print_error is True and the command exits with a non-zero exit code,
    print stdout and stderr to the log. If a KeyboardInterrupt exception
    occurs while waiting for the command to finish, the command is
    terminated.
    """
    
    command_str = command if isinstance(command, str) else " ".join(command)
    
    if os.name == "nt":
        # We need to call bash explicitly as the WU COMMAND assumes POSIX
        # syntax which the Windows command line shell does not implement.
        # Turn command into an array so that "bash -c" gets the WU command
        # as a single parameter, i.e., we defer propery quoting command_str
        # to to the subprocess module.
        command = ["bash",  "-c", command_str]
        command_str = " ".join(command)
        close_fds = False
    else:
        close_fds = True

    logging.info ("Running %s", command_str)

    child = subprocess.Popen(command,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             close_fds=close_fds,
                             **kwargs)

    # If we receive SIGTERM (the default signal for "kill") while a
    # subprocess is running, we want to be able to terminate the
    # subprocess, too, so that the system is not kepy busy with
    # orphaned processes.
    # Python installs by default a signal handler for SIGINT which
    # raises the KeyboardInterrupt exception. This is convenient, as
    # it lets us simply terminate the child in an exception handler.
    # Thus we install the signal handler of SIGINT for SIGTERM as well,
    # so that SIGTERM likewise raises a KeyboardInterrupt exception.

    sigint_handler = signal.getsignal(signal.SIGINT)
    signal.signal(signal.SIGTERM , sigint_handler)  

    # Wait for command to finish executing, capturing stdout and stderr 
    # in output tuple
    try:
        (stdout, stderr) = child.communicate()
    except KeyboardInterrupt:
        logging.critical("KeyboardInterrupt received, killing child "
                         "process with PID %d", child.pid)
        child.terminate()
        raise # Re-raise KeyboardInterrupt to terminate cado-nfs-client.py
    
    # Un-install our handler and revert to the default handler
    signal.signal(signal.SIGTERM , signal.SIG_DFL)
    
    if print_error and child.returncode != 0:
        logging.error("Command resulted in exit code %d", child.returncode)
        if stdout:
            logging.error("Stdout: %s", stdout.rstrip())
        if stderr:
            logging.error("Stderr: %s", stderr.rstrip())
    return (child.returncode, stdout, stderr)


class WorkunitProcessor(object):
    def __init__(self, workunit, settings):
        self.settings = settings
        self.workunit = workunit
        self.errorcode = 0 # will get set if any command exits with code != 0
        self.failedcommand = None # If any command exits with code != 0, this
                                  # get set to the index of the failed command
        self.stdio = {"stdout": [], "stderr": []}

    def  __str__(self):
        return "Processor for Workunit:\n%s" % super(WorkunitProcessor, self)

    def renice(self):
        os.nice(int(self.settings["NICENESS"]))

    @staticmethod
    def is_executable(filename):
        """ Test that the file exists and, if the stat object knows the
        "executable by user" constant, that it is executable
        """
        return os.path.isfile(filename) and not (hasattr(stat, "S_IXUSR") and
                (os.stat(filename).st_mode & stat.S_IXUSR) == 0)

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

    def run_commands(self):
        if self.result_exists():
            if self.settings["KEEPOLDRESULT"]:
                return True
            else:
                self.cleanup()
        files = {}
        
        # To which directory do workunit files map?
        dirs = {"FILE": self.settings["DLDIR"],
                "RESULT": self.settings["WORKDIR"]}
        
        for key in dirs:
            for (index, filename) in enumerate(self.workunit.get(key, [])):
                if not isinstance(filename, str):
                        filename = filename[0] # Drop checksum value
                # index is 0-based, add 1 to make FILE1, FILE2, etc. 1-based
                files["%s%d" % (key, index + 1)] = \
                        os.path.join(dirs[key], filename)

        key = "EXECFILE"
        for (index, filename) in enumerate(self.workunit.get(key, [])):
            if not isinstance(filename, str):
                    filename = filename[0] # Drop checksum value
            if self.settings["BINDIR"]:
                searchpath=self.settings["BINDIR"].split(';')
                suggest=self.workunit.get("SUGGEST_" + key, None)
                if suggest:
                    searchpath+=[os.path.join(x, suggest) for x in searchpath]
                binfile = self.find_binary(filename, searchpath)
                if binfile is None:
                    raise Exception("Binary file %s not found" % filename)
            else:
                binfile = os.path.join(self.settings["DLDIR"], filename)
            files["%s%d" % (key, index + 1)] = binfile
        
        for (counter, command) in enumerate(self.workunit.get("COMMAND", [])):
            command = Template(command).safe_substitute(files)

            # If niceness command line parameter was set, call self.renice() 
            # in child process, before executing command
            if int(self.settings["NICENESS"]) > 0:
                renice_func = self.renice
            else:
                renice_func = None

            if self.settings["override"]:
                mangled=[]
                orig=re.split(' *', command)
                while orig:
                    a=orig.pop(0)
                    repl=None
                    for sub in self.settings["override"]:
                        if re.match('^-{1,2}' + sub[0] + '$', a):
                            repl=sub[1]
                    mangled.append(a)
                    if repl is not None:
                        oldvalue=orig.pop(0)
                        logging.info("Overriding argument %s %s by %s %s in command line (substitution %s %s)" % (a, oldvalue, a, sub[1], sub[0], sub[1]))
                        mangled.append(sub[1])
                command=' '.join(mangled)

            (returncode, stdout, stderr) = run_command(command, shell=True,
                    preexec_fn=renice_func)

            self.stdio["stdout"].append(stdout)
            self.stdio["stderr"].append(stderr)

            if returncode != 0:
                self.failedcommand = counter
                self.errorcode = returncode
                return False
            else:
                logging.debug ("Command exited successfully")
        return True

    def result_exists(self):
        ''' Check whether all result files already exist, returns True of False 
        '''
        # If there is no RESULT line in the workunit, always run commands
        if self.workunit.get("RESULT", None) is None:
            return False
        for filename in self.workunit.get("RESULT", []):
            filepath = os.path.join(self.settings["WORKDIR"], filename)
            if not os.path.isfile(filepath):
                logging.info ("Result file %s does not exist", filepath)
                return False
            logging.info ("Result file %s already exists", filepath)
        logging.info ("All result files already exist")
        return True

    def cleanup(self):
        ''' Delete uploaded result files and files from DELETE lines '''
        logging.info ("Cleaning up for workunit %s", self.workunit.get_id())
        for filename in self.workunit.get("RESULT", []):
            filepath = os.path.join(self.settings["WORKDIR"], filename)
            logging.info ("Removing result file %s", filepath)
            try:
                os.remove(filepath)
            except OSError as err:
                logging.error("Could not remove file: %s", err)
        for filename in self.workunit.get("DELETE", []):
            filepath = os.path.join(self.settings["WORKDIR"], filename)
            logging.info ("Removing file %s", filepath)
            os.remove(filepath)

class WorkunitParseError(ValueError):
    """ Parsing the workunit failed """
    pass

class WorkunitClientToFinish(Exception):
    """ we received a 410 (probably while attempting to download a WU) """
    def __init__(self, explanation):
        self.text = explanation
    def  __str__(self):
        return self.text

class WorkunitClient(object):
    def __init__(self, settings):
        self.settings = settings
        
        self.wu_filename = os.path.join(self.settings["DLDIR"], 
                                        self.settings["WU_FILENAME"])

        force_reload=False

        while True:
            self.download_wu(force_reload=force_reload)

            # Get an exclusive lock to avoid two clients working on the same 
            # workunit
            try:
                self.wu_file = open_exclusive(self.wu_filename)
            except FileLockedException:
                logging.error("File '%s' is already locked. This may "
                              "indicate that two clients with clientid '%s' are "
                              "running. Terminating.", 
                              self.wu_filename, self.settings["CLIENTID"])
                raise

            logging.debug ("Parsing workunit from file %s", self.wu_filename)
            wu_text = self.wu_file.read()
            # WU file stays open so we keep the lock

            try:
                self.workunit = Workunit(wu_text)
            except Exception as err:
                logging.error("Invalid workunit file: %s", err)
                self.cleanup()
                raise WorkunitParseError()
            if not force_reload and self.workunit.get("DEADLINE") and time.time() > float(self.workunit.get("DEADLINE")):
                logging.warn("Old workunit file %s has passed deadline (%s), ignoring",
                        self.wu_filename, 
                        time.asctime(time.localtime(float(self.workunit.get("DEADLINE")))))
                os.remove(self.wu_filename)
                close_exclusive(self.wu_file)
                force_reload=True
            else:
                break

        logging.debug ("Workunit ID is %s", self.workunit.get_id())
    
    def download_wu(self, *args, **kwargs):
        # Download the WU file if none exists
        url = self.settings["GETWUPATH"]
        options = "clientid=" + self.settings["CLIENTID"]
        self.get_missing_file(url, self.wu_filename, options=options, *args, **kwargs)

    def cleanup(self):
        logging.info ("Removing workunit file %s", self.wu_filename)
        close_exclusive(self.wu_file)
        os.remove(self.wu_filename)

    @staticmethod
    def _urlopen_maybe_https(request, cafile=None, check_hostname=True):
        """ Treat requests for HTTPS differently depending on whether we are
        on Python 2 or Python 3.
        """
        if isinstance(request, urllib_request.Request):
            if sys.version_info[0:2] < (3,3):
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
            if sys.version_info[0] == 3:
                # Python 3 implements HTTPS certificate checks, we can just
                # let urllib do the work for us

                # To skip the hostname check under Python 3, I'll have to
                # override # HTTPSConnection.__init__() and register that
                # with the urlopen director
                if not check_hostname:
                    raise Exception("Error, not checking hostname not "
                                    "implemented for Python 3 yet")
                return urllib_request.urlopen(request, cafile=cafile)
            else:
                # For the time being, we just use HTTPS without check.
                # We should never get here, as we use wget or curl as
                # fall-backs under Python 2, and if neither is available,
                # cado-nfs-client.py aborts in the initialisation phase.

                # Ugly hack: urllib2 does not provide for parameter passing
                # to HTTPSConnection, so we modify the class variable
                # default_cert_file of MyHTTPSConnection. YUCK.
                MyHTTPSConnection.ca_file = cafile
                MyHTTPSConnection.check_hostname = check_hostname
                return urllib_request.urlopen(request)
        else:
            # If we are not using HTTPS, we can just let urllib do it, 
            # and there is no need for a cafile parameter (which Python 2
            # urlopen() does not accept)
            return urllib_request.urlopen(request)
    
    def _urlopen(self, request, wait, is_upload=False, cafile=None):
        """ Wrapper around urllib2.urlopen() that retries in case of
        connection failure.
        """
        conn = None
        waiting_since = 0
        # this knowingly mixes http status codes in the 400- 500- range
        # with errno errors. It's ugly.
        last_error = 0
        current_error = 0
        # record the number of connection failures
        connfailed = 0
        maxconnfailed = int(self.settings["MAX_CONNECTION_FAILURES"])
        silent_wait = self.settings["SILENT_WAIT"]
        check_hostname = not self.settings["NO_CN_CHECK"]
        while conn is None:
            try:
                conn = WorkunitClient._urlopen_maybe_https(request, cafile=cafile,
                        check_hostname=check_hostname)
            except urllib_error.HTTPError as error:
                current_error = error.code
                if error.code == 410:
                    # We interpret error code 410 as the work unit server
                    # being gone for good. This instructs us to terminate
                    # the workunit client, which we do by letting an
                    # exception pop up a few levels up (eeek)
                    raise WorkunitClientToFinish("Received 410 from server")
                conn = None
                errorstr = "URL error: %s" % str(error)
            except urllib_error.URLError as error:
                conn = None
                errorstr = "URL error: %s" % str(error)
                current_error = error.errno
            except BadStatusLine as error:
                conn = None
                errorstr = "Bad Status line: %s" % str(error)
            except socket.error as error:
                conn = None
                errorstr = "Connection error: %s" % str(error)
            if not conn:
                givemsg = is_upload or not silent_wait or waiting_since == 0
                if current_error > 0:
                    if current_error != last_error:
                        givemsg = True
                    if current_error == errno.ECONNREFUSED or \
                            current_error == errno.ECONNRESET:
                        connfailed += 1
                    else:
                        connfailed = 0

                if givemsg:
                    logging.error("%s failed, %s", 
                                  'Upload' if is_upload else 'Download', 
                                  errorstr)
                    if waiting_since > 0:
                        logging.error("Waiting %s seconds before retrying (I have been waiting since %s seconds)", wait, waiting_since)
                    else:
                        logging.error("Waiting %s seconds before retrying", wait)
                if current_error > 0:
                    last_error = current_error
                    if connfailed > maxconnfailed:
                        raise WorkunitClientToFinish("Connection failed %s times" % maxconnfailed)
                time.sleep(wait)
                waiting_since+=wait
        if waiting_since > 0:
            logging.info ("Opened URL %s after %s seconds wait", request, waiting_since)
        return conn

    @staticmethod
    def get_content_charset(conn):
        """ Returns the character set of the server's response. 
        
        Defaults to latin-1 if no charset header was sent.
        The encoding may matter if the path names for the uploaded files 
        contain special characters, like accents.
        """
        if sys.version_info[0] == 3:
            charset = conn.info().get_content_charset()
            if charset is None:
                charset = "latin-1"
        else:
            charset = "latin-1" # Default value
            for item in conn.info().getplist():
                pair = item.split("=")
                if len(pair) == 2 and pair[0].strip() == "charset":
                    charset = pair[1].strip()
        return charset

    def external_get_file(self, command, url, wait):
        """ Runs a command to download a file, retrying indefinitely in case
        of error
        """
        silent_wait=self.settings["SILENT_WAIT"]
        waiting_since = 0
        while True:
            (rc, stdout, stderr) = run_command (command)
            if rc == 0:
                return True
            if waiting_since == 0 or not silent_wait:
                logging.error("Download of %s failed. Waiting %s seconds before "
                              "retrying,", url, wait)
            waiting_since+=wait
            time.sleep(wait)

    def wget_file(self, url, wait, dlpath, cafile=None):
        """ Download via wget
        
        This is used as a fall-back for doing HTTPS downloads when running
        under Python 2, whose ssl module does not implement SSL certificate
        checks.
        """
        command = ["wget", "-O", dlpath]
        if cafile:
            command.append("--ca-certificate=%s" % cafile)
        command.append(url)
        return self.external_get_file(command, url, wait)
        
    def curl_get_file(self, url, wait, dlpath, cafile=None):
        """ Download via curl
        
        Like wget_file(), this is used as a fall-back.
        """
        command = ["curl", "--silent", "--show-error", "--fail", "--output", dlpath]
        if cafile:
            command += ["--cacert", cafile]
        command.append(url)
        return self.external_get_file(command, url, wait)
        
    def get_file(self, urlpath, dlpath=None, options=None):
        # print('get_file("' + urlpath + '", "' + dlpath + '")')
        if dlpath == None:
            filename = urlpath.split("/")[-1]
            dlpath = os.path.join(self.settings["DLDIR"], filename)
        url = self.settings["SERVER"].rstrip("/") + "/" + urlpath.lstrip("/")
        if options:
            url = url + "?" + options
        cafile = self.settings.get("CERTFILE", None)
        logging.info ("Downloading %s to %s (cafile = %s)", url, dlpath, cafile)
        wait = float(self.settings["DOWNLOADRETRY"])
        # If we want HTTPS and are running under Python 2, we use wget to do
        # the actual download, as the Python 2 urllib does not implement
        # actually checking the certificate
        # This is a rather ugly hack. It would be nicer to copy the required
        # parts from a fully functional SSL library. TODO.
        if self.settings["USE_EXTERNAL_DL"]:
            if HAVE_WGET:
                return self.wget_file(url, wait, dlpath, cafile=cafile)
            elif HAVE_CURL:
                return self.curl_get_file(url, wait, dlpath, cafile=cafile)

        request = self._urlopen(url, wait, cafile=cafile)
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
                # bytes and return
                logging.warning("Looks like another process already created "
                                "file %s", dlpath)
                self.wait_until_positive_filesize(dlpath)
                return
            else:
                raise
        outfile = os.fdopen(fd, "wb")
        FileLock.lock(outfile, exclusive=True)
        shutil.copyfileobj (request, outfile)
        FileLock.unlock(outfile)
        outfile.close() # This should also close the fd
        request.close()
    
    def get_missing_file(self, urlpath, filename, checksum=None,
                         options=None, force_reload=False):
        """ Downloads a file if it does not exit already.

        Also checks the checksum, if specified; if the file already exists and
        has a wrong checksum, it is deleted an downloaded anew. If the
        downloaded file has the wrong checksum, it is deleted and downloaded
        anew. If the downloaded file has the same, incorrect checksum twice
        in a row, the function returns False. In all other cases, it returns
        True.
        """
        # print('get_missing_file(%s, %s, %s)' % (urlpath, filename, checksum))
        if os.path.isfile(filename):
            if force_reload:
                logging.info ("%s already exists, removing because of force_reload", filename)
                os.remove(filename)
            else:
                logging.info ("%s already exists, not downloading", filename)
                if checksum is None:
                    return True
                filesum = self.do_checksum(filename)
                if filesum.lower() == checksum.lower():
                    return True
                logging.error ("Existing file %s has wrong checksum %s, "
                               "workunit specified %s. Deleting file.", 
                               filename, filesum, checksum)
                os.remove(filename)

        # If checksum is wrong and does not change during two downloads, exit 
        # with failue, as apparently the file on the server and checksum in 
        # workunit do not agree
        last_filesum = None
        while True:
            # we were catching HTTPError here previously. Useless now ?
            self.get_file(urlpath, filename, options)
            if checksum is None:
                return True
            filesum = self.do_checksum(filename)
            if filesum.lower() == checksum.lower():
                return True
            if not last_filesum is None and filesum == last_filesum:
                logging.error ("Downloaded file %s has same wrong checksum %s" 
                               " again. Exiting.", filename, filesum)
                return False
            logging.error ("Downloaded file %s has wrong checksum %s, " 
                           "workunit specified %s. Deleting file.", 
                           filename, filesum, checksum)
            os.remove(filename)
            last_filesum = filesum
        return True

    def get_files(self):
        files_to_download = self.workunit.get("FILE", [])
        if not self.settings["BINDIR"]:
            files_to_download += self.workunit.get("EXECFILE", [])
        for (filename, checksum) in files_to_download:
            templ = Template(filename)
            archname = templ.safe_substitute({"ARCH": self.settings["ARCH"]})
            dlname = templ.safe_substitute({"ARCH": ""})
            dlpath = os.path.join(self.settings["DLDIR"], dlname)
            if self.settings["NOSHA1CHECK"]:
                checksum = None
            if not self.get_missing_file (archname, dlpath, checksum):
                return False
            # Try to lock the file once to be sure that download has finished
            # if another cado-nfs-client is doing the downloading
            with open(dlpath) as file_to_lock:
                FileLock.lock(file_to_lock)
                FileLock.unlock(file_to_lock)
            
            if os.name != "nt" and \
                    filename in dict(self.workunit.get("EXECFILE", [])):
                mode = os.stat(dlpath).st_mode
                if mode & stat.S_IXUSR == 0:
                    logging.info ("Setting executable flag for %s", dlpath)
                    os.chmod(dlpath, mode | stat.S_IXUSR)
        return True

    def attach_result(self, processor, mimedata):
        # Build a multi-part MIME document containing the WU id and result file
        mimedata.attach_key("WUid", self.workunit.get_id())
        mimedata.attach_key("clientid", self.settings["CLIENTID"])
        if processor.errorcode:
            mimedata.attach_key("errorcode", processor.errorcode)
        if not processor.failedcommand is None:
            mimedata.attach_key("failedcommand", processor.failedcommand)
        for filename in self.workunit.get("RESULT", []):
            filepath = os.path.join(self.settings["WORKDIR"], filename)
            logging.info("Attaching file %s to upload", filepath)
            mimedata.attach_file("results", filename, filepath, "RESULT")
        for name in processor.stdio:
            for (counter, data) in enumerate(processor.stdio[name]):
                if data:
                    logging.info ("Attaching %s for command %s to upload", 
                                  name, counter)
                    filename = "%s.%s%d" % (self.workunit.get_id(), name,
                                            counter)
                    mimedata.attach_data("results", filename, data, name,
                                         counter)
        return mimedata

    def upload_result(self, processor):
        # Make POST data
        mimedata = WuMIMEMultipart()
        self.attach_result(processor, mimedata)
        postdata = mimedata.flatten(debug=int(self.settings["DEBUG"]))
        # logging.debug("POST data: %s", mimedata)

        url = self.settings["SERVER"].rstrip("/") + "/" + \
                self.settings["POSTRESULTPATH"].lstrip("/")
        logging.info("Sending result for workunit %s to %s", 
                     self.workunit.get_id(), url)
        request = urllib_request.Request(url, data=postdata, 
                                         headers=dict(mimedata.items()))
        wait = float(self.settings["DOWNLOADRETRY"])
        cafile = self.settings.get("CERTFILE", None)
        conn = self._urlopen(request, wait, is_upload=True, cafile=cafile)
        if not conn:
            return False
        response = conn.read()

        encoding = self.get_content_charset(conn)
        if sys.version_info[0] == 2:
            response_str = unicode(response, encoding=encoding)
        else:
            response_str = response.decode(encoding=encoding)
        logging.debug ("Server response:\n%s", response_str)
        conn.close()
        return True

    @staticmethod
    def wait_until_positive_filesize(filename, timeout = 60):
        slept = 0
        while slept < timeout and os.path.getsize(filename) == 0:
            logging.warning("Sleeping until %s contains data", filename)
            time.sleep(1)
            slept += 1
        if slept == timeout:
            logging.warning("Slept %d seconds, %s still has no data", 
                            timeout, filename)
        return

    @staticmethod
    def do_checksum(filename, checksum = None):
        """ Computes the SHA1 checksum for a file. If checksum is None, returns 
            the computed checksum. If checksum is not None, return whether the
            computed SHA1 sum and checksum agree """
        blocksize = 65536
        sha1hash = hashlib.sha1() # pylint: disable=E1101
        # Like when downloading, we wait until the file has positive size, to
        # avoid getting the shared lock right after the other process created
        # the file but before it gets the exclusive lock
        WorkunitClient.wait_until_positive_filesize(filename)
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
        else:
            return filesum.lower() == checksum.lower()

    def have_terminate_request(self):
        return not self.workunit.get("TERMINATE", None) is None

    def process(self):
        # If all output files exist, send them, return WU as done
        # Otherwise, run commands in WU. If no error and all output 
        #   files exist, send them, return WU as done
        # print(wu)
        if self.have_terminate_request():
            logging.info ("Received TERMINATE, exiting")
            return False

        if self.get_files():
            processor = WorkunitProcessor(self.workunit, self.settings)
            processor.run_commands()
            upload_ok = self.upload_result(processor)
            processor.cleanup()
        else:
            # TODO: notify server of error so it can re-issue immediately?
            logging.error("Could not download a required file, discarding "
                          "workunit %s", self.workunit.get_id())
            upload_ok = True # Client should continue
        self.cleanup()
        return upload_ok


def get_ssl_certificate(server, port=443, retry=False, retrytime=0):
    """ Download the SSL certificate from the server.
    
    In case of connection refused error, if retry is True, retry
    indefinitely waiting retrytime seconds between tries, and if
    retry is False, return None.
    """
    while True:
        try:
            cert = ssl.get_server_certificate((server, int(port)),
                                              ssl_version=ssl.PROTOCOL_SSLv23,
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


def get_missing_certificate(certfilename, netloc, fingerprint, retry=False,
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
        cert = get_ssl_certificate(*address_port, retry=retry,
                retrytime=retrytime)
        if cert is None:
            return False
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

def try_wget():
    try:
        return run_command(["wget", "-V"])[0] == 0
    except OSError:
        return False

def try_curl():
    try:
        (rc, stdout, stderr) = run_command(["curl", "-V"])
    except OSError:
        return False
    if rc != 0:
        return False
    match = False
    version_lines = stdout.splitlines()
    if version_lines:
        match = re.match("curl (?:\d+.\d+)", version_lines[0])
    if not match:
        logging.error("curl did not print its version with -V")
        logging.error("Stdout: %s", stdout)
        return False

    if re.search("SecureTransport", version_lines[0]):
        logging.error("Found Apple version of curl with SecureTransport SSL "
                      "backend") 
        logging.error("SecureTransport doesn't allow self-signed certificates "
                      "provided in files. Please see the SSL section in "
                      "README.Python for details.")
        return False
    return True

# Settings which we require on the command line (no defaults)
REQUIRED_SETTINGS = {"SERVER" : (None, "Base URL for WU server")}

# Optional settings with defaults, overrideable on command line, 
# and a help text
OPTIONAL_SETTINGS = {"WU_FILENAME" : 
                     (None, "Filename under which to store WU files"), 
                     "CLIENTID" : (None, "Unique ID for this client. If not "
                                   "specified, a default of "
                                   "<hostname>.<random hex number> is used"), 
                     "DLDIR" : ('download/', "Directory for downloading files"),
                     "WORKDIR" : (None, "Directory for result files"),
                     "BINDIR" : (None, "Directory with existing executable "
                                       "files to use"),
                     "BASEPATH" : (None, "Base directory for download and work "
                                         "directories"),
                     "GETWUPATH" : 
                     ("/cgi-bin/getwu", 
                      "Path segment of URL for requesting WUs from server"), 
                     "POSTRESULTPATH" : 
                     ("/cgi-bin/upload.py", 
                      "Path segment of URL for reporting results to server"), 
                     "DEBUG" : ("0", "Debugging verbosity"),
                     "ARCH" : ("", "Architecture string for this client"),
                     "DOWNLOADRETRY" : 
                     ("10", "Time to wait before download retries"),
                     "CERTSHA1" : (None, "SHA1 of server SSL certificate"),
                     "SILENT_WAIT": (None, "Discard repeated messages about client waiting for work (does not affect uploads)"),
                     "MAX_CONNECTION_FAILURES" : ("999999", "Maximum number of successive connection failures to tolerate"),
                     "NICENESS" : 
                     ("0", "Run subprocesses under this niceness"),
                     "LOGLEVEL" : ("INFO", "Verbosity of logging"),
                     "LOGFILE" : (None, "File to which to write log output. "
                                  "In demon mode, if no file is specified, a "
                                  "default of <workdir>/<clientid>.log is used")
                     }
# Merge the two, removing help string
SETTINGS = dict([(a, b) for (a, (b, c)) in list(REQUIRED_SETTINGS.items()) + \
                                        list(OPTIONAL_SETTINGS.items())])

BAD_WU_MAX = 3 # Maximum allowed number of bad WUs

if __name__ == '__main__':

    def parse_cmdline():
        # Create command line parser from the keys in SETTINGS
        parser = optparse.OptionParser()
        for (arg, default) in REQUIRED_SETTINGS.items():
            parser.add_option('--' + arg.lower(), help=default[1])
        for (arg, default) in OPTIONAL_SETTINGS.items():
            if not default[0] is None:
                parser.add_option('--' + arg.lower(), default=default[0], 
                    help=default[1] + " (default: " + default[0] + ")")
            else:
                parser.add_option('--' + arg.lower(), help=default[1])
        parser.add_option("-d", "--daemon", action="store_true", dest="daemon",
                          help="Daemonize the client")
        parser.add_option("--keepoldresult", default=False, action="store_true", 
                          help="Keep and upload old results when client starts")
        parser.add_option("--nosha1check", default=False, action="store_true", 
                          help="Skip checking the SHA1 for input files")
        parser.add_option("--single", default=False, action="store_true", 
                          help="process only a single WU, then exit")
        parser.add_option("--nocncheck", default=False, action="store_true", 
                          help="Don't check common name/SAN of certificate. "
                          "Currently works only under Python 2.")
        parser.add_option("--externdl", default=False, action="store_true", 
                          help="Use wget or curl for HTTPS downloads")
        parser.add_option("--override", nargs=2, action='append',
                          metavar=('REGEXP', 'VALUE'),
                          help="Modify command-line arguments which match ^-{1,2}REGEXP$ to take the given VALUE. Note that REGEXP cannot start with a dash")
        # Parse command line
        (options, args) = parser.parse_args()

        if args:
            sys.stderr.write("Did not understand command line arguments %s\n" %
                             " ".join(args))
            raise Exception()
        # Copy values to SETTINGS
        for arg in SETTINGS.keys():
            if hasattr(options, arg.lower()):
                SETTINGS[arg] = getattr(options, arg.lower())
        for arg in REQUIRED_SETTINGS.keys():
            if SETTINGS[arg] is None:
                raise Exception("Command line parameter --%s is required" 
                                % arg.lower())
        return options

    def makedirs(path, mode=None, exist_ok=False):
        # Python 3.2 os.makedirs() has exist_ok, but older Python do not
        if sys.version_info[0:2] >= (3,2):
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
    SETTINGS["USE_EXTERNAL_DL"] = options.externdl
    SETTINGS["NO_CN_CHECK"] = options.nocncheck

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
    logfile = None if logfilename is None else open(logfilename, "a")
    logging.basicConfig(level=loglevel)
    if logfile:
        logging.getLogger().addHandler(logging.StreamHandler(logfile))
    logging.info("Starting client %s", SETTINGS["CLIENTID"])
    logging.info("Python version is %d.%d.%d", *sys.version_info[0:3])

    if ENCODERBUG is BUGGY_MIMEENCODER1:
        logging.info("Using work-around #1 for buggy BytesGenerator")
    elif ENCODERBUG is BUGGY_MIMEENCODER2:
        logging.info("Using work-around #2 for buggy BytesGenerator")

    (scheme, netloc) = urlparse(SETTINGS["SERVER"])[0:2]
    still_need_cert = False # This will be set to True if we need the certi-
                            # ficate, but could not download it right away
    if not SETTINGS["CERTSHA1"] is None and scheme != "https":
        logging.warn("Option --certsha1 makes sense only with an https URL, ignoring it.")
    elif SETTINGS["CERTSHA1"] is None and scheme == "https":
        logging.warn("An https URL was given but no --certsha1 option, NO SSL VALIDATION WILL BE PERFORMED.")
    elif not SETTINGS["CERTSHA1"] is None and scheme == "https":
        certfilename = os.path.join(SETTINGS["DLDIR"], "server.%s.pem" % SETTINGS["CERTSHA1"][0:8])
        SETTINGS["CERTFILE"] = certfilename
        # Try downloading the certificate once. If connection is refused,
        # proceed to daemonizing - hopefully server will come up later
        if not get_missing_certificate(certfilename, netloc, SETTINGS["CERTSHA1"]):
            still_need_cert = True
            logging.info("Could not download SSL certificate: The connection was refused.")
            logging.info("Assuming the server will come up later. Will keep trying%s.",
                         " after daemonizing" if options.daemon else "")
        
        # Can we download with HTTPS at all?
        if sys.version_info[0] == 2:
            HAVE_WGET = try_wget()
            if not HAVE_WGET:
                HAVE_CURL = try_curl()
            if not HAVE_WGET and not HAVE_CURL:
                logging.critical("HTTPS requested, but not implemented in "
                        "Python 2 and can't find working wget or curl as "
                        "fall-back. Aborting.")
                sys.exit(1)

    SETTINGS["override"] = options.override

    if options.daemon:
        create_daemon(keepfd=None if logfile is None else [logfile.fileno()])

    if still_need_cert:
        get_missing_certificate(certfilename, netloc, SETTINGS["CERTSHA1"],
                retry=True, retrytime=SETTINGS["DOWNLOADRETRY"])

    client_ok = True
    bad_wu_counter = 0
    while client_ok:
#        try:
            try:
                client = WorkunitClient(settings = SETTINGS)
            except WorkunitParseError:
                bad_wu_counter += 1
                if bad_wu_counter > BAD_WU_MAX:
                    logging.critical("Had %d bad workunit files. Aborting.", bad_wu_counter)
                    break
                continue
            except WorkunitClientToFinish as e:
                logging.info("Client finishing: %s. Bye." % e)
                break
            client_ok = client.process()
            if options.single:
                logging.info("Client processed its WU. Finishing now as implied by --single")
                sys.exit(0)
#        except Exception:
#            logging.exception("Exception occurred")
#            break
