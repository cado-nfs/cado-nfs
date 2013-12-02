#!/usr/bin/env python

import sys
import os
import fcntl
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
from workunit import Workunit
import ssl

# In Python 3.0, 3.1, 3.2.x < 3.2.4, 3.3.x < 3.3.1, use a fixed BytesGenerator
# which accepts a bytes input. The fact that the BytesGenerator in these Python
# versions doesn't is a bug, see http://bugs.python.org/issue16564
# Update: the first bugfix committed in that bugtracker and shipped in Python
# versions 3.2.4, 3.2.5, 3.3.2 is still buggy, and we have to use a different
# work-around...

# These Python version have bug type #1
BUGGY_MIMEENCODER1 = (
    (3,0,0), (3,0,1),
    (3,1,0), (3,1,1), (3,1,2), (3,1,3), (3,1,4), (3,1,5),
    (3,2,0), (3,2,1), (3,2,2), (3,2,3),
    (3,3,0)
)
# These Python version have bug type #2
BUGGY_MIMEENCODER2 = (
    (3,2,4), (3,2,5),
    (3,3,2), (3,3,3)
)

HAVE_WGET = False
HAVE_CURL = False

if tuple(sys.version_info)[0:3] in BUGGY_MIMEENCODER1:
    class FixedBytesGenerator(email.generator.BytesGenerator):
        bug_type = BUGGY_MIMEENCODER1
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
elif tuple(sys.version_info)[0:3] in BUGGY_MIMEENCODER2:
    import re
    from email.utils import _has_surrogates

    fcre = re.compile(r'^From ', re.MULTILINE)
    class FixedBytesGenerator(email.generator.BytesGenerator):
        bug_type = BUGGY_MIMEENCODER2
        # pylint: disable=W0232
        # pylint: disable=E1101
        # pylint: disable=E1102
        # pylint: disable=E1002
        def _handle_application(self, msg):
            # If the string has surrogates the original source was bytes,
            # so just write it back out.
            if msg._payload is None:
                return
            if _has_surrogates(msg._payload) and \
                    not self.policy.cte_type == '7bit':
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


def create_daemon(workdir=None, umask=None, keepfd=None):
    """Disk And Execution MONitor (Daemon)

    Configurable daemon behaviors:

       1.) The current working directory set to the "/" directory.
       2.) The current file creation mode mask set to 0.
       3.) Close all open files (1024). 
       4.) Redirect standard I/O streams to "/dev/null".

    A failed call to fork() now raises an exception.

    References:
       1) Advanced Programming in the Unix Environment: W. Richard Stevens
       2) Unix Programming Frequently Asked Questions:
             http://www.erlenstar.demon.co.uk/unix/faq_toc.html
    Detach a process from the controlling terminal and run it in the
    background as a daemon.

    Published at
    http://code.activestate.com/recipes/278731-creating-a-daemon-the-python-way/

    Changes: workdir is now a parameter, daemon changes CWD only if workdir 
    parameter is specified. umask is also a parameter, and the process' umask
    is set only if a value is specified.
    """

    __author__ = "Chad J. Schroeder"
    __copyright__ = "Copyright (C) 2005 Chad J. Schroeder"

    __revision__ = "$Id$"
    __version__ = "0.2"
    
    # Use a one-element array to fool Python 2 into not binding a local
    # name in handler(). Python3 has 'nonlocal' for that
    sigusr1_received = [False]
    def handler(signum, frame):
        if signum == signal.SIGUSR1:
            sigusr1_received[0] = True

    # The pid of the original process, used for sending a SIGUSR1 later
    original_pid = os.getpid()

    # We need to install the signal handler before forking, to guarantee that
    # the original process has the signal handler installed at the time its
    # grand-child sends the signal
    old_handler = signal.signal(signal.SIGUSR1, handler)

    # Default daemon parameters.

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

    if (pid == 0):	# The first child.
        # Un-install the signal handler in the child process
        signal.signal(signal.SIGUSR1, signal.SIG_DFL)

        # To become the session leader of this new session and the process group
        # leader of the new process group, we call os.setsid().  The process is
        # also guaranteed not to have a controlling terminal.
        os.setsid()

        # Is ignoring SIGHUP necessary?
        #
        # It's often suggested that the SIGHUP signal should be ignored before
        # the second fork to avoid premature termination of the process.  The
        # reason is that when the first child terminates, all processes, e.g.
        # the second child, in the orphaned group will be sent a SIGHUP.
        #
        # "However, as part of the session management system, there are exactly
        # two cases where SIGHUP is sent on the death of a process:
        #
        #   1) When the process that dies is the session leader of a session
        #      that is attached to a terminal device, SIGHUP is sent to all
        #      processes in the foreground process group of that terminal
        #      device.
        #   2) When the death of a process causes a process group to become
        #      orphaned, and one or more processes in the orphaned group are
        #      stopped, then SIGHUP and SIGCONT are sent to all members of the
        #      orphaned group." [2]
        #
        # The first case can be ignored since the child is guaranteed not to
        # have a controlling terminal.  The second case isn't so easy to
        # dismiss.
        # The process group is orphaned when the first child terminates and
        # POSIX.1 requires that every STOPPED process in an orphaned process
        # group be sent a SIGHUP signal followed by a SIGCONT signal.  Since the
        # second child is not STOPPED though, we can safely forego ignoring the
        # SIGHUP signal.  In any case, there are no ill-effects if it is
        # ignored.
        #
        # signal.signal(signal.SIGHUP, signal.SIG_IGN)

        try:
            # Fork a second child and exit immediately to prevent zombies.  This
            # causes the second child process to be orphaned, making the init
            # process responsible for its cleanup.  And, since the first child
            # is a session leader without a controlling terminal, it's possible
            # for it to acquire one by opening a terminal in the future (System
            #  V-based systems).  This second fork guarantees that the child is
            # no longer a session leader, preventing the daemon from ever
            # acquiring a controlling terminal.
            pid = os.fork()	# Fork a second child.
        except OSError as e:
            raise Exception("%s [%d]" % (e.strerror, e.errno))

        if (pid == 0):	# The second child.
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
        else:
            # exit() or _exit()?  See below.
            os._exit(0)	# Exit parent (the first child) of the second child.
    else:
        # Wait for the child to send a SIGUSR1 signal. This gives the grand-
        # child time to print its PID to stdout, before the original process
        # exits and possibly closes the SSH connection which started the
        # client - we want to read the PID over the SSH connection, after all.
        while not sigusr1_received[0]:
            signal.pause()

        # exit() or _exit()?
        # _exit is like exit(), but it doesn't call any functions registered
        # with atexit (and on_exit) or any registered signal handlers.  It also
        # closes any open file descriptors.  Using exit() may cause all stdio
        # streams to be flushed twice and any temporary files may be 
        # unexpectedly removed.  It's therefore recommended that child branches
        # of a fork() and the parent branch(es) of a daemon use _exit().
        os._exit(0)	# Exit parent of the first child.

    # Print the daemon's PID. We inherited the file descriptors from the
    # original process so this should go over the SSH connection, if we
    # were started via SSH.
    sys.stdout.write("PID: %d\n" % os.getpid())
    sys.stdout.flush()
    
    # Tell the original process that it's ok to terminate now
    os.kill(original_pid, signal.SIGUSR1)

    # Close all open file descriptors.  This prevents the child from keeping
    # open any file descriptors inherited from the parent.  There is a variety
    # of methods to accomplish this task.  Three are listed below.
    #
    # Try the system configuration variable, SC_OPEN_MAX, to obtain the maximum
    # number of open file descriptors to close.  If it doesn't exists, use
    # the default value (configurable).
    #
    # try:
    #    maxfd = os.sysconf("SC_OPEN_MAX")
    # except (AttributeError, ValueError):
    #    maxfd = maxfd_default
    #
    # OR
    #
    # if (os.sysconf_names.has_key("SC_OPEN_MAX")):
    #    maxfd = os.sysconf("SC_OPEN_MAX")
    # else:
    #    maxfd = maxfd_default
    #
    # OR
    #
    # Use the getrlimit method to retrieve the maximum file descriptor number
    # that can be opened by this process.  If there is not limit on the
    # resource, use the default value.
    #
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
        with open(filepath, "rb") as infile:
            filedata = infile.read()
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


class FileLockedException(IOError):
    """ Locking a file for exclusive access failed """
    pass

def open_exclusive(filename):
    """ Open a file and get an exclusive lock on it """
    fileobj = open(filename, "r+")
    try:
        fcntl.flock(fileobj, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except IOError as err:
        if err.errno == errno.EACCES or err.errno == errno.EAGAIN:
            fileobj.close()
            raise FileLockedException(errno.EACCES, "File locked", filename)
        raise
    return fileobj

def close_exclusive(fileobj):
    """ Close a file, releasing any held lock on it """
    fcntl.flock(fileobj, fcntl.LOCK_UN)
    fileobj.close()

def run_command(command):
    """ Run command, wait for it to finish, return exit status, stdout 
    and stderr
    """
    child = subprocess.Popen(command,
                             stdout=subprocess.PIPE, 
                             stderr=subprocess.PIPE, 
                             close_fds=True)
    (stdout, stderr) = child.communicate()
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

    def run_commands(self):
        if self.result_exists():
            if self.settings["KEEPOLDRESULT"]:
                return True
            else:
                self.cleanup()
        for (counter, command) in enumerate(self.workunit.get("COMMAND", [])):
            command = Template(command).safe_substitute(self.settings)
            logging.info ("Running command for workunit %s: %s", 
                          self.workunit.get_id(), command)

            # If niceness command line parameter was set, call self.renice() 
            # in child process, before executing command
            if int(self.settings["NICENESS"]) > 0:
                renice_func = self.renice
            else:
                renice_func = None

            # Run the command
            child = subprocess.Popen(command, shell=True, 
                                     stdout = subprocess.PIPE, 
                                     stderr = subprocess.PIPE, 
                                     close_fds = True,
                                     preexec_fn = renice_func)

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
                (child_stdout, child_stderr) = child.communicate()
            except KeyboardInterrupt:
                logging.critical("KeyboardInterrupt received, killing child "
                                 "process with PID %d", child.pid)
                child.terminate()
                raise # Re-raise KeyboardInterrupt to terminate wuclient.py
            
            # Un-install our handler and revert to the default handler
            signal.signal(signal.SIGTERM , signal.SIG_DFL)
            
            self.stdio["stdout"].append(child_stdout)
            self.stdio["stderr"].append(child_stderr)

            if child.returncode != 0:
                logging.error ("Command exited with exit code %s", 
                               child.returncode) 
                if self.stdio["stdout"][-1]:
                    logging.error ("Stdout: %s", self.stdio["stdout"][-1]) 
                if self.stdio["stderr"][-1]:
                    logging.error ("Stderr: %s", self.stdio["stderr"][-1]) 
                self.failedcommand = counter
                self.errorcode = child.returncode
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
            os.remove(filepath)
        for filename in self.workunit.get("DELETE", []):
            filepath = os.path.join(self.settings["WORKDIR"], filename)
            logging.info ("Removing file %s", filepath)
            os.remove(filepath)

class WorkunitClient(object):
    def __init__(self, settings):
        self.settings = settings
        
        self.wu_filename = os.path.join(self.settings["DLDIR"], 
                                        self.settings["WU_FILENAME"])
        self.download_wu()

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

        self.workunit = Workunit(wu_text)
        logging.debug ("Workunit ID is %s", self.workunit.get_id())
    
    def download_wu(self):
        # Download the WU file if none exists
        url = self.settings["GETWUPATH"]
        options = "clientid=" + self.settings["CLIENTID"]
        while not self.get_missing_file(url, self.wu_filename, options=options):
            logging.error("Error downloading workunit file")
            wait = float(self.settings["DOWNLOADRETRY"])
            logging.error("Waiting %s seconds before retrying", wait)
            time.sleep(wait)

    def cleanup(self):
        logging.info ("Removing workunit file %s", self.wu_filename)
        close_exclusive(self.wu_file)
        os.remove(self.wu_filename)

    @staticmethod
    def _urlopen_maybe_https(request, cafile=None):
        """ Treat requests for HTTPS differently depending on whether we are
        on Python 2 or Python 3.
        """
        if isinstance(request, urllib_request.Request):
            if sys.version_info[0:2] < (3,3):
                # In Python 2, get_type() must be used to get the scheme
                scheme = request.get_type().lower()
            else:
                # The .get_Type() method was deprecated in 3.3 and removed
                # in 3.4, now the scheme is stored in the .type attribute
                scheme = request.type.lower()
        else:
            # Assume it's a URL string
            scheme = request.split(":")[0].lower()
        if scheme == "https":
            if sys.version_info[0] == 3:
                # Python 3 implements HTTPS certificate checks, we can just
                # let urllib do the work for us
                return urllib_request.urlopen(request, cafile=cafile)
            else:
                # Python 2 urllib does not implement certificate checks.
                # For the time being, we just use HTTPS without check.
                # We should never get here, as we use wget or curl as
                # fall-backs under Python 2, and if neither is available,
                # wuclient2.py aborts in the initialisation phase.
                return urllib_request.urlopen(request)
        else:
            # If we are not using HTTPS, we can just let urllib do it, 
            # and there is no need for a cafile parameter (which Python 2
            # urlopen() does not accept)
            return urllib_request.urlopen(request)
    
    @staticmethod
    def _urlopen(request, wait, is_upload=False, cafile=None):
        """ Wrapper around urllib2.urlopen() that retries in case of
        connection failure.
        """
        conn = None
        while conn is None:
            try:
                conn = WorkunitClient._urlopen_maybe_https(request, cafile=cafile)
            except urllib_error.URLError as error:
                conn = None
                errorstr = "URL error: %s" % str(error)
            except BadStatusLine as error:
                conn = None
                errorstr = "Bad Status line: %s" % str(error)
            except socket.error as error:
                if error.errno == errno.ECONNRESET:
                    conn = None
                    errorstr = "Error connecting: %s" % str(error)
                else:
                    raise
            if not conn:
                logging.error("%s failed, %s", 
                              'Upload' if is_upload else 'Download', 
                              errorstr)
                logging.error("Waiting %s seconds before retrying", wait)
                time.sleep(wait)
        return conn

    @staticmethod
    def get_content_charset(conn):
        """ Returns the character set of the server's response. 
        
        Defaults to latin-1 if no charset header was sent.
        The encoding may matter if the path names for the uploaded files 
        contain special characters, like accents.
        """
        if sys.version_info[0] == 3:
            return conn.info().get_content_charset()
        else:
            encoding = "latin-1" # Default value
            for item in conn.info().getplist():
                pair = item.split("=")
                if len(pair) == 2 and pair[0].strip() == "charset":
                    encoding = pair[1].strip()
        return encoding

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
        while True:
            logging.info ("Running %s", " ".join(command))
            (rc, stdout, stderr) = run_command (command)
            if rc == 0:
                return True
            logging.error("Download of %s failed. Stderr:\n%s" % (url, stderr))
            logging.error("Waiting %s seconds before retrying", wait)
            time.sleep(wait)
        
    def curl_get_file(self, url, wait, dlpath, cafile=None):
        """ Download via curl
        
        Like wget_file(), this is used as a fall-back.
        """
        command = ["curl", "--silent", "--show-error", "--fail", "--output", dlpath]
        if cafile:
            command += ["--cacert", cafile]
        command.append(url)
        while True:
            logging.info ("Running %s", " ".join(command))
            (rc, stdout, stderr) = run_command(command)
            if rc == 0:
                return True
            logging.error("Download of %s failed. Stderr:\n%s" % (url, stderr))
            logging.error("Waiting %s seconds before retrying", wait)
            time.sleep(wait)
        
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
        # acutally checking the certificate
        # This is a rather ugly hack. It would be nicer to copy the required
        # parts from a fully functional SSL library. TODO.
        if url.startswith("https:") and sys.version_info[0] == 2:
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
                slept = 0
                timeout = 60
                while slept < timeout and os.path.getsize(dlpath) == 0:
                    logging.warning("Sleeping until %s contains data", dlpath)
                    time.sleep(1)
                    slept += 1
                if slept == timeout:
                    logging.warning("Slept %d seconds, %s still has no data", 
                                    timeout, dlpath)
                return
            else:
                raise
        fcntl.flock(fd, fcntl.LOCK_EX)
        outfile = os.fdopen(fd, "wb")
        shutil.copyfileobj (request, outfile)
        fcntl.flock(fd, fcntl.LOCK_UN)
        outfile.close() # This should also close the fd
        request.close()
    
    def get_missing_file(self, urlpath, filename, checksum=None,
                         options=None):
        # print('get_missing_file(%s, %s, %s)' % (urlpath, filename, checksum))
        if os.path.isfile(filename):
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
            try:
                self.get_file(urlpath, filename, options)
            except urllib_error.HTTPError as error:
                logging.error (str(error))
                return False
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
        for (filename, checksum) in self.workunit.get("FILE", []) + \
                self.workunit.get("EXECFILE", []):
            templ = Template(filename)
            archname = templ.safe_substitute({"ARCH": self.settings["ARCH"]})
            dlname = templ.safe_substitute({"ARCH": ""})
            dlpath = os.path.join(self.settings["DLDIR"], dlname)
            if not self.get_missing_file (archname, dlpath, checksum):
                return False
            # Try to lock the file once to be sure that download has finished
            # if another wuclient is doing the downloading
            with open(dlpath) as file_to_lock:
                fcntl.flock(file_to_lock, fcntl.LOCK_SH)
                fcntl.flock(file_to_lock, fcntl.LOCK_UN)
            
            if filename in dict(self.workunit.get("EXECFILE", [])):
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
        if processor.failedcommand:
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
    def do_checksum(filename, checksum = None):
        """ Computes the SHA1 checksum for a file. If checksum is None, returns 
            the computed checksum. If checksum is not None, return whether the
            computed SHA1 sum and checksum agree """
        blocksize = 65536
        sha1hash = hashlib.sha1() # pylint: disable=E1101
        infile = open(filename, "rb")
        fcntl.flock(infile, fcntl.LOCK_SH)
        
        data = infile.read(blocksize)
        while data:
            sha1hash.update(data)
            data = infile.read(blocksize)
        fcntl.flock(infile, fcntl.LOCK_UN)
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


def get_ssl_certificate(server, port=443):
    cert = ssl.get_server_certificate((server, int(port)),
                                      ssl_version=ssl.PROTOCOL_SSLv23,
                                      ca_certs=None)
    return cert

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

    # Create download and working directories if they don't exist
    if not os.path.isdir(SETTINGS["DLDIR"]):
        os.makedirs(SETTINGS["DLDIR"])
    if not os.path.isdir(SETTINGS["WORKDIR"]):
        os.makedirs(SETTINGS["WORKDIR"])

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

    generator_bug_type = getattr(FixedBytesGenerator, "bug_type", None)
    if generator_bug_type is BUGGY_MIMEENCODER1:
        logging.info("Using work-around #1 for buggy BytesGenerator")
    elif generator_bug_type is BUGGY_MIMEENCODER2:
        logging.info("Using work-around #2 for buggy BytesGenerator")

    (scheme, netloc) = urlparse(SETTINGS["SERVER"])[0:2]
    if not SETTINGS["CERTSHA1"] is None and scheme != "https":
        logging.warn("Option -certsha1 makes sense only with an https URL, ignoring it.")
    elif SETTINGS["CERTSHA1"] is None and scheme == "https":
        logging.warn("An https URL was given but no -certsha1 option, NO SSL VALIDATION WILL BE PERFORMED.")
    elif not SETTINGS["CERTSHA1"] is None and scheme == "https":
        certfilename = os.path.join(SETTINGS["DLDIR"], "server.%s.pem" % SETTINGS["CERTSHA1"][0:8])
        certfile_exists = os.path.isfile(certfilename)
        if certfile_exists:
            logging.info("Using certificate file stored in %s", certfilename)
            with open(certfilename, 'r') as certfile:
                cert = certfile.read()
        else:
            logging.info("Downloading certificate from %s", netloc)
            address_port = netloc.split(":")
            cert = get_ssl_certificate(*address_port)
        bin_cert = ssl.PEM_cert_to_DER_cert(cert)
        sha1hash = hashlib.sha1()
        sha1hash.update(bin_cert)
        cert_sha1 = sha1hash.hexdigest()
        logging.debug("Certificate has SHA1 %s", cert_sha1)
        if not cert_sha1.lower() == SETTINGS["CERTSHA1"].lower():
            logging.critical("Server certificate's SHA1 hash (%s) differs "
                             "from hash specified on command line (%s). "
                             "Aborting.", cert_sha1, SETTINGS["CERTSHA1"])
            logging.critical("Possible reason: several factorizations with "
                             "same download directory")
            sys.exit(1)
        logging.info("Certificate SHA1 hash matches")
        if not certfile_exists:
            logging.info("Writing certificate to file %s", certfilename)
            # FIXME: Set umask first?
            with open(certfilename, 'w') as certfile:
                certfile.write(cert)
        SETTINGS["CERTFILE"] = certfilename
        
        # Can we download with HTTPS at all?
        if sys.version_info[0] == 2:
            try:
                HAVE_WGET = run_command(["wget", "-V"])[0] == 0
            except OSError:
                pass
            if not HAVE_WGET:
                try:
                    HAVE_CURL = run_command(["curl", "-V"])[0] == 0
                except OSError:
                    pass
            if not HAVE_WGET and not HAVE_CURL:
                logging.critical("HTTPS requested, but not implemented in "
                        "Python 2 and can't find wget or curl as fall-back. "
                        "Aborting.")
                sys.exit(1)

    # Daemonize before we open the log file, or we close the log file's 
    # file descriptor
    if options.daemon:
        create_daemon(keepfd=None if logfile is None else [logfile.fileno()])

    client_ok = True
    while client_ok:
        try:
            client = WorkunitClient(settings = SETTINGS)
            client_ok = client.process()
        except BaseException:
            logging.exception("Exception occurred")
            break
