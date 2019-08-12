import os
import sys
import subprocess
import logging
import re
import cadoprograms
import cadologger
import cadotask
from os import environ as os_environ
#t=nb()
#print(t.taskn)
logger = logging.getLogger("Command")

def shellquote(s, first=True):
    ''' Quote a command line argument
    
    Currently does it the hard way: encloses the argument in single
    quotes, and escapes any single quotes that are part of the argument
    '''
    # If only characters that are known to be shell-safe occur, don't quote
    if re.match("^[a-zA-Z0-9_.:@/+-]*$", s):
        return s
    # An equals sign is safe, except in the first word of the command line
    if not first and re.match("^[a-zA-Z0-9_.:@/+-]+=[a-zA-Z0-9_.:@/+-]+$", s):
        return s
    return "'" + s.replace("'", "'\\''") + "'"

class Command(object):
    ''' Represents a running subprocess
    
    The subprocess is started when the instance is initialised
    '''
    @staticmethod
    def is_same(a,b):
        """ Return true if a and b are the same object and not None """
        return (not a is None and a is b)
    
    @staticmethod
    def _open(fn, is_output, appending = False):
        if fn is None:
            assert not appending
            return subprocess.PIPE if is_output else None
        if is_output:
            mode = "a" if appending else "w"
        else:
            assert not appending
            mode = "r"
        return open(fn, mode)
    
    @staticmethod
    def _close_or_not(fd, ref):
        if not fd in [None, subprocess.PIPE, subprocess.STDOUT, ref]:
            fd.close()
    
    def __init__(self, program, *args, **kwargs):
        self.program = program
        progargs = self.program.make_command_array()
        
        # Each of stdin, stdout, stderr is one of these three:
        # - subprocess.PIPE, if program.std* is None
        # - subprocess.STDOUT for stderr if program.stdout and program.stdin
        #   are the same object
        # - a newly opened file handle which needs to be closed later
        (stdin, (stdout, append_out), (stderr, append_err)) = \
            self.program.get_stdio()
        
        self.stdin = self._open(stdin, False)
        if type(stdout) == str or stdout is None:
            self.stdout = self._open(stdout, True, append_out)
            if self.is_same(stdout, stderr):
                self.stderr = subprocess.STDOUT
            else:
                self.stderr = self._open(stderr, True, append_err)
            self.child = subprocess.Popen(progargs, *args, stdin=self.stdin, stdout=self.stdout, stderr=self.stderr, **kwargs)
        else:
            self.stdout = None
            self.stderr = self._open(stderr, True, append_err)
            self.child = subprocess.Popen(progargs, *args, stdin=self.stdin,
                stdout=subprocess.PIPE, stderr=self.stderr, **kwargs)
            for x in self.child.stdout:
                stdout.filter(x.decode("utf-8"))


        cmdline = self.program.make_command_line()
        logger.cmd(cmdline, self.child.pid)
    
    def wait(self):
        ''' Wait for command to finish executing, capturing stdout and stderr
        in output tuple
        '''
        (stdout, stderr) = self.child.communicate()
        if self.child.returncode != 0:
            logger.warning("Process with PID %d finished with return code %u",
                           self.child.pid, self.child.returncode)
        else:
            logger.debug("Process with PID %d finished successfully",
                         self.child.pid)
        if stdout:
            logger.debug("Process with PID %d stdout: %s", 
                         self.child.pid, stdout)
        if stderr:
            logger.debug("Process with PID %d stderr: %s", 
                         self.child.pid, stderr)
        
        self._close_or_not(self.stdin, self.program.stdin)
        self._close_or_not(self.stdout, self.program.stdout)
        self._close_or_not(self.stderr, self.program.stderr)
        
        return (self.child.returncode, stdout, stderr)

PRINTED_SSHAUTH_WARNING = False

class RemoteCommand(Command):
    def __init__(self, program, host, parameters, *args, **kwargs):
        # We use a make_command_line() instead of make_command_array() so that,
        # e.g., stdio redirection to files specified in program can be added
        # to the command line with and the redirection happens on the remote
        # host
        global PRINTED_SSHAUTH_WARNING
        if not PRINTED_SSHAUTH_WARNING and not "SSH_AUTH_SOCK" in os_environ:
            logger.warning("No SSH_AUTH_SOCK shell environment variable found.")
            logger.warning("Make sure to set up an ssh key and ssh-agent to "
                        "avoid ssh asking for a passphrase.")
            logger.warning("See documentation for ssh-keygen and ssh-agent for "
                        "details.")
            PRINTED_SSHAUTH_WARNING = True
        cmdline = program.make_command_line()
        progparams = parameters.myparams(cadoprograms.SSH.get_accepted_keys(),
                                         cadoprograms.SSH.name)
        ssh = cadoprograms.SSH(host,
                               "env", "sh", "-c", shellquote(cmdline),
                               **progparams)
        super().__init__(ssh, *args, **kwargs)

class SendFile(Command):
    rsync_options = []
    def __init__(self, localfile, hostname, hostpath, parameters, *args,
                 **kwargs):
        if hostname == "localhost":
            target = hostpath
        else:
            target = hostname + ":" + hostpath
        rsync = cadoprograms.RSync(localfile, target, **parameters)
        super().__init__(rsync, *args, **kwargs)

if __name__ == '__main__':
    import cadoparams

    cadologger.init_test_logger()
    logger.setLevel(logging.NOTSET)

    ls_program = cadoprograms.Ls("/", long=True)
    c = Command(ls_program)
    (rc, out, err) = c.wait()
    if out:
        print("Stdout: " + str(out, encoding="utf-8"))
    if err:
        print("Stderr: " + str(err, encoding="utf-8"))
    del(c)

    ls_program = cadoprograms.Ls("/", stdout = "ls.out", long=True)
    ssh_parameters = cadoparams.Parameters({"verbose": False})
    c = RemoteCommand(ls_program, "localhost", ssh_parameters, [])
    (rc, out, err) = c.wait()
    print("Stdout: " + str(out, encoding="utf-8"))
    print("Stderr: " + str(err, encoding="utf-8"))
    del(c)

    c = SendFile("/users/caramel/kruppaal/ls.out", "quiche", "/tmp/foo", {})
    (rc, out, err) = c.wait()
    print("Stdout: " + str(out, encoding="utf-8"))
    print("Stderr: " + str(err, encoding="utf-8"))
    del(c)
