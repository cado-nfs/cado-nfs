import os
import platform
import abc
import inspect
import hashlib
import logging
import cadocommand
import cadologger
import cadoparams


class InspectType(type):
    """ Meta-class that adds a class attribute "init_signature" with the
    signature of the __init__() method, as produced by
    inspect.getfullargspec()
    """
    def __init__(cls, name, bases, dct):
        """
        >>> def f(a:"A", b:"B"=1, *args:"*ARGS", c:"C", d:"D"=3, **kwargs:"**KWARGS"):
        ...    pass
        >>> i = inspect.getfullargspec(f)
        >>> i.args
        ['a', 'b']
        >>> i.varargs
        'args'
        >>> i.varkw
        'kwargs'
        >>> i.defaults
        (1,)
        >>> i.kwonlyargs
        ['c', 'd']
        >>> i.kwonlydefaults
        {'d': 3}
        >>> i.annotations == {'a': 'A', 'b': 'B', 'c': 'C', 'd': 'D', 'args': '*ARGS', 'kwargs': '**KWARGS'}
        True
        """
        super().__init__(name, bases, dct)
        # inspect.getfullargspec() produces an object with attributes:
        # args, varargs, kwonlyargs, annotations (among others)
        # where args is a list that contains the names of positional parameters
        # (including those with default values), varargs contains the name of
        # the list of the variable-length positional parameters (i.e., the name
        # of the * catch-all), and kwonlyargs contains a list of keyword-only
        # parameters.
        cls.init_signature = inspect.getfullargspec(cls.__init__)


class Option(object, metaclass=abc.ABCMeta):
    ''' Base class for command line options that may or may not take parameters
    '''
    # TODO: decide whether we need '-' or '/' as command line option prefix
    # character under Windows. Currently: use "-'
    if False and platform.system() == "Windows":
        prefix = '/'
    else:
        prefix = '-'

    def __init__(self, arg=None, prefix=None, is_input_file=False,
                 is_output_file=False, checktype=None):
        """ Define a mapping from a parameter name and value to command line
        parameters.

        arg is the command line parameter to which this should map, e.g.,
            'v' or 'thr'. If not specified, the default name supplied to the
            map() method will be used.
        prefix is the command line parameter prefix to use; the default is the
            class variable which currently is '-'. Some programs, e.g. bwc.pl,
            want some parameters with a different prefix, e.g., ':complete'
        is_input_file must be set to True if this parameter gives the filename
            of an input file to the command. This will be used to generate FILE
            lines in workunits.
        is_output_file must be set to True if this parameter gives the filename
            of an output file of the command. This will be used to generate
            RESULT lines in workunits.
        if checktype is specified, map() will assert that the value is an
            instance of checktype.
        """
        self.arg = arg
        if prefix:
            self.prefix = prefix # hides the class variable
        self.is_input_file = is_input_file
        self.is_output_file = is_output_file
        self.checktype = checktype
        self.defaultname = None

    def set_defaultname(self, defaultname):
        """ Sets default name for the command line parameter. It must be specified
        before calling map().

        The reason for this is that command line parameters can be specified
        in the constructor of Program sub-classes as, e.g.,
          class Las(Program):
            def __init__(lpba : Parameter()=None):
        and the name "lpba" of the Las constructor's parameter should also be
        the default name of command line parameter. However, the Parameter()
        gets instantiated when the ":" annotation gets parsed, i.e., at the 
        time the class definition of Las is parsed, and the fact that this
        Parameter instance will act as an annotation to the "lpba" parameter
        is not known to the Parameter() instance. I.e., at its instantiation,
        the Parameter instance cannot tell to which parameter it will belong.

        This information must be filled in later, via set_defaultname(), which
        is called from Program.__init__(). It uses introspection to find out
        which Option objects belong to which Program constructor parameter, and
        fills in the constructor parameters' name via set_defaultname("lpba").

        If the Option constructor had received an arg parameter, then that is
        used instead of the defaultname, which allows for using different names
        for the Program constructor's parameter and the command line parameter,
        such as in
          __init__(threads: Parameter("t")):
        """
        self.defaultname = defaultname

    def get_arg(self):
        assert not self.defaultname is None
        return self.defaultname if self.arg is None else self.arg

    def get_checktype(self):
        return self.checktype

    def map(self, value):
        """ Public method that converts the Option instance into an array of
        strings with command line parameters. It also checks the type, if
        checktype was specified in the constructor.
        """
        assert not self.defaultname is None
        if not self.checktype is None:
            # If checktype is float, we allow both float and int as the type of
            # value. Some programs accept parameters that are integers (such as
            # admax) in scientific notation which is typed as floating point,
            # thus we want to allow float parameters to be given in both
            # integer and float type, so that passing, e.g., admax as an int
            # does not trip the assertion
            if self.checktype is float and type(value) is int:
                # Write it to command line as an int
                pass
            elif self.checktype is int and type(value) is float:
                # Can we convert this float to an int without loss?
                if float(int(value)) == value:
                    # Yes, convert it and write to command line as an int
                    value = int(value)
                else:
                    raise ValueError("Cannot convert floating-point value %s "
                                     "for parameter %s to an int without loss" %
                                     (value, self.defaultname))
            elif not isinstance(value, self.checktype):
                raise TypeError("Value %s for parameter %s is of type %s, but "
                                "checktype requires %s" %
                                (repr(value), self.defaultname,
                                 type(value).__name__,
                                 self.checktype.__name__))
        return self._map(value)

    @abc.abstractmethod
    def _map(self, value):
        """ Private method that does the actual translation to an array of string
        with command line parameters. Subclasses need to implement it.

        A simple case of the Template Method pattern.
        """
        pass

class PositionalParameter(Option):
    ''' Positional command line parameter '''
    def _map(self, value):
        return [str(value)]

class Parameter(Option):
    ''' Command line option that takes a parameter '''
    def _map(self, value):
        return [self.prefix + self.get_arg(), str(value)]

class ParameterEq(Option):
    ''' Command line option that takes a parameter, used on command line in
    the form of "key=value" '''
    def _map(self, value):
        return ["%s=%s" % (self.get_arg(), value)]

class Toggle(Option):
    ''' Command line option that does not take a parameter.
    value is interpreted as a truth value, the option is either added or not
    '''
    def __init__(self, arg=None, prefix=None, is_input_file=False,
                 is_output_file=False):
        """ Overridden constructor that hard-codes checktype=bool """
        super().__init__(arg=arg, prefix=prefix, is_input_file=is_input_file,
                 is_output_file=is_output_file, checktype=bool)

    def _map(self, value):
        if value is True:
            return [self.prefix + self.get_arg()]
        elif value is False:
            return []
        else:
            raise ValueError("Toggle.map() for %s requires a boolean "
                             "type argument" % self.get_arg())

class Sha1Cache(object):
    """ A class that computes SHA1 sums for files and caches them, so that a
    later request for the SHA1 sum for the same file is not computed again.
    File identity is checked via the file's realpath and the file's inode, size,
    and modification time.
    """
    def __init__(self):
        self._sha1 = {}

    @staticmethod
    def _read_file_in_blocks(file_object):
        blocksize = 65536
        while True:
            data = file_object.read(blocksize)
            if not data:
                break
            yield data

    def get_sha1(self, filename):
        realpath = os.path.realpath(filename)
        stat = os.stat(realpath)
        file_id = (stat.st_ino, stat.st_size, stat.st_mtime)
        # Check whether the file on disk changed
        if realpath in self._sha1 and not self._sha1[realpath][1] == file_id:
            logger = logging.getLogger("Sha1Cache")
            logger.warn("File %s changed! Discarding old SHA1 sum", realpath)
            del(self._sha1[realpath])
        if not realpath in self._sha1:
            logger = logging.getLogger("Sha1Cache")
            logger.debug("Computing SHA1 for file %s", realpath)
            sha1 = hashlib.sha1()
            with open(realpath, "rb") as inputfile:
                for data in self._read_file_in_blocks(inputfile):
                    sha1.update(data)
            self._sha1[realpath] = (sha1.hexdigest(), file_id)
            logger.debug("SHA1 for file %s is %s", realpath,
                         self._sha1[realpath])
        return self._sha1[realpath][0]

sha1cache = Sha1Cache()

if os.name == "nt":
    defaultsuffix = ".exe"
else:
    defaultsuffix = ""

IS_MINGW = "MSYSTEM" in os.environ

def translate_mingw_path(path):
    if path is None:
        return None
    (drive, path) = os.path.splitdrive(path)
    driveletter = drive.rstrip(":")
    # If a drive letter is given, we need to know the absolute path as we
    # don't handle per-drive current working directories
    dirs = path.split(os.sep)
    if driveletter:
        assert path[0] == "\\"
        dirs[1:1] = [driveletter]
    return '/'.join(dirs)

class Program(object, metaclass=InspectType):
    ''' Base class that represents programs of the CADO suite

    The Program base class is oblivious to how programs get input data,
    or how they provide output data. It does, however, handle redirecting
    stdin/stdout/stderr from/to files and accepts file names to/from which
    to redirect. This is done so that workunits can be generated from
    Program instances; in the workunit text, shell redirection syntax needs
    to be used to connect stdio to files, so a Program instance needs to be
    aware of which file names should be used for stdio redirection.

    Sub-classes correspond to individual programs (such as las) and must
    define these class variables:
    binary: a string with the name of the binary executable file
    name: a name used internally in the Python scripts for this program, must
      start with a letter and contain only letters and digits; if the binary
      file name is of this form, then that can be used
    params: a mapping that tells how to translate configuration parameters to
      command line options. The keys of the mapping are the keys as in the
      configuration file, e.g., "verbose" in a configuartion file line
      tasks.polyselect.verbose = 1
      The values of the mapping are instances of subclasses of Option,
      initialised with the command line parameter they should map to, e.g.,
      "verbose" should map to an instance Toggle("v"), and "threads" should
      map to an instance Parameter("t")

    >>> p = Ls()
    >>> p.make_command_line().replace("ls.exe", "/bin/ls")
    '/bin/ls'
    >>> p = Ls(stdout='foo')
    >>> p.make_command_line().replace("ls.exe", "/bin/ls")
    '/bin/ls > foo'
    >>> p = Ls(stderr='foo')
    >>> p.make_command_line().replace("ls.exe", "/bin/ls")
    '/bin/ls 2> foo'
    >>> p = Ls(stdout='foo', stderr='bar')
    >>> p.make_command_line().replace("ls.exe", "/bin/ls")
    '/bin/ls > foo 2> bar'
    >>> p = Ls(stdout='foo', stderr='foo')
    >>> p.make_command_line().replace("ls.exe", "/bin/ls")
    '/bin/ls > foo 2>&1'
    >>> p = Ls(stdout='foo', append_stdout=True)
    >>> p.make_command_line().replace("ls.exe", "/bin/ls")
    '/bin/ls >> foo'
    >>> p = Ls(stderr='foo', append_stderr=True)
    >>> p.make_command_line().replace("ls.exe", "/bin/ls")
    '/bin/ls 2>> foo'
    >>> p = Ls(stdout='foo', append_stdout=True, stderr='bar', append_stderr=True)
    >>> p.make_command_line().replace("ls.exe", "/bin/ls")
    '/bin/ls >> foo 2>> bar'
    >>> p = Ls(stdout='foo', append_stdout=True, stderr='foo', append_stderr=True)
    >>> p.make_command_line().replace("ls.exe", "/bin/ls")
    '/bin/ls >> foo 2>&1'
    >>> p = Ls('foo', 'bar')
    >>> p.make_command_line().replace("ls.exe", "/bin/ls")
    '/bin/ls foo bar'
    >>> p = Ls('foo', 'bar', long = True)
    >>> p.make_command_line().replace("ls.exe", "/bin/ls")
    '/bin/ls -l foo bar'
    '''

    path = '.'
    subdir = ""
    paramnames = ("execpath", "execsubdir", "execbin", "execsuffix",
        "runprefix")
    # These should be abstract properties, but we want to reference them as
    # class attributes, which properties can't. Ergo dummy variables
    binary = None

    # This class variable definition should not be here. It gets overwritten
    # when the InspectType meta-class creates the class object. The only purpose
    # is to make pylint shut up about the class not having an init_signature
    # attribute
    init_signature = None

    def __init__(self, options, stdin=None,
                 stdout=None, append_stdout=False, stderr=None,
                 append_stderr=False, background=False, execpath=None,
                 execsubdir=None, execbin=None, execsuffix=defaultsuffix,
                 runprefix=None):
        ''' Takes a dict of of command line options. Defaults are filled in
        from the cadoaprams.Parameters instance parameters.

        The stdin, stdout, and stderr parameters accept strings. If a string is
        given, it is interpreted as the file name to use for redirecting that
        stdio stream. In workunit generation, the file name is used for shell
        redirection.
        '''

        # The function annotations define how to convert each value in
        # self.parameters to command line arguments, but the annotations are
        # stored in an unordered dict. We would like to put the command line
        # parameters in a particular order, so that, for example, positional
        # parameters are placed correctly. We also need to handle a variable-
        # length positional parameter list.
        # That is, we need three lists:
        # - Individual positional parameters
        # - Variable-length positional parameter (e.g., list of filenames)
        # - Keyword parameters

        # First off, we set the default names for Option subclass instances.
        # Since the Option instances are always the same (the annotation values
        # are evaluated once at class creation time), we'll set the defaultnames
        # on the same Option instances again each time a Program subclass gets
        # instantiated, but that does not hurt.
        for (name, option) in self.__init__.__annotations__.items():
            option.set_defaultname(name)

        # "options" contains all local symbols defined in the sub-class'
        # __init__() method. This includes self, __class__, and a few others.
        # Filter them so that only those that correspond to annotated
        # parameters remain.
        filtered = self._filter_annotated_keys(options)

        # The default value for all __init__() parameters that correspond to
        # command line parameters is always None. Remove those entries.
        # We could in principle look up the default values from the .defaults
        # attribute given by inspect.getfullargspec(), but it's easier to
        # use the convention that the default is always None.
        filtered = {key:options[key] for key in filtered \
                    if not options[key] is None}

        self.parameters = filtered

        self.stdin = stdin
        self.stdout = stdout
        self.append_stdout = append_stdout
        self.stderr = stderr
        self.append_stderr = append_stderr
        self.runprefix = runprefix

        # If we are to run in background, we add " &" to the command line.
        # We require that stdout and stderr are redirected to files.
        self.background = background
        if self.background and (stdout is None or stderr is None):
            raise Exception("Programs to run in background must redirect "
                            "stdout and stderr")

        # Look for location of the binary executable at __init__, to avoid
        # calling os.path.isfile() multiple times
        path = str(execpath or self.path)
        subdir = str(execsubdir or self.subdir)
        binary = str(execbin or self.binary) + execsuffix
        execfile = os.path.normpath(os.sep.join([path, binary]))
        execsubfile = os.path.normpath(os.sep.join([path, subdir, binary]))
        self.suggest_subdir=dict()
        if execsubfile != execfile and os.path.isfile(execsubfile):
            self.execfile = execsubfile
            self.suggest_subdir[self.execfile]=subdir
        elif os.path.isfile(execfile):
            self.execfile = execfile
        else:
            self.execfile = binary
            # Raising an exception here makes it impossible to run doctests
            # for each Program subclass
            # raise Exception("Binary executable file %s not found" % execfile)

    @classmethod
    def _get_option_annotations(cls):
        """ Extract the elements of this class' __init__() annotations where
        the annotation is an Option instance
        """
        return {key:val for (key, val) in cls.__init__.__annotations__.items()
                if isinstance(val, Option)}

    @classmethod
    def _filter_annotated_options(cls, keys):
        """ From the list of keys given in "keys", return those that are
        parameters of the __init__() method of this class and annotated with an
        Option instance. Returns a dictionary of key:Option-instance pairs.
        """
        options = cls._get_option_annotations()
        return {key:options[key] for key in keys if key in options}

    @classmethod
    def _filter_annotated_keys(cls, keys):
        """ From the list of keys given in "keys", return those that are
        parameters of the __init__() method of this class and annotated with an
        Option instance. Returns a list of keys, where order w.r.t. the input
        keys is preserved.
        """
        options = cls._get_option_annotations()
        return [key for key in keys if key in options]

    @classmethod
    def get_accepted_keys(cls):
        """ Return all parameter file keys which can be used by this program,
        including those that don't directly map to command line parameters,
        like those specifying search paths for the binary executable file.

        This list is the union of annotated positional and keyword parameters to
        the constructor, plus the Program class defined parameters. Note that an
        *args parameter (stored in cls.init_signature.varargs) should not be
        included here; there is no way to specify a variable-length set of
        positional parameters in the parameter file.
        """
        # The fact that we want to exclude varargs is why we need the inspect
        # info here.
        
        # Get a map of keys to Option instances
        parameters = cls._filter_annotated_options(
            cls.init_signature.args + cls.init_signature.kwonlyargs)
        
        # Turn it into a map of keys to checktype
        parameters = {key:parameters[key].get_checktype() for key in parameters}
        
        # Turn checktypes that are not None into one-element lists.
        # The one-element list is treated by cadoparams.Parameters.myparams()
        # as a non-mandatory typed parameter.
        parameters = {key:None if checktype is None else [checktype]
            for key,checktype in parameters.items()}
        
        return cadoparams.UseParameters.join_params(parameters,
                                                    Program.paramnames)

    def get_stdout(self):
        return self.stdout

    def get_stderr(self):
        return self.stderr

    def get_stdio(self):
        """ Returns a 3-tuple with information on the stdin, stdout, and stderr
        streams for this program.

        The first entry is for stdin, and is either a file name or None.
        The second entry is for stdout and is a 2-tuple, whose first entry is
            either a file name or None, and whose second entry is True if we
            should append to the file and False if we should not.
        The third entry is like the second, but for stderr.
        """
        return (self.stdin,
                (self.stdout, self.append_stdout),
                (self.stderr, self.append_stderr)
            )

    def _get_files(self, is_output):
        """ Helper method to get a set of input/output files for this Program.
        This method returns the list of filenames without considering files for
        stdio redirection. If "is_output" evaluates to True, the list of output
        files is generated, otherwise the list of input files.
        """
        files = set()
        parameters = self.init_signature.args + self.init_signature.kwonlyargs
        if not self.init_signature.varargs is None:
            parameters.append(self.init_signature.varargs)
        parameters = self._filter_annotated_keys(parameters)
        for param in parameters:
            ann = self.__init__.__annotations__[param]
            if param in self.parameters and \
                    (ann.is_output_file if is_output else ann.is_input_file):
                if param == self.init_signature.varargs:
                    # vararg is a list; we need to convert each entry to string
                    # and append this list of strings to files
                    files |= set(map(str, self.parameters[param]))
                else:
                    files |= {str(self.parameters[param])}
        return files

    def get_input_files(self, with_stdio=True):
        """ Returns a list of input files to this Program instance. If
        with_stdio is True, includes stdin if such redirection is used.
        """
        input_files = self._get_files(is_output=False)
        if with_stdio and isinstance(self.stdin, str):
            input_files |= {self.stdin}
        return list(input_files)

    def get_output_files(self, with_stdio=True):
        """ Returns a list of output files. If with_stdio is True, includes
        files for stdout/stderr redirection if such redirection is used.
        """
        output_files = self._get_files(is_output=True)
        if with_stdio:
            for filename in (self.stdout, self.stderr):
                if isinstance(filename, str):
                    output_files |= {filename}
        return list(output_files)

    def get_exec_file(self):
        return self.execfile

    def get_exec_files(self):
        return [self.get_exec_file()]

    @staticmethod
    def translate_path(filename, filenametrans=None):
        if not filenametrans is None and str(filename) in filenametrans:
            return filenametrans[str(filename)]
        return str(filename)

    def make_command_array(self, filenametrans=None):
        # Begin command line with program to execute
        command = []
        if not self.runprefix is None:
            command=self.runprefix.split(' ')
        command.append(self.translate_path(self.get_exec_file(), filenametrans))

        # Add keyword command line parameters, then positional parameters
        parameters = self.init_signature.kwonlyargs + self.init_signature.args
        parameters = self._filter_annotated_keys(parameters)
        for key in parameters:
            if key in self.parameters:
                ann = self.__init__.__annotations__[key]
                value = self.parameters[key]
                # If this is an input or an output file name, we may have to
                # translate it, e.g., for workunits
                assert not (ann.is_input_file and ann.is_output_file)
                if ann.is_input_file:
                    value = self.translate_path(value, filenametrans)
                elif ann.is_output_file:
                    value = self.translate_path(value, filenametrans)
                command += ann.map(value)

        # Add positional command line parameters
        key = self.init_signature.varargs
        if not key is None:
            paramlist = self.parameters[key]
            ann = self.__init__.__annotations__.get(key, None)
            for value in paramlist:
                command += ann.map(value)
        return command

    def make_command_line(self, filenametrans=None):
        """ Make a shell command line for this program.

        If files are given for stdio redirection, the corresponding redirection
        tokens are added to the command line.
        """
        cmdarr = self.make_command_array(filenametrans=filenametrans)
        cmdline = " ".join(map(cadocommand.shellquote, cmdarr))
        if isinstance(self.stdin, str):
            translated = self.translate_path(self.stdin,
                                             filenametrans=filenametrans)
            cmdline += ' < ' + cadocommand.shellquote(translated)
        if isinstance(self.stdout, str):
            redir = ' >> ' if self.append_stdout else ' > '
            translated = self.translate_path(self.stdout,
                                             filenametrans=filenametrans)
            cmdline += redir + cadocommand.shellquote(translated)
        if not self.stderr is None and self.stderr is self.stdout:
            cmdline += ' 2>&1'
        elif isinstance(self.stderr, str):
            redir = ' 2>> ' if self.append_stderr else ' 2> '
            translated = self.translate_path(self.stderr,
                                             filenametrans=filenametrans)
            cmdline += redir + cadocommand.shellquote(translated)
        if self.background:
            cmdline += " &"
        return cmdline

    def make_wu(self, wuname):
        filenametrans = {}
        counters = {"FILE": 1, "EXECFILE": 1, "RESULT": 1}
        def append_file(wu, key, filename, with_checksum=True):
            assert not filename in filenametrans
            filenametrans[filename] = "${%s%d}" % (key, counters[key])
            counters[key] += 1
            wu.append('%s %s' % (key, os.path.basename(filename)))
            if with_checksum:
                wu.append('CHECKSUM %s' % sha1cache.get_sha1(filename))
            if self.suggest_subdir.get(filename, None):
                wu.append('SUGGEST_%s %s' % (key, self.suggest_subdir[filename]))
        
        workunit = ['WORKUNIT %s' % wuname]
        for filename in self.get_input_files():
            append_file(workunit, 'FILE', str(filename))
        for filename in self.get_exec_files():
            append_file(workunit, 'EXECFILE', str(filename))
        for filename in self.get_output_files():
            append_file(workunit, 'RESULT', str(filename), with_checksum=False)
        cmdline = self.make_command_line(filenametrans=filenametrans)
        workunit.append('COMMAND %s' % cmdline)
        workunit.append("") # Make a trailing newline
        return '\n'.join(workunit)


class Polyselect(Program):
    """
    >>> p = Polyselect(P=5, N=42, degree=4, verbose=True)
    >>> p.make_command_line().replace(defaultsuffix + " ", " ", 1)
    'polyselect -P 5 -N 42 -degree 4 -v'
    >>> p = Polyselect(P=5, N=42, degree=4, verbose=True)
    >>> p.make_command_line().replace(defaultsuffix + " ", " ", 1)
    'polyselect -P 5 -N 42 -degree 4 -v'
    """
    binary = "polyselect"
    name = binary
    subdir = "polyselect"

    def __init__(self, *,
                 P : Parameter(checktype=int)=None,
                 N : Parameter(checktype=int)=None,
                 degree : Parameter(checktype=int)=None,
                 verbose : Toggle("v")=None,
                 quiet : Toggle("q")=None,
                 threads : Parameter("t", checktype=int)=None,
                 admin : Parameter(checktype=int)=None,
                 admax : Parameter(checktype=int)=None,
                 incr : Parameter(checktype=int)=None,
                 nq : Parameter(checktype=int)=None,
                 maxtime : Parameter(checktype=float)=None,
                 out : Parameter(is_output_file=True)=None,
                 printdelay : Parameter("s", checktype=int)=None,
                 keep: Parameter(checktype=int)=None,
                 sopteffort: Parameter(checktype=int)=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)


class PolyselectRopt(Program):
    """
    >>> p = PolyselectRopt(ropteffort=5, inputpolys="foo.polys", verbose=True)
    >>> p.make_command_line().replace(defaultsuffix + " ", " ", 1)
    'polyselect_ropt -v -inputpolys foo.polys -ropteffort 5'
    """
    binary = "polyselect_ropt"
    name = binary
    subdir = "polyselect"

    def __init__(self, *,
                 verbose : Toggle("v")=None,
                 threads : Parameter("t", checktype=int)=None,
                 inputpolys : Parameter(is_input_file=True)=None,
                 ropteffort: Parameter(checktype=int)=None,
                 area : Parameter(checktype=float)=None,
                 Bf : Parameter(checktype=float)=None,
                 Bg : Parameter(checktype=float)=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class PolyselectGFpn(Program):
    binary = "polyselect_gfpn"
    name = binary
    subdir = "polyselect"

    def __init__(self, *,
                 verbose : Toggle("v")=None,
                 p: Parameter(checktype=int)=None,
                 n: Parameter(checktype=int)=None,
                 out: Parameter(is_output_file=True)=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class MakeFB(Program):
    """
    >>> p = MakeFB(poly="foo.poly", lim=1)
    >>> p.make_command_line().replace(defaultsuffix + " ", " ", 1)
    'makefb -poly foo.poly -lim 1'
    >>> p = MakeFB(poly="foo.poly", lim=1, maxbits=5, stdout="foo.roots")
    >>> p.make_command_line().replace(defaultsuffix + " ", " ", 1)
    'makefb -poly foo.poly -lim 1 -maxbits 5 > foo.roots'
    """
    binary = "makefb"
    name = binary
    subdir = "sieve"

    def __init__(self, *,
                 poly: Parameter(is_input_file=True),
                 lim: Parameter(checktype=int),
                 maxbits: Parameter(checktype=int)=None,
                 out: Parameter(is_output_file=True)=None,
                 side: Parameter(checktype=int)=None,
                 threads : Parameter("t", checktype=int)=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)


class FreeRel(Program):
    """
    >>> p = FreeRel(poly="foo.poly", renumber="foo.renumber", lpbr=1, lpba=2, out="foo.freerel")
    >>> p.make_command_line().replace(defaultsuffix + " ", " ", 1)
    'freerel -poly foo.poly -renumber foo.renumber -lpb0 1 -lpb1 2 -out foo.freerel'
    >>> p = FreeRel(poly="foo.poly", renumber="foo.renumber", lpbr=1, lpba=2, out="foo.freerel", badideals="foo.bad", pmin=123, pmax=234)
    >>> p.make_command_line().replace(defaultsuffix + " ", " ", 1)
    'freerel -poly foo.poly -renumber foo.renumber -lpb0 1 -lpb1 2 -out foo.freerel -badideals foo.bad -pmin 123 -pmax 234'
    """
    binary = "freerel"
    name = binary
    subdir = "sieve"
    def __init__(self, *,
                 poly: Parameter(is_input_file=True),
                 renumber: Parameter(is_output_file=True),
                 lpbr: Parameter("lpb0", checktype=int),
                 lpba: Parameter("lpb1", checktype=int),
                 out: Parameter(is_output_file=True),
                 badideals: Parameter(is_output_file=True)=None,
                 pmin: Parameter(checktype=int)=None,
                 pmax: Parameter(checktype=int)=None,
                 lcideals: Toggle() = None,
                 threads: Parameter("t")=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class Las(Program):
    binary = "las"
    name = binary
    subdir = "sieve"
    def __init__(self,
                 I: Parameter(checktype=int),
                 poly: Parameter(is_input_file=True),
                 q0: Parameter(checktype=int),
                 q1: Parameter(checktype=int)=None,
                 rho: Parameter(checktype=int)=None,
                 tdthresh: Parameter(checktype=int)=None,
                 bkthresh: Parameter(checktype=int)=None,
                 rlim: Parameter(checktype=int)=None,
                 alim: Parameter(checktype=int)=None,
                 lpbr: Parameter(checktype=int)=None,
                 lpba: Parameter(checktype=int)=None,
                 mfbr: Parameter(checktype=int)=None,
                 mfba: Parameter(checktype=int)=None,
                 rlambda: Parameter(checktype=float)=None,
                 alambda: Parameter(checktype=float)=None,
                 ncurves0: Parameter(checktype=int)=None,
                 ncurves1: Parameter(checktype=int)=None,
                 skewness: Parameter("S", checktype=float)=None,
                 verbose: Toggle("v")=None,
                 rpowlim: Parameter(checktype=int)=None,
                 apowlim: Parameter(checktype=int)=None,
                 factorbase: Parameter("fb", is_input_file=True)=None,
                 factorbase0: Parameter("fb0", is_input_file=True)=None,
                 factorbase1: Parameter("fb1", is_input_file=True)=None,
                 out: Parameter(is_output_file=True)=None,
                 threads: Parameter("t", checktype=int)=None,
                 ratq: Toggle()=None,
                 batch: Toggle()=None,
                 sqside: Parameter(checktype=int)=None,
                 dup: Toggle()=None,
                 galois: Parameter() = None,
                 allow_largesq: Toggle("allow-largesq")=None,
                 stats_stderr: Toggle("stats-stderr")=None,
                 # We have no checktype for parameters of the form <int>,<int>,
                 # so these are passed just as strings
                 traceab: Parameter() = None,
                 traceij: Parameter() = None,
                 traceNx: Parameter() = None,
                 # Let's make fbcache neither input nor output file. It should
                 # not be distributed to clients, nor sent back to the server.
                 # It's a local temp file, but re-used between different runs.
                 fbcache: Parameter("fbc")=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)


class Duplicates1(Program):
    binary = "dup1"
    name = binary
    subdir = "filter"
    def __init__(self,
                 *args: PositionalParameter(is_input_file=True),
                 prefix : Parameter(),
                 out: Parameter()=None,
                 outfmt: Parameter()=None,
                 bzip: Toggle("bz")=None,
                 only_ab: Toggle("ab")=None,
                 abhexa: Toggle()=None,
                 force_posix_threads: Toggle("force-posix-threads")=None,
                 only: Parameter(checktype=int)=None,
                 nslices_log: Parameter("n", checktype=int)=None,
                 lognrels: Parameter(checktype=int)=None,
                 filelist: Parameter(is_input_file=True)=None,
                 basepath: Parameter()=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)


class Duplicates2(Program):
    binary = "dup2"
    name = binary
    subdir = "filter"
    def __init__(self,
                 *args: PositionalParameter(is_input_file=True),
                 rel_count: Parameter("nrels", checktype=int),
                 renumber: Parameter(is_input_file=True),
                 filelist: Parameter(is_input_file=True)=None,
                 badidealinfo: Parameter(is_input_file=True)=None,
                 force_posix_threads: Toggle("force-posix-threads")=None,
                 dlp: Toggle("dl")=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class GaloisFilter(Program):
    binary = "filter_galois"
    name = binary
    subdir = "filter"
    def __init__(self,
                 *args: PositionalParameter(is_input_file=True),
                 nrels: Parameter(checktype=int),
                 poly: Parameter(is_input_file=True),
                 renumber: Parameter(is_input_file=True),
                 filelist: Parameter(is_input_file=True)=None,
                 basepath: Parameter()=None,
                 galois: Parameter("galois")=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)


class Purge(Program):
    binary = "purge"
    name = binary
    subdir = "filter"
    def __init__(self,
                 *args: PositionalParameter(is_input_file=True),
                 out: Parameter(is_output_file=True),
                 filelist: Parameter(is_input_file=True)=None,
                 basepath: Parameter()=None,
                 subdirlist: Parameter()=None,
                 nrels: Parameter(checktype=int)=None,
                 outdel: Parameter(is_output_file=True)=None,
                 keep: Parameter(checktype=int)=None,
                 col_minindex: Parameter("col-min-index", checktype=int)=None,
                 nprimes: Parameter("col-max-index", checktype=int)=None,
                 threads: Parameter("t", checktype=int)=None,
                 npass: Parameter(checktype=int)=None,
                 force_posix_threads: Toggle("force-posix-threads")=None,
                 required_excess: Parameter(checktype=float)=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class Merge(Program):
    binary = "merge"
    name = binary
    subdir = "filter"
    def __init__(self,
                 purged: Parameter("mat", is_input_file=True),
                 out: Parameter(is_output_file=True),
                 maxlevel: Parameter(checktype=int)=None,
                 keep: Parameter(checktype=int)=None,
                 skip: Parameter(checktype=int)=None,
                 target_density: Parameter(checktype=float)=None,
                 resume: Parameter(is_input_file=True)=None,
                 mkztype: Parameter(checktype=int)=None,
                 wmstmax: Parameter(checktype=int)=None,
                 force_posix_threads: Toggle("force-posix-threads")=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class MergeDLP(Program):
    binary = "merge-dl"
    name = binary
    subdir = "filter"
    def __init__(self,
                 purged: Parameter("mat", is_input_file=True),
                 out: Parameter(is_output_file=True),
                 maxlevel: Parameter(checktype=int)=None,
                 keep: Parameter(checktype=int)=None,
                 skip: Parameter(checktype=int)=None,
                 target_density: Parameter(checktype=float)=None,
                 resume: Parameter(is_input_file=True)=None,
                 mkztype: Parameter(checktype=int)=None,
                 wmstmax: Parameter(checktype=int)=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

# Todo: define is_input_file/is_output_file for remaining programs
class Replay(Program):
    binary = "replay"
    name = binary
    subdir = "filter"
    def __init__(self,
                 purged: Parameter()=None,
                 history: Parameter("his")=None,
                 index: Parameter()=None,
                 out: Parameter()=None,
                 for_msieve: Toggle()=None,
                 skip: Parameter(checktype=int)=None,
                 force_posix_threads: Toggle("force-posix-threads")=None,
                 bwcostmin: Parameter(checktype=int)=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class ReplayDLP(Program):
    binary = "replay-dl"
    name = binary
    subdir = "filter"
    def __init__(self,
                 purged: Parameter()=None,
                 ideals: Parameter()=None,
                 history: Parameter("his")=None,
                 index: Parameter()=None,
                 out: Parameter()=None,
                 skip: Parameter()=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class MagmaNmbrthry(Program):
    binary = "magma-nmbrthry-wrapper.sh"
    name = binary
    subdir = "scripts"
    def __init__(self,
                 poly: Parameter(),
                 N: Parameter("p"),
                 badidealinfo: Parameter("badinfofile"),
                 badideals: Parameter("badfile"),
                 ell: Parameter(),
                 **kwargs):
        super().__init__(locals(), **kwargs)


class MagmaLinalg(Program):
    binary = "magma-linalg-wrapper.sh"
    name = binary
    subdir = "scripts"
    def __init__(self,
                 ell: Parameter(),
                 nmaps: Parameter(),
                 sparsemat: Parameter(),
                 sm: Parameter(),
                 ker: Parameter()=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class BWC(Program):
    binary = "bwc.pl"
    name = "bwc"
    subdir = "linalg/bwc"
    def __init__(self,
                 complete: Toggle(prefix=":")=None,
                 dryrun: Toggle("d")=None,
                 verbose: Toggle("v")=None,
                 mpi: ParameterEq()=None,
                 lingen_mpi: ParameterEq()=None,
                 allow_zero_on_rhs: ParameterEq()=None,
                 threads: ParameterEq("thr")=None,
                 m: ParameterEq()=None,
                 n: ParameterEq()=None,
                 nullspace: ParameterEq()=None,
                 interval: ParameterEq()=None,
                 ys: ParameterEq()=None,
                 matrix: ParameterEq()=None,
                 rhs: ParameterEq()=None,
                 prime: ParameterEq()=None,
                 wdir: ParameterEq()=None,
                 mpiexec: ParameterEq()=None,
                 hosts: ParameterEq()=None,
                 hostfile: ParameterEq()=None,
                 interleaving: ParameterEq()=None,
                 shuffled_product: ParameterEq()=None,
                 bwc_bindir: ParameterEq()=None,
                 mm_impl: ParameterEq()=None,
                 cpubinding: ParameterEq()=None,
                 cantor_threshold: ParameterEq()=None,
                 lingen_threshold: ParameterEq()=None,
                 precmd: ParameterEq()=None,
                 # put None below for a random seed,
                 # or any value (for example 1) for a fixed seed
                 seed: ParameterEq()=None,
                 **kwargs):
        if os.name == "nt":
            kwargs.setdefault("runprefix", "perl.exe")
            kwargs.setdefault("execsuffix", "")
        if IS_MINGW:
            matrix = translate_mingw_path(matrix)
            wdir = translate_mingw_path(wdir)
            mpiexec = translate_mingw_path(mpiexec)
            hostfile = translate_mingw_path(hostfile)
            if bwc_bindir is None and "execpath" in kwargs:
                bwc_bindir = os.path.normpath(os.sep.join([kwargs["execpath"], self.subdir]))
            bwc_bindir = translate_mingw_path(bwc_bindir)
        super().__init__(locals(), **kwargs)

class SM(Program):
    binary = "sm"
    name = binary
    subdir = "filter"
    def __init__(self, *,
                 poly: Parameter(),
                 purged: Parameter(),
                 index: Parameter(),
                 out: Parameter(),
                 ell: Parameter(),
                 nsm: Parameter()=None,
                 threads: Parameter("t")=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)
 
class ReconstructLog(Program):
    binary = "reconstructlog-dl"
    name = binary
    subdir = "filter"
    def __init__(self, *,
                 ell: Parameter(),
                 threads: Parameter("mt")=None,
                 ker: Parameter("log", is_input_file=True),
                 dlog: Parameter("out"),
                 renumber: Parameter(),
                 poly: Parameter(),
                 purged: Parameter(),
                 ideals: Parameter(),
                 relsdel: Parameter(),
                 nrels: Parameter(),
                 partial: Toggle()=None,
                 nsm: Parameter(),
                 **kwargs):
        super().__init__(locals(), **kwargs)

class Descent(Program):
    binary = "descent.py"
    name = "descent"
    subdir = "scripts"
    def __init__(self, *,
                 target: Parameter(prefix="--"),
                 db: Parameter(prefix="--", is_input_file=True),
                 cadobindir: Parameter(prefix="--"),
                 descent_hint: Parameter("descent-hint", prefix="--",
                     is_input_file=True),
                 init_I: Parameter("init-I", prefix="--"),
                 init_ncurves: Parameter("init-ncurves", prefix="--"),
                 init_lpb: Parameter("init-lpb", prefix="--"),
                 init_lim: Parameter("init-lim", prefix="--"),
                 init_mfb: Parameter("init-mfb", prefix="--"),
                 init_tkewness: Parameter("init-tkewness", prefix="--"),
                 I: Parameter(prefix="--"),
                 lpb0: Parameter(prefix="--"),
                 lpb1: Parameter(prefix="--"),
                 mfb0: Parameter(prefix="--"),
                 mfb1: Parameter(prefix="--"),
                 lim0: Parameter(prefix="--"),
                 lim1: Parameter(prefix="--"),
                 ell: Parameter(prefix="--"),
                 **kwargs):
        super().__init__(locals(), **kwargs)


class Characters(Program):
    binary = "characters"
    name = binary
    subdir = "linalg"
    def __init__(self, *,
                 poly: Parameter(),
                 purged: Parameter(),
                 index: Parameter(),
                 heavyblock: Parameter(),
                 out: Parameter(),
                 wfile: Parameter("ker"),
                 lpbr: Parameter("lpb0"),
                 lpba: Parameter("lpb1"),
                 nchar: Parameter()=None,
                 nratchars: Parameter()=None,
                 threads: Parameter("t")=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class Sqrt(Program):
    binary = "sqrt"
    name = binary
    subdir = "sqrt"
    def __init__(self, *,
                 poly: Parameter(),
                 prefix: Parameter(),
                 purged: Parameter(),
                 index: Parameter(),
                 kernel: Parameter("ker"),
                 dep: Parameter()=None,
                 threads : Parameter("t", checktype=int)=None,
                 ab: Toggle()=None,
                 side0: Toggle()=None,
                 side1: Toggle()=None,
                 gcd: Toggle()=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class CadoNFSClient(Program):
    binary = "cado-nfs-client.py"
    name = "cado_nfs_client"
    subdir = ""
    def __init__(self,
                 server: Parameter(prefix='--'),
                 daemon: Toggle(prefix='--')=None,
                 keepoldresult: Toggle(prefix='--')=None,
                 nosha1check: Toggle(prefix='--')=None,
                 dldir: Parameter(prefix='--')=None,
                 workdir: Parameter(prefix='--')=None,
                 bindir: Parameter(prefix='--')=None,
                 clientid: Parameter(prefix='--')=None,
                 basepath: Parameter(prefix='--')=None,
                 getwupath: Parameter(prefix='--')=None,
                 loglevel: Parameter(prefix='--')=None,
                 postresultpath: Parameter(prefix='--')=None,
                 downloadretry: Parameter(prefix='--')=None,
                 logfile: Parameter(prefix='--')=None,
                 debug: Parameter(prefix='--')=None,
                 niceness: Parameter(prefix='--')=None,
                 wu_filename: Parameter(prefix='--')=None,
                 arch: Parameter(prefix='--')=None,
                 certsha1: Parameter(prefix='--')=None,
                 **kwargs):
        if os.name == "nt":
            kwargs.setdefault("runprefix", "python3.exe")
        super().__init__(locals(), **kwargs)

class SSH(Program):
    binary = "ssh"
    name = binary
    path = "/usr/bin"
    def __init__(self,
                 host: PositionalParameter(),
                 *args: PositionalParameter(),
                 compression: Toggle("C")=None,
                 verbose: Toggle("v")=None,
                 cipher: Parameter("c")=None,
                 configfile: Parameter("F")=None,
                 identity_file: Parameter("i")=None,
                 login_name: Parameter("l")=None,
                 port: Parameter("p")=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class RSync(Program):
    binary = "rsync"
    name = binary
    path = "/usr/bin"
    def __init__(self,
                 sourcefile: PositionalParameter(),
                 remotefile: PositionalParameter(),
                 **kwargs):
        super().__init__(locals(), **kwargs)

class Ls(Program):
    binary = "ls"
    name = binary
    path = "/bin"
    def __init__(self,
                 *args : PositionalParameter(),
                 long : Toggle('l')=None,
                 **kwargs):
        super().__init__(locals(), **kwargs)

class Kill(Program):
    binary = "kill"
    name = binary
    path = "/bin"
    def __init__(self,
                 *args: PositionalParameter(),
                 signal: Parameter("s"),
                 **kwargs):
        super().__init__(locals(), **kwargs)


if __name__ == '__main__':
    PROGRAM = Ls('foo', 'bar')
    CMDLINE = PROGRAM.make_command_line()
    print(CMDLINE)
