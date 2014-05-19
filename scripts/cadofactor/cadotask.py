#!/usr/bin/env python3
import re
import os.path
from fractions import gcd
import abc
import random
import time
import datetime
from collections import OrderedDict, defaultdict
from itertools import zip_longest
from math import log, sqrt
import logging
import socket
import gzip
import heapq
import patterns
import wudb
import cadoprograms
import cadoparams
import cadocommand
import wuserver
import workunit
from struct import error as structerror
from workunit import Workunit

# Pattern for floating-point number
RE_FP = r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
CAP_FP = "(%s)" % RE_FP


class Polynomial(list):
    """
    >>> p = Polynomial()
    >>> p.degree < 0
    True
    >>> p[0] = 1
    >>> p.degree == 0
    True
    """
    @property
    def degree(self):
        return len(self) - 1 if len(self) > 0 else float("-inf")

    def __setitem__(self, index, value):
        if index >= len(self):
            self.extend([0]*(index + 1 - len(self)))
        list.__setitem__(self, index, value)

    def __str__(self):
        xpow = ["", "*x"] + ["*x^%d" % i for i in range(2, len(self))]
        arr = ["%+d%s" % (self[idx], xpow[idx]) for idx in range(0, len(self))
               if self[idx]]
        poly = "".join(reversed(arr)).lstrip('+')
        poly = re.sub(r'\b1\*', "", poly)
        return poly

    def eval(self, x):
        """ Evaluate the polynomial at x """
        if len(self) == 0:
            return 0
        deg = self.degree
        value = self[deg]
        for i in range(deg):
            value = value * x + self[deg - i - 1]
        return value


class PolynomialParseException(Exception):
    """ Exception class for signaling errors during polynomial parsing """
    pass


class Polynomials(object):
    r""" A class that represents a polynomial
    
    >>> Polynomials([""])
    Traceback (most recent call last):
    cadotask.PolynomialParseException: No polynomials found
    >>> t="n: 1021\nc0: 1\nc1: -1\nc5: 1\nY0: 4\nY1: -1\nskew: 1.0\n"
    >>> p=Polynomials(t.splitlines())
    >>> str(p)
    'n: 1021\nskew: 1.0\nc0: 1\nc1: -1\nc5: 1\nY0: 4\nY1: -1\n# f(x) = x^5-x+1\n# g(x) = -x+4\n'
    >>> t="n: 1021\nc0: -1\nc1: 1\nc5: -1\nY0: -4\nY1: 1\nskew: 1.0\n"
    >>> p=Polynomials(t.splitlines())
    >>> str(p)
    'n: 1021\nskew: 1.0\nc0: -1\nc1: 1\nc5: -1\nY0: -4\nY1: 1\n# f(x) = -x^5+x-1\n# g(x) = x-4\n'
    """

    re_pol_f = re.compile(r"c(\d+)\s*:\s*(-?\d+)")
    re_pol_g = re.compile(r"Y(\d+)\s*:\s*(-?\d+)")
    re_Murphy = re.compile(r"\s*#\s*MurphyE\s*(?:\(.*\))?=%s$" % CAP_FP)
    re_lognorm = re.compile(r"\s*#\s*lognorm\s+%s" % CAP_FP)
    
    # Keys that can occur in a polynomial file, in their preferred ordering,
    # and whether the key is mandatory or not. The preferred ordering is used
    # when turning a polynomial back into a string.
    keys = OrderedDict(
        (
            ("n", (int, True)),
            ("skew", (float, False)),
            ("type", (str, False))
        ))
    
    def __init__(self, lines):
        """ Parse a polynomial file in the syntax as produced by polyselect2l
        """
        self.MurphyE = 0.
        self.lognorm = 0.
        self.params = {}
        polyf = Polynomial()
        polyg = Polynomial()

        def match_poly(line, poly, regex):
            match = regex.match(line)
            if match:
                (idx, coeff) = map(int, match.groups())
                if idx <= poly.degree and poly[idx]:
                    raise PolynomialParseException(
                        "Line '%s' redefines value %d" % line)
                poly[idx] = coeff
                return True
            return False

        for line in lines:
            # print ("Parsing line: >%s<" % line.strip())
            # If this is a comment line telling the Murphy E value,
            # extract the value and store it
            match = self.re_Murphy.match(line)
            if match:
                self.MurphyE = float(match.group(1))
                continue
            # If this is a comment line telling the lognorm,
            # extract the value and store it
            match = self.re_lognorm.match(line)
            if match:
                self.lognorm = float(match.group(1))
                continue
            # Drop comment, strip whitespace
            line2 = line.split('#', 1)[0].strip()
            # If nothing is left, process next line
            if not line2:
                continue
            # Try to parse polynomial coefficients
            if match_poly(line, polyf, self.re_pol_f) or \
                    match_poly(line, polyg, self.re_pol_g):
                continue
            # All remaining lines must be of the form "x: y"
            array = line2.split(":")
            if not len(array) == 2:
                raise PolynomialParseException("Invalid line '%s'" % line)
            key = array[0].strip()
            value = array[1].strip()
            
            if not key in self.keys:
                raise PolynomialParseException("Invalid key '%s' in line '%s'" %
                                               (key, line))
            if key in self.params:
                raise PolynomialParseException("Key %s in line %s has occurred "
                                               "before" % (key, line))
            (_type, isrequired) = self.keys[key]
            self.params[key] = _type(value)

        # If no polynomial was found at all (not even partial data), assume
        # that polyselect simply did not find anything in this search range
        if polyf.degree < 0 and polyg.degree < 0 and self.params == {} and \
                self.MurphyE == 0.:
            raise PolynomialParseException("No polynomials found")

        # Test that all required keys are there
        for (key, (_type, isrequired)) in self.keys.items():
            if isrequired and not key in self.params:
                raise PolynomialParseException("Key %s missing" % key)
        self.polyf = polyf
        self.polyg = polyg
        return

    def __str__(self):
        arr = ["%s: %s\n" % (key, self.params[key])
               for key in self.keys if key in self.params]
        arr += ["c%d: %d\n" % (idx, coeff) for (idx, coeff)
                in enumerate(self.polyf) if not coeff == 0]
        arr += ["Y%d: %d\n" % (idx, coeff) for (idx, coeff)
                in enumerate(self.polyg) if not coeff == 0]
        if not self.MurphyE == 0.:
            arr.append("# MurphyE = %g\n" % self.MurphyE)
        if not self.lognorm == 0.:
            arr.append("# lognorm %g\n" % self.lognorm)
        arr.append("# f(x) = %s\n" % str(self.polyf))
        arr.append("# g(x) = %s\n" % str(self.polyg))
        return "".join(arr)

    def __eq__(self, other):
        return self.polyf == other.polyf and self.polyg == other.polyg \
            and self.params == other.params

    def __ne__(self, other):
        return not (self == other)

    def create_file(self, filename):
        # Write polynomial to a file
        with open(str(filename), "w") as poly_file:
            poly_file.write(str(self))

    def getN(self):
        return self.params["n"]


class FilePath(object):
    """ A class that represents a path to a file, where the path should be
    somewhat relocateable.
    
    In particular, we separate the path to the working directory, and the file
    path relative to the working directory. For persistent storage in the DB,
    the path relative to the workdir should be used, whereas for any file
    accesses, the full path needs to be used.
    """

    def __init__(self, workdir, filepath):
        self.workdir = workdir.rstrip(os.sep)
        self.filepath = filepath

    def __str__(self):
        return "%s%s%s" % (self.workdir, os.sep, self.filepath)

    def get_wdir_relative(self):
        return self.filepath

    def isfile(self):
        return os.path.isfile(str(self))

    def isdir(self):
        return os.path.isdir(str(self))

    def mkdir(self, *, parent=False, mode=None):
        """ Creates a directory.
        
        parent acts much like the Unix mkdir's '-p' parameter: required parent
        directories are created if they don't exist, and no error is raised
        if the directory to be created already exists.
        If parent==True, a mode for the directory to be created can be specified
        as well.
        """
        if parent:
            # os.makedirs specifies 0o777 as the default value for mode,
            # thus we can't pass None to get the default value. We also
            # want to avoid hard-coding 0x777 as the default in this
            # method's signature, or using **kwargs magic. Thus we use
            # a default of None in this method, and pass the mode value
            # to makedirs only if it is not None.
            if mode is None:
                os.makedirs(str(self), exist_ok=True)
            else:
                os.makedirs(str(self), exist_ok=True, mode=mode)
        else:
            os.mkdir(str(self))
    def realpath(self):
        return os.path.realpath(str(self))
    def open(self, *args, **kwargs):
        return open(str(self), *args, **kwargs)


class WorkDir(object):
    """ A class that allows generating file and directory names under a
    working directory.
    
    The directory layout is as follows:
    The current project (i.e., the factorization) has a jobname, e.g.,
    "RSA512". Each task has a name, e.g., "sieving".
    A task can create various files under
    workdir/jobname.taskname.file
    or put them in a subdirectory
    workdir/jobname.taskname/file
    or, for multiple subdirectories,
    workdir/jobname.taskname/subdir/file
    
    >>> f = WorkDir("/foo/bar", "jobname", "taskname")
    >>> str(f.make_dirname()).replace(os.sep,'/')
    '/foo/bar/jobname.taskname/'
    >>> str(f.make_filename('file')).replace(os.sep,'/')
    '/foo/bar/jobname.taskname.file'
    >>> str(f.make_filename('file', use_subdir=True)).replace(os.sep,'/')
    '/foo/bar/jobname.taskname/file'
    >>> str(f.make_filename('file', use_subdir=True, subdir='subdir')).replace(os.sep,'/')
    '/foo/bar/jobname.taskname/subdir/file'
    """
    def __init__(self, workdir, jobname=None, taskname=None):
        self.workdir = str(workdir).rstrip(os.sep)
        self.jobname = jobname
        self.taskname = taskname
    
    def path_in_workdir(self, filename):
        return FilePath(self.workdir, filename)
    
    def make_filename2(self, jobname=None, taskname=None, filename=None):
        if jobname is None:
            jobname = self.jobname
        if taskname is None:
            taskname = self.taskname
        filename_arr = [s for s in [jobname, taskname, filename] if s]
        return FilePath(self.workdir, ".".join(filename_arr))
    
    def _make_path(self, extra):
        """ Make a path of the form: 'workdir/jobname.taskname''extra' """
        return self.path_in_workdir("%s.%s%s" % (self.jobname, self.taskname,
                                                 extra))
    
    def make_dirname(self):
        """ Make a directory name of the form workdir/jobname.taskname/ """
        return self._make_path(os.sep)
    
    def make_filename(self, name, use_subdir=False, subdir=None):
        """ If use_subdir is False, make a filename of the form
        workdir/jobname.taskname.name
        If use_subdir is True and subdir is None, make a filename of the form
        workdir/jobname.taskname/name
        If use_subdir is True and subdir is a string, make a filename of the
        form
        workdir/jobname.taskname/subdir/name
        """
        if use_subdir:
            if subdir:
                return self._make_path("%s%s%s%s" % (os.sep, subdir, os.sep,
                                                     name))
            else:
                return self._make_path("%s%s" % (os.sep, name))
        else:
            assert subdir is None
            return self._make_path(".%s" % name)


class Statistics(object):
    """ Class that holds statistics on program execution, and can merge two
    such statistics.
    """
    def __init__(self, conversions, formats):
        self.conversions = conversions
        self.stat_formats = formats
        self.stats = {}
    
    @staticmethod
    def typecast(values, types):
        """ Cast the values in values to the types specified in types """
        if type(types) is type:
            return [types(v) for v in values]
        else:
            return [t(v) for (v, t) in zip(values, types)]
    
    @staticmethod
    def _to_str(stat):
        """ Convert one statistic to a string """
        return " ".join(map(str, stat))
    
    @staticmethod
    def _from_str(string, types):
        """ Convert a string (probably from a state dict) to a statistic """
        return Statistics.typecast(string.split(), types)
    
    def from_dict(self, stats):
        """ Initialise values in self from the strings in the "stats"
        dictionary
        """
        for (key, types, defaults, combine, regex) in self.conversions:
            if key in stats:
                assert not key in self.stats
                self.stats[key] = self._from_str(stats.get(key, defaults),
                                                 types)
                assert not self.stats[key] is None
    
    def parse_line(self, line):
        """ Parse one line of program output and look for statistics.
        
        If they are found, they are added to self.stats.
        """
        for (key, types, defaults, combine, regex) in self.conversions:
            match = regex.match(line)
            if match:
                assert not key in self.stats
                # print (pattern.pattern, match.groups())
                # Optional groups that did not match are returned as None.
                # Skip over those so typecast doesn't raise TypeError
                groups = [group for group in match.groups() if not group is None]
                self.stats[key] = self.typecast(groups, types)
                assert not self.stats[key] is None
    
    def merge_one_stat(self, key, new_val, combine):
        if key in self.stats:
            self.stats[key] = combine(self.stats[key], new_val)
        else:
            self.stats[key] = new_val
        assert not self.stats[key] is None
        # print(self.stats)
    
    def merge_stats(self, new_stats):
        """ Merge the stats currently in self with the Statistics in
        "new_stats"
        """
        
        assert self.conversions == new_stats.conversions
        for (key, types, defaults, combine, regex) in self.conversions:
            if key in new_stats.stats:
                self.merge_one_stat(key, new_stats.stats[key], combine)
    
    def as_dict(self):
        return {key: self._to_str(self.stats[key]) for key in self.stats}
    
    def as_strings(self):
        """ Convert statistics to lines of output
        
        The self.stat_formats is an array, with each entry corresponding to
        a line that should be output.
        Each such entry is again an array, containing the format strings that
        should be used for the conversion of statistics. If a conversion
        fails with a KeyError or an IndexError, it is silently skipped over.
        This is to allow producing lines on which some statistics are not
        printed if the value is not known.
        """
        result = []
        for format_arr in self.stat_formats:
            line = []
            for format_str in format_arr:
                try:
                    line.append(format_str.format(**self.stats))
                except (KeyError, IndexError):
                    pass
            if line:
                result.append("".join(line))
        return result
    
    # Helper functions for processing statistics.
    # We can't make them @staticmethod or references are not callable
    def add_list(*lists):
        """ Add zero or more lists elementwise.
        
        Short lists are handled as if padded with zeroes.
        
        >>> Statistics.add_list([])
        []
        >>> Statistics.add_list([1])
        [1]
        >>> Statistics.add_list([1,2], [3,7])
        [4, 9]
        >>> Statistics.add_list([1,2], [3,7], [5], [3,1,4,1,5])
        [12, 10, 4, 1, 5]
        """
        return [sum(items) for items in zip_longest(*lists, fillvalue=0)]
    
    def weigh(samples, weights):
        return [sample * weight for (sample, weight) in zip(samples, weights)]
    
    def combine_mean(means, samples):
        """ From two lists, one containing values and the other containing
        the respective sample sizes (i.e., weights of the values), compute
        the combined mean (i.e. the weighted mean of the values).
        The two lists must have equal length.
        """
        assert len(means) == len(samples)
        total_samples = sum(samples)
        weighted_sum = sum(Statistics.weigh(means, samples))
        return [weighted_sum / total_samples, total_samples]
    
    def zip_combine_mean(*lists):
        """ From a list of 2-tuples, each tuple containing a value and a
        weight, compute the weighted average of the values.
        """
        for l in lists:
            assert len(l) == 2
        (means, samples) = zip(*lists)
        return Statistics.combine_mean(means, samples)
    
    def combine_stats(*stats):
        """ Computes the combined mean and std.dev. for the stats
        
        stats is a list of 3-tuples, each containing number of sample points,
        mean, and std.dev.
        Returns a 3-tuple with the combined number of sample points, mean,
        and std. dev.
        """
        
        # Samples is a list containing the first item (number of samples) of
        # each item of stats, means is list of means, stddevs is list of
        # std. dev.s
        for s in stats:
            assert len(s) == 3
        
        (samples, means, stddevs) = zip(*stats)
        
        (total_mean, total_samples) = Statistics.combine_mean(means, samples)
        # t is the E[X^2] part of V(X)=E(X^2) - (E[X])^2
        t = [mean**2 + stddev**2 for (mean, stddev) in zip(means, stddevs)]
        # Compute combined variance
        total_var = Statistics.combine_mean(t, samples)[0] - total_mean**2
        return [total_samples, total_mean, sqrt(total_var)]
    
    def test_combine_stats():
        """ Test function for combine_stats()
        
        >>> Statistics.test_combine_stats()
        True
        """
        
        from random import randrange
        
        def mean(x):
            return float(sum(x))/float(len(x))
        def var(x):
            E = mean(x)
            return mean([(a-E)**2 for a in x])
        def stddev(x):
            return sqrt(var(x))
        
        # Generate between 1 and 5 random integers in [1,100]
        lengths = [randrange(100) + 1 for i in range(randrange(5) + 1)]
        lengths = [1, 10]
        # Generate lists of random integers in [1,100]
        lists = [[randrange(100) for i in range(l)] for l in lengths]
        stats = [(length, mean(l), stddev(l))
            for (length, l) in zip(lengths, lists)]
        
        combined = []
        for l in lists:
            combined += l
        
        combined1 = Statistics.combine_stats(*stats)
        combined2 = [len(combined), mean(combined), stddev(combined)]
        if abs(combined1[2] - combined2[2]) > 0.2 * combined2[2]:
            print("lists = %r" % lists)
            print("combineds = %r" % combined)
            print("stats = %r" % stats)
            print("combined1 = %r" % combined1)
            print("combined2 = %r" % combined2)
            print(combined1[2], combined2[2])
            print(abs(combined1[2] / combined2[2] - 1))
        return combined1[0] == combined2[0] and \
                abs(combined1[1] / combined2[1] - 1) < 1e-10 and \
                abs(combined1[2] - combined2[2]) <= 1e-10 * combined2[2]
    
    def smallest_n(*lists, n=10):
        concat = []
        for l in lists:
            concat += l
        concat.sort()
        return concat[0:n]

    def parse_stats(self, filename):
        """ Parse statistics from the file with name "filename" and merge them
        into self
        
        Returns the newly parsed stats as a dictionary
        """
        new_stats = Statistics(self.conversions, self.stat_formats)
        with open(str(filename), "r") as inputfile:
            for line in inputfile:
                new_stats.parse_line(line)
        self.merge_stats(new_stats)
        return new_stats.as_dict()


class HasName(object, metaclass=abc.ABCMeta):
    @abc.abstractproperty
    def name(self):
        # The name of the task in a simple form that can be used as
        # a Python dictionary key, a directory name, part of a file name,
        # part of an SQL table name, etc. That pretty much limits it to
        # alphabetic first letter, and alphanumeric rest.
        pass

class HasTitle(object, metaclass=abc.ABCMeta):
    @abc.abstractproperty
    def title(self):
        # A pretty name for the task, will be used in screen output
        pass

class DoesLogging(HasTitle, metaclass=abc.ABCMeta):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.logger = logging.getLogger(self.title)

class MakesTablenames(HasName):
    @property
    def tablename_prefix(self):
        """ Prefix string for table names
        
        By default, the table name prefix is the name attribute, but this can
        be overridden
        """
        return self.name
    
    def make_tablename(self, extra=None):
        """ Return a name for a DB table """
        # Maybe replace SQL-disallowed characters here, like digits and '.' ?
        # Could be tricky to avoid collisions
        name = self.tablename_prefix
        if extra:
            name = name + '_' + extra
        wudb.check_tablename(name)
        return name

class HasState(MakesTablenames, wudb.HasDbConnection):
    """ Declares that the class has a DB-backed dictionary in which the class
    can store state information.
    
    The dicatonary is available as an instance attribute "state".
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        name = self.make_tablename()
        self.state = self.make_db_dict(name, connection=self.db_connection)


class FilesCreator(MakesTablenames, wudb.HasDbConnection, metaclass=abc.ABCMeta):
    """ A base class for classes that produce a list of output files, with
    some auxiliary information stored with each file (e.g., nr. of relations).
    This info is stored in the form of a DB-backed dictionary, with the file
    name as the key and the auxiliary data as the value.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        tablename = self.make_tablename("outputfiles")
        self.output_files = self.make_db_dict(tablename,
                                              connection=self.db_connection)
    
    def add_output_files(self, filenames, *, commit):
        """ Adds a dict of files to the list of existing output files """
        for filename in filenames:
            if filename in self.output_files:
                raise KeyError("%s already in output files table" % filename)
        self.output_files.update(filenames, commit=commit)
    
    def get_output_filenames(self, condition=None):
        """ Return output file names, optionally those that match a condition
        
        If a condition is given, it must be callable with 1 parameter and
        boolean return type; then only those filenames are returned where
        for the auxiliary data s (i.e., the value stored in the dictionary
        with the file name as key) satisfies condition(s) == True.
        """
        if condition is None:
            return list(self.output_files.keys())
        else:
            return [f for (f, s) in self.output_files.items() if condition(s)]
    
    def forget_output_filenames(self, filenames, *, commit):
        self.output_files.clear(filenames, commit=commit)


class BaseStatistics(object):
    """ Base class for HasStatistics and SimpleStatistics that terminates the
    print_stats() call chain.
    """
    def print_stats(self):
        pass


class HasStatistics(BaseStatistics, HasState, DoesLogging, metaclass=abc.ABCMeta):
    @property
    def stat_conversions(self):
        """ Sub-classes should override """
        return []

    @property
    def stat_formats(self):
        """ Sub-classes should override """
        return []

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.statistics = Statistics(self.stat_conversions, self.stat_formats)
        self.statistics.from_dict(self.state)

    def get_statistics_as_strings(self):
        """ Return the statistics collected so far as a List of strings.
        
        Sub-classes can override to add/remove/change strings.
        """
        return self.statistics.as_strings()

    def print_stats(self):
        stat_msgs = self.get_statistics_as_strings()
        if stat_msgs:
            self.logger.info("Aggregate statistics:")
            for msg in stat_msgs:
                self.logger.info(msg)
        super().print_stats()
    
    def parse_stats(self, filename, *, commit):
        new_stats = self.statistics.parse_stats(filename)
        self.logger.debug("Newly arrived stats: %s", new_stats)
        update = self.statistics.as_dict()
        self.logger.debug("Combined stats: %s", update)
        self.state.update(update, commit=commit)


class SimpleStatistics(BaseStatistics, HasState, DoesLogging, 
        metaclass=abc.ABCMeta):

    @abc.abstractproperty
    def programs(self):
        # A list of classes of Programs which this tasks uses
        pass

    def print_cpu_real_time(self, cputotal, realtotal, program):
        """ Print cpu and/or real time to logger """
        # Uses self only for access to the logger
        pairs = zip((cputotal, realtotal), ("cpu", "real"))
        usepairs = [pair for pair in pairs if pair[0]]
        if usepairs:
            printformat = "/".join(["%g"] * len(usepairs))
            usepairs = tuple(zip(*usepairs))
            timestr = '/'.join(usepairs[1])
            self.logger.info("Total %s time for %s: " + printformat,
                    timestr, program, *usepairs[0])
    
    def update_cpu_or_real_time(self, is_cpu, program, seconds, commit=True):
        """ Add seconds to the statistics of cpu time spent by program,
        and return the new total.
        """
        assert isinstance(program, str)
        key = "%s_%s" % ("cputime" if is_cpu else "realtime", program)
        total = self.state.get(key, 0.) + seconds
        # Update only if the value changed, i.e., if seconds != 0.
        if seconds:
            self.state.update({key: total}, commit=commit)
        return total

    def get_total_cpu_or_real_time(self, is_cpu):
        """ Return number of seconds of cpu time spent by all programs of
        this Task
        """
        total = 0.
        for program in self.programs:
            total += self.update_cpu_or_real_time(is_cpu, program.name, 0.,
                    commit=False)
        return total

    def print_stats(self):
        for program in self.programs:
            cputotal = self.update_cpu_or_real_time(True, program.name, 0.)
            realtotal = self.update_cpu_or_real_time(False, program.name, 0.)
            self.print_cpu_real_time(cputotal, realtotal, program.name)
        super().print_stats()


class Task(patterns.Colleague, SimpleStatistics, HasState, DoesLogging,
           cadoparams.UseParameters, metaclass=abc.ABCMeta):
    """ A base class that represents one task that needs to be processed.
    
    Sub-classes must define class variables:
    """
    
    # Properties that subclasses need to define
    @abc.abstractproperty
    def programs(self):
        # A list of classes of Programs which this tasks uses
        pass
    @abc.abstractproperty
    def progparam_override(self):
        # A list of lists of program parameters that get overridden by the
        # Task, such as admin for polyselect2l. Specifying these overridden
        # parameters here lets cadofactor print a warning message when the
        # user tries to supply them in the parameter file. Not specifying
        # them here may also lead to a program abort, when a parameter is
        # specified both in the progparams dictionary (as extracted via
        # .myparams()) and as a direct parameter as supplied by the Task.
        pass
    @abc.abstractproperty
    def paramnames(self):
        # Parameters that all tasks use
        return self.join_params(super().paramnames, 
            {"name": str, "workdir": str, "run": True})
    @property
    def param_nodename(self):
        return self.name
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        ''' Sets up a database connection and a DB-backed dictionary for
        parameters. Reads parameters from DB, and merges with hierarchical
        parameters in the parameters argument. Parameters passed in by
        parameters argument do not override values in the DB-backed
        parameter dictionary.
        '''
        
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.logger.debug("Enter Task.__init__(%s)",
                          self.name)
        self.logger.debug("state = %s", self.state)
        # Set default parameters for this task, if any are given
        self.params = self.parameters.myparams(self.paramnames)
        self.logger.debug("self.parameters = %s", self.parameters)
        self.logger.debug("params = %s", self.params)
        # Set default parameters for our programs
        self.progparams = []
        for prog, override in zip(self.programs, self.progparam_override):
            progparams = self.parameters.myparams(prog.get_accepted_keys(),
                                                  prog.name)
            for param in set(override) & set(progparams):
                self.logger.warn('Parameter "%s" for program "%s" is '
                                 'generated at run time and cannot be '
                                 'supplied through the parameter file',
                                 param, prog.name)
                del(progparams[param])
            
            self.progparams.append(progparams)
        # FIXME: whether to init workdir or not should not be controlled via
        # presence of a "workdir" parameter, but by class definition
        if "workdir" in self.params:
            self.workdir = WorkDir(self.params["workdir"], self.params["name"],
                               self.name)
        # Request mediator to run this task, unless "run" parameter is set
        # to false
        if self.params["run"]:
            self.send_notification(Notification.WANT_TO_RUN, None)
        self.logger.debug("Exit Task.__init__(%s)", self.name)
        return
    
    def translate_input_filename(self, filename):
        return filename

    def test_outputfile_exists(self, filename):
        return filename.isfile()
    
    @staticmethod
    def check_files_exist(filenames, filedesc, shouldexist):
        """ Check that the output files in "filenames" exist or don't exist,
        according to shouldexist.
        
        Raise IOError if any check fails, return None
        """
        for filename in filenames:
            if isinstance(filename, FilePath):
                exists = filename.isfile()
            else:
                exists = os.path.isfile(filename)
            if shouldexist and not exists:
                raise IOError("%s file %s does not exist" % (filedesc, filename))
            elif not shouldexist and exists:
                raise IOError("%s file %s already exists" % (filedesc, filename))
        return
    
    # These two function go together, one produces a workunit name from the
    # name of the factorization, the task name, and a task-provided identifier,
    # and the other function splits them again
    wu_paste_char = '_'
    wu_attempt_char = '#'
    def make_wuname(self, identifier, attempt=None):
        """ Generates a wuname from project name, task name, identifier, and
        attempt number.
        """
        assert not self.wu_paste_char in self.name # self.name is task name
        assert not self.wu_paste_char in identifier # identifier is, e.g., range string
        assert not self.wu_attempt_char in identifier
        wuname = self.wu_paste_char.join([self.params["name"], self.name,
                                          identifier])
        if not attempt is None:
            wuname += "%s%d" % (self.wu_attempt_char, attempt)
        return wuname
    
    def split_wuname(self, wuname):
        """ Splits a wuname into project name, task name, identifier, and
        attempt number.
        
        Always returns a list of length 4; if there is not attempt given in
        the wuname, then the last array entry is None

        >>> # Test many possible combinations of "_" and "#" occuring in names
        >>> # where these characters are allowed
        >>> class Klass():
        ...     params = {"name": None}
        ...     wu_paste_char = '_'
        ...     wu_attempt_char = '#'
        >>> inst = Klass()
        >>> from itertools import product
        >>> prod = product(*[["", "_", "#"]] * 4 + [["", "#"]]*2)
        >>> for sep in prod:
        ...     inst.params["name"] = "%s%sprojectname%s%s" % sep[0:4]
        ...     inst.name = "%staskname%s" % sep[4:6]
        ...     for attempt in [None, 2, 3]:
        ...         identifier = "identifier"
        ...         wuname = Task.make_wuname(inst, "identifier", attempt=attempt)
        ...         wu_split = Task.split_wuname(inst, wuname)
        ...         assert wu_split == [inst.params["name"], inst.name, identifier, attempt]
        """
        arr = wuname.rsplit(self.wu_paste_char, 2)
        assert len(arr) == 3
        attempt = None
        # Split off attempt number, if available
        if "#" in arr[2]:
            (arr[2], attempt) = arr[2].split('#')
            attempt = int(attempt)
        arr.append(attempt)
        return arr
    
    class ResultInfo(wudb.WuResultMessage):
        def __init__(self, wuid, rc, stdout, stderr, program, cmd_line, host):
            self.wuid = wuid
            self.rc = rc
            self.stdout = stdout if stdout else None
            self.stdoutfile = program.get_stdout()
            # stdout must be either in a string or in a file, but not both
            assert self.stdout is None or not self.stdoutfile
            self.stderr = stderr if stderr else None
            self.stderrfile = program.get_stderr()
            # stderr must be either in a string or in a file, but not both
            assert self.stderr is None or not self.stderrfile
            self.output_files = program.get_output_files(with_stdio=False)
            self.cmd_line = cmd_line
            self.host = host
        def get_wu_id(self):
            return self.wuid
        def get_output_files(self):
            return self.output_files
        def get_stdout(self, command_nr):
            assert command_nr == 0
            return self.stdout
        def get_stdoutfile(self, command_nr):
            assert command_nr == 0
            return self.stdoutfile
        def get_stderr(self, command_nr):
            assert command_nr == 0
            return self.stderr
        def get_stderrfile(self, command_nr):
            assert command_nr == 0
            return self.stderrfile
        def get_exitcode(self, command_nr):
            assert command_nr == 0
            return self.rc
        def get_command_line(self, command_nr):
            assert command_nr == 0
            return self.cmd_line
        def get_host(self):
            return self.host

    def log_failed_command_error(self, message, command_nr):
        host = message.get_host()
        host_msg = " run on %s" % host if host else ""
        self.logger.error("Program%s failed with exit code %d", 
                          host_msg, message.get_exitcode(command_nr))
        cmd_line = message.get_command_line(command_nr)
        if cmd_line:
            self.logger.error("Command line was: %s", cmd_line)
        stderr = message.read_stderr(command_nr)
        stderrfilename = message.get_stderrfile(command_nr)
        if stderrfilename:
            stderrmsg = " (stored in file %s)" % stderrfilename
        else:
            stderrmsg = ""
        if stderr:
            self.logger.error("Stderr output follows%s:\n%s", stderrmsg, stderr)

    def submit_command(self, command, identifier, commit=True, log_errors=False):
        ''' Run a command.
        Return the result tuple. If the caller is an Observer, also send
        result to updateObserver().
        '''
        wuname = self.make_wuname(identifier)
        process = cadocommand.Command(command)
        cputime_used = os.times()[2] # CPU time of child processes
        realtime_used = time.time()
        (rc, stdout, stderr) = process.wait()
        cputime_used = os.times()[2] - cputime_used
        realtime_used = time.time() - realtime_used
        self.update_cpu_or_real_time(True, command.name, cputime_used, False)
        self.update_cpu_or_real_time(False, command.name, realtime_used, commit)
        message = Task.ResultInfo(wuname, rc, stdout, stderr, command, 
                                  command.make_command_line(), "server")
        if rc != 0 and log_errors:
            self.log_failed_command_error(message, 0)

        if isinstance(self, patterns.Observer):
            # pylint: disable=E1101
            self.updateObserver(message)
        return message
    
    def filter_notification(self, message):
        wuid = message.get_wu_id()
        rc = message.get_exitcode(0)
        stdout = message.read_stdout(0)
        stderr = message.read_stderr(0)
        output_files = message.get_output_files()
        self.logger.message("%s: Received notification for wuid=%s, rc=%d, "
                            "output_files=[%s]",
                            self.name, wuid, rc, ", ".join(output_files))
        (name, task, identifier, attempt) = self.split_wuname(wuid)
        if name != self.params["name"] or task != self.name:
            # This notification is not for me
            self.logger.message("Notification %s is not for me", wuid)
            return
        self.logger.message("Notification %s is for me", wuid)
        if rc != 0:
            self.logger.debug("Return code is: %d", rc)
        if stdout:
            self.logger.debug("stdout is: %s", stdout)
        if stderr:
            self.logger.debug("stderr is: %s", stderr)
        if output_files:
            self.logger.message("Output files are: %s", ", ".join(output_files))
        return identifier
    
    def send_notification(self, key, value):
        """ Wrapper around Colleague.send_notification() that instantiates a
        Notification with self as the sender
        """
        notification = Notification(self, key, value)
        super().send_notification(notification)
    
    def send_request(self, key, *args):
        """ Wrapper around Colleague.send_request() that instantiates a
        Request with self as the sender
        """
        request = Request(self, key, *args)
        return super().send_request(request)
    
    def get_number_outstanding_wus(self):
        return 0
    
    def verification(self, wuid, ok, *, commit):
        pass

    def get_state_filename(self, key):
        if not key in self.state:
            return None
        return self.workdir.path_in_workdir(self.state[key])

    def make_std_paths(self, progname, do_increment=True):
        count = self.state.get("stdiocount", 0)
        if do_increment:
            count += 1
        did_increment = do_increment
        while True:
            try:
                stdoutname = "%s.stdout.%d" % (progname, count)
                stderrname = "%s.stderr.%d" % (progname, count)
                self.check_files_exist((stdoutname, stderrname), "stdio",
                                       shouldexist=False)
            except IOError:
                count += 1
                did_increment = True
                self.logger.warning("Stdout or stderr files with index %d "
                                    "already exist", count)
            else:
                break
        stdoutpath = self.workdir.make_filename(stdoutname)
        stderrpath = self.workdir.make_filename(stderrname)
        if did_increment:
            self.state["stdiocount"] = count
        return (stdoutpath, stderrpath)


class ClientServerTask(Task, wudb.UsesWorkunitDb, patterns.Observer):
    @abc.abstractproperty
    def paramnames(self):
        return self.join_params(super().paramnames,  
            {"maxwu": 10, 
             "wutimeout": 10800,  # Default: 3h
             "maxresubmit": 5, 
             "maxtimedout": 100, 
             "maxfailed": 100})
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.state.setdefault("wu_submitted", 0)
        self.state.setdefault("wu_received", 0)
        self.state.setdefault("wu_timedout", 0)
        self.state.setdefault("wu_failed", 0)
        assert self.get_number_outstanding_wus() >= 0
        self.send_notification(Notification.SUBSCRIBE_WU_NOTIFICATIONS, None)
    
    def submit_wu(self, wu, commit=True):
        """ Submit a WU and update wu_submitted counter """
        key = "wu_submitted"
        self.state.update({key: self.state[key] + 1}, commit=False)
        self.wuar.create(str(wu), commit=commit)
    
    def cancel_wu(self, wuid, commit=True):
        """ Cancel a WU and update wu_timedout counter """
        self.logger.debug("Cancelling: %s", wuid)
        key = "wu_timedout"
        maxtimedout = self.params["maxtimedout"]
        if not self.state[key] < maxtimedout:
            self.logger.error("Exceeded maximum number of timed out "
                              "workunits, maxtimedout=%d ", maxtimedout)
            raise Exception("Too many timed out work units")
        self.state.update({key: self.state[key] + 1}, commit=False)
        self.wuar.cancel(wuid, commit=commit)
    
    def submit_command(self, command, identifier, commit=True, log_errors=False):
        ''' Submit a workunit to the database. '''
        
        while self.get_number_available_wus() >= self.params["maxwu"]:
            self.wait()
        wuid = self.make_wuname(identifier)
        wutext = command.make_wu(wuid)
        for filename in command.get_exec_files() + command.get_input_files():
            basename = os.path.basename(filename)
            self.send_notification(Notification.REGISTER_FILENAME,
                                   {basename:filename})
        
        self.logger.info("Adding workunit %s to database", wuid)
        # self.logger.debug("WU:\n%s" % wutext)
        self.submit_wu(wutext, commit=commit)
        # Write command line to a file
        cmdline = command.make_command_line()
        client_cmd_filename = self.workdir.make_filename2(taskname="",
                                                          filename="wucmd")
        with client_cmd_filename.open("a") as client_cmd_file:
            client_cmd_file.write("# Command for work unit: %s\n%s\n" %
                                  (wuid, cmdline))
    
    def verification(self, wuid, ok, *, commit):
        """ Mark a workunit as received and update wu_received counter """
        ok_str = "ok" if ok else "not ok"
        self.logger.info("Marking workunit %s as %s", wuid, ok_str)
        assert self.get_number_outstanding_wus() >= 1
        key = "wu_received"
        self.state.update({key: self.state[key] + 1}, commit=False)
        self.wuar.verification(wuid, ok, commit=commit)
    
    def cancel_available_wus(self):
        self.logger.info("Cancelling remaining workunits")
        self.wuar.cancel_all_available()
    
    def get_number_outstanding_wus(self):
        return self.state["wu_submitted"] - self.state["wu_received"] \
                - self.state["wu_timedout"]

    def get_number_available_wus(self):
        return self.wuar.count_available()

    def test_outputfile_exists(self, filename):
        # Can't test
        return False
    
    def wait(self):
        # Ask the mediator to check for workunits of status Received,
        # and if there are any, to send WU result notifications to the
        # subscribed listeners.
        # If we get notification on new results reliably from the HTTP server,
        # we might not need this poll. But they probably won't be totally
        # reliable
        if not self.send_request(Request.GET_WU_RESULT):
            self.resubmit_timed_out_wus()
            time.sleep(1)
    
    def resubmit_one_wu(self, wu, commit=True, maxresubmit=None):
        """ Takes a Workunit instance and adds it to workunits table under
        a modified name.
        """
        wuid = wu.get_id()
        (name, task, identifier, attempt) = self.split_wuname(wuid)
        attempt = 2 if attempt is None else attempt + 1
        # Don't do "if not maxresubmit:" as 0 is legit value
        if maxresubmit is None:
            maxresubmit = self.params["maxresubmit"]
        if attempt > maxresubmit:
            self.logger.info("Not resubmitting workunit %s, failed %d times",
                             wuid, attempt - 1)
            self.wuar.commit(commit)
            return
        new_wuid = self.make_wuname(identifier, attempt)
        wu.set_id(new_wuid)
        self.logger.info("Resubmitting workunit %s as %s", wuid, new_wuid)
        self.submit_wu(wu, commit=commit)
        
    def resubmit_timed_out_wus(self):
        """ Check for any timed out workunits and resubmit them """
        # We don't store the lastcheck in state as we do *not* want to check
        # instantly when we start up - clients should get a chance to upload
        # results first
        now = time.time()
        if not hasattr(self, "last_timeout_check"):
            self.logger.debug("Setting last timeout check to %f", now)
            self.last_timeout_check = now
            return
        
        check_every = 60 # Check every xx seconds
        if self.last_timeout_check + check_every >= now:
            # self.logger.info("It's not time to check yet, now = %f", now)
            return
        self.last_timeout_check = now
        
        timeout = self.params["wutimeout"]
        delta = datetime.timedelta(seconds=timeout)
        cutoff = str(datetime.datetime.utcnow() - delta)
        # self.logger.debug("Doing timeout check, cutoff=%s, and setting last check to %f",
        #                   cutoff, now)
        results = self.wuar.query(eq={"status":1}, lt={"timeassigned": cutoff})
        if not results:
            # self.logger.debug("Found no timed-out workunits")
            pass
        for entry in results:
            self.cancel_wu(entry["wuid"], commit=False)
            self.resubmit_one_wu(Workunit(entry["wu"]), commit=True)

    def handle_error_result(self, message):
        """ Handle workunit with non-zero exit code
        
        If the result message indicates a failed command, log an error
        message, set the workunit to VERIFIED_ERROR in the DB, resubmit
        the work unit (but no more than once) and return True.
        If it indicates no error, return False. """
        if message.get_exitcode(0) == 0:
            return False
        self.log_failed_command_error(message, 0)
        key = "wu_failed"
        maxfailed = self.params["maxfailed"]
        if not self.state[key] < maxfailed:
            self.logger.error("Exceeded maximum number of failed "
                              "workunits, maxfailed=%d ", maxfailed)
            raise Exception("Too many failed work units")
        results = self.wuar.query(eq={"wuid":message.get_wu_id()})
        assert len(results) == 1 # There must be exactly 1 WU
        assert results[0]["status"] == wudb.WuStatus.RECEIVED_ERROR
        wu = workunit.Workunit(results[0]["wu"])
        self.state.update({key: self.state[key] + 1}, commit=False)
        self.verification(message.get_wu_id(), False, commit=False)
        self.resubmit_one_wu(wu, commit=True, maxresubmit=2)
        return True


class PolyselTask(ClientServerTask, HasStatistics, patterns.Observer):
    """ Finds a polynomial, uses client/server """
    @property
    def name(self):
        return "polyselect"
    @property
    def title(self):
        return "Polynomial Selection"
    @property
    def programs(self):
        return (cadoprograms.Polyselect2l,)
    @property
    def progparam_override(self):
        # admin and admax are special, which is a bit ugly: these parameters
        # to the Polyselect2l constructor are supplied by the task, but the
        # task has itself admin, admax parameters, which specify the total
        # size of the search range. Thus we don't include admin, admax here,
        # or PolyselTask would incorrectly warn about them not being used.
        return [[]]
    @property
    def paramnames(self):
        return self.join_params(super().paramnames, {
            "adrange": int, "admin": 0, "admax": int, "import": None,
            "I": int, "alim": int, "rlim": int})
    def update_lognorms(old_lognorm, new_lognorm):
        lognorm = [0, 0, 0, 0, 0]
        # print("update_lognorms: old_lognorm: %s" % old_lognorm)
        # print("update_lognorms: new_lognorm: %s" % new_lognorm)
        # New minimum. Don't use default value of 0 for minimum
        lognorm[1] = min(old_lognorm[1] or new_lognorm[1], new_lognorm[1])
        # New maximum
        lognorm[3] = max(old_lognorm[3], new_lognorm[3])
        # Rest is done by combine_stats(). [0::2] selects indices 0,2,4
        lognorm[0::2] = Statistics.combine_stats(old_lognorm[0::2],
                                                 new_lognorm[0::2])
        return lognorm
    
    # Stat: potential collisions=124.92 (2.25e+00/s)
    # Stat: raw lognorm (nr/min/av/max/std): 132/18.87/21.83/24.31/0.48
    # Stat: optimized lognorm (nr/min/av/max/std): 125/20.10/22.73/24.42/0.69
    # Stat: av. g0/adm2 ratio: 8.594e+04
    # Stat: tried 83 ad-value(s), found 132 polynomial(s), 125 below maxnorm
    # Stat: best logmu: 20.10 21.05 21.41 21.48 21.51 21.57 21.71 21.74 21.76 21.76
    # Stat: total phase took 55.47s
    # Stat: rootsieve took 54.54s
    @property
    def stat_conversions(self):
        return (
        (
            "stats_collisions",
            float,
            "0",
            Statistics.add_list,
            re.compile(r"# Stat: potential collisions=%s" % CAP_FP)
        ),
        (
            "stats_rawlognorm",
            (int, float, float, float, float),
            "0 0 0 0 0",
            update_lognorms,
            re.compile(r"# Stat: raw lognorm \(nr/min/av/max/std\): (\d+)/%s/%s/%s/%s" % ((CAP_FP,) * 4))
        ),
        (
            "stats_optlognorm",
            (int, float, float, float, float),
            "0 0 0 0 0",
            update_lognorms,
            re.compile(r"# Stat: optimized lognorm \(nr/min/av/max/std\): (\d+)/%s/%s/%s/%s" % ((CAP_FP,) * 4))
        ),
        (
            "stats_tries",
            int,
            "0 0 0",
            Statistics.add_list,
            re.compile(r"# Stat: tried (\d+) ad-value\(s\), found (\d+) polynomial\(s\), (\d+) below maxnorm")
        ),
        # Note for "best logmu" pattern: a regex like (%s )* does not work;
        # the number of the capture group is determined by the parentheses
        # in the regex string, so trying to repeat a group like this will
        # always capture to the *same* group, overwriting previous matches,
        # so that in the end, only the last match is in the capture group.
        (
            "stats_logmu",
            float,
            "",
            Statistics.smallest_n,
            re.compile(r"# Stat: best logmu:" + (" " + CAP_FP)*10)
        ),
        (
            "stats_total_time",
            float,
            "0",
            Statistics.add_list,
            re.compile(r"# Stat: total phase took %ss" % CAP_FP)
        ),
        (
            "stats_rootsieve_time",
            float,
            "0",
            Statistics.add_list,
            re.compile(r"# Stat: rootsieve took %ss" % CAP_FP)
        )
    )
    @property
    def stat_formats(self):
        return (
            ["potential collisions: {stats_collisions[0]:g}"],
            ["raw lognorm (nr/min/av/max/std): {stats_rawlognorm[0]:d}"] + 
                ["/{stats_rawlognorm[%d]:.3f}" % i for i in range(1, 5)],
            ["optimized lognorm (nr/min/av/max/std): {stats_optlognorm[0]:d}"] +
                ["/{stats_optlognorm[%d]:.3f}" % i for i in range(1, 5)],
            ["10 best logmu: {stats_logmu[0]}"] +
                [" {stats_logmu[%d]}" % i for i in range(1, 10)],
            ["Total time: {stats_total_time[0]:g}",
             ", rootsieve time: {stats_rootsieve_time[0]:g}"]
            )
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.state["adnext"] = \
            max(self.state.get("adnext", 0), self.params["admin"])
        self.bestpoly = None
        if "bestpoly" in self.state:
            self.bestpoly = Polynomials(self.state["bestpoly"].splitlines())
        self.did_import = False
        self.progparams[0].setdefault("area", 2.**(2*self.params["I"]-1) \
                * self.params["alim"])
        self.progparams[0].setdefault("Bf", float(self.params["alim"]))
        self.progparams[0].setdefault("Bg", float(self.params["rlim"]))
    
    def run(self):
        self.logger.info("Starting")
        self.logger.debug("%s.run(): Task state: %s", self.name, self.state)
        
        if not self.bestpoly is None:
            self.logger.info("Best polynomial previously found in %s has "
                             "Murphy_E = %g",
                             self.state["bestfile"], self.bestpoly.MurphyE)
        else:
            self.logger.info("No polynomial was previously found")
        
        if "import" in self.params and not self.did_import:
            self.process_polyfile(self.params["import"])
            if not self.bestpoly is None:
                self.write_poly_file()
            self.did_import = True
        
        if self.is_done():
            self.logger.info("Already finished - nothing to do")
            # If the poly file got lost somehow, write it again
            filename = self.get_state_filename("polyfilename")
            if filename is None or not filename.isfile():
                self.logger.warn("Polynomial file disappeared, writing again")
                self.write_poly_file()
            return True
        
        # Submit all the WUs we need to reach admax
        while self.need_more_wus():
            self.submit_one_wu()
        
        # Wait for all the WUs to finish
        while self.get_number_outstanding_wus() > 0:
            self.wait()
        
        if self.bestpoly is None:
            self.logger.error ("No polynomial found. Consider increasing the "
                               "search range bound admax, or maxnorm")
            return False
        self.logger.info("Finished, best polynomial from file %s has Murphy_E "
                         "= %g", self.state["bestfile"] , self.bestpoly.MurphyE)
        self.write_poly_file()
        return True
    
    def is_done(self):
        return not self.bestpoly is None and not self.need_more_wus() and \
            self.get_number_outstanding_wus() == 0
    
    def updateObserver(self, message):
        identifier = self.filter_notification(message)
        if not identifier:
            # This notification was not for me
            return False
        if self.handle_error_result(message):
            return True
        (filename, ) = message.get_output_files()
        self.process_polyfile(filename, commit=False)
        self.parse_stats(filename, commit=False)
        # Always mark ok to avoid warning messages about WUs that did not
        # find a poly
        self.verification(message.get_wu_id(), True, commit=True)
        return True
    
    def process_polyfile(self, filename, commit=True):
        poly = self.parse_poly(filename)
        if not poly is None:
            self.bestpoly = poly
            update = {"bestpoly": str(poly), "bestfile": filename}
            self.state.update(update, commit=commit)
    
    def read_log_warning(self, filename):
        """ Read lines from file. If a "# WARNING" line occurs, log it.
        """
        re_warning = re.compile("# WARNING")
        with open(filename, "r") as inputfile:
            for line in inputfile:
                if re_warning.match(line):
                    self.logger.warn("File %s contains: %s",
                                     filename, line.strip())
                yield line

    def parse_poly(self, filename):
        poly = None
        try:
            poly = Polynomials(self.read_log_warning(filename))
        except (OSError, IOError) as e:
            if e.errno == 2: # No such file or directory
                self.logger.error("File '%s' does not exist", filename)
                return None
            else:
                raise
        except PolynomialParseException as e:
            if str(e) != "No polynomials found":
                self.logger.warn("Invalid polyselect file '%s': %s",
                                  filename, e)
                return None
        except UnicodeDecodeError as e:
            self.logger.error("Error reading '%s' (corrupted?): %s", filename, e)
            return None
        
        if not poly:
            self.logger.info('No polynomial found in %s', filename)
            return None
        if not poly.MurphyE:
            self.logger.warn("Polynomial in file %s has no Murphy E value",
                             filename)
        if self.bestpoly is None or poly.MurphyE > self.bestpoly.MurphyE:
            self.logger.info("New best polynomial from file %s:"
                             " Murphy E = %g" % (filename, poly.MurphyE))
            self.logger.debug("New best polynomial is:\n%s", poly)
            return poly
        else:
            self.logger.info("Best polynomial from file %s with E=%g is "
                             "no better than current best with E=%g",
                             filename, poly.MurphyE, self.bestpoly.MurphyE)
        return None
    
    def write_poly_file(self):
        filename = self.workdir.make_filename("poly")
        self.bestpoly.create_file(filename)
        self.state["polyfilename"] = filename.get_wdir_relative()
    
    def get_poly(self):
        if not "bestpoly" in self.state:
            return None
        return Polynomials(self.state["bestpoly"].splitlines())
    
    def get_poly_filename(self):
        return self.get_state_filename("polyfilename")

    def need_more_wus(self):
        return self.state["adnext"] < self.params["admax"]
    
    def submit_one_wu(self):
        adstart = self.state["adnext"]
        adend = adstart + self.params["adrange"]
        adend = min(adend, self.params["admax"])
        outputfile = self.workdir.make_filename("%d-%d" % (adstart, adend))
        if self.test_outputfile_exists(outputfile):
            self.logger.info("%s already exists, won't generate again",
                             outputfile)
        else:
            # Remove admin and admax from the parameter-file-supplied
            # parameters as those would conflict with the computed values
            self.progparams[0].pop("admin", None)
            self.progparams[0].pop("admax", None)
            p = cadoprograms.Polyselect2l(admin=adstart, admax=adend,
                                          stdout=str(outputfile),
                                          **self.progparams[0])
            self.submit_command(p, "%d-%d" % (adstart, adend), commit=False)
        self.state.update({"adnext": adend}, commit=True)

    def get_total_cpu_or_real_time(self, is_cpu):
        """ Return number of seconds of cpu time spent by polyselect2l """
        return float(self.state.get("stats_total_time", 0.)) if is_cpu else 0.

class Polysel1Task(ClientServerTask, HasStatistics, patterns.Observer):
    """ Finds a number of size-optimized polynomial, uses client/server """
    @property
    def name(self):
        return "polyselect1"
    @property
    def title(self):
        return "Polynomial Selection (size optimized)"
    @property
    def programs(self):
        return (cadoprograms.Polyselect2l,)
    @property
    def progparam_override(self):
        # admin and admax are special, which is a bit ugly: these parameters
        # to the Polyselect2l constructor are supplied by the task, but the
        # task has itself admin, admax parameters, which specify the total
        # size of the search range. Thus we don't include admin, admax here,
        # or PolyselTask would incorrectly warn about them not being used.
        return [[]]
    @property
    def paramnames(self):
        return self.join_params(super().paramnames, {
            "adrange": int, "admin": 0, "admax": int, "import": None,
            "I": int, "alim": int, "rlim": int, "nrkeep": 20})
    @staticmethod
    def update_lognorms(old_lognorm, new_lognorm):
        lognorm = [0, 0, 0, 0, 0]
        # print("update_lognorms: old_lognorm: %s" % old_lognorm)
        # print("update_lognorms: new_lognorm: %s" % new_lognorm)
        # New minimum. Don't use default value of 0 for minimum
        lognorm[1] = min(old_lognorm[1] or new_lognorm[1], new_lognorm[1])
        # New maximum
        lognorm[3] = max(old_lognorm[3], new_lognorm[3])
        # Rest is done by combine_stats(). [0::2] selects indices 0,2,4
        lognorm[0::2] = Statistics.combine_stats(old_lognorm[0::2],
                                                 new_lognorm[0::2])
        return lognorm
    
    # Stat: potential collisions=124.92 (2.25e+00/s)
    # Stat: raw lognorm (nr/min/av/max/std): 132/18.87/21.83/24.31/0.48
    # Stat: optimized lognorm (nr/min/av/max/std): 125/20.10/22.73/24.42/0.69
    # Stat: av. g0/adm2 ratio: 8.594e+04
    # Stat: tried 83 ad-value(s), found 132 polynomial(s), 125 below maxnorm
    # Stat: best raw logmu: 20.25 20.88
    # Stat: best opt logmu: 19.68 20.25
    # Stat: best logmu: 19.68 20.25
    # Stat: total phase took 55.47s
    @property
    def stat_conversions(self):
        return (
        (
            "stats_collisions",
            float,
            "0",
            Statistics.add_list,
            re.compile(r"# Stat: potential collisions=%s" % CAP_FP)
        ),
        (
            "stats_rawlognorm",
            (int, float, float, float, float),
            "0 0 0 0 0",
            self.update_lognorms,
            re.compile(r"# Stat: raw lognorm \(nr/min/av/max/std\): (\d+)/%s/%s/%s/%s" % ((CAP_FP,) * 4))
        ),
        (
            "stats_optlognorm",
            (int, float, float, float, float),
            "0 0 0 0 0",
            self.update_lognorms,
            re.compile(r"# Stat: optimized lognorm \(nr/min/av/max/std\): (\d+)/%s/%s/%s/%s" % ((CAP_FP,) * 4))
        ),
        (
            "stats_tries",
            int,
            "0 0 0",
            Statistics.add_list,
            re.compile(r"# Stat: tried (\d+) ad-value\(s\), found (\d+) polynomial\(s\), (\d+) below maxnorm")
        ),
        # Note for "best raw logmu" pattern: a regex like (%s )* does not work;
        # the number of the capture group is determined by the parentheses
        # in the regex string, so trying to repeat a group like this will
        # always capture to the *same* group, overwriting previous matches,
        # so that in the end, only the last match is in the capture group.
        (
            "stats_raw_logmu",
            float,
            "",
            Statistics.smallest_n,
            re.compile(r"# Stat: best raw logmu:" + ("(?: " + CAP_FP + ")?")*10)
        ),
        (
            "stats_opt_logmu",
            float,
            "",
            Statistics.smallest_n,
            re.compile(r"# Stat: best opt logmu:" + ("(?: " + CAP_FP + ")?")*10)
        ),
        (
            "stats_total_time",
            float,
            "0",
            Statistics.add_list,
            re.compile(r"# Stat: total phase took %ss" % CAP_FP)
        ),
    )
    @property
    def stat_formats(self):
        return (
            ["potential collisions: {stats_collisions[0]:g}"],
            ["raw lognorm (nr/min/av/max/std): {stats_rawlognorm[0]:d}"] + 
                ["/{stats_rawlognorm[%d]:.3f}" % i for i in range(1, 5)],
            ["optimized lognorm (nr/min/av/max/std): {stats_optlognorm[0]:d}"] +
                ["/{stats_optlognorm[%d]:.3f}" % i for i in range(1, 5)],
            ["10 best raw logmu:"] +
                [" {stats_raw_logmu[%d]}" % i for i in range(10)],
            ["10 best opt logmu:"] +
                [" {stats_opt_logmu[%d]}" % i for i in range(10)],
            ["Total time: {stats_total_time[0]:g}",
             ", rootsieve time: {stats_rootsieve_time[0]:g}"]
            )
    
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        assert self.params["nrkeep"] > 0
        self.state["adnext"] = \
            max(self.state.get("adnext", 0), self.params["admin"])
        self.did_import = False
        self.progparams[0].setdefault("area", 2.**(2*self.params["I"]-1) \
                * self.params["alim"])
        self.progparams[0].setdefault("Bf", float(self.params["alim"]))
        self.progparams[0].setdefault("Bg", float(self.params["rlim"]))
        # Remove admin and admax from the parameter-file-supplied program
        # parameters as those would conflict with the computed values
        self.progparams[0].pop("admin", None)
        self.progparams[0].pop("admax", None)

        tablename = self.make_tablename("bestpolynomials")
        self.best_polynomials = self.make_db_dict(
                tablename, connection=self.db_connection)
        self._check_best_polynomials()

        self.poly_heap = []
        self.import_existing_polynomials()
        self._check_best_polynomials()
        self._compare_heap_db()
    
    def _check_best_polynomials(self):
        # Check that the keys form a sequence of consecutive non-negative
        # integers
        oldkeys = list(self.best_polynomials.keys())
        oldkeys.sort(key=int)
        assert oldkeys == list(map(str, range(len(self.best_polynomials))))

    def _compare_heap_db(self):
        """ Compare that the polynomials in the heap and in the DB agree
        
        They must contain an equal number of entries, and each polynomial
        stored in the heap must be at the specified index in the DB.
        """
        assert len(self.poly_heap) == len(self.best_polynomials)
        for lognorm, (key, poly) in self.poly_heap:
            assert self.best_polynomials[key] == str(poly)

    def import_existing_polynomials(self):
        debug = False
        oldkeys = list(self.best_polynomials.keys())
        oldkeys.sort(key=int) # Sort by numerical value
        for oldkey in oldkeys:
            if debug:
                print("Adding old polynomial at DB index %s: %s" %
                      (oldkey, self.best_polynomials[oldkey]))
            poly = Polynomials(self.best_polynomials[oldkey].splitlines())
            if not poly.lognorm:
                self.logger.error("Polynomial at DB index %s has no lognorm", oldkey)
                continue
            newkey = self._add_poly_heap(poly)
            if newkey is None:
                # Heap is full, and the poly was worse than the worst one on
                # the heap. Thus it did not get added and must be removed from
                # the DB
                if debug:
                    print("Deleting polynomial lognorm=%f, key=%s" % 
                          (poly.lognorm, oldkey))
                del(self.best_polynomials[oldkey])
            elif newkey != oldkey:
                # Heap is full, worst one in heap (with key=newkey) was
                # overwritten and its DB entry gets replaced with poly from
                # key=oldkey
                if debug:
                    print("Overwriting poly lognorm=%f, key=%s with poly "
                          "lognorm=%f, key=%s" %
                          (self.poly_heap[0][0], newkey, poly, oldkey))
                self.best_polynomials.clear(oldkey, commit=False)
                self.best_polynomials.update({newkey: poly}, commit=True)
            else:
                # Last case newkey == oldkey: nothing to do
                if debug:
                    print("Adding lognorm=%f, key=%s" % (poly.lognorm, oldkey))

    def run(self):
        self.logger.info("Starting")
        self.logger.debug("%s.run(): Task state: %s", self.name, self.state)
        
        worstmsg = ", worst lognorm %f" % -self.poly_heap[0][0] \
                if self.poly_heap else ""
        self.logger.info("%d polynomials in queue from previous run%s", 
                         len(self.poly_heap), worstmsg)
        
        if "import" in self.params and not self.did_import:
            self.process_polyfile(self.params["import"])
            self.did_import = True
        
        if self.is_done():
            self.logger.info("Already finished - nothing to do")
            return True
        
        # Submit all the WUs we need to reach admax
        while self.need_more_wus():
            self.submit_one_wu()
        
        # Wait for all the WUs to finish
        while self.get_number_outstanding_wus() > 0:
            self.wait()
        
        self._compare_heap_db()
        self.logger.info("Finished")
        return True
    
    def is_done(self):
        return not self.need_more_wus() and \
            self.get_number_outstanding_wus() == 0
    
    def updateObserver(self, message):
        identifier = self.filter_notification(message)
        if not identifier:
            # This notification was not for me
            return False
        if self.handle_error_result(message):
            return True
        (filename, ) = message.get_output_files()
        self.process_polyfile(filename, commit=False)
        self.parse_stats(filename, commit=False)
        # Always mark ok to avoid warning messages about WUs that did not
        # find a poly
        self.verification(message.get_wu_id(), True, commit=True)
        return True
    
    @staticmethod
    def read_blocks(input):
        """ Return blocks of consecutive non-empty lines from input
        
        Whitespace is stripped; a line containing only whitespace is
        considered empty. An empty block is never returned.
        
        >>> list(Polysel1Task.read_blocks(['', 'a', 'b', '', 'c', '', '', 'd', 'e', '']))
        [['a', 'b'], ['c'], ['d', 'e']]
        """
        block = []
        for line in input:
            line = line.strip()
            if line:
                block.append(line)
            else:
                if block:
                    yield block
                block = []
        if block:
            yield block

    def process_polyfile(self, filename, commit=True):
        """ Read all size-optimized polynomials in a file and add them to the
        DB and priority queue if worthwhile.
        
        Different polynomials must be separated by a blank line.
        """
        try:
            polyfile = self.read_log_warning(filename)
        except (OSError, IOError) as e:
            if e.errno == 2: # No such file or directory
                self.logger.error("File '%s' does not exist", filename)
                return None
            else:
                raise
        totalparsed, totaladded = 0, 0
        for block in self.read_blocks(polyfile):
            parsed, added = self.parse_and_add_poly(block, filename)
            totalparsed += parsed
            totaladded += added
        have = len(self.poly_heap)
        nrkeep = self.params["nrkeep"]
        fullmsg = ("%d/%d" % (have, nrkeep)) if have < nrkeep else "%d" % nrkeep
        self.logger.info("Parsed %d polynomials, added %d to priority queue (has %s)",
                         totalparsed, totaladded, fullmsg)
        if totaladded:
            self.logger.info("Worst polynomial in queue now has lognorm %f",
                             -self.poly_heap[0][0])
                                     
    
    def read_log_warning(self, filename):
        """ Read lines from file. If a "# WARNING" line occurs, log it.
        """
        re_warning = re.compile("# WARNING")
        with open(filename, "r") as inputfile:
            for line in inputfile:
                if re_warning.match(line):
                    self.logger.warn("File %s contains: %s",
                                     filename, line.strip())
                yield line

    def parse_and_add_poly(self, text, filename):
        """ Parse a polynomial from an iterable of lines and add it to the
        priority queue and DB. Return a two-element list with the number of
        polynomials parsed and added, i.e., (0,0) or (1,0) or (1,1).
        """
        poly = self.parse_poly(text, filename)
        if poly is None:
            return (0, 0)
        if not poly.lognorm:
            self.logger.warn("Polynomial in file %s has no lognorm, skipping it",
                             filename)
            return (0, 0)
        if self._add_poly_heap_db(poly):
            return (1, 1)
        else:
            return (1, 0)

    def _add_poly_heap_db(self, poly):
        """ Add a polynomial to the heap and DB, if it's good enough.
        
        Returns True if the poly was added, False if not. """
        key = self._add_poly_heap(poly)
        if key is None:
            return False
        self.best_polynomials[key] = str(poly)
        return True

    def _add_poly_heap(self, poly):
        """ Add a polynomial to the heap
        
        If the heap is full (nrkeep), the worst polynomial (i.e., with the
        largest lognorm) is replaced if the new one is better.
        Returns the key (as a str) under which the polynomial was added,
        or None if it was not added.
        """
        assert len(self.poly_heap) <= self.params["nrkeep"]
        debug = False

        # Find DB index under which to store this new poly. If the heap
        # is not full, use the next bigger index.
        key = len(self.poly_heap)
        # Is the heap full?
        if key == self.params["nrkeep"]:
            # Should we store this poly at all, i.e., is it better than
            # the worst one in the heap?
            worstnorm = -self.poly_heap[0][0]
            if worstnorm <= poly.lognorm:
                if debug:
                    self.logger.debug("_add_poly_heap(): new poly lognorm %f, "
                          "worst in heap has %f. Not adding",
                          poly.lognorm, worstnorm)
                return None
            # Pop the worst poly from heap and re-use its DB index
            key = heapq.heappop(self.poly_heap)[1][0]
            if debug:
                self.logger.debug("_add_poly_heap(): new poly lognorm %f, "
                    "worst in heap has %f. Replacing DB index %s",
                     poly.lognorm, worstnorm, key)
        else:
            # Heap was not full
            if debug:
                self.logger.debug("_add_poly_heap(): heap was not full, adding "
                    "poly with lognorm %f at DB index %s", poly.lognorm, key)

        # The DB requires the key to be a string. In order to have
        # identical data in DB and heap, we store key as str everywhere.
        key = str(key)

        # Python heapq stores a minheap, so in order to have the worst
        # polynomial (with largest norm) easily accessible, we use
        # -lognorm as the heap key
        new_entry = (-poly.lognorm, (key, poly))
        heapq.heappush(self.poly_heap, new_entry)
        return key

    def parse_poly(self, text, filename):
        poly = None
        try:
            poly = Polynomials(text)
        except PolynomialParseException as e:
            if str(e) != "No polynomials found":
                self.logger.warn("Invalid polyselect file '%s': %s",
                                  filename, e)
                return None
        except UnicodeDecodeError as e:
            self.logger.error("Error reading '%s' (corrupted?): %s", filename, e)
            return None
        
        if not poly:
            return None
        return poly
    
    def get_raw_polynomials(self):
        # Extract polynomials from heap and return as list
        return [entry[1][1] for entry in self.poly_heap]

    def get_poly_rank(self, search_poly):
        """ Return how many polynomnials with lognorm less than the lognorm
        of the size-optimized version of search_poly there are in the
        priority queue.
        
        The size-optimized version of search_poly is identified by comparing
        the leading coefficients of both polynomials.
        """
        df = search_poly.polyf.degree
        dg = search_poly.polyg.degree
        # Search for the raw polynomial pair by comparing the leading
        # coefficients of both polynomials
        found = None
        for (index, (lognorm, (key, poly))) in enumerate(self.poly_heap):
            if poly.polyf.degree == df and poly.polyg.degree == dg and \
                   poly.polyf[df] == search_poly.polyf[df] and \
                   poly.polyg[dg] == search_poly.polyg[dg]:
               if not found is None:
                   self.logger.error("Found more than one match for %s", search_poly)
               else:
                   found = index
        if found is None:
            self.logger.error("Could not find polynomial rank for %s", search_poly)
            return None
        # print("search_poly: %s" % search_poly)
        # print("Poly found in heap: %s" % self.poly_heap[found][1][1])
        search_lognorm = -self.poly_heap[found][0]
        rank = 0
        for (lognorm, (key, poly)) in self.poly_heap:
            if -lognorm < search_lognorm:
                rank += 1
        return rank

    def need_more_wus(self):
        return self.state["adnext"] < self.params["admax"]
    
    def submit_one_wu(self):
        adstart = self.state["adnext"]
        adend = adstart + self.params["adrange"]
        adend = min(adend, self.params["admax"])
        outputfile = self.workdir.make_filename("%d-%d" % (adstart, adend))
        if self.test_outputfile_exists(outputfile):
            self.logger.info("%s already exists, won't generate again",
                             outputfile)
        else:
            p = cadoprograms.Polyselect2l(admin=adstart, admax=adend,
                                          stdout=str(outputfile),
                                          sizeonly=True,
                                          **self.progparams[0])
            self.submit_command(p, "%d-%d" % (adstart, adend), commit=False)
        self.state.update({"adnext": adend}, commit=True)

    def get_total_cpu_or_real_time(self, is_cpu):
        """ Return number of seconds of cpu time spent by polyselect2l """
        return float(self.state.get("stats_total_time", 0.)) if is_cpu else 0.


class Polysel2Task(ClientServerTask, HasStatistics, patterns.Observer):
    """ Finds a polynomial, uses client/server """
    @property
    def name(self):
        return "polyselect2"
    @property
    def title(self):
        return "Polynomial Selection (root optimized)"
    @property
    def programs(self):
        return (cadoprograms.Polyselect2l,)
    @property
    def progparam_override(self):
        return [[]]
    @property
    def paramnames(self):
        return self.join_params(super().paramnames, {
            "import": None, "I": int, "alim": int, "rlim": int, "batch": 5})
    @property
    def stat_conversions(self):
        return (
        (
            "stats_total_time",
            float,
            "0",
            Statistics.add_list,
            re.compile(r"# Stat: total phase took %ss" % CAP_FP)
        ),
        (
            "stats_rootsieve_time",
            float,
            "0",
            Statistics.add_list,
            re.compile(r"# Stat: rootsieve took %ss" % CAP_FP)
        )
    )
    @property
    def stat_formats(self):
        return (
            ["Total time: {stats_total_time[0]:g}",
             ", rootsieve time: {stats_rootsieve_time[0]:g}"],
            )

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.bestpoly = None
        if "bestpoly" in self.state:
            self.bestpoly = Polynomials(self.state["bestpoly"].splitlines())
        self.did_import = False
        self.state.setdefault("nr_poly_submitted", 0)
        self.progparams[0].setdefault("area", 2.**(2*self.params["I"]-1) \
                * self.params["alim"])
        self.progparams[0].setdefault("Bf", float(self.params["alim"]))
        self.progparams[0].setdefault("Bg", float(self.params["rlim"]))
        self.poly_to_submit = None

    def run(self):
        self.logger.info("Starting")
        self.logger.debug("%s.run(): Task state: %s", self.name, self.state)
        
        if self.bestpoly is None:
            self.logger.info("No polynomial was previously found")
        else:
            self.logger.info("Best polynomial previously found in %s has "
                             "Murphy_E = %g",
                             self.state["bestfile"], self.bestpoly.MurphyE)
        
        if "import" in self.params and not self.did_import:
            self.process_polyfile(self.params["import"])
            if not self.bestpoly is None:
                self.write_poly_file()
            self.did_import = True
        
        # Get the list of polynomials to submit
        self.poly_to_submit = self.send_request(Request.GET_RAW_POLYNOMIALS)
        
        if self.is_done():
            self.logger.info("Already finished - nothing to do")
            rank = self.send_request(Request.GET_POLY_RANK, self.bestpoly)
            if not rank is None:
                self.logger.info("Best overall polynomial was %d-th in list "
                                 "after size optimization", rank)
            # If the poly file got lost somehow, write it again
            filename = self.get_state_filename("polyfilename")
            if filename is None or not filename.isfile():
                self.logger.warn("Polynomial file disappeared, writing again")
                self.write_poly_file()
            return True
        
        # Submit all the WUs we need
        while self.need_more_wus():
            self.submit_one_wu()
        
        # Wait for all the WUs to finish
        while self.get_number_outstanding_wus() > 0:
            self.wait()
        
        if self.bestpoly is None:
            self.logger.error ("No polynomial found. Consider increasing the "
                               "search range bound admax, or maxnorm")
            return False
        self.logger.info("Finished, best polynomial from file %s has Murphy_E "
                         "= %g", self.state["bestfile"] , self.bestpoly.MurphyE)
        rank = self.send_request(Request.GET_POLY_RANK, self.bestpoly)
        if not rank is None:
            self.logger.info("Best overall polynomial was %d-th in list "
                             "after size optimization", rank)
        self.write_poly_file()
        return True
    
    def is_done(self):
        return not self.bestpoly is None and not self.need_more_wus() and \
            self.get_number_outstanding_wus() == 0
    
    def need_more_wus(self):
        return self.state["nr_poly_submitted"] < len(self.poly_to_submit)
    
    def updateObserver(self, message):
        identifier = self.filter_notification(message)
        if not identifier:
            # This notification was not for me
            return False
        if self.handle_error_result(message):
            return False
        (filename, ) = message.get_output_files()
        self.process_polyfile(filename, commit=False)
        self.parse_stats(filename, commit=False)
        # Always mark ok to avoid warning messages about WUs that did not
        # find a poly
        # FIXME: wrong, we should always get an optimized poly for a raw one
        self.verification(message.get_wu_id(), True, commit=True)
        return True
    
    def process_polyfile(self, filename, commit=True):
        poly = self.parse_poly(filename)
        if not poly is None:
            self.bestpoly = poly
            update = {"bestpoly": str(poly), "bestfile": filename}
            self.state.update(update, commit=commit)
    
    def read_log_warning(self, filename):
        """ Read lines from file. If a "# WARNING" line occurs, log it.
        """
        re_warning = re.compile("# WARNING")
        with open(filename, "r") as inputfile:
            for line in inputfile:
                if re_warning.match(line):
                    self.logger.warn("File %s contains: %s",
                                     filename, line.strip())
                yield line

    def parse_poly(self, filename):
        poly = None
        try:
            poly = Polynomials(self.read_log_warning(filename))
        except (OSError, IOError) as e:
            if e.errno == 2: # No such file or directory
                self.logger.error("File '%s' does not exist", filename)
                return None
            else:
                raise
        except PolynomialParseException as e:
            if str(e) != "No polynomials found":
                self.logger.warn("Invalid polyselect file '%s': %s",
                                  filename, e)
                return None
        except UnicodeDecodeError as e:
            self.logger.error("Error reading '%s' (corrupted?): %s", filename, e)
            return None
        
        if not poly:
            self.logger.info('No polynomial found in %s', filename)
            return None
        if not poly.MurphyE:
            self.logger.warn("Polynomial in file %s has no Murphy E value",
                             filename)
        if self.bestpoly is None or poly.MurphyE > self.bestpoly.MurphyE:
            self.logger.info("New best polynomial from file %s:"
                             " Murphy E = %g" % (filename, poly.MurphyE))
            self.logger.debug("New best polynomial is:\n%s", poly)
            return poly
        else:
            self.logger.info("Best polynomial from file %s with E=%g is "
                             "no better than current best with E=%g",
                             filename, poly.MurphyE, self.bestpoly.MurphyE)
        return None
    
    def write_poly_file(self):
        filename = self.workdir.make_filename("poly")
        self.bestpoly.create_file(filename)
        self.state["polyfilename"] = filename.get_wdir_relative()
    
    def get_poly(self):
        if not "bestpoly" in self.state:
            return None
        return Polynomials(self.state["bestpoly"].splitlines())
    
    def get_poly_filename(self):
        return self.get_state_filename("polyfilename")

    def submit_one_wu(self):
        assert self.need_more_wus()
        to_submit = len(self.poly_to_submit)
        nr = self.state["nr_poly_submitted"]
        inputfilename = self.workdir.make_filename("raw_%d" % nr)
        # Write one raw polynomial to inputfile
        batchsize = min(to_submit - nr, self.params["batch"])
        with inputfilename.open("w") as inputfile:
            for i in range(batchsize):
                inputfile.write(str(self.poly_to_submit[nr + i]))
                inputfile.write("\n")
        outputfile = self.workdir.make_filename("opt_%d" % nr)
        if self.test_outputfile_exists(outputfile):
            self.logger.info("%s already exists, won't generate again",
                             outputfile)
        else:
            p = cadoprograms.Polyselect2l(rootsieve=str(inputfilename),
                                          stdout=str(outputfile),
                                          **self.progparams[0])
            self.submit_command(p, "%d" % nr, commit=False)
        self.state.update({"nr_poly_submitted": nr + batchsize}, commit=True)

    def get_total_cpu_or_real_time(self, is_cpu):
        """ Return number of seconds of cpu time spent by polyselect2l """
        return float(self.state.get("stats_total_time", 0.)) if is_cpu else 0.

class FactorBaseTask(Task):
    """ Generates the factor base for the polynomial(s) """
    @property
    def name(self):
        return "factorbase"
    @property
    def title(self):
        return "Generate Factor Base"
    @property
    def programs(self):
        return (cadoprograms.MakeFB,)
    @property
    def progparam_override(self):
        return [["poly", "out"]]
    @property
    def paramnames(self):
        return self.join_params(super().paramnames, {"gzip": True, "I": int})

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        # Invariant: if we have a result (in self.state["outputfile"]) then we
        # must also have a polynomial (in self.state["poly"])
        if "outputfile" in self.state:
            assert "poly" in self.state
            # The target file must correspond to the polynomial "poly"
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.name, self.state)
        
        # Get best polynomial found by polyselect
        poly = self.send_request(Request.GET_POLYNOMIAL)
        if not poly:
            raise Exception("FactorBaseTask(): no polynomial "
                            "received from PolyselTask")
        
        # Check if we have already computed the target file for this polynomial
        if "poly" in self.state:
            prevpoly = Polynomials(self.state["poly"].splitlines())
            if poly != prevpoly:
                if "outputfile" in self.state:
                    self.logger.info("Received different polynomial, "
                                     "discarding old one")
                    del(self.state["outputfile"])
                self.state["poly"] = str(poly)
        else:
            self.state["poly"] = str(poly)
        
        if not "outputfile" in self.state:
            self.logger.info("Starting")
            polyfilename = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
            
            # Make file name for factor base/free relations file
            # We use .gzip by default, unless set to no in parameters
            use_gz = ".gz" if self.params["gzip"] else ""
            outputfilename = self.workdir.make_filename("roots" + use_gz)

            # Run command to generate factor base/free relations file
            self.progparams[0].setdefault("maxbits", self.params["I"] - 1)
            (stdoutpath, stderrpath) = \
                    self.make_std_paths(cadoprograms.MakeFB.name)
            p = cadoprograms.MakeFB(poly=polyfilename,
                                    out=str(outputfilename),
                                    stdout=str(stdoutpath),
                                    stderr=str(stderrpath),
                                    **self.progparams[0])
            message = self.submit_command(p, "", log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            
            self.state["outputfile"] = outputfilename.get_wdir_relative()
            self.logger.info("Finished")

        self.check_files_exist([self.get_filename()], "output",
                               shouldexist=True)
        return True
    
    def get_filename(self):
        return self.get_state_filename("outputfile")


class FreeRelTask(Task):
    """ Generates free relations for the polynomial(s) """
    @property
    def name(self):
        return "freerel"
    @property
    def title(self):
        return "Generate Free Relations"
    @property
    def programs(self):
        return (cadoprograms.FreeRel,)
    @property
    def progparam_override(self):
        return [["poly", "renumber", "badideals", "out"]]
    @property
    def paramnames(self):
        return self.join_params(super().paramnames,
                {"dlp": False, "gzip": True, "addfullcol": None})
    wanted_regex = {
        'nfree': (r'# Free relations: (\d+)', int),
        'nprimes': (r'Renumbering struct: nprimes=(\d+)', int)
    }
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        # Invariant: if we have a result (in self.state["freerelfilename"])
        # then we must also have a polynomial (in self.state["poly"])
        if "freerelfilename" in self.state:
            assert "poly" in self.state
            # The target file must correspond to the polynomial "poly"
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.name, self.state)
        
        # Get best polynomial found by polyselect
        poly = self.send_request(Request.GET_POLYNOMIAL)
        if not poly:
            raise Exception("FreerelTask(): no polynomial "
                            "received from PolyselTask")
        
        # Check if we have already computed the target file for this polynomial
        if "poly" in self.state:
            prevpoly = Polynomials(self.state["poly"].splitlines())
            if poly != prevpoly:
                if "freerelfilename" in self.state:
                    self.logger.info("Received different polynomial, "
                                     "discarding old one")
                    del(self.state["freerelfilename"])
                    del(self.state["renumberfilename"])
                self.state["poly"] = str(poly)
        else:
            self.state["poly"] = str(poly)
        
        if not "freerelfilename" in self.state:
            self.logger.info("Starting")
            # Write polynomial to a file
            polyfilename = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
            
            # Make file name for factor base/free relations file
            # We use .gzip by default, unless set to no in parameters
            use_gz = ".gz" if self.params["gzip"] else ""
            freerelfilename = self.workdir.make_filename("freerel" + use_gz)
            renumberfilename = self.workdir.make_filename("renumber" + use_gz)
            (stdoutpath, stderrpath) = \
                    self.make_std_paths(cadoprograms.FreeRel.name)
            if self.params["dlp"]:
                badidealfilename = self.send_request(Request.GET_BADIDEAL_FILENAME)
                # default for dlp is addfullcol
                if not "addfullcol" in self.progparams[0]:
                    self.progparams[0]["addfullcol"] = True;
                p = cadoprograms.FreeRel(poly=polyfilename,
                                         renumber=renumberfilename,
                                         badideals=badidealfilename,
                                         out=str(freerelfilename),
                                         stdout=str(stdoutpath),
                                         stderr=str(stderrpath),
                                         **self.progparams[0])
            else:
                # Run command to generate factor base/free relations file
                p = cadoprograms.FreeRel(poly=polyfilename,
                                         renumber=renumberfilename,
                                         out=str(freerelfilename),
                                         stdout=str(stdoutpath),
                                         stderr=str(stderrpath),
                                         **self.progparams[0])
            message = self.submit_command(p, "", log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            stderr = message.read_stderr(0).decode("utf-8")
            update = self.parse_file(stderr.splitlines())
            update["freerelfilename"] = freerelfilename.get_wdir_relative()
            update["renumberfilename"] = renumberfilename.get_wdir_relative()
            self.state.update(update) 
            self.logger.info("Found %d free relations" % self.state["nfree"])
            self.logger.info("Finished")

        self.check_files_exist([self.get_freerel_filename(),
                                self.get_renumber_filename()], "output",
                               shouldexist=True)
        return True

    def parse_file(self, text):
        found = {}
        for line in text:
            for (key, (regex, datatype)) in self.wanted_regex.items():
                match = re.match(regex, line)
                if match:
                    if key in found:
                        raise Exception("Received two values for %s" % key)
                    found[key] = datatype(match.group(1))
        
        for key in self.wanted_regex:
            if not key in found:
                raise Exception("Received no value for %s" % key)
        return found
    
    def get_freerel_filename(self):
        return self.get_state_filename("freerelfilename")
    
    def get_renumber_filename(self):
        return self.get_state_filename("renumberfilename")
    
    def get_nrels(self):
        return self.state.get("nfree", None)
    
    def get_nprimes(self):
        return self.state.get("nprimes", None)


class SievingTask(ClientServerTask, FilesCreator, HasStatistics,
                  patterns.Observer):
    """ Does the sieving, uses client/server """
    @property
    def name(self):
        return "sieving"
    @property
    def title(self):
        return "Lattice Sieving"
    @property
    def programs(self):
        return (cadoprograms.Las,)
    @property
    def progparam_override(self):
        return [["q0", "q1", "poly", "factorbase", "out", "stats_stderr"]]
    @property
    def paramnames(self):
        return self.join_params(super().paramnames, {
            "qmin": 0, "qrange": int, "rels_wanted": 0, "alim": int, 
            "import": None, "gzip": True})
    @property
    def stat_conversions(self):
        # Average J=1017 for 168 special-q's, max bucket fill 0.737035
        # Total cpu time 7.0s [precise timings available only for mono-thread]
        # Total 26198 reports [0.000267s/r, 155.9r/sq]
        return (
        (
            "stats_avg_J",
            (float, int),
            "0 0",
            Statistics.zip_combine_mean,
            re.compile(r"# Average J=%s for (\d+) special-q's" % CAP_FP)
        ),
        (
            "stats_max_bucket_fill",
            float,
            "0",
            max,
            re.compile(r"#.*max bucket fill %s" % CAP_FP)
        ),
        (
            "stats_total_cpu_time",
            float,
            "0",
            Statistics.add_list,
            re.compile(r"# Total cpu time %ss" % CAP_FP)
        ),
        (
            "stats_total_time",
            (float, ),
            "0",
            Statistics.add_list,
            re.compile(r"# Total time %ss" % CAP_FP)
        )
    )
    @property
    def stat_formats(self):
        return (
            ["Average J: {stats_avg_J[0]:g} for {stats_avg_J[1]:d} special-q",
                ", max bucket fill: {stats_max_bucket_fill[0]:g}"],
            ["Total CPU time: {stats_total_cpu_time[0]:g}s"],
            ["Total time: {stats_total_time[0]:g}s"],
        )
    # We seek to this many bytes before the EOF to look for the
    # "Total xxx reports" message
    file_end_offset = 1000
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        qmin = self.params["qmin"]
        if "qnext" in self.state:
            self.state["qnext"] = max(self.state["qnext"], qmin)
        else:
            # qmin = 0 is magic value that uses alim instead. Not pretty.
            self.state["qnext"] = qmin if qmin > 0 else self.params["alim"]
        
        self.state.setdefault("rels_found", 0)
        self.state["rels_wanted"] = self.params["rels_wanted"]
        if self.state["rels_wanted"] == 0:
            # TODO: Choose sensible default value
            pass
        self.did_import = False
    
    def run(self):
        self.logger.info("Starting")
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        if "import" in self.params and not self.did_import:
            self.import_files(self.params["import"])
            self.did_import = True
        
        while self.get_nrels() < self.state["rels_wanted"]:
            q0 = self.state["qnext"]
            q1 = q0 + self.params["qrange"]
            # We use .gzip by default, unless set to no in parameters
            use_gz = ".gz" if self.params["gzip"] else ""
            outputfilename = \
                self.workdir.make_filename("%d-%d%s" % (q0, q1, use_gz))
            self.check_files_exist([outputfilename], "output",
                                   shouldexist=False)
            polyfilename = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
            factorbase = self.send_request(Request.GET_FACTORBASE_FILENAME)
            p = cadoprograms.Las(q0=q0, q1=q1,
                                 poly=polyfilename, factorbase=factorbase,
                                 out=outputfilename, stats_stderr=True,
                                 **self.progparams[0])
            self.submit_command(p, "%d-%d" % (q0, q1), commit=False)
            self.state.update({"qnext": q1}, commit=True)
        self.logger.info("Reached target of %d relations, now have %d",
                         self.state["rels_wanted"], self.get_nrels())
        self.logger.debug("Exit SievingTask.run(" + self.name + ")")
        return True
    
    def updateObserver(self, message):
        identifier = self.filter_notification(message)
        if not identifier:
            # This notification was not for me
            return False
        if self.handle_error_result(message):
            return False
        output_files = message.get_output_files()
        assert len(output_files) == 1
        stderrfilename = message.get_stderrfile(0)
        ok = self.add_file(output_files[0], stderrfilename, commit=False)
        self.verification(message.get_wu_id(), ok, commit=True)
        return True

    def add_file(self, filename, stats_filename=None, commit=True):
        use_stats_filename = stats_filename
        if stats_filename is None:
            self.logger.info("No statistics output received for file '%s', "
                             "have to scan file", filename)
            use_stats_filename = filename
        rels = self.parse_rel_count(use_stats_filename)
        if rels is None:
            return False
        update = {"rels_found": self.get_nrels() + rels}
        self.state.update(update, commit=False)
        self.add_output_files({filename: rels}, commit=commit)
        if stats_filename:
            self.parse_stats(stats_filename, commit=commit)
        self.logger.info("Found %d relations in '%s', total is now %d/%d",
                         rels, filename, self.get_nrels(),
                         self.state["rels_wanted"])
        return True
    
    def parse_rel_count(self, filename):
        (name, ext) = os.path.splitext(filename)
        try:
            if ext == ".gz":
                f = gzip.open(filename, "rb")
            else:
                size = os.path.getsize(filename)
                f = open(filename, "rb")
                f.seek(max(size - self.file_end_offset, 0))
        except (OSError, IOError) as e:
            if e.errno == 2: # No such file or directory
                self.logger.error("File '%s' does not exist", filename)
                return None
            else:
                raise
        try:
            for line in f:
                match = re.match(br"# Total (\d+) reports ", line)
                if match:
                    rels = int(match.group(1))
                    f.close()
                    return rels
        except (IOError, TypeError, structerror, EOFError) as e:
            if isinstance(e, IOError) and str(e) == "Not a gzipped file" or \
                    isinstance(e, TypeError) and str(e).startswith("ord() expected a character") or \
                    isinstance(e, structerror) and str(e).startswith("unpack requires a bytes object") or \
                    isinstance(e, EOFError) and str(e).startswith("Compressed file ended before"):
                self.logger.error("Error reading '%s' (corrupted?): %s", filename, e)
                return None
            else:
                raise
        f.close()
        self.logger.error("Number of relations not found in file '%s'", 
                          filename)
        return None
    
    def import_files(self, input_filename):
        if input_filename.startswith('@'):
            with open(input_filename[1:], "r") as f:
                filenames = f.read().splitlines()
        else:
            filenames = [input_filename]
        for filename in filenames:
            if filename in self.get_output_filenames():
                self.logger.info("Re-scanning file %s", filename)
                nrels = self.get_nrels() - self.get_nrels(filename)
                self.state.update({"rels_found": nrels}, commit=False)
                self.forget_output_filenames([filename], commit=True)
            self.add_file(filename)

    def get_statistics_as_strings(self):
        strings = ["Total number of relations: %d" % self.get_nrels()]
        strings += super().get_statistics_as_strings()
        return strings
    
    def get_nrels(self, filename=None):
        """ Return the number of relations found, either the total so far or
        for a given file
        """
        if filename is None:
            return self.state["rels_found"]
        else:
            # Fixme: don't access self.output_files directly
            return self.output_files[filename]
    
    def request_more_relations(self, target):
        wanted = self.state["rels_wanted"]
        if target > wanted:
            self.state["rels_wanted"] = target
            wanted = target
        nrels = self.get_nrels()
        if wanted > nrels:
            self.send_notification(Notification.WANT_TO_RUN, None)
            self.logger.info("New goal for number of relations is %d, "
                             "currently have %d. Need to sieve more",
                             wanted, nrels)
        else:
            self.logger.info("New goal for number of relations is %d, but "
                             "already have %d. No need to sieve more",
                             wanted, nrels)

    def get_total_cpu_or_real_time(self, is_cpu):
        """ Return number of seconds of cpu time spent by las """
        if is_cpu:
            return float(self.statistics.stats.get("stats_total_cpu_time", [0.])[0])
        else:
            return float(self.statistics.stats.get("stats_total_time", [0.])[0])


class Duplicates1Task(Task, FilesCreator, HasStatistics):
    """ Removes duplicate relations """
    @property
    def name(self):
        return "duplicates1"
    @property
    def title(self):
        return "Filtering - Duplicate Removal, splitting pass"
    @property
    def programs(self):
        return (cadoprograms.Duplicates1,)
    @property
    def progparam_override(self):
        return [["filelist", "prefix", "out"]]
    @property
    def paramnames(self):
        return self.join_params(super().paramnames, {"nslices_log": 1})
    @property
    def stat_conversions(self):
        # "End of read: 229176 relations in 0.9s -- 21.0 MB/s -- 253905.7 rels/s"
        # Without leading "# " !
        return (
        (
            "stats_dup1_time",
            float,
            "0",
            Statistics.add_list,
            re.compile(r"End of read: \d+ relations in %ss" % CAP_FP)
        ),
    )
    @property
    def stat_formats(self):
        return (
            ["CPU time for dup1: {stats_dup1_time[0]}s"],
        )
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.nr_slices = 2**self.params["nslices_log"]
        tablename = self.make_tablename("infiles")
        self.already_split_input = self.make_db_dict(tablename,
                                                     connection=self.db_connection)
        self.slice_relcounts = self.make_db_dict(self.make_tablename("counts"),
                                                 connection=self.db_connection)
        # Default slice counts to 0, in single DB commit
        self.slice_relcounts.setdefault(
            None, {str(i): 0 for i in range(0, self.nr_slices)})
    
    def run(self):
        self.logger.info("Starting")
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        # Check that previously split files were split into the same number
        # of pieces that we want now
        for (infile, parts) in self.already_split_input.items():
            if parts != self.nr_slices:
                # TODO: ask interactively (or by -recover) whether to delete
                # old output files and generate new ones, if input file is
                # still available
                # If input file is not available but the previously split
                # parts are, we could join them again... not sure if want
                raise Exception("%s was previously split into %d parts, "
                                "now %d parts requested",
                                infile, parts, self.nr_slices)
        
        # Check that previously split files do, in fact, exist.
        # FIXME: Do we want this? It may be slow when there are many files.
        # Reading the directory and comparing the lists would probably be
        # faster than individual lookups.
        self.check_files_exist(self.get_output_filenames(), "output",
                               shouldexist=True)
        
        siever_files = self.send_request(Request.GET_SIEVER_FILENAMES)
        newfiles = [f for f in siever_files
                    if not f in self.already_split_input]
        self.logger.debug ("new files to split are: %s", newfiles)
        
        if not newfiles:
            self.logger.info("No new files to split")
        else:
            self.logger.info("Splitting %d new files", len(newfiles))
            # TODO: can we recover from missing input files? Ask Sieving to
            # generate them again? Just ignore the missing ones?
            self.check_files_exist(newfiles, "input", shouldexist=True)
            # Split the new files
            if self.nr_slices == 1:
                # If we should split into only 1 part, we don't actually
                # split at all. We simply write the input file name
                # to the table of output files, so the next stages will
                # read the original siever output file, thus avoiding
                # having another copy of the data on disk. Since we don't
                # process the file at all, we need to ask the Siever task
                # for the relation count in the files
                # TODO: pass a list or generator expression in the request
                # here?
                total = self.slice_relcounts["0"]
                for f in newfiles:
                    total += self.send_request(Request.GET_SIEVER_RELCOUNT, f)
                self.slice_relcounts.update({"0": total}, commit=False)
                update1 = dict.fromkeys(newfiles, self.nr_slices)
                self.already_split_input.update(update1, commit=False)
                update2 = dict.fromkeys(newfiles, 0)
                self.add_output_files(update2, commit=True)
            else:
                # Make a task-specific subdirectory name under out working
                # directory
                outputdir = self.workdir.make_dirname()
                # Create this directory if it does not exist
                # self.logger.info("Creating directory %s", outputdir)
                outputdir.mkdir(parent=True)
                # Create a WorkDir instance for this subdirectory
                dup1_subdir = WorkDir(outputdir)
                # For each slice index [0, ..., nr_slices-1], create under the
                # subdirectory another directory for that slice's output files
                for i in range(0, self.nr_slices):
                    subdir = dup1_subdir.make_filename2(filename=str(i))
                    # self.logger.info("Creating directory %s", subdir)
                    subdir.mkdir(parent=True)
                run_counter = self.state.get("run_counter", 0)
                prefix = "dup1.%s" % run_counter
                (stdoutpath, stderrpath) = \
                        self.make_std_paths(cadoprograms.Duplicates1.name)
                if len(newfiles) <= 10:
                    p = cadoprograms.Duplicates1(*newfiles,
                                                 prefix=prefix,
                                                 out=outputdir,
                                                 stdout=str(stdoutpath),
                                                 stderr=str(stderrpath),
                                                 **self.progparams[0])
                else:
                    filelistname = self.workdir.make_filename("filelist")
                    with filelistname.open("w") as filelistfile:
                        filelistfile.write("\n".join(newfiles) + "\n")
                    p = cadoprograms.Duplicates1(filelist=filelistname,
                                                 prefix=prefix,
                                                 out=outputdir,
                                                 stdout=str(stdoutpath),
                                                 stderr=str(stderrpath),
                                                 **self.progparams[0])
                message = self.submit_command(p, "", log_errors=True)
                if message.get_exitcode(0) != 0:
                    raise Exception("Program failed")
                    # Check that the output files exist now
                    # TODO: How to recover from error? Presumably a dup1
                    # process failed, but that should raise a return code
                    # exception
                with stderrpath.open("r") as stderrfile:
                    stderr = stderrfile.read()
                outfilenames = self.parse_output_files(stderr)
                if not outfilenames:
                    self.logger.critical("No output files produced by %s",
                                         p.name)
                    return False
                self.logger.debug("Output file names: %s", outfilenames)
                self.check_files_exist(outfilenames.keys(), "output",
                                       shouldexist=True)
                self.state.update({"run_counter": run_counter + 1},
                                  commit=False)
                current_counts = self.parse_slice_counts(stderr)
                self.parse_stats(stderrpath, commit=False)
                # Add relation count from the newly processed files to the
                # relations-per-slice dict
                update1 = {str(idx): self.slice_relcounts[str(idx)] + 
                                     current_counts[idx]
                           for idx in range(self.nr_slices)}
                self.slice_relcounts.update(update1, commit=False)
                # Add the newly processed input files to the list of already
                # processed input files
                update2 = dict.fromkeys(newfiles, self.nr_slices)
                self.already_split_input.update(update2, commit=False)
                # Add the newly produced output files and commit everything
                self.add_output_files(outfilenames, commit=True)
        totals = ["%d: %d" % (i, self.slice_relcounts[str(i)])
                  for i in range(0, self.nr_slices)]
        self.logger.info("Relations per slice: %s", ", ".join(totals))
        self.logger.debug("Exit Duplicates1Task.run(" + self.name + ")")
        return True
    
    @staticmethod
    def parse_output_files(stderr):
        files = {}
        for line in stderr.splitlines():
            match = re.match(r'# Opening output file for slice (\d+) : (.+)$',
                             line)
            if match:
                (slicenr, filename) = match.groups()
                files[filename] = int(slicenr)
        return files
    
    def parse_slice_counts(self, stderr):
        """ Takes lines of text and looks for slice counts as printed by dup1
        """
        counts = [None] * self.nr_slices
        for line in stderr.splitlines():
            match = re.match(r'# slice (\d+) received (\d+) relations', line)
            if match:
                (slicenr, nrrels) = map(int, match.groups())
                if not counts[slicenr] is None:
                    raise Exception("Received two values for relation count "
                                    "in slice %d" % slicenr)
                counts[slicenr] = nrrels
        for (slicenr, nrrels) in enumerate(counts):
            if nrrels is None:
                raise Exception("Received no value for relation count in "
                                "slice %d" % slicenr)
        return counts
    
    def get_nr_slices(self):
        return self.nr_slices
    
    def get_nrels(self, idx):
        return self.slice_relcounts[str(idx)]
    
    def request_more_relations(self, target):
        self.send_notification(Notification.WANT_MORE_RELATIONS, target)
        self.send_notification(Notification.WANT_TO_RUN, None)

class Duplicates2Task(Task, FilesCreator, HasStatistics):
    """ Removes duplicate relations """
    @property
    def name(self):
        return "duplicates2"
    @property
    def title(self):
        return "Filtering - Duplicate Removal, removal pass"
    @property
    def programs(self):
        return (cadoprograms.Duplicates2,)
    @property
    def progparam_override(self):
        return [["poly", "rel_count", "badidealinfo", "renumber", "filelist"]]
    @property
    def paramnames(self):
        return self.join_params(super().paramnames, 
            {"dlp": False, "nslices_log": 1})
    @property
    def stat_conversions(self):
        # "End of read: 229176 relations in 0.9s -- 21.0 MB/s -- 253905.7 rels/s"
        # Without leading "# " !
        return (
        (
            "stats_dup2_time",
            float,
            "0",
            Statistics.add_list,
            re.compile(r"End of read: \d+ relations in %ss" % CAP_FP)
        ),
    )
    @property
    def stat_formats(self):
        return (
            ["CPU time for dup2: {stats_dup2_time[0]}s"],
        )
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.nr_slices = 2**self.params["nslices_log"]
        tablename = self.make_tablename("infiles")
        self.already_done_input = self.make_db_dict(tablename, connection=self.db_connection)
        self.slice_relcounts = self.make_db_dict(self.make_tablename("counts"), connection=self.db_connection)
        self.slice_relcounts.setdefault(
            None, {str(i): 0 for i in range(0, self.nr_slices)})
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        self.logger.info("Starting")
        input_nrel = 0
        for i in range(0, self.nr_slices):
            polyfilename = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
            renumber_filename = self.send_request(Request.GET_RENUMBER_FILENAME)
            files = self.send_request(Request.GET_DUP1_FILENAMES, i.__eq__)
            rel_count = self.send_request(Request.GET_DUP1_RELCOUNT, i)
            input_nrel += rel_count
            newfiles = [f for f in files if not f in self.already_done_input]
            if not newfiles:
                self.logger.info("No new files for slice %d, nothing to do", i)
                continue
            # If there are any new files in a slice, we remove duplicates on
            # the whole file set of the slice, as we currently cannot store
            # the duplicate removal state to be able to add more relations
            # in another pass
            # Forget about the previous output filenames of this slice
            # FIXME: Should we delete the files, too?
            self.forget_output_filenames(self.get_output_filenames(i.__eq__),
                                         commit=True)
            del(self.slice_relcounts[str(i)])
            name = "%s.slice%d" % (cadoprograms.Duplicates2.name, i)
            (stdoutpath, stderrpath) = \
                self.make_std_paths(name, do_increment=(i == 0))
             
            if self.params["dlp"]:
                badinfofilename = self.send_request(Request.GET_BADIDEALINFO_FILENAME)
            else:
                badinfofilename = None

            if len(files) <= 10:
                p = cadoprograms.Duplicates2(*files,
                                             poly=polyfilename,
                                             rel_count=rel_count,
                                             badidealinfo=badinfofilename,
                                             renumber=renumber_filename,
                                             stdout=str(stdoutpath),
                                             stderr=str(stderrpath),
                                             **self.progparams[0])
            else:
                filelistname = self.workdir.make_filename("filelist")
                with filelistname.open("w") as filelistfile:
                    filelistfile.write("\n".join(files) + "\n")
                p = cadoprograms.Duplicates2(poly=polyfilename,
                                             rel_count=rel_count,
                                             badidealinfo=badinfofilename,
                                             renumber=renumber_filename,
                                             filelist=filelistname,
                                             stdout=str(stdoutpath),
                                             stderr=str(stderrpath),
                                             **self.progparams[0])
            message = self.submit_command(p, "", log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            with stderrpath.open("r") as stderrfile:
                nr_rels = self.parse_remaining(stderrfile)
            # Mark input file names and output file names
            for f in files:
                self.already_done_input[f] = True
            outfilenames = {f:i for f in files}
            self.add_output_files(outfilenames, commit=True)
            # Disabled for now, there are multiple lines of the same format
            # which we can't parse atm.
            # self.parse_stats(stderrpath)
            self.logger.info("%d unique relations remain on slice %d",
                             nr_rels, i)
            self.slice_relcounts[str(i)] = nr_rels
        self.update_ratio(input_nrel, self.get_nrels())
        self.logger.info("%d unique relations remain in total",
                         self.get_nrels())
        self.logger.debug("Exit Duplicates2Task.run(" + self.name + ")")
        return True
    
    @staticmethod
    def parse_remaining(text):
        # "     112889 remaining relations"
        for line in text:
            match = re.match(r'At the end:\s*(\d+) remaining relations', line)
            if match:
                remaining = int(match.group(1))
                return remaining
        raise Exception("Received no value for remaining relation count")

    def get_nrels(self):
        nrels = 0
        for i in range(0, self.nr_slices):
            nrels += self.slice_relcounts[str(i)]
        return nrels
    
    def update_ratio(self, input_nrel, output_nrel):
        last_input_nrel = self.state.get("last_input_nrel", 0)
        last_output_nrel = self.state.get("last_output_nrel", 0)
        new_in = input_nrel - last_input_nrel
        new_out = output_nrel - last_output_nrel
        if new_in < 0:
            self.logger.error("Negative number %d of new relations?", new_in)
            return
        if new_in == 0:
            return
        if new_out > new_in:
            self.logger.error("More new output relations (%d) than input (%d)?",
                              new_out, new_in)
            return
        ratio = new_out / new_in
        self.logger.info("Of %d newly added relations %d were unique "
                         "(ratio %f)", new_in, new_out, ratio)
        self.state.update({"last_input_nrel": input_nrel,
            "last_output_nrel": output_nrel, "unique_ratio": ratio})
    
    def request_more_relations(self, target):
        nrels = self.get_nrels()
        if target <= nrels:
            return
        additional_out = target - nrels
        ratio = self.state.get("unique_ratio", 1.)
        additional_in = int(additional_out / ratio)
        newtarget = self.state["last_input_nrel"] + additional_in
        self.logger.info("Got request for %d (%d additional) output relations, "
                         "estimate %d (%d additional) needed in input",
                         target, additional_out, newtarget, additional_in)
        self.send_notification(Notification.WANT_MORE_RELATIONS, newtarget)
        self.send_notification(Notification.WANT_TO_RUN, None)


class PurgeTask(Task):
    """ Removes singletons and computes excess """
    @property
    def name(self):
        return "purge"
    @property
    def title(self):
        return "Filtering - Singleton removal"
    @property
    def programs(self):
        return (cadoprograms.Purge,)
    @property
    def progparam_override(self):
        return [["nrels", "out", "minindex", "outdel", "nprimes", "filelist"]]
    @property
    def paramnames(self):
        return self.join_params(super().paramnames, 
            {"dlp": False, "alim": int, "rlim": int, "gzip": True})
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.state.setdefault("input_nrels", 0)
    
    def run(self):
        self.logger.info("Starting")
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        nfree = self.send_request(Request.GET_FREEREL_RELCOUNT)
        nunique = self.send_request(Request.GET_UNIQUE_RELCOUNT)
        minlim = min(self.params["alim"], self.params["rlim"])
        minindex = int(2. * minlim / (log(minlim) - 1))
        nprimes = self.send_request(Request.GET_RENUMBER_PRIMECOUNT)
        if not nunique:
            self.logger.critical("No unique relation count received")
            return False
        input_nrels = nfree + nunique
        
        if "purgedfile" in self.state and \
                input_nrels == self.state["input_nrels"]:
            self.logger.info("Already have a purged file, and no new input "
                             "relations available. Nothing to do")
            self.send_notification(Notification.HAVE_ENOUGH_RELATIONS, None)
            return True
        
        self.state.pop("purgedfile", None)
        self.state.pop("input_nrels", None)
        
        self.logger.info("Reading %d unique and %d free relations, total %d"
                         % (nunique, nfree, input_nrels))
        use_gz = ".gz" if self.params["gzip"] else ""
        purgedfile = self.workdir.make_filename("purged" + use_gz)
        if self.params["dlp"]:
            relsdelfile = self.workdir.make_filename("relsdel" + use_gz)
        else:
            relsdelfile = None
        freerel_filename = self.send_request(Request.GET_FREEREL_FILENAME)
        unique_filenames = self.send_request(Request.GET_UNIQUE_FILENAMES)
        files = unique_filenames + [str(freerel_filename)]
        (stdoutpath, stderrpath) = self.make_std_paths(cadoprograms.Purge.name)
        
        if self.params["dlp"]:
            nmaps = self.send_request(Request.GET_NMAPS)
            self.progparams[0]["keep"] = nmaps
        
        if len(files) <= 10:
            p = cadoprograms.Purge(*files,
                                   nrels=input_nrels, out=purgedfile,
                                   outdel=relsdelfile,
                                   minindex=minindex, nprimes=nprimes,
                                   stdout=str(stdoutpath),
                                   stderr=str(stderrpath),
                                   **self.progparams[0])
        else:
            filelistname = self.workdir.make_filename("filelist")
            with filelistname.open("w") as filelistfile:
                filelistfile.write("\n".join(files))
            p = cadoprograms.Purge(nrels=input_nrels,
                                   out=purgedfile, minindex=minindex,
                                   outdel=relsdelfile,
                                   nprimes=nprimes, filelist=filelistname,
                                   stdout=str(stdoutpath),
                                   stderr=str(stderrpath),
                                   **self.progparams[0])
        message = self.submit_command(p, "")
        stdout = message.read_stdout(0).decode('utf-8')
        stderr = message.read_stderr(0).decode('utf-8')
        if self.parse_stderr(stderr, input_nrels):
            stats = self.parse_stdout(stdout)
            self.logger.info("After purge, %d relations with %d primes remain "
                             "with weight %s and excess %s", *stats)
            update = {"purgedfile": purgedfile.get_wdir_relative(),
                      "input_nrels": input_nrels }
            if self.params["dlp"]:
                update["relsdelfile"] = relsdelfile.get_wdir_relative()
            self.state.update(update)
            self.logger.info("Have enough relations")
            self.send_notification(Notification.HAVE_ENOUGH_RELATIONS, None)
        else:
            self.logger.info("Not enough relations")
            self.request_more_relations(nunique)
        self.logger.debug("Exit PurgeTask.run(" + self.name + ")")
        return True
    
    def request_more_relations(self, nunique):
        r"""
        We want an excess of $e = 0.1 n_p$,
        with $n_p$ primes among $n_r$ relations in the output file.
        
        We have $\Delta_r$ additional relations in the output file per
        additional input relation, and $\Delta_p$ additional primes in
        the output file per additional input relation.
        
        We need $a$ additional relations so that
        \begin{eqnarray*}
            n_r + a  \Delta_r - (n_p + a  \Delta_p) & = & 0.1  (n_p + a  \Delta_p) \\
            n_r + a  \Delta_r & = & 1.1  (n_p + a  \Delta_p) \\
            n_r - 1.1 n_p & = & 1.1 a \Delta_p - a \Delta_r \\
            \frac{n_r - 1.1 n_p}{1.1 \Delta_p - \Delta_r} & = & a \\
        \end{eqnarray*}
        """
        
        additional = nunique * 0.1
        if "delta_r" in self.state:
            excess = self.state["last_input_nrels"] - \
                self.state["last_input_nprimes"]
            n_r = self.state["last_output_nrels"]
            n_p = self.state["last_output_nprimes"]
            delta_r = self.state["delta_r"]
            delta_p = self.state["delta_p"]
            if abs(1.1 * delta_p - delta_r) > 0.01:
                a = (n_r - 1.1 * n_p) / (1.1 * delta_p - delta_r)
                self.logger.info("a = %f, excess = %d", a, excess)
                # Use the negative excess among the input relations as an
                # upper limit on the additional relations needed, to keep
                # a small donominator from causing a huge value
                if excess < 0 and a > -excess:
                    a = -excess
                if a > 10000. and a < nunique * 0.5:
                    additional = int(a)
        # Always request at least 10k more
        additional = max(additional, 10000)
        
        self.logger.info("Requesting %d additional relations", additional)
        self.send_notification(Notification.WANT_MORE_RELATIONS,
                               nunique + additional)
        self.send_notification(Notification.WANT_TO_RUN, None)
    
    def get_purged_filename(self):
        return self.get_state_filename("purgedfile")
    
    def get_relsdel_filename(self):
        return self.get_state_filename("relsdelfile")
    
    def parse_stderr(self, stderr, input_nrels):
        # If stderr ends with
        # b'excess < 0.10 * #primes. See -required_excess argument.'
        # then we need more relations from filtering and return False
        input_nprimes = None
        have_enough = True
        # not_enough1 = re.compile(r"excess < (\d+.\d+) \* #primes")
        not_enough1 = re.compile(r"\(excess / nprimes\) = \d+.?\d* < \d+.?\d*. "
                                 r"See -required_excess argument.")
        not_enough2 = re.compile(r"number of relations <= number of ideals")
        nrels_nprimes = re.compile(r"\s*nrels=(\d+), nprimes=(\d+); "
                                   r"excess=(-?\d+)")
        for line in stderr.splitlines():
            match = not_enough1.match(line)
            if match:
                have_enough = False
                break
            if not_enough2.match(line):
                have_enough = False
                break
            match = nrels_nprimes.match(line)
            if not match is None:
                (nrels, nprimes, excess) = map(int, match.groups())
                assert nrels - nprimes == excess
                # The first occurrence of the message counts input relations
                if input_nprimes is None:
                    assert input_nrels == nrels
                    input_nprimes = nprimes
        
        # At this point we shoud have:
        # input_nrels, input_nprimes: rels and primes among input
        # nrels, nprimes, excess: rels and primes when purging stopped
        if not input_nprimes is None:
            self.update_excess_per_input(input_nrels, nrels, nprimes)
        return have_enough
    
    def update_excess_per_input(self, input_nrels, nrels, nprimes):
        if input_nrels == 0:
            return # Nothing sensible that we can do
        last_input_nrels = self.state.get("last_input_nrels", 0)
        last_nrels = self.state.get("last_output_nrels", 0)
        last_nprimes = self.state.get("last_output_nprimes", 0)
        if last_input_nrels >= input_nrels:
            self.logger.warn("Previously stored input nrels (%d) is no "
                             "smaller than value from new run (%d)",
                             last_input_nrels, input_nrels)
            return
        if nrels <= last_nrels:
            self.logger.warn("Previously stored nrels (%d) is no "
                             "smaller than value from new run (%d)",
                             last_nrels, nrels)
            return
        if nprimes <= last_nprimes:
            self.logger.warn("Previously stored nprimes (%d) is no "
                             "smaller than value from new run (%d)",
                             last_nprimes, nprimes)
            return
        self.logger.info("Previous run had %d input relations and ended with "
                         "%d relations and %d primes, new run had %d input "
                         "relations and ended with %d relations and %d primes",
                         last_input_nrels, last_nrels, last_nprimes,
                         input_nrels, nrels, nprimes)
        delta_r = (nrels - last_nrels) / (input_nrels - last_input_nrels)
        delta_p = (nprimes - last_nprimes) / (input_nrels - last_input_nrels)
        self.logger.info("Gained %f output relations and %f primes per input "
                         "relation", delta_r, delta_p)
        update = {"last_output_nrels": nrels, "last_output_nprimes": nprimes,
                  "last_input_nrels": input_nrels, "delta_r": delta_r,
                  "delta_p": delta_p}
        self.state.update(update)
    
    def parse_stdout(self, stdout):
        # Program stdout is expected in the form:
        #   Final values:
        #   nrels=23105 nprimes=22945 excess=160
        #   weight=382433 weight*nrels=8.84e+09
        # but we allow some extra whitespace
        r = {}
        keys = ("nrels", "nprimes", "weight", "excess")
        for line in stdout.splitlines():
            for key in keys:
                # Match the key at the start of a line, or after a whitespace
                # Note: (?:) is non-capturing group
                match = re.search(r"(?:^|\s)%s\s*=\s*(\d+)" % key, line)
                if match:
                    if key in r:
                        raise Exception("Found multiple values for %s" % key)
                    r[key] = int(match.group(1))
        for key in keys:
            if not key in r:
                raise Exception("%s: output of %s did not contain value for %s: %s"
                                % (self.title, self.programs[0].name, key, stdout))
        return [r[key] for key in keys]

class MergeDLPTask(Task):
    """ Merges relations """
    @property
    def name(self):
        return "merge"
    @property
    def title(self):
        return "Filtering - Merging"
    @property
    def programs(self):
        return (cadoprograms.MergeDLP, cadoprograms.ReplayDLP)
    @property
    def progparam_override(self):
        return [["mat", "out", "keep"],
            ["purged", "ideals", "history", "index", "out"]]
    @property
    def paramnames(self):
        return self.join_params(super().paramnames, {"gzip": True})
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.progparams[0]["skip"] = 0

    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        if not "mergedfile" in self.state:
            self.logger.info("Starting")
            if "idealfile" in self.state:
                del(self.state["idealfile"])
            if "indexfile" in self.state:
                del(self.state["indexfile"])
            if "mergedfile" in self.state:
                del(self.state["mergedfile"])
            if "densefile" in self.state:
                del(self.state["densefile"])
            
            purged_filename = self.send_request(Request.GET_PURGED_FILENAME)
            keep = self.send_request(Request.GET_NMAPS)
            # We use .gzip by default, unless set to no in parameters
            use_gz = ".gz" if self.params["gzip"] else ""
            historyfile = self.workdir.make_filename("history" + use_gz)
            (stdoutpath, stderrpath) = self.make_std_paths(cadoprograms.MergeDLP.name)
            p = cadoprograms.MergeDLP(mat=purged_filename,
                                   out=historyfile,
                                   keep=keep,
                                   stdout=str(stdoutpath),
                                   stderr=str(stderrpath),
                                   **self.progparams[0])
            message = self.submit_command(p, "", log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            
            indexfile = self.workdir.make_filename("index" + use_gz)
            mergedfile = self.workdir.make_filename("small.txt")
            (stdoutpath, stderrpath) = self.make_std_paths(cadoprograms.Replay.name)
            idealfile = self.workdir.make_filename("ideal")
            p = cadoprograms.ReplayDLP(purged=purged_filename,
                                    ideals=idealfile,
                                    history=historyfile, index=indexfile,
                                    out=mergedfile, stdout=str(stdoutpath),
                                    stderr=str(stderrpath),
                                    **self.progparams[1])
            message = self.submit_command(p, "", log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            
            if not idealfile.isfile():
                raise Exception("Output file %s does not exist" % idealfile)
            if not indexfile.isfile():
                raise Exception("Output file %s does not exist" % indexfile)
            if not mergedfile.isfile():
                raise Exception("Output file %s does not exist" % mergedfile)
            update = {"indexfile": indexfile.get_wdir_relative(),
                      "mergedfile": mergedfile.get_wdir_relative(),
                      "idealfile": idealfile.get_wdir_relative()}
            densefilename = self.workdir.make_filename("small.dense.bin")
            if densefilename.isfile():
                update["densefile"] = densefilename.get_wdir_relative()
            self.state.update(update)
            
        self.logger.debug("Exit MergeTask.run(" + self.name + ")")
        return True
    
    def get_index_filename(self):
        return self.get_state_filename("indexfile")
    
    def get_ideal_filename(self):
        return self.get_state_filename("idealfile")
    
    def get_merged_filename(self):
        return self.get_state_filename("mergedfile")
    
    def get_dense_filename(self):
        return self.get_state_filename("densefile")


class MergeTask(Task):
    """ Merges relations """
    @property
    def name(self):
        return "merge"
    @property
    def title(self):
        return "Filtering - Merging"
    @property
    def programs(self):
        return (cadoprograms.Merge, cadoprograms.Replay)
    @property
    def progparam_override(self):
        return [["mat", "out"], ["purged", "history", "index"]]
    @property
    def paramnames(self):
        return self.join_params(super().paramnames,  \
            {"skip": None, "keep": None, "gzip": True})
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        skip = int(self.progparams[0].get("skip", 32))
        self.progparams[0].setdefault("skip", skip)
        self.progparams[0].setdefault("keep", skip + 128)

    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        if not "mergedfile" in self.state:
            self.logger.info("Starting")
            if "indexfile" in self.state:
                del(self.state["indexfile"])
            if "mergedfile" in self.state:
                del(self.state["mergedfile"])
            if "densefile" in self.state:
                del(self.state["densefile"])
            
            purged_filename = self.send_request(Request.GET_PURGED_FILENAME)
            # We use .gzip by default, unless set to no in parameters
            use_gz = ".gz" if self.params["gzip"] else ""
            historyfile = self.workdir.make_filename("history" + use_gz)
            (stdoutpath, stderrpath) = self.make_std_paths(cadoprograms.Merge.name)
            p = cadoprograms.Merge(mat=purged_filename,
                                   out=historyfile,
                                   stdout=str(stdoutpath),
                                   stderr=str(stderrpath),
                                   **self.progparams[0])
            message = self.submit_command(p, "", log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            
            indexfile = self.workdir.make_filename("index" + use_gz)
            mergedfile = self.workdir.make_filename("small.bin")
            (stdoutpath, stderrpath) = self.make_std_paths(cadoprograms.Replay.name)
            p = cadoprograms.Replay(purged=purged_filename,
                                    history=historyfile, index=indexfile,
                                    out=mergedfile, stdout=str(stdoutpath),
                                    stderr=str(stderrpath),
                                    **self.progparams[1])
            message = self.submit_command(p, "", log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            
            if not indexfile.isfile():
                raise Exception("Output file %s does not exist" % indexfile)
            if not mergedfile.isfile():
                raise Exception("Output file %s does not exist" % mergedfile)
            update = {"indexfile": indexfile.get_wdir_relative(),
                      "mergedfile": mergedfile.get_wdir_relative()}
            densefilename = self.workdir.make_filename("small.dense.bin")
            if densefilename.isfile():
                update["densefile"] = densefilename.get_wdir_relative()
            self.state.update(update)
            
        self.logger.debug("Exit MergeTask.run(" + self.name + ")")
        return True
    
    def get_index_filename(self):
        return self.get_state_filename("indexfile")
    
    def get_merged_filename(self):
        return self.get_state_filename("mergedfile")
    
    def get_dense_filename(self):
        return self.get_state_filename("densefile")

class NmbrthryTask(Task):
    """ Number theory tasks for dlp"""
    @property
    def name(self):
        return "magmanmbrthry"
    @property
    def title(self):
        return "Number Theory for DLP"
    @property
    def programs(self):
        return (cadoprograms.MagmaNmbrthry,)
    @property
    def progparam_override(self):
        return [["poly", "badidealinfo", "badideals"]]
    @property
    def paramnames(self):
        return super().paramnames
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)

    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        self.logger.info("Starting")
        
        badfile = self.workdir.make_filename("badideals")
        badinfofile = self.workdir.make_filename("badidealinfo")
        polyfile = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
        (stdoutpath, stderrpath) = self.make_std_paths(cadoprograms.MagmaNmbrthry.name)
        p = cadoprograms.MagmaNmbrthry(poly=polyfile,
                               badidealinfo=badinfofile,
                               badideals=badfile,
                               stdout=str(stdoutpath),
                               stderr=str(stderrpath),
                               **self.progparams[0])
        message = self.submit_command(p, "", log_errors=True)
        if message.get_exitcode(0) != 0:
            raise Exception("Program failed")

        stdout = message.read_stdout(0).decode("utf-8")
        update = {}
        for line in stdout.splitlines():
            match = re.match(r'ell (\d+)', line)
            if match:
                update["ell"] = int(match.group(1))
            match = re.match(r'smexp (\d+)', line)
            if match:
                update["smexp"] = int(match.group(1))
            match = re.match(r'nmaps (\d+)', line)
            if match:
                update["nmaps"] = int(match.group(1))
        update["badinfofile"] = badinfofile.get_wdir_relative()
        update["badfile"] = badfile.get_wdir_relative()
        
        if not "ell" in update:
            raise Exception("Stdout does not give ell")
        if not "smexp" in update:
            raise Exception("Stdout does not give smexp")
        if not "nmaps" in update:
            raise Exception("Stdout does not give nmaps")
        if not badfile.isfile():
            raise Exception("Output file %s does not exist" % badfile)
        if not badinfofile.isfile():
            raise Exception("Output file %s does not exist" % badinfofile)
        # Update the state entries atomically
        self.state.update(update)

        self.logger.info("Will compute Dlog modulo %s", self.state["ell"])
        self.logger.debug("Exit NmbrthryTask.run(" + self.name + ")")
        return True

    def get_badinfo_filename(self):
        return self.get_state_filename("badinfofile")
    
    def get_bad_filename(self):
        return self.get_state_filename("badfile")
    
    def get_ell(self):
        return self.state["ell"]
    
    def get_smexp(self):
        return self.state["smexp"]
    
    def get_nmaps(self):
        return self.state["nmaps"]


class LinAlgDLPTask(Task):
    """ Runs the linear algebra step for dlp"""
    @property
    def name(self):
        return "magmalinalg"
    @property
    def title(self):
        return "Linear Algebra for DLP"
    @property
    def programs(self):
        return (cadoprograms.MagmaLinalg,)
    @property
    def progparam_override(self):
        return [["sparsemat", "ker", "sm", "ell", "nmaps"]]
    @property
    def paramnames(self):
        return super().paramnames
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)

    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        if not "kerfile" in self.state: 
            self.logger.info("Starting")
            
            kerfile = self.workdir.make_filename("ker")
            mergedfile = self.send_request(Request.GET_MERGED_FILENAME)
            smfile = self.send_request(Request.GET_SM_FILENAME)
            gorder = self.send_request(Request.GET_ELL)
            nmaps = self.send_request(Request.GET_NMAPS)
            (stdoutpath, stderrpath) = self.make_std_paths(cadoprograms.MagmaLinalg.name)
            p = cadoprograms.MagmaLinalg(sparsemat=mergedfile,
                                   ker=kerfile,
                                   sm=smfile,
                                   ell=gorder,
                                   nmaps=nmaps,
                                   stdout=str(stdoutpath),
                                   stderr=str(stderrpath),
                                   **self.progparams[0])
            message = self.submit_command(p, "", log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            
            if not kerfile.isfile():
                raise Exception("Output file %s does not exist" % kerfile)
            update = {"kerfile": kerfile.get_wdir_relative()}
            self.state.update(update)
            
        self.logger.debug("Exit LinAlgDLPTask.run(" + self.name + ")")
        return True
    
    def get_kernel_filename(self):
        return self.get_state_filename("kerfile")


class LinAlgTask(Task, HasStatistics):
    """ Runs the linear algebra step """
    @property
    def name(self):
        return "bwc"
    @property
    def title(self):
        return "Linear Algebra"
    @property
    def programs(self):
        return (cadoprograms.BWC,)
    @property
    def progparam_override(self):
        return [["complete", "matrix",  "wdir", "nullspace"]]
    @property
    def paramnames(self):
        return super().paramnames
    @property
    def stat_conversions(self):
        return (
        (
            "krylov_time",
            float,
            "0",
            Statistics.add_list,
            re.compile(r"krylov done, N=\d+ ; CPU: %s" % CAP_FP)
        ),
        (
            "krylov_comm",
            float,
            "0",
            Statistics.add_list,
            re.compile(r"krylov done, N=\d+ ; COMM: %s" % CAP_FP)
        ),
        (
            "lingen_time",
            float,
            "0",
            Statistics.add_list,
            re.compile(r"Total computation took %s" % CAP_FP)
        ),
        (
            "mksol_time",
            float,
            "0",
            Statistics.add_list,
            re.compile(r"mksol done, N=\d+ ; CPU: %s" % CAP_FP)
        ),
        (
            "mksol_comm",
            float,
            "0",
            Statistics.add_list,
            re.compile(r"mksol done, N=\d+ ; COMM: %s" % CAP_FP)
        ),
    )
    @property
    def stat_formats(self):
        return (
            ["Krylov: CPU time {krylov_time[0]}", ", COMM time {krylov_comm[0]}"],
            ["Lingen CPU time {lingen_time[0]}"],
            ["Mksol: CPU time {mksol_time[0]}", ", COMM time {mksol_comm[0]}"],
        )
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        if not "dependency" in self.state:
            self.logger.info("Starting")
            workdir = self.workdir.make_dirname()
            workdir.mkdir(parent=True)
            mergedfile = self.send_request(Request.GET_MERGED_FILENAME)
            if mergedfile is None:
                self.logger.critical("No merged file received.")
                return False
            (stdoutpath, stderrpath) = self.make_std_paths(cadoprograms.BWC.name)
            matrix = mergedfile.realpath()
            wdir = workdir.realpath()
            p = cadoprograms.BWC(complete=True,
                                 matrix=matrix,  wdir=wdir, nullspace="left",
                                 stdout=str(stdoutpath),
                                 stderr=str(stderrpath),
                                 **self.progparams[0])
            message = self.submit_command(p, "", log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            dependencyfilename = self.workdir.make_filename("W", use_subdir=True)
            if not dependencyfilename.isfile():
                raise Exception("Kernel file %s does not exist" % dependencyfilename)
            self.parse_stats(stdoutpath, commit=False)
            update = {"dependency": dependencyfilename.get_wdir_relative()}
            self.state.update(update, commit=True)
        self.logger.debug("Exit LinAlgTask.run(" + self.name + ")")
        return True

    def get_dependency_filename(self):
        return self.get_state_filename("dependency")
    
    def get_prefix(self):
        return "%s%s%s.%s" % (self.params["workdir"].rstrip(os.sep), os.sep,
                              self.params["name"], "dep")


class CharactersTask(Task):
    """ Computes Quadratic Characters """
    @property
    def name(self):
        return "characters"
    @property
    def title(self):
        return "Quadratic Characters"
    @property
    def programs(self):
        return (cadoprograms.Characters,)
    @property
    def progparam_override(self):
        return [["poly", "purged", "index", "wfile", "out", "heavyblock"]]
    @property
    def paramnames(self):
        return super().paramnames

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        if not "kernel" in self.state:
            self.logger.info("Starting")
            polyfilename = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
            kernelfilename = self.workdir.make_filename("kernel")
            
            purgedfilename = self.send_request(Request.GET_PURGED_FILENAME)
            indexfilename = self.send_request(Request.GET_INDEX_FILENAME)
            densefilename = self.send_request(Request.GET_DENSE_FILENAME)
            dependencyfilename = self.send_request(Request.GET_DEPENDENCY_FILENAME)
            
            (stdoutpath, stderrpath) = \
                    self.make_std_paths(cadoprograms.Characters.name)
            p = cadoprograms.Characters(poly=polyfilename,
                    purged=purgedfilename, index=indexfilename,
                    wfile=dependencyfilename, out=kernelfilename,
                    heavyblock=densefilename, stdout=str(stdoutpath),
                    stderr=str(stderrpath),
                    **self.progparams[0])
            message = self.submit_command(p, "", log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            if not kernelfilename.isfile():
                raise Exception("Output file %s does not exist" % kernelfilename)
            self.state["kernel"] = kernelfilename.get_wdir_relative()
        self.logger.debug("Exit CharactersTask.run(" + self.name + ")")
        return True
    
    def get_kernel_filename(self):
        return self.get_state_filename("kernel")


class SqrtTask(Task):
    """ Runs the square root """
    @property
    def name(self):
        return "sqrt"
    @property
    def title(self):
        return "Square Root"
    @property
    def programs(self):
        return (cadoprograms.Sqrt,)
    @property
    def progparam_override(self):
        return [["ab", "poly", "purged", "index", "kernel", "prefix", "rat",
                 "alg", "gcd", "dep"]]
    @property
    def paramnames(self):
        return self.join_params(super().paramnames, {"N": int, "gzip": True})
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.factors = self.make_db_dict(self.make_tablename("factors"), connection=self.db_connection)
        self.add_factor(self.params["N"])
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        if not self.is_done():
            self.logger.info("Starting")
            polyfilename = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
            purgedfilename = self.send_request(Request.GET_PURGED_FILENAME)
            indexfilename = self.send_request(Request.GET_INDEX_FILENAME)
            kernelfilename = self.send_request(Request.GET_KERNEL_FILENAME)
            prefix = self.send_request(Request.GET_LINALG_PREFIX)
            if self.params["gzip"]:
                prefix += ".gz"
            (stdoutpath, stderrpath) = \
                self.make_std_paths(cadoprograms.Sqrt.name)
            self.logger.info("Creating file of (a,b) values")
            p = cadoprograms.Sqrt(ab=True,
                    poly=polyfilename, purged=purgedfilename,
                    index=indexfilename, kernel=kernelfilename,
                    prefix=prefix, stdout=str(stdoutpath),
                    stderr=str(stderrpath), **self.progparams[0])
            message = self.submit_command(p, "", log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            
            while not self.is_done():
                dep = self.state.get("next_dep", 0)
                self.logger.info("Trying dependency %d", dep)
                (stdoutpath, stderrpath) = \
                    self.make_std_paths(cadoprograms.Sqrt.name)
                p = cadoprograms.Sqrt(ab=False, rat=True,
                        alg=True, gcd=True, dep=dep, poly=polyfilename,
                        purged=purgedfilename, index=indexfilename,
                        kernel=kernelfilename, prefix=prefix,
                        stdout=str(stdoutpath), stderr=str(stderrpath), 
                        **self.progparams[0])
                message = self.submit_command(p, "dep%d" % dep, log_errors=True)
                if message.get_exitcode(0) != 0:
                    raise Exception("Program failed")
                with stdoutpath.open("r") as stdoutfile:
                    stdout = stdoutfile.read()
                lines = stdout.splitlines()
                # Skip last factor which cannot produce a new split on top
                # of what the smaller factors did
                for line in lines[:-1]:
                    if line == "Failed":
                        break
                    self.add_factor(int(line))
                self.state["next_dep"] = dep+1
            self.logger.info("finished")
        self.logger.info("Factors: %s" % " ".join(self.get_factors()))
        self.logger.debug("Exit SqrtTask.run(" + self.name + ")")
        return True
    
    def is_done(self):
        for (factor, isprime) in self.factors.items():
            if not isprime:
                return False
        return True
    
    def add_factor(self, factor):
        assert factor > 0
        if str(factor) in self.factors:
            return
        for oldfac in list(map(int, self.factors.keys())):
            g = gcd(factor, oldfac)
            if 1 < g and g < factor:
                self.add_factor(g)
                self.add_factor(factor // g)
                break
            if 1 < g and g < oldfac:
                # We get here only if newfac is a proper factor of oldfac
                assert factor == g
                del(self.factors[str(oldfac)])
                self.add_factor(g)
                self.add_factor(oldfac // g)
                break
        else:
            # We get here if the new factor is coprime to all previously
            # known factors
            isprime = SqrtTask.miller_rabin_tests(factor, 10)
            self.factors[str(factor)] = isprime
    
    def get_factors(self):
        return self.factors.keys()
    
    @staticmethod
    def miller_rabin_pass(number, base):
        """
        >>> SqrtTask.miller_rabin_pass(3, 2)
        True
        >>> SqrtTask.miller_rabin_pass(9, 2)
        False
        >>> SqrtTask.miller_rabin_pass(91, 2)
        False
        >>> SqrtTask.miller_rabin_pass(1009, 2)
        True
        >>> SqrtTask.miller_rabin_pass(10000000019, 2)
        True
        >>> SqrtTask.miller_rabin_pass(10000000019*10000000021, 2)
        False
        
        # Check some pseudoprimes. First a few Fermat pseudoprimes which
        # Miller-Rabin should recognize as composite
        >>> SqrtTask.miller_rabin_pass(341, 2)
        False
        >>> SqrtTask.miller_rabin_pass(561, 2)
        False
        >>> SqrtTask.miller_rabin_pass(645, 2)
        False
        
        # Now some strong pseudo-primes
        >>> SqrtTask.miller_rabin_pass(2047, 2)
        True
        >>> SqrtTask.miller_rabin_pass(703, 3)
        True
        >>> SqrtTask.miller_rabin_pass(781, 5)
        True
        """
        po2 = 0
        exponent = number - 1
        while exponent % 2 == 0:
            exponent >>= 1
            po2 += 1
        
        result = pow(base, exponent, number)
        if result == 1:
            return True
        for i in range(0, po2 - 1):
            if result == number - 1:
                return True
            result = pow(result, 2, number)
        return result == number - 1
    
    @staticmethod
    def miller_rabin_tests(number, passes):
        if number <= 3:
            return number >= 2
        if number % 2 == 0:
            return False
        for i in range(0, passes):
            # random.randrange(n) produces random integer in [0, n-1].
            # We want [2, n-2]
            base = random.randrange(number - 3) + 2
            if not SqrtTask.miller_rabin_pass(number, base):
                return False
        return True

    @staticmethod
    def nextprime(N):
        """ Return the smallest strong probable prime no smaller than N
        >>> prps = [SqrtTask.nextprime(i) for i in range(30)]
        >>> prps == [2, 2, 2, 3, 5, 5, 7, 7, 11, 11, 11, 11, 13, 13, 17, 17, \
                     17, 17, 19, 19, 23, 23, 23, 23, 29, 29, 29, 29, 29, 29]
        True
        """
        if N <= 2:
            return 2
        if N % 2 == 0:
            N += 1
        while not SqrtTask.miller_rabin_tests(N, 5):
            N += 2     
        return N

class SMTask(Task):
    """ Computes Schirokauer Maps """
    @property
    def name(self):
        return "sm"
    @property
    def title(self):
        return "Schirokauer Maps"
    @property
    def programs(self):
        return (cadoprograms.SM,)
    @property
    def progparam_override(self):
        return [["poly", "renumber", "purged", "index", "ell", "smexp",
                 "nmaps", "out"]]
    @property
    def paramnames(self):
        return super().paramnames

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)
        
        if not "sm" in self.state:
            self.logger.info("Starting")
            polyfilename = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
            renumberfilename = self.send_request(Request.GET_RENUMBER_FILENAME)
            badidealinfofilename = self.send_request(Request.GET_BADIDEALINFO_FILENAME)
            purgedfilename = self.send_request(Request.GET_PURGED_FILENAME)
            indexfilename = self.send_request(Request.GET_INDEX_FILENAME)
            nmaps = self.send_request(Request.GET_NMAPS)
            smfilename = self.workdir.make_filename("sm")

            gorder = self.send_request(Request.GET_ELL)
            smexp = self.send_request(Request.GET_SMEXP)
            
            (stdoutpath, stderrpath) = \
                    self.make_std_paths(cadoprograms.SM.name)
            p = cadoprograms.SM(poly=polyfilename, renumber=renumberfilename,
	      	    badidealinfo=badidealinfofilename,
                    purged=purgedfilename, index=indexfilename,
                    ell=gorder, smexp=smexp,
                    nmaps=nmaps,
                    out=smfilename,
                    stdout=str(stdoutpath),
                    stderr=str(stderrpath),
                    **self.progparams[0])
            message = self.submit_command(p, "", log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            if not smfilename.isfile():
                raise Exception("Output file %s does not exist" % smfilename)
            self.state["sm"] = smfilename.get_wdir_relative()
        self.logger.debug("Exit SMTask.run(" + self.name + ")")
        return True
    
    def get_sm_filename(self):
        return self.get_state_filename("sm")



class ReconstructLogTask(Task):
    """ Logarithms Reconstruction Task """
    @property
    def name(self):
        return "reconstructlog"
    @property
    def title(self):
        return "Logarithms Reconstruction"
    @property
    def programs(self):
        return (cadoprograms.ReconstructLog,)
    @property
    def progparam_override(self):
        return [["poly", "purged", "renumber", "dlog",  "ell", "smexp",
                 "nmaps", "partial", "ker", "ideals", "relsdel", "nrels"]]
    @property
    def paramnames(self):
        return super().paramnames

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
    
    def run(self):
        self.logger.debug("%s.run(): Task state: %s", self.__class__.name,
                          self.state)

        if not "dlog" in self.state:
            self.logger.info("Starting")
            polyfilename = self.send_request(Request.GET_POLYNOMIAL_FILENAME)
            renumberfilename = self.send_request(Request.GET_RENUMBER_FILENAME)
            purgedfilename = self.send_request(Request.GET_PURGED_FILENAME)
            kerfilename = self.send_request(Request.GET_KERNEL_FILENAME)
            idealfilename = self.send_request(Request.GET_IDEAL_FILENAME)
            relsdelfilename = self.send_request(Request.GET_RELSDEL_FILENAME)
            dlogfilename = self.workdir.make_filename("dlog")
            gorder = self.send_request(Request.GET_ELL)
            smexp = self.send_request(Request.GET_SMEXP)
            nmaps = self.send_request(Request.GET_NMAPS)

            nfree = self.send_request(Request.GET_FREEREL_RELCOUNT)
            nunique = self.send_request(Request.GET_UNIQUE_RELCOUNT)
            nrels = nfree+nunique
            
            (stdoutpath, stderrpath) = \
                    self.make_std_paths(cadoprograms.SM.name)
            p = cadoprograms.ReconstructLog(poly=polyfilename,
                    purged=purgedfilename,
                    renumber=renumberfilename,
                    dlog=dlogfilename,
                    ell=gorder, smexp=smexp,
                    nmaps=nmaps,
                    partial=True,
                    ker=kerfilename,
                    ideals=idealfilename,
                    relsdel=relsdelfilename,
                    nrels=nrels,
                    stdout=str(stdoutpath),
                    stderr=str(stderrpath),
                    **self.progparams[0])
            message = self.submit_command(p, "", log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            if not dlogfilename.isfile():
                raise Exception("Output file %s does not exist" % dlogfilename)
            self.state["dlog"] = dlogfilename.get_wdir_relative()
        self.logger.debug("Exit ReconstructLogTask.run(" + self.name + ")")
        return True
    
    def get_dlog_filename(self):
        return self.get_state_filename("dlog")
    
    def get_log2log3(self):
        filename = self.get_state_filename("dlog").get_wdir_relative()
        fullfile = self.params["workdir"].rstrip(os.sep) + os.sep + filename
        log2 = None
        log3 = None
        myfile = open(fullfile, "rb")
        data = myfile.read()
        for line in data.splitlines():
            match = re.match(br'(\d+) 2 0 rat (\d+)', line)
            if match:
                log2 = match.group(2)
            match = re.match(br'(\d+) 3 0 rat (\d+)', line)
            if match:
                log3 = match.group(2)
            if log2 != None and log3 != None:
                myfile.close()
                return [ log2, log3 ]
        raise Exception("Could not find log2 and log3 in %s" % filename)

class StartServerTask(DoesLogging, cadoparams.UseParameters, wudb.HasDbConnection):
    """ Starts HTTP server """
    @property
    def name(self):
        return "server"
    @property
    def title(self):
        return "Server Launcher"
    @property
    def paramnames(self):
        return {"name": str, "workdir": None, "address": None, "port": 8001,
                "threaded": False, "ssl": True, "whitelist": None,
                "only_registered": True}
    @property
    def param_nodename(self):
        return self.name

    # The whitelist parameter here is an iterable of strings in CIDR notation.
    # The whitelist parameter we get from params (i.e., from the parameter file)
    # is a string with comma-separated CIDR strings. The two are concatenated
    # to form the server whitelist.
    def __init__(self, *, default_workdir, parameters, path_prefix, db, whitelist=None):
        super().__init__(db=db, parameters=parameters, path_prefix=path_prefix)
        # self.logger.info("path_prefix = %s, parameters = %s", path_prefix, parameters)
        self.params = self.parameters.myparams(self.paramnames)
        serveraddress = self.params.get("address", None)
        serverport = self.params["port"]
        basedir = self.params.get("workdir", default_workdir).rstrip(os.sep) + os.sep
        uploaddir = basedir + self.params["name"] + ".upload/"
        threaded = self.params["threaded"]
        # By default, allow access only to files explicitly registered by tasks,
        # i.e., those files required by clients when downloading input files for
        # their workunits. By setting only_registered=False, access to all files
        # under the server working directory is allowed.
        only_registered = self.params["only_registered"]
        if self.params["ssl"]:
            cafilename = basedir + self.params["name"] + ".server.cert"
        else:
            cafilename = None

        server_whitelist = []
        if not whitelist is None:
            server_whitelist += whitelist
        if "whitelist" in self.params:
            server_whitelist += [h.strip() for h in self.params["whitelist"].split(",")]

        # If (1) any clients are to be started on localhost, but (2) the server
        # is listening on a network-visible address, then we need to whitelist
        # the network-visible address(es) of the current host as well, because
        # the client's connection will come from (one of) the network-visible
        # addresses of the current host.
        # For test (1), it should suffice to look for "localhost" or
        # "127.0.0.1" in the existing whitelist, as the host names on which to
        # start clients are inserted verbatim.
        # For (2), we check that serveraddress is either None (i.e., the
        # wildcard address which is network-visible), or anything other than
        # "localhost"
        if (serveraddress is None or socket.getfqdn(serveraddress) != "localhost") \
                and set(server_whitelist) & {"localhost", "127.0.0.1"}:
            hostname = socket.gethostname()
            if not hostname in server_whitelist:
                self.logger.info("Adding %s to whitelist to allow clients on localhost to connect", hostname)
                server_whitelist.append(hostname)
        if not server_whitelist:
            server_whitelist = None

        self.registered_filenames = self.make_db_dict('server_registered_filenames')
        self.server = wuserver.ServerLauncher(serveraddress, serverport,
            threaded, self.get_db_filename(), self.registered_filenames,
            uploaddir, bg=True, only_registered=only_registered, cafile=cafilename,
            whitelist=server_whitelist)

    def run(self):
        self.server.serve()

    def shutdown(self):
        self.server.shutdown()        

    def stop_serving_wus(self):
        self.server.stop_serving_wus()

    def get_url(self):
        return self.server.get_url()

    def get_cert_sha1(self):
        return self.server.get_cert_sha1()

    def register_filename(self, d):
        for key in d:
            if not key in self.registered_filenames:
                self.logger.debug("Registering file name %s with target %s",
                                  key, d[key])
                self.registered_filenames[key] = d[key]
            elif d[key] != self.registered_filenames[key]:
                # It was already registered with a different target. This will
                # happen if, e.g., the user chooses a different build directory
                # between runs. It's still fragile, as the server will try to
                # serve the old file(s) until a new workunit is generated and
                # overrides the target.
                # The proper solution would be to make Program classes
                # Templates, so Tasks can instantiate them at __init__() and
                # register the resolved paths once. The registered_filenames
                # dict could then be memory-backed. Tasks would also have to
                # register their own files which may need to be served in
                # __init__().
                self.logger.warning("Filename %s, to be registered for "
                                    "target %s, was previously registered for "
                                    "target %s. Overriding with new target.",
                                    key, d[key], self.registered_filenames[key])
                self.registered_filenames[key] = d[key]
            else:
                # Was already registered with the same target. Nothing to do
                pass

class StartClientsTask(Task):
    """ Starts clients on slave machines """
    @property
    def name(self):
        return "slaves"
    @property
    def title(self):
        return "Client Launcher"
    @property
    def programs(self):
        return (cadoprograms.WuClient,)
    @property
    def progparam_override(self):
        return [["clientid", "certsha1"]]
    @property
    def paramnames(self):
        return {'hostnames': str, 'scriptpath': None, "nrclients": [int], "run": True}
    @property
    def param_nodename(self):
        return None
    
    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.used_ids = {}
        self.pids = self.make_db_dict(self.make_tablename("client_pids"), connection=self.db_connection)
        self.hosts = self.make_db_dict(self.make_tablename("client_hosts"), connection=self.db_connection)
        assert set(self.pids) == set(self.hosts)
        # Invariants: the keys of self.pids and of self.hosts are the same set.
        # The keys of self.used_ids are a subset of the keys of self.pids.
        # A clientid is in self.used_ids if we know that clientid to be
        # currently running.
        
        if 'scriptpath' in self.params:
            self.progparams[0]['execpath'] = self.params['scriptpath']
        
        # If hostnames are of the form @file, read host names from file,
        # one host name per line
        match = re.match(r"@(.*)", self.params["hostnames"])
        if match:
            with open(match.group(1)) as f:
                self.hosts_to_launch = [line.strip() for line in f]
        else:
            self.hosts_to_launch = [host.strip() for host in
                    self.params["hostnames"].split(",")]

        if "nrclients" in self.params:
            self.hosts_to_launch = self.make_multiplicity(self.hosts_to_launch,
                   self.params["nrclients"])

    @staticmethod
    def make_multiplicity(names, multi):
        """ Produce a list in which each unique entry of the list "names"
        occurs "multi" times. The order of elements in names is preserved.
        
        >>> names = ['a', 'b', 'a', 'c', 'c', 'a', 'a']
        >>> StartClientsTask.make_multiplicity(names, 1)
        ['a', 'b', 'c']
        >>> StartClientsTask.make_multiplicity(names, 2)
        ['a', 'a', 'b', 'b', 'c', 'c']
        """
        result = []
        # Use OrderedDict to get unique names, preserving order
        for name in OrderedDict.fromkeys(names, None):
            result.extend([name] * multi)
        return result

    def get_hosts_to_launch(self):
        """ Get list host names on which clients should run """
        return self.hosts_to_launch

    def is_alive(self, clientid):
        # Simplistic: just test if process with that pid exists and accepts
        # signals from us. TODO: better testing here, probably with ps|grep
        # or some such
        (rc, stdout, stderr) = self.kill_client(clientid, signal=0)
        return (rc == 0)
    
    def _add_cid(self, clientid, pid, host):
        """ Add a client id atomically to both the "pids" and "hosts"
        dictionaries
        """
        self.pids.update({clientid: pid}, commit=False)
        self.hosts.update({clientid: host}, commit=True)

    def _del_cid(self, clientid):
        """ Remove a client id atomically from both the "pids" and "hosts"
        dictionaries
        """
        self.pids.clear([clientid], commit=False)
        self.hosts.clear([clientid], commit=True)
    
    def launch_clients(self, server, certsha1=None):
        for host in self.hosts_to_launch:
            self.launch_one_client(host.strip(), server, certsha1=certsha1)
        running_clients = [(cid, self.hosts[cid], pid) for (cid, pid) in
            self.pids.items()]
        s = ", ".join(["%s (Host %s, PID %d)" % t for t in running_clients])
        self.logger.info("Running clients: %s" % s)
        # Check for old clients which we did not mean to start this run
        for cid in set(self.pids) - set(self.used_ids):
            if self.is_alive(cid):
                self.logger.warn("Client id %s (Host %s, PID %d), launched "
                                 "in a previous run and not meant to be "
                                 "launched this time, is still running",
                                 cid, self.hosts[cid], self.pids[cid])
            else:
                self.logger.warn("Client id %s (Host %s, PID %d), launched "
                                 "in a previous run and not meant to be "
                                 "launched this time, seems to have died. "
                                 "I'll forget about this client.",
                                 cid, self.hosts[cid], self.pids[cid])
                self._del_cid(cid)
    
    def make_unique_id(self, host):
        # Make a unique client id for host
        clientid = host
        i = 1
        while clientid in self.used_ids:
            assert clientid in self.pids
            assert clientid in self.hosts
            i += 1
            clientid = "%s+%d" % (host, i)
        return clientid
    
    # Cases:
    # Client was never started. Start it, add to state
    # Client was started, but does not exist any more. Remove from state,
    #   then start and add again
    # Client was started, and does still exists. Nothing to do.
    
    def launch_one_client(self, host, server, *, clientid=None, certsha1=None):
        if clientid is None:
            clientid = self.make_unique_id(host)
        # Check if client is already running
        if clientid in self.pids:
            assert self.hosts[clientid] == host
            if self.is_alive(clientid):
                self.logger.info("Client %s on host %s with PID %d already "
                                 "running",
                                 clientid, host, self.pids[clientid])
                self.used_ids[clientid] = True
                return
            else:
                self.logger.info("Client %s on host %s with PID %d seems to have died",
                                 clientid, host, self.pids[clientid])
                self._del_cid(clientid)
        
        self.logger.info("Starting client id %s on host %s", clientid, host)
        wuclient = cadoprograms.WuClient(server=server,
                                         clientid=clientid, daemon=True,
                                         certsha1=certsha1,
                                         **self.progparams[0])
        if host == "localhost":
            process = cadocommand.Command(wuclient)
        else:
            process = cadocommand.RemoteCommand(wuclient, host, self.parameters)
        (rc, stdout, stderr) = process.wait()
        if rc != 0:
            self.logger.warning("Starting client on host %s failed.", host)
            if stdout:
                self.logger.warning("Stdout: %s", stdout.decode("utf-8").strip())
            if stderr:
                self.logger.warning("Stderr: %s", stderr.decode("utf-8").strip())
            return
        match = None
        if not stdout is None:
            match = re.match(r"PID: (\d+)", stdout.decode("utf-8"))
        if not match:
            self.logger.warning("Client did not print PID")
            if not stdout is None:
                self.logger.warning("Stdout: %s", stdout.decode("utf-8").strip())
            if not stderr is None:
                self.logger.warning("Stderr: %s", stderr.decode("utf-8").strip())
            return
        self.used_ids[clientid] = True
        self._add_cid(clientid, int(match.group(1)), host)

    def kill_all_clients(self):
        # Need the list() to make a copy as dict will change in loop body
        for clientid in list(self.pids):
            (rc, stdout, stderr) = self.kill_client(clientid)
            if rc == 0:
                self.logger.info("Stopped client %s (Host %s, PID %d)",
                                 clientid, self.hosts[clientid], self.pids[clientid])
                self._del_cid(clientid)
            else:
                self.logger.warning("Stopping client %s (Host %s, PID %d) failed",
                                    clientid, self.hosts[clientid], self.pids[clientid])
                if stdout:
                    self.logger.warning("Stdout: %s", stdout.decode("utf-8").strip())
                if stderr:
                    self.logger.warning("Stderr: %s", stderr.decode("utf-8").strip())
                # Assume that the client is already dead and remove it from
                # the list of running clients
                self._del_cid(clientid)
    
    def kill_client(self, clientid, signal=None):
        pid = self.pids[clientid]
        host = self.hosts[clientid]
        kill = cadoprograms.Kill(pid, signal=signal)
        if host == "localhost":
            process = cadocommand.Command(kill)
        else:
            process = cadocommand.RemoteCommand(kill, host, self.parameters)
        return process.wait()

class Message(object):
    def __init__(self, sender, key, value=None):
        self.sender = sender
        self.key = key
        self.value = value
    def get_sender(self):
        return self.sender
    def get_key(self):
        return self.key
    def get_value(self):
        return self.value
    @classmethod
    def reverse_lookup(cls, reference):
        for key in dir(cls):
            if getattr(cls, key) == reference:
                return key


class Notification(Message):
    FINISHED_POLYNOMIAL_SELECTION = object()
    WANT_MORE_RELATIONS = object()
    HAVE_ENOUGH_RELATIONS = object()
    REGISTER_FILENAME = object()
    UNREGISTER_FILENAME = object()
    WANT_TO_RUN = object()
    SUBSCRIBE_WU_NOTIFICATIONS = object()

class Request(Message):
    # Lacking a proper enum before Python 3.4, we generate dummy objects
    # which have separate identity and can be used as dict keys
    GET_RAW_POLYNOMIALS = object()
    GET_POLYNOMIAL = object()
    GET_POLYNOMIAL_FILENAME = object()
    GET_POLY_RANK = object()
    GET_FACTORBASE_FILENAME = object()
    GET_FREEREL_FILENAME = object()
    GET_RENUMBER_FILENAME = object()
    GET_FREEREL_RELCOUNT = object()
    GET_RENUMBER_PRIMECOUNT = object()
    GET_SIEVER_FILENAMES = object()
    GET_SIEVER_RELCOUNT = object()
    GET_DUP1_FILENAMES = object()
    GET_DUP1_RELCOUNT = object()
    GET_UNIQUE_RELCOUNT = object()
    GET_UNIQUE_FILENAMES = object()
    GET_PURGED_FILENAME = object()
    GET_MERGED_FILENAME = object()
    GET_INDEX_FILENAME = object()
    GET_IDEAL_FILENAME = object()
    GET_DENSE_FILENAME = object()
    GET_DEPENDENCY_FILENAME = object()
    GET_LINALG_PREFIX = object()
    GET_KERNEL_FILENAME = object()
    GET_RELSDEL_FILENAME = object()
    GET_SM_FILENAME = object()
    GET_BADIDEAL_FILENAME = object()
    GET_BADIDEALINFO_FILENAME = object()
    GET_SMEXP = object()
    GET_ELL = object()
    GET_NMAPS = object()
    GET_WU_RESULT = object()

class CompleteFactorization(HasState, wudb.DbAccess, 
        DoesLogging, cadoparams.UseParameters, patterns.Mediator):
    """ The complete factorization, aggregate of the individual tasks """
    @property
    def name(self):
        return "tasks"
    @property
    def param_nodename(self):
        return self.name
    @property
    def paramnames(self):
        # This isn't a Task subclass so we don't really need to define
        # paramnames, but we do it out of habit
        return {"name": str, "workdir": str, "N": int, "dlp": False,
                "trybadwu": False}
    @property
    def title(self):
        return "Complete Factorization"
    @property
    def programs(self):
        return []
    
    def __init__(self, db, parameters, path_prefix):
        super().__init__(db=db, parameters=parameters, path_prefix=path_prefix)
        self.params = self.parameters.myparams(self.paramnames)
        self.db_listener = self.make_db_listener()
        
        # Init WU BD
        self.wuar = self.make_wu_access()
        self.wuar.create_tables()
        if self.params["trybadwu"]:
            # Test behaviour when a WU is in the DB that does not belong to
            # any task. It should get cancelled with an error message.
            self.wuar.create(["WORKUNIT FAKE_WU_%s\nCOMMAND true" % time.time()])

        # Start with an empty list of tasks that want to run. Tasks will add
        # themselves during __init__().
        self.tasks_that_want_to_run = list()

        # Init client lists
        self.clients = []
        whitelist = set()
        for (path, key) in self.parameters.get_parameters().find(['slaves'], 'hostnames'):
            self.clients.append(StartClientsTask(mediator=self,
                                                 db=db,
                                                 parameters=self.parameters,
                                                 path_prefix=path))
            hostnames = self.clients[-1].get_hosts_to_launch()
            whitelist |= set(hostnames)

        whitelist = list(whitelist) if whitelist else None
        # Init server task
        self.servertask = StartServerTask(default_workdir=self.params["workdir"],
                parameters=parameters, path_prefix=path_prefix, db=db,
                whitelist=whitelist)

        parampath = self.parameters.get_param_path()
        polyselpath = parampath + ['polyselect']
        sievepath = parampath + ['sieve']
        filterpath = parampath + ['filter']
        linalgpath = parampath + ['linalg']
        
        ## tasks that are common to factorization and dlp
        self.polysel1 = Polysel1Task(mediator=self,
                                   db=db,
                                   parameters=self.parameters,
                                   path_prefix=polyselpath)
        self.polysel2 = Polysel2Task(mediator=self,
                                   db=db,
                                   parameters=self.parameters,
                                   path_prefix=polyselpath)
        self.fb = FactorBaseTask(mediator=self,
                                 db=db,
                                 parameters=self.parameters,
                                 path_prefix=sievepath)
        self.freerel = FreeRelTask(mediator=self,
                                   db=db,
                                   parameters=self.parameters,
                                   path_prefix=sievepath)
        self.sieving = SievingTask(mediator=self,
                                   db=db,
                                   parameters=self.parameters,
                                   path_prefix=sievepath)
        self.dup1 = Duplicates1Task(mediator=self,
                                    db=db,
                                    parameters=self.parameters,
                                    path_prefix=filterpath)
        self.dup2 = Duplicates2Task(mediator=self,
                                    db=db,
                                    parameters=self.parameters,
                                    path_prefix=filterpath)
        self.purge = PurgeTask(mediator=self,
                               db=db,
                               parameters=self.parameters,
                               path_prefix=filterpath)
        if self.params["dlp"]:
            ## Tasks specific to dlp
            self.nmbrthry = NmbrthryTask(mediator=self,
                             db=db,
                             parameters=self.parameters,
                             path_prefix=parampath)
            self.sm = SMTask(mediator=self,
                             db=db,
                             parameters=self.parameters,
                             path_prefix=filterpath)
            self.merge = MergeDLPTask(mediator=self,
                                   db=db,
                                   parameters=self.parameters,
                                   path_prefix=filterpath)
            self.linalg = LinAlgDLPTask(mediator=self,
                                     db=db,
                                     parameters=self.parameters,
                                     path_prefix=linalgpath)
            self.reconstructlog = ReconstructLogTask(mediator=self,
                                     db=db,
                                     parameters=self.parameters,
                                     path_prefix=filterpath)
        else:
            ## Tasks specific to factorization
            self.merge = MergeTask(mediator=self,
                                   db=db,
                                   parameters=self.parameters,
                                   path_prefix=filterpath)
            self.linalg = LinAlgTask(mediator=self,
                                     db=db,
                                     parameters=self.parameters,
                                     path_prefix=linalgpath)
            self.characters = CharactersTask(mediator=self,
                                             db=db,
                                             parameters=self.parameters,
                                             path_prefix=linalgpath)
            self.sqrt = SqrtTask(mediator=self,
                                 db=db,
                                 parameters=self.parameters,
                                 path_prefix=parampath)
        
        # Defines an order on tasks in which tasks that want to run should be
        # run
        if self.params["dlp"]:
            self.tasks = (self.polysel1, self.polysel2, self.nmbrthry, self.fb,
                          self.freerel, self.sieving,
                          self.dup1, self.dup2, self.purge, self.merge,
                          self.sm, self.linalg, self.reconstructlog)
        else:
            self.tasks = (self.polysel1, self.polysel2, self.fb, self.freerel, self.sieving,
                          self.dup1, self.dup2, self.purge, self.merge,
                          self.linalg, self.characters, self.sqrt)

        for (path, key, value) in parameters.get_unused_parameters():
            self.logger.warning("Parameter %s = %s was not used anywhere",
                                ".".join(path + [key]), value)

        self.request_map = {
            Request.GET_RAW_POLYNOMIALS: self.polysel1.get_raw_polynomials,
            Request.GET_POLY_RANK: self.polysel1.get_poly_rank,
            Request.GET_POLYNOMIAL: self.polysel2.get_poly,
            Request.GET_POLYNOMIAL_FILENAME: self.polysel2.get_poly_filename,
            Request.GET_FACTORBASE_FILENAME: self.fb.get_filename,
            Request.GET_FREEREL_FILENAME: self.freerel.get_freerel_filename,
            Request.GET_RENUMBER_FILENAME: self.freerel.get_renumber_filename,
            Request.GET_FREEREL_RELCOUNT: self.freerel.get_nrels,
            Request.GET_RENUMBER_PRIMECOUNT: self.freerel.get_nprimes,
            Request.GET_SIEVER_FILENAMES: self.sieving.get_output_filenames,
            Request.GET_SIEVER_RELCOUNT: self.sieving.get_nrels,
            Request.GET_DUP1_FILENAMES: self.dup1.get_output_filenames,
            Request.GET_DUP1_RELCOUNT: self.dup1.get_nrels,
            Request.GET_UNIQUE_RELCOUNT: self.dup2.get_nrels,
            Request.GET_UNIQUE_FILENAMES: self.dup2.get_output_filenames,
            Request.GET_PURGED_FILENAME: self.purge.get_purged_filename,
            Request.GET_MERGED_FILENAME: self.merge.get_merged_filename,
            Request.GET_INDEX_FILENAME: self.merge.get_index_filename,
            Request.GET_DENSE_FILENAME: self.merge.get_dense_filename,
            Request.GET_WU_RESULT: self.db_listener.send_result
        }
        ## add requests specific to dlp or factoring
        if self.params["dlp"]:
            self.request_map[Request.GET_IDEAL_FILENAME] = self.merge.get_ideal_filename
            self.request_map[Request.GET_BADIDEAL_FILENAME] = self.nmbrthry.get_bad_filename
            self.request_map[Request.GET_BADIDEALINFO_FILENAME] = self.nmbrthry.get_badinfo_filename
            self.request_map[Request.GET_NMAPS] = self.nmbrthry.get_nmaps
            self.request_map[Request.GET_SMEXP] = self.nmbrthry.get_smexp
            self.request_map[Request.GET_ELL] = self.nmbrthry.get_ell
            self.request_map[Request.GET_SM_FILENAME] = self.sm.get_sm_filename
            self.request_map[Request.GET_KERNEL_FILENAME] = self.linalg.get_kernel_filename
            self.request_map[Request.GET_RELSDEL_FILENAME] = self.purge.get_relsdel_filename
        else:
            self.request_map[Request.GET_KERNEL_FILENAME] = self.characters.get_kernel_filename
            self.request_map[Request.GET_DEPENDENCY_FILENAME] = self.linalg.get_dependency_filename
            self.request_map[Request.GET_LINALG_PREFIX] = self.linalg.get_prefix

    def run(self):
        had_interrupt = False
        if self.params["dlp"]:
            self.logger.info("Computing Discrete Logs in GF(%s)", self.params["N"])
        else:
            self.logger.info("Factoring %s", self.params["N"])
        self.start_elapsed_time()

        self.servertask.run()
        
        last_task = None
        last_status = True
        try:
            self.start_all_clients()
            
            while last_status:
                last_status, last_task = self.run_next_task()
            
            for task in self.tasks:
                task.print_stats()
        
        except KeyboardInterrupt:
            self.logger.fatal("Received KeyboardInterrupt. Terminating")
            had_interrupt = True

        self.stop_all_clients()
        self.servertask.shutdown()
        elapsed = self.end_elapsed_time()

        if had_interrupt:
            return None

        cputotal = self.get_sum_of_cpu_or_real_time(True)
        # Do we want the sum of real times over all sub-processes for
        # something?
        # realtotal = self.get_sum_of_cpu_or_real_time(False)
        self.logger.info("Total cpu/elapsed time for entire factorization: %g/%g",
                         cputotal, elapsed)

        if last_task and not last_status:
            self.logger.fatal("Premature exit within %s. Bye.", last_task)
            return None

        if self.params["dlp"]:
            return [ self.params["N"], self.nmbrthry.get_ell() ] + self.reconstructlog.get_log2log3()
        else:
            return self.sqrt.get_factors()
    
    def start_all_clients(self):
        url = self.servertask.get_url()
        certsha1 = self.servertask.get_cert_sha1()
        for clients in self.clients:
            clients.launch_clients(url, certsha1)
    
    def stop_all_clients(self):
        for clients in self.clients:
            clients.kill_all_clients()

    def start_elapsed_time(self):
        if "starttime" in self.state:
            self.logger.warning("The start time of the last cadofactor.py "
                                "run was recorded, but not its end time, "
                                "maybe because it died unexpectedly.")
            self.logger.warning("Elapsed time of last run is not known and "
                                "will not be counted towards total.")
        self.state["starttime"] = time.time()
        
    def end_elapsed_time(self):
        if not "starttime" in self.state:
            self.logger.error("Missing starttime in end_elapsed_time(). "
                              "This should not have happened.")
            return
        elapsed = time.time() - self.state["starttime"]
        elapsed += self.state.get("elapsed", 0)
        self.state.__delitem__("starttime", commit=False)
        self.state.update({"elapsed": elapsed}, commit=True)
        return elapsed
    
    def run_next_task(self):
        for task in self.tasks:
            if task in self.tasks_that_want_to_run:
                # self.logger.info("Next task that wants to run: %s",
                #                  task.title)
                self.tasks_that_want_to_run.remove(task)
                return [task.run(), task.title]
        return [False, None]
    
    def get_sum_of_cpu_or_real_time(self, is_cpu):
        total = 0
        for task in self.tasks:
            task_time = task.get_total_cpu_or_real_time(is_cpu)
            total += task_time
            # self.logger.info("Task %s reports %s time of %g, new total: %g",
            #         task.name, "cpu" if is_cpu else "real", task_time, total)
        return total

    def register_filename(self, d):
        return self.servertask.register_filename(d)
    
    def relay_notification(self, message):
        """ The relay for letting Tasks talk to us and each other """
        assert isinstance(message, Notification)
        sender = message.get_sender()
        key = message.get_key()
        value = message.get_value()
        self.logger.message("Received notification from %s, key = %s, value = %s",
                            sender, Notification.reverse_lookup(key), value)
        if key is Notification.WANT_MORE_RELATIONS:
            if sender is self.purge:
                self.dup2.request_more_relations(value)
            elif sender is self.dup2:
                self.dup1.request_more_relations(value)
            elif sender is self.dup1:
                self.sieving.request_more_relations(value)
            else:
                raise Exception("Got WANT_MORE_RELATIONS from unknown sender")
        elif key is Notification.HAVE_ENOUGH_RELATIONS:
            if sender is self.purge:
                self.servertask.stop_serving_wus()
                self.sieving.cancel_available_wus()
                self.stop_all_clients()
            else:
                raise Exception("Got HAVE_ENOUGH_RELATIONS from unknown sender")
        elif key is Notification.REGISTER_FILENAME:
            if isinstance(sender, ClientServerTask):
                self.register_filename(value)
            else:
                raise Exception("Got REGISTER_FILENAME, but not from a ClientServerTask")
        elif key is Notification.WANT_TO_RUN:
            if sender in self.tasks_that_want_to_run:
                raise Exception("Got request from %s to run, but it was in run queue already",
                                sender)
            else:
                self.tasks_that_want_to_run.append(sender)
        elif key is Notification.SUBSCRIBE_WU_NOTIFICATIONS:
            return self.db_listener.subscribeObserver(sender)
        else:
            raise KeyError("Notification from %s has unknown key %s" % (sender, key))
    
    def answer_request(self, request):
        assert isinstance(request, Request)
        sender = request.get_sender()
        key = request.get_key()
        value = request.get_value()
        self.logger.message("Received request from %s, key = %s, values = %s",
                            sender, Request.reverse_lookup(key), value)
        if not key in self.request_map:
            raise KeyError("Unknown Request key %s from sender %s" %
                           (key, sender))
        if value is None:
            result = self.request_map[key]()
        else:
            result = self.request_map[key](value)
        self.logger.message("Completed request from %s, key = %s, values = %s, result = %s",
                            sender, Request.reverse_lookup(key), value, result)
        return result
    
    def handle_message(self, message):
        if isinstance(message, Notification):
            self.relay_notification(Notification)
        elif isinstance(message, Request):
            return self.answer_request(message)
        else:
            raise TypeError("Message is neither Notification nor Request")
