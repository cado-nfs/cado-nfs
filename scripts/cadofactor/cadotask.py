#!/usr/bin/env python3

import re
import os
from math import gcd
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
import errno
from cadofactor import patterns, wudb, cadoprograms
from cadofactor import cadoparams, cadocommand, workunit
from cadofactor.workunit import Workunit
from struct import error as structerror
from shutil import rmtree
from cadofactor.api_server import ApiServer
from cadofactor.cadoutils import Algorithm, Computation
from cadofactor.cadoutils import xgcd, CRT, primes_above

# Pattern for floating-point numbers
RE_FP = r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?"


def re_fp_compile(pattern, *args, **kwargs):
    """
    This is just like re.compile(), except that any string "{fp}" within the
    pattern is subsituted by a regular expression that matches floating-point
    numbers (incl. in scientific notation).
    "{a}" or "{a,b}" with a, b numeric are repetition operators in regular
    expressions, but if a or b are not numeric text, they match just as literal
    text, so using "{fp}" shouldn't introduce conflicts with valid regexes.
    """

    # We can't use str.format() here or it will try to interpret regex
    # repetition operators {n} or {n,m} and complain
    fp_pattern = pattern.replace("{fp}", RE_FP)
    return re.compile(fp_pattern, *args, **kwargs)


# convert the time s in [day,]hours,minutes,seconds
# (patch from Hermann Stamm-Wilbrandt)
def tstr(s):
    i = round(s)
    if i < 3600 or i >= 366 * 24 * 3600:
        return ""
    fmt = ""
    if i >= 24 * 3600:
        fmt = "%-jd "
        i -= 24 * 3600
    return " [" + time.strftime(fmt + "%H:%M:%S", time.gmtime(i)) + "]"


class Polynomial(list):
    r"""
    >>> # This docstring needs to be a raw string (r prefix) because
    >>> # backslashes occur in it
    >>> p = Polynomial()
    >>> p.degree < 0
    True
    >>> p[0] = 1
    >>> p.degree == 0
    True
    >>> p[42] = 1
    >>> p.degree == 42
    True
    >>> p[42] = 0
    >>> p.degree == 0
    True
    >>> p = Polynomial([1,-5,0,0,0,7])
    >>> str(p)
    '7*x^5-5*x+1'
    >>> p.as_lines("c")
    ['c0: 1\n', 'c1: -5\n', 'c5: 7\n']
    >>> p = Polynomial([3,2,1]) # x^2 + 2*x + 3
    >>> p.eval(0)
    3
    >>> p.eval(1)
    6
    >>> p.eval(2)
    11
    >>> p.eval(-3)
    6
    >>> p.eval_h(2,7)
    179
    >>> p.eval_h(-3,5)
    54
    """
    @property
    def degree(self):
        return len(self) - 1 if len(self) > 0 else float("-inf")

    def __setitem__(self, index, value):
        if index >= len(self):
            self.extend([0]*(index + 1 - len(self)))
        list.__setitem__(self, index, value)
        # Remove leading zeroes
        while len(self) > 0 and self[-1] == 0:
            self.pop()

    def __str__(self):
        xpow = ["", "*x"] + ["*x^%d" % i for i in range(2, len(self))]
        arr = ["%+d%s" % (self[idx], xpow[idx]) for idx in range(0, len(self))
               if self[idx]]
        poly = "".join(reversed(arr)).lstrip('+')
        poly = re.sub(r'\b1\*', "", poly)
        return poly

    # Returns the coefficients as an array of strings in a format like we use
    # in polynomial files. See doctest example.
    # Coefficients which are 0 get omitted, each line has a trailing newline
    def as_lines(self, prefix):
        return ["%s%d: %d\n" % (prefix, idx, coeff)
                for (idx, coeff) in enumerate(self)
                if not coeff == 0]

    def eval(self, x):
        """
        Evaluate the polynomial at x
        """
        if len(self) == 0:
            return 0
        deg = self.degree
        value = self[deg]
        for i in range(deg):
            value = value * x + self[deg - i - 1]
        return value

    def eval_h(self, a, b):
        """
        Evaluate homogenized bi-variate polynomial at a,b
        """
        if len(self) == 0:
            return 0
        powers_a = [a**i for i in range(self.degree + 1)]
        powers_b = [b**i for i in range(self.degree + 1)]
        return sum([coeff * pow_a * pow_b for (coeff, pow_a, pow_b)
                    in zip(self, powers_a, reversed(powers_b))])

    def same_lc(self, other):
        """
        Return true if the two polynomials have the same degree
        and leading coefficient
        """
        return self.degree == other.degree and \
            self[self.degree] == other[other.degree]


class PolynomialParseException(Exception):
    """
    Exception class for signaling errors during polynomial parsing
    """
    pass


class TaskException(Exception):
    """
    Exception class for signaling errors during task execution
    """
    pass


class EarlyStopException(Exception):
    """
    Exception class for cases like tasks.sieve.run=false
    """
    pass


class Polynomials(object):
    r"""
    A class that represents a polynomial

    >>> Polynomials([""])
    Traceback (most recent call last):
    cadotask.PolynomialParseException: No polynomials found
    >>> t="n: 1021\nc0: 1\nc5: -1\nc5: 1\nY0: 4\nY1: -1\nskew: 1.0\n"
    >>> p=Polynomials(t.splitlines())  # doctest: +NORMALIZE_WHITESPACE
    Traceback (most recent call last):
    cadotask.PolynomialParseException:
    Line 'c5: 1' redefines coefficient of x^5
    >>> t="n: 1021\n"
    >>> p=Polynomials(t.splitlines())  # doctest: +NORMALIZE_WHITESPACE
    Traceback (most recent call last):
    cadotask.PolynomialParseException:
    No polynomial f specified (c: lines, or poly1)
    >>> t="n: 1021\nc0: 1\nc1: -1\nc5: 1\nY0: 4\nY1: -1\nskew: 1.0\n"
    >>> p=Polynomials(t.splitlines())
    >>> str(p)                         # doctest: +NORMALIZE_WHITESPACE
    'n: 1021\nskew: 1.0\nc0: 1\nc1: -1\nc5: 1\nY0:
    4\nY1: -1\n# f(x) = x^5-x+1\n# g(x) = -x+4\n'
    >>> t="n: 1021\nc0: -1\nc1: 1\nc5: -1\nY0: -4\nY1: 1\nskew: 1.0\n"
    >>> p=Polynomials(t.splitlines())
    >>> str(p)                         # doctest: +NORMALIZE_WHITESPACE
    'n: 1021\nskew: 1.0\nc0: -1\nc1: 1\nc5: -1\nY0:
    -4\nY1: 1\n# f(x) = -x^5+x-1\n# g(x) = x-4\n'
    >>> # Without skew
    >>> t="n: 1021\nc0: 1\nc1: -1\nc5: 1\nY0: 4\nY1: -1\n\n"
    >>> p1=Polynomials(t.splitlines())
    >>> str(p1)                        # doctest: +NORMALIZE_WHITESPACE
    'n: 1021\nc0: 1\nc1: -1\nc5: 1\nY0:
    4\nY1: -1\n# f(x) = x^5-x+1\n# g(x) = -x+4\n'
    >>> # With all optional lines
    >>> t="n: 1021\nc0: 1\nc1: -1\nc5: 1\n" \
    ... "Y0: 4\nY1: -1\nskew: 1.0\ntype: gnfs\nm: 123\n"
    >>> p2=Polynomials(t.splitlines())
    >>> str(p2)                        # doctest: +NORMALIZE_WHITESPACE
    'n: 1021\nskew: 1.0\ntype: gnfs\nc0: 1\nc1: -1\nc5: 1\nY0:
    4\nY1: -1\n# f(x) = x^5-x+1\n# g(x) = -x+4\n'
    >>> p1.same_lc(p2)
    True
    >>> t="n: 1021\npoly0: 1, 2, 3\npoly1: 4, 5, 6\nskew: 1.0\n"
    >>> p=Polynomials(t.splitlines())
    >>> str(p)                         # doctest: +NORMALIZE_WHITESPACE
    'n: 1021\nskew: 1.0\npoly0: 1,2,3\npoly1: 4,5,6\n#
    poly0 = 3*x^2+2*x+1\n# poly1 = 6*x^2+5*x+4\n'
    """

    re_pol_f = re.compile(r"c(\d+)\s*:\s*(-?\d+)")
    re_pol_g = re.compile(r"Y(\d+)\s*:\s*(-?\d+)")
    re_polys = re.compile(r"poly(\d+)\s*:")  # FIXME: do better?
    re_Murphy = re_fp_compile(r"\s*#\s*(side\s+\d+\s+)?"
                              r"MurphyE\s*\((.*)\)\s*=\s*({fp})")
    # MurphyF is the refined value of MurphyE produced by polyselect3
    re_MurphyF = re_fp_compile(r"\s*#\s*MurphyF\s*\((.*)\)\s*=\s*({fp})")
    re_n = re.compile(r"n\s*:\s* (-?\d+)")  # Ex. "n: 1234567"
    re_skew = re_fp_compile(r"skew:\s*({fp})")  # Ex. "skew: 1.3e5"
    re_type = re.compile(r"type\s*:\s*(snfs|gnfs)")  # Ex. "type: snfs"
    # Note, the value of m is ignored by CADO-NFS, but should not trigger an
    # error if it occurs in a polynomial file.
    re_m = re.compile(r"m\s*:\s*(\d+)")  # Ex. "m: 123"
    re_best = re.compile(r"# Best polynomial found \(revision (.*)\):")
    re_exp_E = re_fp_compile(r"\s*#\s*exp_E\s*({fp})")

    def __init__(self, lines, allow_only_one_poly=False):
        """
        Parse a polynomial file in the syntax as produced by polyselect
        and polyselect_ropt
        """
        self.MurphyE = 0.
        self.MurphyF = 0.
        self.n = None
        self.skew = None
        self.type = None
        self.MurphyParams = None
        self.revision = None
        self.exp_E = 0.
        polyf = Polynomial()
        polyg = Polynomial()
        # in case of multiple fields
        tabpoly = {}

        def match_poly(line, poly, regex):
            match = regex.match(line)
            if match:
                (idx, coeff) = map(int, match.groups())
                if idx <= poly.degree and poly[idx]:
                    raise PolynomialParseException(
                        "Line '%s' redefines coefficient of x^%d"
                        % (line, idx))
                poly[idx] = coeff
                return True
            return False

        def match_poly_all(line, regex):
            # line = "poly0: 1, 2, 3" => poly[0] = [1, 2, 3] = 1+2*X+3*X^2
            match = regex.match(line)
            if match:
                line2 = line.split(":")
                # get index of poly
                ip = int(line2[0].split("poly")[1])
                # get coeffs of 1+2*X+3*X^2
                line3 = line2[1].split(",")
                pol = Polynomial()
                for idx in range(len(line3)):
                    pol[idx] = int(line3[idx])
                return ip, pol
            return -1, []

        for line in lines:
            # print ("Parsing line: >%s<" % line.strip())

            # First, match comment lines that contain useful info

            # If this is a comment line telling the Murphy E value,
            # extract the value and store it
            match = self.re_Murphy.match(line)
            if match:
                if self.MurphyParams or self.MurphyE:
                    raise PolynomialParseException(
                        "Line '%s' redefines Murphy E value" % line)
                self.MurphyParams = match.group(2)
                self.MurphyE = float(match.group(3))
                continue
            match = self.re_MurphyF.match(line)
            if match:
                if self.MurphyF:
                    raise PolynomialParseException(
                        "Line '%s' redefines Murphy F value" % line)
                self.MurphyF = float(match.group(2))
                continue
            match = self.re_best.match(line)
            if match:
                self.revision = match.group(1)
                continue
            # If this is a comment line telling the expected E-value,
            # extract the value and store it
            match = self.re_exp_E.match(line)
            if match:
                if self.exp_E != 0:
                    raise PolynomialParseException(
                        "Line '%s' redefines exp_E value" % line)
                self.exp_E = float(match.group(1))
                continue

            # If it's a comment that doesn't contain interesting data,
            # drop it

            line2 = line.split('#', 1)[0].strip()
            # If nothing is left, process next line
            if not line2:
                continue

            # Try to parse polynomial coefficients
            if match_poly(line, polyf, self.re_pol_f) or \
                    match_poly(line, polyg, self.re_pol_g):
                continue
            # is it in format "poly*: ..."
            ip, tip = match_poly_all(line, self.re_polys)
            if ip != -1:
                tabpoly[ip] = tip
                continue

            match = self.re_n.match(line)
            if match:
                if self.n is not None:
                    raise PolynomialParseException(
                        "Value of n redefined in line %s" % line)
                self.n = int(match.group(1))
                continue

            match = self.re_skew.match(line)
            if match:
                if self.skew is not None:
                    raise PolynomialParseException(
                        "Value of skewness redefined in line %s" % line)
                self.skew = float(match.group(1))
                continue

            match = self.re_type.match(line)
            if match:
                if self.type is not None:
                    raise PolynomialParseException(
                        "Type of factorization redefined in line %s" % line)
                self.type = match.group(1)
                continue

            match = self.re_m.match(line)
            if match:
                # Simply ignore "m:" lines
                continue

            # If nothing matches, complain
            raise PolynomialParseException("Invalid line '%s'" % line)

        # If no polynomial was found at all (not even partial data), assume
        # that polyselect simply did not find anything in this search range
        if polyf.degree < 0 and \
                polyg.degree < 0 and \
                self.n is None and \
                self.skew is None and \
                self.MurphyE == 0.:
            raise PolynomialParseException("No polynomials found")

        # Test that all required keys are there. Currently, only n is required
        if self.n is None:
            raise PolynomialParseException("Value of n missing")

        if len(tabpoly) > 0:
            polyg = tabpoly[0]
        if len(tabpoly) > 1:
            polyf = tabpoly[1]

        # Check that the polynomials were specified
        if not allow_only_one_poly and polyf.degree < 0:
            raise PolynomialParseException(
                    "No polynomial f specified (c: lines, or poly1)")
        if polyg.degree < 0:
            raise PolynomialParseException(
                    "No polynomial g specified (Y: lines, or poly0)")

        self.polyf = polyf
        self.polyg = polyg
        self.tabpoly = tabpoly
        return

    def __str__(self):
        arr = ["n: %d\n" % self.n]
        if self.skew is not None:
            arr += ["skew: %s\n" % self.skew]
        if self.type is not None:
            arr += ["type: %s\n" % self.type]
        if len(self.tabpoly) > 0:
            for i in range(len(self.tabpoly)):
                poltmp = self.tabpoly[i]
                arr += ["poly%d: %s" % (i, poltmp[0])]
                arr += [","+str(poltmp[j]) for j in range(1, len(poltmp))]
                arr += "\n"
        else:
            arr += self.polyf.as_lines("c")
            arr += self.polyg.as_lines("Y")
        if not self.MurphyE == 0.:
            if self.MurphyParams:
                arr.append("# MurphyE (%s) = %.3e\n" %
                           (self.MurphyParams, self.MurphyE))
            else:
                arr.append("# MurphyE = %.3e\n" % self.MurphyE)
        if self.revision is not None:
            arr.append("# found by revision %s\n" % self.revision)
        if not self.exp_E == 0.:
            arr.append("# exp_E %g\n" % self.exp_E)
        if len(self.tabpoly) > 0:
            for i in range(len(self.tabpoly)):
                arr.append("# poly%d = %s\n" % (i, str(self.tabpoly[i])))
        else:
            arr.append("# f(x) = %s\n" % str(self.polyf))
            arr.append("# g(x) = %s\n" % str(self.polyg))
        return "".join(arr)

    def __eq__(self, other):
        return self.polyf == other.polyf and \
                self.polyg == other.polyg and \
                self.n == other.n

    def __ne__(self, other):
        return not (self == other)

    def same_lc(self, other):
        """
        Returns true if both polynomial pairs have same degree and
        leading coefficient
        """
        return self.polyf.same_lc(other.polyf) \
            and self.polyg.same_lc(other.polyg)

    @property
    def nsides(self):
        if self.tabpoly:
            return len(self.tabpoly)
        else:
            return (self.polyf.degree >= 0) + (self.polyg.degree >= 0)

    def get_all_nonlinear_sides(self):
        return [s for s in range(self.nsides)
                if self.get_polynomial(s).degree > 1]

    def get_polynomial(self, side):
        """
        Returns one of the two polynomial as indexed by side
        """
        assert 0 <= side < self.nsides
        # Welp, f is side 1 and g is side 0 :(
        if side == 0:
            return self.polyg
        elif side == 1:
            return self.polyf
        else:
            return self.tabpoly[side]


class FilePath(object):
    """
    A class that represents a path to a file, where the path should be
    somewhat relocateable.

    In particular, we separate the path to the working directory, and the file
    path relative to the working directory. For persistent storage in the DB,
    the path relative to the workdir should be used, whereas for any file
    accesses, the full path needs to be used.
    It also piggy-backs a version information field.
    """

    def __init__(self, workdir, filepath, version=None):
        self.workdir = workdir.rstrip(os.sep)
        self.filepath = filepath
        self.version = version

    def __str__(self):
        return "%s%s%s" % (self.workdir, os.sep, self.filepath)

    def get_wdir_relative(self):
        return self.filepath

    def isfile(self):
        return os.path.isfile(str(self))

    def isdir(self):
        return os.path.isdir(str(self))

    def get_version(self):
        return self.version

    def mkdir(self, *, parent=False, mode=None):
        """
        Creates a directory.

        parent acts much like the Unix mkdir's '-p' parameter: required
        parent directories are created if they don't exist, and no error
        is raised if the directory to be created already exists.  If
        parent==True, a mode for the directory to be created can be
        specified as well.
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

    def rmtree(self, ignore_errors=False):
        rmtree(str(self), ignore_errors)


class WorkDir(object):
    """
    A class that allows generating file and directory names under a
    working directory.

    The directory layout is as follows:
    The current project (i.e., the factorization) has a jobname, e.g.,
    "RSA512". Each task may have a name, e.g., "sieving".
    A task can create various files under
    workdir/jobname.taskname.file
    or put them in a subdirectory
    workdir/jobname.taskname/file
    or, for multiple subdirectories,
    workdir/jobname.taskname/subdir/file

    It is also ok for tasks to have no particular name that is
    reflected in the filename hierarchy.

    >>> f = WorkDir("/foo/bar", "jobname", "taskname")
    >>> str(f.make_dirname("foo")).replace(os.sep,'/')
    '/foo/bar/jobname.foo/'
    >>> str(f.make_filename('file')).replace(os.sep,'/')
    '/foo/bar/jobname.file'
    >>> str(f.make_filename('file', subdir="foo")).replace(os.sep,'/')
    '/foo/bar/jobname.foo/file'
    >>> str(f.make_filename('file', prefix="bar", subdir='foo')).replace(
    ...     os.sep,'/')
    '/foo/bar/jobname.foo/jobname.bar.file'
    """

    def __init__(self, workdir, jobname=None, taskname=None):
        self.workdir = str(workdir).rstrip(os.sep)
        self.jobname = jobname
        self.taskname = taskname

    def path_in_workdir(self, filename, version=None):
        return FilePath(self.workdir, filename, version=version)

    def make_filename2(self, jobname=None, taskname=None, filename=None):
        if jobname is None:
            jobname = self.jobname
        if taskname is None:
            taskname = self.taskname
        filename_arr = [s for s in [jobname, taskname, filename] if s]
        return FilePath(self.workdir, ".".join(filename_arr))

    def make_dirname(self, subdir):
        """
        Make a directory name of the form workdir/jobname.prefix/
        """
        return self.path_in_workdir("".join([self.jobname,
                                             ".",
                                             subdir,
                                             os.sep]))

    def make_filename(self, name, prefix=None, subdir=None):
        """
        If subdir is None, make a filename of the form
        workdir/jobname.prefix.name or workdir/jobname.name depending on
        whether prefix is None or not.
        If subdir is not None, make a filename of the form
        workdir/jobname.subdir/jobname.prefix.name
        or workdir/jobname.subdir/name
        """
        components = [self.jobname]
        if subdir is not None:
            components += [".", subdir, os.sep]
            if prefix is not None:
                components += [self.jobname, ".", prefix, "."]
            components += [name]
        else:
            if prefix is not None:
                components += [".", prefix]
            components += [".", name]
        return self.path_in_workdir("".join(components))

    def get_workdir_jobname(self):
        return self.jobname

    def get_workdir_path(self):
        return self.workdir


class Statistics(object):
    """
    Class that holds statistics on program execution, and can merge two
    such statistics.
    """

    def __init__(self, conversions, formats):
        self.conversions = conversions
        self.stat_formats = formats
        self.stats = {}

    @staticmethod
    def typecast(values, types):
        """
        Cast the values in values to the types specified in types
        """
        if type(types) is type:
            return [types(v) for v in values]
        else:
            return [t(v) for (v, t) in zip(values, types)]

    @staticmethod
    def _to_str(stat):
        """
        Convert one statistic to a string
        """
        return " ".join(map(str, stat))

    @staticmethod
    def _from_str(string, types):
        """
        Convert a string (probably from a state dict) to a statistic
        """
        return Statistics.typecast(string.split(), types)

    def from_dict(self, stats):
        """
        Initialise values in self from the strings in the "stats"
        dictionary
        """
        for (key, types, defaults,
             combine, regex, allow_several) in self.conversions:
            if key in stats:
                if key in self.stats:
                    print("duplicate %s\n" % key)
                assert key not in self.stats
                self.stats[key] = self._from_str(stats.get(key, defaults),
                                                 types)
                assert self.stats[key] is not None

    def parse_line(self, line):
        """
        Parse one line of program output and look for statistics.

        If they are found, they are added to self.stats.
        """
        for (key, types, defaults,
             combine, regex, allow_several) in self.conversions:
            match = regex.match(line)
            if match:
                # print (pattern.pattern, match.groups())
                # Optional groups that did not match are returned as None.
                # Skip over those so typecast doesn't raise TypeError
                groups = [group
                          for group in match.groups()
                          if group is not None]
                new_val = self.typecast(groups, types)
                if not allow_several:
                    assert key not in self.stats
                    self.stats[key] = new_val
                else:
                    # Some output files inherently have several values.
                    # This is the case of bwc output files if we use
                    # multiple sequences.
                    if key in self.stats:
                        self.stats[key] = combine(self.stats[key], new_val)
                    else:
                        self.stats[key] = new_val
                assert self.stats[key] is not None

    def merge_one_stat(self, key, new_val, combine):
        if key in self.stats:
            self.stats[key] = combine(self.stats[key], new_val)
        else:
            self.stats[key] = new_val
        assert self.stats[key] is not None
        # print(self.stats)

    def merge_stats(self, new_stats):
        """
        Merge the stats currently in self with the Statistics in
        "new_stats"
        """

        assert self.conversions == new_stats.conversions
        for (key, types, defaults,
             combine, regex, allow_several) in self.conversions:
            if key in new_stats.stats:
                self.merge_one_stat(key, new_stats.stats[key], combine)

    def as_dict(self):
        return {key: self._to_str(self.stats[key]) for key in self.stats}

    def as_strings(self):
        """
        Convert statistics to lines of output

        The self.stat_formats is an array, with each entry corresponding to
        a line that should be output.
        Each such entry is again an array, containing the format strings that
        should be used for the conversion of statistics. If a conversion
        fails with a KeyError or an IndexError, it is silently skipped over.
        This is to allow producing lines on which some statistics are not
        printed if the value is not known.
        """
        result = []
        errors = []
        for format_arr in self.stat_formats:
            line = []
            for format_str in format_arr:
                try:
                    line.append(format_str.format(**self.stats))
                except KeyError:
                    errors.append("KeyError with \"%s\"" % format_str)
                except IndexError:
                    errors.append("IndexError with \"%s\"" % format_str)
            if line:
                result.append("".join(line))
        if len(errors) > 0:
            errors.append("(registered stats: %s)" % self.stats)
            return result, errors
        else:
            return result, None

    # Helper functions for processing statistics.
    # We can't make them @staticmethod or references are not callable
    def add_list(*lists):
        """
        Add zero or more lists elementwise.

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
        """
        From two lists, one containing values and the other containing
        the respective sample sizes (i.e., weights of the values), compute
        the combined mean (i.e. the weighted mean of the values).
        The two lists must have equal length.
        """
        assert len(means) == len(samples)
        total_samples = sum(samples)
        weighted_sum = sum(Statistics.weigh(means, samples))
        return [weighted_sum / total_samples, total_samples]

    def zip_combine_mean(*lists):
        """
        From a list of 2-tuples, each tuple containing a value and a
        weight, compute the weighted average of the values.
        """
        for ell in lists:
            assert len(ell) == 2
        (means, samples) = zip(*lists)
        return Statistics.combine_mean(means, samples)

    def combine_stats(*stats):
        """
        Computes the combined mean and std.dev. for the stats

        stats is a list of 3-tuples, each containing number of sample points,
        mean, and standard deviations.
        Returns a 3-tuple with the combined number of sample points, mean,
        and standard deviations.
        """

        # Samples is a list containing the first item (number of samples)
        # of each item of stats, means is list of means, stddevs is list
        # of standard deviations.
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
        """
        Test function for combine_stats()

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
        lists = [[randrange(100) for i in range(ell)] for ell in lengths]
        stats = [(length, mean(ell), stddev(ell))
                 for (length, ell) in zip(lengths, lists)]

        combined = sum(lists, [])

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
        concat = sum(lists, [])
        concat.sort()
        return concat[0:n]

    def parse_stats(self, filename):
        """
        Parse statistics from the file with name "filename" and merge them
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
        # The name of the task in a simple form that can be used as a
        # Python dictionary key, a directory name, part of a file name,
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
    def database_state_table_name(self):
        """
        Prefix string for table names

        By default, the table name prefix is the name attribute, but this
        can be overridden
        """
        return self.name

    def make_tablename(self, extra=None):
        """
        Return a name for a DB table
        """
        # Maybe replace SQL-disallowed characters here, like digits and '.' ?
        # Could be tricky to avoid collisions
        name = self.database_state_table_name
        if extra:
            name = name + '_' + extra
        wudb.check_tablename(name)
        return name


class HasState(MakesTablenames, wudb.HasDbConnection):
    """
    Declares that the class has a DB-backed dictionary in which the class
    can store state information.

    The dictionary is available as an instance attribute "state".
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        name = self.make_tablename()
        self.state = self.make_db_dict(name, connection=self.db_connection)


class FilesCreator(MakesTablenames,
                   wudb.HasDbConnection,
                   metaclass=abc.ABCMeta):
    """
    A base class for classes that produce a list of output files, with
    some auxiliary information stored with each file (e.g., nr. of
    relations).  This info is stored in the form of a DB-backed
    dictionary, with the file name as the key and the auxiliary data as
    the value.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        tablename = self.make_tablename("outputfiles")
        self.output_files = self.make_db_dict(tablename,
                                              connection=self.db_connection)

    def add_output_files(self, filenames, *, commit):
        """
        Adds a dict of files to the list of existing output files
        """
        final_files = {}
        for filename in filenames:
            if filename in self.output_files:
                self.logger.warning(f"{filename} already in"
                                    " output files table")
            else:
                final_files[filename] = filenames[filename]
        self.output_files.update(final_files, commit=commit)

    def get_output_filenames(self, condition=None):
        """
        Return output file names, optionally those that match a condition

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
    """
    Base class for HasStatistics and SimpleStatistics that terminates the
    print_stats() call chain.
    """
    def print_stats(self):
        pass


class HasStatistics(BaseStatistics,
                    HasState,
                    DoesLogging,
                    metaclass=abc.ABCMeta):
    @property
    def stat_conversions(self):
        """
        Sub-classes should override
        """
        return []

    @property
    def stat_formats(self):
        """
        Sub-classes should override
        """
        return []

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.statistics = Statistics(self.stat_conversions, self.stat_formats)
        self.statistics.from_dict(self.state)

    def get_statistics_as_strings(self):
        """
        Return the statistics collected so far as a List of strings.

        Sub-classes can override to add/remove/change strings.

        Typically, classes that subclass HasStatistics *AND* DoesImport
        may wish to report stats differently if import has happened.
        """
        result, errors = self.statistics.as_strings()
        return result, errors

    def print_stats(self):
        stat_msgs, errors = self.get_statistics_as_strings()
        if stat_msgs:
            self.logger.info("Aggregate statistics:")
            for msg in stat_msgs:
                self.logger.info(msg)
        if errors is not None:
            self.logger.warning("some stats could not be displayed"
                                f" for {self.name}"
                                " (see log file for debug info)")
            for e in errors:
                self.logger.debug(e)
            if "STATS_PARSING_ERRORS_ARE_FATAL" in os.environ:
                raise RuntimeError("Aborting now, since"
                                   " STATS_PARSING_ERRORS_ARE_FATAL is set")
        super().print_stats()

    def parse_stats(self, filename, *, commit):
        # self.logger.info("Parsing filename %s\n", filename)
        new_stats = self.statistics.parse_stats(filename)
        self.logger.debug("Newly arrived stats: %s", new_stats)
        update = self.statistics.as_dict()
        self.logger.debug("Combined stats: %s", update)
        self.state.update(update, commit=commit)


class SimpleStatistics(BaseStatistics,
                       HasState,
                       DoesLogging,
                       metaclass=abc.ABCMeta):

    @abc.abstractproperty
    def programs(self):
        # A list of classes of Programs which this tasks uses
        pass

    def print_cpu_real_time(self, cputotal, realtotal, program):
        """
        Print cpu and/or real time to logger
        """
        # Uses self only for access to the logger
        pairs = zip((cputotal, realtotal), ("cpu", "real"))
        usepairs = [pair for pair in pairs if pair[0]]
        if usepairs:
            printformat = "/".join(["%g"] * len(usepairs))
            usepairs = tuple(zip(*usepairs))
            timestr = '/'.join(usepairs[1])
            self.logger.info("Total %s time for %s: " + printformat,
                             timestr, program, *usepairs[0])

    @staticmethod
    def keyname(is_cpu, programname):
        if is_cpu:
            return "cputime_%s" % programname
        else:
            return "realtime_%s" % programname

    def update_cpu_real_time(self, programname,
                             cpu=None, real=None, commit=True):
        """
        Add seconds to the statistics of cpu time spent by program,
        and return the new total.
        """
        assert isinstance(programname, str)
        update = {}
        for (is_cpu, t) in ((True, cpu), (False, real)):
            if t is not None:
                key = self.keyname(is_cpu, programname)
                update[key] = self.state.get(key, 0.) + t
        if update:
            self.state.update(update, commit=commit)

    def get_cpu_real_time(self, program):
        """
        Return list of cpu and real time spent by program
        """
        return [self.state.get(self.keyname(is_cpu, program.name), 0.)
                for is_cpu in (True, False)]

    def get_total_cpu_or_real_time(self, is_cpu):
        """
        Return tuple with number of seconds of cpu and real time spent
        by all programs of this Task
        """
        if not self.programs:
            return 0
        times = [self.get_cpu_real_time(p) for p, o, i in self.programs]
        times = tuple(map(sum, zip(*times)))
        return times[0 if is_cpu else 1]

    def print_stats(self):
        for program, o, i in self.programs:
            cputotal, realtotal = self.get_cpu_real_time(program)
            self.print_cpu_real_time(cputotal, realtotal, program.name)
        super().print_stats()


class Runnable(object):
    @abc.abstractmethod
    def run(self):
        pass


class DoesImport(DoesLogging, cadoparams.UseParameters, Runnable,
                 metaclass=abc.ABCMeta):
    @abc.abstractproperty
    def paramnames(self):
        return self.join_params(super().paramnames, {"import": None})

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._did_import = False

    def run(self):
        super().run()

    def do_import(self):
        if "import" in self.params and not self._did_import:
            self.import_files(self.params["import"])
            self._did_import = True

    def import_files(self, input_filename):
        if input_filename.startswith('@'):
            self.logger.info("Importing files listed in %s",
                             input_filename[1:])
            with open(input_filename[1:], "r") as f:
                filenames = f.read().splitlines()
        else:
            self.logger.info("Importing file %s", input_filename)
            filenames = [input_filename]
        for filename in filenames:
            self.import_one_file(filename)

    def did_import(self):
        return self._did_import

    @abc.abstractmethod
    def import_one_file(self, filename):
        pass


def chain_dict(d1, d2):
    """
    Chain two mappings.

    If d[x] == y and e[y] == z, then chain_dict(d, e)[x] == z.
    >>> chain_dict({1: 17}, {17: 42})
    {1: 42}
    """
    return {key: d2[value] for key, value in d1.items()}


class RealTimeOutputFilter:

    def __init__(self, logger, filename):
        # the buffering=1 option is needed so that the stdout files are
        # not buffered, cf
        # docs.python.org/3.7/library/functions.html?highlight=open#open
        self.stdout = open(filename, mode="w", buffering=1)
        self.logger = logger

    def filter(self, data):
        self.stdout.write(data)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.stdout.close()


class Task(patterns.Colleague, SimpleStatistics, HasState, DoesLogging,
           cadoparams.UseParameters, Runnable, metaclass=abc.ABCMeta):
    """
    A base class that represents one task that needs to be processed.

    Sub-classes must define class variables:
    """

    # Properties that subclasses need to define

    @abc.abstractproperty
    def programs(self):
        # A tuple of 3-tuples, with each 3-tuple containing
        # 1. the class of Program which this tasks uses
        # 2. a tuple of parameters to the Program which the Task computes and
        #    which therefore should not be filled in from the Parameters file
        # 3. a dict of parameters which are file names and which should be
        #    filled in by sending Requests to other Tasks. This also enables
        #    testing whether input files have been changed by the other Task.
        pass

    @abc.abstractproperty
    def paramnames(self):
        # Parameters that all tasks use
        return self.join_params(super().paramnames,
                                {"name": str, "workdir": str, "run": True})

    @property
    def param_nodename(self):
        # avoid segregating our parameters, which are user-visible
        # things, underneath tree nodes whose name depends on some
        # implementation detail which is the task name. Except in
        # specific cases, a "task" does not (no longer) define a nesting
        # level in the parameter hierarchy.
        #
        # return self.name
        return None

    def __init__(self, *, mediator, db, parameters, path_prefix):
        """
        Sets up a database connection and a DB-backed dictionary for
        parameters. Reads parameters from DB, and merges with hierarchical

        parameters in the parameters argument. Parameters passed in by
        parameters argument do not override values in the DB-backed
        parameter dictionary.
        """

        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.logger.debug("Enter Task.__init__(%s)",
                          self.name)
        self.logger.debug("state = %s", self.state)
        # Set default parameters for this task, if any are given
        # The parameters that are relative to the _task_ itself ("run",
        # "wdir"), and that affect all programs in a given task (e.g.,
        # below "mergetask": the "merge" and "replay" programs) are
        # prefixed with the task name, simply because there is no program
        # to play with.
        self.params = self.parameters.myparams(self.paramnames, self.name)
        self.logger.debug("self.parameters = %s", self.parameters)
        self.logger.debug("params = %s", self.params)
        # Set default parameters for our programs
        # The progparams entries should not be modified after a class'
        # constuctor (within __init__() is fine tho)
        self.progparams = []
        maindict = self.parameters.parameters
        for prog, override, needed_input in self.programs:
            # Parameters listed in needed_input are assumed to be overridden
            for key in (set(override) & set(needed_input)):
                self.logger.warning("Parameter %s listed in both overridden "
                                    "parameters and in input files for %s, "
                                    "only one is needed", key, prog.name)
            prog_param_path = self.parameters.get_param_path() + [prog.name]
            progparams = self.parameters.myparams(prog.get_accepted_keys(),
                                                  prog.name)
            for c in progparams:
                finergrain = '.'.join(prog_param_path+[c])
                coarsegrain = maindict.locate(finergrain)
                self.logger.debug("%s found from %s", finergrain, coarsegrain)

            for param in (set(needed_input) | set(override)) & set(progparams):
                finergrain = '.'.join(prog_param_path+[param])
                coarsegrain = maindict.locate(finergrain)
                # Whenever we see a parameter that is marked as override,
                # we will discard it and let the task level fill data for
                # this parameter. There are cases where this really is a
                # user error and we want to complain:
                #  - when the parameter file *explicitly* sets this
                #    parameter at this level. This does not make sense
                #    and is a troubling no-op. Typical example is
                #    specifying tasks.linalg.bwc.m in dlp mode.
                #  - when the parameter file sets it at a level above,
                #    but the task level does *not* know about this
                #    parameter anyway. This is ignored as well, and the
                #    task level will fill that parameter based on data it
                #    knows. But leaving the user with the feeling that he
                #    might be able to control that parameter is
                #    inelegant. A typical example is
                #    tasks.sieve.makefb.lim (which the tasks.sieve level
                #    sets based on the lim0 and lim1 parameters it knows
                #    about). Likewise, many "out" parameters behave
                #    similarly.
                if finergrain == coarsegrain or \
                        param not in set(self.paramnames):
                    self.logger.error('Parameter "%s" for program "%s" is '
                                      'generated at run time and cannot be '
                                      'supplied through the parameter file',
                                      param, prog.name)
                    self.logger.error('Ignoring %s, we rely on '
                                      '%s to compute it '
                                      'based on parameters at level %s only',
                                      '.'.join(path_prefix+[prog.name, param]),
                                      self.__class__,
                                      '.'.join(path_prefix))
                # We'll anyway discard it, but it's normal if we
                # inherited the parameter from a level above.
                del progparams[param]

            self.progparams.append(progparams)

        # FIXME: whether to init workdir or not should not be controlled via
        # presence of a "workdir" parameter, but by class definition
        if "workdir" in self.params:
            self.workdir = WorkDir(self.params["workdir"],
                                   self.params["name"],
                                   self.name)
        # Request mediator to run this task. It the "run" parameter is set
        # to false, then run() below will abort.
        self.send_notification(Notification.WANT_TO_RUN, None)
        self.logger.debug("Exit Task.__init__(%s)", self.name)
        return

    def run(self):
        if not self.params["run"]:
            self.logger.info("Stopping at %s", self.name)
            raise EarlyStopException("Job stops here"
                                     " because of a forcibly disabled task"
                                     " -- stopped at " + self.name)
        self.logger.info("Starting")
        self.logger.debug("%s.run(): Task state: %s", self.name, self.state)
        super().run()
        # Make set of requests so multiply listed requests are sent only once
        # The input_file dict maps key -> Request. Make set union of requests
        if self.programs:
            requests = set.union(*[set(i.values())
                                   for p, o, i in self.programs])
        else:
            requests = set()
        # Make dict mapping Request -> answer (i.e., FileName object)
        answers = self.batch_request(dict(zip(requests, requests)))
        # Make list of dicts mapping key -> answer
        self._input_files = [chain_dict(i, answers)
                             for p, o, i in self.programs]
        # Since merged_args is re-generated in each run(), the subclass can
        # modify it as it pleases (unlike progparams)
        # For each program, merge the progparams and the input_files dicts
        self.merged_args = [dict(p.items() | i.items())
                            for p, i in zip(self.progparams,
                                            self._input_files)]

    def translate_input_filename(self, filename):
        return filename

    def test_outputfile_exists(self, filename):
        return filename.isfile()

    def cmp_input_version(self, state_key, version):
        """
        Compare the version of the input file with the version we have
        processed before

        Returns None if filename does not include file version information,
        returns -2 if we have never processed the file before,
        returns -1 if the previously processed file is older,
        returns 0 if they are the same version,
        throws an exception if the processed version is newer than the
        current one.
        """
        if version is None:
            return None
        if state_key not in self.state:
            return -2
        if self.state[state_key] < version:
            return -1
        if self.state[state_key] == version:
            return 0
        raise ValueError("Previously processed version is newer than current")

    @staticmethod
    def _input_file_to_state_dict(key, index, filename):
        return ("processed_version_%d_%s" % (index, key),
                filename.get_version())

    def have_new_input_files(self):
        # Change this to "self.logger.info" for showing each check on screen
        log = self.logger.debug
        result = False
        for index, input_files in enumerate(self._input_files):
            for key, filename in input_files.items():
                (state_key, version) = \
                    self._input_file_to_state_dict(key, index, filename)
                c = self.cmp_input_version(state_key, version)
                if c == -2:
                    log("File %s was not processed before", filename)
                    result = True
                if c == -1:
                    log("File %s is newer than last time", filename)
                    result = True
        # Collapse programs from all dict into one set
        all = set.union(*[set(i.values()) for i in self._input_files])
        if result is False and all:
            (n, v) = (("", "is"), ("s", "are"))[len(all) > 1]
            log("Input file%s %s %s unchanged since last time",
                n, ", ".join(map(str, all)), v)
        return result

    def remember_input_versions(self, commit=True):
        update = {}
        for index, input_files in enumerate(self._input_files):
            for (key, filename) in input_files.items():
                (state_key, version) = \
                    self._input_file_to_state_dict(key, index, filename)
                if version is not None:
                    update[state_key] = version
        self.state.update(update, commit)

    @staticmethod
    def check_files_exist(filenames, filedesc, shouldexist):
        """
        Check that the output files in "filenames" exist or don't exist,
        according to shouldexist.

        Raise IOError if any check fails, return None
        """
        for filename in filenames:
            if isinstance(filename, FilePath):
                exists = filename.isfile()
            else:
                exists = os.path.isfile(filename)
            if shouldexist and not exists:
                raise IOError(f"{filedesc} file {filename} does not exist")
            elif not shouldexist and exists:
                raise IOError(f"{filedesc} file {filename} already exists")
        return

    # The next two function (make_wuname and split_wuname) go together,
    # one produces a workunit name from the name of the factorization,
    # the task name, and a task-provided identifier, and the other
    # function splits them again

    wu_paste_char = '_'
    wu_attempt_char = '__R'

    def make_wuname(self, identifier, attempt=None):
        """
        Generates a wuname from project name, task name, identifier, and
        attempt number.
        """
        arr = [self.params["name"], self.name]
        assert self.wu_attempt_char not in self.name
        if identifier:
            arr.append(identifier)
        wuname = self.wu_paste_char.join(arr)
        if attempt is not None:
            wuname += "%s%d" % (self.wu_attempt_char, attempt)
        return wuname

    def split_wuname(self, wuname):
        """
        Splits a wuname into project name, task name, identifier, and
        attempt number.

        Always returns a list of length 4; if there is not attempt given in
        the wuname, then the last array entry is None

        >>> # Test how various names are split. Note that if
        >>> # wu_paste_char is a substring of wu_attempt_char, we cannot
        >>> # have wu_attempt_char in the identifier name since that
        >>> # would lead to ambiguous parsing.
        >>> class Klass():
        ...     params = {"name": None}
        ...     wu_paste_char = '_'
        ...     wu_attempt_char = '__R'
        >>> inst = Klass()
        >>> from itertools import product
        >>> ranges = [["", Klass.wu_paste_char]] * 4
        >>> ranges += [[""]]*2
        >>> prod = product(*ranges)
        >>> for sep in prod:
        ...     inst.params["name"] = "%s%sprojectname%s%s" % sep[0:4]
        ...     inst.name = "%staskname%s" % sep[4:6]
        ...     for attempt in [None, 2, 3]:
        ...         identifier = "identifier"
        ...         wuname = Task.make_wuname(inst, "identifier",
        ...                                   attempt=attempt)
        ...         wu_split = Task.split_wuname(inst, wuname)
        ...         assert wu_split == [inst.params["name"],
        ...                             inst.name,
        ...                             identifier,
        ...                             attempt]
        """
        attempt = None
        # Split off attempt number, if available. We do it first so that
        # things work correctly even if wu_paste_char is a substring of
        # wu_attempt_char
        m = re.match(r'(.*)' + self.wu_attempt_char + r'(\d+)', wuname)
        if m:
            wuname = m.group(1)
            attempt = int(m.group(2))
        arr = wuname.rsplit(self.wu_paste_char, 2)
        assert len(arr) == 3
        if arr[-1] == "":
            arr[-1] = None
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
            self.logger.error("Stderr output (last 10 lines only) follows%s:",
                              stderrmsg)
            for ell in stderr.decode().split('\n')[-10:]:
                self.logger.error("\t" + ell)

    def submit_command(self, command, identifier,
                       commit=True, log_errors=False):
        """
        Run a command.
        Return the result tuple. If the caller is an Observer, also send
        result to updateObserver().
        """

        # Task objects may submit commands with an identifier, but in
        # most cases it's unnecessary, and anyway the identifier isn't
        # meaningful to any function here (in contrast with the
        # ClientServerTask situation)

        wuname = self.make_wuname(identifier)
        process = cadocommand.Command(command)
        cputime_used = os.times()[2]  # CPU time of child processes
        realtime_used = time.time()
        (rc, stdout, stderr) = process.wait()
        cputime_used = os.times()[2] - cputime_used
        realtime_used = time.time() - realtime_used
        self.update_cpu_real_time(command.name,
                                  cputime_used, realtime_used, commit)
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
            self.logger.message("Output files are: %s",
                                ", ".join(output_files))
        return identifier

    def send_notification(self, key, value):
        """
        Wrapper around Colleague.send_notification() that instantiates a
        Notification with self as the sender
        """
        notification = Notification(self, key, value)
        super().send_notification(notification)

    def send_request(self, key, *args):
        """
        Wrapper around Colleague.send_request() that instantiates a
        Request with self as the sender
        """
        request = Request(self, key, *args)
        return super().send_request(request)

    def batch_request(self, requests):
        """
        Given a dict from keys to Request objects, return a dict with the
        same keys to the results of the requests.
        """
        return {key: self.send_request(request)
                for key, request in requests.items()}

    def get_number_outstanding_wus(self):
        return 0

    def verification(self, wuid, ok, *, commit):
        pass

    def get_state_filename(self, key, version=None):
        """
        Return a file name stored in self.state as a FilePath object

        If a version parameter is passed, then this version is set as the
        version field of the FilePath object. If no parameter is passed,
        but our state includes an "output_version" key, then that is
        used.
        """
        if key not in self.state:
            return None
        if version is None:
            version = self.state.get("output_version", None)
        return self.workdir.path_in_workdir(self.state[key], version)

    def make_std_paths(self, progname, do_increment=True, prefix=None):
        count = self.state.get("stdiocount", 0)
        if do_increment:
            count += 1
        did_increment = do_increment
        while True:
            try:
                stdoutname = "%s.stdout.%d" % (progname, count)
                stderrname = "%s.stderr.%d" % (progname, count)
                self.check_files_exist((stdoutname, stderrname),
                                       "stdio",
                                       shouldexist=False)
            except IOError:
                count += 1
                did_increment = True
                self.logger.warning("Stdout or stderr files with index %d "
                                    "already exist", count)
            else:
                break
        stdoutpath = self.workdir.make_filename(stdoutname, prefix=prefix)
        stderrpath = self.workdir.make_filename(stderrname, prefix=prefix)
        if did_increment:
            self.state["stdiocount"] = count
        return (stdoutpath, stderrpath)

    def make_filelist(self, files, prefix=None):
        """
        Create file file containing a list of files, one per line
        """
        filelist_idx = self.state.get("filelist_idx", 0) + 1
        self.state["filelist_idx"] = filelist_idx
        filelistname = self.workdir.make_filename(f"filelist.{filelist_idx}",
                                                  prefix=prefix)
        with filelistname.open("w") as filelistfile:
            filelistfile.write("\n".join(files) + "\n")
        return filelistname

    def collect_usable_parameters(self, rl):
        message = []
        message.append("Parameters used by Task %s" % self.name)
        prefix = '.'.join(self.parameters.get_param_path())
        for p in self.paramnames:
            message.append("  %s.%s.%s" % (prefix, self.name, p))
            rl[p].append(prefix)
        for prog, override, needed_input in self.programs:
            message.append("  Parameters for program %s"
                           " (general form %s.%s.*)" % (
                               prog.name, prefix, prog.name))
            for p in sorted(prog.get_accepted_keys()):
                t = "%s.%s.%s" % (prefix, prog.name, p)
                rl[p].append("%s.%s" % (prefix, prog.name))
                if p in set(override):
                    message.append("    [excluding internal parameter %s]" % t)
                elif p in set(needed_input):
                    message.append("    [excluding internal file name %s]" % t)
                else:
                    message.append("    %s" % t)
        message.append("")
        return "\n".join(message)


class ClientServerTask(Task, wudb.UsesWorkunitDb, patterns.Observer):
    # Note that we get self.wuar = self.make_wu_access(db.connect()) via
    # inheritance of wudb.UsesWorkunitDb

    @abc.abstractproperty
    def paramnames(self):
        # maxwuerror can be increased if jobs are often killed badly.
        return self.join_params(super().paramnames,
                                {"maxwu": 10,
                                 "wutimeout": 10800,  # Default: 3h
                                 "wutimeoutcheck": 60,  # Default: every minute
                                 "maxresubmit": 5,
                                 "maxwuerror": 2,
                                 "maxtimedout": 100,
                                 "maxfailed": 100})

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.state.setdefault("wu_submitted", 0)
        # wu_received only counts the WUs as individual units. It is not
        # robust to range parameter changing in the course of the
        # computation.
        # -> wu_range_received is a better measure of what we've done thus far.
        self.state.setdefault("wu_received", 0)
        self.state.setdefault("wu_range_received", 0)
        self.state.setdefault("wu_timedout", 0)
        self.state.setdefault("wu_failed", 0)
        assert self.get_number_outstanding_wus() >= 0
        # start_real_time will be a float giving the number of seconds since
        # Jan 1 1900 at the beginning of the task
        self.state.update({"start_real_time": 0})
        # start_achievement is a variable that tells us how far we were at
        # the beginning of this run (for example if a factorization is
        # restarted in the middle of a polyselect or sieve task.)
        # It should be in [0,1], and if not initialized yet it is -1.
        self.state.update({"start_achievement": -1})
        self.send_notification(Notification.SUBSCRIBE_WU_NOTIFICATIONS, None)
        self.clients = None

    def submit_wu(self, wu, commit=True):
        """
        Submit a WU and update wu_submitted counter
        """
        # at beginning of the task, set "start_real_time" to the number of
        # seconds since Jan 1 1900
        if self.state["start_real_time"] == 0:
            delta = datetime.datetime.now() - datetime.datetime(1900, 1, 1)
            self.state.update({"start_real_time": delta.total_seconds()})
        key = "wu_submitted"
        self.state.update({key: self.state[key] + 1}, commit=False)
        self.wuar.create(wu, commit=commit)

    def cancel_wu(self, wuid, commit=True):
        """
        Cancel a WU and update wu_timedout counter
        """
        self.logger.debug("Cancelling: %s", wuid)
        key = "wu_timedout"
        maxtimedout = self.params["maxtimedout"]
        if not self.state[key] < maxtimedout:
            self.logger.error("Exceeded maximum number of timed out "
                              "workunits, maxtimedout=%d ", maxtimedout)
            raise Exception("Too many timed out work units."
                            " Please increase tasks.maxtimedout"
                            " (current value is %d)" % maxtimedout)
        self.state.update({key: self.state[key] + 1}, commit=False)
        self.wuar.cancel(wuid, commit=commit)

    def enough_work_received(self):
        """
        This is independent of the loop in the run() function. Here,
        we're only telling whether it makes sense to drop outstanding WUs
        because the amount of work that was collected by the server is
        already enough to proceed with the next steps of the computation.
        To disable this early abort check, the default behaviour is
        simply to return False, in which case the client-server task will
        finish once all the WUs that were scheduled from the run()
        command are received.
        """
        return False

    def submit_command(self, command, identifier,
                       commit=True, log_errors=False):
        """
        Submit a workunit to the database.

        This is the main entry point of a ClientServerTask, which gets
        called by subclasses.
        """

        # client-server tasks *must* have identifiers...
        assert identifier is not None

        if re.match(r'\d+-\d+', identifier):
            pass
        elif re.match(r'\d+', identifier):
            identifier = "%d-%d" % (int(identifier), int(identifier)+1)
        else:
            raise ValueError("Bad WU identifer %s in %s" % (identifier,
                                                            self.name))

        # XXX design trap here. If clients are pulling hard on the server,
        # it may never be able to keep the number of available WUs
        # afloat. It's pulled low. And then, we never enter self.wait(),
        # which also mean that we never request wu results!
        #
        # XXX There's also a second flaw, perhaps even more annoying: the
        # server will only terminate if it's done with _all_ the WUs it
        # has scheduled, which might not be desirable in all cases. If
        # enough WUs have been received to obtain enough relations (in
        # the SievingTask case), then we should probably proceed and not
        # wait for the WUs that are still lagging.
        while self.get_number_available_wus() >= self.params["maxwu"]:
            if self.enough_work_received():
                # This fixes the second issue above. But really, it's a
                # telling sign that the control flow should be improved.
                o = self.get_number_outstanding_wus()
                a = self.get_number_available_wus()
                self.logger.info(f"{self.title} finishes early, since"
                                 " enough work has been collected already."
                                 " Dropping %d outstanding WUs, %d or which"
                                 " were not even assigned yet.", o, a)
                return
            self.wait()

        # To fix the first issue above, let's ensure that we call
        # GET_WU_RESULT at least once. This should ensure some minimal
        # fairness overall, but it's still true that backlog may
        # accumulate in cases where we fetch only one wu result for each
        # wu we produce.

        # It's also quite unsatisfactory that we do time.sleep() only
        # based on a condition that is the absence of result (from within
        # wait), and lose the occasion for a fast-path recreation of WUs.

        self.send_request(Request.GET_WU_RESULT)

        wuid = self.make_wuname(identifier)

        # ...and we want to be sure that the range size can be extracted
        # from the identifier.
        assert self.get_wusize(wuid) > 0

        wutext = command.make_wu(wuid)

        for filename in command.get_exec_files() + command.get_input_files():
            basename = os.path.basename(filename)
            self.send_notification(Notification.REGISTER_FILENAME,
                                   {basename: filename})

        self.logger.info("Adding workunit %s to database", wuid)
        # Note that submit_wu parses the workunit to deduce the wuid
        # (which was created by make_wu in the first place)
        self.submit_wu(wutext, commit=commit)
        # Write command line to a file
        cmdline = command.make_command_line()
        client_cmd_filename = self.workdir.make_filename2(taskname="",
                                                          filename="wucmd")
        with client_cmd_filename.open("a") as client_cmd_file:
            client_cmd_file.write("# Command for work unit: %s\n%s\n" %
                                  (wuid, cmdline))

    def get_eta(self):
        delta = datetime.datetime.now() - datetime.datetime(1900, 1, 1)
        seconds = delta.total_seconds() - self.state["start_real_time"]
        a = self.get_achievement()
        a0 = self.state["start_achievement"]
        if a0 == -1:
            self.state["start_achievement"] = a
            a0 = a
        elif a0 > a:
            # if a0 > a, it means we had a failing filtering try, which means
            # a had attained 100%, and then decreased to say 95% for example,
            # thus we need to update a0 by multiplying it by a
            a0 = a0 * a
        try:
            remaining_time = seconds / (a - a0) * (1.0 - a)
            now = datetime.datetime.now()
            arrival = now + datetime.timedelta(seconds=remaining_time)
            return arrival.ctime()
        except (OverflowError, ZeroDivisionError):
            return "Unknown"

    def get_wusize(self, wuid):
        """
        parses a wuid that is relevant for the current task, and
        return the size of the attached workunit
        """
        (name, task, identifier, attempt) = self.split_wuname(wuid)
        m = re.match(r'(\d+)-(\d+)', identifier)
        if not m:
            raise ValueError(wuid)
        return int(m.group(2))-int(m.group(1))

    def verification(self, wuid, ok, *, commit):
        """
        Mark a workunit as verified ok or verified with error and update
        wu_received counter
        """
        ok_str = "ok" if ok else "not ok"
        assert self.get_number_outstanding_wus() >= 1
        key = "wu_received"
        self.state.update({key: self.state[key] + 1}, commit=False)
        # Like wu_received, we count here the failed workunits as well as
        # the good ones. So the range in wu_range_received might be
        # different from what has been "done".
        key = "wu_range_received"
        z = self.get_wusize(wuid)
        self.state.update({key: self.state[key] + z}, commit=False)
        # only print ETA when achievement > 0 to avoid division by zero
        a = self.get_achievement()
        if a > 0:
            self.logger.info("Marking workunit %s as %s (%.1f%% => ETA %s)",
                             wuid, ok_str, 100.0 * a, self.get_eta())
        self.wuar.verification(wuid, ok, commit=commit)

    def cancel_available_wus(self):
        self.logger.info("Cancelling remaining workunits")
        self.wuar.cancel_all_available()

    def get_number_outstanding_wus(self):
        return self.state["wu_submitted"] \
                - self.state["wu_received"] \
                - self.state["wu_timedout"]

    def get_number_available_wus(self):
        return self.wuar.count_available()

    def test_outputfile_exists(self, filename):
        # Can't test
        return False

    def wait(self):
        # Ask the mediator to check for workunits of status Received, and
        # if there are any, to send WU result notifications to the
        # subscribed listeners.
        # If we get notification on new results reliably from the HTTP
        # server, we might not need this poll. But they probably won't be
        # totally reliable
        if self.clients is None:
            self.clients = self.send_request(Request.GET_CLIENTS)
        if not self.send_request(Request.GET_WU_RESULT):
            self.resubmit_timed_out_wus()
            # Watch at least any local clients we may have.
            for c in self.clients:
                c.check_health()
            time.sleep(1)

    def resubmit_one_wu(self, wu, commit=True, maxresubmit=None):
        """
        Takes a Workunit instance and adds it to workunits table under
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
            return
        new_wuid = self.make_wuname(identifier, attempt)
        wu.set_id(new_wuid)
        self.logger.info("Resubmitting workunit %s as %s", wuid, new_wuid)
        self.submit_wu(wu, commit=commit)

    def resubmit_timed_out_wus(self):
        """
        Check for any timed out workunits and resubmit them
        """
        # We don't store the lastcheck in state as we do *not* want to check
        # instantly when we start up - clients should get a chance to upload
        # results first
        now = time.time()
        if not hasattr(self, "last_timeout_check"):
            self.logger.debug("Setting last timeout check to %f", now)
            self.last_timeout_check = now
            return

        check_every = self.params["wutimeoutcheck"]  # Check every xx seconds
        if self.last_timeout_check + check_every >= now:
            # self.logger.info("It's not time to check yet, now = %f", now)
            return
        self.last_timeout_check = now

        timeout = self.params["wutimeout"]
        delta = datetime.timedelta(seconds=timeout)
        cutoff = str(datetime.datetime.utcnow() - delta)
        # self.logger.debug("Doing timeout check,"
        #                   " cutoff=%s, and setting last check to %f",
        #                   cutoff, now)
        results = self.wuar.query(eq={"status": wudb.WuStatus.ASSIGNED},
                                  lt={"timeassigned": cutoff})
        results += self.wuar.query(eq={"status": wudb.WuStatus.NEED_RESUBMIT})
        if not results:
            # self.logger.debug("Found no timed-out workunits")
            pass
        self.logger.debug("Timeout check took %f s, found %d WUs",
                          time.time() - now, len(results))
        for entry in results:
            # XXX we would like to do this atomically, but the current db
            # access methods ignore the commit argument, which is a huge
            # pain. In fact, it's not that much of a problem, though. We
            # can live with the wu disappearing for a jiffy or two.
            self.cancel_wu(entry["wuid"], commit=False)
            self.resubmit_one_wu(Workunit(entry["wu"]), commit=True)

    def handle_error_result(self, message):
        """
        Handle workunit with non-zero exit code

        If the result message indicates a failed command, log an error
        message, set the workunit to VERIFIED_ERROR in the DB, resubmit
        the work unit (but no more than once) and return True.
        If it indicates no error, return False.
        """
        if message.get_exitcode(0) == 0:
            return False
        self.log_failed_command_error(message, 0)
        key = "wu_failed"
        maxfailed = self.params["maxfailed"]
        maxwuerror = self.params["maxwuerror"]
        if not self.state[key] < maxfailed:
            self.logger.error("Exceeded maximum number of failed "
                              "workunits, maxfailed=%d ", maxfailed)
            raise Exception("Too many failed work units")
        results = self.wuar.query(eq={"wuid": message.get_wu_id()})
        assert len(results) == 1  # There must be exactly 1 WU
        assert results[0]["status"] == wudb.WuStatus.RECEIVED_ERROR
        wu = workunit.Workunit(results[0]["wu"])
        self.state.update({key: self.state[key] + 1}, commit=False)
        self.verification(message.get_wu_id(), False, commit=False)
        self.resubmit_one_wu(wu, commit=True, maxresubmit=maxwuerror)
        return True


class Polysel1Task(ClientServerTask,
                   DoesImport,
                   HasStatistics,
                   patterns.Observer):
    """
    Finds a number of size-optimized polynomial, uses client/server
    """

    @property
    def name(self):
        return "polyselect1"

    @property
    def title(self):
        return "Polynomial Selection (size optimized)"

    @property
    def programs(self):
        # admin and admax are special, which is a bit ugly: these parameters
        # to the Polyselect constructor are supplied by the task, but the
        # task has itself admin, admax parameters, which specify the total
        # size of the search range. Thus we don't include admin, admax here,
        # or PolyselTask would incorrectly warn about them not being used.
        return ((cadoprograms.Polyselect, ("out"), {}),)

    @property
    def paramnames(self):
        return self.join_params(super().paramnames, {
            "N": int, "adrange": int, "admin": 0, "admax": int,
            "nrkeep": 20, "import_sopt": [str]})

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
    # Stat: total phase took 55.47s
    @property
    def stat_conversions(self):
        return (
                (
                    "stats_collisions",
                    float,
                    "0",
                    Statistics.add_list,
                    re_fp_compile(r'# Stat: potential collisions=({fp})'),
                    False
                    ),
                (
                    "stats_rawlognorm",
                    (int, float, float, float, float),
                    "0 0 0 0 0",
                    self.update_lognorms,
                    re_fp_compile(r"# Stat: raw lognorm \(nr/min/av/max/std\):"
                                  r" (\d+)/({fp})/({fp})/({fp})/({fp})"),
                    False
                    ),
                (
                    "stats_optlognorm",
                    (int, float, float, float, float),
                    "0 0 0 0 0",
                    self.update_lognorms,
                    re_fp_compile(r"# Stat: optimized lognorm"
                                  r" \(nr/min/av/max/std\):"
                                  r" (\d+)/({fp})/({fp})/({fp})/({fp})"),
                    False
                    ),
                (
                    "stats_total_time",
                    float,
                    "0",
                    Statistics.add_list,
                    re_fp_compile(r'# Stat: total phase took ({fp})s'),
                    False
                    ),
                )

    @property
    def stat_formats(self):
        L = "lognorm (nr/min/av/max/std)"
        return (
            ["potential collisions: {stats_collisions[0]:g}"],
            [f"raw {L}" + ": {stats_rawlognorm[0]:d}"] +
            ["/{stats_rawlognorm[%d]:.3f}" % i for i in range(1, 5)],
            [f"optimized {L}" + ": {stats_optlognorm[0]:d}"] +
            ["/{stats_optlognorm[%d]:.3f}" % i for i in range(1, 5)],
            ["Total time: {stats_total_time[0]:g}"],
            )

    def get_statistics_as_strings(self):
        # technically, polyselect1 does not import anything: it just
        # doesn't run. So we can't check self.did_import, as it will
        # remain false. The (final) import thing happens in polyselect2,
        # while import that happens here is not exclusive with the fact
        # of actually running.
        if self.send_request(Request.GET_WILL_IMPORT_FINAL_POLYNOMIAL):
            return [], None
        else:
            return super().get_statistics_as_strings()

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        assert self.params["nrkeep"] > 0
        self.state["adnext"] = \
            max(self.state.get("adnext", 0), self.params["admin"])
        # Remove admin and admax from the parameter-file-supplied program
        # parameters as those would conflict with the computed values
        self.progparams[0].pop("admin", None)
        self.progparams[0].pop("admax", None)

        tablename = self.make_tablename("bestpolynomials")
        self.best_polynomials = self.make_db_dict(
                tablename, connection=self.db_connection)
        self._check_best_polynomials()

        self.poly_heap = []
        # If we have "import", discard any existing polynomials
        if "import" in self.params and self.best_polynomials:
            self.logger.warning('Have "import" parameter, discarding '
                                'previously found polynomials')
            self.best_polynomials.clear()
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
        """
        Compare that the polynomials in the heap and in the DB agree

        They must contain an equal number of entries, and each polynomial
        stored in the heap must be at the specified index in the DB.
        """
        assert len(self.poly_heap) == len(self.best_polynomials)
        for exp_E, (key, poly) in self.poly_heap:
            assert self.best_polynomials[key] == str(poly)

    def import_existing_polynomials(self):
        debug = False
        oldkeys = list(self.best_polynomials.keys())
        oldkeys.sort(key=int)  # Sort by numerical value
        for oldkey in oldkeys:
            if debug:
                print("Adding old polynomial at DB index %s: %s" %
                      (oldkey, self.best_polynomials[oldkey]))
            poly = Polynomials(self.best_polynomials[oldkey].splitlines())
            if not poly.exp_E:
                self.logger.error("Polynomial at DB index %s has no exp_E",
                                  oldkey)
                continue
            newkey = self._add_poly_heap(poly)
            if newkey is None:
                # Heap is full, and the poly was worse than the worst one on
                # the heap. Thus it did not get added and must be removed from
                # the DB
                if debug:
                    print("Deleting polynomial exp_E=%f, key=%s" %
                          (poly.exp_E, oldkey))
                del self.best_polynomials[oldkey]
            elif newkey != oldkey:
                # Heap is full, worst one in heap (with key=newkey) was
                # overwritten and its DB entry gets replaced with poly from
                # key=oldkey
                if debug:
                    print("Overwriting poly exp_E=%f, key=%s with poly "
                          "exp_E=%f, key=%s" %
                          (self.poly_heap[0][0], newkey, poly, oldkey))
                self.best_polynomials.clear(oldkey, commit=False)
                self.best_polynomials.update({newkey: poly}, commit=True)
            else:
                # Last case newkey == oldkey: nothing to do
                if debug:
                    print("Adding exp_E=%f, key=%s" % (poly.exp_E, oldkey))

    def run(self):
        if self.send_request(Request.GET_WILL_IMPORT_FINAL_POLYNOMIAL):
            self.logger.info("Skipping this phase,"
                             " as we will import the final polynomial")
            return True

        super().run()

        self.do_import()
        if self.did_import() and "import_sopt" in self.params:
            self.logger.critical("The import and import_sopt parameters "
                                 "are mutually exclusive")
            return False

        if self.did_import():
            self.logger.info("Imported polynomial(s), skipping this phase")
            return True

        if "import_sopt" in self.params:
            self.import_files(self.params["import_sopt"])

        worstmsg = ""
        if self.poly_heap:
            worstmsg = ", worst exp_E %f" % -self.poly_heap[0][0]
        self.logger.info("%d polynomials in queue from previous run%s",
                         len(self.poly_heap), worstmsg)

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

    def get_achievement(self):
        # Note that wu_range_received (like wu_received, by the way)
        # counts ERROR'd workunits as well !
        adspan = self.params["admax"] - self.params["admin"]
        return self.state["wu_range_received"] / adspan

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
        """
        Return blocks of consecutive non-empty lines from input

        Whitespace is stripped; a line containing only whitespace is
        considered empty. An empty block is never returned.

        >>> list(Polysel1Task.read_blocks(['', 'a', 'b', '',
        ...                                'c', '', '', 'd', 'e', '']))
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

    def import_one_file(self, filename):
        self.process_polyfile(filename)

    def process_polyfile(self, filename, commit=True):
        """
        Read all size-optimized polynomials in a file and add them to the
        DB and priority queue if worthwhile.

        Different polynomials must be separated by a blank line.
        """
        try:
            polyfile = self.read_log_warning(filename)
        except (OSError, IOError) as e:
            if e.errno == errno.ENOENT:
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
        if have < nrkeep:
            fullmsg = "%d/%d" % (have, nrkeep)
        else:
            fullmsg = "%d" % nrkeep
        self.logger.info("Parsed %d polynomials,"
                         " added %d to priority queue (has %s)",
                         totalparsed, totaladded, fullmsg)
        if totaladded:
            self.logger.info("Worst polynomial in queue now has exp_E %f",
                             -self.poly_heap[0][0])

    def read_log_warning(self, filename):
        """
        Read lines from file. If a "# WARNING" line occurs, log it.
        """
        re_warning = re.compile("# WARNING")
        with open(filename, "r") as inputfile:
            for line in inputfile:
                if re_warning.match(line):
                    self.logger.warning("File %s contains: %s",
                                        filename, line.strip())
                yield line

    def parse_and_add_poly(self, text, filename):
        """
        Parse a polynomial from an iterable of lines and add it to the
        priority queue and DB. Return a two-element list with the number of
        polynomials parsed and added, i.e., (0,0) or (1,0) or (1,1).
        """
        poly = self.parse_poly(text, filename)
        if poly is None:
            return (0, 0)
        if poly.n != self.params["N"]:
            self.logger.error("Polynomial is for the wrong number"
                              " to be factored:\n%s",
                              poly)
            return (0, 0)
        if not poly.exp_E:
            self.logger.warning("Polynomial in file %s has no exp_E,"
                                " skipping it",
                                filename)
            return (0, 0)
        if self._add_poly_heap_db(poly):
            return (1, 1)
        else:
            return (1, 0)

    def _add_poly_heap_db(self, poly):
        """
        Add a polynomial to the heap and DB, if it's good enough.

        Returns True if the poly was added, False if not.
        """
        key = self._add_poly_heap(poly)
        if key is None:
            return False
        self.best_polynomials[key] = str(poly)
        return True

    def _add_poly_heap(self, poly):
        """
        Add a polynomial to the heap

        If the heap is full (nrkeep), the worst polynomial (i.e., with
        the largest exp_E) is replaced if the new one is better.  Returns
        the key (as a str) under which the polynomial was added, or None
        if it was not added.
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
            worst_exp_E = -self.poly_heap[0][0]
            if worst_exp_E <= poly.exp_E:
                if debug:
                    self.logger.debug("_add_poly_heap(): new poly exp_E %f,"
                                      " worst in heap has %f. Not adding",
                                      poly.exp_E, worst_exp_E)
                return None
            # Pop the worst poly from heap and re-use its DB index
            key = heapq.heappop(self.poly_heap)[1][0]
            if debug:
                self.logger.debug("_add_poly_heap(): new poly exp_E %f,"
                                  " worst in heap has %f."
                                  " Replacing DB index %s",
                                  poly.exp_E, worst_exp_E, key)
        else:
            # Heap was not full
            if debug:
                self.logger.debug("_add_poly_heap(): heap was not full,"
                                  " adding poly with exp_E %f"
                                  " at DB index %s", poly.exp_E, key)

        # The DB requires the key to be a string. In order to have
        # identical data in DB and heap, we store key as str everywhere.
        key = str(key)

        # Python heapq stores a minheap, so in order to have the worst
        # polynomial (with largest exp_E) easily accessible, we use
        # -exp_E as the heap key
        new_entry = (-poly.exp_E, (key, poly))
        heapq.heappush(self.poly_heap, new_entry)
        return key

    def parse_poly(self, text, filename):
        poly = None
        try:
            poly = Polynomials(text)
        except PolynomialParseException as e:
            if str(e) != "No polynomials found":
                self.logger.warning("Invalid polyselect file '%s': %s",
                                    filename, e)
                return None
        except UnicodeDecodeError as e:
            self.logger.error("Error reading '%s' (corrupted?): %s",
                              filename, e)
            return None

        if not poly:
            return None
        return poly

    def get_raw_polynomials(self):
        # Extract polynomials from heap and return as list
        return [entry[1][1] for entry in self.poly_heap]

    def get_poly_rank(self, search_poly):
        """
        Return how many polynomnials with exp_E less than the exp_E
        of the size-optimized version of search_poly there are in the
        priority queue.

        The size-optimized version of search_poly is identified by comparing
        the leading coefficients of both polynomials.
        """

        # Search for the raw polynomial pair by comparing the leading
        # coefficients of both polynomials
        found = None
        for (index, (exp_E, (key, poly))) in enumerate(self.poly_heap):
            if search_poly.polyg.same_lc(poly.polyg):
                if found is not None:
                    self.logger.warning("Found more than one match for:\n%s",
                                        search_poly)
                else:
                    found = index
        if found is None:
            self.logger.warning("Could not find polynomial rank for %s",
                                search_poly)
            return None
        # print("search_poly: %s" % search_poly)
        # print("Poly found in heap: %s" % self.poly_heap[found][1][1])
        search_exp_E = -self.poly_heap[found][0]
        rank = 0
        for (exp_E, (key, poly)) in self.poly_heap:
            if -exp_E < search_exp_E:
                rank += 1
        return rank

    def need_more_wus(self):
        return self.state["adnext"] < self.params["admax"]

    def submit_one_wu(self):
        adstart = self.state["adnext"]
        adend = adstart + self.params["adrange"]
        adend = adend - (adend % self.params["adrange"])
        assert adend > adstart
        adend = min(adend, self.params["admax"])
        outputfile = self.workdir.make_filename("%d-%d" % (adstart, adend),
                                                prefix=self.name)
        if self.test_outputfile_exists(outputfile):
            self.logger.info("%s already exists, won't generate again",
                             outputfile)
        else:
            p = cadoprograms.Polyselect(admin=adstart, admax=adend,
                                        stdout=str(outputfile),
                                        skip_check_binary_exists=True,
                                        **self.progparams[0])
            self.submit_command(p, "%d-%d" % (adstart, adend), commit=False)
        self.state.update({"adnext": adend}, commit=True)

    def get_total_cpu_or_real_time(self, is_cpu):
        """
        Return number of seconds of cpu time spent by polyselect
        """
        return float(self.state.get("stats_total_time", 0.)) if is_cpu else 0.


class Polysel2Task(ClientServerTask,
                   HasStatistics,
                   DoesImport,
                   patterns.Observer):
    """
    Finds a polynomial, uses client/server
    """

    @property
    def name(self):
        return "polyselect2"

    @property
    def title(self):
        return "Polynomial Selection (root optimized)"

    @property
    def programs(self):
        return ((cadoprograms.PolyselectRopt, (), {}),
                (cadoprograms.Skewness, (), {}))

    @property
    def paramnames(self):
        return self.join_params(super().paramnames, {
            "N": int,
            "I": [int],
            "A": [int],
            "qmin": int,
            "lpb1": int,
            "lpb0": int,
            "batch": [int],
            "import_ropt": [str]
            })

    @property
    def stat_conversions(self):
        return (
            (
                "stats_total_time",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r'# Stat: total phase took ({fp})s'),
                False
            ),
            (
                "stats_rootsieve_time",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r'# Stat: rootsieve took ({fp})s'),
                False
            )
            )

    @property
    def stat_formats(self):
        return (
            ["Total time: {stats_total_time[0]:g}"],
            ["Rootsieve time: {stats_rootsieve_time[0]:g}"],
            )

    def get_statistics_as_strings(self):
        if self.did_import():
            return [], None
        else:
            return super().get_statistics_as_strings()

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.bestpoly = None
        self.best_polys = []
        if "import" in self.params and "bestpoly" in self.state:
            self.logger.warning('Have "import" parameter, discarding '
                                'previously found best polynomial')
            self.state.clear(["bestpoly"])
        if "bestpoly" in self.state:
            self.bestpoly = Polynomials(self.state["bestpoly"].splitlines())
            self.best_polys = [self.bestpoly]
        self.state.setdefault("nr_poly_submitted", 0)
        # I don't understand why the area is based on one particular side.
        if "I" in self.params:
            areabits = 2*self.params["I"]-1
        elif "A" in self.params:
            areabits = self.params["A"]
        else:
            msg = "Required parameter I or A not found " \
                  "for polyselects's area, under %s ; " \
                  "consider setting tasks.I or tasks.A" % path_prefix
            self.logger.critical(msg)
            raise KeyError(msg)
        self.progparams[0].setdefault("area",
                                      2.**areabits * self.params["qmin"])
        # on Sep 26, 2018, changed Bf, Bg from lim1/lim0 to 2^lpb1/2^lpb0
        self.progparams[0].setdefault("Bf", float(2**self.params["lpb1"]))
        self.progparams[0].setdefault("Bg", float(2**self.params["lpb0"]))
        if "batch" not in self.params:
            t = self.progparams[0].get("threads", 1)
            # batch = 5 rounded up to a multiple of t
            self.params["batch"] = (4 // t + 1) * t
        self.poly_to_submit = None

    def run(self):
        super().run()

        self.do_import()
        if self.bestpoly is None:
            self.logger.info("No polynomial was previously found")
        else:
            self.logger.info("Best polynomial previously has Murphy_E = %.3e",
                             self.bestpoly.MurphyE)

        if self.did_import() and "import_ropt" in self.params:
            self.logger.critical("The import and import_ropt parameters "
                                 "are mutually exclusive")
            return False

        if self.did_import():
            self.logger.info("Imported polynomial, skipping this phase")
            return True

        if "import_ropt" in self.params:
            self.import_files(self.params["import_ropt"])

        # Get the list of polynomials to submit
        self.poly_to_submit = self.send_request(Request.GET_RAW_POLYNOMIALS)

        if self.is_done():
            self.logger.info("Already finished - nothing to do")
            self.print_rank()
            # If the poly file got lost somehow, write it again
            filename = self.get_state_filename("polyfilename")
            if filename is None or not filename.isfile():
                self.logger.warning("Polynomial file disappeared,"
                                    " writing again")
                self.write_poly_file()
            return True

        # Submit all the WUs we need
        while self.need_more_wus():
            self.submit_one_wu()

        # Wait for all the WUs to finish
        while self.get_number_outstanding_wus() > 0:
            self.wait()

        # Print best polynomials found
        if False:
            for i in range(len(self.best_polys)):
                self.logger.info(f"Best polynomial {i} is\n"
                                 f"{self.best_polys[i]}")
        else:
            n = len(self.best_polys)
            self.logger.info("Kept %d polynomials with"
                             " MurphyE from %.3e to %.3e",
                             n,
                             self.best_polys[0].MurphyE,
                             self.best_polys[n-1].MurphyE)

        # save best polynomials in cxxx.poly.0, cxxx.poly.1, ...
        # and run polyselect3
        n = len(self.best_polys)
        for i in range(n):
            filename = self.workdir.make_filename("poly." + str(i))
            with open(str(filename), "w") as poly_file:
                poly_file.write(str(self.best_polys[i]))
        # remove ropteffort since polyselect3 does not use it
        self.progparams[0].pop("ropteffort", None)
        filename = self.workdir.make_filename("poly")
        p = cadoprograms.Polyselect3(poly=filename, num=n,
                                     **self.progparams[0])
        process = cadocommand.Command(p)
        (rc, stdout, stderr) = process.wait()

        # select the polynomial pair with the largest Murphy-E
        best_i = -1
        best_MurphyF = 0.0
        for i in range(n):
            filename = self.workdir.make_filename("poly." + str(i))
            f = open(str(filename), "r")
            poly = Polynomials(f.read().splitlines())
            f.close()
            self.logger.info("Polynomial %s had MurphyE %.3e, refined to %.3e",
                             filename, poly.MurphyE, poly.MurphyF)
            if best_i == -1 or poly.MurphyF > best_MurphyF:
                best_i = i
                best_MurphyF = poly.MurphyF
                self.bestpoly = poly

        filename = self.workdir.make_filename("poly." + str(best_i))
        self.logger.info("Best polynomial is %s", filename)

        if self.bestpoly is None:
            self.logger.error("No polynomial found. Consider increasing the "
                              "search range bound admax")
            return False

        self.logger.info("Finished, best polynomial has Murphy_E = %.3e",
                         self.bestpoly.MurphyE)
        self.logger.info("Best polynomial is:\n%s", str(self.bestpoly))
        self.print_rank()
        self.write_poly_file()

        return True

    def is_done(self):
        return self.bestpoly is not None and \
                not self.need_more_wus() and \
                self.get_number_outstanding_wus() == 0

    def need_more_wus(self):
        return self.state["nr_poly_submitted"] < len(self.poly_to_submit)

    def get_achievement(self):
        departed = self.state["nr_poly_submitted"]
        in_flight = self.params["batch"] * self.get_number_outstanding_wus()
        total = len(self.poly_to_submit)
        return (departed - in_flight) / total

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
        # FIXME: wrong, we should always get an optimized poly for a raw one
        self.verification(message.get_wu_id(), True, commit=True)
        return True

    def import_one_file(self, filename):
        old_bestpoly = self.bestpoly
        self.process_polyfile(filename, allow_no_skewness=True)

        if self.bestpoly is not old_bestpoly:
            self.write_poly_file()

    @staticmethod
    def read_blocks(input):
        """
        Return blocks of consecutive non-empty lines from input

        Whitespace is stripped; a line containing only whitespace is
        considered empty. An empty block is never returned.

        >>> list(Polysel1Task.read_blocks(['', 'a', 'b', '',
        ...                                'c', '', '', 'd', 'e', '']))
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

    def process_polyfile(self, filename, *,
                         allow_no_skewness=False,
                         commit=True):
        try:
            polyfile = self.read_log_warning(filename)
        except (OSError, IOError) as e:
            if e.errno == errno.ENOENT:
                self.logger.error("File '%s' does not exist", filename)
                return None
            else:
                raise
        for block in self.read_blocks(polyfile):
            poly = self.parse_poly(block, filename,
                                   allow_no_skewness=allow_no_skewness)
            if poly is not None:
                self.bestpoly = poly
                update = {"bestpoly": str(poly)}
                self.state.update(update, commit=commit)

    def read_log_warning(self, filename):
        """
        Read lines from file. If a "# WARNING" line occurs, log it.
        """
        re_warning = re.compile("# WARNING")
        with open(filename, "r") as inputfile:
            for line in inputfile:
                if re_warning.match(line):
                    self.logger.warning("File %s contains: %s",
                                        filename, line.strip())
                yield line

    # If allow_no_skewness is True and the polynomial file does not contain
    # a skew: line, compute the skewness and use it.
    def parse_poly(self, text, filename, *, allow_no_skewness=False):
        poly = None
        try:
            poly = Polynomials(text)
        except (OSError, IOError) as e:
            if e.errno == errno.ENOENT:
                self.logger.error("File '%s' does not exist", filename)
                return None
            else:
                raise
        except PolynomialParseException as e:
            if str(e) != "No polynomials found":
                self.logger.warning("Invalid polyselect file '%s': %s",
                                    filename, e)
                return None
        except UnicodeDecodeError as e:
            self.logger.error("Error reading '%s' (corrupted?): %s",
                              filename, e)
            return None

        if not poly:
            # This happens when all polynomials are parsed
            return None
        if poly.n != self.params["N"]:
            self.logger.error("Polynomial is for the wrong number"
                              " to be factored:\n%s",
                              poly)
            return None
        if not poly.MurphyE:
            self.logger.info("Polynomial in file %s has no Murphy E value, "
                             "computing it", filename)
            p = cadoprograms.Score(inputpoly=filename,
                                   **self.merged_args[1])
            process = cadocommand.Command(p)
            (rc, stdout, stderr) = process.wait()
            if rc != 0:
                self.logger.error("Computing murphyE failed with exit code %d",
                                  rc)
            poly.MurphyE = float(stdout)
            self.logger.info("Computed murphyE is %.5g", poly.MurphyE)
        if poly.skew is None:
            if not allow_no_skewness:
                self.logger.error("Polynomial in file %s has no skew value",
                                  filename)
                return None
            else:
                self.logger.info("Polynomial in file %s has no skew value, "
                                 "computing it", filename)
                p = cadoprograms.Skewness(inputpoly=filename,
                                          **self.merged_args[1])
                process = cadocommand.Command(p)
                (rc, stdout, stderr) = process.wait()
                if rc != 0:
                    self.logger.error("Computing skewness failed"
                                      " with exit code %d",
                                      rc)
                    return None
                poly.skew = float(stdout)
                self.logger.info("Computed skewness is %.5g", poly.skew)

        margin = 0.80  # we keep polys with MurphyE >= margin*bestMurphyE
        if not self.best_polys or \
                poly.MurphyE >= margin * self.bestpoly.MurphyE:
            self.best_polys.append(poly)

        i = len(self.best_polys) - 1
        while i > 0 and \
                self.best_polys[i].MurphyE > self.best_polys[i-1].MurphyE:
            t = self.best_polys[i]
            self.best_polys[i] = self.best_polys[i-1]
            self.best_polys[i-1] = t
            i = i - 1

        # remove polynomials with MurphyE < margin*bestMurphyE
        i = len(self.best_polys) - 1
        while self.best_polys[i].MurphyE < margin * self.best_polys[0].MurphyE:
            del self.best_polys[i]
            i = i - 1
        # in case poly.MurphyE = self.bestpoly.MurphyE (MurphyE is printed
        # only with 3 digits in the cxxx.poly file), we choose the polynomial
        # with the smallest skewness, to avoid non-determinism in polynomial
        # selection. Indeed, if we have two polynomials with the same MurphyE
        # value (when rounded on 3 digits), the one found first was chosen
        # before that change.
        if self.bestpoly is None \
           or poly.MurphyE > self.bestpoly.MurphyE \
           or (poly.MurphyE == self.bestpoly.MurphyE
               and poly.skew < self.bestpoly.skew):
            self.logger.info("New best polynomial from file %s:"
                             " Murphy E = %.3e" % (filename, poly.MurphyE))
            self.logger.debug("New best polynomial is:\n%s", poly)
            return poly
        else:
            self.logger.info("Best polynomial from file %s with E=%.3e is "
                             "no better than current best with E=%.3e",
                             filename, poly.MurphyE, self.bestpoly.MurphyE)
        return None

    def write_poly_file(self):
        filename = self.workdir.make_filename("poly")
        with open(str(filename), "w") as poly_file:
            poly_file.write(str(self.bestpoly))
        self.state["polyfilename"] = filename.get_wdir_relative()

    def get_poly(self):
        if "bestpoly" not in self.state:
            return None
        return Polynomials(self.state["bestpoly"].splitlines())

    def get_poly_filename(self):
        return self.get_state_filename("polyfilename")

    def submit_one_wu(self):
        assert self.need_more_wus()
        to_submit = len(self.poly_to_submit)
        nr = self.state["nr_poly_submitted"]
        inputfilename = self.workdir.make_filename("raw_%d" % nr,
                                                   prefix=self.name)
        # Write one raw polynomial to inputfile
        batchsize = min(to_submit - nr, self.params["batch"])
        with inputfilename.open("w") as inputfile:
            for i in range(batchsize):
                inputfile.write(str(self.poly_to_submit[nr + i]))
                inputfile.write("\n")
        outputfile = self.workdir.make_filename("opt_%d" % nr,
                                                prefix=self.name)
        if self.test_outputfile_exists(outputfile):
            self.logger.info("%s already exists, won't generate again",
                             outputfile)
        else:
            p = cadoprograms.PolyselectRopt(inputpolys=str(inputfilename),
                                            stdout=str(outputfile),
                                            skip_check_binary_exists=True,
                                            **self.progparams[0])
            self.submit_command(p, "%d" % nr, commit=False)
        self.state.update({"nr_poly_submitted": nr + batchsize}, commit=True)

    def get_total_cpu_or_real_time(self, is_cpu):
        """
        Return number of seconds of cpu time spent by polyselect_ropt
        """
        return float(self.state.get("stats_total_time", 0.)) if is_cpu else 0.

    def print_rank(self):
        if not self.did_import():
            rank = self.send_request(Request.GET_POLY_RANK, self.bestpoly)
            if rank is not None:
                self.logger.info("Best overall polynomial was %d-th in list "
                                 "after size optimization", rank)

    def get_will_import(self):
        return "import" in self.params


# TODO: add HasStatistics
class PolyselJLTask(ClientServerTask, DoesImport, patterns.Observer):
    """
    Find a polynomial pair using Joux-Lercier for DL in GF(p), uses
    client/server
    """

    @property
    def name(self):
        return "polyselectJL"

    @property
    def title(self):
        return "Polynomial Selection (Joux-Lercier)"

    @property
    def programs(self):
        return ((cadoprograms.PolyselectJL, (), {}),)

    @property
    def paramnames(self):
        return self.join_params(super().paramnames, {
            "N": int, "modr": 0, "nrkeep": 20,
            "bound": int, "modm": int, "degree": int,
            "I": [int],
            "A": [int],
            "lim1": int, "lim0": int,
            "lpb0": int, "lpb1": int,
            "qmin": 0,
            "ell": int,
            "fastSM": False
            })

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        if "import" in self.params:
            self.logger.warning('Have "import" parameter')
        assert self.params["nrkeep"] > 0
        self.state["rnext"] = self.state.get("rnext", 0)
        qmin = self.params["qmin"]
        if qmin == 0:
            qmin = max(self.params["lim0"], self.params["lim1"])
        if "I" in self.params:
            areabits = 2*self.params["I"]-1
        elif "A" in self.params:
            areabits = self.params["A"]
        else:
            msg = "Required parameter I or A not found " \
                  "for polyselects's area, under %s ; " \
                  "consider setting tasks.I or tasks.A" % path_prefix
            self.logger.critical(msg)
            raise KeyError(msg)
        self.progparams[0].setdefault("area", 2.**areabits * qmin)
        self.progparams[0].setdefault("Bf", float(2**self.params["lpb1"]))
        self.progparams[0].setdefault("Bg", float(2**self.params["lpb0"]))
        if self.params["fastSM"]:
            self.progparams[0].setdefault("easySM", self.params["ell"])
        self.bestpoly = None
        if "bestpoly" in self.state:
            self.bestpoly = Polynomials(self.state["bestpoly"].splitlines())

    def run(self):
        super().run()

        self.do_import()
        if "import" not in self.params:
            if self.is_done():
                self.logger.info("Already finished - nothing to do")
                return True

            # Submit all the WUs we need to reach the final modr
            while self.need_more_wus():
                self.submit_one_wu()

            # Wait for all the WUs to finish
            while self.get_number_outstanding_wus() > 0:
                self.wait()

            self.logger.info("Finished")

        filename = self.workdir.make_filename("poly")
        with open(str(filename), "w") as poly_file:
            poly_file.write(str(self.bestpoly))
        self.state["polyfilename"] = filename.get_wdir_relative()
        self.state["bestpoly"] = str(self.bestpoly)
        self.logger.info("Selected polynomial has MurphyE %f",
                         self.bestpoly.MurphyE)
        return True

    def is_done(self):
        return not self.need_more_wus() and \
            self.get_number_outstanding_wus() == 0

    def get_achievement(self):
        # Note that wu_range_received (like wu_received, by the way)
        # counts ERROR'd workunits as well !
        return self.state["wu_range_received"] / self.params["modm"]

    def updateObserver(self, message):
        identifier = self.filter_notification(message)
        if not identifier:
            # This notification was not for me
            return False
        if self.handle_error_result(message):
            return True
        (filename, ) = message.get_output_files()
        self.process_polyfile(filename, commit=False)
        # Always mark ok to avoid warning messages about WUs that did not
        # find a poly
        self.verification(message.get_wu_id(), True, commit=True)
        return True

    @staticmethod
    def read_blocks(input):
        """
        Return blocks of consecutive non-empty lines from input

        Whitespace is stripped; a line containing only whitespace is
        considered empty. An empty block is never returned.

        >>> list(Polysel1Task.read_blocks(['', 'a', 'b', '',
        ...                                'c', '', '', 'd', 'e', '']))
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

    def import_one_file(self, filename):
        self.process_polyfile(filename)

    def process_polyfile(self, filename, commit=True):
        try:
            polyfile = self.read_log_warning(filename)
        except (OSError, IOError) as e:
            if e.errno == errno.ENOENT:
                self.logger.error("File '%s' does not exist", filename)
                return None
            else:
                raise

        totalparsed = 0
        for block in self.read_blocks(polyfile):
            parsed = self.parse_and_add_poly(block, filename)
            totalparsed += parsed
        self.logger.info("Parsed %d polynomials", totalparsed)

    def read_log_warning(self, filename):
        """
        Read lines from file. If a "# WARNING" line occurs, log it.
        """
        re_warning = re.compile("# WARNING")
        with open(filename, "r") as inputfile:
            for line in inputfile:
                if re_warning.match(line):
                    self.logger.warning("File %s contains: %s",
                                        filename, line.strip())
                yield line

    def parse_and_add_poly(self, text, filename):
        poly = self.parse_poly(text, filename)
        if poly is None:
            return 0
        if poly.n != self.params["N"]:
            self.logger.error("Polynomial is for the wrong prime:\n%s",
                              poly)
            return 0
        if not poly.MurphyE:
            self.logger.warning("Polynomial in file %s has no MurphyE,"
                                " skipping it",
                                filename)
            return 0
        if self.bestpoly is None \
           or self.bestpoly.MurphyE < poly.MurphyE \
           or (self.bestpoly.MurphyE == poly.MurphyE
               and self.bestpoly.skew > poly.skew):
            self.bestpoly = poly
            self.logger.info("Best polynomial so far has MurphyE %f",
                             poly.MurphyE)
        return 1

    def parse_poly(self, text, filename):
        poly = None
        try:
            poly = Polynomials(text)
        except PolynomialParseException as e:
            if str(e) != "No polynomials found":
                self.logger.warning("Invalid polyselect file '%s': %s",
                                    filename, e)
                return None
        except UnicodeDecodeError as e:
            self.logger.error("Error reading '%s' (corrupted?): %s",
                              filename, e)
            return None

        if not poly:
            return None
        return poly

    def get_poly_filename(self):
        return self.get_state_filename("polyfilename")

    def get_poly(self):
        if "bestpoly" not in self.state:
            return None
        return Polynomials(self.state["bestpoly"].splitlines())

    def need_more_wus(self):
        return 1 + self.state["rnext"] < self.params["modm"]

    def submit_one_wu(self):
        modr = self.state["rnext"]
        df = self.params["degree"]
        dg = self.params["degree"]-1
        outputfile = self.workdir.make_filename("%d" % (modr,),
                                                prefix=self.name)
        if self.test_outputfile_exists(outputfile):
            self.logger.info("%s already exists, won't generate again",
                             outputfile)
        else:
            p = cadoprograms.PolyselectJL(modr=modr, df=df, dg=dg,
                                          stdout=str(outputfile),
                                          skip_check_binary_exists=True,
                                          **self.progparams[0])
            self.submit_command(p, "%d" % (modr,), commit=False)
        self.state.update({"rnext": modr+1}, commit=True)


# TODO: add HasStatistics
class PolyselGFpnTask(Task, DoesImport):
    """
    Polynomial selection for DL in extension fields
    """

    @property
    def name(self):
        return "polyselgfpn"

    @property
    def title(self):
        return "Polynomial Selection (for GF(p^n))"

    @property
    def programs(self):
        override = ("p", "n", "out")
        return ((cadoprograms.PolyselectGFpn, (override), {}),)

    @property
    def paramnames(self):
        return self.join_params(super().paramnames,
                                {"N": int,
                                 "gfpext": int,
                                 "import": None})

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)

    def run(self):
        super().run()

        self.do_import()
        if "polyfilename" not in self.state:
            polyfilename = self.workdir.make_filename("poly")
            # Import mode
            if self.did_import():
                if not self.state["imported_poly"]:
                    raise Exception("Import failed?")
                with open(str(polyfilename), "w") as outfile:
                    outfile.write(self.state["imported_poly"])
                update = {
                        "poly": self.state["imported_poly"],
                        "polyfilename": polyfilename.get_wdir_relative()
                        }
                self.state.update(update)
                return True

            # Check that user does not ask for degree > 2
            if self.params["gfpext"] != 2:
                raise Exception("Polynomial selection not implemented "
                                + "for extension degree > 2")

            # Call binary
            (stdoutpath, stderrpath) = self.make_std_paths(
                cadoprograms.PolyselectGFpn.name)
            p = cadoprograms.PolyselectGFpn(p=self.params["N"],
                                            n=self.params["gfpext"],
                                            out=polyfilename,
                                            stdout=str(stdoutpath),
                                            stderr=str(stderrpath),
                                            **self.merged_args[0])
            message = self.submit_command(p, None, log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            with open(str(polyfilename), "r") as inputfile:
                poly = Polynomials(list(inputfile))
            update = {"poly": str(poly),
                      "polyfilename": polyfilename.get_wdir_relative()
                      }
            self.state.update(update)
        return True

    def import_one_file(self, filename):
        with open(filename, "r") as inputfile:
            poly = Polynomials(list(inputfile))
        update = {"imported_poly": str(poly)}
        self.state.update(update)

    def get_poly(self):
        return Polynomials(self.state["poly"].splitlines())

    def get_poly_filename(self):
        return self.get_state_filename("polyfilename")


class PolyselQSTask(Task):
    """ Fake polynomial selection for QS algorithm """
    @property
    def name(self):
        return "polysel"

    @property
    def title(self):
        return "Polynomial Selection (for QS algorithm)"

    @property
    def programs(self):
        return ()

    @property
    def paramnames(self):
        return self.join_params(super().paramnames,
                                {"N": int, "computation": str, "import": None})

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)

    def run(self):
        super().run()
        if "import" in self.params:
            self.logger.critical("It is not possible to import polynomial for "
                                 "the Quadratic Sieve algorithm")
            return False

        if "polyfilename" not in self.state:
            polyfilename = self.workdir.make_filename("poly")
            # Create the polynomial file.
            N = self.params["N"]
            if self.params["computation"] == Computation.CL:
                if N % 4 == 0:
                    N = N//4
                    if N % 4 not in (2, 3):
                        raise ValueError("Discriminant should be 1, 2 or 3 "
                                         "modulo 4 or divisible by 4 with "
                                         "quotient equal to 2 or 3 modulo 4")

            poly = str(Polynomials([f"n: {N}",
                                    f"skew: {int(sqrt(abs(N)))}",
                                    f"poly0:{-N},0,1"],
                                   allow_only_one_poly=True))

            with open(str(polyfilename), "w") as outfile:
                outfile.write(poly)
            update = {
                    "poly": poly,
                    "polyfilename": polyfilename.get_wdir_relative()
                    }
            self.state.update(update)
        return True

    def get_poly(self):
        if "poly" not in self.state:
            return None
        return Polynomials(self.state["poly"].splitlines(),
                           allow_only_one_poly=True)

    def get_poly_filename(self):
        return self.get_state_filename("polyfilename")


class FactorBaseTask(Task):
    """
    Generates the factor base for the polynomial(s)
    """

    @property
    def name(self):
        return "factorbase"

    @property
    def title(self):
        return "Generate Factor Base"

    @property
    def programs(self):
        return ((cadoprograms.MakeFB,
                 ("out", "side", "lim"),
                 {"poly": Request.GET_POLYNOMIAL_FILENAME}),)

    @property
    def paramnames(self):
        lims = {f"lim{side}": int for side in range(self._nsides)}
        return self.join_params(super().paramnames,
                                {"gzip": True,
                                 "algo": str,
                                 "I": [int],
                                 "A": [int],
                                 **lims})

    def __init__(self, nsides, *, mediator, db, parameters, path_prefix):
        self._nsides = int(nsides)
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        if "I" in self.params:
            self.progparams[0].setdefault("maxbits", self.params["I"])
        elif "A" in self.params:
            self.progparams[0].setdefault("maxbits", (self.params["A"]+1)//2)
        else:
            path = ".".join(self.parameters.get_param_path())
            msg = "Required parameter I or A not found " \
                  "for makefb's maxbits under %s ; " \
                  "consider setting tasks.I or tasks.A" % path
            self.logger.critical(msg)
            raise KeyError(msg)

    def run(self):
        super().run()

        # Get best polynomial found by polyselect
        poly = self.send_request(Request.GET_POLYNOMIAL)
        if not poly:
            raise Exception("FactorBaseTask(): no polynomial "
                            "received from PolyselTask")
        if poly.nsides != self._nsides:
            raise Exception("FactorBaseTask(): expecting a polynomial with "
                            f"{self._nsides} side(s), got one with "
                            f"{poly.nsides} side(s) instead")

        # If there is already a poly in state but it is different from the
        # current one, we discard all previous computation (by cleaning state)
        if "poly" in self.state:
            one_sided = self.params["algo"] == Algorithm.QS
            prevpoly = Polynomials(self.state["poly"].splitlines(),
                                   allow_only_one_poly=one_sided)
            if poly != prevpoly:
                self.logger.warning("Received different polynomial, "
                                    "discarding old factor base file")
                # remove previous state + store new poly
                self.state = {"poly": str(poly)}
        else:
            self.state = {"poly": str(poly)}

        # For each side corresponding to a nonlinear polynomial, we compute the
        # factor base if
        #   - the output file does not exist;
        #   - the parameter lim has changed.
        if not self.has_all_outputfiles_in_state(poly):
            # Make file name for factor base/free relations file
            # We use .gzip by default, unless set to no in parameters
            use_gz = ".gz" if self.params["gzip"] else ""

            # Run command to generate factor base file
            (stdoutpath, stderrpath) = self.make_std_paths(
                cadoprograms.MakeFB.name)
            for side in poly.get_all_nonlinear_sides():
                keylim = f"lim{side}"
                if self.state.get(keylim, None) != self.params[keylim]:
                    if keylim in self.state:  # params was changed
                        self.logger.warning("Parameter %s changed,"
                                            " discarding old"
                                            " factor base file", keylim)
                    outputfilename = self.workdir.make_filename(
                                                    f"roots{side}" + use_gz)
                    p = cadoprograms.MakeFB(out=str(outputfilename),
                                            side=side,
                                            lim=self.params[keylim],
                                            stdout=str(stdoutpath),
                                            stderr=str(stderrpath),
                                            **self.merged_args[0])
                    message = self.submit_command(p, None, log_errors=True)
                    if message.get_exitcode(0) != 0:
                        raise Exception("Program failed")

                    self.state.update({
                        f"outputfile{side}":
                            outputfilename.get_wdir_relative(),
                        keylim: self.params[keylim],
                    })
                else:
                    assert f"outputfile{side}" in self.state

            self.logger.info("Finished")

        for side in poly.get_all_nonlinear_sides():
            self.check_files_exist([self.get_filename(side)],
                                   "output",
                                   shouldexist=True)
        return True

    def has_all_outputfiles_in_state(self, poly):
        return all(f"outputfile{side}" in self.state
                   for side in poly.get_all_nonlinear_sides())

    def get_filename(self, side):
        return self.get_state_filename(f"outputfile{side}")


class CheckDiscriminantTask(Task):
    """
    Check that the discriminant has no square factor that belongs to the factor
    base.
    """

    @property
    def name(self):
        return "checkdisc"

    @property
    def title(self):
        return "Check Discriminant"

    @property
    def programs(self):
        return ()

    @property
    def paramnames(self):
        return self.join_params(super().paramnames,
                                {"N": int, "computation": str, "gzip": True})

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        if self.params["computation"] != Computation.CL:
            msg = "CheckDiscriminantTask(): this task is only for class "\
                  "groups computation"
            self.logger.critical(msg)
            raise Exception(msg)

    def run(self):
        super().run()

        renumfile = self.send_request(Request.GET_RENUMBER_FILENAME)
        if not renumfile:
            raise Exception("CheckDiscriminantTask(): no renumber file "
                            "received from FreeRelTask")

        open_fun = gzip.open if self.params["gzip"] else open

        with open_fun(str(renumfile), "rb") as f:
            for line in f:
                if line.startswith(b"#"):
                    continue
                m = re.fullmatch(rb'(\d+) 0\n', line)
                if m is not None:
                    p = int(m.group(1))
                    if p > 2 and self.params["N"] % (p*p) == 0:
                        msg = f"The discriminant {self.params['N']} has a " \
                              f"square factor {p}^2 belonging to the factor " \
                              "base"
                        self.logger.critical(msg)
                        return False
                    # case p=2 is check in __init__ of CompleteFactorization
        return True


class FreeRelTask(Task):
    """
    Generates free relations for the polynomial(s)
    """

    @property
    def name(self):
        return "freerel"

    @property
    def title(self):
        return "Generate Free Relations"

    @property
    def programs(self):
        input = {"poly": Request.GET_POLYNOMIAL_FILENAME}
        return ((cadoprograms.FreeRel, ("renumber", "out", "dl"), input),)

    @property
    def paramnames(self):
        return self.join_params(super().paramnames,
                                {"computation": str,
                                 "algo": str,
                                 "gzip": True})

    wanted_regex = {
        'nfree': (r'# Free relations: (\d+)', int),
        'nprimes': (r'Renumbering struct: nprimes=(\d+)', int)
    }

    def __init__(self, nsides, *, mediator, db, parameters, path_prefix):
        self._nsides = int(nsides)
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.progparams[0]["dl"] = \
            self.params["computation"] in (Computation.DLP, Computation.CL)
        # Invariant: if we have a result (in self.state["freerelfilename"])
        # then we must also have a polynomial (in self.state["poly"]) and
        # the lpb values used in self.state[f"lpb{side}"]
        if "freerelfilename" in self.state:
            assert "poly" in self.state
            for side in range(self._nsides):
                assert f"lpb{side}" in self.state, f"lpb{side} not in state"
            # The target file must correspond to the polynomial "poly"

    def run(self):
        super().run()

        # Get best polynomial found by polyselect
        poly = self.send_request(Request.GET_POLYNOMIAL)
        if not poly:
            raise Exception("FreeRelTask(): no polynomial "
                            "received from PolyselTask")

        if poly.nsides != self._nsides:
            raise Exception("FreeRelTask(): expecting a polynomial with "
                            f"{self._nsides} side(s), got one with "
                            f"{poly.nsides} side(s) instead")

        # Check if we have already computed the freerelfile for this polynomial
        # and lpb. If any of the inputs mismatch, we remove freerelfilename
        # from state
        if "freerelfilename" in self.state:
            discard = False
            one_sided = self.params["algo"] == Algorithm.QS
            prevpoly = Polynomials(self.state["poly"].splitlines(),
                                   allow_only_one_poly=one_sided)
            if poly != prevpoly:
                self.logger.warning("Received different polynomial,"
                                    " discarding old free relations file")
                discard = True
            else:
                for side in range(self._nsides):
                    keylpb = f"lpb{side}"
                    if self.state[keylpb] != self.progparams[0][keylpb]:
                        self.logger.warning(f"Parameter {keylpb} changed,"
                                            " discarding old free relations"
                                            " file")
                        discard = True
                        break

            if discard:
                del self.state["freerelfilename"]
                del self.state["renumberfilename"]
        # If outputfile is not in state, because we never produced it or
        # because input parameters changed, we remember our current input
        # parameters
        if "freerelfilename" not in self.state:
            lpbs = {f"lpb{side}": self.progparams[0][f"lpb{side}"]
                    for side in range(self._nsides)}
            self.state.update({"poly": str(poly), **lpbs})

        if "freerelfilename" not in self.state or self.have_new_input_files():
            # Make file name for factor base/free relations file
            # We use .gzip by default, unless set to no in parameters
            use_gz = ".gz" if self.params["gzip"] else ""
            freerelfilename = self.workdir.make_filename("freerel" + use_gz)
            renumberfilename = self.workdir.make_filename("renumber" + use_gz)
            (stdoutpath, stderrpath) = self.make_std_paths(
                cadoprograms.FreeRel.name)
            # Run command to generate factor base/free relations file
            p = cadoprograms.FreeRel(renumber=renumberfilename,
                                     out=str(freerelfilename),
                                     stdout=str(stdoutpath),
                                     stderr=str(stderrpath),
                                     **self.merged_args[0])
            message = self.submit_command(p, None, log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            stderr = message.read_stderr(0).decode()
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
            if key not in found:
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


class SievingTask(ClientServerTask, DoesImport, FilesCreator, HasStatistics):
    """
    Does the sieving, uses client/server
    """

    @property
    def name(self):
        return "sieving"

    @property
    def title(self):
        return "Lattice Sieving"

    @property
    def programs(self):
        override = ("q0", "q1", "factorbase0", "factorbase1", "out",
                    "stats_stderr")
        input = {"poly": Request.GET_POLYNOMIAL_FILENAME}
        return ((cadoprograms.Las, override, input),)

    @property
    def paramnames(self):
        lims = {f"lim{side}": int for side in range(self._nsides)}
        return self.join_params(super().paramnames, {
            "qmin": 0, "qmax": [int], "qrange": int, "rels_wanted": 0,
            "gzip": True, "sqside": 1 if self._nsides > 1 else 0,
            "adjust_strategy": 0, **lims})

    def combine_bkmult(*lists):
        d = {}
        dv = 0
        for ell in lists:
            assert len(ell) == 1
            p = ell[0].split(",")
            s = 0
            if p[0].find(":") < 0:
                f = float(p[0])
                if f > dv:
                    dv = f
                s = 1
            for pp in p[s:]:
                r = re.compile(r"(\d+[a-z]+):(.*)")
                m = r.match(pp)
                assert m is not None
                assert len(m.groups()) == 2
                key, m = m.groups()
                m = float(m)
                if key not in d or m > d[key]:
                    d[key] = m
        return [",".join([str(dv)] + ["%s:%f" % v for v in d.items()])]

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
                re_fp_compile(r"# Average J=({fp})\s*for (\d+) special-q's"),
                False
            ),
            (
                "stats_max_bucket_fill",
                (str,),
                "0",
                SievingTask.combine_bkmult,
                re.compile(r"#.*max bucket fill -bkmult "
                           r"([\d\.]+(?:,\d[sl]:[\d\.]+)*)"),
                False
            ),
            (
                "stats_total_cpu_time",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r'# Total cpu time ({fp})s'),
                False
            ),
            (
                "stats_total_time",
                (float, ),
                "0",
                Statistics.add_list,
                re_fp_compile(r'# Total elapsed time ({fp})s'),
                False
            )
        )

    @property
    def stat_formats(self):
        return (
            ["Average J: {stats_avg_J[0]:g} for {stats_avg_J[1]:d} special-q",
             ", max bucket fill -bkmult {stats_max_bucket_fill[0]}"],
            ["Total time: {stats_total_time[0]:g}s"],
        )

    def __init__(self, nsides, *, mediator, db, parameters, path_prefix):
        self._nsides = int(nsides)
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        qmin = self.params["qmin"]
        if "qnext" in self.state:
            self.state["qnext"] = max(self.state["qnext"], qmin)
        else:
            # qmin = 0 is a magic value (undefined)
            if qmin > 0:
                self.state["qnext"] = qmin
            else:
                sqside = self.params["sqside"]
                self.state["qnext"] = int(self.params[f"lim{sqside}"]/2)

        self.state.setdefault("rels_found", 0)
        self.logger.info("param rels_wanted is %d", self.params["rels_wanted"])
        self.state["rels_wanted"] = self.params["rels_wanted"]
        if self.state["rels_wanted"] == 0:
            # For rsa140 with lpb[01]=29, using guess_factor = 0.91
            # gives an initial excess of about 13%, we estimate
            # guess_factor = 0.85 would be enough to get a square matrix.
            # We use 0.81 which might yield some extra filtering attempts
            # in some cases, but which should avoid too much initial excess
            # (the cado-nfs linear algebra is so fast that we should avoid
            # oversieving if we want to optimize the total wall-clock time,
            # assuming we use as many cores for sieving and linear algebra).
            guess_factor = 0.81
            s = 0.0
            for side in range(self._nsides):
                n = 2 ** self.progparams[0][f"lpb{side}"]
                s += n / log(n)
            self.state["rels_wanted"] = int(guess_factor * s)

    def enough_work_received(self):
        return self.get_nrels() >= self.state["rels_wanted"]

    def reached_qmax(self):
        if "qmax" not in self.params:
            return False
        qnext = self.state["qnext"]
        qmax = self.params["qmax"]
        return qnext >= qmax

    def should_schedule_more_work(self):
        return not self.enough_work_received() and not self.reached_qmax()

    def run(self):
        super().run()
        self.do_import()
        poly = self.send_request(Request.GET_POLYNOMIAL)

        if not poly:
            raise Exception("SievingTask(): no polynomial "
                            "received from PolyselTask")

        if poly.nsides != self._nsides:
            raise Exception("SievingTask(): expecting a polynomial with "
                            f"{self._nsides} side(s), got one with "
                            f"{poly.nsides} side(s) instead")
        fb = {}
        for side in poly.get_all_nonlinear_sides():
            fb[f"factorbase{side}"] = self.send_request(
                                            Request.GET_FACTORBASE_FILENAME,
                                            side)

        self.logger.info("We want %d relation(s)", self.state["rels_wanted"])
        qrange = self.params["qrange"]
        while self.should_schedule_more_work():
            q0 = self.state["qnext"]
            q1 = q0 + qrange
            q1 = q1 - (q1 % qrange)
            assert q1 > q0
            # We use .gzip by default, unless set to no in parameters
            use_gz = ".gz" if self.params["gzip"] else ""
            outputfilename = \
                self.workdir.make_filename("%d-%d%s" % (q0, q1, use_gz))
            self.check_files_exist([outputfilename], "output",
                                   shouldexist=False)
            p = self.programs[0][0](q0=q0, q1=q1,
                                    out=outputfilename,
                                    stats_stderr=True,
                                    skip_check_binary_exists=True,
                                    **fb,
                                    **self.merged_args[0])
            # Note that submit_command may call wait() !
            self.submit_command(p, "%d-%d" % (q0, q1), commit=False)
            self.state.update({"qnext": q1}, commit=True)

        if self.reached_qmax():
            self.logger.info("Reached maximum q-range value %d"
                             " -- not scheduling any further WUs.",
                             self.params["qmax"])
            waitloop = 0
            nextwaitmsg = 1
            while self.get_number_outstanding_wus() > 0:
                self.logger.info("check outstanding wus")
                if self.enough_work_received():
                    break
                waitloop = waitloop + 1
                if waitloop >= nextwaitmsg:
                    self.logger.info("There are still %d outstanding WUs."
                                     "  Server waiting"
                                     " (relations: %d/%d --"
                                     " we have been waiting for %d seconds)",
                                     self.get_number_outstanding_wus(),
                                     self.get_nrels(),
                                     self.state["rels_wanted"],
                                     waitloop)
                    if nextwaitmsg < 30:
                        nextwaitmsg = nextwaitmsg + 1
                    elif nextwaitmsg < 300:
                        nextwaitmsg = nextwaitmsg + 10
                    elif nextwaitmsg < 600:
                        nextwaitmsg = nextwaitmsg + 60
                    elif nextwaitmsg < 3600:
                        nextwaitmsg = nextwaitmsg + 300
                    else:
                        nextwaitmsg = nextwaitmsg + 3600
                self.wait()
            if self.get_number_outstanding_wus() == 0:
                self.logger.info("No outstanding WUs. Sieving stops")
        if self.enough_work_received():
            self.logger.info("Reached target of %d relations, now have %d",
                             self.state["rels_wanted"], self.get_nrels())
        else:
            self.logger.critical("Could not reach target of %d relations,"
                                 " we only have %d -- sieving has failed",
                                 self.state["rels_wanted"], self.get_nrels())
            return False
        self.logger.debug("Exit SievingTask.run(" + self.name + ")")
        return True

    def get_achievement(self):
        return self.state["rels_found"] / self.state["rels_wanted"]

    def updateObserver(self, message):
        identifier = self.filter_notification(message)
        if not identifier:
            # This notification was not for me
            return False
        if self.handle_error_result(message):
            return True
        output_files = message.get_output_files()
        if len(output_files) != 1:
            self.logger.warning("Received output with %d files: %s" % (
                len(output_files), ", ".join(output_files)))
            return False
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

    re_rel = re.compile(bytes(r"(-?\d*),(\d*):(.*)", "ascii"))

    def verify_relation(self, line, poly):
        """
        Check that the primes listed for a relation divide the value of
        the polynomials
        """
        match = self.re_rel.match(line)
        if match:
            a, b, rest = match.groups()
            a, b = int(a), int(b)
            sides = rest.split(b":")
            assert len(sides) == 2
            for side, primes_as_str in enumerate(sides):
                if not primes_as_str:  # empty string
                    continue
                value = poly.get_polynomial(side).eval_h(a, b)
                primes = [int(s, 16) for s in primes_as_str.split(b",")]
                for prime in primes:
                    if value % prime != 0:
                        self.logger.error("Relation %d,%d invalid:"
                                          " %d does not divide %d",
                                          a, b, prime, value)
                        return False
            return True
        return None

    def parse_rel_count(self, filename):
        (name, ext) = os.path.splitext(filename)
        try:
            if ext == ".gz":
                f = gzip.open(filename, "rb")
            else:
                f = open(filename, "rb")
        except (OSError, IOError) as e:
            if e.errno == errno.ENOENT:
                self.logger.error("File '%s' does not exist", filename)
                return None
            else:
                raise
        relations_to_check = 10
        poly = self.send_request(Request.GET_POLYNOMIAL)
        try:
            for line in f:
                if relations_to_check > 0:
                    result = self.verify_relation(line, poly)
                    if result is True:
                        relations_to_check -= 1
                    elif result is False:
                        f.close()
                        return None
                    else:  # Did not match: try again
                        pass
                match = re.match(br"# Total (\d+) reports ", line)
                if match:
                    rels = int(match.group(1))
                    f.close()
                    return rels
        except (IOError, TypeError, structerror, EOFError) as e:
            filtered = [
                    (IOError, r"Not a gzipped file"),
                    (TypeError, r"^ord\(\) expected a character"),
                    (structerror, "^unpack requires a bytes object"),
                    (EOFError, r"^Compressed file ended before")
                    ]
            for t, r in filtered:
                if isinstance(e, t) and re.search(r, str(e)):
                    self.logger.error(f"Error reading '{filename}'"
                                      f"(corrupted?): {e}")
                    return None
            else:
                raise
        f.close()
        self.logger.error(f"Number of relations not found in file {filename}")
        return None

    def import_one_file(self, filename):
        if filename in self.get_output_filenames():
            self.logger.info("Re-scanning file %s", filename)
            nrels = self.get_nrels() - self.get_nrels(filename)
            self.state.update({"rels_found": nrels}, commit=False)
            self.forget_output_filenames([filename], commit=True)
        filename_with_stats_extension = filename + ".stats.txt"
        if os.path.isfile(filename_with_stats_extension):
            self.add_file(filename, filename_with_stats_extension)
        else:
            self.add_file(filename)

    def get_statistics_as_strings(self):
        strings = ["Total number of relations: %d" % self.get_nrels()]
        s1, errors = super().get_statistics_as_strings()
        return strings + s1, errors

    def get_nrels(self, filename=None):
        """
        Return the number of relations found, either the total so far or
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
        """
        Return number of seconds of cpu time spent by las
        """
        s = "stats_total_cpu_time" if is_cpu else "stats_total_time"
        return float(self.statistics.stats.get(s, [0.])[0])


class QuadraticSievingTask(SievingTask):
    @property
    def title(self):
        return "Quadratic Sieving"

    @property
    def programs(self):
        ((_, override, input), ) = super().programs
        return ((cadoprograms.Siqs, override, input),)

    @property
    def paramnames(self):
        return self.join_params(super().paramnames,
                                {"qfac_nfac": int,
                                 "qfac_min": int,
                                 "qfac_max": [int]})


class Duplicates1Task(Task, FilesCreator, HasStatistics):
    """
    Removes duplicate relations
    """

    @property
    def name(self):
        return "duplicates1"

    @property
    def title(self):
        return "Filtering - Duplicate Removal, splitting pass"

    @property
    def programs(self):
        return ((cadoprograms.Duplicates1,
                 ("filelist", "prefix", "out", "nslices_log", "large_ab"),
                 {}),)

    @property
    def paramnames(self):
        return self.join_params(super().paramnames,
                                {"nslices_log": 1, "algo": str})

    @property
    def stat_conversions(self):
        return (
            (
                "stats_dup1_time",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r"# Done: Read \d+ relations in ({fp})s"),
                False
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
        # Enforce the fact that our children *MUST* use the same
        # nslices_log value as the one we have.
        self.progparams[0]["nslices_log"] = self.params["nslices_log"]
        self.progparams[0]["large_ab"] = self.params["algo"] == Algorithm.QS
        tablename = self.make_tablename("infiles")
        self.already_split_input = \
            self.make_db_dict(tablename,
                              connection=self.db_connection)
        self.slice_relcounts = self.make_db_dict(self.make_tablename("counts"),
                                                 connection=self.db_connection)
        # Default slice counts to 0, in single DB commit
        self.slice_relcounts.setdefault(
            None, {str(i): 0 for i in range(0, self.nr_slices)})

    def run(self):
        super().run()

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
                    if f not in self.already_split_input]
        self.logger.debug("new files to split are: %s", newfiles)

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
                outputdir = self.workdir.make_dirname(subdir="dup1")
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
                (stdoutpath, stderrpath) = self.make_std_paths(
                    cadoprograms.Duplicates1.name)

                if len(newfiles) <= 10:
                    p = cadoprograms.Duplicates1(*newfiles,
                                                 prefix=prefix,
                                                 out=outputdir,
                                                 stdout=str(stdoutpath),
                                                 stderr=str(stderrpath),
                                                 **self.progparams[0])
                else:
                    filelistname = self.make_filelist(newfiles, prefix="dup1")
                    p = cadoprograms.Duplicates1(filelist=filelistname,
                                                 prefix=prefix,
                                                 out=outputdir,
                                                 stdout=str(stdoutpath),
                                                 stderr=str(stderrpath),
                                                 **self.progparams[0])
                message = self.submit_command(p, None, log_errors=True)
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
                self.parse_stats(stdoutpath, commit=False)
                # Add relation count from the newly processed files to the
                # relations-per-slice dict
                update1 = {str(idx):
                           self.slice_relcounts[str(idx)] + current_counts[idx]
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
        """
        Takes lines of text and looks for slice counts as printed by dup1
        """
        counts = [None] * self.nr_slices
        for line in stderr.splitlines():
            match = re.match(r'# slice (\d+) received (\d+) relations', line)
            if match:
                (slicenr, nrrels) = map(int, match.groups())
                if counts[slicenr] is not None:
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
    """
    Removes duplicate relations
    """

    @property
    def name(self):
        return "duplicates2"

    @property
    def title(self):
        return "Filtering - Duplicate Removal, removal pass"

    @property
    def programs(self):
        input = {"poly": Request.GET_POLYNOMIAL_FILENAME,
                 "renumber": Request.GET_RENUMBER_FILENAME}
        return ((cadoprograms.Duplicates2,
                ("dlp", "rel_count", "filelist", "large_ab"),
                input),)

    @property
    def paramnames(self):
        return self.join_params(super().paramnames,
                                {"computation": str,
                                 "algo": str,
                                 "nslices_log": 1})

    @property
    def stat_conversions(self):
        return (
            (
                "stats_dup2_time",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r"# Done: Read \d+ relations in ({fp})s"),
                True    # allow several !
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
        self.progparams[0]["large_ab"] = self.params["algo"] == Algorithm.QS
        self.progparams[0]["dlp"] = \
            self.params["computation"] in (Computation.DLP, Computation.CL)
        self.nr_slices = 2**self.params["nslices_log"]
        tablename = self.make_tablename("infiles")
        self.already_done_input = \
            self.make_db_dict(tablename,
                              connection=self.db_connection)
        self.slice_relcounts =  \
            self.make_db_dict(self.make_tablename("counts"),
                              connection=self.db_connection)
        self.slice_relcounts.setdefault(
            None, {str(i): 0 for i in range(0, self.nr_slices)})

    def run(self):
        super().run()

        input_nrel = 0
        for i in range(0, self.nr_slices):
            files = self.send_request(Request.GET_DUP1_FILENAMES, i.__eq__)
            rel_count = self.send_request(Request.GET_DUP1_RELCOUNT, i)
            input_nrel += rel_count
            newfiles = [f for f in files if f not in self.already_done_input]
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
            del self.slice_relcounts[str(i)]
            name = "%s.slice%d" % (cadoprograms.Duplicates2.name, i)
            (stdoutpath, stderrpath) = self.make_std_paths(
                name, do_increment=(i == 0))

            if len(files) <= 10:
                p = cadoprograms.Duplicates2(*files,
                                             rel_count=rel_count,
                                             stdout=str(stdoutpath),
                                             stderr=str(stderrpath),
                                             **self.merged_args[0])
            else:
                filelistname = self.make_filelist(files, prefix="dup1")
                p = cadoprograms.Duplicates2(rel_count=rel_count,
                                             filelist=filelistname,
                                             stdout=str(stdoutpath),
                                             stderr=str(stderrpath),
                                             **self.merged_args[0])
            message = self.submit_command(p, None, log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            with stderrpath.open("r") as stderrfile:
                nr_rels = self.parse_remaining(stderrfile)
            # Mark input file names and output file names
            for f in files:
                self.already_done_input[f] = True
            outfilenames = {f: i for f in files}
            self.add_output_files(outfilenames, commit=False)
            # XXX How do we add the timings ?
            self.parse_stats(stdoutpath, commit=False)
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
            self.logger.error("More new output relations"
                              " (%d) than input (%d)?",
                              new_out, new_in)
            return
        ratio = new_out / new_in
        self.logger.info("Of %d newly added relations %d were unique "
                         "(ratio %f)", new_in, new_out, ratio)
        self.state.update({"last_input_nrel": input_nrel,
                           "last_output_nrel": output_nrel,
                           "unique_ratio": ratio})

    def request_more_relations(self, target):
        nrels = self.get_nrels()
        if target <= nrels:
            return
        additional_out = target - nrels
        ratio = self.state.get("unique_ratio", 1.)
        ratio = max(0.5, ratio)  # avoid stupidly large values of rels_wanted
        additional_in = int(additional_out / ratio)
        newtarget = self.state["last_input_nrel"] + additional_in
        self.logger.info("Got request for"
                         " %d (%d additional) output relations, "
                         "estimate %d (%d additional) needed in input",
                         target, additional_out, newtarget, additional_in)
        self.send_notification(Notification.WANT_MORE_RELATIONS, newtarget)
        self.send_notification(Notification.WANT_TO_RUN, None)


class PurgeTask(Task):
    """
    Removes singletons and computes excess
    """

    @property
    def name(self):
        return "purge"

    @property
    def title(self):
        return "Filtering - Singleton removal"

    @property
    def programs(self):
        override = ("nrels", "out", "outdel",
                    "nprimes", "filelist", "required_excess")
        return ((cadoprograms.Purge, override,
                 {"_freerel": Request.GET_FREEREL_FILENAME}),)

    @property
    def paramnames(self):
        return self.join_params(super().paramnames,
                                {"computation": str,
                                 "galois": "none",
                                 "gzip": True,
                                 "add_ratio": 0.01,
                                 "required_excess": 0.0,
                                 "nmatrices": 0})

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.state.setdefault("input_nrels", 0)
        # We use a computed keep value for DLP and CL
        self.keep = self.progparams[0].pop("keep", None)

    def run(self):
        super().run()

        # Enforce the fact that our children *MUST* use the same
        # required_excess value as the one we have.
        self.progparams[0]["required_excess"] = self.params["required_excess"]

        after_filter_galois = \
            self.params["galois"] in FilterGaloisTask.known_parameters
        if not after_filter_galois:
            nfree = self.send_request(Request.GET_FREEREL_RELCOUNT)
            nunique = self.send_request(Request.GET_UNIQUE_RELCOUNT)
            if not nunique:
                self.logger.critical("No unique relation count received")
                return False
            input_nrels = nfree + nunique
        else:
            # Freerels and Galois are not yet fully compatible.
            input_nrels = self.send_request(Request.GET_GAL_UNIQUE_RELCOUNT)
            if not input_nrels:
                self.logger.critical(
                    "No Galois unique relation count received")
                return False

        nprimes = self.send_request(Request.GET_RENUMBER_PRIMECOUNT)
        # If the user didn't give col_min_index, let's compute it.
        col_min_index = int(self.progparams[0].get("col_min_index", -1))
        if col_min_index == -1:
            # Note: on RSA-120, reducing col_min_index from nprimes/20 to
            # nprimes/40 decreases the matrix size after purge by 3.1%,
            # and the final matrix by 1.4%, while not increasing the memory
            # usage of purge, and increasing the cpu time of purge by only 13%.
            col_min_index = int(nprimes / 40.0)
            # For small cases, we want to avoid degenerated cases, so let's
            # keep most of the ideals: memory is not an issue in that case.
            if (col_min_index < 10000):
                col_min_index = min(500, nprimes-1)
            self.progparams[0].setdefault("col_min_index", col_min_index)

        if "purgedfile" in self.state and not self.have_new_input_files() and \
                input_nrels == self.state["input_nrels"]:
            self.logger.info("Already have a purged file, and no new input "
                             "relations available. Nothing to do")
            self.send_notification(Notification.HAVE_ENOUGH_RELATIONS, None)
            return True

        self.state.pop("purgedfile", None)
        self.state.pop("input_nrels", None)

        if not after_filter_galois:
            self.logger.info("Reading %d unique and %d free relations,"
                             " total %d"
                             % (nunique, nfree, input_nrels))
        else:
            self.logger.info("Reading %d Galois unique" % input_nrels)

        use_gz = ".gz" if self.params["gzip"] else ""
        purgedfile = self.workdir.make_filename("purged" + use_gz)
        if self.params["computation"] == Computation.DLP:
            relsdelfile = self.workdir.make_filename("relsdel" + use_gz)
            nmaps = self.send_request(Request.GET_NMAPS)
            keep = sum(nmaps)
        elif self.params["computation"] == Computation.CL:
            relsdelfile = None
            keep = self.params["nmatrices"]-1
            if self.keep is not None:
                keep = max(self.keep, keep)
                if keep > self.keep:
                    self.logger.warn("Increasing keep from %d to %d because "
                                     "tasks.nmatrices is %d",
                                     self.keep,
                                     keep,
                                     self.params['nmatrices'])
        else:
            relsdelfile = None
            keep = self.keep
        freerel_filename = self.merged_args[0].pop("_freerel", None)
        # Remark: "Galois unique" and "unique" are in the same files
        # because filter_galois works in place. Same request.
        unique_filenames = self.send_request(Request.GET_UNIQUE_FILENAMES)
        if not after_filter_galois and freerel_filename is not None:
            files = unique_filenames + [str(freerel_filename)]
        else:
            files = unique_filenames
        (stdoutpath, stderrpath) = self.make_std_paths(
            cadoprograms.Purge.name)

        if len(files) <= 10:
            p = cadoprograms.Purge(*files,
                                   nrels=input_nrels, out=purgedfile,
                                   outdel=relsdelfile, keep=keep,
                                   nprimes=nprimes,
                                   stdout=str(stdoutpath),
                                   stderr=str(stderrpath),
                                   **self.progparams[0])
        else:
            filelistname = self.make_filelist(files, prefix=self.name)
            p = cadoprograms.Purge(nrels=input_nrels,
                                   out=purgedfile,
                                   outdel=relsdelfile, keep=keep,
                                   nprimes=nprimes,
                                   filelist=filelistname,
                                   stdout=str(stdoutpath),
                                   stderr=str(stderrpath),
                                   **self.progparams[0])
        message = self.submit_command(p, None)
        stdout = message.read_stdout(0).decode()
        # stderr = message.read_stderr(0).decode()

        if self.parse_output(stdout, input_nrels):
            stats = self.parse_stdout(stdout)
            self.logger.info("After purge, %d relations with %d primes remain "
                             "with weight %s and excess %s", *stats)
            output_version = self.state.get("output_version", 0) + 1
            update = {"purgedfile": purgedfile.get_wdir_relative(),
                      "input_nrels": input_nrels,
                      "output_version": output_version}
            if self.params["computation"] == Computation.DLP:
                update["relsdelfile"] = relsdelfile.get_wdir_relative()
            self.state.update(update)
            self.logger.info("Have enough relations")
            self.send_notification(Notification.HAVE_ENOUGH_RELATIONS, None)
        else:
            stats = self.parse_stdout(stdout)
            self.logger.info("After purge, %d relations with %d primes remain "
                             "with excess %d", stats[0], stats[1], stats[3])
            excess = stats[3]
            self.logger.info("Not enough relations")
            if not after_filter_galois:
                self.request_more_relations(nunique, excess)
            else:
                self.request_more_relations(input_nrels, excess)
        self.logger.debug("Exit PurgeTask.run(" + self.name + ")")
        return True

    def request_more_relations(self, nunique, excess):
        r"""
        Given 'nunique' relations and an excess of 'excess',
        estimate how many new (unique) relations we need.
        """

        additional = nunique * self.params["add_ratio"]
        # if the excess is negative, we need at least -excess new relations
        if excess < 0:
            additional = max(additional, -excess)
        required_excess = self.params["required_excess"]
        if required_excess > 0.0:
            # we might need more relations due to required_excess
            nprimes = nunique - excess  # correct whatever the sign of excess
            # we have nprimes ideals, and we want an excess of at least
            # required_excess * nprimes
            target_excess = int(required_excess * nprimes)
            additional = max(additional, target_excess - excess)
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

    def parse_output(self, stdout, input_nrels):
        # If stdout ends with
        # (excess / ncols) = ... < .... See -required_excess argument.'
        # then we need more relations from filtering and return False
        input_nprimes = None
        have_enough = True
        # not_enough1 = re.compile(r"excess < (\d+.\d+) \* #primes")
        not_enough1 = re.compile(r"\(excess / ncols\) = \d+.?\d* < \d+.?\d*. "
                                 r"See -required_excess argument.")
        not_enough2 = re.compile(r"number of rows < number of columns \+ keep")
        nrels_nprimes = re.compile(r"\s*nrows=(\d+), ncols=(\d+); "
                                   r"excess=(-?\d+)")
        for line in stdout.splitlines():
            match = not_enough1.match(line)
            if match:
                have_enough = False
                break
            if not_enough2.match(line):
                have_enough = False
                break
            match = nrels_nprimes.match(line)
            if match is not None:
                (nrels, nprimes, excess) = map(int, match.groups())
                assert nrels - nprimes == excess
                # The first occurrence of the message counts input relations
                if input_nprimes is None:
                    assert input_nrels == nrels
                    input_nprimes = nprimes

        # At this point we should have:
        # input_nrels, input_nprimes: rels and primes among input
        # nrels, nprimes, excess: rels and primes when purging stopped
        if input_nprimes is not None:
            self.update_excess_per_input(input_nrels, nrels, nprimes)
        return have_enough

    def update_excess_per_input(self, input_nrels, nrels, nprimes):
        if input_nrels == 0:
            return  # Nothing sensible that we can do
        last_input_nrels = self.state.get("last_input_nrels", 0)
        last_nrels = self.state.get("last_output_nrels", 0)
        last_nprimes = self.state.get("last_output_nprimes", 0)
        if last_input_nrels >= input_nrels:
            self.logger.warning("Previously stored input nrels (%d) is no "
                                "smaller than value from new run (%d)",
                                last_input_nrels, input_nrels)
            return
        if nrels <= last_nrels:
            self.logger.warning("Previously stored nrels (%d) is no "
                                "smaller than value from new run (%d)",
                                last_nrels, nrels)
            return
        if nprimes <= last_nprimes:
            self.logger.warning("Previously stored nprimes (%d) is no "
                                "smaller than value from new run (%d)",
                                last_nprimes, nprimes)
            return
        self.logger.info("Previous run had %d input relations and ended with "
                         "%d relations and %d primes, new run had %d input "
                         "relations and ended with %d relations and %d primes",
                         last_input_nrels, last_nrels, last_nprimes,
                         input_nrels, nrels, nprimes)
        update = {"last_output_nrels": nrels, "last_output_nprimes": nprimes,
                  "last_input_nrels": input_nrels}
        self.state.update(update)

    def parse_stdout(self, stdout):
        # Program stdout is expected in the form:
        #   Final values:
        #   nrows=23105 ncols=22945 excess=160
        #   weight=382433 weight*nows=8.84e+09
        # but we allow some extra whitespace
        r = {}
        keys = ("nrows", "ncols", "weight", "excess")
        final_values_line = "Final values:"
        had_final_values = False
        for line in stdout.splitlines():
            # Look for values only after we saw "Final values:" line
            if not had_final_values:
                if re.match(final_values_line, line):
                    had_final_values = True
                else:
                    continue
            for key in keys:
                # Match the key at the start of a line, or after a whitespace
                # Note: (?:) is non-capturing group
                match = re.search(r"(?:^|\s)%s\s*=\s*([-]?\d+)" % key, line)
                if match:
                    if key in r:
                        raise Exception("Found multiple values for %s" % key)
                    r[key] = int(match.group(1))
        if not had_final_values:
            raise Exception("%s: output of %s did not contain '%s'" %
                            (self.title,
                             self.programs[0][0].name,
                             final_values_line))
        for key in keys:
            if key not in r:
                raise Exception("%s: output of %s"
                                " did not contain value for %s: %s"
                                % (self.title,
                                   self.programs[0][0].name,
                                   key,
                                   stdout))
        return [r[key] for key in keys]


class FilterGaloisTask(Task):
    """
    Galois Filtering
    """

    known_parameters = ["1/y", "_y", "autom3.1g", "autom3.2g"]

    @property
    def name(self):
        return "filtergalois"

    @property
    def title(self):
        return "Filtering - Galois"

    @property
    def programs(self):
        input = {"poly": Request.GET_POLYNOMIAL_FILENAME,
                 "renumber": Request.GET_RENUMBER_FILENAME}
        return ((cadoprograms.GaloisFilter,
                 ("dl", "nrels", "large_ab"),
                 input),)

    @property
    def paramnames(self):
        return self.join_params(super().paramnames,
                                {"computation": str,
                                 "algo": str,
                                 "galois": "none"})

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.progparams[0]["large_ab"] = self.params["algo"] == Algorithm.QS
        self.progparams[0]["dl"] = \
            self.params["computation"] in (Computation.DLP, Computation.CL)

    def run(self):
        # This task must be run only if galois is recognized by filter_galois
        if self.params["galois"] not in self.known_parameters:
            return True

        super().run()

        files = self.send_request(Request.GET_UNIQUE_FILENAMES)
        nrels = self.send_request(Request.GET_UNIQUE_RELCOUNT)

        (stdoutpath, stderrpath) = self.make_std_paths(
            cadoprograms.GaloisFilter.name)

        # TODO: if there are too many files, pass a filelist (cf PurgeTask)
        p = cadoprograms.GaloisFilter(*files,
                                      nrels=nrels,
                                      stdout=str(stdoutpath),
                                      stderr=str(stderrpath),
                                      **self.merged_args[0])
        message = self.submit_command(p, None, log_errors=True)
        if message.get_exitcode(0) != 0:
            raise Exception("Program failed")
        with stderrpath.open("r") as stderrfile:
            noutrels = self.parse_remaining(stderrfile)
        self.state["noutrels"] = noutrels

        self.logger.debug("Exit FilterGaloisTask.run(" + self.name + ")")
        return True

    def get_nrels(self):
        return self.state["noutrels"]

    def parse_remaining(self, text):
        # "Number of output relations: 303571"
        for line in text:
            match = re.match(r'Number of output relations:\s*(\d+)', line)
            if match:
                remaining = int(match.group(1))
                return remaining
        raise Exception("Received no value for output relation count")

    def request_more_relations(self, target):
        self.send_notification(Notification.WANT_TO_RUN, None)


class MergeDLPTask(Task):
    """
    Merges relations
    """

    @property
    def name(self):
        return "mergedlp"

    @property
    def title(self):
        return "Filtering - Merging"

    @property
    def programs(self):
        input = {"purged": Request.GET_PURGED_FILENAME}
        return ((cadoprograms.MergeDLP, ("out"), input),
                (cadoprograms.ReplayDLP, ("ideals",
                                          "history",
                                          "index",
                                          "out"), input))

    @property
    def paramnames(self):
        return self.join_params(super().paramnames, {"gzip": True})

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.progparams[0]["skip"] = 0

    def run(self):
        super().run()

        if "mergedfile" not in self.state or self.have_new_input_files():
            if "idealfile" in self.state:
                del self.state["idealfile"]
            if "indexfile" in self.state:
                del self.state["indexfile"]
            if "mergedfile" in self.state:
                del self.state["mergedfile"]
            if "densefile" in self.state:
                del self.state["densefile"]

            # nmaps = self.send_request(Request.GET_NMAPS)
            # keep = sum(nmaps)
            # We use .gzip by default, unless set to no in parameters
            use_gz = ".gz" if self.params["gzip"] else ""
            historyfile = self.workdir.make_filename("history" + use_gz)
            (stdoutpath, stderrpath) = self.make_std_paths(
                cadoprograms.MergeDLP.name)
            p = cadoprograms.MergeDLP(out=historyfile,
                                      stdout=str(stdoutpath),
                                      stderr=str(stderrpath),
                                      **self.merged_args[0])
            message = self.submit_command(p, None, log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            stdout = message.read_stdout(0).decode()
            matsize = 0
            matweight = 0
            for line in stdout.splitlines():
                match = re.match(r'Final matrix has'
                                 r' N=(\d+) nc=\d+ \(\d+\) W=(\d+)', line)
                if match:
                    matsize = int(match.group(1))
                    matweight = int(match.group(2))
            if (matsize == 0) or (matweight == 0):
                raise Exception("Could not read matrix size and weight")
            self.logger.info("Merged matrix has"
                             " %d rows and total weight %d"
                             " (%.1f entries per row on average)"
                             % (matsize, matweight,
                                float(matweight)/float(matsize)))

            indexfile = self.workdir.make_filename("index" + use_gz)
            mergedfile = self.workdir.make_filename("sparse.bin")
            (stdoutpath, stderrpath) = self.make_std_paths(
                cadoprograms.Replay.name)
            idealfile = self.workdir.make_filename("ideal")
            p = cadoprograms.ReplayDLP(ideals=idealfile,
                                       history=historyfile, index=indexfile,
                                       out=mergedfile, stdout=str(stdoutpath),
                                       stderr=str(stderrpath),
                                       **self.merged_args[1])
            message = self.submit_command(p, None, log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")

            if not idealfile.isfile():
                raise Exception("Output file %s does not exist" % idealfile)
            if not indexfile.isfile():
                raise Exception("Output file %s does not exist" % indexfile)
            if not mergedfile.isfile():
                raise Exception("Output file %s does not exist" % mergedfile)
            output_version = self.state.get("output_version", 0) + 1
            self.remember_input_versions(commit=False)
            update = {"indexfile": indexfile.get_wdir_relative(),
                      "mergedfile": mergedfile.get_wdir_relative(),
                      "idealfile": idealfile.get_wdir_relative(),
                      "output_version": output_version}
            densefilename = self.workdir.make_filename("dense.bin")
            if densefilename.isfile():
                update["densefile"] = densefilename.get_wdir_relative()
            self.state.update(update)

        self.logger.debug("Exit MergeDLPTask.run(" + self.name + ")")
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
    """
    Merges relations
    """
    @property
    def name(self):
        return "mergetask"

    @property
    def title(self):
        return "Filtering - Merging"

    @property
    def programs(self):
        input = {"purged": Request.GET_PURGED_FILENAME}
        return ((cadoprograms.Merge, ("out",), input),
                (cadoprograms.Replay, ("out", "history", "index"), input))

    @property
    def paramnames(self):
        return self.join_params(super().paramnames,
                                {"skip": None, "gzip": True})

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        skip = int(self.progparams[0].get("skip", 32))
        self.progparams[0].setdefault("skip", skip)

    def run(self):
        super().run()

        if "mergedfile" in self.state:
            fn = self.workdir.path_in_workdir(self.state["mergedfile"])
            if not fn.isfile():
                self.logger.warning("Output file %s disappeared,"
                                    " generating it again", fn)
                del self.state["mergedfile"]
        if "mergedfile" not in self.state or self.have_new_input_files():
            if "indexfile" in self.state:
                del self.state["indexfile"]
            if "densefile" in self.state:
                del self.state["densefile"]

            # We use .gzip by default, unless set to no in parameters
            use_gz = ".gz" if self.params["gzip"] else ""
            historyfile = self.workdir.make_filename("history" + use_gz)
            (stdoutpath, stderrpath) = self.make_std_paths(
                cadoprograms.Merge.name)
            p = cadoprograms.Merge(out=historyfile,
                                   stdout=str(stdoutpath),
                                   stderr=str(stderrpath),
                                   **self.merged_args[0])
            message = self.submit_command(p, None, log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            stdout = message.read_stdout(0).decode()
            matsize = 0
            matweight = 0
            for line in stdout.splitlines():
                match = re.match(r'Final matrix has'
                                 r' N=(\d+) nc=\d+ \(\d+\) W=(\d+)', line)
                if match:
                    matsize = int(match.group(1))
                    matweight = int(match.group(2))
            if (matsize == 0) or (matweight == 0):
                raise Exception("Could not read matrix size and weight")
            self.logger.info("Merged matrix has"
                             " %d rows and total weight %d"
                             " (%.1f entries per row on average)"
                             % (matsize,
                                matweight,
                                float(matweight)/float(matsize)))

            indexfile = self.workdir.make_filename("index" + use_gz)
            mergedfile = self.workdir.make_filename("sparse.bin")
            (stdoutpath, stderrpath) = self.make_std_paths(
                cadoprograms.Replay.name)
            p = cadoprograms.Replay(history=historyfile, index=indexfile,
                                    out=mergedfile, stdout=str(stdoutpath),
                                    stderr=str(stderrpath),
                                    **self.merged_args[1])
            message = self.submit_command(p, None, log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")

            if not indexfile.isfile():
                raise Exception("Output file %s does not exist" % indexfile)
            if not mergedfile.isfile():
                raise Exception("Output file %s does not exist" % mergedfile)
            self.remember_input_versions(commit=False)
            output_version = self.state.get("output_version", 0) + 1
            update = {"indexfile": indexfile.get_wdir_relative(),
                      "mergedfile": mergedfile.get_wdir_relative(),
                      "output_version": output_version}
            densefilename = self.workdir.make_filename("dense.bin")
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


class NumberTheoryTask(Task):
    """
    Number theory tasks for dlp
    """
    @property
    def name(self):
        return "numbertheory"

    @property
    def title(self):
        return "Number Theory for DLP"

    @property
    def programs(self):
        return ((cadoprograms.NumberTheory,
                (),
                 {"poly": Request.GET_POLYNOMIAL_FILENAME}),)

    @property
    def paramnames(self):
        return self.join_params(super().paramnames, {"nsm0": -1, "nsm1": -1})

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)

    def run(self):
        super().run()

        # Check if the numbertheory program was run already.
        if "nmaps0" in self.state:
            self.logger.info("NumberTheory task has already run,"
                             " reusing the result.")
            return True

        # Create output files and start the computation
        (stdoutpath, stderrpath) = self.make_std_paths(
            cadoprograms.NumberTheory.name)
        p = cadoprograms.NumberTheory(stdout=str(stdoutpath),
                                      stderr=str(stderrpath),
                                      **self.merged_args[0])
        message = self.submit_command(p, None, log_errors=True)
        if message.get_exitcode(0) != 0:
            raise Exception("Program failed")

        stdout = message.read_stdout(0).decode()
        update = {}
        for line in stdout.splitlines():
            match = re.match(r'# nmaps0 (\d+)', line)
            if match:
                update["nmaps0"] = int(match.group(1))
            match = re.match(r'# nmaps1 (\d+)', line)
            if match:
                update["nmaps1"] = int(match.group(1))
        # Allow user-given parameter to override what we compute:
        if self.params["nsm0"] != -1:
            update["nmaps0"] = self.params["nsm0"]
        if self.params["nsm1"] != -1:
            update["nmaps1"] = self.params["nsm1"]

        if "nmaps0" not in update:
            raise Exception("Stdout does not give nmaps0")
        if "nmaps1" not in update:
            raise Exception("Stdout does not give nmaps1")
        # Update the state entries atomically
        self.state.update(update)

        self.logger.debug("Exit NumberTheoryTask.run(" + self.name + ")")
        return True

    def get_nmaps(self):
        return (self.state["nmaps0"], self.state["nmaps1"])


class bwc_output_filter(RealTimeOutputFilter):
    def filter(self, data):
        super().filter(data)
        if ("ETA" or "Timings") in data:
            self.logger.info(data.rstrip())


# I've just ditched the statistics bit, cause I don't know to make its
# despair cry a little bit more useful.
class LinAlgDLPTask(Task):
    """
    Runs the linear algebra step for DLP
    """
    @property
    def name(self):
        return "linalgdlp"

    @property
    def title(self):
        return "Linear Algebra"

    @property
    def programs(self):
        override = ("complete", "rhs", "matrix",  "wdir",
                    "nullspace", "m", "n")
        return ((cadoprograms.BWC, override,
                 {"merged": Request.GET_MERGED_FILENAME,
                  "sm": Request.GET_SM_FILENAME}),)

    @property
    def paramnames(self):
        # the default value for m and n is to use the number of SMs for
        # n, and then m=2*n
        return self.join_params(super().paramnames,
                                {"m": 0, "n": 0,
                                 "ell": int,
                                 "force_wipeout": False})

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.state.setdefault("ran_already", False)

    def run(self):
        super().run()

        if self.state["ran_already"] and self.params["force_wipeout"]:
            self.logger.warning("Ran before, but force_wipeout is set. "
                                "Wiping out working directory.")
            self.workdir.make_dirname(subdir="bwc").rmtree()
            self.state["ran_already"] = False
            self.state.pop("virtual_logs", None)

        if "virtual_logs" not in self.state or self.have_new_input_files():
            workdir = self.workdir.make_dirname(subdir="bwc")
            workdir.mkdir(parent=True)
            mergedfile = self.merged_args[0].pop("merged")
            smfile = self.merged_args[0].pop("sm")
            if mergedfile is None:
                self.logger.critical("No merged file received.")
                return False
            (stdoutpath, stderrpath) = self.make_std_paths(
                cadoprograms.BWC.name)
            matrix = mergedfile.realpath()
            wdir = workdir.realpath()
            self.state["ran_already"] = True
            nmaps = self.send_request(Request.GET_NMAPS)
            nsm = sum(nmaps)
            if self.params["n"] == 0:
                self.logger.info("Using %d as default value for n"
                                 " to account for Schirokauer maps"
                                 % nsm)
                n = nsm
            else:
                n = self.params["n"]
                if n < nsm:
                    self.logger.critical("n must be greater than or equal to"
                                         " the number of Schirokauer maps,"
                                         " which is %d (got n=%d)" % (nsm, n))
                    raise Exception("Program failed")

            if self.params["m"] == 0:
                m = 2*n
                self.logger.info("Using 2*n=%d as default value for m" % m)
            else:
                m = self.params["m"]

            if n == 0:
                self.logger.error("Error: homogeneous Linalg"
                                  " is not implemented")
                raise Exception("Program failed")

            p = cadoprograms.BWC(complete=True,
                                 matrix=matrix,
                                 wdir=wdir,
                                 rhs=smfile,
                                 prime=self.params["ell"],
                                 nullspace="right",
                                 stdout=str(stdoutpath),
                                 stderr=str(stderrpath),
                                 m=m,
                                 n=n,
                                 **self.progparams[0])
            message = self.submit_command(p, None, log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            virtual_logs_filename = \
                self.workdir.make_filename("K.sols0-1.0.txt", subdir="bwc")
            if not virtual_logs_filename.isfile():
                raise Exception("Kernel file %s does not exist"
                                % virtual_logs_filename)
            self.remember_input_versions(commit=False)
            output_version = self.state.get("output_version", 0) + 1
            update = {"virtual_logs":
                      virtual_logs_filename.get_wdir_relative(),
                      "output_version":
                      output_version}
            self.state.update(update, commit=True)
        self.logger.debug("Exit LinAlgTask.run(" + self.name + ")")
        return True

    def get_virtual_logs_filename(self):
        return self.get_state_filename("virtual_logs")

    def get_prefix(self):
        return "%s%s%s.%s" % (self.params["workdir"].rstrip(os.sep), os.sep,
                              self.params["name"], "dep")


class LinAlgClTask(Task):
    """ Runs the linear algebra step for Cl """
    @property
    def name(self):
        return "linalgcl"

    @property
    def title(self):
        return "Linear Algebra"

    @property
    def programs(self):
        override = ("determinant", "matrix", "wdir", "nullspace", "m", "n")
        return ((cadoprograms.BWC, override,
                 {"merged": Request.GET_MERGED_FILENAME}),)

    @property
    def paramnames(self):
        # the default value for m and n is n=1, and then m=2*n
        return self.join_params(super().paramnames,
                                {"m": [int], "n": [int],
                                 "force_wipeout": False,
                                 "nmatrices": 1})

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.state.setdefault("ran_already", "")

    def run(self):
        super().run()
        if self.state["ran_already"] and self.params["force_wipeout"]:
            self.logger.warning("Ran before, but force_wipeout is set. "
                                "Wiping out working directories.")
            for dirname in self.state["ran_already"].split(","):
                self.workdir.make_dirname(subdir=dirname).rmtree()
            self.state["ran_already"] = ""
            self.state.pop("classnumber", None)

        if "classnumber" not in self.state:
            # First build the list of matrices
            matrices = []
            mergedfile = self.merged_args[0].pop("merged")
            if mergedfile is None:
                self.logger.critical("No merged file received.")
                return False
            basefilepath = mergedfile.filepath[:-4]
            for i in range(self.params["nmatrices"]):
                if i:
                    mergedfile.filepath = "%s.%d.bin" % (basefilepath, i)
                if not mergedfile.isfile():
                    self.logger.critical(f"Missing merged file {mergedfile}.")
                    return False
                matrices.append(mergedfile.realpath())

            if "n" not in self.params:
                self.logger.info("Using 1 as default value for n")
                n = 1
            else:
                n = self.params["n"]
                if n < 1:
                    self.logger.critical("n must be greater than 0 (got n=%d)"
                                         % n)
                    raise Exception("Program failed")

            if "m" not in self.params:
                m = 2*n
                self.logger.info("Using 2*n=%d as default value for m" % m)
            else:
                m = self.params["m"]

            h = None
            for i, matrix in enumerate(matrices):
                self.logger.info("Starting linalg for matrix #%d" % i)
                det, M = 0, 1
                break_on_next_eq = False
                for j, ell in enumerate(primes_above(2**63)):
                    self.logger.info(f"Performing linalg modulo {ell}")
                    dirname = f"bwc.h{i:02d}.p{j:05d}"
                    workdir = self.workdir.make_dirname(subdir=dirname)
                    workdir.mkdir(parent=True)
                    wdir = workdir.realpath()
                    if not self.state["ran_already"]:
                        self.state["ran_already"] = dirname
                    else:
                        self.state["ran_already"] += "," + dirname

                    (stdoutpath, stderrpath) = self.make_std_paths(
                                                    dirname,
                                                    do_increment=False)
                    p = cadoprograms.BWC(determinant=True,
                                         matrix=matrix,
                                         wdir=wdir,
                                         prime=ell,
                                         nullspace="right",
                                         stdout=str(stdoutpath),
                                         stderr=str(stderrpath),
                                         m=m,
                                         n=n,
                                         **self.progparams[0])
                    message = self.submit_command(p, None)
                    stdout = message.read_stdout(0).decode("utf-8")
                    match = re.search(r"^determinant modulo p: (.*)$", stdout,
                                      re.MULTILINE)
                    if message.get_exitcode(0) != 0:
                        if match and match.group(1) == "inconclusive":
                            self.logger.warn("Skipping prime %d, could not "
                                             "compute the determinant", ell)
                        else:
                            self.log_failed_command_error(message, 0)
                            raise Exception("Program failed")
                    else:
                        if not match:
                            raise Exception("Could not parse output")
                        detell = int(match.group(1)) % ell
                        if 2*detell > ell:
                            detell -= ell
                        old_det = det
                        det, M = CRT(det, M, detell, ell, signed=True)
                        if det == old_det:
                            if break_on_next_eq:
                                break
                            else:
                                break_on_next_eq = True
                        else:
                            break_on_next_eq = False
                else:
                    self.logger.error("No more primes for CRT, probably due "
                                      "to an error: current value is "
                                      "%d modulo %d", det, M)
                    raise Exception("Could not finish CRT computation")

                self.logger.info(f"h multiple #{i} = {det}")
                if h is None:
                    h = abs(det)
                else:
                    h = xgcd(h, det)[0]
                self.logger.info(f"Current multiple of class number: {h}")

            self.state["classnumber"] = h

        h = self.state["classnumber"]
        self.logger.info(f"Candidate class number: {h}")

        self.logger.debug("Exit LinAlgClTask.run(" + self.name + ")")
        return True

    def get_candidate_class_number(self):
        return self.state.get("classnumber", None)


class LinAlgTask(Task, HasStatistics):
    """
    Runs the linear algebra step
    """
    @property
    def name(self):
        return "linalg"

    @property
    def title(self):
        return "Linear Algebra"

    @property
    def programs(self):
        return ((cadoprograms.BWC,
                 ("complete", "matrix",  "wdir", "nullspace"),
                 {"merged": Request.GET_MERGED_FILENAME}),)

    @property
    def paramnames(self):
        return self.join_params(super().paramnames, {"force_wipeout": False})

    @property
    def stat_conversions(self):
        return (
            (
                "prep_wct",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r'Timings for prep: .wct. ({fp})'),
                False
            ),
            (
                "prep_cpu",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r'Timings for prep: .cpu. ({fp})'),
                False
            ),
            (
                "secure_wct",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r'Timings for secure: .wct. ({fp})'),
                False
            ),
            (
                "secure_cpu",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r'Timings for secure: .cpu. ({fp})'),
                False
            ),
            (
                "gather_wct",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r'Timings for gather: .wct. ({fp})'),
                False
            ),
            (
                "gather_cpu",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r'Timings for gather: .cpu. ({fp})'),
                False
            ),
            (
                "krylov_wct",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r"Timings for krylov: .wct. ({fp})"),
                True
            ),
            (
                "krylov_cpu",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r"Timings for krylov: .cpu. ({fp})"),
                True
            ),
            (
                "krylov_iteration_cpu",
                (int, float),
                "0",
                Statistics.add_list,
                re_fp_compile(r"krylov done N=(\d+) ; CPU\d*: ({fp})"),
                True
            ),
            (
                "krylov_iteration_cpu_wait",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r"krylov done N=\d+ ; cpu-wait\d*: ({fp})"),
                True
            ),
            (
                "krylov_iteration_comm",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r"krylov done N=\d+ ; COMM\d*: ({fp})"),
                True
            ),
            (
                "krylov_iteration_comm_wait",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r"krylov done N=\d+ ; comm-wait\d*: ({fp})"),
                True
            ),
            (
                "lingen_wct",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r'Timings for lingen_\w+: .wct. ({fp})'),
                False
            ),
            (
                "lingen_cpu",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r'Timings for lingen_\w+: .cpu. ({fp})'),
                False
            ),
            (
                "mksol_wct",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r"Timings for mksol: .wct. ({fp})"),
                True
            ),
            (
                "mksol_cpu",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r"Timings for mksol: .cpu. ({fp})"),
                True
            ),
            (
                "mksol_iteration_cpu",
                (int, float),
                "0",
                Statistics.add_list,
                re_fp_compile(r"mksol done N=(\d+) ; CPU\d*: ({fp})"),
                True
            ),
            (
                "mksol_iteration_cpu_wait",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r"mksol done N=\d+ ; cpu-wait\d*: ({fp})"),
                True
            ),
            (
                "mksol_iteration_comm",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r"mksol done N=\d+ ; COMM\d*: ({fp})"),
                True
            ),
            (
                "mksol_iteration_comm_wait",
                float,
                "0",
                Statistics.add_list,
                re_fp_compile(r"mksol done N=\d+ ; comm-wait\d*: ({fp})"),
                True
            ),
        )

    @property
    def stat_formats(self):
        return (
            ["Krylov: CPU time {krylov_cpu[0]}",
                ", WCT time {krylov_wct[0]}",
                ", iteration CPU time {krylov_iteration_cpu[1]:g}",
                ", COMM {krylov_iteration_comm[0]}",
                ", cpu-wait {krylov_iteration_cpu_wait[0]}",
                ", comm-wait {krylov_iteration_comm_wait[0]}",
                " ({krylov_iteration_cpu[0]:d} iterations)"
             ],
            ["Lingen CPU time {lingen_cpu[0]}", ", WCT time {lingen_wct[0]}"],
            ["Mksol: CPU time {mksol_cpu[0]}",
                ",  WCT time {mksol_wct[0]}",
                ", iteration CPU time {mksol_iteration_cpu[1]:g}",
                ", COMM {mksol_iteration_comm[0]}",
                ", cpu-wait {mksol_iteration_cpu_wait[0]}",
                ", comm-wait {mksol_iteration_comm_wait[0]}",
                " ({mksol_iteration_cpu[0]:d} iterations)"
             ],
        )

    def get_total_cpu_or_real_time(self, is_cpu):
        """
        Return number of seconds of cpu time spent by bwc
        """
        what = "_cpu" if is_cpu else "_wct"
        s = 0
        for step in ["prep", "secure", "gather", "lingen", "krylov", "mksol"]:
            s += float(self.statistics.stats.get(step + what, [0.])[0])
        return s

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.state.setdefault("ran_already", False)

    def run(self):
        super().run()

        if self.state["ran_already"] and self.params["force_wipeout"]:
            self.logger.warning("Ran before, but force_wipeout is set. "
                                "Wiping out working directory.")
            self.workdir.make_dirname(subdir="bwc").rmtree()
            self.state["ran_already"] = False
            self.state.pop("dependency", None)

        if "dependency" not in self.state or self.have_new_input_files():
            workdir = self.workdir.make_dirname(subdir="bwc")
            workdir.mkdir(parent=True)
            mergedfile = self.merged_args[0].pop("merged")
            if mergedfile is None:
                self.logger.critical("No merged file received.")
                return False
            (stdoutpath, stderrpath) = self.make_std_paths(
                cadoprograms.BWC.name)
            matrix = mergedfile.realpath()
            wdir = workdir.realpath()
            self.state["ran_already"] = True
            self.remember_input_versions(commit=True)
            with bwc_output_filter(self.logger, str(stdoutpath)) as outfilter:
                p = cadoprograms.BWC(complete=True,
                                     matrix=matrix,
                                     wdir=wdir,
                                     nullspace="left",
                                     stdout=outfilter,
                                     stderr=str(stderrpath),
                                     **self.progparams[0])
                message = self.submit_command(p, None, log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            dependencyfilename = self.workdir.make_filename("W", subdir="bwc")
            if not dependencyfilename.isfile():
                raise Exception("Kernel file %s does not exist"
                                % dependencyfilename)
            self.logger.debug("Parsing stats from %s" % stdoutpath)
            self.parse_stats(stdoutpath, commit=False)

            # for some reason, submit_command does not properly catch the
            # time taken by bwc. Is it related to bwc being a perl
            # script, or maybe is it related to the output filter ? At
            # any rate, we're better off resetting the real time to what
            # we've read from the output.
            real = self.get_total_cpu_or_real_time(False)
            self.state.update({"realtime_bwc": real}, commit=True)

            output_version = self.state.get("output_version", 0) + 1
            update = {"dependency": dependencyfilename.get_wdir_relative(),
                      "output_version": output_version}
            self.state.update(update, commit=True)
        self.logger.debug("Exit LinAlgTask.run(" + self.name + ")")
        return True

    def get_dependency_filename(self):
        return self.get_state_filename("dependency")

    def get_prefix(self):
        return "%s%s%s.%s" % (self.params["workdir"].rstrip(os.sep), os.sep,
                              self.params["name"], "dep")


class CharactersTask(Task):
    """
    Computes Quadratic Characters
    """

    @property
    def name(self):
        return "characters"

    @property
    def title(self):
        return "Quadratic Characters"

    @property
    def programs(self):
        input = {"poly": Request.GET_POLYNOMIAL_FILENAME,
                 "wfile": Request.GET_DEPENDENCY_FILENAME,
                 "purged": Request.GET_PURGED_FILENAME,
                 "index": Request.GET_INDEX_FILENAME,
                 "heavyblock": Request.GET_DENSE_FILENAME}
        override = ("out", "large_ab", "only_sign_chars")
        return ((cadoprograms.Characters, override, input),)

    @property
    def paramnames(self):
        return self.join_params(super().paramnames, {"algo": str})

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.progparams[0]["large_ab"] = self.params["algo"] == Algorithm.QS
        self.progparams[0]["only_sign_chars"] = \
            self.params["algo"] == Algorithm.QS

    def run(self):
        super().run()

        if "kernel" not in self.state or self.have_new_input_files():
            kernelfilename = self.workdir.make_filename("kernel")
            (stdoutpath, stderrpath) = self.make_std_paths(
                cadoprograms.Characters.name)
            p = cadoprograms.Characters(out=kernelfilename,
                                        stdout=str(stdoutpath),
                                        stderr=str(stderrpath),
                                        **self.merged_args[0])
            message = self.submit_command(p, None, log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            if not kernelfilename.isfile():
                raise Exception("Output file %s does not exist"
                                % kernelfilename)
            self.remember_input_versions(commit=False)
            update = {"kernel": kernelfilename.get_wdir_relative()}
            self.state.update(update)
        self.logger.debug("Exit CharactersTask.run(" + self.name + ")")
        return True

    def get_kernel_filename(self):
        return self.get_state_filename("kernel")


class SqrtTask(Task):
    """
    Runs the square root
    """

    @property
    def name(self):
        return "sqrt"

    @property
    def title(self):
        return "Square Root"

    @property
    def programs(self):
        input = {"poly": Request.GET_POLYNOMIAL_FILENAME,
                 "purged": Request.GET_PURGED_FILENAME,
                 "index": Request.GET_INDEX_FILENAME,
                 "kernel": Request.GET_KERNEL_FILENAME}
        return ((cadoprograms.Sqrt,
                 ("ab", "prefix", "side0", "side1", "gcd", "dep", "large_ab",
                  "qs"),
                 input), )

    @property
    def paramnames(self):
        return self.join_params(super().paramnames,
                                {"N": int, "gzip": True, "first_dep": [int],
                                 "algo": str})

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.progparams[0]["large_ab"] = self.params["algo"] == Algorithm.QS
        self.factors = self.make_db_dict(self.make_tablename("factors"),
                                         connection=self.db_connection)
        self.add_factor(self.params["N"])
        if "first_dep" in self.params:
            self.state["next_dep"] = self.params["first_dep"]

    def run(self):
        super().run()

        if not self.is_done() or self.have_new_input_files():
            prefix = self.send_request(Request.GET_LINALG_PREFIX)
            if self.params["gzip"]:
                prefix += ".gz"
            (stdoutpath, stderrpath) = self.make_std_paths(
                cadoprograms.Sqrt.name)
            self.logger.info("Creating file of (a,b) values")
            purged = self.merged_args[0].pop("purged", None)
            index = self.merged_args[0].pop("index", None)
            kernel = self.merged_args[0].pop("kernel", None)
            p = cadoprograms.Sqrt(ab=True,
                                  prefix=prefix, purged=purged,
                                  index=index, kernel=kernel,
                                  stdout=str(stdoutpath),
                                  stderr=str(stderrpath),
                                  **self.merged_args[0])
            message = self.submit_command(p, None, log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")

            t = self.progparams[0].get("threads", 1)
            while not self.is_done():
                dep = self.state.get("next_dep", 0)
                # if t == 1:
                #    self.logger.info("Trying dependency %d", dep)
                # else:
                #    self.logger.info("Trying dependencies %d to %d",
                #                     dep, dep+t-1)
                (stdoutpath, stderrpath) = self.make_std_paths(
                    cadoprograms.Sqrt.name)
                args = {"ab": False, "dep": dep, "prefix": prefix, "gcd": True,
                        "stdout": str(stdoutpath), "stderr": str(stderrpath)}
                if self.params["algo"] == Algorithm.NFS:
                    args["side0"] = args["side1"] = True
                else:
                    args["qs"] = True
                p = cadoprograms.Sqrt(**args, **self.merged_args[0])
                message = self.submit_command(p, f"dep{dep}", log_errors=True)
                if message.get_exitcode(0) != 0:
                    raise Exception("Program failed")
                with stdoutpath.open("r") as stdoutfile:
                    stdout = stdoutfile.read()
                lines = stdout.splitlines()
                for line in lines:
                    if line == "Failed":
                        # try next lines (if any) in multi-thread mode
                        continue
                    self.add_factor(int(line))
                self.state.update({"next_dep": dep+t})
            self.remember_input_versions(commit=True)
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
                del self.factors[str(oldfac)]
                self.add_factor(g)
                self.add_factor(oldfac // g)
                break
        else:
            # We get here if the new factor is coprime to all previously
            # known factors
            isprime = SqrtTask.miller_rabin_tests(factor, 10)
            self.factors[str(factor)] = isprime

    def get_factors(self):
        N = self.params["N"]
        F = []
        for p in sorted([int(p) for p in self.factors.keys()]):
            assert N % p == 0
            while N % p == 0:
                F.append(str(p))
                N = N // p
        return F

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
        """
        Return the smallest strong probable prime no smaller than N
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
    """
    Computes Schirokauer Maps
    """

    @property
    def name(self):
        return "sm"

    @property
    def title(self):
        return "Schirokauer Maps"

    @property
    def programs(self):
        override = ("nsm", "out")
        input = {"poly": Request.GET_POLYNOMIAL_FILENAME,
                 "purged": Request.GET_PURGED_FILENAME,
                 "index": Request.GET_INDEX_FILENAME}
        return ((cadoprograms.SM, override, input),)

    @property
    def paramnames(self):
        return super().paramnames

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)

    def run(self):
        super().run()

        if "sm" not in self.state or self.have_new_input_files():
            nmaps = self.send_request(Request.GET_NMAPS)
            if sum(nmaps) == 0:
                self.logger.info("Number of SM is 0: skipping this part.")
                return True
            smfilename = self.workdir.make_filename("sm")

            (stdoutpath, stderrpath) = self.make_std_paths(
                cadoprograms.SM.name)
            p = cadoprograms.SM(nsms=",".join(str(nmap) for nmap in nmaps),
                                out=smfilename,
                                stdout=str(stdoutpath),
                                stderr=str(stderrpath),
                                **self.merged_args[0])
            message = self.submit_command(p, None, log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            if not smfilename.isfile():
                raise Exception("Output file %s does not exist" % smfilename)
            self.state["sm"] = smfilename.get_wdir_relative()
            self.remember_input_versions()
        self.logger.debug("Exit SMTask.run(" + self.name + ")")
        return True

    def get_sm_filename(self):
        return self.get_state_filename("sm")


class ReconstructLogTask(Task):
    """
    Logarithms Reconstruction Task
    """

    @property
    def name(self):
        return "reconstructlog"

    @property
    def title(self):
        return "Logarithms Reconstruction"

    @property
    def programs(self):
        input = {
                "ker": Request.GET_KERNEL_FILENAME,
                "poly": Request.GET_POLYNOMIAL_FILENAME,
                "renumber": Request.GET_RENUMBER_FILENAME,
                "purged": Request.GET_PURGED_FILENAME,
                "ideals": Request.GET_IDEAL_FILENAME,
                "relsdel": Request.GET_RELSDEL_FILENAME,
                }
        override = ("dlog", "nrels")
        return ((cadoprograms.ReconstructLog, override, input),)

    @property
    def paramnames(self):
        return self.join_params(super().paramnames,
                                {"checkdlp": True, "jlpoly": False})

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)

    def run(self):
        super().run()

        if "dlog" not in self.state or self.have_new_input_files():
            dlogfilename = self.workdir.make_filename("dlog")
            nmaps = self.send_request(Request.GET_NMAPS)
            nsms = ",".join(str(nmap) for nmap in nmaps)

            nfree = self.send_request(Request.GET_FREEREL_RELCOUNT)
            nunique = self.send_request(Request.GET_UNIQUE_RELCOUNT)

            (stdoutpath, stderrpath) = self.make_std_paths(
                cadoprograms.ReconstructLog.name)
            p = cadoprograms.ReconstructLog(dlog=dlogfilename,
                                            nsms=nsms,
                                            nrels=nfree+nunique,
                                            stdout=str(stdoutpath),
                                            stderr=str(stderrpath),
                                            **self.merged_args[0])
            message = self.submit_command(p, None, log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            if not dlogfilename.isfile():
                raise Exception("Output file %s does not exist" % dlogfilename)
            self.state["dlog"] = dlogfilename.get_wdir_relative()
            self.remember_input_versions()
        self.logger.debug("Exit ReconstructLogTask.run(" + self.name + ")")
        return True

    def get_dlog_filename(self):
        return self.get_state_filename("dlog")


# TODO: This is a bit ugly. We're leaning on the functionality that
# descent.py infers the complete set of file names from the prefix (or
# from the database, if it so wishes). However, the cadofactor way would
# be to pass each and every needed file name as provided by the mediator.
class DescentTask(Task):
    """
    Individual logarithm Task
    """

    @property
    def name(self):
        return "descent"

    @property
    def title(self):
        return "Individual logarithm"

    @property
    def programs(self):
        input = {
                "prefix": Request.GET_WORKDIR_JOBNAME,
                "datadir": Request.GET_WORKDIR_PATH,
                }
        override = ("cadobindir", "target")
        return ((cadoprograms.Descent, override, input),)

    @property
    def paramnames(self):
        return self.join_params(super().paramnames,
                                {"target": [str],
                                 "gfpext": [int, 1],
                                 "execpath": str})

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.logtargets = []

    def run(self):
        super().run()

        if self.params["target"] is None:
            self.logger.info("Skipping descent,"
                             " as no target= argument was passed")
            return

        for ts in self.params["target"].split(","):
            target = int(ts)

            self.logger.info("Now doing descent for target=%s" % target)

            (stdoutpath, stderrpath) = self.make_std_paths(
                cadoprograms.Descent.name)
            p = cadoprograms.Descent(cadobindir=self.params["execpath"],
                                     stdout=str(stdoutpath),
                                     stderr=str(stderrpath),
                                     target=target,
                                     **self.merged_args[0])
            message = self.submit_command(p, None, log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")

            stdout = message.read_stdout(0).decode()
            for line in stdout.splitlines():
                match = re.match(r'log\(target\)=(\d+)', line)
                if match:
                    logtarget = int(match.group(1))
                    self.logger.info("Descent yields log(%s)=%s"
                                     % (target, logtarget))
                    self.logtargets.append(logtarget)

                    checker = self.send_request(Request.GET_LOGQUERY_CHECKER)
                    checker.check_new_log(target, logtarget)
                    break
        return True

    def get_logtargets(self):
        return self.logtargets


# This simple task is here to collect in one single file the logarithms
# that we've queried, and check their consistency.
class LogQueryTask(Task):
    """
    Log Query Task
    """

    @property
    def name(self):
        return "logquery"

    @property
    def title(self):
        return "Log queries"

    @property
    def programs(self):
        return ()

    @property
    def paramnames(self):
        return self.join_params(super().paramnames,
                                {"checkdlp": True,
                                 "jlpoly": False,
                                 "gfpext": [int, 1],
                                 "N": int,
                                 "ell": 0})

    def get_logquery_filename(self):
        return self.get_state_filename("logquery")

    def get_logquery_checker(self):
        return self

    def get_logbase(self):
        return self.logbase

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)

        self.history = dict()
        if "logbase" in self.state:
            self.logbase = self.state["logbase"]
        else:
            self.logbase = None
        self.p = self.params["N"]
        self.ell = self.params["ell"]
        self.gfpext = self.params["gfpext"]
        # TODO gfpext
        assert (self.p ** self.gfpext - 1) % self.ell == 0
        self.cof = (self.p ** self.gfpext - 1) // self.ell

        if "logquery" not in self.state:
            logquery_filename = self.workdir.make_filename("logquery")
            update = {"logquery": logquery_filename.get_wdir_relative()}
            self.state.update(update, commit=True)

    def get_small_logs(self, bound=20):
        small_logs = dict()
        with open(str(self.send_request(Request.GET_DLOG_FILENAME)), "r") as f:
            for line in f:
                mm = re.match(r'(\w+) (\w+) \d+ rat (\d+)', line)
                if mm:
                    target = int(mm.group(2), 16)
                    logtarget = int(mm.group(3))
                    if target > bound:
                        break
                    small_logs[target] = logtarget
        return small_logs

    def commit_logs(self, *args):
        with open(str(self.get_logquery_filename()), "a") as f:
            for x in args:
                target, logtarget = x
                f.write("%d %d\n" % (target, logtarget))

    def check_new_log(self, target, logtarget, commit=True):
        if target in self.history:
            return
        if logtarget == 0:
            msg = "Checking that log of %d is zero..." % target
            check = pow(target, self.cof, self.p) == 1
            if check:
                self.logger.info(msg + " passed")
            else:
                self.logger.critical(msg + " FAILED")
                raise ValueError("Failed log check, log(%d)=0 seems wrong\n"
                                 % target)
            return
        just_deduced_gen = False
        if self.logbase is None:
            gt, ilogt, foo = xgcd(logtarget * self.cof, self.ell)
            ilogt = ilogt % self.ell
            if gt == 1:
                # then target^((p-1)/ell * ilogt) is a generator
                self.logbase = pow(int(target), ilogt*self.cof, self.p)
                self.state.update({"logbase": self.logbase}, commit=True)
                context = "log(%d)=%d" % (target, logtarget)
                conclusion = "logarithms are given in base %d" % self.logbase

                self.logger.info("Based on %s, we expect that %s"
                                 % (context, conclusion))
                just_deduced_gen = True
        if not just_deduced_gen:
            msg = "Checking that log(%d)=%d is correct in base %d..." \
                  % (target, logtarget, self.logbase)
            self.logger.info(msg)

            left = pow(self.logbase, logtarget*self.cof, self.p)
            right = pow(target, self.cof, self.p)

            if left == right:
                self.logger.info(msg + " passed")
            else:
                self.logger.critical(msg + " FAILED")
                raise ValueError("Failed log check,"
                                 " log(%d)=%d seems wrong in base %d\n"
                                 % (target, logtarget, self.logbase))
        else:
            self.logger.info("Check skipped for this log,"
                             " as it is the only known log value"
                             " at this point")

        self.history[target] = logtarget
        if commit:
            self.commit_logs((target, logtarget))

    def run(self):
        super().run()

        # read the logs that we already queried.
        if self.get_logquery_filename().isfile():
            with open(str(self.get_logquery_filename()), "r") as f:
                for line in f:
                    if re.match(r'^#', line):
                        continue
                    mm = re.match(r"(\d+) (\d+)", line)
                    if not mm:
                        self.logger.warning("Unparsed line in %s: "
                                            % self.get_logquery_filename(),
                                            line)
                    target = int(mm.group(1))
                    logtarget = int(mm.group(2))
                    self.check_new_log(target, logtarget, commit=False)

        for target, logtarget in self.get_small_logs().items():
            self.check_new_log(target, logtarget)

        return True


class ClGroupStructureTask(Task):
    """ Compute the group structure for Cl """
    @property
    def name(self):
        return "groupstructure"

    @property
    def title(self):
        return "CL Group Structure"

    @property
    def programs(self):
        return ((cadoprograms.ClStructure, (),
                 {"poly": Request.GET_POLYNOMIAL_FILENAME,
                  "order": Request.GET_CANDIDATE_CLASSNUMBER_FACTORED}),)

    @property
    def paramnames(self):
        return self.join_params(super().paramnames,
                                {"verbose": [bool],
                                 "generators_bound": [int],
                                 "pSylow_bound": [int]})

    def run(self):
        super().run()
        order = self.merged_args[0]['order']
        if not self.has_all_data_in_state(order):
            self.state.clear()
            self.state["order"] = order
            (stdoutpath, stderrpath) = self.make_std_paths(
                cadoprograms.ClStructure.name)
            p = cadoprograms.ClStructure(stdout=str(stdoutpath),
                                         stderr=str(stderrpath),
                                         **self.merged_args[0])
            message = self.submit_command(p, None, log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")
            # TODO if program failed due to too large problematic p^e in the
            # group order, try to recover using linear algebra modulo p

            self.state["classgroup_desc"] = ""
            stdout = message.read_stdout(0).decode("utf-8")
            for line in stdout.splitlines():
                if line.startswith("structure: "):
                    self.state["classgroup_desc"] += line[11:] + "\n"

            if not self.state["classgroup_desc"]:
                self.logger.error("No information on the class group "
                                  "structure was found in the output file")
                raise Exception("Parsing failed or uncaught error")

        self.logger.debug("Exit ClGroupStructureTask.run(" + self.name + ")")
        return True

    def has_all_data_in_state(self, order):
        return "order" in self.state and self.state["order"] == order \
            and "classgroup_desc" in self.state

    def get_class_group_structure(self):
        return self.state.get("classgroup_desc", None)


class FactorTask(Task):
    """
    Use the binary misc/factor to **try** to factor completly a integer
    Beware: it is not the main task (which is CompleteFactorization) but a task
    used to factor potential class number.
    """

    @property
    def name(self):
        return "factor"

    @property
    def title(self):
        return "Factor"

    @property
    def programs(self):
        return ((cadoprograms.Factor, (), {"input": self._input_req}),)

    @property
    def paramnames(self):
        return self.join_params(super().paramnames, {"verbose": [bool]})

    def __init__(self, input_req, *, mediator, db, parameters, path_prefix):
        self._input_req = input_req
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)

    def run(self):
        super().run()
        input = self.merged_args[0]['input']
        if not self.has_all_data_in_state(input):
            self.state.clear()
            self.state['input'] = input
            (stdoutpath, stderrpath) = self.make_std_paths(
                cadoprograms.Factor.name)
            p = cadoprograms.Factor(stdout=str(stdoutpath),
                                    stderr=str(stderrpath),
                                    **self.merged_args[0])
            message = self.submit_command(p, None, log_errors=True)
            if message.get_exitcode(0) != 0:
                raise Exception("Program failed")

            stdout = message.read_stdout(0).decode("utf-8")
            for line in stdout.splitlines():
                if line.startswith("factorization: "):
                    self.state["fact_str"] = line[15:].replace(" ", "")
                    break
            else:
                raise Exception("Could not find pattern '^factorization: '")

        self.logger.info(f"factorization of {input}: {self.state['fact_str']}")
        self.logger.debug("Exit Factor.run(" + self.name + ")")
        return True

    def has_all_data_in_state(self, n):
        return 'input' in self.state and 'fact_str' in self.state \
                and self.state['input'] == n

    def get_fact_str(self):
        return self.state.get('fact_str', None)


class StartServerTask(DoesLogging, cadoparams.UseParameters, HasState):
    """
    Starts HTTP server
    """

    @property
    def name(self):
        return "server"

    @property
    def title(self):
        return "Server Launcher"

    @property
    def paramnames(self):
        return {"name": str, "workdir": None, "address": None, "port": 0,
                "threaded": False, "ssl": True, "whitelist": None,
                "only_registered": True, "forgetport": False,
                "timeout_hint": None, "nrsubdir": 0,
                "linger_before_quit": 0}

    @property
    def param_nodename(self):
        return self.name

    # The whitelist parameter here is an iterable of strings in CIDR
    # notation.
    # The whitelist parameter we get from params (i.e., from the
    # parameter file) is a string with comma-separated CIDR strings. The
    # two are concatenated to form the server whitelist.
    def __init__(self, *,
                 default_workdir,
                 parameters,
                 path_prefix,
                 db,
                 whitelist=None):
        super().__init__(db=db,
                         parameters=parameters,
                         path_prefix=path_prefix)
        self.params = self.parameters.myparams(self.paramnames)
        serveraddress = self.params.get("address", None)
        serverport = self.params["port"]
        basedir = self.params.get("workdir", default_workdir)
        basedir = basedir.rstrip(os.sep) + os.sep
        uploaddir = basedir + self.params["name"] + ".upload/"
        threaded = self.params["threaded"]
        nrsubdir = self.params["nrsubdir"]
        # By default, allow access only to files explicitly registered by
        # tasks, i.e., those files required by clients when downloading
        # input files for their workunits. By setting
        # only_registered=False, access to all files under the server
        # working directory is allowed.
        only_registered = self.params["only_registered"]
        if self.params["ssl"]:
            cafilename = basedir + self.params["name"] + ".server.cert"
        else:
            cafilename = None

        servertimeout_hint = self.params.get("timeout_hint")

        server_whitelist = []
        if whitelist is not None:
            server_whitelist += whitelist
        if "whitelist" in self.params:
            server_whitelist += [h.strip()
                                 for h in self.params["whitelist"].split(",")]

        # If we should auto-assign an available port, try to use the same one
        # as last time, if we had run before. This can be overridden with
        # server.forgetport=yes
        if self.params["forgetport"] and "port" in self.state:
            del self.state["port"]

        if serverport == 0 and "port" in self.state:
            serverport = self.state["port"]

        # If (1) any clients are to be started on localhost, but (2) the
        # server is listening on a network-visible address, then we need
        # to whitelist the network-visible address(es) of the current
        # host as well, because the client's connection will come from
        # (one of) the network-visible addresses of the current host.
        # For test (1), it should suffice to look for "localhost" or
        # "127.0.0.1" in the existing whitelist, as the host names on
        # which to start clients are inserted verbatim.
        # For (2), we check that serveraddress is either None (i.e., the
        # wildcard address which is network-visible), or anything other
        # than "localhost"
        if (serveraddress is None
            or socket.getfqdn(serveraddress) != "localhost") \
           and set(server_whitelist) & {"localhost", "127.0.0.1"}:
            hostname = socket.gethostname()
            if hostname not in server_whitelist:
                try:
                    # foo=socket.gethostbyname(hostname)
                    self.logger.info("Adding %s to whitelist"
                                     " to allow clients on localhost"
                                     " to connect", hostname)
                    server_whitelist.append(hostname)
                except socket.gaierror:
                    self.logger.info("Not adding %s to whitelist"
                                     " (cannot be resolved),"
                                     " clients will only be allowed to"
                                     " connect on 127.0.0.1", hostname)
        if not server_whitelist:
            server_whitelist = None

        self.registered_filenames = \
            self.make_db_dict('server_registered_filenames')

        # lbq = self.params["linger_before_quit"]

        self.server = ApiServer(
            serveraddress, serverport, db,
            threaded=threaded,
            uploaddir=uploaddir,
            nrsubdir=nrsubdir,
            only_registered=only_registered,
            cafile=cafilename,
            whitelist=server_whitelist,
            timeout_hint=servertimeout_hint,
            # linger_before_quit=lbq
            )
        self.state["port"] = self.server.get_port()

    def run(self):
        self.server.serve()

    def shutdown(self, *args):
        self.server.shutdown(*args)

    def stop_serving_wus(self):
        self.server.stop_serving_wus()

    def get_url(self, **kwargs):
        return self.server.get_url(**kwargs)

    def get_cert_sha1(self):
        return self.server.get_cert_sha1()

    def register_filename(self, d):
        for key in d:
            if key not in self.registered_filenames:
                self.logger.debug("Registering file name %s with target %s",
                                  key, d[key])
                self.registered_filenames[key] = d[key]
            elif d[key] != self.registered_filenames[key]:
                # It was already registered with a different target. This
                # will happen if, e.g., the user chooses a different
                # build directory between runs. It's still fragile, as
                # the server will try to serve the old file(s) until a
                # new workunit is generated and overrides the target.
                # The proper solution would be to make Program classes
                # Templates, so Tasks can instantiate them at __init__()
                # and register the resolved paths once. The
                # registered_filenames dict could then be memory-backed.
                # Tasks would also have to register their own files which
                # may need to be served in __init__().
                self.logger.warning("Filename %s,"
                                    " to be registered for target %s,"
                                    " was previously registered for target %s."
                                    " Overriding with new target.",
                                    key,
                                    d[key],
                                    self.registered_filenames[key])
                self.registered_filenames[key] = d[key]
            else:
                # Was already registered with the same target. Nothing to do
                pass


class StartClientsTask(Task):
    """
    Starts clients on slave machines
    """

    @property
    def name(self):
        return "slaves"

    @property
    def title(self):
        return "Client Launcher"

    @property
    def programs(self):
        return ((cadoprograms.CadoNFSClient, ("clientid", "certsha1"), {}),)

    @property
    def paramnames(self):
        return {'hostnames': str,
                'scriptpath': None,
                "nrclients": [int],
                "run": True}

    @property
    def param_nodename(self):
        return None

    def __init__(self, *, mediator, db, parameters, path_prefix):
        super().__init__(mediator=mediator, db=db, parameters=parameters,
                         path_prefix=path_prefix)
        self.used_ids = {}
        self.pids = self.make_db_dict(self.make_tablename("client_pids"),
                                      connection=self.db_connection)
        self.hosts = self.make_db_dict(self.make_tablename("client_hosts"),
                                       connection=self.db_connection)
        self.servertask = None
        assert set(self.pids) == set(self.hosts)
        # Invariants: the keys of self.pids and of self.hosts are the
        # same set.  The keys of self.used_ids are a subset of the keys
        # of self.pids.  A clientid is in self.used_ids if we know that
        # clientid to be currently running.

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
            self.hosts_to_launch = \
                self.make_multiplicity(self.hosts_to_launch,
                                       self.params["nrclients"])

    @staticmethod
    def make_multiplicity(names, multi):
        """
        Produce a list in which each unique entry of the list "names"
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
        """
        Get list host names on which clients should run
        """
        return self.hosts_to_launch

    def is_alive(self, clientid):
        # Simplistic: just test if process with that pid exists and
        # accepts signals from us. TODO: better testing here, probably
        # with ps|grep or some such
        # (rc, stdout, stderr) = self.kill_client(clientid, signal=0)
        (rc, stdout, stderr) = self.ping_client(clientid)
        return (rc == 0, stdout, stderr)

    def _add_cid(self, clientid, pid, host):
        """
        Add a client id atomically to both the "pids" and "hosts"
        dictionaries
        """
        self.pids.update({clientid: pid}, commit=False)
        self.hosts.update({clientid: host}, commit=True)

    def _del_cid(self, clientid):
        """
        Remove a client id atomically from both the "pids" and "hosts"
        dictionaries
        """
        self.pids.clear([clientid], commit=False)
        self.hosts.clear([clientid], commit=True)

    def launch_clients(self, servertask):
        """
        This now takes server as a servertask object, so that we can
        get an URL which is special-cased for localhost
        """
        self.servertask = servertask
        url = servertask.get_url()
        url_loc = servertask.get_url(origin="localhost")
        certsha1 = servertask.get_cert_sha1()
        for host in self.hosts_to_launch:
            if host == "localhost":
                self.launch_one_client(host.strip(),
                                       url_loc,
                                       certsha1=certsha1)
            else:
                self.launch_one_client(host.strip(),
                                       url,
                                       certsha1=certsha1)
        running_clients = [(cid, self.hosts[cid], pid)
                           for (cid, pid) in self.pids.items()]
        s = ", ".join(["%s (Host %s, PID %d)" % t
                       for t in running_clients])
        self.logger.info("Running clients: %s" % s)
        # Check for old clients which we did not mean to start this run
        for cid in set(self.pids) - set(self.used_ids):
            if self.is_alive(cid)[0]:
                self.logger.warning("Client id %s (Host %s, PID %d), launched"
                                    " in a previous run and not meant to be"
                                    " launched this time, is still running",
                                    cid, self.hosts[cid], self.pids[cid])
            else:
                self.logger.warning("Client id %s (Host %s, PID %d), launched"
                                    " in a previous run and not meant to be"
                                    " launched this time, seems to have died."
                                    " I'll forget about this client.",
                                    cid, self.hosts[cid], self.pids[cid])
                self._del_cid(cid)

    def check_health(self, localhost=True):
        """
        This checks that clients are alive, and raises an exception
        if they aren't. A priori, it is reasonable to do this only on
        localhost, but optionally we support the idea of making this sort
        of check on all remote hosts as well.
        """
        look = [k
                for k, v in self.hosts.items()
                if not localhost or v == 'localhost']
        for cid in look:
            alive, stdout, stderr = self.is_alive(cid)
            h = self.hosts[cid]
            p = self.pids[cid]
            who = "client id %s (Host %s, PID %d)" % (cid, h, p)
            # self.logger.debug("check %s: %s", who, alive)
            if not alive:
                self.logger.critical("DEAD: " + who)
                if stdout:
                    self.logger.critical("Stdout: %s",
                                         stdout.decode().strip())
                if stderr:
                    self.logger.critical("Stderr: %s",
                                         stderr.decode().strip())
                if localhost:
                    raise RuntimeError("localhost clients should never die")
                else:
                    raise RuntimeError("clients should never die")

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
            if self.is_alive(clientid)[0]:
                self.logger.info("Client %s on host %s with PID %d already "
                                 "running",
                                 clientid, host, self.pids[clientid])
                self.used_ids[clientid] = True
                return
            else:
                self.logger.info("Client %s on host %s"
                                 " with PID %d seems to have died",
                                 clientid, host, self.pids[clientid])
                self._del_cid(clientid)

        self.logger.info("Starting client id %s on host %s", clientid, host)
        cado_nfs_client = cadoprograms.CadoNFSClient(server=server,
                                                     clientid=clientid,
                                                     daemon=True,
                                                     certsha1=certsha1,
                                                     **self.progparams[0])
        if host == "localhost":
            process = cadocommand.Command(cado_nfs_client)
        else:
            process = cadocommand.RemoteCommand(cado_nfs_client,
                                                host,
                                                self.parameters)
        (rc, stdout, stderr) = process.wait()
        if rc != 0:
            self.logger.warning("Starting client on host %s failed.", host)
            if stdout:
                self.logger.warning("Stdout: %s", stdout.decode().strip())
            if stderr:
                self.logger.warning("Stderr: %s", stderr.decode().strip())
            return
        match = None
        if stdout is not None:
            # The client doesn't print much in daemon mode, but still
            # does print a little. re.search is better than re.match!
            match = re.search(r"PID: (\d+)", stdout.decode())
        if not match:
            self.logger.warning("Client did not print PID")
            if stdout is not None:
                self.logger.warning("Stdout: %s", stdout.decode().strip())
            if stderr is not None:
                self.logger.warning("Stderr: %s", stderr.decode().strip())
            return
        self.used_ids[clientid] = True
        self._add_cid(clientid, int(match.group(1)), host)

    def kill_all_clients(self):
        # Need the list() to make a copy as dict will change in loop body
        for clientid in list(self.pids):
            (rc, stdout, stderr) = self.kill_client(clientid)
            if rc == 0:
                self.logger.info("Stopped client %s (Host %s, PID %d)",
                                 clientid,
                                 self.hosts[clientid],
                                 self.pids[clientid])
                self._del_cid(clientid)
            else:
                self.logger.warning("Stopping client %s"
                                    " (Host %s, PID %d) failed",
                                    clientid,
                                    self.hosts[clientid],
                                    self.pids[clientid])
                if stdout:
                    self.logger.warning("Stdout: %s", stdout.decode().strip())
                if stderr:
                    self.logger.warning("Stderr: %s", stderr.decode().strip())
                # Assume that the client is already dead and remove it from
                # the list of running clients
                self._del_cid(clientid)

    def ping_client(self, clientid):
        """
        runs cado-nfs-client to see if the client that we started is
        still running, and try to recover its output if it died.

        Old behaviour was to kill -0 the pid, but of course that doesn't
        give us the error log.
        """
        pid = self.pids[clientid]
        host = self.hosts[clientid]
        url = self.servertask.get_url()
        url_loc = self.servertask.get_url(origin="localhost")
        certsha1 = self.servertask.get_cert_sha1()
        if host == "localhost":
            server = url_loc
            ping = cadoprograms.CadoNFSClient(server=server,
                                              clientid=clientid,
                                              ping=pid,
                                              daemon=True,
                                              certsha1=certsha1,
                                              **self.progparams[0])
            process = cadocommand.Command(ping)
        else:
            server = url
            ping = cadoprograms.CadoNFSClient(server=server,
                                              clientid=clientid,
                                              ping=pid,
                                              daemon=True,
                                              certsha1=certsha1,
                                              **self.progparams[0])
            process = cadocommand.RemoteCommand(ping, host, self.parameters)
        return process.wait()

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
    GET_WILL_IMPORT_FINAL_POLYNOMIAL = object()
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
    GET_GAL_UNIQUE_RELCOUNT = object()
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
    GET_VIRTUAL_LOGS_FILENAME = object()
    GET_RELSDEL_FILENAME = object()
    GET_SM_FILENAME = object()
    GET_UNITS_DIRNAME = object()
    GET_NMAPS = object()
    GET_WU_RESULT = object()
    GET_WORKDIR_JOBNAME = object()
    GET_WORKDIR_PATH = object()
    GET_DLOG_FILENAME = object()
    GET_LOGQUERY_FILENAME = object()
    GET_LOGQUERY_CHECKER = object()
    GET_CLIENTS = object()
    GET_CANDIDATE_CLASSNUMBER = object()
    GET_CANDIDATE_CLASSNUMBER_FACTORED = object()


class CompleteFactorization(HasState,
                            # (implicit) wudb.UsesWorkunitDb,
                            wudb.ListensToWorkunitDb,
                            DoesLogging,
                            cadoparams.UseParameters,
                            patterns.Mediator):
    """
    The complete factorization / dlp, aggregate of the individual tasks
    """

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
        return {"name": str,
                "workdir": str,
                "N": int,
                "ell": 0,
                "algo": Algorithm.NFS.name,
                "computation": Computation.FACT.name,
                "gfpext": 1,
                "jlpoly": False,
                "trybadwu": False,
                "target": ""}

    @property
    def title(self):
        try:
            if self.params["computation"] == Computation.DLP:
                return "Discrete logarithm"
            elif self.params["computation"] == Computation.CL:
                return "Class group"
            else:
                return "Complete Factorization"
        except AttributeError:
            return "Complete Factorization / Discrete logarithm / Class group"

    @property
    def programs(self):
        return []

    def get_clients(self):
        return self.clients

    def __init__(self, db, parameters, path_prefix):
        self.db = db
        super().__init__(db=db, parameters=parameters, path_prefix=path_prefix)
        self.params = self.parameters.myparams(self.paramnames)

        # check compatibility
        Incompatible = (
            (Computation.CL, Algorithm.NFS),
            (Computation.DLP, Algorithm.QS),
            )
        if (self.params["computation"], self.params["algo"]) in Incompatible:
            raise ValueError(f"Algorithm {self.params['computation']} is "
                             f"incompatible with {self.params['algo']}")

        if self.params["computation"] == Computation.DLP:
            p = self.params["N"]
            k = self.params["gfpext"]
            ell = self.params["ell"]
            # Check that ell divides Phi_k(p) (or simply p^k-1 for large k)
            if k == 1:
                if (p - 1) % ell != 0:
                    raise ValueError("ell must divide p-1")
            elif k == 2:
                if (p + 1) % ell != 0:
                    raise ValueError("ell must divide p+1")
            elif k == 3:
                if (p*p + p + 1) % ell != 0:
                    raise ValueError("ell must divide p^2+p+1")
            elif k == 4:
                if (p*p + 1) % ell != 0:
                    raise ValueError("ell must divide p^2+1")
            elif k == 5:
                if (p**4 + p**3 + p**2 + p + 1) % ell != 0:
                    raise ValueError("ell must divide (p^5-1)/(p-1)")
            elif k == 6:
                if (p**2 - p + 1) % ell != 0:
                    raise ValueError("ell must divide p^2-p+1")
            elif (p**k - 1) % ell != 0:
                raise ValueError("ell must divide p^%d-1" % k)
        if self.params["computation"] == Computation.CL:
            disc = self.params["N"]
            if disc >= 0:
                raise ValueError("discriminant must be negative")
            elif disc % 16 in (0, 4):
                raise ValueError("discriminant divisible by 4 must be 8 or 12 "
                                 "mod 16 (discriminant/4 must be 2 or 3)")
            # Set tasks.linalg.nmatrices to 2 by default
            parameters.set_if_unset("tasks.nmatrices", 2)

        if self.params["algo"] == Algorithm.QS:
            # set galois to _y; it must be unset or already set to _y
            if parameters.set_if_unset("tasks.galois", "_y") != "_y":
                raise Exception("For QS algorithm, parameter \"galois\" "
                                "must be unset or set to _y")

        # Init WU BD
        # Note that we get self.wuar = self.make_wu_access(db.connect()) via
        # inheritance of wudb.UsesWorkunitDb
        # self.wuar = self.make_wu_access()

        self.wuar.create_tables()
        if self.params["trybadwu"]:
            # Test behaviour when a WU is in the DB that does not belong to
            # any task. It should get cancelled with an error message.
            self.wuar.create(["WORKUNIT FAKE_WU_%s\n"
                              "COMMAND true"
                              % time.time()])

        # Start with an empty list of tasks that want to run. Tasks will add
        # themselves during __init__().
        self.tasks_that_want_to_run = list()

        # Init client lists
        self.clients = []
        whitelist = set()
        start_list = self.parameters.get_parameters().find(['slaves'],
                                                           'hostnames')
        for (path, key) in start_list:
            self.clients.append(StartClientsTask(mediator=self,
                                                 db=db,
                                                 parameters=self.parameters,
                                                 path_prefix=path))
            hostnames = self.clients[-1].get_hosts_to_launch()
            whitelist |= set(hostnames)

        whitelist = list(whitelist) if whitelist else None
        # Init server task
        self.servertask = \
            StartServerTask(default_workdir=self.params["workdir"],
                            parameters=parameters,
                            path_prefix=path_prefix,
                            db=db,
                            whitelist=whitelist)

        # Create tasks
        parampath = self.parameters.get_param_path()
        computation = self.params["computation"]
        algo = self.params["algo"]
        tasks_kwargs = {'mediator': self, 'db': db,
                        'parameters': self.parameters}

        # tasks.polyselect
        tasks_kwargs['path_prefix'] = parampath + ['polyselect']
        if algo == Algorithm.QS:
            self.polysel = (PolyselQSTask(**tasks_kwargs), )
        elif computation == Computation.DLP and self.params["gfpext"] != 1:
            self.polysel = (PolyselGFpnTask(**tasks_kwargs), )
        elif computation == Computation.DLP and self.params["jlpoly"]:
            self.polysel = (PolyselJLTask(**tasks_kwargs), )
        else:  # default dlp with gfpext = 1 and factorization
            self.polysel = (Polysel1Task(**tasks_kwargs),
                            Polysel2Task(**tasks_kwargs))

        # tasks.numbertheory
        tasks_kwargs['path_prefix'] = parampath + ['numbertheory']
        if computation == Computation.DLP:
            self.numbertheory = NumberTheoryTask(**tasks_kwargs)

        # tasks.sieve
        tasks_kwargs['path_prefix'] = parampath + ['sieve']
        nsides = 2 if not algo == Algorithm.QS else 1
        self.fb = FactorBaseTask(nsides, **tasks_kwargs)
        if computation == Computation.CL:
            self.checkdisc = CheckDiscriminantTask(**tasks_kwargs)
        self.freerel = FreeRelTask(nsides, **tasks_kwargs)
        if algo == Algorithm.QS:
            self.sieving = QuadraticSievingTask(nsides, **tasks_kwargs)
        else:
            self.sieving = SievingTask(nsides, **tasks_kwargs)

        # tasks.filter
        tasks_kwargs['path_prefix'] = parampath + ['filter']
        self.dup1 = Duplicates1Task(**tasks_kwargs)
        self.dup2 = Duplicates2Task(**tasks_kwargs)
        if computation == Computation.DLP or algo == Algorithm.QS:
            self.filtergalois = FilterGaloisTask(**tasks_kwargs)
        self.purge = PurgeTask(**tasks_kwargs)
        if computation == Computation.FACT:
            self.merge = MergeTask(**tasks_kwargs)
        else:
            self.merge = MergeDLPTask(**tasks_kwargs)
        if computation == Computation.DLP:
            self.sm = SMTask(**tasks_kwargs)

        # tasks.linalg
        tasks_kwargs['path_prefix'] = parampath + ['linalg']
        if computation == Computation.DLP:
            self.linalg = LinAlgDLPTask(**tasks_kwargs)
        elif computation == Computation.CL:
            self.linalg = LinAlgClTask(**tasks_kwargs)
        else:  # factorization
            self.linalg = LinAlgTask(**tasks_kwargs)
            self.characters = CharactersTask(**tasks_kwargs)

        # tasks.sqrt
        tasks_kwargs['path_prefix'] = parampath + ['sqrt']
        if computation == Computation.FACT:
            self.sqrt = SqrtTask(**tasks_kwargs)

        # tasks.reconstructlog
        tasks_kwargs['path_prefix'] = parampath + ['reconstructlog']
        if computation == Computation.DLP:
            self.reconstructlog = ReconstructLogTask(**tasks_kwargs)

        # tasks.logquery
        tasks_kwargs['path_prefix'] = parampath + ['logquery']
        if computation == Computation.DLP:
            # By specifying a single endpoint, we make it easier to do
            # consistency checks.
            self.logquery = LogQueryTask(**tasks_kwargs)

        # tasks.descent
        tasks_kwargs['path_prefix'] = parampath + ['descent']
        if computation == Computation.DLP:
            self.descent = DescentTask(**tasks_kwargs)

        # tasks.groupstructure
        tasks_kwargs['path_prefix'] = parampath + ['groupstructure']
        if computation == Computation.CL:
            self.hfactor = FactorTask(Request.GET_CANDIDATE_CLASSNUMBER,
                                      **tasks_kwargs)
            self.grstruct = ClGroupStructureTask(**tasks_kwargs)

        # Defines an order on tasks in which tasks that want to run should be
        # run
        if computation == Computation.DLP:
            self.tasks = self.polysel \
                + (self.numbertheory, self.fb, self.freerel, self.sieving,
                   self.dup1, self.dup2, self.filtergalois, self.purge,
                   self.merge, self.sm, self.linalg, self.reconstructlog,
                   self.logquery)
            if self.params["target"]:
                self.tasks = self.tasks + (self.descent,)
        elif computation == Computation.CL:
            self.tasks = self.polysel \
                + (self.fb, self.freerel, self.checkdisc, self.sieving,
                   self.dup1, self.dup2, self.filtergalois, self.purge,
                   self.merge, self.linalg, self.hfactor, self.grstruct)
        else:
            fg = () if algo == Algorithm.NFS else (self.filtergalois, )
            self.tasks = self.polysel \
                + (self.fb, self.freerel, self.sieving, self.dup1, self.dup2) \
                + fg + (self.purge, self.merge, self.linalg) \
                + (self.characters, self.sqrt)

        reverse_lookup = defaultdict(list)
        self.parameter_help = ""
        for t in self.tasks:
            self.parameter_help += t.collect_usable_parameters(reverse_lookup)

        for (path, key, value) in parameters.get_unused_parameters():
            self.logger.warning("Parameter %s = %s was not used anywhere",
                                ".".join(path + [key]), value)
            if key in reverse_lookup.keys():
                ell = reverse_lookup[key]
                if len(ell) == 1:
                    self.logger.warning("Perhaps you meant %s.%s ?" %
                                        (ell[0], key))
                else:
                    self.logger.warning("Perhaps you meant"
                                        " one of the following ?")
                    for x in ell:
                        self.logger.warning("  %s.%s ?" % (x, key))
                    prefix = ".".join(os.path.commonprefix([x.split(".")
                                                            for x in ell]))
                    self.logger.warning("(If you wish to set all of these"
                                        " consistently, you may set %s.%s)"
                                        % (prefix, key))

        self.request_map = {
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
            Request.GET_WU_RESULT: self.send_result,
            Request.GET_WORKDIR_JOBNAME: self.fb.workdir.get_workdir_jobname,
            Request.GET_WORKDIR_PATH: self.fb.workdir.get_workdir_path,
            Request.GET_CLIENTS: self.get_clients,
        }

        # Set requests related to polynomial selection
        if len(self.polysel) == 1:
            self.request_map[Request.GET_POLYNOMIAL] = \
                self.polysel[0].get_poly
            self.request_map[Request.GET_POLYNOMIAL_FILENAME] = \
                self.polysel[0].get_poly_filename
        else:
            self.request_map[Request.GET_RAW_POLYNOMIALS] = \
                self.polysel[0].get_raw_polynomials
            self.request_map[Request.GET_POLY_RANK] = \
                self.polysel[0].get_poly_rank
            self.request_map[Request.GET_POLYNOMIAL] = \
                self.polysel[1].get_poly
            self.request_map[Request.GET_POLYNOMIAL_FILENAME] = \
                self.polysel[1].get_poly_filename
            self.request_map[Request.GET_WILL_IMPORT_FINAL_POLYNOMIAL] = \
                self.polysel[1].get_will_import

        # Add requests specific the type of computation
        if computation == Computation.DLP or algo == Algorithm.QS:
            self.request_map[Request.GET_GAL_UNIQUE_RELCOUNT] = \
                self.filtergalois.get_nrels

        if computation == Computation.DLP:
            self.request_map[Request.GET_IDEAL_FILENAME] = \
                self.merge.get_ideal_filename
            self.request_map[Request.GET_NMAPS] = \
                self.numbertheory.get_nmaps
            self.request_map[Request.GET_SM_FILENAME] = \
                self.sm.get_sm_filename
            self.request_map[Request.GET_RELSDEL_FILENAME] = \
                self.purge.get_relsdel_filename
            self.request_map[Request.GET_DLOG_FILENAME] = \
                self.reconstructlog.get_dlog_filename
            self.request_map[Request.GET_LOGQUERY_FILENAME] = \
                self.logquery.get_logquery_filename
            self.request_map[Request.GET_LOGQUERY_CHECKER] = \
                self.logquery.get_logquery_checker
            self.request_map[Request.GET_KERNEL_FILENAME] = \
                self.linalg.get_virtual_logs_filename
            self.request_map[Request.GET_VIRTUAL_LOGS_FILENAME] = \
                self.linalg.get_virtual_logs_filename
        elif computation == Computation.FACT:
            self.request_map[Request.GET_KERNEL_FILENAME] = \
                self.characters.get_kernel_filename
            self.request_map[Request.GET_DEPENDENCY_FILENAME] = \
                self.linalg.get_dependency_filename
            self.request_map[Request.GET_LINALG_PREFIX] = \
                self.linalg.get_prefix
        elif computation == Computation.CL:
            self.request_map[Request.GET_CANDIDATE_CLASSNUMBER] = \
                self.linalg.get_candidate_class_number
            self.request_map[Request.GET_CANDIDATE_CLASSNUMBER_FACTORED] = \
                self.hfactor.get_fact_str

    def enter_subtask_chain(self):
        self.start_elapsed_time()
        self.servertask.run()
        self.start_all_clients()

    def exit_subtask_chain(self, exc):
        self.servertask.stop_serving_wus()
        # print everybody's stats before we exit.
        for task in self.tasks_that_have_run:
            task.print_stats()
        self.stop_all_clients()
        self.elapsed = self.end_elapsed_time()
        self.cputotal = self.get_sum_of_cpu_or_real_time(True)
        self.servertask.shutdown(exc)

    def run(self):
        if self.params["computation"] == Computation.DLP:
            self.logger.info("Computing Discrete Logs in GF(%s)",
                             self.params["N"])
        elif self.params["computation"] == Computation.CL:
            self.logger.info("Computing class number for %s", self.params["N"])
        else:
            if self.params["algo"] == Algorithm.QS:
                self.logger.info("Factoring (with QS) %s", self.params["N"])
            else:
                self.logger.info("Factoring %s", self.params["N"])

        class wrapme(object):
            def __init__(self, s):
                self.s = s

            def __enter__(self):
                self.s.enter_subtask_chain()
                return self

            def __exit__(self, e_type, e_value, traceback):
                self.s.exit_subtask_chain(e_value)

        last_task = None
        last_status = True
        # we rely here on Task not having a weird comparison operator
        self.tasks_that_have_run = set()
        try:
            with wrapme(self):
                while last_status:
                    task = self.next_task()
                    if task is None:
                        break
                    last_task = task.title
                    last_status = task.run()
                    self.tasks_that_have_run.add(task)
                    self.logger.info(task.title)
                    task.print_stats()

        except KeyboardInterrupt:
            self.logger.fatal("Received KeyboardInterrupt. Terminating")
            return None

        except EarlyStopException as e:
            self.logger.info("Total cpu/elapsed time"
                             " for incomplete %s: %g/%g",
                             self.title, self.cputotal, self.elapsed)
            self.logger.info("Finishing early: " + str(e))
            raise e

        # Do we want the sum of real times over all sub-processes for
        # something?
        # realtotal = self.get_sum_of_cpu_or_real_time(False)
        if self.params["computation"] == Computation.DLP:
            self.logger.info("Total cpu/elapsed time"
                             " for entire %s: %g/%g",
                             self.title, self.cputotal, self.elapsed)
        else:
            self.logger.info("Total cpu/elapsed time"
                             " for entire %s %g/%g" + tstr(self.elapsed),
                             self.title, self.cputotal, self.elapsed)

        if last_task and not last_status:
            self.logger.fatal("Premature exit within %s. Bye.", last_task)
            return None

        # Print the defining polynomial of the finite field used for
        # representing elements.
        # This assumes that the last line of the poly file contains this
        # information. This is currently the case for polyselect_gfpn.c
        # but of course, this won't be the case for a user-defined poly
        # file that has been imported (anyway, in that case, the user
        # should know what she is doing).
        if self.params["computation"] == Computation.DLP and \
                self.params["gfpext"] > 1:
            polyfile = self.request_map[Request.GET_POLYNOMIAL_FILENAME]()
            with open(str(polyfile), "r") as ff:
                s = ff.read().splitlines()[-1].split()[-1]
                self.logger.info("The polynomial defining"
                                 " the finite field is %s", s)

        if self.params["computation"] == Computation.DLP:
            if self.params["target"]:
                logt = self.descent.get_logtargets()
                base = self.logquery.get_logbase()
                return [base] + logt
            else:
                return [0]
        elif self.params["computation"] == Computation.CL:
            return self.grstruct.get_class_group_structure()
        else:
            return self.sqrt.get_factors()

    def start_all_clients(self):
        for clients in self.clients:
            clients.launch_clients(self.servertask)

    def stop_all_clients(self):
        for clients in self.clients:
            clients.kill_all_clients()

    def start_elapsed_time(self):
        if "starttime" in self.state:
            self.logger.warning("The start time of the last cado-nfs.py "
                                "run was recorded, but not its end time, "
                                "maybe because it died unexpectedly.")
            self.logger.warning("Elapsed time of last run is not known and "
                                "will not be counted towards total.")
        self.state["starttime"] = time.time()

    def end_elapsed_time(self):
        if "starttime" not in self.state:
            self.logger.error("Missing starttime in end_elapsed_time(). "
                              "This should not have happened.")
            return
        elapsed = time.time() - self.state["starttime"]
        elapsed += self.state.get("elapsed", 0)
        self.state.__delitem__("starttime", commit=False)
        self.state.update({"elapsed": elapsed}, commit=True)
        return elapsed

    def next_task(self):
        for task in self.tasks:
            if task in self.tasks_that_want_to_run:
                # self.logger.info("Next task that wants to run: %s",
                #                  task.title)
                self.tasks_that_want_to_run.remove(task)
                return task
        return None

    def get_sum_of_cpu_or_real_time(self, is_cpu):
        total = 0
        for task in self.tasks:
            task_time = task.get_total_cpu_or_real_time(is_cpu)
            total += task_time
            # self.logger.info("Task %s reports %s time of %g,"
            #                  " new total: %g",
            #                  task.name,
            #                  "cpu" if is_cpu else "real", task_time, total)
        return total

    def register_filename(self, d):
        return self.servertask.register_filename(d)

    def relay_notification(self, message):
        """
        The relay for letting Tasks talk to us and each other
        """
        assert isinstance(message, Notification)
        sender = message.get_sender()
        key = message.get_key()
        value = message.get_value()
        self.logger.message("Received notification from %s,"
                            " key = %s, value = %s",
                            sender, Notification.reverse_lookup(key), value)
        if key is Notification.WANT_MORE_RELATIONS:
            if sender is self.purge:
                self.dup2.request_more_relations(value)
                if hasattr(self, "filtergalois"):
                    self.filtergalois.request_more_relations(value)
            elif sender is self.dup2:
                self.dup1.request_more_relations(value)
            elif sender is self.dup1:
                self.sieving.request_more_relations(value)
            else:
                raise Exception("Got WANT_MORE_RELATIONS"
                                " from unknown sender")
        elif key is Notification.HAVE_ENOUGH_RELATIONS:
            if sender is self.purge:
                self.servertask.stop_serving_wus()
                self.sieving.cancel_available_wus()
                self.stop_all_clients()
            else:
                raise Exception("Got HAVE_ENOUGH_RELATIONS"
                                " from unknown sender")
        elif key is Notification.REGISTER_FILENAME:
            if isinstance(sender, ClientServerTask):
                self.register_filename(value)
            else:
                raise Exception("Got REGISTER_FILENAME,"
                                " but not from a ClientServerTask")
        elif key is Notification.WANT_TO_RUN:
            if sender in self.tasks_that_want_to_run:
                raise Exception("Got request from %s to run,"
                                " but it was in run queue already",
                                sender)
            else:
                self.tasks_that_want_to_run.append(sender)
        elif key is Notification.SUBSCRIBE_WU_NOTIFICATIONS:
            return self.subscribeObserver(sender)
        else:
            raise KeyError("Notification from %s has unknown key %s"
                           % (sender, key))

    def answer_request(self, request):
        assert isinstance(request, Request)
        sender = request.get_sender()
        key = request.get_key()
        value = request.get_value()
        self.logger.message("Received request from %s, key = %s, values = %s",
                            sender, Request.reverse_lookup(key), value)
        if key not in self.request_map:
            raise KeyError("Unknown Request key %s from sender %s" %
                           (key, sender))
        if value is None:
            result = self.request_map[key]()
        else:
            result = self.request_map[key](value)
        self.logger.message("Completed request from"
                            " %s, key = %s, values = %s,"
                            " result = %s",
                            sender, Request.reverse_lookup(key), value,
                            result)
        return result

    # Colleagues now call either send_notification or send_request.
    # So handle_message is no more.
    #
    # def handle_message(self, message):
    #     if isinstance(message, Notification):
    #         self.relay_notification(Notification)
    #     elif isinstance(message, Request):
    #         return self.answer_request(message)
    #     else:
    #         raise TypeError("Message is neither Notification nor Request")
