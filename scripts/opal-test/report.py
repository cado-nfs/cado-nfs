#!/usr/bin/env python3

import sys
import re
from math import log, sqrt
from collections import Iterable

RE_FP = r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?"
CAP_FP = "(%s)" % RE_FP
REGEXES = {"cap_fp" : CAP_FP}

PATTERN_SIEVE_SQ = re.compile(r"# Now sieving (algebraic|rational|side-0|side-1) q=(\d+);")
PATTERN_SQ = re.compile(r"# Average (.*) for (\d+) special-q's, max bucket fill (.*)")
PATTERN_CPUTIME = re.compile(r"# Total cpu time {cap_fp}s, useful {cap_fp}s \[norm {cap_fp}\+{cap_fp}, sieving {cap_fp} \({cap_fp}\+{cap_fp} \+ {cap_fp}\), factor {cap_fp} \({cap_fp}\+{cap_fp} \+ {cap_fp}\), rest {cap_fp}\], wasted\+waited {cap_fp}s, rest {cap_fp}s".format(**REGEXES))
PATTERN_REPORTS = re.compile(r"# Total (\d+) reports".format(**REGEXES))
PATTERN_DUPE = re.compile("# DUPE ")

class NumInt(object):
    """
    >>> n = NumInt()
    >>> n.add(0,1)
    >>> n.add(1,2)
    >>> n.get_value()
    1.5
    >>> n.add(2,3)
    >>> n.get_value()
    4.0
    >>> x = n.interpolate_for_value(3.); 1.6457 < x < 1.6458
    True
    >>> 2.9999 < n.interpolate_at_coord(x) < 3.0001
    True
    """
    def __init__(self):
        # At index 0, the most recent value,coordinate pair.
        # At index 1, the value,coordinate pair before that.
        self.lastvalue = [None, None]
        self.lastcoord = [None, None]
        self.sum = 0
    def trapez_area(self):
        """ Return area of the trapezoidal defined by the last two pairs of coordinates """
        return (self.lastcoord[0] - self.lastcoord[1])*(self.lastvalue[0] + self.lastvalue[1]) / 2.
    def add(self, coord, value):
        self.lastvalue[1], self.lastcoord[1] = self.lastvalue[0], self.lastcoord[0]
        self.lastvalue[0], self.lastcoord[0] = value, coord
        if not self.lastcoord[1] is None:
            assert coord > self.lastcoord[1]
            self.sum += self.trapez_area()
    def get_value(self):
        return self.sum
    def interpolate_for_value(self, value):
        """ Find a coordinate c, greater than the second-last one, such that
        cutting or extending the last trapezoidal up to abscissa c results in
        sum = "value"
        """
        prev_sum = self.sum - self.trapez_area()
        diff = value - prev_sum
        t = (self.lastvalue[0] - self.lastvalue[1]) / (self.lastcoord[0] - self.lastcoord[1])
        if t <= 0: # can happen due to special-q correction
            # we estimate the function is constant to (v0+v1)/2
            v = (self.lastvalue[0] + self.lastvalue[1]) / 2
            return self.lastcoord[1] + diff / v
        # Choose offset x such that, with c = c1 + x and v = v1 + x * t,
        # prev_sum + (c - c1)*(v + v1) / 2 = value
        # thus (c - c1)*(v + v1) / 2 = diff
        # c - c1 = x, v + v1 = 2*v1 + x*t
        # x*(2*v1 + x*t) / 2 = diff
        # t/2*x^2 + v1*x - diff = 0
        # x = (-v1 +- sqrt(v1^2 + 2*t*diff)) / t
        # We need only the positive solution
        v1 = self.lastvalue[1]
        disc = v1**2 + 2*t*diff
        if disc < 0:
            sys.stderr.write("discriminant = %f < 0! t = %f, diff = %d, lv=%s, lc=%s, prev_sum = %f, value = %f\n" % (disc, t, diff, self.lastvalue, self.lastcoord, prev_sum, value))
        x = (-v1 + sqrt(disc)) / t
        return self.lastcoord[1] + x
    def interpolate_at_coord(self, coord):
        """ Return the sum that would result if the last trapezoidal had been
        cut or extended to abscissa "coord".
        """
        # due to the special-q correction, coord <= self.lastcoord[0] might
        # not hold, thus we disable the following assertion
        # assert self.lastcoord[1] <= coord <= self.lastcoord[0]
        x = coord - self.lastcoord[1]
        prev_sum = self.sum - self.trapez_area()
        return prev_sum + x * self.lastvalue[1]
        

class ListArith(list):
    """
    >>> a = ListArith([1,2,3])
    >>> b = ListArith([3,4,5])
    >>> a + 1
    [2, 3, 4]
    >>> a - 1
    [0, 1, 2]
    >>> a * 2
    [2, 4, 6]
    >>> a + b
    [4, 6, 8]
    >>> b - a
    [2, 2, 2]
    >>> a * b
    [3, 8, 15]
    """
    def __add__(self, other):
        if isinstance(other, Iterable):
            return ListArith([a + b for a,b in zip(self, other)])
        else:
            return ListArith([a + other for a in self])

    def __sub__(self, other):
        if isinstance(other, Iterable):
            return ListArith([a - b for a,b in zip(self, other)])
        else:
            return ListArith([a - other for a in self])

    def __mul__(self, other):
        if isinstance(other, Iterable):
            return ListArith([a * b for a,b in zip(self, other)])
        else:
            return ListArith([a * other for a in self])

    def to_str(self):
        formats = ["%s"] * len(self)
        pat = " ".join(formats)
        return pat % tuple(self)


class LasStats(object):
    def __init__(self):
        self.dupes = 0
        self.nr_sq = 0
        self.cputimes = ListArith([0.] * 8)
        self.reports = 0
        self.relations_int = NumInt()
        self.dupes_int = NumInt()
        self.elapsed_int = NumInt()

    def parse_one_input(self, lines, verbose=False):
        nr_sq = None
        cputimes = None
        reports = None
        new_dupes = 0
        first_sq = None
        for line in lines:
            if PATTERN_DUPE.match(line):
                new_dupes += 1
            match = PATTERN_SIEVE_SQ.match(line)
            if match:
                last_sq = int(match.group(2))
                if first_sq is None:
                    first_sq = last_sq
                else:
                    assert first_sq <= last_sq
            match = PATTERN_SQ.match(line)
            if match:
                nr_sq = int(match.group(2))
            match = PATTERN_CPUTIME.match(line)
            if match:
                cputimes = list(map(float, match.groups()))
            match = PATTERN_REPORTS.match(line)
            if match:
                reports = int(match.group(1))
        if cputimes is None:
            sys.stderr.write("Did not receive value for cputimes\n")
            return False
        if reports is None:
            sys.stderr.write("Did not receive value for reports\n")
            return False
        # check number of relations before number of special-q's
        if reports == 0:
            sys.stderr.write("No relation found in sample run, please increase q_range in las_run.py\n")
            return False
        if nr_sq is None:
            sys.stderr.write("Did not receive value for nr_sq\n")
            return False
        self.nr_sq += nr_sq
        self.dupes += new_dupes
        self.reports += reports
        self.cputimes += cputimes
        cputimes_str = self.cputimes.to_str()
        sq = (last_sq + first_sq) / 2
        sq_correction = 1./nr_sq/log(sq)
        self.relations_int.add(sq, reports * sq_correction)
        self.dupes_int.add(sq, new_dupes * sq_correction)
        self.elapsed_int.add(sq, cputimes[0] * sq_correction)
        if verbose:
            names = ("sq", "nr_sq", "sq_sum", "cputimes_str", "elapsed", "elapsed/sq", "elapsed/rel", "reports", "reports/nr_sq", "reports/sqrange", "dupes")
            values = (sq, nr_sq, self.nr_sq, cputimes_str, reports, reports/nr_sq, reports * sq_correction, self.dupes)
            print(", ".join( (":".join(map(str, x)) for x in zip(names, values)) ))
        return True

    def parse_one_file(self, filename, verbose=False):
        with open(filename, "r") as f:
            return self.parse_one_input(f, verbose)
     
    def get_rels(self):
        return self.relations_int.get_value()

    def get_dupes(self):
        return self.dupes_int.get_value()

    def get_qmax(self, nr_relations):
            return self.relations_int.interpolate_for_value(nr_relations)

    def get_time(self, nr_relations=None):
        if nr_relations is None:
            return self.elapsed_int.get_value()
        else:
            qmax = self.get_qmax(nr_relations)
            return self.elapsed_int.interpolate_at_coord(qmax)

    def print_stats(self):
        print("Estimated total relations: %f" % self.get_rels())
        print("Estimated total dupes: %f" % self.get_dupes())
        print("Estimated total elapsed time: %f" % self.get_time())
        print("Estimated relations/second: %f" % (self.get_rels() / self.get_time()))


def run():
    stats = LasStats()
    for filename in sys.argv[1:]:
        print("Parsing file %s" % filename)
        stats.parse_one_file(filename, verbose=True)
    stats.print_stats()

if __name__ == '__main__':
    run()
