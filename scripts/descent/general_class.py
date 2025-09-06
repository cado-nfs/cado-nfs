import os
import re
import shutil
import subprocess
import sqlite3
import tempfile

from .log_base import LogBase


# This gives the boring info about the file names for everything, and
# the boilerplate arguments to be fed to binaries.
class GeneralClass(object):

    def declare_args(parser):
        # parser.add_argument("--no-wipe",
        #         help="Keep working files",
        #         action="store_true")
        parser.add_argument("--sm-mode",
                            help="Select SM mode",
                            type=str)
        parser.add_argument("--datadir",
                            help="cadofactor working directory",
                            type=str)
        parser.add_argument("--prefix",
                            help="project prefix",
                            type=str)
        parser.add_argument("--db",
                            help="SQLite db file name",
                            type=str)
        # the next few are optional file names
        parser.add_argument("--tmpdir", help="Temporary working directory")
        parser.add_argument("--cadobindir",
                            help="Cado build directory",
                            required=True,
                            type=str)
#        parser.add_argument("--todofile",
#                help="Output primes to this toto file",
#                type=str)
        parser.add_argument("--poly",
                            help="Polynomial file",
                            type=str)
        parser.add_argument("--renumber",
                            help="Renumber file",
                            type=str)
        parser.add_argument("--fb1",
                            help="Factor base file for the algebraic side",
                            type=str)
        parser.add_argument("--log",
                            help="File with known logs",
                            type=str)
        parser.add_argument("--no-logs",
                            help="Work without logs, only compute relations",
                            default=False,
                            action='store_true')
        parser.add_argument("--gfpext",
                            help="Degree of extension (default 1)",
                            type=int)
        # This one applies to both las in the initial step
        parser.add_argument("--threads",
                            help="Number of threads to use",
                            type=int, default=4)
        # the arguments below are really better fetched from the
        # database.
        parser.add_argument("--ell", help="Group order (a.k.a. ell)")
        parser.add_argument("--nsm0",
                            help="Number of Schirokauer maps on side 0")
        parser.add_argument("--nsm1",
                            help="Number of Schirokauer maps on side 1")
        # Those are used both for the middle and lower levels of the
        # descent.
        parser.add_argument("--lpb0",
                            help="Default large prime bound on side 0",
                            required=True,
                            type=int)
        parser.add_argument("--lpb1",
                            help="Default large prime bound on side 1",
                            required=True,
                            type=int)
        parser.add_argument("--memory-margin",
                            help="Keep this amount of RAM free for"
                            " the rest of the world (see -t help)"
                            " (in gigabytes, floating point values allowed)",
                            type=int)
        parser.add_argument("--descent-max-increase-A", type=int, default=4)
        parser.add_argument("--descent-max-increase-lpb", type=int, default=0)

    def __init__(self, args):
        self._conn = None
        self.args = args
        if bool(self.args.db) == bool(self.args.prefix and self.args.datadir):
            raise ValueError(
                "Either --db (with an sqlite db) "
                "or the combo --prefix + --datadir must be specified")
        if self.args.tmpdir:
            self._tmpdir = self.args.tmpdir
            # do mkdir ???
        else:
            self._tmpdir = tempfile.mkdtemp(prefix="cado-nfs.")
        self.hello()
        self.__load_badidealdata()
        if self.args.no_logs:
            self.logDB = None
        else:
            self.logDB = LogBase(self)
        self.initrandomizer = 1

    def __connect(self):
        if self.args.db and not self._conn:
            self._conn = sqlite3.connect(self.args.db)

    def __getdb(self, query):
        if not self.args.db:
            return None
        self.__connect()
        self._cursor = self._conn.cursor()
        self._cursor.execute(query)
        v = self._cursor.fetchone()
        self._cursor.close()
        return v

    def __getfile(self, shortname, typical, table, key):
        try:
            v = self.args.__dict__[shortname]
            if v:
                return v
        except KeyError:
            pass
        if self.args.db and table:
            v = self.__getdb(f"select value from {table} where kkey='{key}'")
            if v is not None and len(v) > 0:
                return os.path.join(os.path.dirname(self.args.db), v[0])
        elif self.args.datadir and self.args.prefix:
            return os.path.join(self.args.datadir,
                                self.args.prefix + "." + typical)
        raise ValueError("no %s file known" % shortname)

    def __getarg(self, shortname, table, key):
        try:
            v = self.args.__dict__[shortname]
            if v:
                return v
        except KeyError:
            pass
        if self.args.db:
            v = self.__getdb(f"select value from {table} where kkey='{key}'")
            if v is not None and len(v) > 0:
                return v[0]
        raise ValueError("no %s parameter known" % shortname)

    def prefix(self):
        if self.args.prefix:
            return self.args.prefix
        else:
            return os.path.basename(self.args.db).split('.')[0]

    def datadir(self):
        if self.args.datadir:
            return self.args.datadir
        elif self.args.db:
            return os.path.dirname(self.args.db)
        else:
            raise ValueError("Need --datadir or --db with an sqlite db")

    def poly(self):
        return self.__getfile("poly", "poly", "polyselect2", "polyfilename")

    def renumber(self):
        return self.__getfile("renumber",
                              "renumber.gz",
                              "freerel",
                              "renumberfilename")

    def log(self):
        if self.args.no_logs:
            return None
        else:
            return self.__getfile("log", "dlog", "reconstructlog", "dlog")

    def badideals(self):
        return os.path.join(self.datadir(), self.prefix() + ".badideals")

    def badidealinfo(self):
        return os.path.join(self.datadir(), self.prefix() + ".badidealinfo")

    def fb1(self):
        return self.__getfile("fb1", "roots1.gz", "factorbase", "outputfile")

    def fb0(self):
        return self.__getfile("fb0", "roots0.gz", "factorbase", "outputfile")

    def ell(self):
        return int(self.args.ell)

    def lpb0(self):
        return self.args.lpb0

    def lpb1(self):
        return self.args.lpb1

    def tmpdir(self):
        return self._tmpdir

    def threads(self):
        return int(self.args.threads)

    def poly_data(self):
        d = {}
        with open(self.poly(), "r") as file:
            for line in file:
                if re.match(r"^\s*#", line):
                    continue
                if re.match(r"^\s*$", line):
                    continue
                key, value = line.split(":")
                key = key.strip()
                foo = re.match(r"^([cY])(\d+)$", key)
                if foo:
                    s, i = foo.groups()
                    if s not in d:
                        d[s] = []
                    while int(i) >= len(d[s]):
                        d[s] += [None]
                    d[s][int(i)] = value.strip()
                else:
                    d[key] = value.strip()
            # after all, we don't need to parse the polynomials at all.
            # only has_rational_side is ever used, and all we care about
            # is its boolean value.

            # If we _must_ parse the polynomials, it's a bit annoying now
            # that they can be written in algebraic form...
            if 'poly0' in d:
                assert 'Y' not in d
                v = d["poly0"]
                if 'x' in v:
                    raise ValueError("How do I parse a polynomial?")
                    # d['Y'] = ZZ['x'](v).list()
                else:
                    d['Y'] = [int(x) for x in v.split(',')]
            if 'poly1' in d:
                assert 'c' not in d
                v = d["poly1"]
                if 'x' in v:
                    raise ValueError("How do I parse a polynomial?")
                    # d['c'] = ZZ['x'](v).list()
                else:
                    d['c'] = [int(x) for x in v.split(',')]
        return d

    def p(self):
        d = self.poly_data()
        return int(d["n"])

    def extdeg(self):
        if self.args.gfpext:
            return self.args.gfpext
        else:
            return 1

    def target(self):
        if self.extdeg() == 1:
            return int(self.args.target)
        else:
            return [int(x) for x in self.args.target.split(",")]

    # short name for the target, to be used in filenames
    def short_target(self):
        target = str(self.args.target)
        if len(target) <= 20:
            return target
        else:
            return target[:10] + "..." + target[-10:]

    def has_rational_side(self):
        d = self.poly_data()
        return len(d["Y"]) == 2

    def cleanup(self):
        # if not self.args.tmpdir and not self.args.no_wipe:
        if False:
            shutil.rmtree(self.tmpdir())

    def __del__(self):
        if self._conn:
            self._conn.close()

    def descentinit_bin(self):
        return os.path.join(self.args.cadobindir, "misc", "descent_init_Fp")

    def las_bin(self):
        return os.path.join(self.args.cadobindir, "sieve", "las")

    def sm_simple_bin(self):
        return os.path.join(self.args.cadobindir, "filter", "sm_simple")

    def numbertheory_bin(self):
        return os.path.join(self.args.cadobindir, "utils", "numbertheory_tool")

    def lasMiddle_base_args(self):
        # TODO add threads once it's fixed.
        s = [self.las_bin() + "_descent",
             "--recursive-descent",
             "--allow-largesq",
             "--never-discard",  # useful for small computations.
             "--adjust-strategy", 2,  # avoids shrinking.
             "--fb1", self.fb1(),
             "--poly", self.poly(),
             "--descent-max-increase-A", self.args.descent_max_increase_A,
             "--descent-max-increase-lpb", self.args.descent_max_increase_lpb,
             ]
        if not self.args.no_logs:
            s += ["--renumber", self.renumber()]
            s += ["--log", self.log()]
        if not self.has_rational_side():
            s += ["--fb0", self.fb0()]
        return [str(x) for x in s]

    # There's no las_init_base_args, since DescentUpperClass uses only
    # its very own arguments.

    def hello(self):
        print("Working in GF(p), p=%d" % self.p())
        print("Subgroup considered in GF(p)^* has size %d" % self.ell())
        print("prefix is %s" % self.prefix())
        errors = []
        if not os.path.exists(self.las_bin()):
            errors.append("las not found (make las ?)")
        if not os.path.exists(self.las_bin() + "_descent"):
            errors.append("las_descent not found"
                          " (make las_descent ?)")
        if not os.path.exists(self.sm_simple_bin()):
            errors.append("sm_simple not found"
                          " (make sm_simple ?)")
        if not os.path.exists(self.numbertheory_bin()):
            errors.append("numbertheory_tool not found"
                          " (make numbertheory_tool ?)")

        for f in [self.log(),
                  self.poly(),
                  self.renumber(),
                  self.fb1()
                  ]:
            if f is not None and not os.path.exists(f):
                errors.append("%s missing" % f)

        if len(errors):
            msg = "Some data files and/or binaries missing:\n"
            msg += "\n".join(["\t" + x for x in errors])
            raise RuntimeError(msg)

    # self.list_badideals will contain a list of (p,r,side)
    # self.list_bad_ncols will contain a list of the corresponding nb of cols
    # self.badidealdata will contain a list of
    #     (p, k, rk, side, [exp1, exp2, ..., expi])
    def __load_badidealdata(self):
        # Note that we can as well get this information from the renumber
        # file.
        if not os.path.exists(self.badideals()) \
                or not os.path.exists(self.badidealinfo()):
            call_that = [self.numbertheory_bin(),
                         "-poly", self.poly(),
                         "-badideals", self.badideals(),
                         "-badidealinfo", self.badidealinfo(),
                         "-ell", self.ell()
                         ]
            call_that = [str(x) for x in call_that]
            print("command line:\n" + " ".join(call_that))
            with open(os.devnull, 'w') as devnull:
                subprocess.check_call(call_that, stderr=devnull)

        self.list_badideals = []
        self.list_ncols = []
        with open(self.badideals(), 'r') as bad:
            for line in bad:
                if line[0] == '#':
                    continue
                foo = re.match(r"^(\d+),(\d+):(\d+): (\d+)$", line)
                if foo:
                    self.list_badideals.append((int(foo.groups()[0]),
                                                int(foo.groups()[1]),
                                                int(foo.groups()[2])))
                    self.list_ncols.append(int(foo.groups()[3]))
                else:
                    raise ValueError(f"Error while reading {self.badideal()}")

        self.badidealdata = []
        with open(self.badidealinfo(), 'r') as bad:
            for line in bad:
                if line[0] == '#':
                    continue
                pattern = r"^(\d+) (\d+) (\d+) (\d+) (.+)$"
                foo = re.match(pattern, line)
                if foo:
                    self.badidealdata.append(
                        (int(foo.groups()[0]),  # p
                         int(foo.groups()[1]),  # k
                         int(foo.groups()[2]),  # rk
                         int(foo.groups()[3]),  # side
                         [int(x) for x in foo.groups()[4].split()]  # exp
                         ))
                else:
                    raise ValueError("Error while reading %s" %
                                     self.badidealinfo())

        print("Bad ideal information loaded: %s bad ideals,"
              " and %s lines in badidealinfo"
              % (str(len(self.list_badideals)),
                 str(len(self.badidealdata))))
        print("badideal data: %s" % str(self.badidealdata))
