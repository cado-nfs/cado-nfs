import os
import re
import sys
import math
import random
import subprocess
import shutil
import functools

from .ideals_above_p import ideals_above_p
from .descent_utils import feature_get_hwloc
from .descent_utils import object_holder
from .descent_utils import a_over_b_mod_p

if sys.version_info >= (3, 8):
    from .descent_helper_asyncio import monitor_important_files
else:
    from .descent_helper_fallback import monitor_important_files


class DescentUpperClass(object):
    def declare_args(parser):
        c = " (specific for the descent bootstrap)"
        parser.add_argument("--init-tkewness",
                            help="Tkewness" + c,
                            type=int,
                            default=2**30)
        parser.add_argument("--init-lim",
                            help="Factor base bound" + c,
                            default=2**26)
        parser.add_argument("--init-lpb",
                            help="Large prime bound" + c,
                            default=64)
        parser.add_argument("--init-mfb",
                            help="Cofactor bound" + c,
                            default=100)
        parser.add_argument("--init-ncurves",
                            help="ECM effort in cofactorization" + c,
                            default=80)
        parser.add_argument("--init-I",
                            help="Sieving range" + c,
                            default=14)
        parser.add_argument("--init-minB1",
                            help="ECM first B1" + c,
                            default=200)
        parser.add_argument("--init-mineff",
                            help="ECM minimal effort" + c,
                            default=1000)
        parser.add_argument("--init-maxeff",
                            help="ECM maximal effort" + c,
                            default=100000)
        parser.add_argument("--init-side",
                            help="Side of the bootstrap"
                            " (when there is no rational side)",
                            default=1)
        # Slave las processes in the initial step.
        parser.add_argument("--slaves",
                            help="Number of slaves to use",
                            type=int, default=1)
        # In case we used an external process
        parser.add_argument("--external-init",
                            help="Use precomputed external data"
                            " for the descent bootstrap",
                            type=str,
                            default=None)

    def has_log(self, *args):
        if self.general.logDB is None:
            p, r, side = args
            return p < self.lim
        return self.general.logDB.has(*args)

    def __init__(self, general, args):
        self.general = general
        self.args = args

        if self.args.external_init is not None:
            self.external = self.args.external_init
            if not os.path.exists(self.external):
                raise NameError("Given external file for init does not exist")
        else:
            self.external = None
            for a in ["tkewness", "lim", "lpb", "mfb",
                      "ncurves", "I", "side",
                      "mineff", "maxeff", "minB1"]:
                setattr(self, a, int(getattr(self.args, "init_" + a)))
            self.slaves = int(self.args.slaves)
            # the final step needs to know the init side as well.
            general.init_side = int(self.args.init_side)

    def __isqrt(self, n):
        x = n
        y = (x + 1) // 2
        while y < x:
            x = y
            y = (x + n // x) // 2
        return x

    def __myxgcd(self, a, b, T):
        assert type(a) is int
        assert type(b) is int
        assert type(T) is int
        # ainit = a
        # binit = b
        bound = self.__isqrt(b * T)
        x = 0
        lastx = 1
        y = 1
        lasty = 0
        while abs(b) > bound:
            q = a // b
            r = a % b
            a = b
            b = r
            newx = lastx - q * x
            lastx = x
            x = newx
            newy = lasty - q * y
            lasty = y
            y = newy
        return [[b, x], [a, lastx]]

    def use_external_data(self, z):
        general = self.general
        fil = open(self.external, "r")
        rrr = fil.read()
        fil.close()
        lines = rrr.splitlines()
        e = int(lines[0])
        Num = int(lines[1])
        Den = int(lines[2])
        # check that we are talking about the same z!
        p = self.general.p()
        zz = pow(z, e, p)
        assert (zz * Den - Num) % p == 0
        general.initrandomizer = e       # for later use
        fnum = [int(x) for x in lines[3].split()]
        fden = [int(x) for x in lines[4].split()]
        large_q = [int(x) for x in lines[5].split()]
        descrelfile = lines[6]

        # create todolist from fnum and fden, skipping primes of the
        # large_q list
        prefix = f"{general.prefix()}.descent.{general.short_target()}.init."
        todofilename = os.path.join(general.datadir(), prefix + "todo")
        with open(todofilename, "w") as f:
            for q in fnum + fden:
                if q in large_q:
                    continue
                if self.has_log(q, -1, 0):
                    continue
                logq = math.ceil(math.log(q, 2))
                print("Will do further descent"
                      " for %d-bit rational prime %d"
                      % (logq, q))
                # las can understand when the rational root is missing
                f.write("0 %d\n" % q)
        fil = open(descrelfile, "r")
        rrr = fil.read()
        fil.close()
        lines = rrr.splitlines()
        with open(todofilename, "a") as f:
            for line in lines:
                foo = re.match(r"^Taken: ([0-9\-]+),([0-9\-]+):"
                               r"([0-9a-fA-F,]+):([0-9a-fA-F,]+)",
                               line)
                assert foo
                foog = foo.groups()
                a = int(foog[0])
                b = int(foog[1])
                list_p = [[int(x, 16) for x in foog[i].split(",")]
                          for i in [2, 3]]

                for side in range(2):
                    for p in list_p[side]:
                        if p in large_q:
                            continue
                        if side == 0:
                            if not self.has_log(p, -1, 0):
                                f.write("0 %d\n" % p)
                        else:
                            if b % p == 0:
                                continue
                            ideal = ideals_above_p(p, 1, a, b, side, general)
                            if ideal.get_log() is not None:
                                continue
                            else:
                                r = a_over_b_mod_p(a, b, p)
                                f.write("1 %d %d\n" % (p, r))
        return todofilename, [Num, Den, fnum, fden], descrelfile

    def do_descent_for_real(self, z, seed,
                            randomize_multiplicatively=None):
        """
        This is a convenient entry point, limited to the first step of
        the descent. The goal is to find an initial split for the target
        z, i.e. an expression z=u/v with reasonaly smooth u and v. The
        "T-kewness" approach is used.

        This step typically requires randomization. We have two types of
        randomization.
         - by default, we do "exponential randomization", which is
           accessed by `randomize_multiplicatively=False`. Here, a random
           exponent h is chosen, and we look for a split of `z^h`. Of
           course multiple such exponents h are tried, and this depends
           on the seed. This is only doable if we know how to invert h
           modulo the group order, and this is adapted to discrete
           logarithms. The integer h is returned as a fourth return
           value.
         - The "multiplicative randomization", for which
           `randomize_multiplicatively` is set to a non-zero integer e. In
           this mode, a random number is picked, and we look for a split
           of `m^e*z`. This is adapted to the computation of e-th roots.
           The integer m is returned as a fourth return value. Note that
           this mode also *implies* that we call this function directly,
           and it is not compatible with what is done in
           descent_lower_class.
        """
        general = self.general
        p = general.p()
        bound = p.bit_length() // 2 + 20
        # make the randomness deterministic to be able to replay
        # interrupted computations.
        random.seed(seed)
        while True:
            if randomize_multiplicatively is None:
                h = random.randrange(p)
                zz = pow(z, h, p)
                general.initrandomizer = h
            else:
                e = randomize_multiplicatively
                m = random.randrange(p)
                zz = (pow(m, e, p) * z) % p
                general.initrandomizer = m

            gg = self.__myxgcd(zz, p, self.tkewness)
            if (gg[0][0].bit_length() < bound and
                    gg[1][0].bit_length() < bound and
                    gg[0][1].bit_length() < bound and
                    gg[1][1].bit_length() < bound):
                # we're happy
                break
            print("Skewed reconstruction. Let's randomize the input.")

        tmpdir = general.tmpdir()
        prefix = f"{general.prefix()}.descent" \
                 f".{general.short_target()}" \
                 f".{seed}.init."

        polyfilename = os.path.join(tmpdir, prefix + "poly")
        with open(polyfilename, 'w') as f:
            f.write("n: %d\n" % p)
            f.write("skew: 1\n")
            f.write("c1: %d\n" % gg[0][0])
            f.write("c0: %d\n" % gg[1][0])
            f.write("Y1: %d\n" % gg[0][1])
            f.write("Y0: %d\n" % gg[1][1])

        print("--- Sieving (initial) ---")

        def relation_filter(data):
            line = data.strip()
            if re.match(r"^[^#]-?\d+,\d+:(\w+(,\w+)*)?:(\w+(,\w+)*)?$", line):
                return line

        relsfilename = os.path.join(general.datadir(), prefix + "rels")

        if os.path.exists(relsfilename):
            sources = [(relsfilename, [])]
        else:
            fbcfilename = os.path.join(tmpdir, prefix + "fbc")
            call_common = [general.las_bin(),
                           "-poly", polyfilename,
                           "-lim0", self.lim,
                           "-lim1", self.lim,
                           "-lpb0", self.lpb,
                           "-lpb1", self.lpb,
                           "-mfb0", self.mfb,
                           "-mfb1", self.mfb,
                           "-ncurves0", self.ncurves,
                           "-ncurves1", self.ncurves,
                           "-fbc", fbcfilename,
                           "-I", self.I
                           ]

            def fbc_call():
                call_that = call_common + [
                    "-q0", self.tkewness,
                    "-q1", self.tkewness,
                    "-nq", 0]
                if os.environ.get("CADO_NFS_MAX_THREADS") is not None:
                    call_that += ["-t", os.environ["CADO_NFS_MAX_THREADS"]]
                elif feature_get_hwloc():
                    call_that += ["-t", "machine,1,pu"]
                else:
                    call_that += ["-t", 4]
                call_that = [str(x) for x in call_that]
                return call_that

            def construct_call(q0, q1):
                call_that = call_common + [
                    "-q0", q0,
                    "-q1", q1,
                    "--exit-early", 2]
                if os.environ.get("CADO_NFS_MAX_THREADS") is not None:
                    call_that += ["-t", os.environ["CADO_NFS_MAX_THREADS"]]
                elif feature_get_hwloc():
                    call_that += ["-t", "auto"]
                else:
                    call_that += ["-t", 4]
                call_that = [str(x) for x in call_that]
                return call_that

            # outfile
            call_params = [(os.path.join(relsfilename + "." + str(i)),
                            self.tkewness + 100000 * i,  # q0
                            self.tkewness + 100000 * (i + 1))
                           for i in range(self.slaves)]  # q1

            if not os.path.exists(fbcfilename):
                all_ok = True
                for t in call_params:
                    if not os.path.exists(t[0]):
                        all_ok = False
                        break

                if all_ok:
                    print(" - Using %s" % fbcfilename)
                else:
                    print(" - Factor base cache -")
                    with open(os.devnull, 'w') as devnull:
                        subprocess.check_call(fbc_call(), stdout=devnull)
                    print(" - done -")

            # Whether or not the output files are already present, this
            # will do the right thing and run the new processes only if
            # needed.
            sources = [(outfile, construct_call(q0, q1))
                       for (outfile, q0, q1) in call_params]

        rel_holder = object_holder()

        def consume(rel_holder, idx, line):
            if line[0] != '#':
                rel_holder.set((idx, line))
                return True
            if re.match(r"^# (Now sieving.*q=|\d+ relation)", line):
                sys.stdout.write('\n')
                print(line.rstrip())
                sys.stdout.flush()

        monitor_important_files(sources,
                                consume,
                                (rel_holder,),
                                temporary_is_reusable=True)

        sys.stdout.write('\n')
        if rel_holder.v is None:
            print("No relation found!")
            print("Trying again with another random seed...")
            return None, None, None, None

        idx, rel = rel_holder.v

        if not os.path.exists(relsfilename):
            shutil.copyfile(sources[idx][0], relsfilename)

        print("Taking relation %s\n" % rel)
        rel = rel.split(':')
        a, b = [int(x) for x in rel[0].split(',')]

        Num = a * gg[0][0] + b * gg[1][0]
        Den = a * gg[0][1] + b * gg[1][1]
        assert (zz * Den - Num) % p == 0

        factNum = [int(x, 16) for x in rel[2].split(',')]
        factDen = [int(x, 16) for x in rel[1].split(',')]
        print(Num, Den, factNum, factDen)

        assert abs(Num) == functools.reduce(lambda x, y: x * y, factNum, 1)
        assert abs(Den) == functools.reduce(lambda x, y: x * y, factDen, 1)

        lc_ratpol = int(general.poly_data()["Y"][1])
        for q in factNum + factDen:
            if not self.has_log(q, -1, 0):
                if lc_ratpol % q == 0:
                    print("Would need to descend %s" % q,
                          "which divides the lc of the rational poly.")
                    print("Trying again with a new seed.")
                    return None, None, None, None

        todofilename = os.path.join(general.datadir(), prefix + "todo")

        if not os.path.exists(todofilename):
            with open(todofilename, "w") as f:
                for q in factNum + factDen:
                    if self.has_log(q, -1, 0):
                        continue
                    logq = math.ceil(math.log(q, 2))
                    print("Will do further descent",
                          "for %d-bit rational prime %d" % (logq, q))
                    # las can understand when the rational root is missing
                    f.write("0 %d\n" % q)
        else:
            with open(todofilename, "r") as f:
                for line in f:
                    side, q = line.strip().split(' ')
                    q = int(q)
                    if self.has_log(q, -1, 0):
                        continue
                    logq = math.ceil(math.log(q, 2))
                    print("Will do further descent",
                          "for %d-bit rational prime %d" % (logq, q))

        return (todofilename,
                [Num, Den, factNum, factDen],
                None,
                general.initrandomizer)

    def do_descent_nonlinear(self, z):
        general = self.general
        p = general.p()
        # tmpdir = general.tmpdir()
        prefix = f"{general.prefix()}.descent.{general.short_target()}.upper."
        # polyfilename = os.path.join(tmpdir, prefix + "poly")
        if general.extdeg() == 1:
            zz = [z]
        else:
            zz = z
        call_that = [general.descentinit_bin(),
                     "-poly", general.poly(),
                     "-mt", 4,
                     "-minB1", self.minB1,
                     "-mineff", self.mineff,
                     "-maxeff", self.maxeff,
                     "-side", self.side,
                     "-extdeg", general.extdeg(),
                     "-lpb", self.lpb,
                     "-seed", 42,
                     "-jl",
                     p
                     ] + zz
        call_that = [str(x) for x in call_that]
        initfilename = os.path.join(general.datadir(), prefix + "init")

        has_winner = object_holder(True)

        def consume(has_winner, general, idx, line):
            foo = re.match(r"^Youpi: e = (\d+) is a winner", line)
            if foo:
                has_winner.set(True)
                general.initrandomizer = int(foo.groups()[0])
            foo = re.match(r"^U = ([0-9\-,]+)", line)
            if foo:
                general.initU = [int(x) for x in foo.groups()[0].split(',')]
            foo = re.match(r"^V = ([0-9\-,]+)", line)
            if foo:
                general.initV = [int(x) for x in foo.groups()[0].split(',')]
            foo = re.match(r"^u = ([0-9]+)", line)
            if foo:
                general.initu = int(foo.groups()[0])
            foo = re.match(r"^v = ([0-9]+)", line)
            if foo:
                general.initv = int(foo.groups()[0])
            foo = re.match(r"^fac_u = ([, 0-9]+)", line)
            if foo:
                general.initfacu = [[int(y) for y in x.split(',')]
                                    for x in foo.groups()[0].split(' ')]
            foo = re.match(r"^fac_v = ([, 0-9]+)", line)
            if foo:
                general.initfacv = [[int(y) for y in x.split(',')]
                                    for x in foo.groups()[0].split(' ')]

        monitor_important_files([(initfilename, call_that)],
                                consume,
                                (has_winner, general))

        if not has_winner.v:
            raise ValueError("initial descent failed for target %s" % zz)

        todofilename = os.path.join(general.datadir(), prefix + "todo")
        print(general.initfacu)
        print(general.initfacv)

        if not os.path.exists(todofilename):
            with open(todofilename, "w") as f:
                for ideal in general.initfacu + general.initfacv:
                    q = ideal[0]
                    r = ideal[1]
                    if self.has_log(q, r, self.side):
                        continue
                    logq = math.ceil(math.log(q, 2))
                    print("Will do further descent",
                          "for %d-bit prime %d" % (logq, q))
                    f.write("%d %d %d\n" % (self.side, q, r))
        else:
            with open(todofilename, "r") as f:
                for line in f:
                    ll = line.strip().split(' ')
                    # side = ll[0]
                    q = ll[1]
                    q = int(q)
                    logq = math.ceil(math.log(q, 2))
                    print("Will do further descent",
                          "for %d-bit prime %d" % (logq, q))

        return todofilename, [general.initU,
                              general.initV,
                              general.initfacu,
                              general.initfacv], None

    def do_descent(self, z):
        if self.external:
            return self.use_external_data(z)

        if not self.general.has_rational_side():
            return self.do_descent_nonlinear(z)

        seed = 42

        while True:
            tdf, spl, frf, *tail = self.do_descent_for_real(z, seed)
            if tdf is not None:
                return tdf, spl, frf
            else:
                seed += 1
