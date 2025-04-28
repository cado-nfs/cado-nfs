import os
import re
import sys

from .descent_utils import feature_get_hwloc
from .descent_utils import object_holder
from .descent_utils import FailedDescent

if sys.version_info >= (3, 8):
    from .descent_helper_asyncio import monitor_important_files
else:
    from .descent_helper_fallback import monitor_important_files


class DescentMiddleClass(object):
    def declare_args(parser):
        # TODO: import default values from the sieving parameters.
        parser.add_argument("--descent-hint",
                            help="Hintfile for the descent",
                            required=True,
                            # TODO: fall back on default values.
                            )
        parser.add_argument("--I",
                            help="Default value for I (must match hint file)",
                            required=True,
                            type=int)
        parser.add_argument("--B",
                            help="set bucket region size in las calls",
                            required=False,
                            type=int)
        for side in range(2):
            parser.add_argument(f"--mfb{side}",
                                help=f"Default cofactor bound on side {side}",
                                required=True,
                                type=int)
            parser.add_argument(f"--lim{side}",
                                help="Default factor base bound"
                                     f" on side {side} (must match hint file)",
                                required=True,
                                type=int)

    def __init__(self, general, args):
        self.general = general
        self.args = args
        # We need to do some safety checking
        values_I = set()
        values_lim0 = set()
        values_lim1 = set()
        values_I.add(self.args.I)
        values_lim0.add(self.args.lim0)
        values_lim1.add(self.args.lim1)
        with open(self.args.descent_hint, 'r') as file:
            for line in file:
                if re.match(r"^\s*#", line):
                    continue
                if re.match(r"^\s*$", line):
                    continue
                line = line.strip()
                foo = re.match(r"^.*I=(\d+)\s+(\d+),[\d.,]+"
                               r"\s+(\d+),[\d.,]+$",
                               line)
                if not foo:
                    print("Warning, parse error in hint file",
                          "on line:\n" + line)
                    continue
                I, lim0, lim1 = foo.groups()
                values_I.add(int(I))
                values_lim0.add(int(lim0))
                values_lim1.add(int(lim1))
        if False:
            if len(values_lim0) > 1:
                raise ValueError("lim0 values should match"
                                 " between cmdline and hint file")
            if len(values_lim1) > 1:
                raise ValueError("lim1 values should match"
                                 " between cmdline and hint file")
            if len(values_I) > 1:
                raise ValueError("I values should match"
                                 " between cmdline and hint file")
        print("Consistency check for las_descent passed")
        print("\tI=%d" % values_I.pop())
        print("\tlim0=%d" % values_lim0.pop())
        print("\tlim1=%d" % values_lim1.pop())

    def do_descent(self, todofile):
        # tmpdir = self.general.tmpdir()
        prefix = f"{self.general.prefix()}.descent." \
                 f"{self.general.short_target()}.middle."

        f = open(todofile, 'r')
        ntodo = len(list(f))
        f.close()
        print("--- Sieving (middle, %d rational primes) ---" % ntodo)
        s = self.general.lasMiddle_base_args()
        if self.args.descent_hint:
            s += ["--descent-hint-table", self.args.descent_hint]
        s += ["--I", self.args.I,
              "--lim0", self.args.lim0,
              "--lim1", self.args.lim1,
              "--lpb0", self.general.lpb0(),
              "--mfb0", self.args.mfb0,
              "--lpb1", self.general.lpb1(),
              "--mfb1", self.args.mfb1,
              ]
        if os.environ.get("CADO_NFS_MAX_THREADS") is not None:
            s += ["-t", os.environ["CADO_NFS_MAX_THREADS"]]
        elif feature_get_hwloc():
            s += ["-t", "machine,1,pu"]
        else:
            s += ["-t", 4]
        if self.args.B is not None:
            s += ["--B", self.args.B]
        s += ["--todo", todofile]
        call_that = [str(x) for x in s]
        relsfilename = os.path.join(self.general.datadir(), prefix + "rels")

        printing = object_holder(False)
        failed = []

        def consume(printing, failed, idx, line):
            if re.match(r"^# taking path", line):
                print(line.rstrip())
            elif re.match(r"^# END TREE", line):
                print("")
                printing.unset()
            elif printing.v:
                print(line.rstrip())
                foo = re.match(r"# FAILED (\d+\@\d+)", line)
                if foo:
                    failed.append(foo.groups()[0])
            elif re.match(r"^# BEGIN TREE", line):
                print("")
                printing.set(True)

        monitor_important_files([(relsfilename, call_that)],
                                consume,
                                (printing, failed)
                                )

        if failed:
            raise FailedDescent(failed)

        return relsfilename
