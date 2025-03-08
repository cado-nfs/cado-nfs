import os
import subprocess
import sys
import argparse
import re

from .general_class import GeneralClass
from .descent_upper_class import DescentUpperClass
from .descent_middle_class import DescentMiddleClass
from .descent_lower_class import DescentLowerClass
from .descent_utils import drive_me_crazy, feature_set_hwloc, FailedDescent


def descent_main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        prog='descent',
        description="Descent initialization for DLP")

    # Required
    parser.add_argument("--target",
                        help="Element whose DL is wanted",
                        type=str,
                        required=True)
    parser.add_argument("--timestamp",
                        help="Prefix all lines with a time stamp",
                        action="store_true")

    GeneralClass.declare_args(parser)
    DescentUpperClass.declare_args(parser)
    DescentMiddleClass.declare_args(parser)
    DescentLowerClass.declare_args(parser)

    args = parser.parse_args()

    sys.stdout = drive_me_crazy(sys.stdout, args.timestamp)

    las_bin = os.path.join(args.cadobindir, "sieve", "las")
    cp = subprocess.Popen([ las_bin, "-help" ],
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)

    feature_set_hwloc(not bool(re.search(r"unused, needs hwloc",
                                         cp.stderr.read().decode())))

    general = GeneralClass(args)

    if general.target() == 1:
        # the re-randomization does not work for target=1
        print("# p=%d" % general.p())
        print("# ell=%d" % general.ell())
        print("# target=%s" % args.target)
        print("log(target)=0")
    else:
        init = DescentUpperClass(general, args)
        middle = DescentMiddleClass(general, args)
        lower = DescentLowerClass(general, args)

        todofile, initial_split, firstrelsfile = \
            init.do_descent(general.target())
        relsfile = middle.do_descent(todofile)
        if firstrelsfile:
            lower.do_descent([firstrelsfile, relsfile], initial_split)
        else:
            lower.do_descent([relsfile], initial_split)

    general.cleanup()
