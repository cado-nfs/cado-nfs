#!/usr/bin/env python3
import re
import argparse
import subprocess
from cadofactor import wudb


def do_check(acc, wu, programs):
    result = None
    if wu["files"] is None:
        return result
    for file in wu["files"]:
        (filename, filepath) = (file["filename"], file["path"])
        for (progname, pattern) in programs:
            m = re.match(pattern, filename)
            if m is None:
                continue
            assert m.start() == 0
            # Check that the whole string matched the pattern,
            # not just a prefix
            if not m.end() == len(filename):
                continue
            command = progname + " " + filepath
            print("Running: " + command)
            child = subprocess.Popen(command, shell=True,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
            (child_stdout, child_stderr) = child.communicate()
            print("Return code is " + str(child.returncode))
            if child.returncode == 0:
                result = True
            else:
                return False
    return result


if __name__ == "__main__":
    description = "Run a verification tool on files of received work units"
    epilog = """
    For each WU processed, run PROGRAM on those files of the WU whose
    filenames match REGEX. The client side filename (as specified in a
    RESULT line of a WU) is used for the REGEX match. Multiple programs
    can be specified; each program is run for any file(s) that match the
    respective REGEX. If the exit status of the PROGRAMs is zero for all
    tested files of a WU, the WU is considered verified ok. If the exit
    status is non-zero for any tested file, the WU is considered verified
    with error. The database is updated accordingly only if --update is
    given; otherwise only output to screen (and badidealsfile if
    specified) happens.
    """.strip().replace('    ', '').replace('\n', ' ')

    wudb_name = "wudb"
    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument('-f', '--file', action='store', help="Does nothing")
    parser.add_argument('-p', '--program', nargs=2, required=True,
                        action='append', metavar=('PROGRAM', 'REGEX'),
                        help="run PROGRAM on files whose names match REGEX")
    parser.add_argument('--wudb', action='store', default=wudb_name,
                        metavar="FILE", help="use FILE as Workunit DB")
    parser.add_argument('--badidealsfile', action='store',
                        metavar="FILE",
                        help="dump text of failed workunits to FILE")
    parser.add_argument('-r', '--received', action='store_true',
                        help='Process all workunits'
                             ' currently marked as received without error')
    parser.add_argument('--update',
                        action='store_true',
                        help='Update status of processed WUs'
                             ' to verified with/without error'
                             ' according to test result')
    args = parser.parse_args()

    if "wudb" in args:
        wudb_name = args.wudb

    acc = wudb.WuAccess(wudb_name)

    if "file" in args:
        # find wu for this file
        # do_check(acc, wu)
        pass

    if args.received:
        # Test all workunits whose status is RECEIVED_OK
        received = acc.query(eq={"status": wudb.WuStatus.RECEIVED_OK})
        badidealsfile = None
        print("Found " + str(len(received)) + " workunits")
        for wu in received:
            result = do_check(acc, wu, list(args.program))
            if result is None:
                print("Warning: none of the supplied programs"
                      " matched for workunit " + wu["wuid"])
            else:
                if result:
                    print("Workunit " + wu["wuid"] + " verified ok")
                else:
                    print("Workunit " + wu["wuid"] + " verified with error")
                    if args.badidealsfile:
                        if badidealsfile is None:
                            badidealsfile = open(args.badidealsfile, "a")
                        badidealsfile.write(wu["wu"] + '\n')
            if args.update:
                acc.verification(wu["wuid"], result)
        if badidealsfile is not None:
            badidealsfile.close()
