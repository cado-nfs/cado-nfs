#!/usr/bin/env python3

# Example for the p59 in tests/test_full_p59
# scripts/hintfile-helper.py
#   --cadobindir /tmp/cado-nfs-build-master
#   --ntrials 100
#   --descent-hint-table ~/Local/p59.hint
#   "27a 0 0 I=15 400000,47,65,1.4 400000,48,94,2.0"

# where the last line has exactly the same format as what should be put
# in a descent hint table, except that the second and third items are zero.


import os
import subprocess
import argparse
import re
import tempfile

wdir=None
ntrials = 100
las_binary=None
las_args={}
tmphintfilename=None
hintfile_configs={}
last_active_config=[None,None]
verbose=False

def make_selector(side, bitsize):
    return "%d@%d" % (bitsize,side)

def parse_selector(qq):
    foo=re.match("^(\d+)([ra]|@\d+)$",qq)
    if not foo:
        return
    bitsize, side = foo.groups()
    bitsize = int(bitsize)
    if side == 'r':
        side = 0
    elif side == 'a':
        side = 1
    else:
        side = int(side[1:])
    return side, bitsize

def write_hintfile_inner(f):
    selectors=list(hintfile_configs.keys())
    selectors.sort()
    for sel in selectors:
        hc=hintfile_configs[sel]
        implied,config=hc[:2]
        time=0
        if len(hc) > 2:
            time=hc[2]
        suxs=0
        if len(hc) > 3:
            suxs=hc[3]
        f.write("%s %.4f %.4f %s\n" % (sel, time, suxs, config))

def write_temporary_hintfile():
    # write the current status of the hint file to a temp file.
    fd,tmphintfilename = tempfile.mkstemp('.txt', wdir + '/')
    tmphintfile=os.fdopen(fd, 'w')
    write_hintfile_inner(tmphintfile)
    tmphintfile.close()
    return tmphintfilename

def write_hintfile(filename):
    f=open(filename, 'w')
    write_hintfile_inner(f)
    f.close()

def parse_config_line(configline):
    match = re.match(
            "^((\d+)([ra]|@\d+))\s+([\d\.]+)\s+([\d\.]+)\s+(I=(\d+)\s+(\d+),(\d+),([\d.]+)\s+(\d+),(\d+),([\d.]+))\s*$",
            configline)
    if not match:
        return
    selector, bitsize, side, time, suxs, confstring, I, lim0, lpb0, mix0, lim1, lpb1, mix1 = match.groups()

    time = float(time)
    suxs = float(suxs)

    bitsize = int(bitsize)
    if side == 'r':
        side = 0
    elif side == 'a':
        side = 1
    else:
        side = int(side[1:])
    # canonicalize the selector
    selector=make_selector(side,bitsize)

    implied={}

    implied['-I'] = int(I)

    # FIXME: this is really a problem. We give the impression that these
    # parameters are important, even though they are not. In truth, las
    # will look inside the hint table for the parameters to use.
    implied['--lim0'] = int(lim0)
    implied['--lpb0'] = int(lpb0)

    implied['--lim1'] = int(lim1)
    implied['--lpb1'] = int(lpb1)

    if float(mix0) >= 10:
        implied['--mfb0'] = int(float(mix0))
    else:
        implied['--mfb0'] = int(float(mix0)*int(lpb0))
        implied['--lambda0'] = float(mix0)
    if float(mix1) >= 10:
        implied['--mfb1'] = int(float(mix1))
    else:
        implied['--mfb1'] = int(float(mix1)*int(lpb1))
        implied['--lambda1'] = float(mix1)

    return implied, selector, confstring, time, suxs


def testoneline(selector, implied, *args, **kwargs):

    argdict=las_args.copy()

    tmphintfilename = write_temporary_hintfile()
    tmprelfilename = tempfile.mktemp('.rels', wdir + '/')

    argdict['-out'] = tmprelfilename

    for key,value in implied.items():
        argdict[key]=str(value)


    side,bitsize = parse_selector(selector)
    argdict['--q0'] = int(2**(bitsize-1))
    argdict['--q1'] = int(2**bitsize)
    argdict['--sqside'] = str(side)

    argseq=[las_binary]
    for key,value in argdict.items():
        argseq += [key, str(value)]
    argseq += [
        "--descent-hint-table", tmphintfilename,
        "--prepend-relation-time",
        "--exit-early", "1",
        "--never-discard",  # this is a kludge, we're risking a segfault.
                            # see bug 15617
        "-allow-largesq",
        ]

    if verbose:
        print("las command line:\n" + " ".join(argseq))
    subprocess.check_call(argseq)

    # now parse the output
    f=open(tmprelfilename, 'r')
    nok=0
    nfail=0
    tok=0
    current=None
    for line in f.readlines():
        foo=re.match("^# Sieving\s+(.*q=\d+;\s+rho=\d+)", line)
        if foo:
            if current:
                nfail+=1
            current=[foo.groups()[0]]
            continue
        elif current:
            # maybe we've foud a relation !
            foo=re.match("^\(([\d\.]+)\)", line)
            if foo:
                tok += float(foo.groups()[0])
                nok+=1
                current = None
    if current:
        nfail+=1
    f.close()

    time = tok/nok
    suxs = 1.0*nok/ntrials

    print("%s %.4f %.4f %s" % (selector, time, suxs, confstring))

    os.unlink(tmprelfilename)
    os.unlink(tmphintfilename)


    return time, suxs


if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Descent hint file setup for DLP")
    # Required
    parser.add_argument("--cadobindir", help="Cado build directory",
            required=True, type=str)
    parser.add_argument("--datadir", help="cadofactor working directory")
    parser.add_argument("--prefix", help="prefix for data files")
    parser.add_argument("--poly", help="poly file name")
    parser.add_argument("--fb", help="factor base file name")
    parser.add_argument("--workdir", help="Temporary working directory")
    parser.add_argument("--descent-hint-table", help="Preliminary hint file")
    parser.add_argument("--ntrials", type=int, help="Number of special-q's to generate and test")
    parser.add_argument("--qrange", help="Comma-separated list of special-q sizes to optimize")
    parser.add_argument("--hintline", type=str, help="Example hint line to be evaluated")
    parser.add_argument("--replace", action="store_true")
    parser.add_argument("--uselogtable", action="store_true")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    wdir = args.workdir
    if wdir is None:
        wdir = tempfile.mkdtemp(dir="/tmp")
    if args.ntrials:
        ntrials=args.ntrials
    las_args['--random-sample'] = str(ntrials)
    las_binary = args.cadobindir + "/sieve/las_descent"
    if args.verbose:
        verbose=True

    if args.datadir and args.prefix:
        if (args.poly or args.fb):
            raise ValueError("Please specify either datadir&prefix or poly&fb")
        las_args["--poly"] = args.datadir+"/"+args.prefix+".poly"
        las_args["--fb1"]  = args.datadir+"/"+args.prefix+".roots.gz"
        if args.uselogtable:
            las_args["--renumber"] = args.datadir+"/"+args.prefix+".renumber.gz"
            las_args["--log"] = args.datadir+"/"+args.prefix+".dlog"
    elif args.poly and args.fb:
        if (args.datadir or args.prefix):
            raise ValueError("Please specify either datadir&prefix or poly&fb")
        las_args["--poly"] = args.poly
        las_args["--fb1"] = args.fb1
        if args.uselogtable:
            raise ValueError("--uselogtable requires --datadir & --prefix")
    else:
        raise ValueError("Please specify either datadir&prefix or poly&fb")


    qranges=[]        # list of (side,bitsize)

    if args.qrange:
        if args.hintline:
            raise ValueError("--qrange and hintline are incompatible")
        qranges=[]
        # parsing of the qrange is complicated.
        # it's a comma-separated list of ranges
        last_bitsizes=[0,0]
        for qq in sum([a.split() for a in args.qrange.split(',')],[]):
            side,bitsize = parse_selector(qq)
            if bitsize <= last_bitsizes[side]:
                raise ValueError("Arguments for --qrange not in increasing order for side %d" % side)
            last_bitsizes[side] = bitsize
            qranges.append((side,bitsize))


    # read the preliminary hint file. It may contain some config defaults
    # which we are going to use.
    if args.descent_hint_table:
        oldhintfile = open(args.descent_hint_table, 'r')
        for l in oldhintfile.readlines():
            res = parse_config_line(l)
            if not res:
                continue
            implied, selector, confstring, time, suxs = res
            hintfile_configs[selector] = (implied, confstring, time, suxs)
        oldhintfile.close()

    if args.hintline:
        if args.qrange is not None:
            raise ValueError("--qrange and hintline are incompatible")
        res = parse_config_line(args.hintline)
        if not res:
            raise ValueError("Error parsing %s"%args.hintline)
        implied, selector, confstring = res[:3]
        hintfile_configs[selector] = (implied, confstring)
        side,bitsize = parse_selector(selector)
        qranges.append((side,bitsize))


    for qsqs in qranges:
        side,bitsize = qsqs
        config=None
        selector = make_selector(side,bitsize)
        try:
            config=hintfile_configs[selector]
            last_active_config[side]=config
        except KeyError:
            config=last_active_config[side]
            if not config:
                b=bitsize
                while b>0:
                    try:
                        config=hintfile_configs[make_selector(side, b)]
                        break
                    except KeyError:
                        b -= 1
                if not config:
                    raise ValueError("cannot find configuration for %s"% selector)
        implied, confstring = config[:2]

        time, suxs = testoneline(selector, implied)

        hintfile_configs[selector] = (implied, confstring, time, suxs)


    if args.workdir is None:
        os.rmdir(wdir)

    if args.replace:
        if not args.descent_hint_table:
            raise ValueError("--replace requires --descent-hint-table")
        write_hintfile(args.descent_hint_table)


