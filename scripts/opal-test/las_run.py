from opal.core.io import *
import sys
import subprocess
import re
from math import log
import os

from report import LasStats

def update_existing(a, b):
    """ Update those keys in "a" which also exist in "b" with the values
    from "b"
    """
    for key in set(a) & set(b):
        a[key] = b[key]

# return prime_pi(2^x)
def primepi(x):
    # precomputed values up to x=40 (obtained with Sage)
    l = [0, 1, 2, 4, 6, 11, 18, 31, 54, 97, 172, 309, 564, 1028, 1900, 3512, 6542, 12251, 23000, 43390, 82025, 155611, 295947, 564163, 1077871, 2063689, 3957809, 7603553, 14630843, 28192750, 54400028, 105097565, 203280221, 393615806, 762939111, 1480206279, 2874398515, 5586502348, 10866266172, 21151907950, 41203088796]
    return l[x]

# return the best time so far (or 1e308 at the beginning)
def get_best_time(file):
    try:
        f = open(file, "r")
        best_time = float(f.read())
        f.close()
    except (OSError, IOError):
        best_time = 1e308
    return best_time

def set_best_time(file, t):
    f = open(file, "w")
    f.write(str(t) + "\n")
    f.close()

def run(param_file, problem):
    "Run las with given parameters until the required number of relations is found."

    build_dir = os.environ["CADO_BUILD"]
    cado_poly_type = os.environ["OPAL_CADO_TYPE"].split(',')

    makefb = "%s/sieve/makefb" % build_dir
    las = "%s/sieve/las" % build_dir
    sqside = os.environ["OPAL_CADO_SQSIDE"].strip()

    las_params = {
        "I": 11,
        "poly": "c59.polyselect2.poly",
        "lim0": 50000,
        "lim1": 100000,
        "lpb0": 22,
        "lpb1": 22,
        "mfb0": 22,
        "mfb1": 22,
        "ncurves0": 6,
        "ncurves1": 6,
        "qmin": 100000,
        "bkthresh1": 1000000,
        "t": 2 # number of threads for las
    }
    if sqside != "":
      las_params["sqside"] = sqside


    params = read_params_from_file(param_file)

    # Update parameters for las
    update_existing(las_params, params)
    las_params["mfb0"] = max(las_params["mfb0"], las_params["lpb0"])
    las_params["mfb1"] = max(las_params["mfb1"], las_params["lpb1"])

    las_params["lim0"] = min(las_params["lim0"], 2 ** las_params["lpb0"])
    las_params["lim1"] = min(las_params["lim1"], 2 ** las_params["lpb1"])
    
    to_print = ["I", "lim1", "lpb1", "mfb1", "lim0", "lpb0", "mfb0", "ncurves0", "ncurves1", "qmin", "bkthresh1"]
    sys.stderr.write("Using parameters %s\n" % " ".join(["%s:%s" % (key, las_params[key]) for key in to_print]))

    ### Call makefb on all algebraic polynomials
    makefb_params = {
      "poly" : las_params["poly"],
      "maxbits": las_params["I"],
      "t": las_params["t"]
    }
    for side, t in enumerate(cado_poly_type):
      if t == "alg":
        makefb_params["lim"] = las_params["lim%d"%side]
        makefb_params["out"] = "factorbase.roots%d.gz" % side
        makefb_params["side"] = side
        las_params["fb%s"%side] = makefb_params["out"]

        makefb_cmd_line = [makefb]
        for (key, value) in makefb_params.items():
          makefb_cmd_line += ["-%s" % key, str(value)]
        # sys.stderr.write("Running: %s\n" % " ".join(makefb_cmd_line))
        subprocess.check_call(makefb_cmd_line)

    stats = LasStats()
    q0 = las_params["qmin"]
    q_range = 1000
    q_inc = 0
    # since we remove duplicates, 80% of the total ideals should be enough
    rels_wanted = int(0.8 * primepi(las_params["lpb1"]) + 0.8 * primepi(las_params["lpb0"]))
    # read the best time so far (if any)
    best_time = get_best_time ("las.best")
    if best_time >= 1e308:
       sys.stderr.write("Estimate %u unique relations needed\n" % rels_wanted)
    else:
       sys.stderr.write("Estimate %u relations needed, best time so far is %.0f seconds\n" % (rels_wanted, best_time))

    while stats.get_rels() < rels_wanted:
        # Set q0, q1
        outputfile = "las.%d.out" % q0
        las_params.update({"q0": q0, "q1": q0 + q_range, "out": outputfile})

        las_cmd_line = [las, "-allow-largesq", "-never-discard", "-dup", "-dup-qmin"]
        if sqside == 0:
            las_cmd_line += ["%s,0" % str(las_params["qmin"])]
        else:
            las_cmd_line += ["0,%s" % str(las_params["qmin"])]
        for (key, value) in las_params.items():
            if key == "qmin":
               continue
            las_cmd_line += ["-%s" % key, str(value)]
        # sys.stderr.write("Running: %s\n" % " ".join(las_cmd_line))
        las_process = subprocess.Popen(las_cmd_line)
        (stdout, stderr) = las_process.communicate()
        if las_process.returncode != 0:
            raise Exception("las exited with status %d" % las_process.returncode)
        if not stats.parse_one_file(outputfile):
            raise Exception("Could not read statistics from %s" % outputfile)
        os.unlink(outputfile)
        cur_rels = stats.get_rels()
        if cur_rels == 0:
           cur_time = 0
        else:
           cur_time = stats.get_time(cur_rels)
        # sys.stderr.write("   Up to q=%u, estimate %u/%u relations in %.0f seconds\n" % (q0, cur_rels, rels_wanted, cur_time))

        # if the total time so far exceeds best_time, and we don't have enough
        # relations, we can stop here
        if cur_time > best_time and cur_rels < rels_wanted:
           break

        # set q_inc so that we do about 10 sieving tests of length q_range
        if q_inc == 0:
           v = stats.relations_int.lastvalue[0]
           q_inc = int(rels_wanted / (10.0 * v))
           # but make sure q_inc is not too large: q_inc = 0.45*q0 gives
           # about 20 sieving tests for [q0,10*q0]
           if q_inc > int(0.45 * q0):
               q_inc = int(0.45 * q0)
        q0 += q_inc
    qmax = stats.get_qmax(rels_wanted)
    sievetime = stats.get_time(rels_wanted)
    if sievetime < best_time: # store the new best time
        set_best_time ("las.best", sievetime)
    sys.stderr.write("Estimated special-q up to %.0f and %.0f seconds for %d relations\n"
                      % (qmax, sievetime, rels_wanted))
    # Estimate how far the sieving should have gone to get the desired number
    # of relations
    return {'SIEVETIME': sievetime, 'RELATIONS': rels_wanted}

if __name__ == '__main__':
    param_file  = sys.argv[1]
    problem     = sys.argv[2]
    output_file = sys.argv[3]

    # Solve, gather measures and write to file.
    measures = run(param_file, problem)
    write_measures_to_file(output_file, measures)
