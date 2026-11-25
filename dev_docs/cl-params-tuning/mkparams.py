#!/usr/bin/env python3
"""Emit parameters/cl/params.qs.dNNN from the chosen operating points (choices.json)."""
import json, math, sys, os
from math import comb, log

HEADER = """###########################################################################
#     Parameter file for Cado-NFS
###########################################################################
# See parameters/factor/params.c90 for an example which contains some
# documentation.
# See parameters/cl/params.qs.d90 for an example which contains some
# documentation specific to the QS algorithm for class groups computation.
#
{notes}
###########################################################################
# General parameters
###########################################################################

name = d{d}

###########################################################################
# Sieve
###########################################################################

tasks.I = {I}
tasks.lim0 = {lim}
tasks.lpb0 = {lpb}
tasks.sieve.mfb0 = {mfb}

# The special-q's are products of qfac_nfac primes taken in increasing order
# from qfac_min on.  The range below suits |D| around 10^{dm}; for a given
# discriminant the best yield is obtained when the primes actually used are
# centred on (sqrt(|D|)/2^(I-1))^(1/qfac_nfac).
tasks.sieve.qfac_nfac = {nfac}
tasks.sieve.qfac_min = {qmin}
tasks.sieve.qfac_max = {qmax}

tasks.sieve.qrange = {qrange}

tasks.sieve.ncurves0 = {ncurves}

tasks.sieve.rels_wanted = {rels_wanted}

###########################################################################
# Filtering
###########################################################################

tasks.filter.target_density = 170.0
tasks.nmatrices = 2

###########################################################################
# Linalg
###########################################################################

# The determinants are large (about 1.9 bits per row of the matrix), so many
# moduli are needed; wide moduli are cheaper per bit of determinant.  Do not
# go above 704: the prime is part of the workunit file names, and those
# become too long for most filesystems.
tasks.linalg.moduli_nbits = 704
"""

def qfac_range(d, I, nfac, nsq):
    """special-q factors: split primes around g = (sqrt|D|/2^(I-1))^(1/nfac) for |D| ~ 10^(d-1/2).
    The special-q's are enumerated from the smallest primes on; c primes are enough for 4*nsq
    special-q's.  About half of the primes split, so they are 2 ln g apart on average."""
    g = math.exp((0.5 * (d - 0.5) * log(10) - (I - 1) * log(2)) / nfac)
    c = nfac
    while comb(c, nfac) < 4 * nsq:
        c += 1
    gap = 2 * log(g)
    qmin = int(g - 0.5 * c * gap)
    qmax = int(g + 2.0 * c * gap)
    return max(qmin, 50), qmax, g, c

def emit(c, outdir):
    d = c["d"]
    qmin, qmax, g, cc = qfac_range(d, c["I"], c["nfac"], c["nsq"])
    rels_wanted = int(round(c["R"] * c.get("margin", 1.13), -3))
    qrange = c.get("qrange") or max(10, min(1000, int(round(400 / c["t_sq"], -1)) or 10))
    notes = "\n".join("# " + l if l else "#" for l in c["notes"].split("\n")) + "\n"
    txt = HEADER.format(dm=f"{d - 0.5:g}", d=d, I=c["I"], lim=c["lim"], lpb=c["lpb"], mfb=2 * c["lpb"], nfac=c["nfac"], qmin=qmin,
                        qmax=qmax, qrange=qrange, ncurves=c.get("ncurves", 20), rels_wanted=rels_wanted, notes=notes)
    fn = f"{outdir}/params.qs.d{d}"
    open(fn, "w").write(txt)
    return fn, dict(qmin=qmin, qmax=qmax, g=g, c=cc, rels_wanted=rels_wanted, qrange=qrange)

if __name__ == "__main__":
    choices = json.load(open(sys.argv[1]))
    for c in choices:
        print(emit(c, sys.argv[2]))
