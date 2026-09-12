#!/usr/bin/env bash

# This checks that mm_bucket_heatmap produces a json file that is
# consistent with the matrix it describes. The viewer page is not
# checked.

N=100000
dens=20
seed=1
niter=10
python="${PYTHON_EXECUTABLE:-python3}"

set -e

: ${bindir:=$PROJECT_BINARY_DIR}

usage() {
    echo "Usage: $0 [--matrix-size <N>] [--density <d>] [--seed <s>]" \
         "[--bindir <dir>] [--python <exe>]" >&2
    exit 1
}

while [ $# -gt 0 ] ; do
    if [ "$1" = "--matrix-size" ] ; then
        shift
        N=$1
        shift
    elif [ "$1" = "--density" ] ; then
        shift
        dens=$1
        shift
    elif [ "$1" = "--seed" ] ; then
        shift
        seed=$1
        shift
    elif [ "$1" = "--bindir" ] ; then
        shift
        bindir=$1
        shift
    elif [ "$1" = "--python" ] ; then
        shift
        python=$1
        shift
    else
        usage
    fi
done

: ${bindir:?missing variable}

wdir=$(mktemp -d  ${TMPDIR-/tmp}/cado-nfs.XXXXXXXX)

cleanup() { if ! [ "$CADO_DEBUG" ] ; then rm -rf $wdir ; fi ; }
trap cleanup EXIT

cat > $wdir/check.py <<'ENDOFPYTHON'
import json
import sys

out, nrows, ncols, niter = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])

with open(out) as f:
    J = json.load(f)


def check(cond, msg):
    if not cond:
        raise SystemExit("%s: %s" % (out, msg))


check(J["nrows"] == nrows, "nrows is %s, expected %d" % (J["nrows"], nrows))
check(J["ncols"] == ncols, "ncols is %s, expected %d" % (J["ncols"], ncols))
check(sum(J["iterations"]) == niter,
      "iterations %s do not add up to %d" % (J["iterations"], niter))
check(J["blocks"], "no block at all")
check(J["total"] > 0, "the cpu-bound loop took no time at all")

ncoeffs = 0
attributed = 0.0
for b in J["blocks"]:
    i0, i1, j0, j1, nc, t, dispatch, combine = b
    check(0 <= i0 < i1 <= J["nrows"], "bad row range %d..%d" % (i0, i1))
    check(0 <= j0 < j1 <= J["ncols"], "bad column range %d..%d" % (j0, j1))
    check(nc <= (i1 - i0) * (j1 - j0),
          "block %d..%d x %d..%d holds more coefficients than it has cells"
          % (i0, i1, j0, j1))
    check(0 <= t < len(J["types"]), "bad type index %d" % t)
    check(dispatch >= 0 and combine >= 0, "negative time")
    ncoeffs += nc
    attributed += dispatch + combine

check(ncoeffs == J["ncoeffs"],
      "blocks hold %d coefficients, header says %d" % (ncoeffs, J["ncoeffs"]))
# the per-block times are part of the loop that "total" measures. Leave
# room for the timer noise of a very short run.
check(attributed <= J["total"] * 1.05 + 1e-6,
      "blocks account for %g s, more than the %g s of the whole loop"
      % (attributed, J["total"]))

sys.stderr.write("%s: %d blocks, %d coefficients, %.1f%% of the loop attributed\n"
                 % (out, len(J["blocks"]), ncoeffs, 100.0 * attributed / J["total"]))
print(ncoeffs)
ENDOFPYTHON

$bindir/linalg/bwc/random_matrix -nrows $N -d $dens --binary \
    -o $wdir/mat.bin --freq -s $seed

# The three ways of cutting the matrix into blocks give three different
# shapes of the slice header tree, and the heat map walks all of them.
# small2 only ever appears as a child of small1, hence the pair.
reference=
for methods in small1,small2 large vsc ; do
    rm -f $wdir/mat.bin-bucket.bin
    $bindir/linalg/bwc/bench_matcache -r --nmax $niter --nchecks 1 \
        -impl bucket \
        -matmul_bucket_methods $methods \
        -mm_bucket_heatmap $wdir/heat.json \
        $wdir/mat.bin > $wdir/bench.$methods.out 2>&1

    # the name of the local matrix is inserted in the file name, so that
    # the instances of an mpi/thread grid do not share one file
    out=$wdir/heat.mat.bin.json

    if ! [ -f "$out" ] ; then
        echo "$methods: $out was not created" >&2
        cat $wdir/bench.$methods.out >&2
        exit 1
    fi

    ncoeffs=$("$python" $wdir/check.py "$out" $N $N $niter)

    # the blocks must cover the matrix exactly once, so their
    # coefficients must add up to what the matrix file holds
    expected=$(sed -n 's/^total \([0-9]*\) coeffs$/\1/p' $wdir/bench.$methods.out)
    if ! [ "$expected" ] ; then
        echo "$methods: could not read the coefficient count from" \
             "$wdir/bench.$methods.out" >&2
        exit 1
    fi
    if [ "$ncoeffs" != "$expected" ] ; then
        echo "$methods: blocks hold $ncoeffs coefficients, the matrix has" \
             "$expected" >&2
        exit 1
    fi

    # and that must not depend on how the matrix was cut
    if [ "$reference" ] && [ "$ncoeffs" != "$reference" ] ; then
        echo "$methods: $ncoeffs coefficients, but $reference with the" \
             "other cuttings" >&2
        exit 1
    fi
    reference=$ncoeffs

    rm -f "$out"
done

echo "matcache heatmap ok ($reference coefficients)"
