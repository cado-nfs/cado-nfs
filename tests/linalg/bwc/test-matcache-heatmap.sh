#!/usr/bin/env bash

# This checks that mm_bucket_heatmap produces a json file that is
# consistent with the matrix it describes. The viewer page is not
# checked.

N=100000
M=140000
dens=20
seed=1
niter=10
python="${PYTHON_EXECUTABLE:-python3}"

set -e

: ${bindir:=$PROJECT_BINARY_DIR}

usage() {
    echo "Usage: $0 [--matrix-size <N>] [--matrix-cols <M>] [--density <d>] [--seed <s>]" \
         "[--bindir <dir>] [--python <exe>]" >&2
    exit 1
}

while [ $# -gt 0 ] ; do
    if [ "$1" = "--matrix-size" ] ; then
        shift
        N=$1
        shift
    elif [ "$1" = "--matrix-cols" ] ; then
        shift
        M=$1
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
# the blocks must cover the matrix, not its transpose. The bucket code
# indexes its slices along dim[store_transposed], so a rectangular matrix
# in the -t direction catches a missing swap.
check(max(b[1] for b in J["blocks"]) == nrows,
      "blocks stop at row %d, the matrix has %d"
      % (max(b[1] for b in J["blocks"]), nrows))
check(max(b[3] for b in J["blocks"]) == ncols,
      "blocks stop at column %d, the matrix has %d"
      % (max(b[3] for b in J["blocks"]), ncols))
# the per-block times are part of the loop that "total" measures. Leave
# room for the timer noise of a very short run.
check(attributed <= J["total"] * 1.05 + 1e-6,
      "blocks account for %g s, more than the %g s of the whole loop"
      % (attributed, J["total"]))

# The vertical staircase sees every coefficient twice, and the second
# pass is attributed from the flush batches. If that walk goes wrong the
# combine time silently vanishes, so insist that it is there.
dis = J["types"].index("d-dis") if "d-dis" in J["types"] else None
if dis is not None and any(b[5] == dis for b in J["blocks"]):
    combine = sum(b[7] for b in J["blocks"])
    check(combine > 0, "the staircase blocks were given no combine time")
    check(all(b[7] > 0 for b in J["blocks"] if b[5] == dis and b[4]),
          "some staircase block was given no combine time")

sys.stderr.write("%s: %d blocks, %d coefficients, %.1f%% of the loop attributed\n"
                 % (out, len(J["blocks"]), ncoeffs, 100.0 * attributed / J["total"]))
print(ncoeffs)
ENDOFPYTHON

# deliberately not square, so that a confusion between the rows and the
# columns of a block cannot go unnoticed
$bindir/linalg/bwc/random_matrix -nrows $N -ncols $M -d $dens --binary \
    -o $wdir/mat.bin --freq -s $seed

# The three ways of cutting the matrix into blocks give three different
# shapes of the slice header tree, and the heat map walks all of them.
# small2 only ever appears as a child of small1, hence the pair. Each is
# tried in both directions, since the bucket code stores the matrix
# column-major in one of them.
reference=
for methods in small1,small2 large vsc ; do
  for direction in plain transposed ; do
    if [ "$direction" = transposed ] ; then
        # bench_matcache -t reads the transpose, so the dimensions swap
        t_arg=-t
        e_nrows=$M
        e_ncols=$N
    else
        t_arg=
        e_nrows=$N
        e_ncols=$M
    fi
    what="$methods/$direction"

    rm -f $wdir/mat.bin-bucket.bin $wdir/mat.bin-bucketT.bin
    $bindir/linalg/bwc/bench_matcache -r --nmax $niter --nchecks 1 \
        -impl bucket $t_arg \
        -matmul_bucket_methods $methods \
        -mm_bucket_heatmap $wdir/heat.json \
        $wdir/mat.bin > $wdir/bench.$methods.$direction.out 2>&1

    # the name of the local matrix is inserted in the file name, so that
    # the instances of an mpi/thread grid do not share one file
    out=$wdir/heat.mat.bin.json

    if ! [ -f "$out" ] ; then
        echo "$what: $out was not created" >&2
        cat $wdir/bench.$methods.$direction.out >&2
        exit 1
    fi

    ncoeffs=$("$python" $wdir/check.py "$out" $e_nrows $e_ncols $niter)

    # the blocks must cover the matrix exactly once, so their
    # coefficients must add up to what the matrix file holds
    expected=$(sed -n 's/^total \([0-9]*\) coeffs$/\1/p' \
        $wdir/bench.$methods.$direction.out)
    if ! [ "$expected" ] ; then
        echo "$what: could not read the coefficient count from" \
             "$wdir/bench.$methods.$direction.out" >&2
        exit 1
    fi
    if [ "$ncoeffs" != "$expected" ] ; then
        echo "$what: blocks hold $ncoeffs coefficients, the matrix has" \
             "$expected" >&2
        exit 1
    fi

    # and that must not depend on how the matrix was cut
    if [ "$reference" ] && [ "$ncoeffs" != "$reference" ] ; then
        echo "$what: $ncoeffs coefficients, but $reference with the" \
             "other cuttings" >&2
        exit 1
    fi
    reference=$ncoeffs

    rm -f "$out"
  done
done

# Now the other half: matmul_top collects the heat maps of all the
# submatrices of the mpi/thread grid into the file that the command line
# actually named. A 2x2 thread grid is enough to tell whether the pieces
# land where they belong.
if [ -x "$bindir/linalg/bwc/krylov" ] ; then
    ( cd $wdir && $bindir/linalg/bwc/krylov wdir=. thr=2x2 nullspace=left \
        interval=100 mn=64 prime=2 ys=0..64 start=0 end=100 \
        skip_online_checks=1 rebuild_cache=1 seed=1 \
        sequential_cache_build=1 no_save_cache=1 \
        random_matrix=nrows=$N,density=$dens,seed=$seed \
        mm_impl=bucket mm_bucket_heatmap=./collected.json ) \
        > $wdir/krylov.out 2>&1

    if ! [ -f "$wdir/collected.json" ] ; then
        echo "the collected heat map was not created" >&2
        tail -20 $wdir/krylov.out >&2
        exit 1
    fi

    "${python}" - "$wdir/collected.json" <<-'ENDOFPYTHON'
import json
import sys

out = sys.argv[1]
with open(out) as f:
    J = json.load(f)


def check(cond, msg):
    if not cond:
        raise SystemExit("%s: %s" % (out, msg))


check(J.get("grid"), "no grid in the collected file")
nh, nv = J["grid"]
snr, snc = J["submatrix"]
check(nh * nv == 4, "grid is %dx%d, expected 2x2" % (nh, nv))
check(nh * snr == J["nrows"],
      "%d row blocks of %d do not make %d rows" % (nh, snr, J["nrows"]))
check(nv * snc == J["ncols"],
      "%d column blocks of %d do not make %d columns" % (nv, snc, J["ncols"]))

# every submatrix of the grid must have contributed, and no block may
# straddle a boundary -- the submatrices are handled independently
seen = set()
ncoeffs = 0
for b in J["blocks"]:
    i0, i1, j0, j1, nc, t, dispatch, combine = b
    check(0 <= i0 < i1 <= J["nrows"], "bad row range %d..%d" % (i0, i1))
    check(0 <= j0 < j1 <= J["ncols"], "bad column range %d..%d" % (j0, j1))
    check(i0 // snr == (i1 - 1) // snr and j0 // snc == (j1 - 1) // snc,
          "block %d..%d x %d..%d straddles a submatrix boundary"
          % (i0, i1, j0, j1))
    seen.add((i0 // snr, j0 // snc))
    ncoeffs += nc

check(len(seen) == nh * nv,
      "only %d of the %d submatrices contributed" % (len(seen), nh * nv))
check(ncoeffs == J["ncoeffs"],
      "blocks hold %d coefficients, header says %d" % (ncoeffs, J["ncoeffs"]))
check(J["total_max"] <= J["total"] + 1e-9,
      "the busiest submatrix took more than all of them together")

sys.stderr.write("%s: %dx%d grid, %d blocks, %d coefficients\n"
                 % (out, nh, nv, len(J["blocks"]), ncoeffs))
ENDOFPYTHON
    echo "collected heatmap ok"
fi

echo "matcache heatmap ok ($reference coefficients)"
