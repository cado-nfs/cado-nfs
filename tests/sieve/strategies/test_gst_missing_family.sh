#!/usr/bin/env bash

# gst splits its input file into the six factoring method families it
# knows about, and used to read the first element of each without ever
# checking that the family is there at all. Depending on which family was
# missing, that was a SIGSEGV, or an empty point set handed to the convex
# hull computation. Check that each of them is now diagnosed instead.
#
# The input is a pinned sample file rather than something produced by the
# gfm/benchfm chain: benchfm selects methods from wall-clock measurements
# taken on the spot, so that chain does not give the same output twice.

set -e

: ${CADO_NFS_SOURCE_DIR:?missing}
: ${CADO_NFS_BINARY_DIR:?missing}
: ${wdir:?missing}

GST="$CADO_NFS_BINARY_DIR/sieve/strategies/gst"
sample="$CADO_NFS_SOURCE_DIR/tests/sieve/strategies/sample_methods.txt"

lim0=1024
lpb0=16
r0=22

cd "$wdir"

$GST -gdc -lim0 $lim0 -mfb0 $r0 -out decomp_${lim0}_${r0}

run_gst_r() {
    # the -gst_r mode of gst, on the method file given as $1
    mkdir -p out
    $GST -gst_r -lim0 $lim0 -lpb0 $lpb0 -r0 $r0 -ncurves 3 \
        -in "$1" -decomp decomp_${lim0}_${r0} -out out
}

# a complete file must still go through, and produce the output file.
run_gst_r "$sample"
if ! [ -f out/strategies${lim0}_${r0} ] ; then
    echo "gst produced no strategies file on a complete input" >&2
    exit 1
fi

# Each case below is an extended regexp matching the lines of one family,
# and the name gst is expected to use for it. "4 1" is ECM-B12, "4 2" is
# ECM-M12, "4 4" is ECM-M16. Removing PP1-27 ("2 ") or PP1-65 ("3 ")
# alone is fine, so the two go together: gst only uses their union.
cases=(
    "1 :PM1"
    "2 |3 :PP1-27 or PP1-65"
    "4 1:ECM-B12"
    "4 2:ECM-M12"
    "4 4:ECM-M16"
)

rc=0
for c in "${cases[@]}" ; do
    family="${c%%:*}"
    name="${c#*:}"
    grep -E -v "^($family)" "$sample" > amputated.txt
    rm -rf out
    status=0
    run_gst_r amputated.txt > stdout.txt 2> stderr.txt || status=$?
    if [ $status = 0 ] ; then
        echo "gst accepted an input file with no $name method" >&2
        rc=1
    elif [ $status -ge 128 ] ; then
        # not merely cosmetic: a runtime that kills the process on an
        # uncaught exception need not print what() first, and libc++
        # does not.
        echo "gst died on a signal on a file with no $name method" >&2
        rc=1
    elif ! grep -q "contains no $name method" stderr.txt ; then
        echo "gst failed on a file with no $name method," \
            "but without saying so:" >&2
        tail -3 stderr.txt >&2
        rc=1
    fi
done

exit $rc
