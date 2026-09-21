#!/usr/bin/env bash

# This checks that build_matcache and bench_matcache have consistent
# behaviour.

N=100
dens=10
seed=1

set -e

: ${bindir:=$PROJECT_BINARY_DIR}

usage() {
    echo "Usage: $0 <N>" >&2
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
    elif [ "$1" = "--variable-width" ] ; then
        # a simd width that no fixed-width arithmetic backend serves, so
        # that arith_generic::instance() falls back to the
        # variable-width layer. Empty if that layer is not compiled in.
        shift
        variable_width=$1
        shift
    else
        usage
    fi
done

: ${bindir:?missing variable}

if ! [ "$N" ] ; then usage ; fi

SHA1BIN=sha1sum
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=sha1 ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=shasum ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then
    echo "Could not find a SHA-1 checksumming binary !" >&2
    exit 1
fi


wdir=$(mktemp -d  ${TMPDIR-/tmp}/cado-nfs.XXXXXXXX)

cleanup() { if ! [ "$CADO_DEBUG" ] ; then rm -rf $wdir ; fi ; }
trap cleanup EXIT

$bindir/linalg/bwc/random_matrix -nrows $N -d $dens  --binary -o $wdir/mat.bin --freq -s $seed

bench_arg_left=-t
cachesuffix_left=T

bench_arg_right=
cachesuffix_right=

for impl in basic sliced bucket ; do
    for direction in left right ; do
        eval cachefiledirection=\$cachesuffix_${direction}
        cachefile=$wdir/mat-${impl}${cachefiledirection}.bin

        $bindir/linalg/bwc/build_matcache --matrix-file $wdir/mat.bin -impl $impl -direction $direction -tmpdir $wdir > $wdir/build.$impl.$direction.out 2>&1
        SHA1=$($SHA1BIN ${cachefile})
        SHA1="${SHA1%% *}"
        REFSHA1=$SHA1
        unset SHA1


        cachefile2=$wdir/mat.bin-${impl}${cachefiledirection}.bin

        eval argtail=\(\$bench_arg_${direction}\)
        $bindir/linalg/bwc/bench_matcache -r  --nmax 100 -impl $impl "${argtail[@]}" $wdir/mat.bin > $wdir/bench$impl.$direction.out 2>&1

        SHA1=$($SHA1BIN ${cachefile2})
        SHA1="${SHA1%% *}"

        if [ "$SHA1" != "$REFSHA1" ] ; then
            echo "weird, build_matcache and bench_matcache computed different things ?" >&2
            exit 1
        fi

        echo "matcache $impl $direction ok"
    done
done

# The bucket implementation cuts the matrix in several ways, and each of
# them has its own multiplication code. The default cutting of a small
# matrix only ever produces small slices, so ask for the other ones
# explicitly. Widths are swept as well: the variable-width layer used to
# get all of its pointer arithmetic and all of its cache sizing wrong,
# and nothing noticed because the tests only ever asked for widths that
# a fixed-width backend serves.
for nbys in 64 $variable_width ; do
    for methods in small1 small1,small2 large vsc ; do
        for direction in left right ; do
            eval argtail=\(\$bench_arg_${direction}\)
            rm -f $wdir/mat.bin-bucket.bin $wdir/mat.bin-bucketT.bin
            out=$wdir/bench.bucket.$nbys.$methods.$direction.out
            $bindir/linalg/bwc/bench_matcache -r --nmax 2 --nchecks 2 \
                -impl bucket -nbys $nbys \
                matmul_bucket_methods=$methods \
                "${argtail[@]}" $wdir/mat.bin > $out 2>&1 || :
            if ! grep -q "^All 2 checks passed" $out ; then
                echo "bucket nbys=$nbys methods=$methods $direction: check failed" >&2
                cat $out >&2
                exit 1
            fi
            echo "matcache bucket nbys=$nbys methods=$methods $direction ok"
        done
    done
done

