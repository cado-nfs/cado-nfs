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
    else
        usage
    fi
done

: ${bindir?missing variable}

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

$bindir/linalg/bwc/random_matrix $N -d $dens  --binary -o $wdir/mat.bin -s $seed

# $bindir/linalg/bwc/mf_scan mfile=$wdir/mat.bin --binary-in --freq

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
        $bindir/linalg/bwc/bench_matcache -r $wdir/mat.bin  --nmax 100 -impl $impl "${argtail[@]}" > $wdir/bench$impl.$direction.out 2>&1

        SHA1=$($SHA1BIN ${cachefile2})
        SHA1="${SHA1%% *}"

        if [ "$SHA1" != "$REFSHA1" ] ; then
            echo "weird, build_matcache and bench_matcache computed different things ?" >&2
            exit 1
        fi

        echo "matcache $impl $direction ok"
    done
done

