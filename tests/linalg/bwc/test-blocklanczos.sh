#!/usr/bin/env bash


nrows=100
ncols=100
density=10
seed=1
nullspace=left
bwc_extra=()
mpithr_args=()
mf_bal_extra=()
nh=1
nv=1

set -e

if [ "$CADO_DEBUG" ] ; then set -x ; fi

usage() {
    echo "Usage: $0 <N>" >&2
    exit 1
}

while [ $# -gt 0 ] ; do
    if [ "$1" = "-b" ] ; then
        shift
        CADO_NFS_BINARY_DIR="$1"
        shift
    elif [ "$1" = "--matrix-size" ] ; then
        shift
        nrows=`echo $1 | cut -d, -f1`
        ncols=`echo $1 | cut -d, -f2`   # = nrows if no comma
        shift
    elif [ "$1" = "--density" ] ; then
        shift
        density=$1
        shift
    elif [ "$1" = "--seed" ] ; then
        shift
        seed=$1
        shift
    elif [ "$1" = "--nullspace" ] ; then
        shift
        bwc_extra+=("nullspace=$1")
        nullspace=$1
        shift
    else
        case "$1" in
            thr=*|mpi=*)
                x=`echo $1 | cut -d= -f2 | cut -dx -f1`
                nh=$((x*nh))
                x=`echo $1 | cut -d= -f2 | cut -dx -f2`
                nv=$((x*nv))
                mpithr_args+=("$1")
                bwc_extra+=("$1")
                shift
                ;;
            *) usage;;
        esac
    fi
done

: "${CADO_NFS_BINARY_DIR:?missing}"

redirect_unless_debug() {
    file="$1"
    shift
    if [ "$CADO_DEBUG" ] ; then
        "$@" 2>&1 | tee "$file"
    else
        "$@" > $file 2>&1
    fi
}

if ! [ "$nrows" ] ; then usage ; fi

wdir=$(mktemp -d  ${TMPDIR-/tmp}/cado-nfs.XXXXXXXX)
cleanup() { if ! [ "$CADO_DEBUG" ] ; then rm -rf $wdir ; fi ; }
trap cleanup EXIT

kleft=128
if [ "$kleft" -ge "$((nrows/4))" ] ; then
    kleft=$((nrows/4))
fi

$CADO_NFS_BINARY_DIR/linalg/bwc/random_matrix $nrows $ncols $density seed=1 kleft=$kleft --binary --output $wdir/mat.bin --freq

file_is_zero() {
    tt=$(mktemp "$wdir/XXXXXXXXXXXXXX")
    tr -d '\000' < "$1" > "$tt"
    ! [ -s "$tt" ]
    rc=$?
    rm -f "$tt"
    return $rc
}

# We pass --rectangular because that is the way we say that we do not
# care about replicating permutations.
$CADO_NFS_BINARY_DIR/linalg/bwc/mf_bal $nh $nv $wdir/mat.bin skip_decorrelating_permutation=1 out=$wdir/bal.bin --rectangular

$CADO_NFS_BINARY_DIR/linalg/bwc/bwc.pl :mpirun "${mpithr_args[@]}" -- $CADO_NFS_BINARY_DIR/linalg/bwc/blocklanczos m=64 n=64 ys=0..64 matrix=mat.bin wdir=$wdir balancing=bal.bin seed=1 interval=$((nrows/10+1)) no_save_cache=1 "${bwc_extra[@]}" skip_bw_early_rank_check=1 wdir=$wdir

if [ $nullspace = left ] ; then
    $CADO_NFS_BINARY_DIR/tests/linalg/bwc/short_matmul -t $wdir/mat.bin $wdir/blsolution0-64.0 $wdir/blsolution0-64.0.image
    if ! file_is_zero $wdir/blsolution0-64.0.image ; then
        echo "Apparently we have v*M != 0 ; BAD !" >&2
        exit 1
    elif file_is_zero $wdir/blsolution0-64.0 ; then
        if [ $nrows -gt $ncols ] ; then
            echo "Found trivial kernel although M is singular ; BAD !" >&2
            echo "Or it could be that Ker M \\cap {\\langle y \\rangle \\cup Im M} = 0" >&2
            exit 1
        else
            echo "Found trivial kernel, put perhaps M ($nrows*$ncols) has full rank" >&2
            exit 1
        fi
    fi
else
    $CADO_NFS_BINARY_DIR/tests/linalg/bwc/short_matmul $wdir/mat.bin $wdir/blsolution0-64.0 $wdir/blsolution0-64.0.image
    if ! file_is_zero $wdir/blsolution0-64.0.image ; then
        echo "Apparently we have M*v != 0 ; BAD !" >&2
        exit 1
    elif file_is_zero $wdir/blsolution0-64.0 ; then
        if [ $nrows -lt $ncols ] ; then
            echo "Found trivial kernel although M is singular ; BAD !" >&2
            echo "Or it could be that Ker M \\cap {\\langle y \\rangle \\cup Im M} = 0" >&2
            exit 1
        else
            echo "Found trivial kernel, put perhaps M ($nrows*$ncols) has full rank" >&2
            exit 1
        fi
    fi
fi
