#!/usr/bin/env bash

# This tests that the matrix times vector product does what it is
# expected to do, and also checks that the permutations do what they are
# expected to do as well.

: ${TMPDIR=/tmp}
nrows=100
ncols=100
density=10
seed=1
: ${bindir:=$PROJECT_BINARY_DIR}
bwc_extra=()
mf_bal_extra=()
nh=1
nv=1
prime=2

# inject the variables that were provided by guess_mpi_configs
if [ "$mpi" ] ; then
    eval "$exporter_mpirun"
    eval "$exporter_mpi_extra_args"
fi

set -e

if [ "$CADO_DEBUG" ] ; then set -x ; fi

usage() {
    echo "Usage: $0 <N>" >&2
    exit 1
}

while [ $# -gt 0 ] ; do
    if [ "$1" = "--matrix-size" ] ; then
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
    elif [ "$1" = "--bindir" ] ; then
        shift
        bindir=$1
        shift
    elif [ "$1" = "--arith-layer" ] ; then
        shift
        arith_layer=$1
        shift
    elif [ "$1" = "--backends" ] ; then
        shift
        backends=($1)
        shift
    elif [ "$1" = "--prime" ] ; then
        shift
        prime=$1
        shift
    else
        case "$1" in
            thr=*|mpi=*)
                x=`echo $1 | cut -d= -f2 | cut -dx -f1`
                nh=$((x*nh))
                x=`echo $1 | cut -d= -f2 | cut -dx -f2`
                nv=$((x*nv))
                bwc_extra+=("$1")
                shift
                ;;
            *) echo "Unexpected arg: $1" >&2 ; usage;;
        esac
    fi
done

redirect_unless_debug() {
    file="$1"
    shift
    if [ "$CADO_DEBUG" ] ; then
        if ! "$@" > >(tee "$file") 2>&1 ; then
            echo "Failed command: $@" >&2
            exit 1
        fi
    else
        if ! "$@" > $file 2>&1 ; then
            echo "Failed command: $@" >&2
            exit 1
        fi
    fi
}

if ! [ "$nrows" ] ; then usage ; fi

wdir=$(mktemp -d  $TMPDIR/cado-nfs.XXXXXXXX)
cleanup() { if ! [ "$CADO_DEBUG" ] ; then rm -rf $wdir ; fi ; }
argh() { echo "Failed on command error" >&2 ; cleanup ; }
trap cleanup EXIT
trap argh ERR

extra_args_random_matrix=()
extra_args_mf=()
if [ $prime != 2 ] ; then 
    extra_args_random_matrix=(-c 2)
    extra_args_mf=(--withcoeffs)
fi

redirect_unless_debug $wdir/random_matrix.out $bindir/linalg/bwc/random_matrix $nrows $ncols -d $density  --binary -o $wdir/mat.bin -s $seed --freq "${extra_args_random_matrix[@]}"

if [ $nrows != $ncols ] ; then
    mf_bal_extra=(--rectangular --reorder both)
else
    mf_bal_extra=(--reorder columns)
fi
# mf_bal_extra+=(skip_decorrelating_permutation=true)
redirect_unless_debug $wdir/bal.out $bindir/linalg/bwc/mf_bal $nh $nv $wdir/mat.bin out=$wdir/bal.bin "${mf_bal_extra[@]}" "${extra_args_mf[@]}"
B=$wdir/bal.bin

echo "## arith layer is $arith_layer"
echo "## backends to test: ${backends[*]}"

if [ "$prime" = 2 ] ; then
    case "$arith_layer" in
        bz) m=1024; n=256;;
        b*) n=`echo $arith_layer | cut -c2-`; m=$n;;
        *) echo "unknown arithmetic layer $arith_layer" >&2; exit 1;;
    esac
    bwc_common=(m=$m n=$n)
    nullspace_values=(left RIGHT)
    splitwidth=64
else
    bwc_common=(m=1 n=1)
    nullspace_values=(LEFT right)
    splitwidth=1
fi


for impl in "${backends[@]}" ; do
    rm -f $wdir/Xa0-$splitwidth.0 $wdir/Xb0-$splitwidth.0
    rm -f $wdir/Xa0-$splitwidth.1 $wdir/Xb0-$splitwidth.1
    rm -f $wdir/XTa0-$splitwidth.0 $wdir/XTb0-$splitwidth.0
    rm -f $wdir/XTa0-$splitwidth.1 $wdir/XTb0-$splitwidth.1
    rm -f $wdir/MY0-$splitwidth.0 $wdir/sMY0-$splitwidth.0
    rm -f $wdir/WM0-$splitwidth.0 $wdir/sWM0-$splitwidth.0

    # The nullspace argument is just for selecting the preferred
    # direction for the matrix times vector product. But in the spmv_test
    # file, we unconditionally compute M times Y in the file MY, and W
    # times M in the file WM. Which of these products uses the "fast" or
    # "slow" code depdens on $ns.
    for ns in ${nullspace_values[@]} ; do
        redirect_unless_debug $wdir/spmv-$impl-left.out $bindir/linalg/bwc/bwc.pl :mpirun ${bwc_extra} -- $bindir/tests/linalg/bwc/spmv_test wdir=$wdir "${bwc_common[@]}" prime=$prime balancing=$B matrix=$wdir/mat.bin nullspace=$ns mm_impl=$impl no_save_cache=1 "${bwc_extra[@]}" skip_bw_early_rank_check=1
        # check done within the C code.
        # diff -q $wdir/Z0-$splitwidth.0 $wdir/ZI0-$splitwidth.0
        # diff -q $wdir/Z0-$splitwidth.0 $wdir/ZII0-$splitwidth.0
        diff -q $wdir/Xa0-$splitwidth.0 $wdir/Xb0-$splitwidth.0
        diff -q $wdir/Xa0-$splitwidth.1 $wdir/Xb0-$splitwidth.1
        diff -q $wdir/XTa0-$splitwidth.0 $wdir/XTb0-$splitwidth.0
        diff -q $wdir/XTa0-$splitwidth.1 $wdir/XTb0-$splitwidth.1
    done
    $bindir/tests/linalg/bwc/short_matmul -p $prime $wdir/mat.bin  $wdir/Y0-$splitwidth.0  $wdir/sMY0-$splitwidth.0 > /dev/null 2>&1
    diff -q $wdir/MY0-$splitwidth.0 $wdir/sMY0-$splitwidth.0
    $bindir/tests/linalg/bwc/short_matmul -p $prime -t $wdir/mat.bin  $wdir/W0-$splitwidth.0  $wdir/sWM0-$splitwidth.0 > /dev/null 2>&1
    diff -q $wdir/WM0-$splitwidth.0 $wdir/sWM0-$splitwidth.0
    echo "spmv ${nrows}x${ncols} $impl left ok ${bwc_extra[@]}"
done
