#!/usr/bin/env bash 

set -e
if [ "$CADO_DEBUG" ] ; then set -x ; fi

if ! [ "$WDIR" ] ; then
    echo "Want \$WDIR" >&2
    exit 1
fi


echo "Using mpi=$mpi; $mpirun; $mpi_extra_args"

# scan for arguments that are useful to pass to lingen_verify_checkpoints
#       $prime
#       $mpi    (we already have it)
m="$2"
n="$3"
# kmax="$4"
prime="$5"
# seed="$6"
if [ "$mpi" ] ; then
    eval "$exporter_mpirun"
    eval "$exporter_mpi_extra_args"
    # don't append mpi=$mpi to args just now, as this will also be done
    # by test-plingen.sh
    # set -- "$@" mpi="$mpi"
fi

# This is not accurate, unfortunately. There seems to be no way to do
# long integer arithmetic in pure bash.
sizeinbase2() {
    a=$1
    if [ "${#a}" -le 10 ] ; then 
        x=$(printf '%x' $1)
        len=$((4*${#x}))
        case $x in
            0) echo $len; return;;
            1*) echo $((len-3)); return;;
            [23]*) echo $((len-2)); return;;
            [4567]*) echo $((len-1)); return;;
            [89a-fA-F]*) echo $len; return;;
        esac
        echo "sizeinbase2 error $1 -> $x -> uh ?" >&2; exit 1
    else
        # a is M * 10^e, so log2(a) = log2(M) + e*log2(10), and too bad
        # for the inaccuracy we get...
        logM=$(sizeinbase2 "${a:0:10}")
        elog10=$((332*(${#a}-10)/100))
        loga=$((logM+elog10))
        echo $loga
        return
    fi
}
# this is used only for GF(p) code, which is resolutely 64-bit only
wordsize=64
nbits_prime=$(sizeinbase2 $p)
nwords=$((1+nbits_prime/wordsize))

verify() {
    set -- "${mpirun[@]}"
    mpirun_single=()
    while [ $# -gt 0 ] ; do
        mpirun_single+=("$1")
        if [ "$1" = "-n" ] ; then
            shift
            mpirun_single+=(1)
        fi
        shift
    done
    if [ "$prime" ] && [ "$prime" != 2 ] ; then
        list_tool="`dirname $0`"/../../../linalg/bwc/list_lingen_checkpoints.sh
        checks=()
        while read x ; do checks+=("$x") ; done < <("$list_tool" "$cpdir")
        for c in "${checks[@]}" ; do
            aa=($c)
            bb=()
            for x in "${aa[@]}" ; do
                bb+=("$cpdir/$x")
            done
            verifier=("$PROJECT_BINARY_DIR/linalg/bwc/lingen_verify_checkpoints_p_$nwords" "prime=$prime" m=$m n=$n)
            if [ "$mpi" ] ; then
                verifier+=(mpi=$mpi)
            fi
            "${mpirun_single[@]}" "${verifier[@]}" -- "${bb[@]}"
        done
    fi
}


args=("$@")

set -- "${args[@]}" tuning_schedule_filename="$WDIR/ts.txt"

# First decide once and for all on the schedule for the multiplications.
"`dirname $0`"/test-plingen.sh "$@" --tune

# Do a first run, and save a series of checkpoints
cpdir="$WDIR/cp"
mkdir "$cpdir"
set -- "$@" tuning_quiet=1 checkpoint_directory="$cpdir"
"`dirname $0`"/test-plingen.sh "$@"
verify

# intentionally destroy part of the checkpoints
# osx does not have -r option to xargs.
ls -l "$cpdir"/pi*aux
find "$cpdir" -name 'pi.0.*' | xargs rm -f || :
nfiles=$(ls -rt "$cpdir"/pi*aux | wc -l)
ls -t "$cpdir"/pi*aux | head -n $((2*nfiles/3)) | xargs rm -vf || :
ls -l "$cpdir"/pi*aux

# now run again, 
"`dirname $0`"/test-plingen.sh "$@"
verify
