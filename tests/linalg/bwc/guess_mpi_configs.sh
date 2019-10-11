#!/bin/bash

# This file can be sourced from shell scripts, and defines a few shell
# variables and arrays.

# The $bindir must be set to the top of the build_tree

# not sure that it makes sense to have this script distinct from bwc.pl

detect_mpi_family()
{
    output="`$mpicc -v 2>&1`"
    family=
    suffix=
    if [[ $output =~ ^mpi.*for.MVAPICH2.version[[:space:]](.*) ]] ; then
        family=mvapich2
        suffix="${BASH_REMATCH[1]}"
    elif [[ $output =~ ^mpi.*for.MPICH.*version[[:space:]](.*) ]] ; then
        family=mpich
        suffix="${BASH_REMATCH[1]}"
    elif [[ $output =~ ^mpi.*for.*Intel.*MPI.*Library[[:space:]]([0-9].*)[[:space:]]for.* ]] ; then
        family=impi
        suffix="${BASH_REMATCH[1]}"
    else
        output0="$output"
        output="`$mpicc -showme:version 2>&1`"
        if [[ $output =~ ^.*Open[[:space:]]MPI[[:space:]]([^[:space:]]*) ]] ; then
            family=openmpi
            suffix="${BASH_REMATCH[1]}"
        else
            echo "MPI C compiler front-end not recognized, proceeding anyway (output0 was: $output0 ; later was $output)" >&2
        fi
    fi
}


set_choices_from_n()
{
    jobsize="$1"
    # it is also possible to set things such as _mpi_rect1_mpi_args=(foo bar)
    # in order to adjust the per-config parameters.
    if [ "$jobsize" -ge 6 ] ;   then _mpi_rect1=2x3 ; _mpi_rect2=3x2 ;
    elif [ "$jobsize" -ge 2 ] ; then _mpi_rect1=2x1 ; _mpi_rect2=1x2 ; fi
    if [ "$jobsize" -ge 9 ] ;   then _mpi_square2=3x3 ; fi
    if [ "$jobsize" -ge 4 ] ;   then _mpi_square1=2x2 ; fi
}

set_mpi_derived_variables()
{
    if ! [ "$mpi_magic" ] ; then
        # We typically expect some of our magic strings such as "_mpi_rect1",
        # "_mpi_rect2", "_mpi_square1", "_mpi_square2"
        return
    fi

    case "$mpi_magic" in
        _mpi_rect[12]|_mpi_square[12]) : ;;
        *) echo "Unkown magic for mpi auto choice: $mpi_magic" >&2 ; return ;;
    esac

    if ! [ "$bindir" ] ; then
        echo "guess_mpi_configs.sh must be sourced with \$bindir set" >&2
    fi

    mpi_bindir=$(perl -ne '/HAVE_MPI\s*"(.*)"\s*$/ && print "$1\n";' $bindir/cado_mpi_config.h)
    mpiexec=$(perl -ne '/MPIEXEC\s*"(.*)"\s*$/ && print "$1\n";' $bindir/cado_mpi_config.h)
    mpicc=$(perl -ne '/MPI_C_COMPILER\s*"(.*)"\s*$/ && print "$1\n";' $bindir/cado_mpi_config.h)

    if ! [ "$mpi_bindir" ] ; then
        echo "Not an MPI build, mpi checks are disabled"
        return
    fi

    detect_mpi_family

    if [ "$OAR_NODEFILE" ] ; then
        nnodes=$(wc -l < $OAR_NODEFILE)
        mpi_args_common=(-machinefile $OAR_NODEFILE --map-by node --mca plm_rsh_agent oarsh)
        # maybe auto-detect some hardware and decide on the right config
        # options ? should we really have to do that ?
    elif [ "$SLURM_NPROCS" ] ; then
        nnodes=$SLURM_NPROCS
    else
        nnodes=1
    fi

    ncores=$(egrep -c '^processor[[:space:]]+:' /proc/cpuinfo)

    case "$nnodes,$ncores,$family" in
        1,*,openmpi) 
            set_choices_from_n $ncores
            ;;
        *,openmpi) 
            set_choices_from_n $nnodes
            ;;
        *)
            echo "Script does not know which mpi tests to enable"
            ;;
    esac
    mpi="${!mpi_magic}"
    if ! [ "$mpi" ] ; then
        echo "No MPI run possible for magic choice $mpi_magic ; no-op exit" >&2
        exit 0
        return
    fi
    _t="${mpi_magic}_mpi_args"[@]
    set `echo $mpi | tr 'x' ' '`
    if ! [ "$1" ] || ! [ "$2" ] ; then
        # should never happen.
        echo "Bad test configuration, mpi should be of the form \d+x\d+ for MPI test" >&2
        exit 1
    fi
    njobs=$(($1*$2))
    mpi_args_common+=(-n $njobs)
    mpirun=("$mpiexec" "${mpi_args_common[@]}" "${!_t}" "${mpi_extra_args[@]}")
    # pass on to subcalls
    export mpi
    export mpirun
}
