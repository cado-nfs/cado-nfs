#!/usr/bin/env bash

# This file can be sourced from shell scripts, and defines a few shell
# variables and arrays.
#
# Input: we want $1 to be one of:
#       "" (empty string)
#       mpi_rect1
#       mpi_rect2
#       mpi_square1
#       mpi_square2
#       We also want to have either $bindir or $PROJECT_BINARY_DIR point
#       to the top of the build tree.
# Output: 
#       $mpi is maked as exported, and set to the mpi geometry, or an
#       empty string if running without mpi.
#       $exporter_mpirun is an exported variable that contain bash-quoted contents
#       meant to be de-serialized by eval "$mpirun".
#       ditto for exporter_mpi_extra_args and mpi_extra_args

: ${bindir:=$PROJECT_BINARY_DIR}
: ${bindir?missing variable}

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
    nompi=1x1
    overcommit_openmpi=
    # it is also possible to set things such as mpi_rect1_mpi_args=(foo bar)
    # in order to adjust the per-config parameters.
    if [ "$jobsize" -ge 6 ] ;   then
        mpi_rect1=2x3 ; mpi_rect2=3x2 ;
    elif [ "$nnodes" -eq 1 ] && [ "$family" = openmpi ] ; then
        # we know how to overcommit.
        mpi_rect1=2x3 ; mpi_rect1_mpi_args=(--host localhost:6)
        mpi_rect2=3x2 ; mpi_rect2_mpi_args=(--host localhost:6)
        if ! [ "$OPENMPI_INSIDE_DOCKER" ] ; then
            # make sure we don't add the same stanza twice.
            mpi_rect1_mpi_args+=(--bind-to none)
            mpi_rect2_mpi_args+=(--bind-to none)
        fi
        overcommit_openmpi=1
    elif [ "$jobsize" -ge 2 ] ; then
        mpi_rect1=2x1 ; mpi_rect2=1x2 ;
    fi
    if [ "$jobsize" -ge 9 ] ; then
        mpi_square2=3x3
    elif [ "$nnodes" -eq 1 ] && [ "$family" = openmpi ] ; then
        mpi_square2=3x3 ; mpi_square2_mpi_args=(--host localhost:9)
        if ! [ "$OPENMPI_INSIDE_DOCKER" ] ; then
            # make sure we don't add the same stanza twice.
            mpi_square2_mpi_args+=(--bind-to none)
        fi
        overcommit_openmpi=1
    fi
    if [ "$jobsize" -ge 4 ] ; then
        mpi_square1=2x2
    elif [ "$nnodes" -eq 1 ] && [ "$family" = openmpi ] ; then
        mpi_square1=2x2
        mpi_square1_mpi_args=(--host localhost:4)
        if ! [ "$OPENMPI_INSIDE_DOCKER" ] ; then
            # make sure we don't add the same stanza twice.
            mpi_square1_mpi_args+=(--bind-to none)
        fi
        overcommit_openmpi=1
    fi
}

create_exporters() {
    export mpi
    # only bash5 can do that.
    if (set -e ; a=(false) ; eval "${a[@]@A}") ; then
        export exporter_mpirun="${mpirun[@]@A}"
        export exporter_mpi_extra_args="${mpi_extra_args[@]@A}"
    else
        export exporter_mpirun="mpirun=(${mpirun[*]})"
        export exporter_mpi_extra_args="mpi_extra_args=(${mpi_extra_args[*]})"
    fi
}

set_mpi_derived_variables()
{
    case "$mpi_magic" in
        nompi) ;;
        mpi_rect[12]|mpi_square[12]) : ;;
        *) echo "Unkown magic for mpi auto choice: $mpi_magic" >&2 ; return ;;
    esac

    if ! [ "$bindir" ] ; then
        echo "guess_mpi_configs.sh must be sourced with \$bindir set" >&2
    fi

    mpi_bindir=$(perl -ne '/HAVE_MPI\s*"(.*)"\s*$/ && print "$1\n";' $bindir/cado_mpi_config.h)
    mpiexec=$(perl -ne '/MPIEXEC\s*"(.*)"\s*$/ && print "$1\n";' $bindir/cado_mpi_config.h)
    mpicc=$(perl -ne '/MPI_C_COMPILER\s*"(.*)"\s*$/ && print "$1\n";' $bindir/cado_mpi_config.h)

    if [ "$mpi_magic" != nompi ] && ! [ "$mpi_bindir" ] ; then
        mpi=skip
        echo "Not an MPI build, mpi checks are disabled"
        return
    fi

    # Note that even if we don't have an mpi build, we do intend to run the
    # test if it's a nompi one !
    if [ "$mpi_magic" = nompi ] && ! [ "$mpi_bindir" ] ; then
        mpi=
        mpirun=()
        mpi_extra_args=()
        create_exporters
        return
    fi

    # A nompi test with an mpi build will need mpirun detection anyway.

    detect_mpi_family

    unset mpi_args_common

    if [ "$OAR_NODEFILE" ] ; then
        nnodes=$(wc -l < $OAR_NODEFILE)
        if [ "$family" = openmpi ] ; then
            mpi_args_common=(-machinefile $OAR_NODEFILE)
            mpi_extra_args+=(--map-by node --bind-to none --mca plm_rsh_agent oarsh)
        fi
        # maybe auto-detect some hardware and decide on the right config
        # options ? should we really have to do that ?
    elif [ "$SLURM_NPROCS" ] ; then
        nnodes=$SLURM_NPROCS
    else
        nnodes=1
    fi

    ncores=$(egrep '^core[[:space:]]+id[[:space:]]+:' /proc/cpuinfo | sort -u | wc -l)

    if [ "$CI_JOB_NAME" ] && [ "$family" = openmpi ] ; then
        # we get failures similar to what is reported there
        # https://github.com/open-mpi/ompi/issues/4948
        mpi_extra_args+=(--mca btl_vader_single_copy_mechanism none)

        # See https://github.com/horovod/horovod/issues/1985 for the
        # rationale of the fix below.
        OPENMPI_INSIDE_DOCKER=1
        mpi_extra_args+=(--bind-to none)
    fi

    case "$nnodes,$ncores,$family" in
        1,*,openmpi) 
            mpi_extra_args+=(-mca mtl ^psm2,ofi,cm --mca btl ^openib)
            set_choices_from_n $ncores
            if ! [ "$overcommit_openmpi" ] && ! [ "$OPENMPI_INSIDE_DOCKER" ] ; then
                mpi_extra_args+=(--bind-to core)
            fi
            ;;
        *,openmpi) 
            set_choices_from_n $nnodes
            ;;
        # this works only when impi is configured to used mpiexec. If we
        # are in a configuration that requires PMI / PMIx (and in fact,
        # the same question holds for openmpi), then we must do something
        # different.
        *,impi) 
            # see https://gitlab.inria.fr/cado-nfs/cado-nfs/-/merge_requests/123#note_898921
            set_choices_from_n $nnodes
            ;;
        *)
            echo "Script does not know which mpi tests to enable (nnode=$nnodes ncores=$ncores mpi_family=$family)"
            ;;
    esac
    mpi="${!mpi_magic}"
    if [ "$mpi_magic" != nompi ] && ! [ "$mpi" ] ; then
        echo "No MPI run possible for magic choice $mpi_magic ; no-op exit" >&2
        mpi=skip
        return
    fi
    if [ "$mpi" ] ; then
        set `echo $mpi | tr 'x' ' '`
        if ! [ "$1" ] || ! [ "$2" ] ; then
            # should never happen.
            echo "Bad test configuration, mpi should be of the form \d+x\d+ for MPI test" >&2
            exit 1
        fi
        njobs=$(($1*$2))
    else
        njobs=1
    fi
    mpi_args_common+=(-n $njobs)
    _t="${mpi_magic}_mpi_args"[@]
    mpi_extra_args+=("${!_t}")
    # mpirun=("$mpiexec" "${mpi_args_common[@]}" "${!_t}" "${mpi_extra_args[@]}")
    mpirun=("$mpiexec" "${mpi_args_common[@]}" "${mpi_extra_args[@]}")
    # pass on to subcalls
    create_exporters
}
