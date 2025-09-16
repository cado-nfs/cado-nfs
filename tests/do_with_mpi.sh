#!/usr/bin/env bash

if [ "$CADO_DEBUG" ] ; then
    set -x
fi

set -x

if [ "$1" = --handle-mpirun-directly ] ; then
    handle_mpirun=1
    shift
fi

mpiconfs="$1"
shift

. "`dirname $0`"/guess_mpi_configs.sh

if [ "${#mpiconfs[@]}" -gt 1 ] ; then
    echo "script does not like multi configs" >&2
    # main question is how do we deal with return codes...
    exit 1
fi

for mpi_magic in "${mpiconfs[@]}" ; do
    mpi=
    mpirun=()
    unset exporter_mpirun
    unset exporter_mpi_extra_args
    . "`dirname $0`/guess_mpi_configs.sh"
    set_mpi_derived_variables
    # hmmm. maybe it's the job of the callee to use these variables. It
    # very much depends on what we're calling, after all.
    # eval "$exporter_mpirun"
    # eval "$exporter_mpi_extra_args"
    export mpirun
    if [ "$mpi" = skip ] ; then
        exit 125
    fi
    if [[ "$1" =~ bwc.pl ]] ; then
        if [ "$mpi" ] ; then
            set -- "$@" mpi="$mpi"
        fi
        if [[ "$*" =~ :mpirun ]] ; then
            if [ "$mpi" ] ; then
                # If :mpirun is found, we add the mpi= thing just
                # after it _as well_.
                nargs=()
                for x in "$@" ; do
                    nargs+=("$x")
                    if [ "$x" = :mpirun ] ; then
                        nargs+=(mpi="$mpi")
                        nargs+=(mpi_extra_args="${mpi_extra_args[*]}")
                    fi
                done
                set -- "${nargs[@]}"
            fi
        else
            # we probably shouldn't pass mpi_extra_args at all,
            # in fact.
            set -- "$@" mpi_extra_args="${mpi_extra_args[*]}"
        fi
    fi
    if [ "$handle_mpirun" ] ; then
        "${mpirun[@]}" "$@"
    else
        "$@"
    fi
done
