#!/usr/bin/env bash

mpiconfs="$1"
shift

. "`dirname $0`"/guess_mpi_configs.sh

for mpi_magic in "${mpiconfs[@]}" ; do
    if [ "$mpi_magic" ] && [ "$mpi_magic" != nompi ] ; then
        mpi=
        mpirun=()
        unset exporter_mpirun
        unset exporter_mpi_extra_args
        . "`dirname $0`/guess_mpi_configs.sh"
        set_mpi_derived_variables
        if [ "$mpi" ] ; then
            if [[ "$1" =~ bwc.pl ]] ; then
                set -- "$@" mpi="$mpi"
                if [[ "$*" =~ :mpirun ]] ; then
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
                else
                    # we probably shouldn't pass mpi_extra_args at all,
                    # in fact.
                    set -- "$@"mpi_extra_args="${mpi_extra_args[*]}"
                fi
            fi
            "$@"
        fi
        # if mpi_magic was set, but mpi wasn't, then we don't run.
    else
        "$@"
    fi
done
