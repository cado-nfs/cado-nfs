#!/usr/bin/env bash

# This shell script will do the following.
#
# - create a temp directory.
# - run the given command line, by providing the information about the
#   location of the temp directory in any desired way (see below for
#   typical uses)
# - clean up the temp directory.  If CADO_DEBUG is set to something, the
#   temp directory is not cleaned.
#
# Typical uses:
#
# The provided temp dir can be given to the subcommand via a command line
# option (anything that starts with a dash is interpreted as such).
#
#       ./provide-wdir.sh --arg --wdir ./foo.sh [args]
#
#       creates a temp dir (say $T) and run ./foo.sh [args] --wdir $T
#
# If the token that follows --arg does not start with a dash, the command
# line addition looks like an assignment instead. Note that several
# arguments can be set to the created temp directory.
#
#       ./provide-wdir.sh --arg wdir --arg x ./foo.sh [args]
#
#       creates a temp dir (say $T) and run ./foo.sh [args] wdir=$T x=$T
#
# The provided temp dir can be given to the subcommand via an environment
# variable.
#
#       ./provide-wdir.sh --env wdir ./foo.sh [args]
#
#       creates a temp dir (say $T) and run env wdir=$T ./foo.sh [args]
#
#
# It is also possible to create several temp directories.
#
#       ./provide-wdir.sh --env wdir --other --env wdir2 ./foo.sh [args]
#
#       creates two temp dirs (say $T and $T2) and run env wdir=$T wdir2=$T2 ./foo.sh [args]
#

: ${TMPDIR:=/tmp}
t=`mktemp -d $TMPDIR/cado-nfs.XXXXXXXXXXXXXX`
temps=("$t")

cleanup() { rm -rf "${temps[@]}" ; }

if [ "$CADO_DEBUG" ] ; then
    echo "debug mode, data will be left in $t"
    set -x
else
    trap cleanup EXIT
fi

extra=()

while [ $# -gt 0 ] ; do
    if [ "$1" = "--" ] ; then
        shift
        break;
    elif [ "$1" = "--env" ] ; then
        shift
        var="$1"
        eval "$var=$t"
        eval "export $var"
        shift
    elif [ "$1" = "--arg" ] ; then
        shift
        case "$1" in
            -*) extra+=("$1" "$t");;
            *) extra+=("$1=$t");;
        esac
        shift
    elif [ "$1" = "--other" ] ; then
        t=`mktemp -d $TMPDIR/cado-nfs.XXXXXXXXXXXXXX`
        temps+=("$t")
        if [ "$CADO_DEBUG" ] ; then
            echo "debug mode, data will be left in $t (in addition to other temp directories above)"
        fi
        shift
    else
        # Then we finish processing, presumably the user omitted the --
        # separator.
        break
    fi
done

# add our extra parameters before the first occurrence of the -- separator, if we happen to find any.
main=()

while [ $# -gt 0 ] ; do
    if [ "$1" = "--" ] ; then
        break
    else
        main+=("$1")
        shift
    fi
done

"${main[@]}" "${extra[@]}" "$@"

rc=$?

if [ "$CADO_DEBUG" ] ; then
    echo "debug mode, data left in ${temps[@]}"
    exit $rc
fi

if [ $rc != 0 ] ; then
    trap - EXIT
    # what do we do when the script failed ?
    rm -rf "${temps[@]}"
fi

exit $rc
