#!/usr/bin/env bash

# This shell script will do the following.
#
# - create a temp directory.
# - run the given command line, by providing the information about the
#   location of the temp directory in any desired way.
# - clean up the temp directory

# If CADO_DEBUG is set to something, the temp directory is not cleaned.


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
        extra=("${extra[@]}" "$1=$t")
        shift
    elif [ "$1" = "--other" ] ; then
        t=`mktemp -d $TMPDIR/XXXXXXXXXXXXXX`
        temps=("${temps[@]}" "$t")
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
        main=("${main[@]}" "$1")
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
