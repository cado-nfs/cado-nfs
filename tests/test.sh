#!/usr/bin/env bash

# We need this silly script inside the test definition just in order to
# catch potential precommands...

set -eo pipefail

if [ "$CADO_DEBUG" ] ; then set -x ; fi

usage() {
    echo "Usage: $0 [--filter-output <regexp>] [--expect-sha1 <hash>] -- <program> [<args> ...]" >&2
    exit 1
}

tests=()

filter_regex=.
expect_sha1=
stdin=/dev/null
abort_to_fail=

while [ $# -gt 0 ] ; do
    if [ "$1" = "--filter-output" ] ; then
        shift
        filter_regex="$1"
        shift
    elif [ "$1" = "--sed-output" ] ; then
        shift
        filter_sed="$1"
        shift
    elif [ "$1" = "--expect-sha1" ] ; then
        shift
        expect_sha1="$1"
        shift
    elif [ "$1" = "--stdin" ] ; then
        shift
        stdin="$1"
        shift
    elif [ "$1" = "--abort-to-fail" ] ; then
        shift
        abort_to_fail="$1"
    elif [ "$1" = "--" ] ; then
        shift
        break
    else
        usage
    fi
done

SHA1BIN=sha1sum
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=sha1 ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=shasum ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then
    echo "Could not find a SHA-1 checksumming binary !" >&2
    exit 1
fi

exec 3>&1

dotest() {
    if [ "$abort_to_fail" ] ; then
        if ! "${TEST_PRECOMMAND[@]}" "$@" ; then exit 1 ; else true ; fi
    else
        "${TEST_PRECOMMAND[@]}" "$@"
    fi
}


if [ "$filter_sed" ] ; then
    do_filter() { (grep "$filter_regex" || : ) | sed -e "$filter_sed" ; }
else
    do_filter() { grep "$filter_regex" || : ; }
fi

got_sha1=$(dotest "$@" <"$stdin" | tee >(cat >&3) | do_filter | $SHA1BIN)

if ! [ "$expect_sha1" ] ; then
    # echo "========= $got_sha1 ========"
    exit 0
fi

got_sha1="${got_sha1%% *}"

match=
reftab=(`echo $expect_sha1 | tr ',' ' '`)
for ref in "${reftab[@]}" ; do
    if [ "${got_sha1}" = "$ref" ] ; then
        match="$ref"
    fi
done

if ! [ "$match" ] ; then
    echo "Got sha1 checksum $got_sha1 different from expected $expect_sha1" >&2
    exit 1
fi

exit 0
