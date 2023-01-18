#!/usr/bin/env bash

set -eo pipefail

if [ "$CADO_DEBUG" ] ; then set -x ; fi

usage() {
    echo "Usage: $0 [--filter-output <regexp>] [--expect-sha1 <hash>] -- <program> [<args> ...]" >&2
    exit 1
}

tests=()

filter_regex=.
expect_sha1=

while [ $# -gt 0 ] ; do
    if [ "$1" = "--filter-output" ] ; then
        shift
        filter_regex=$1
        shift
    elif [ "$1" = "--expect-sha1" ] ; then
        shift
        expect_sha1=$1
        shift
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

got_sha1=$("$@" | grep "$filter_regex" | $SHA1BIN)

if ! [ "$expect_sha1" ] ; then
    exit 0
fi

got_sha1="${got_sha1%% *}"

if [ "$got_sha1" != "$expect_sha1" ] ; then
    echo "Got sha1 checksum=$got_sha1, different from expected $expect_sha1" >&2
    exit 1
fi
