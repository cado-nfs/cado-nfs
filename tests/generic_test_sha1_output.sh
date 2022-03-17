#!/usr/bin/env bash

: ${REFERENCE_SHA1=unknown}

if [ "$STDIN" ] ; then
    exec < "$STDIN"
fi

SHA1BIN=sha1sum
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=sha1 ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=shasum ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then
    echo "Could not find a SHA-1 checksumming binary !" >&2
    exit 1
fi

got=$("$@" | "$SHA1BIN")
got="${got%% *}"

if [ "$got" != "$REFERENCE_SHA1" ] ; then
    echo "Wrong SHA-1 obtained from $@" >&2
    echo "Got $got, expected $REFERENCE_SHA1" >&2
    exit 1
fi
