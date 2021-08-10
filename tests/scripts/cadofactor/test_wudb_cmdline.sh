#!/usr/bin/env bash

SCRIPTDIR="$1"

# The correct string has a sneaky whitespace at the end of the first line
CORRECT_OUTPUT="All existing workunits: 
0"

if test -z "$SCRIPTDIR"
then
    echo "SCRIPTDIR parameter must be supplied"
    exit 1
fi

rm -f test.db
if ! "$SCRIPTDIR"/wudb.py -dbfile test.db -create || \
   ! OUTPUT=`"$SCRIPTDIR"/wudb.py -dbfile test.db -all -dump` ||
   [ "$OUTPUT" != "$CORRECT_OUTPUT" ]
then
    rm -f test.db
    exit 1
fi

rm -f test.db
