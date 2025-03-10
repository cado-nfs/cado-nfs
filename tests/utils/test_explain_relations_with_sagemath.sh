#!/usr/bin/env bash

set -e

if [ "$CADO_DEBUG" ] ; then
    set -x
fi

build_tree="$PROJECT_BINARY_DIR"
DL=

while [ $# -gt 0 ] ; do
    if [ "$1" = "-b" ] ; then shift ; build_tree="$1" ; shift
    elif [ "$1" = "-poly" ] ; then shift ; POLY="$1" ; shift
    elif [ "$1" = "-lim0" ] ; then shift ; LIM0="$1" ; shift
    elif [ "$1" = "-lim1" ] ; then shift ; LIM1="$1" ; shift
    elif [ "$1" = "-lpb0" ] ; then shift ; LPB0="$1" ; shift
    elif [ "$1" = "-lpb1" ] ; then shift ; LPB1="$1" ; shift
    elif [ "$1" = "-mfb0" ] ; then shift ; MFB0="$1" ; shift
    elif [ "$1" = "-mfb1" ] ; then shift ; MFB1="$1" ; shift
    elif [ "$1" = "-q0" ] ; then shift ; Q0="$1" ; shift
    elif [ "$1" = "-q1" ] ; then shift ; Q1="$1" ; shift
    elif [ "$1" = "-I" ] ; then shift ; I="$1" ; shift
    elif [ "$1" = "-dl" ] ; then shift ; DL="-dl"
    else
        echo "bad arg: $1" >&2
        exit 1
    fi
done

: ${SAGE:?missing}
: ${build_tree:?missing}
: ${wdir:?missing}
: ${POLY:?missing}
: ${LIM0:?missing}
: ${LIM1:?missing}
: ${LPB0:?missing}
: ${LPB1:?missing}
: ${MFB0:?missing}
: ${MFB1:?missing}
: ${Q0:?missing}
: ${Q1:?missing}
: ${I:?missing}

name=$(basename $POLY .poly)

"$build_tree/sieve/makefb" -poly $POLY -lim $LIM0 -maxbits $I -out "$wdir/$name.roots0.gz"
"$build_tree/sieve/makefb" -poly $POLY -lim $LIM1 -maxbits $I -out "$wdir/$name.roots1.gz"

las_args=(
    -poly $POLY 
    -q0 $Q0 
    -q1 $Q1 
    -I $I 
    -lim0 $LIM0 -lpb0 $LPB0 -mfb0 $MFB0 -fb0 "$wdir/$name.roots0.gz"
    -lim1 $LIM1 -lpb1 $LPB1 -mfb1 $MFB1 -fb1 "$wdir/$name.roots1.gz"
)

"$build_tree/sieve/las" "${las_args[@]}" > "$wdir/$name.rels_raw.txt"

RENUMBER="${wdir}/$name.renumber.gz"

# it's okay to test in DL mode.
$build_tree/sieve/freerel -poly $POLY -renumber $RENUMBER -lpbs "$LPB0,$LPB1" $DL

RELATION_FILES=("$wdir/$name.rels_raw.txt")

mkdir -p $wdir/$name.dup1/0
"$build_tree/filter/dup1" -prefix "$name.dup1" -out "$wdir/$name.dup1" -n 0 "${RELATION_FILES[@]}"

FILTERED="$wdir/$name.dup1/0/$name.dup1.0000"
nrels=$(wc -l < "$FILTERED")
"$build_tree/filter/dup2" -poly $POLY -nrels $nrels -renumber $RENUMBER $DL "$FILTERED"

"$build_tree/misc/explain_indexed_relation" -renumber "$RENUMBER" -poly "$POLY" $DL -python -relations "$FILTERED" > "$wdir/check.py"

export PYTHONUNBUFFERED=true
"$SAGE" "$wdir/check.py"


