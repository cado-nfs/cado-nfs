#!/usr/bin/env bash

CHECK_RELS="$1"
SOURCE_TEST_DIR="`dirname "$0"`"
: ${wdir?missing}


poly="${SOURCE_TEST_DIR}/c60.poly"
rels1="${SOURCE_TEST_DIR}/c60.rels1"
rels2="${SOURCE_TEST_DIR}/c60.rels2"

out="$wdir/out"
fixed="$wdir/fixed"

set -e

if [ "$CADO_DEBUG" ] ; then set -x ; fi

############## Check with correct relations
echo "$0: recognizing correct relations"
"$CHECK_RELS" -poly $poly -lpb0 18 -lpb1 19 -check_primality $rels1 > $out 2>&1
echo "$0: recognizing correct relations: binary ran OK"

nrels=`awk '/read relations/ {print $5}' $out`
correctrels=`awk '/correct relations/ {print $5}' $out`
if [ $nrels != $correctrels ] ; then
    echo "$0: check_rels binary failed to recognize correct relations"
    exit 1
fi
echo "$0: recognizing correct relations: OK"


############## Check with wrong relations
echo "$0: testing wrong relations"
! "$CHECK_RELS" -poly $poly -v -lpb0 18 -lpb1 18 -check_primality -fixit -out $fixed $rels2 > $out 2>&1
echo "$0: testing wrong relations: binary ran OK"

nrels=`awk '/read relations/ {print $5}' $out`
correctrels=`awk '/correct relations/ {print $5}' $out`
fixedrels=`awk '/fixed relations/ {print $5}' $out`
wrongrels=`awk '/of wrong relations/ {print $5}' $out`
composite=`awk '/1 non-prime factor/ {print $2}' $out`
lpb=`awk '/ideal larger than a lpb/ {print $2}' $out`
if [ $nrels != "5" -o $correctrels != "0" -o $fixedrels != "3" ] ; then
    echo "$0: check_rels binary failed to fix wrong relations"
    exit 1
fi
if [ $wrongrels != "2" ] ; then
    echo "$0: check_rels binary failed to detect wrong relations"
    exit 1
fi
if [ $composite != "1" ] ; then
    echo "$0: check_rels binary failed to detect composite factor"
    exit 1
fi
if [ $lpb != "1" ] ; then
    echo "$0: check_rels binary failed to detect larger than lpb factor"
    exit 1
fi
