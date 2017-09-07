#!/usr/bin/env bash

CHECK_RELS="$1"
SOURCE_TEST_DIR="`dirname "$0"`"
WORKDIR=`mktemp -d ${TMPDIR-/tmp}/cadotest.XXXXXXXX`
chmod a+rx "${WORKDIR}"

poly="${SOURCE_TEST_DIR}/c60.poly"
rels1="${SOURCE_TEST_DIR}/c60.rels1"
rels2="${SOURCE_TEST_DIR}/c60.rels2"

out="$WORKDIR/out"
fixed="$WORKDIR/fixed"


############## Check with correct relations
if ! ("$CHECK_RELS" -poly $poly -lpb0 18 -lpb1 19 -check_primality $rels1 > $out 2>&1) ; then
    echo "$0: check_rels binary failed to recognize correct relations"
    exit 1
fi

nrels=`awk '/read relations/ {print $5}' $out`
correctrels=`awk '/correct relations/ {print $5}' $out`
if [ $nrels != $correctrels ]; then
    echo "$0: check_rels binary failed to recognize correct relations"
    exit 1
fi


############## Check with wrong relations
if ("$CHECK_RELS" -poly $poly -v -lpb0 18 -lpb1 18 -check_primality -fixit -out $fixed $rels2 > $out 2>&1) ; then
    echo "$0: check_rels binary failed when testing wrong relations"
    exit 1
fi

nrels=`awk '/read relations/ {print $5}' $out`
correctrels=`awk '/correct relations/ {print $5}' $out`
fixedrels=`awk '/fixed relations/ {print $5}' $out`
wrongrels=`awk '/of wrong relations/ {print $5}' $out`
composite=`awk '/1 non-prime factor/ {print $2}' $out`
lpb=`awk '/ideal larger than a lpb/ {print $2}' $out`
if [ $nrels != "5" -o $correctrels != "0" -o $fixedrels != "3" ]; then
    echo "$0: check_rels binary failed to fix wrong relations"
    exit 1
fi
if [ $wrongrels != "2" ]; then
    echo "$0: check_rels binary failed to detect wrong relations"
    exit 1
fi
if [ $composite != "1" ]; then
    echo "$0: check_rels binary failed to detect composite factor"
    exit 1
fi
if [ $lpb != "1" ]; then
    echo "$0: check_rels binary failed to detect larger than lpb factor"
    exit 1
fi

rm -rf $WORKDIR
