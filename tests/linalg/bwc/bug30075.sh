#!/usr/bin/env bash

if [ "$CADO_DEBUG" ] ; then
    set -x
fi

: ${wdir?missing}
: ${PROJECT_BINARY_DIR?missing}

test_matrix() {
    one="\x01\x00\x00\x00"
    zero="\x00\x00\x00\x00"
    f="\x00\x00\x00"
    echo -ne "${one}A$f${one}B$f${one}D$f${one}E$f${one}G$f"
    echo -ne "${one}H$f${one}J$f${one}K$f${one}M$f${one}N$f"
    for i in {0..85} ; do echo -ne "$zero" ; done
}

test_matrix > $wdir/a.bin
$PROJECT_BINARY_DIR/linalg/bwc/mf_scan mfile=$wdir/a.bin --freq
echo -ne "$zero" >> $wdir/a.cw.bin

$PROJECT_BINARY_DIR/linalg/bwc/mf_bal 1 1 $wdir/a.bin out=$wdir/bal.bin reorder=none skip_decorrelating_permutations=1

for impl in basic sliced bucket ; do
    $PROJECT_BINARY_DIR/tests/linalg/bwc/spmv_test wdir=$wdir/$impl prime=2 matrix=$wdir/a.bin nullspace=left mm_impl=$impl no_save_cache=1 skip_bw_early_rank_check=1 m=64 n=64 balancing=$wdir/bal.bin seed=1
done

### $PROJECT_BINARY_DIR/tests/linalg/bwc/short_matmul -p 2 $wdir/a.bin $wdir/basic/Y0-64.0 $wdir/MY0-64.ref.bin

set -e
for impl in basic sliced bucket ; do
    # make sure that we generated the same vector
    diff $wdir/$impl/Y0-64.0 $wdir/basic/Y0-64.0
    diff $wdir/$impl/MY0-64.0 $wdir/basic/MY0-64.0
done
