#!/usr/bin/env bash

set -e

if [ "$CADO_DEBUG" ] ; then set -x ; fi

: ${wdir?missing}
: ${CADO_NFS_BINARY_DIR?missing}
: ${CADO_NFS_SOURCE_DIR?missing}
: ${test_id:=$1}
: ${test_id:=sample1}

case "$test_id" in
    # By default, this script runs with a very small test matrix.
    sample1)
# Note that the d.bin file is output in sorted form by construction, so
# that comparison only makes sense if this one is already sorted.
cat > $wdir/a.txt <<EOF
1 0
3 1 3 32768
1 2
6 0 1 2 3 4 65569
EOF
"$CADO_NFS_SOURCE_DIR/tests/linalg/bwc/matrix_txt2bin.pl"  $wdir/a.txt
    ;;

    # Another way to use this script is with the sample2 matrix, which is
    # a little bit bigger.
sample2)
(cat <<EOF
2 0 131072
3 1 3 32768
1 2
EOF
for i in `seq 3 65536` ; do echo 1 $i ; done
echo 6 0 1 2 3 4 65569
) > $wdir/a.txt
"$CADO_NFS_SOURCE_DIR/tests/linalg/bwc/matrix_txt2bin.pl"  $wdir/a.txt
    ;;

    # Last, we can work with an existing .bin file, which we transform to
    # our desired format. The conversion from .bin to .rows in perl is
    # dog slow (less than 2MB/sec), but for a small test it's still ok. We
    # don't have any such test in the test suite, and don't intend to run
    # one routinely anyway.
*.bin)
    base=${test_id%%.bin}
    if ! [ -f "$base.cw.bin" ] || ! [ -f "$base.bin" ] ; then
        echo "Cannot test $test_id: not all files are present (.cw.bin .bin)" >&2
        exit 1
    fi
    # create a (sorted) .rows files from the .bin file, and the
    # .row_offsets and .col_offsets. The .rw file is not needed.
    read -d  -r perlcode1 <<'EOF'
binmode(STDIN);
binmode(STDOUT);
my $r = 0;
open(R, ">&3"); binmode(R);
open(B, ">&4"); binmode(B);
while (sysread(STDIN, my $x, 4)) {
    my $w = unpack("L", $x);
    syswrite(B, $x, 4);
    $r += $w;
    syswrite(R, pack("L", $r), 4);
    my @s = ();
    for(my $i = 0 ; $i < $w ; $i++) {
        sysread(STDIN, my $x, 4);
        push @s, unpack("L", $x);
    }
    for my $j (sort { $a <=> $b } @s) {
        $x = pack("L", $j);
        syswrite(STDOUT, $x, 4);
        syswrite(B, $x, 4);
    }
}

EOF
    read -d  -r perlcode2 <<'EOF'
open(CW, "<&STDIN");
open(C, ">&STDOUT");
binmode(CW);
binmode(C);
my $c = 0;
while (sysread(CW, my $x, 4)) {
    my $w = unpack("L", $x);
    $c += $w;
    syswrite(C, pack("L", $c), 4);
}

EOF
    perl -e "$perlcode1" < "$base.bin" > "$wdir/a.rows" 3> "$wdir/a.row_offsets" 4> "$wdir/a.bin"
    perl -e "$perlcode2" < "$base.cw.bin" > "$wdir/a.col_offsets"
    ;;
*)
    echo "Unexpected test sample format: $1" >&2
    exit 1
    ;;
esac




SHA1BIN=sha1sum
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=sha1 ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=shasum ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then
    echo "Could not find a SHA-1 checksumming binary !" >&2
    exit 1
fi

checkit() {
    b="$1"
    b=$(basename "$b")
    b=${b%%.*}
    "$CADO_NFS_BINARY_DIR/tests/linalg/bwc/test_matrix_file" read "$@" | while read f sha ; do
    if [ -f "$wdir/$b.$f" ] ; then
        read filesha filename <<<$($SHA1BIN "$wdir/$b.$f")
        if [ "$filesha" != "$sha" ] ; then
            echo "$filename has wrong sha1" >&2
            exit 1
        else
            echo "$filename sha1 ok"
        fi
    fi
done
}

"$CADO_NFS_BINARY_DIR/tests/linalg/bwc/test_matrix_file" read $wdir/a.rows --binary-mixed $wdir/b.bin --binary-data $wdir/b.rows --binary-offsets $wdir/b.row_offsets

diff -q $wdir/a.rows $wdir/b.rows
diff -q $wdir/a.row_offsets $wdir/b.row_offsets
diff -q $wdir/a.bin $wdir/b.bin
ln -s a.col_offsets $wdir/b.col_offsets

checkit $wdir/a.rows 

"$CADO_NFS_BINARY_DIR/tests/linalg/bwc/test_matrix_file" read $wdir/a.rows --columns --binary-mixed $wdir/c.bin --binary-data $wdir/c.cols --binary-offsets $wdir/c.col_offsets

diff -q $wdir/a.col_offsets $wdir/c.col_offsets
ln -s a.row_offsets $wdir/c.row_offsets

checkit $wdir/c.cols --columns

"$CADO_NFS_BINARY_DIR/tests/linalg/bwc/test_matrix_file" read $wdir/c.cols --binary-mixed $wdir/d.bin --binary-data $wdir/d.rows --binary-offsets $wdir/d.row_offsets

diff -q $wdir/a.bin $wdir/d.bin
diff -q $wdir/a.row_offsets $wdir/d.row_offsets
ln -s a.col_offsets $wdir/d.col_offsets
diff -q $wdir/a.rows $wdir/d.rows

checkit $wdir/d.rows

