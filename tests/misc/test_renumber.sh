#!/usr/bin/env bash

set -e

: ${wdir:?missing}

while [ $# -gt 0 ] ; do
    if [ "$1" = "-b" ] ; then
        shift
        BUILD_DIR="$1"
        shift
    elif [ "$1" = "-poly" ] ; then
        shift
        POLY="$1"
        shift
    elif [ "$1" = "-lpbs" ] ; then
        shift
        LPBS="$1"
        shift
    else
        echo "bad arg: $1" >&2
        exit 1
    fi
done

if ! [ "$BUILD_DIR" ] ; then
    BUILD_DIR=${PROJECT_BINARY_DIR?missing}
fi
: ${POLY:?missing}
: ${LPBS:?missing}

FREEREL="${BUILD_DIR}/sieve/freerel"
DEBUG_RENUMBER="${BUILD_DIR}/misc/debug_renumber"
python="${PYTHON_EXECUTABLE:-python3}"

# Three ways of saving the same table: the binary format in a plain file
# (which the reader mmaps), the binary format through a compressor
# (which it cannot), and the legacy text format.
${FREEREL} -poly ${POLY} -lpbs "$LPBS" -renumber ${wdir}/renumber.bin
${FREEREL} -poly ${POLY} -lpbs "$LPBS" -renumber ${wdir}/renumber.gz
${FREEREL} -poly ${POLY} -lpbs "$LPBS" -renumber ${wdir}/renumber.flat.gz \
           -renumber_format flat
${FREEREL} -poly ${POLY} -lpbs "$LPBS" -renumber ${wdir}/renumber.flat \
           -renumber_format flat

# decompressing a binary table must give back something that can be
# mmapped: the header is padded in that case too
gzip -dc ${wdir}/renumber.gz > ${wdir}/renumber.ungz

for f in renumber.bin renumber.gz renumber.flat.gz renumber.flat renumber.ungz ; do
    ${DEBUG_RENUMBER} -poly ${POLY} -renumber ${wdir}/$f -check \
        > ${wdir}/$f.dump
    grep -v '^#' ${wdir}/$f.dump > ${wdir}/$f.data
done

# the plain binary file must go through mmap, the plain text file must
# be parsed by all threads, and the compressed ones can do neither
grep -q "entries mmapped from" ${wdir}/renumber.bin.dump
grep -q "entries mmapped from" ${wdir}/renumber.ungz.dump
grep -q "entries parsed from" ${wdir}/renumber.flat.dump
! grep -q "entries mmapped from" ${wdir}/renumber.gz.dump
! grep -q "entries parsed from" ${wdir}/renumber.flat.gz.dump

# and all three must describe the very same table
diff ${wdir}/renumber.bin.data ${wdir}/renumber.gz.data
diff ${wdir}/renumber.bin.data ${wdir}/renumber.flat.gz.data
diff ${wdir}/renumber.bin.data ${wdir}/renumber.flat.data
diff ${wdir}/renumber.bin.data ${wdir}/renumber.ungz.data
cmp ${wdir}/renumber.bin ${wdir}/renumber.ungz

${DEBUG_RENUMBER} -poly ${POLY} -renumber ${wdir}/renumber.bin -check -quiet

# The last line of the header says where the data begins. A header that
# lies about it, or that announces an offset that does not suit the
# entries, must be rejected rather than believed. (Finding that line
# means counting bytes in a file that is part text and part binary,
# which no portable shell tool does.)
patch_offset() {        # $1 = what to add to the offset, $2 = output file
    "${python}" - "${wdir}/renumber.bin" "$2" "$1" <<-'ENDOFPYTHON'
	import sys
	d = open(sys.argv[1], "rb").read()
	# the placement line is the one that ends where it says the data
	# begins, which is what identifies it without any marker
	pos = 0
	while True:
	    end = d.index(b"\n", pos) + 1
	    if d[pos:pos+16].isdigit() and int(d[pos:pos+16]) == end:
	        break
	    pos = end
	off = int(d[pos:pos+16]) + int(sys.argv[3])
	open(sys.argv[2], "wb").write(d[:pos] + b"%016d" % off + d[pos+16:])
	ENDOFPYTHON
}

expect_failure() {      # $1 = file, $2 = what the complaint should say
    # the complaint comes out as an uncaught exception, so don't let it
    # leave a core behind
    if ( ulimit -c 0
         ${DEBUG_RENUMBER} -poly ${POLY} -renumber $1 -quiet > $1.out 2>&1 )
    then
        echo "reading $1 should have failed" >&2
        exit 1
    fi
    if ! grep -q "$2" $1.out ; then
        echo "expected the complaint about $1 to mention \"$2\", got:" >&2
        cat $1.out >&2
        exit 1
    fi
}

patch_offset 8 ${wdir}/renumber.lying
expect_failure ${wdir}/renumber.lying "but the header ends at offset"

patch_offset 1 ${wdir}/renumber.odd
expect_failure ${wdir}/renumber.odd "is not a multiple of"

# A table whose entries are not the width this build uses must be read
# and converted. Rewrite the blob with the other width, and patch the
# elementsize field of the header to match. (The test tables are small,
# so their values fit in 32 bits either way.)
"${python}" - "${wdir}/renumber.bin" "${wdir}/renumber.other" <<-'ENDOFPYTHON'
	import struct, sys
	d = open(sys.argv[1], "rb").read()
	# the placement line is the one that ends where it says the data
	# begins; the entry width is the field just after the offset
	pos = 0
	while True:
	    end = d.index(b"\n", pos) + 1
	    if d[pos:pos+16].isdigit() and int(d[pos:pos+16]) == end:
	        break
	    pos = end
	off = int(d[pos:pos+16])
	w = int(d[pos+17:pos+18])
	other = 8 if w == 4 else 4
	blob = d[off:]
	fmt = "I" if w == 4 else "Q"
	vals = struct.unpack("<%d%s" % (len(blob) // w, fmt), blob)
	out = struct.pack("<%d%s" % (len(vals), "I" if other == 4 else "Q"), *vals)
	head = d[:pos+17] + str(other).encode() + d[pos+18:off]
	open(sys.argv[2], "wb").write(head + out)
	ENDOFPYTHON
${DEBUG_RENUMBER} -poly ${POLY} -renumber ${wdir}/renumber.other -check \
    > ${wdir}/renumber.other.dump
grep -q "entries converted from" ${wdir}/renumber.other.dump
grep -v '^#' ${wdir}/renumber.other.dump > ${wdir}/renumber.other.data
diff ${wdir}/renumber.bin.data ${wdir}/renumber.other.data
