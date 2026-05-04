#!/usr/bin/env bash

# This test is here to check that the filter_galois binary rewrite relations
# correctly.

set -e

build_tree="$PROJECT_BINARY_DIR"

# the sha1 hash is changed by the normalization of the relations.
# REFERENCE_SHA1="9a1705ed3505fd6b7ae68df51870440bb225eeee"
REFERENCE_SHA1="aa9060aebbeaf0062af0e5ee9c9f79bb2213415b"

while [ $# -gt 0 ] ; do
    if [ "$1" = "-b" ] ; then
        shift
        build_tree="$1"
        shift
    elif [ "$1" = "-S" ] ; then
        shift
        REFERENCE_SHA1="$1"
        shift
    else
        echo "bad arg: $1" >&2
        exit 1
    fi
done

FG="$build_tree/filter/filter_galois"
SOURCE_TEST_DIR="`dirname "$0"`"

: ${wdir:?missing}

poly="${SOURCE_TEST_DIR}/test_filter_galois.d20.poly"
renumber="${SOURCE_TEST_DIR}/test_filter_galois.d20.renumber"
rels="${SOURCE_TEST_DIR}/test_filter_galois.d20.rels"
outrels="${wdir}/test_filter_galois.d20.rels"
args="-poly ${poly} -nrels 29 -renumber ${renumber} -galois _y -dl -large-ab\
       -outdir ${wdir} ${rels}"

SHA1BIN=sha1sum
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=sha1 ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=shasum ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then
    echo "Could not find a SHA-1 checksumming binary !" >&2
    exit 1
fi

if ! "${FG}" ${args}; then
  echo "$0: filter_galois binary failed. Files remain in ${wdir}"
  exit 1
fi

# XXX It used to be not strictly necessary to sort before computing the
# hash because filter_galois wasn't multi-threaded and the computation
# was deterministic. It's no longer the case.  Another source of
# harmless variation also exists with the ordering of the primes, which
# the implementation can choose to modify. Finally, it is important to
# use "sort" and not "sort -n", because the latter won't behave in a
# reproducible manner if passed integers that exceed the machine word
# size (at the very least, on alpine linux, our test_filter_galois
# testcase exhibits such a discrepancy).
sort_rels() {
    read -s -r -d '' perl_code <<- 'EOF'
        /^[^#]/ or next;
        chomp($_);
        my ($ab,@sides) = split(":", $_, 3);
        for (@sides) {
            $_ = join(",", sort({ hex($a) <=> hex($b) } split(",",$_)));
        }
        print join(":", ($ab, @sides)), "\n";
        
EOF
    perl -ne "$perl_code" "$@" | sort
}
SHA1=$(grep -v "^#" ${outrels} | sort_rels | ${SHA1BIN})
SHA1="${SHA1%% *}"

echo ""
echo "Got SHA1 of ${SHA1}"
echo "expected ${REFERENCE_SHA1}"
if [ "${SHA1}" != "${REFERENCE_SHA1}" ] ; then
  if [ "$CADO_DEBUG" ] ; then
      REFMSG=". Files remain in ${wdir}"
  else
      REFMSG=". Set CADO_DEBUG=1 to examine log output"
  fi
  echo "Got SHA1(${outrels})=${SHA1} but expected ${REFERENCE_SHA1}${REFMSG}"
  grep -v "^#" ${outrels}
  exit 1
fi
