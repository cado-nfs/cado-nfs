#!/usr/bin/env bash

set -e

if [ "$CADO_DEBUG" ] ; then set -x ; fi

dup1="${PROJECT_BINARY_DIR?missing}/filter/dup1"
c60="${CADO_NFS_SOURCE_DIR?missing}/tests/misc/c60.rels1"

# try when we have only one slice (it should just be a copy!)
c60tmp="${wdir?missing}/c60"
mkdir "$c60tmp"
mkdir "$c60tmp"/0
$dup1 -n 0 --out "$c60tmp" --prefix $(basename "$c60") "$c60"
diff -q "$c60" "$c60tmp"/0/$(basename "$c60").0000
rm -rf $c60tmp

# two slices
c60tmp="${wdir?missing}/c60"
mkdir "$c60tmp"
mkdir "$c60tmp"/{0,1}
$dup1 -n 1 --out "$c60tmp" --prefix $(basename "$c60") "$c60"
diff -q <(sort -u <"$c60") <(cat "$c60tmp"/*/$(basename "$c60").* | sort -u)
rm -rf $c60tmp

# two slices, many files
c60tmp="${wdir?missing}/c60"
mkdir "$c60tmp"
mkdir "$c60tmp"/{0,1}
$dup1 -n 1 --lognrels 1 --out "$c60tmp" --prefix $(basename "$c60") "$c60"
diff -q <(sort -u <"$c60") <(cat "$c60tmp"/*/$(basename "$c60").* | sort -u)
rm -rf $c60tmp

# single-slice mode, but do it several times and make sure we lose nothing
c60tmp="${wdir?missing}/c60"
mkdir "$c60tmp"
mkdir "$c60tmp"/{0,1,2,3}
for i in 0 1 2 3 ; do
$dup1 -n 2 --only $i --lognrels 2 --out "$c60tmp" --prefix $(basename "$c60") "$c60"
done
diff -q <(sort -u <"$c60") <(cat "$c60tmp"/*/$(basename "$c60").* | sort -u)

# now test filelist mode, based on the output of the previous test
ls "$c60tmp"/*/$(basename "$c60").* > "$c60tmp/filelist.txt"

# one slice
c60tmp2="${wdir?missing}/c60b"
mkdir "$c60tmp2"
mkdir "$c60tmp2"/0
$dup1 -n 0 --out "$c60tmp2" --prefix $(basename "$c60") --filelist "$c60tmp/filelist.txt"
diff -q <(sort -u <"$c60") <(cat "$c60tmp2"/*/$(basename "$c60").* | sort -u)
rm -rf "$c60tmp2"

# two slices
c60tmp2="${wdir?missing}/c60b"
mkdir "$c60tmp2"
mkdir "$c60tmp2"/{0,1}
$dup1 -n 1 --out "$c60tmp2" --prefix $(basename "$c60") --filelist "$c60tmp/filelist.txt"
diff -q <(sort -u <"$c60") <(cat "$c60tmp2"/*/$(basename "$c60").* | sort -u)
rm -rf "$c60tmp2"

rm -rf $c60tmp


