#!/usr/bin/env bash

# This installs cmake if it is not found in the current system

# minimal functionality is, it seems, provided by cmake 3.4+.

# I wouldn't be extremely surprised to realize that cmake 3.7+ is needed
# somehow, but I'm not sur.e

# cmake 3.9 or later is needed in order to properly see the
# "skipped tests"

name=cmake
version=3.4.3
package=${name}-${version}.tar.gz
url=http://www.cmake.org/files/v3.4/${package}

prefix="$1"
shift

echo "installing ${name} ${version} in ${prefix}"

tmpdir=`mktemp -d ${TMPDIR-/tmp}/cado-nfs.${name}-build.XXXXXXXX`
cd $tmpdir
rm -f ${package}
wget="`which wget 2>/dev/null`"
if [ "$?" = "0" ] ; then
    wget --no-check-certificate $url
else
    curl="`which curl 2>/dev/null`"
    if [ "$?" = "0" ] ; then
        curl -kL $url > $package
    else
        echo "Need either wget or curl to get $url" >&2
        exit 1
    fi
fi
tar xzf ${package}
# Better unset these, since we expect cmake to do the right thing without
# our custom flags !
unset CC
unset CXX
unset CFLAGS
unset CXXFLAGS
cd ${name}-${version}
./configure --prefix=$prefix
make -j 4
make install
rm -rf $tmpdir
