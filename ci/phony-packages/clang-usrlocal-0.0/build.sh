#!/bin/bash

if ! [ -d "$1" ] ; then
    echo "Need directory to place the .deb file" >&2
    exit 1
fi

T=$(mktemp -d /tmp/XXXXXXXXX)
trap "rm -rf $T" EXIT
F=clang-usrlocal_0.0_all.deb

echo 2.0 > $T/debian-binary
touch $T/md5sums
(cd $(dirname $0) ; cp -p control postinst prerm $T/)
(cd $T ; tar czf control.tar.gz control postinst prerm md5sums)
tar czf $T/data.tar.gz --files-from=/dev/null
(cd $T ; ar q $F  debian-binary control.tar.gz data.tar.gz)
cp "$T/$F" "$1/$F"
echo "$1/$F"
