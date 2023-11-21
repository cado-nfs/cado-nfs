#!/bin/sh


# Detect Musl C library.
libc=$(ldd /bin/ls 2>/dev/null | grep 'musl' | head -1 | cut -d ' ' -f1)
if [ -z $libc ] ; then
    # This is not Musl.
    exit 1
fi
T=$(mktemp -d /tmp/XXXXXX)
trap "rm -rf $T" EXIT
tmpf=$T/log
$libc >${tmpf} 2>&1
vstr=$(cat ${tmpf} | grep "Version" | cut -d ' ' -f2)

v_major=$(echo $vstr | cut -d '.' -f1)
v_minor=$(echo $vstr | cut -d '.' -f2)
v_patch=$(echo $vstr | cut -d '.' -f3)

rm -f ${tmpf}

printf "0x%02x%03x%03x\n" $v_major $v_minor $v_patch
