#!/usr/bin/env bash

set -e

if [ "$NCPUS_FAKE" ] ; then
    echo $NCPUS_FAKE
    exit 0
fi

if [ -f /proc/cpuinfo ] ; then
    # grep -c fails with no match, while |wc -l is happy
    nphysical=$(sort -u /proc/cpuinfo  | grep '^physical' | wc -l)
    if [ "$nphysical" -eq 0 ] ; then
        grep -c ^processor /proc/cpuinfo
        exit 0
    fi
    variants="$(grep '^cpu cores' /proc/cpuinfo | uniq | wc -l)"
    if [ "$variants" = 0 ] ; then
        grep -c ^processor /proc/cpuinfo
        exit 0
    fi
    if [ "$variants" != 1 ] ; then
        echo "inhomogeneous platform ?" >&2
        exit 1
    fi
    cores_per_cpu=$(grep '^cpu cores' /proc/cpuinfo | head -1 | cut -d: -f2)
    if [ $((nphysical*cores_per_cpu)) -gt 1 ] ; then
        # make sure we're not restricted by some stubborn job manager
        allowed="$(awk '/^Cpus_allowed:/ { print $2; }' /proc/self/status)"
        if ! [ "$allowed" ] ; then
            # we do not want to bother.
            exit 0
        fi
        if [[ $allowed =~ ^[0,]*[1248,][0,]*$ ]] ; then
            echo "Uh oh. We're bound to only 1 cpu out of $((nphysical*cores_per_cpu)), that does not seem right" >&2
            exit 1
        fi
    fi
    echo $((nphysical*cores_per_cpu))
elif [ "$(uname -s)" = Darwin ] ; then
    # does this count hyperthreading or not ?
    /usr/sbin/sysctl -n hw.ncpu
elif [ "$(uname -s)" = OpenBSD ] ; then
    # does this count hyperthreading or not ?
    sysctl -n hw.ncpu
elif [ "$(uname -s)" = FreeBSD ] ; then
    # does this count hyperthreading or not ?
    sysctl -n hw.ncpu
elif [ "$(uname -s)" = MINGW32_NT-6.1 ] ; then
    # not clear whether it's physical or logical.
    # wmic cpu get Caption | tail -n +2 | grep -c .
    # anyway we don't believe mingw will support multithreading well. For
    # sure make -j2 seems to just not work...
    echo 1
else
    # this would work as well on linux and darwin, but pretty surely does
    # not count hyperthreading. Does not work on openbsd5.3
    getconf _NPROCESSORS_ONLN || getconf NPROCESSORS_ONLN
fi

