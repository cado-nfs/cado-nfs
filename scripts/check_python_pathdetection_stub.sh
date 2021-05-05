#!/usr/bin/env bash

if ! [ -f cado.h ] ; then
    echo "Please call $0 from the top of the source tree" >&2
    exit 1
fi

getstub() {
    perl -ne '$x=1 if /THIS PART MUST BE EXACTLY IDENTICAL/; print if $x; $x=0 if /END OF THE PART THAT MUST BE EXACTLY IDENTICAL/;'  "$@"
}

if ! diff -u <(getstub cado-nfs.py) <(getstub cado-nfs-client.py) ; then
    echo
    echo "### Please fix the differences above in cado-nfs.py and cado-nfs-client.py before committing" >&2
    exit 1
fi

