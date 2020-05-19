#!/bin/bash

# all headers must have some sort of ifndef/define/endif construct, and C
# headers must have the typical extern "C" construct as well. (unless
# they have no prototypes)

has_double_include_protection() {
    # make it a very simple check
    grep -q '\(^#ifndef\|pragma multi include\)' "$1"
    # && grep -q '^#define' "$1" && grep -q '^#endif' "$1"
}

has_cxx_prototype_protection() {
    if [[ "$1" =~ \.hpp$ ]] ; then
        true
    else
        # assume that automatically generated file care about mentioning
        # that they've been generated automatically, and possibly less
        # about sticking in a pragma that we like.
        grep -q '\(__cplusplus\|pragma no prototypes\|generated automatically\|Automatically generated\)' "$1"
    fi
}

is_ok() {
    has_double_include_protection "$1" && has_cxx_prototype_protection "$1"
}


git ls-files '*.h' '*.hpp' | \
    grep -v '^cado\.h$' | \
    grep -v '^gf2x/' | \
    grep -v '^config/' | \
    grep -v '^utils/fmt/' | \
    grep -v '^linalg/bwc/flint-fft/' | \
    grep -v '^linalg/bwc/mpfq/' | \
    (rc=0 ; while read x ; do 
        if ! is_ok "$x" ; then
            if [ $rc = 0 ] ; then
                echo "# The following header files miss some of the mandatory protections (#ifndef guards, c++ guards; please add them)." >&2
                rc=1
            fi
            echo "$x" >&2
        fi
    done ; exit $rc)
