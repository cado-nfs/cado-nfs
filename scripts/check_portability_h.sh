#!/bin/bash

# portability.h is necessary for a range of functions that we happen to
# use.

git ls-files '*.[ch]' '*.[ch]pp' | \
    xargs grep -lw '\(strdup\|strndup\|strlcat\|strlcpy\|asprintf\|lrand48\|realpath\|pagesize\|sleep\)' | \
    grep -v '^cado\.h$' | \
    grep -v '^portability\.h$' | \
    grep -v '^config/' | \
    grep -v '^gf2x/' | \
    grep -v '^linalg/bwc/flint-fft/' | \
    (rc=0 ; while read x ; do 
        if ! grep -q portability\.h "$x" ; then
            if [ $rc = 0 ] ; then
                echo "# The following files use a function whose prototype is, in certain cases, provided by portability.h ; you must #include it." >&2
                rc=1
            fi
            echo "$x" >&2
        fi
    done ; exit $rc)
