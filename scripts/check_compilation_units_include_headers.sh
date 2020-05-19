#!/bin/bash

# check that .c and .cpp files follow good practice of always including
# the header file of the same name. This could avoid terrible errors with
# inconsistent prototype and implementation.
git ls-files '*.c' '*.cpp' | \
    grep -v linalg/bwc/flint-fft | \
    (rc=0; while read x ; do
    x0="$x"
    x0=${x0%%.cpp}
    x0=${x0%%.c} 
    b=`basename $x0`
    if [ -e "$x0.h" ] || [ -e "$x0.hpp" ] ; then
        if ! grep -q '"'"$b"'.h\(pp\)\?"' "$x" ; then
            if [ $rc = 0 ] ; then
                echo "# The following compilation units should include the header file of the same name." >&2
                rc=1
            fi
            echo "$x" >&2
        fi
    fi
done ; exit $rc)

