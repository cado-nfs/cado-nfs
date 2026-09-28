#!/usr/bin/env bash

# Driver for tests/sieve/test-unsieve-patterns, which decides whether the
# survivor search should use its bound patterns for the primes 3 and 5. The
# answer is a compile-time constant in sieve/las-unsieve-sse2.cpp, and it
# depends on the microarchitecture, so it has to be measured per target.
#
# Usage:
#   tests/sieve/tune-unsieve-patterns.sh [build tree] [-- extra tuner args]
#
# With no build tree, the one that local.sh would pick is used, and the tuner
# is built if need be. Note that the tuner is only configured when
# CHECKS_EXPENSIVE is set, which this script therefore does.
#
# It sweeps the line length and the survivor density around the values las
# really works at, and prints the setting to compile with. A sweep that does
# not agree with itself is reported as such rather than averaged away: it
# means the effect is too small to act on for this target.

set -eu

export CHECKS_EXPENSIVE=1

srcdir="$(cd "$(dirname "$0")/../.." && pwd)"
build=
args=()
while [ $# -gt 0 ] ; do
    case "$1" in
        --) shift ; args=("$@") ; break ;;
        -*) args+=("$1") ;;
        *)  build="$1" ;;
    esac
    shift
done

if [ -z "$build" ] ; then
    build="$(cd "$srcdir" && "$srcdir"/scripts/build_environment.sh --show \
        | sed -n 's/^build_tree=//p' | tr -d '"')"
fi
if [ -z "$build" ] || ! [ -d "$build" ] ; then
    echo "no build tree found; pass one as an argument" >&2
    exit 1
fi

tuner="$build/tests/sieve/test-unsieve-patterns"
if ! [ -x "$tuner" ] ; then
    echo "# building the tuner in $build" >&2
    (cd "$srcdir" && force_build_tree="$build" make -j "$(nproc)" \
        test-unsieve-patterns) >&2
fi

echo "# $(uname -m), $(sed -n 's/^model name[^:]*: //p' /proc/cpuinfo | head -1)"

verdicts=()
for I in 14 15 16 ; do
    for density in 3e-5 1e-4 2.5e-4 ; do
        echo "## -I $I -density $density"
        out="$("$tuner" -I "$I" -density "$density" -report-only "${args[@]}")"
        echo "$out" | sed -n '/configuration/,/nanoseconds/p'
        verdicts+=("$(echo "$out" | sed -n 's/^UNSIEVE_PATTERNS=//p')")
    done
done

uniq_verdicts="$(printf '%s\n' "${verdicts[@]}" | sort -u | tr '\n' ' ')"
echo "#"
case "$uniq_verdicts" in
    "0 ") echo "# verdict: compile with -DUNSIEVE_PATTERNS=0 (no bound patterns)" ;;
    "1 ") echo "# verdict: compile with -DUNSIEVE_PATTERNS=1 (bound patterns)" ;;
    *)   echo "# verdict: none. The sweep does not agree with itself"
         echo "# ($uniq_verdicts), so the effect is too small to act on here." ;;
esac
echo "# The default for this target is in sieve/las-unsieve-sse2.cpp."
