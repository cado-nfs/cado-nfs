#!/usr/bin/env bash

# Run the test suite against a coverage build, and render a report.
#
# See dev_docs/README.coverage.md. In short:
#
#     COV=1 make -j8 all all_test_dependencies
#     COV=1 scripts/cov.sh -R numbertheory
#
# Any argument is passed on to ctest, so -R / -E / -j and friends work as
# usual. With no argument at all, the whole test suite runs, which takes
# a while under coverage instrumentation.
#
# This uses the very same code as the ci does (ci/ci.bash), so what you
# see here is what the pipeline will show.

set -e

cd "$(dirname "$0")/.."

: ${COV=1}
export COV

source_tree="$PWD"
export source_tree

# local.sh appends ".cov" to the build tree when COV is set
eval "$(make show)"
export build_tree

if ! [ -d "$build_tree" ] ; then
    echo "No $build_tree. Run: COV=1 make -j8 all all_test_dependencies" >&2
    exit 1
fi

for prog in grcov genhtml ; do
    if ! type -p "$prog" > /dev/null ; then
        echo "$prog is not in \$PATH." >&2
        case "$prog" in
            grcov) echo "  It is a single static binary, see" >&2
                   echo "  https://github.com/mozilla/grcov/releases" >&2
                   echo "  (debian and friends also have a grcov package)" >&2;;
            genhtml) echo "  It comes with lcov (apt install lcov)." >&2;;
        esac
        exit 1
    fi
done

NCPUS="$(nproc 2> /dev/null || echo 4)"
export NCPUS
COMMIT_SHORT_SHA="$(git rev-parse --short=8 HEAD 2> /dev/null || echo unknown)"
PROJECT_URL="https://gitlab.inria.fr/cado-nfs/cado-nfs"

# ci/ci.bash is only function and variable definitions, so sourcing it is
# harmless, and it saves us from keeping a second copy of the exclusion
# list and of the report recipe.
enter_section() { shift 2 ; echo "### $*" ; }
leave_section() { : ; }
yellow_message() { echo "$@" >&2 ; }
major_message() { echo "$@" ; }
fatal_error() { echo "$@" >&2 ; exit 1 ; }
. ./ci/ci.bash

# Keep the report out of the checkout: everything below the build tree is
# throwaway anyway.
coverage_report_workdir="$build_tree/cov-work"
coverage_report_dir="$build_tree/coverage"
coverage_report_xml="$build_tree/coverage.xml"

# Start from a clean slate, so that what we measure is this run and not
# an accumulation of the previous ones.
find "$build_tree" -name '*.gcda' -delete
purge_unused_coverage_files

rc=0
OMP_DYNAMIC=true "${MAKE:-make}" check ARGS="-E ^builddep $*" || rc=$?

coverage_tree="$build_tree"
step_coverage_report

index="$coverage_report_dir/index.html"
echo
echo "Coverage report: file://$index"
if type -p gio > /dev/null ; then
    gio open "$index" || :
fi

exit $rc
