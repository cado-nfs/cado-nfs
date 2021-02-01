#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

NCPUS=`"$(dirname $0)/utilities/ncpus.sh"`
export NCPUS
if [ "$CHECKS_EXPENSIVE" ] ; then
    enter_section "xtest" "Running expensive tests"
else
    enter_section "test" "Running tests"
fi
export OMP_DYNAMIC=true STATS_PARSING_ERRORS_ARE_FATAL=1
eval $(make show)

if [ "$coverage" ] ; then
    # The "base" coverage file has zero coverage for every instrumented
    # line of the project. At a later stage, we will combine this data
    # file with coverage data files captured after the test run. This way
    # the percentage of total lines covered will always be correct, even
    # when not all source code files were loaded during the test(s).
    enter_section "coverage" "preparing base coverage data"
    C=coverage-$CI_COMMIT_SHORT_SHA-$CI_JOB_ID
    set -x
    lcov -q -c -i -d $build_tree -b . -o ${C}-base.info --no-external
    $(dirname $0)/utilities/coverage_local_infofile_modifications.pl -d $build_tree ${C}-base.info
    set +x
    leave_section
fi

set +e
# --no-compress-output is perhaps better for test uploading, as ctest
# likes to store as zlib but headerless, which is a bit of a pain
make check ARGS="-j$NCPUS -T Test --no-compress-output --test-output-size-passed 4096 --test-output-size-failed 262144"
rc=$?
set -e

if [ "$coverage" ] ; then
    enter_section "coverage" "extracting coverage data"
    gcovr --json . > ${C}-app.json
    set -x
    geninfo --ignore-errors gcov,source -q --output-filename ${C}-app.info -b . $build_tree --no-external
    $(dirname $0)/utilities/coverage_local_infofile_modifications.pl -d $build_tree ${C}-app.info
    set +x
    leave_section
fi

STATUS=
if ! [ "$rc" = 0 ] ; then
    STATUS="FAILED "
fi
enter_section "postprocessing" "Post-processing ${STATUS}ctest result into JUnit format"
xmls=(`find "$build_tree/Testing" -name Test.xml`)
if [ "${#xmls[@]}" != 1 ] ; then
    echo "Error, we expected one test xml file, we got ${#xmls[@]}" >&2
    for f in "${xmls[@]}" ; do ls -l "$f" >&2 ; done
    exit 1
fi
xsltproc "$(dirname $0)/utilities/ctest-to-junit.xsl" "${xmls[0]}" > junit.xml
$ECHO_E "${CSI_BLUE}20 most expensive tests (real time):${CSI_RESET}"
perl -ne '/testcase.*" name="([^"]+)" time="([\d\.]+)"/ && print "$2 $1\n";' junit.xml  |sort -n | tail -n 20
leave_section

exit $rc
