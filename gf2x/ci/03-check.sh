#!/usr/bin/env bash

# This is actually __ONLY__ called from ci/40-testsuite.sh now

if ! [ "$sourced_from_testsuite" ] ; then
    . "${CI_PATH:-$(dirname $0)}/000-functions.sh"
    . "${CI_PATH:-$(dirname $0)}/001-environment.sh"
fi

NCPUS=`"${CI_PATH:-$(dirname $0)}/utilities/ncpus.sh"`
export NCPUS
export OMP_DYNAMIC=true STATS_PARSING_ERRORS_ARE_FATAL=1

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
    ${CI_PATH:-$(dirname $0)}/utilities/coverage_local_infofile_modifications.pl -d $build_tree ${C}-base.info
    set +x
    leave_section
fi

if [ "$CHECKS_EXPENSIVE" ] ; then
    enter_section "xtest" "Running expensive tests"
else
    enter_section "test" "Running tests"
fi
set +e
# --no-compress-output is perhaps better for test uploading, as ctest
# likes to store as zlib but headerless, which is a bit of a pain

if [ "$out_of_source" ] ; then
    set -o pipefail
    (cd "$build_tree" ; make -j$NCPUS check)
else
    "${MAKE}" -j$NCPUS check
fi
rc=$?
set -e
leave_section # test (or xtest)

if [ "$coverage" ] ; then
    enter_section "coverage" "extracting coverage data"
    gcovr --json . > ${C}-app.json
    set -x
    geninfo --ignore-errors gcov,source -q --output-filename ${C}-app.info -b . $build_tree --no-external
    ${CI_PATH:-$(dirname $0)}/utilities/coverage_local_infofile_modifications.pl -d $build_tree ${C}-app.info
    set +x
    leave_section
fi

STATUS=
if ! [ "$rc" = 0 ] ; then
    STATUS="FAILED "
fi

## TODO: process the test results, convert to junit. How do we do that?
# (below for ctest, as used in cado-nfs)
#enter_section "postprocessing" "Post-processing ${STATUS}ctest result into JUnit format"
#xmls=(`find "$build_tree/Testing" -name Test.xml`)
#if [ "${#xmls[@]}" = 1 ] ; then
#    xsltproc "${CI_PATH:-$(dirname $0)}/utilities/ctest-to-junit.xsl" "${xmls[0]}" > junit.xml
#    $ECHO_E "${CSI_BLUE}20 most expensive tests (real time):${CSI_RESET}"
#    perl -ne '/testcase.*" name="([^"]+)" time="([\d\.]+)"/ && print "$2 $1\n";' junit.xml  |sort -n | tail -n 20
#else
#    echo "Error, we expected one test xml file, we got ${#xmls[@]} (setting job as failed)" >&2
#    for f in "${xmls[@]}" ; do ls -l "$f" >&2 ; done
#    rc=1
#fi
#leave_section # postprocessing

exit $rc
