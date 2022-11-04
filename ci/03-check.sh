#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

NCPUS=`"$(dirname $0)/utilities/ncpus.sh"`
export NCPUS
export OMP_DYNAMIC=true STATS_PARSING_ERRORS_ARE_FATAL=1

if ! [ "$using_cmake_directly" ] ; then
    eval $("${MAKE}" show)
fi

if [ "$coverage" ] ; then
    # The "base" coverage file has zero coverage for every instrumented
    # line of the project. At a later stage, we will combine this data
    # file with coverage data files captured after the test run. This way
    # the percentage of total lines covered will always be correct, even
    # when not all source code files were loaded during the test(s).
    enter_section "coverage" "preparing base coverage data"

    # might be useful.
    find $build_tree -name '*conftest.gcno' | xargs -r rm
    find $build_tree -name 'CMake*.gcno' | xargs -r rm

    C=coverage-$CI_COMMIT_SHORT_SHA-$CI_JOB_ID
    set -x
    # avoid -b . --no-external ; -b is also used to rewrite symlinks,
    # which we clearly don't want, and without it, --no-external is
    # incapable of properly filtering file names.
    # lcov -q -c -i -d $build_tree -b . -o ${C}-base.info --no-external
    # These two are equivalent.
    lcov -q -c -i -d $build_tree -o ${C}-base.info
    # geninfo -q -i $build_tree -o ${C}-base.info
    $(dirname $0)/utilities/coverage_local_infofile_modifications.pl -d $build_tree ${C}-base.info
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

ctest_args=(
    -T Test
    --no-compress-output
    --test-output-size-passed 4096
    --test-output-size-failed 262144
)

if [ "$using_cmake_directly" ] ; then
    set -o pipefail
    (cd "$build_tree" ; ctest -j$NCPUS "${ctest_args[@]}") | "$source_tree"/scripts/filter-ctest.pl
else
    "${MAKE}" check ARGS="-j$NCPUS ${ctest_args[*]}"
fi
rc=$?
set -e
leave_section # test (or xtest)

if [ "$coverage" ] ; then
    enter_section "coverage" "extracting coverage data"

    # gcovr is terribly picky. There are countless ways to make it lose
    # track of symlinked sources in an out-of-source build.
    # Both versions below work in a fresh checkout with $build_tree
    # somewhere below ./ ; the latter is simpler.
    # gcovr --json ${C}-app.json ./ $build_tree/
    # gcovr --json ${C}-app.json ./
    # when $build_tree is somewhere else, it seems that the following is
    # a more robust way to proceed.
    # TODO: I doubt that gcovr correctly tracks the symlinks in the build
    # tree.
    gcovr --gcov-ignore-parse-errors --json ${C}-app.json $build_tree/ -f . -f $build_tree
    set -x
    # It _seems_ that in fact, we do **NOT** want --no-external, and -b
    # is actually doing more harm than good.
    # geninfo --ignore-errors gcov,source -q --output-filename ${C}-app.info -b . $build_tree --no-external
    # These two are equivalent
    lcov -q -c -d $build_tree -o ${C}-app.info
    # geninfo --ignore-errors gcov,source -q --output-filename ${C}-app.info $build_tree
    $(dirname $0)/utilities/coverage_local_infofile_modifications.pl -d $build_tree ${C}-app.info
    set +x
    # well, no. Let's rather rewrite the references to the build tree
    # sources in the info file directly.
    # tar czf ${C}-generated-sources.tar.gz $(perl -ne "m,^SF:${build_tree#$PWD/}/, && s,^SF:,, && print;" ${C}-base.info  | sort -u)
    leave_section
fi

STATUS=
if ! [ "$rc" = 0 ] ; then
    STATUS="FAILED "
fi
enter_section "postprocessing" "Post-processing ${STATUS}ctest result into JUnit format"
xmls=(`find "$build_tree/Testing" -name Test.xml`)
if [ "${#xmls[@]}" = 1 ] ; then
    xsltproc "$(dirname $0)/utilities/ctest-to-junit.xsl" "${xmls[0]}" > junit.xml
    $ECHO_E "${CSI_BLUE}20 most expensive tests (real time):${CSI_RESET}"
    perl -ne '/testcase.*" name="([^"]+)" time="([\d\.]+)"/ && print "$2 $1\n";' junit.xml  |sort -n | tail -n 20
else
    echo "Error, we expected one test xml file, we got ${#xmls[@]} (setting job as failed)" >&2
    for f in "${xmls[@]}" ; do ls -l "$f" >&2 ; done
    rc=1
fi
leave_section # postprocessing

exit $rc
