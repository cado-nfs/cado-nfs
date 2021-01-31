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

set +e
make check ARGS="-j$NCPUS -T Test"
rc=$?
leave_section

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
xsltproc "$(dirname $0)/ctest-to-junit.xsl" "${xmls[0]}" > junit.xml
leave_section
exit $rc
