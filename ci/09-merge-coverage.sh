#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

find . -name 'coverage-*.info'

set -x

coverage_info_reports=(`ls coverage-${CI_COMMIT_SHORT_SHA}-*.info`)
coverage_json_reports=(`ls coverage-${CI_COMMIT_SHORT_SHA}-*.json`)
# coverage_generated_sources=(`ls coverage-${CI_COMMIT_SHORT_SHA}-*generated-sources.tar.gz`)

tracefiles=()
for c in "${coverage_json_reports[@]}" ; do
    tracefiles+=(--add-tracefile "$c")
done

# for f in "${coverage_generated_sources[@]}" ; do
#     tar xf "$f"
# done

gcovr --gcov-ignore-parse-errors "${tracefiles[@]}" --xml coverage.xml --xml-pretty -s

$(dirname $0)/utilities/coverage-postprocess.sh -o coverage "${coverage_info_reports[@]}"

tar czf coverage.tar.gz coverage

