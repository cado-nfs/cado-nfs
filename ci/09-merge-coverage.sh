#!/usr/bin/env bash

. "$(dirname $0)/000-functions.sh"
. "$(dirname $0)/001-environment.sh"

find . -name 'coverage-*'

rm -rf coverage || :
mkdir coverage

set -x

coverage_json_reports=(`ls coverage-${CI_COMMIT_SHORT_SHA}-*.json`)
# coverage_generated_sources=(`ls coverage-${CI_COMMIT_SHORT_SHA}-*generated-sources.tar.gz`)

gcovr_args=()
gcovr_args+=(--merge-mode-functions=separate)

for c in "${coverage_json_reports[@]}" ; do
    gcovr_args+=(-a "$c")
done

# Do we still need --gcov-ignore-parse-errors ?

commit=$CI_COMMIT_SHORT_SHA
commit_ref="$CI_PROJECT_URL/-/commit/$commit"
title="Coverage for commit <a href=\"$commit_ref\">$commit</a>"

echo "Generating html output"
gcovr "${gcovr_args[@]}"        \
    --html-title "$title"       \
    --html-nested coverage/coverage.html
ln -s coverage.html coverage/index.html

echo "Generating cobertura output"
gcovr "${gcovr_args[@]}"        \
    --cobertura coverage.xml    \
    --cobertura-pretty          \
    --print-summary

# do we want to add coverage.xml to the artifacts ?
tar czf coverage.tar.gz coverage
