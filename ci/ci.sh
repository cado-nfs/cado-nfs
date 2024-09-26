
export project_name=cado-nfs
: ${CI_PROJECT_URL=https://gitlab.inria.fr/cado-nfs/cado-nfs}
export build_system=cmake

needs_bc=1
needs_python=1
needs_perl=1
needs_optional_hwloc=1
needs_optional_ecm=1
needs_optional_fmt=1
needs_gmp=1

tweak_tree_before_configure() { : ; }

step_configure() {
    if [ "$using_cmake_directly" ] ; then
        (cd "$build_tree" ; cmake "$source_tree")
    else
        "${MAKE}" cmake
    fi
    eval `${MAKE} show`
}

build_steps="build1 build2"
build_step_name_build1="Building"

step_build1() {
    if [ "$using_cmake_directly" ] ; then
        SOURCEDIR="$PWD"
        (cd "$build_tree" ; "${MAKE}" -j$NCPUS)
    else
        "${MAKE}" -j$NCPUS
    fi
}

build_step_name_build2="Building test dependencies"
step_build2() {
    if [ "$using_cmake_directly" ] ; then
        SOURCEDIR="$PWD"
        (cd "$build_tree" ; "${MAKE}" -j$NCPUS all_test_dependencies)
    else
        "${MAKE}" -j$NCPUS all_test_dependencies
    fi
}

check_environment() {
    export OMP_DYNAMIC=true
    export STATS_PARSING_ERRORS_ARE_FATAL=1
}


step_coverage() {
    # This takes a coverage file prefix as $1, and info for the kind of
    # file (whether it's "base" or "app", for instance) in $2.
   
    # build_tree and source_tree are used

    outfile="$source_tree/$1-$2.json"

    # Rationale for having this "base" coverage step (and I'm not sure it
    # is still relevant):
    # The "base" coverage file has zero coverage for every instrumented
    # line of the project. At a later stage, we will combine this data
    # file with coverage data files captured after the test run. This way
    # the percentage of total lines covered will always be correct, even
    # when not all source code files were loaded during the test(s).

    # might be useful. We don't want to bother with traces of config
    # checks.
    find "$build_tree" -name '*conftest.gcno' -o -name 'CMake*.gcno' -o -name '?-CMake*.gcno' | xargs -r rm -v

    # ci/ci/001-environment.sh sets build_tree to "./generated" for coverage
    # jobs. Therefore, all files, gcno and gcda, are found under
    # $PWD==$src_tree .
    # This is done so because we have to ship the files that are
    # generated by the build process, and expose them to the merged
    # coverage report.
    # If $build_tree is outside $src_tree, we should add "." before
    # --json
    (set -x ; cd "$build_tree" ; time gcovr --merge-mode-functions=separate -r "$source_tree" --json "$outfile")
}

step_coverage_more_artifacts() {
    prefix="$1"
    if [ "$build_tree" != generated ] ; then
        echo "This part of the script assumes that build_tree=generated"
    fi

    # because of /bin/sh, we can't do arrays.
    find "$build_tree" -name '*.[ch]' -o -name '*.[ch]pp' | xargs -x tar czf ${prefix}-generated-sources.tar.gz
}


coverage_expunge_paths="utils/embedded:gf2x:generated:linalg/bwc/flint-fft:linalg/bwc/mpfq"

step_check() {
    # --no-compress-output is perhaps better for test uploading, as ctest
    # likes to store as zlib but headerless, which is a bit of a pain

    ctest_args="-T Test --no-compress-output --test-output-size-passed 4096 --test-output-size-failed 262144"

    if [ "$using_cmake_directly" ] ; then
        set -o pipefail
        (cd "$build_tree" ; ctest -j$NCPUS $ctest_args ) | "$source_tree"/scripts/filter-ctest.pl
    else
        "${MAKE}" check ARGS="-j$NCPUS $ctest_args"
    fi
}

step_doc() { : ; }
