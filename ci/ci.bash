
tweak_tree_before_configure() { : ; }

user_variables() {
    if [ "$using_cmake_directly" ] ; then
        # use build_tree in this case, which matches the variable that
        # call_cmake.sh uses, by the way.
        if [ "$BASH_VERSION" ] ; then
            if ! [ "$build_tree" ] ; then
                build_tree="/tmp/$CI_JOB_NAME"
                # spaces in dir names don't work, mostly because of libtool
                # (look at gf2x/fft/libgf2x-fft.la)
                # This substitution is bash-only, but this should be fine to 
                # have in a conditional that non-bash skips over
                build_tree="${build_tree// /_}"
                export build_tree
            fi
            if ! [ -d "$build_tree" ] ; then
                mkdir -p "$build_tree"
            fi
        else
            # just a safeguard
            build_tree=/no/build_tree/set/because/we/require/bash/for/that
            export build_tree
        fi
    fi
}


step_configure() {
    # now that we're confident that we've made the bwc checks specific to
    # a "with_sagemath" suffix, there's no risk in missing the sagemath
    # code by inadvertence.
    # if [ "$specific_checks" = "bwc.sagemath" ] ; then
    #     export FORCE_BWC_EXTERNAL_CHECKS_OUTPUT_ON_FD3=1
    # fi
    if [ "$specific_checks" = "including_mpi" ] ; then
        export MPI=1
        # sigh. when we run in containers, running as root isn't much of
        # a problem
        export OMPI_ALLOW_RUN_AS_ROOT=1
        export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
    elif [ "$specific_checks" = "only_mpi" ] ; then
        export MPI=1
        # sigh. when we run in containers, running as root isn't much of
        # a problem
        export OMPI_ALLOW_RUN_AS_ROOT=1
        export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
    fi
    if [ "$using_cmake_directly" ] ; then
        (cd "$build_tree" ; cmake "$source_tree" $pass_flags_to_cmake)
        # Ignore local.sh if we're building directly from cmake
    else
        "${MAKE}" cmake
        eval `${MAKE} show`
    fi
}

build_steps="build1 build2"
build_step_name_build1="Building"

step_build1() {
    target=all
    if [ "$specific_checks" = "bwc.sagemath" ] ; then
        target=all_sagemath_test_dependencies
    fi
    if [ "$using_cmake_directly" ] ; then
        SOURCEDIR="$PWD"
        (cd "$build_tree" ; "${MAKE}" -j$NCPUS $target)
    else
        "${MAKE}" -j$NCPUS $target
    fi
}

build_step_name_build2="Building test dependencies"
step_build2() {
    target=all_test_dependencies
    if [ "$specific_checks" = "bwc.sagemath" ] ; then
        # already covered in build1 anyway
        return
    fi
    if [ "$using_cmake_directly" ] ; then
        SOURCEDIR="$PWD"
        (cd "$build_tree" ; "${MAKE}" -j$NCPUS $target)
    else
        "${MAKE}" -j$NCPUS $target
    fi
}

prepare_valgrind_environment() {
    vdir=$PWD/valgrind.$CI_COMMIT_SHORT_SHA-$CI_JOB_ID
    export vdir
    mkdir -p $vdir
    cat > $vdir/v.sh <<EOF
#!/usr/bin/env bash
cado=$PWD
vdir=$vdir
EOF
    cat >> $vdir/v.sh <<'EOF'
if [ -x "./$1" ] ; then
    prg="$1"
    shift
    args=("./$prg" "$@")
    set -- "${args[@]}"
fi

# Don't use --error-exitcode, so that we get a chance to be notified of all potential errors at once.
valgrind --suppressions=$cado/cado-nfs.supp --gen-suppressions=all --trace-children=yes --trace-children-skip=gdb,gzip,libtool,gcc,g++ "--log-file=$vdir/pid-%p" --leak-check=full "$@"
EOF

    VALGRIND="$vdir/v.sh"
    chmod 755 $VALGRIND

    export PYTHONDONTWRITEBYTECODE=1
    test_precommand+=(env TEST_PRECOMMAND=$VALGRIND)
    # valgrind tests can take _ages_ if we run them with openmp
    export OMP_NUM_THREADS=1
    # on the other hand we do want at least 2 threads for things that
    # touch mf_scan2
    export CADO_NFS_MAX_THREADS=2
}


check_environment() {
    Nmax=16
    if [ -x "$build_tree/tests/omp_get_max_threads" ] ; then
        N=$("$build_tree/tests/omp_get_max_threads")
        # originally we sensed a need to do so only on 32-bit machines,
        # but after all it makes sense more generally.
        if [ "$N" -gt "$Nmax" ] ; then
            major_message "reducing the max number of openmp threads to only $Nmax"
            export OMP_NUM_THREADS=$Nmax
            export OMP_THREAD_LIMIT=$Nmax
        fi
    fi
    if [ -f "$build_tree/hwloc-`hostname`.xml" ] ; then
        export HWLOC_XMLFILE="$build_tree/hwloc-`hostname`.xml"
    elif [ -x "$build_tree/tests/hwloc_cado_helper" ] ; then
        export HWLOC_XMLFILE="$build_tree/hwloc-`hostname`.xml"
        "$build_tree/tests/hwloc_cado_helper" -o "$HWLOC_XMLFILE"
    else
        major_message "Forcing a fake hwloc file. This might conflict with some tests"
        export HWLOC_XMLFILE="$PWD/ci/placeholder-machine-for-tests.xml"
    fi
    export CADO_NFS_MAX_THREADS=$Nmax
    export OMP_DYNAMIC=true
    # See https://stackoverflow.com/questions/70126350/openmp-incredibly-slow-when-another-process-is-running
    # It's not totally clear to me if it somewhere specified that
    # lowercase "passive" implies GOMP_SPINCOUNT=0 for gcc. If it's not
    # specified, it may change in the future, so let's force the setting
    # ourselves.
    export OMP_DISPLAY_ENV=verbose
    export OMP_WAIT_POLICY=passive
    export GOMP_SPINCOUNT=0
    # OMP_PROC_BIND helps in certain cases, and is a disaster in other
    # cases. We can't afford it.
    # export OMP_PROC_BIND=true
    export STATS_PARSING_ERRORS_ARE_FATAL=1

    test_precommand=()
    if [ "$valgrind" ] ; then
        prepare_valgrind_environment
    fi
    
}

purge_unused_coverage_files() {
    # might be useful. We don't want to bother with traces of config
    # checks.
    (
        find "$build_tree" -name '*conftest*gcno' -o -name 'CMake*.gcno' -o -name '?-CMake*.gcno'
        find "$build_tree" -name '*conftest*gcda' -o -name 'CMake*.gcda' -o -name '?-CMake*.gcda'
        find "$build_tree"/utils/embedded -name \*.gcda -o -name \*.gcno 2>/dev/null || :
        find "$build_tree"/gf2x -name \*.gcda -o -name \*.gcno 2>/dev/null || :
        find "$build_tree"/linalg/bwc/flint-fft -name \*.gcda -o -name \*.gcno 2>/dev/null || :
    ) | sort -u | xargs -r rm -v
}


step_coverage() {
    # This takes a coverage file prefix as $1, and info for the kind of
    # file (whether it's "base" or "app", for instance) in $2.
   
    # build_tree and source_tree are used


    # Rationale for having this "base" coverage step (and I'm not sure it
    # is still relevant):
    # The "base" coverage file has zero coverage for every instrumented
    # line of the project. At a later stage, we will combine this data
    # file with coverage data files captured after the test run. This way
    # the percentage of total lines covered will always be correct, even
    # when not all source code files were loaded during the test(s).

    if true ; then
        purge_unused_coverage_files

        outfile="$source_tree/$1-$2.json"
        # ci/ci/001-environment.sh sets build_tree to "./generated" for coverage
        # jobs. Therefore, all files, gcno and gcda, are found under
        # $PWD==$src_tree .
        # This is done so because we have to ship the files that are
        # generated by the build process, and expose them to the merged
        # coverage report.
        # If $build_tree is outside $src_tree, we should add "." before
        # --json
        #
        gcovr_args=()
        gcovr_args+=(-j0)
        gcovr_args+=(--exclude-lines-by-pattern '.*ASSERT_ALWAYS')
        gcovr_args+=(--merge-mode-functions=separate)
        gcovr_args+=(-r "$source_tree")
        gcovr_args+=(--json "$outfile")
        gcovr_args+=(--gcov-ignore-parse-errors=suspicious_hits.warn_once_per_file)
        gcovr_args+=(--exclude-throw-branches)
        (set -x ; cd "$build_tree" ; time gcovr "${gcovr_args[@]}")
    fi

    if true ; then
        purge_unused_coverage_files

        # lcov's --no-external is hopeless as long as
        # https://github.com/linux-test-project/lcov/commit/d73281a15 is
        # not checked in.
        outfile_pre0="$source_tree/$1-$2-pre0.info"
        outfile_pre1="$source_tree/$1-$2-pre1.info"
        outfile="$source_tree/$1-$2.info"

        lcov_capture_args=(
            -q
            -c
            --directory "$build_tree"
            --exclude "$build_tree"
            --exclude "/usr/*"
            --ignore-errors unused,inconsistent
            --external
        )
        if [ "$2" = "base" ] ; then
            lcov_capture_args+=(-i)
        fi
        lcov_capture_args+=(-o "$outfile_pre0")
        lcov "${lcov_capture_args[@]}"

        # do removals before path simplification, because the matches are
        # shell globs but not anchored...
        
        lcov_removal_args=(
            --exclude "$PWD/utils/embedded/*"
            --exclude "$PWD/linalg/bwc/flint-fft/*"
            --exclude "$PWD/gf2x/*"
            --ignore-errors unused,inconsistent
        )
        lcov -a "$outfile_pre0" "${lcov_removal_args[@]}" -o "$outfile_pre1"

        # now simplify the paths
        lcov -a "$outfile_pre1" --substitute "s#^$PWD/##" -o "$outfile"

        rm -f "$outfile_pre0"
        rm -f "$outfile_pre1"
    fi
}

step_coverage_more_artifacts() {
    return
    # prefix="$1"
    # if [ "$build_tree" != generated ] ; then
    #     fatal_error "This part of the script assumes that build_tree=generated"
    # fi

    # # because of /bin/sh, we can't do arrays.
    # find "$build_tree" -name '*.[ch]' -o -name '*.[ch]pp' | xargs -x tar czf ${prefix}-generated-sources.tar.gz
}


coverage_expunge_paths="utils/embedded:gf2x:generated:linalg/bwc/flint-fft:linalg/bwc/mpfq"

dispatch_valgrind_files() {
    cd $vdir
    mkdir ok ok-signal nok system
    find . -type f -a -name 'pid-*' | xargs egrep -l "Command: (/usr/bin|/bin|python|perl|env|[^ ]*\.sh)" | xargs -r mv --target-directory system
    # the rm -rf step could be considered an option
    rm -rf system
    # SEGV is something to worry about, but there are cases where we
    # terminate with SIGTERM / SIGINT / SIGHUP and this is just normal
    # business (e.g., cado-nfs-client.py can do that). It is possible
    # that vlagrind report leaks in such cases, but we're not super
    # interested in them
    # SIGABRT is also what we get when an expect-fail test aborts on an
    # exception. Likewise, there is little to worry about _in the
    # valgrind setting_ about aborts in general. (If a SIGABRT error
    # happens for a reason that is not an expect-fail, then the other
    # tests should catch it!)
    ls | grep pid | xargs -r egrep -l 'ERROR SUMMARY: 0' | xargs -r mv -t ok
    ls | grep pid | xargs -r egrep -l 'Process terminating.*signal.*SIG(TERM|INT|HUP|ABRT)' pid-* | xargs -r mv -t ok-signal
    ls | grep pid | xargs -r egrep -l 'ERROR SUMMARY: [^0]' pid-* | xargs -r mv -t nok
}

postprocess_valgrind() {
    (dispatch_valgrind_files)

    set +e
    nok_files=($(find "$vdir/nok" -type f))
    ok_files=($(find "$vdir/ok" -type f))
    ok_signal_files=($(find "$vdir/ok-signal" -type f))
    set -e

    if [ ${#nok_files[@]} -gt 0 ] ; then
        red_message "Found valgrind errors (${#nok_files[@]} different executions)" >&2
    fi

    for f in "${nok_files[@]}" ; do
        cmd=$(perl -ne 'm{Command: \S*/([^/\s]+)} && print "$1\n";' "$f")
        nerr=$(perl -ne 'm{ERROR SUMMARY: (\d+) errors from (\d+) contexts} && print "$1 from $2\n";' "$f")
        enter_section collapsed errors "Errors in $cmd ($nerr)"
        cat "$f"
        leave_section
    done
    tar czf $vdir.tar.gz $vdir/
    rm -rf $vdir
    if [ $rc != 0 ] ; then
        red_message "exit code was $rc" >&2
        exit $rc
    fi
    if [ ${#nok_files[@]} -gt 0 ] ; then
        fatal_error "Found valgrind errors (${#nok_files[@]} different executions)" "See archive of log files in $vdir.tar.gz" 
    else
        green_message "valgrind passed successfully (${#ok_files[@]} different executions)"
        if [ "${#ok_signal_files[@]}" ] ; then
            yellow_message "NOTE: valgrind reported (possibly spurious) errors on ${#ok_signal_files[@]} executions that were terminated by SIG{TERM,INT,HUP}"
        fi
        green_message "See archive of log files in $vdir.tar.gz"
    fi
}

step_check() {
    # --no-compress-output is perhaps better for test uploading, as ctest
    # likes to store as zlib but headerless, which is a bit of a pain
    #
    # -V is to get the output of tests. We want it, since anyway for
    # practical purposes our ctest filter does the required filtering.

    ctest_args=(-V -T Test --no-compress-output --test-output-size-passed 4096 --test-output-size-failed 262144)

    if [ "$specific_checks" = "bwc.sagemath" ] ; then
        ctest_args+=(-R with_sagemath)
        # it's only for our sage-in-docker script, but we really want
        # this in order to avoid long pulls from runners.
        # Note that we'll pull anyway if the image is not there.
        export DOCKER_SAGEMATH_NO_PULL=1
    elif [ "$specific_checks" = "including_mpi" ] ; then
        # nothing to do
        :
    elif [ "$specific_checks" = "only_mpi" ] ; then
        ctest_args+=(-R mpi)
    elif [ "$specific_checks" = "mysql" ] ; then
        ctest_args+=(-R mysql)
    fi

    ctest_args+=(-E ^builddep)

    if [[ $CI_JOB_NAME =~ ([[:digit:]]+)/([[:digit:]]+) ]] ; then
        stride_args+=(-I ${BASH_REMATCH[1]},,${BASH_REMATCH[2]})
        ctest_args+=("${stride_args[@]}")
        enter_section collapsed all_tests "List of tests to run"
        (cd "$build_tree" ; ctest -N "${ctest_args[@]}")
        leave_section
    fi

    if [ "$using_cmake_directly" ] ; then
        set -o pipefail
        (cd "$build_tree" ; "${test_precommand[@]}" ctest -j$NCPUS "${ctest_args[@]}" ) | "$source_tree"/scripts/filter-ctest.pl
    else
        "${test_precommand[@]}" "${MAKE}" check ARGS="-j$NCPUS ${ctest_args[*]}"
    fi
    rc=$?
    export rc

    if [ "$valgrind" ] ; then
        (set +x ; postprocess_valgrind)
    else
        return $rc
    fi
}

step_doc() { : ; }
