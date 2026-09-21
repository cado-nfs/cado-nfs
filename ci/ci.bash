
tweak_tree_before_configure() { : ; }

user_variables() {
    if [ "$ninja_build" ] ; then
        export CMAKE_GENERATOR=Ninja
    fi
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
    if is_ubuntu ; then
        . /etc/lsb-release
        # force newer compiler on ubuntu 24
        # case "$DISTRIB_RELEASE" in
        #     24*)
        #         export CC=gcc-14
        #         export CXX=g++-14
        # esac
    fi
    if [ "$using_cmake_directly" ] ; then
        (cd "$build_tree" ; if [ "$ninja_build" ] ; then unset MAKE; fi ; cmake "$source_tree" $pass_flags_to_cmake)
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
        (cd "$build_tree" ; B=${MAKE}; if [ "$ninja_build" ] ; then B=ninja; fi; $B -j$NCPUS $target)
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
        (cd "$build_tree" ; B=${MAKE}; if [ "$ninja_build" ] ; then B=ninja; fi; $B -j$NCPUS $target)
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

# Third-party or generated code that has no business showing up in our
# coverage reports. This list lives here and nowhere else: it drives both
# the pruning of counter files in the test jobs and the exclusions that
# the report job passes to grcov.
coverage_excluded_paths="utils/embedded gf2x linalg/bwc/flint-fft linalg/bwc/mpfq"

# Where step_coverage_report writes. The defaults are what .gitlab-ci.yml
# collects; scripts/cov.sh points them inside the build tree instead, so
# that a local run does not litter the checkout.
coverage_report_dir="coverage"
coverage_report_xml="coverage.xml"
coverage_report_workdir="cov-work"

purge_unused_coverage_files() {
    # We don't want to bother with traces of config checks, nor with code
    # that we exclude from the reports anyway -- dropping the counter
    # files right away is cheaper than filtering them out afterwards.
    (
        find "$build_tree" -name '*conftest*gcno' -o -name 'CMake*.gcno' -o -name '?-CMake*.gcno'
        find "$build_tree" -name '*conftest*gcda' -o -name 'CMake*.gcda' -o -name '?-CMake*.gcda'
        for d in $coverage_excluded_paths ; do
            find "$build_tree/$d" -name \*.gcda -o -name \*.gcno 2>/dev/null || :
        done
    ) | sort -u | xargs -r rm -f
}

# tar(1) on alpine is busybox tar, which understands -T but none of the
# fancier GNU spellings. Going through a file list works with both.
tar_files_under() {
    # $1 = directory to look under ; $2 = find pattern ; $3 = output tarball
    local list="$3.filelist"
    (cd "$1" && find . -name "$2" > "$list" && tar czf "$3" -T "$list")
    rm -f "$list"
    ls -l "$3"
}

step_coverage() {
    # $1 is a prefix for the artifact file names.
    #
    # All we do here is save the raw gcov counters. Merging them and
    # rendering a report is the business of ci/ci/09-merge-coverage.sh,
    # which does it once for the whole pipeline rather than once per job.
    # This takes about a second, where rendering here used to take
    # several minutes in every single job.
    #
    # ci/ci/005-build-environment.sh sets build_tree to "generated" for
    # coverage jobs, so both the .gcno and the .gcda live under
    # $PWD == $source_tree.

    purge_unused_coverage_files

    tar_files_under "$build_tree" '*.gcda' "$source_tree/$1-gcda.tar.gz"

    # The .gcno files describe the shape of the code. They are identical
    # in every job of the pipeline, and two orders of magnitude bigger
    # than the counters, so exactly one job ships them. Which one is
    # decided in ci/ci.sh.
    if [ "$coverage_ship_gcno" ] ; then
        tar_files_under "$build_tree" '*.gcno' \
            "$source_tree/coverage-$COMMIT_SHORT_SHA-gcno.tar.gz"
    fi
}

step_coverage_report() {
    # Called by ci/ci/09-merge-coverage.sh, once the counters of every
    # coverage job have been merged. $coverage_tree is a directory that
    # holds the merged .gcda next to the .gcno they belong to; it is not
    # a build tree, and it holds no object file.
    #
    # We produce $coverage_report_dir (browsable html, plus the lcov
    # tracefile that it was built from) and $coverage_report_xml
    # (cobertura, which is what gitlab reads).

    local work="$coverage_report_workdir/report"

    # Check this here rather than let the first missing program blow up
    # half way through, inside a collapsed section, leaving the job to
    # be diagnosed from "coverage.xml is missing from the artifacts".
    local missing= prog
    for prog in grcov genhtml perl ; do
        type -p "$prog" > /dev/null 2>&1 || missing="$missing $prog"
    done
    if [ "$missing" ] ; then
        fatal_error "Cannot build the coverage report, missing:$missing"
    fi

    # grcov insists on the output directory existing beforehand when it
    # is asked for more than one output type.
    rm -rf "$work" ; mkdir -p "$work"

    # grcov reports paths relative to -s, so the build tree has to be
    # named that way too. In the ci jobs it already is ("generated"); in
    # a local run it is an absolute path under the checkout.
    local bt_rel="$build_tree"
    case "$bt_rel" in "$source_tree"/*) bt_rel="${bt_rel#"$source_tree"/}";; esac

    grcov_args=(
        "$coverage_tree"
        -s "$source_tree"
        -t lcov,cobertura
        -o "$work"
        --branch
        --ignore-not-existing
        # Lines that only report a fatal error and leave. They are a
        # fifth of everything the report shows as uncovered, and no test
        # is ever going to reach them.
        --excl-line '(ASSERT_ALWAYS|\bthrow\b|\.fail\(|exit\s*\(\s*EXIT_FAILURE|abort\s*\(\s*\))'
        --ignore '/usr/*'
        --ignore "$bt_rel/*"
    )
    for d in $coverage_excluded_paths ; do
        # both spellings: the code as it sits in the checkout, and the
        # copies that the build system generates under the build tree.
        grcov_args+=(--ignore "$d/*")
        grcov_args+=(--ignore "*/$d/*")
    done

    enter_section collapsed grcov "Reading the merged coverage data"
    local gcov_log="$coverage_report_workdir/grcov.log"
    local rc=0
    (set -x ; time grcov "${grcov_args[@]}") > "$gcov_log" 2>&1 || rc=$?
    # A stamp mismatch means the .gcno and the .gcda were produced by
    # different builds. gcov does not fail on it -- it reports 0% and
    # carries on -- and grcov dumps the whole of gcov's output for every
    # file, which is megabytes. Catch it here, where we can say why.
    if grep -q 'stamp mismatch' "$gcov_log" ; then
        grep -m 2 -B 1 'stamp mismatch' "$gcov_log" >&2
        fatal_error "The .gcda and the .gcno come from different builds."   \
            "Every coverage job must compile with the same -frandom-seed;"  \
            "see ci/ci/005-build-environment.sh. Without it the report"     \
            "would silently come out at 0%."
    fi
    tail -c 200000 "$gcov_log"
    [ "$rc" = 0 ] || fatal_error "grcov failed (exit $rc)"
    leave_section

    for f in "$work/lcov" "$work/cobertura.xml" ; do
        [ -s "$f" ] || fatal_error "grcov did not produce $f"
    done

    mkdir -p "$coverage_report_dir"
    cp "$work/lcov" "$coverage_report_dir/coverage.info"

    # gitlab reads coverage.xml to annotate merge request diffs, and
    # silently ignores the file when it is above 10MiB. What it renders
    # is whether a line is covered; the <methods> and <conditions>
    # blocks, which it does not look at, are most of the weight. On one
    # pipeline: 49MB raw, 11MB without <methods>, 5MB without either.
    perl -0777 -pi                                      \
        -e 's{<methods>.*?</methods>}{<methods/>}gs;'    \
        -e 's{<conditions>.*?</conditions>}{}gs;'        \
        "$work/cobertura.xml"
    cp "$work/cobertura.xml" "$coverage_report_xml"
    ls -l "$coverage_report_xml"
    sz=$(wc -c < "$coverage_report_xml")
    if [ "$sz" -gt $((10 * 1024 * 1024)) ] ; then
        yellow_message "$coverage_report_xml is past gitlab's 10MiB limit ($sz bytes)," \
            "so the merge request diff annotations will be missing"
    fi

    # genhtml gives us something browsable. grcov can write html on its
    # own, but its top-level index lists only the files that sit at the
    # root of the source tree, and offers no way down into the
    # subdirectories.
    # genhtml puts this where it would otherwise print a date. Link to
    # the commit when we know where the project lives, and settle for the
    # bare sha when we do not -- a relative href would resolve against the
    # pages host and dangle.
    local commit_html="$COMMIT_SHORT_SHA"
    if [ "$PROJECT_URL" ] ; then
        commit_html="<a href=\"$PROJECT_URL/-/commit/$COMMIT_SHORT_SHA\">$COMMIT_SHORT_SHA</a>"
    fi
    enter_section collapsed genhtml "Rendering the coverage report"
    # About --prefix: the tracefile has paths relative to the top of the
    # source tree, and genhtml resolves those against $PWD. We strip
    # everything *above* the source tree but not the source tree itself,
    # because genhtml 2.5 cannot place a file that ends up at the very
    # root of the prefix -- with --prefix "$source_tree" it dies with
    # "file error for portability.h" and takes its parallel workers down
    # with it. So the report has one extra level, named after the
    # checkout directory, and that is a fair price.
    genhtml -q --parallel "$NCPUS"                                          \
        --prefix "$(dirname "$source_tree")"                                \
        --hierarchical                                                      \
        --ignore-errors inconsistent,category,unmapped,source,corrupt,count \
        --current-date "$commit_html"                                       \
        -o "$coverage_report_dir" "$coverage_report_dir/coverage.info"
    leave_section

    # Because we cannot strip the last component of the prefix (see the
    # comment above), the page genhtml calls "top level" holds a single
    # row, and it is *not* the page you want to land on: the sort-by-line
    # and sort-by-function links only exist from the level below. Point
    # the reader there. README.md links to $coverage_report_dir/index.html
    # and should keep working without knowing the checkout's name.
    local subs=() x
    for x in "$coverage_report_dir"/*/ ; do
        [ -d "$x" ] || continue
        subs+=("$(basename "$x")")
    done
    if [ "${#subs[@]}" = 1 ] && [ -f "$coverage_report_dir/${subs[0]}/index.html" ] ; then
        cat > "$coverage_report_dir/index.html" <<EOF
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<meta http-equiv="refresh" content="0; url=${subs[0]}/index.html">
<title>cado-nfs coverage report</title>
</head>
<body><a href="${subs[0]}/index.html">cado-nfs coverage report</a></body>
</html>
EOF
    fi

    # gitlab scrapes the job log for the number that goes on the coverage
    # badge (see the "coverage:" key in .gitlab-ci.yml). Print a line of
    # our own for it to latch onto, rather than depend on the exact
    # wording of whichever tool happens to print a summary. Plain echo,
    # because major_message would wrap it in escape sequences.
    perl -0777 -ne '
        if (m{<coverage\s+([^>]*?)>}) {
            my %a = ($1 =~ m{([\w-]+)="([^"]*)"}g);
            printf "coverage: lines %.2f%% (%d of %d), branches %.2f%% (%d of %d)\n",
                100 * $a{"line-rate"}, $a{"lines-covered"}, $a{"lines-valid"},
                100 * $a{"branch-rate"}, $a{"branches-covered"}, $a{"branches-valid"};
        }' "$coverage_report_xml"
}


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
