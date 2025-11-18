#!/bin/bash

set -e
set -x

# Usage: (see also dev_docs/README.coverage.md)
#
# COV=1 make -j8 all all_test_dependencies && COV=1 scripts/cov.sh -R numbertheory
#
# (it's not absolutely necessary to compile everything)
#
# The COV=1 CLANG=1 magic assumes that you have this code inside `local.sh`:
#
# if [ "$COV" ] ; then
#     build_tree="${build_tree}.cov"
#     DEBUG=1
#     if [ "$CLANG" ] ; then
#         : ${CFLAGS="-O0 -g -fprofile-instr-generate -fcoverage-mapping"}
#         : ${CXXFLAGS="-O0 -g -fprofile-instr-generate -fcoverage-mapping"}
#         : ${LDFLAGS="-O0 -g -fprofile-instr-generate -fcoverage-mapping"}
#     else
#         : ${CFLAGS="-O0 -g --coverage -fprofile-update=atomic"}
#         : ${CXXFLAGS="-O0 -g --coverage -fprofile-update=atomic"}
#         : ${LDFLAGS="--coverage"}
#     fi
# fi
# 
# Note that currently --llvm does not work as well as I thought. It's
# fast, but it doesn't display an aggregated lcov report

case "$1" in
    --llvm) llvm=1; lcov=1; export CLANG=1; shift;;
    --lcov) lcov=1; shift;;
    --gcovr) gcovr=1; shift;;
    *) lcov=1;;
esac

if [ "$CLANG" ] && ! [ "$llvm" ] ; then
    echo "\$CLANG requires --llvm" >&2
    exit 1
fi

# the llvm package must be installed
export COV=1
eval $(make show)

# This must be an absolute path
: ${COV_BASE=/tmp/cov}

TIMESTAMP=$(date +%Y%m%d%H%M%S%N)

DATA=$COV_BASE/$TIMESTAMP.data
mkdir -p "$DATA"

RENDER=$COV_BASE/$TIMESTAMP.render

if [ "$llvm" ] ; then
    export LLVM_PROFILE_FILE="$DATA/%32m.profraw"
fi

commit=$(./scripts/version.sh)
commit_ref="https://gitlab.inria.fr/cado-nfs/cado-nfs/-/commit/$commit"

export src_tree=$PWD

run() {
    OMP_DYNAMIC=true make check ARGS="-E builddep $*"
    # COV=1 ./cado-nfs.py 90377629292003121684002147101760858109247336549001090677693
}

list_cado_binaries() {
    find "$build_tree"  \
        \( -name CMakeFiles -a -prune \) -o \( \
            -type f     \
            -a -perm -111       \
            -a \! -name '*.so*' \
            -a \! -name '*.status'      \
            -a \! -name 'libtool'      \
            -a \! -name 'gen_bb_mul_code'      \
            -a \! -name '*.p[yl]'       \
            -a -print   \
        \)
}

if [ "$gcovr" ] ; then
    gcovr_try_parallel() {
        gcovr -j$(nproc) "$@" || gcovr -j0 "$@" || gcovr "$@"
    }

    gcovr_args_common=(--exclude-throw-branches --merge-mode-functions=separate --exclude-lines-by-pattern 'ASSERT_ALWAYS')

    software_dependencies() {
        if ! [ -e "cov/bin/activate" ] ; then
            python3 -m venv --system-site-packages cov
        fi
        . cov/bin/activate
        if ! [ -x "cov/bin/gcovr" ] ; then
            pip3 install gcovr
        fi
    }

    before_run() {
        # remove coverage results of the previous run
        (cd $build_tree ; find . -name '*.gcda' | xargs -r rm -f)
        (cd $build_tree ; time gcovr_try_parallel "${gcovr_args_common[@]}" -r $src_tree --json $DATA/base.json)
    }

    after_run() {
        (cd $build_tree ; time gcovr_try_parallel "${gcovr_args_common[@]}" -r $src_tree --json $DATA/app.json)
    }

    final_processing() {
        rm -rf "$RENDER" || :
        mkdir $RENDER
        gcovr_try_parallel "${gcovr_args_common[@]}" -a $DATA/base.json -a $DATA/app.json --html-title "Coverage for commit <a href=\"$commit_ref\">$commit</a>" --html-nested $RENDER/coverage.html --print-summary
    }

    open_url() { gio open $RENDER/coverage.html ; }
elif [ "$lcov" ] && ! [ "$llvm" ] ; then
    (grep '^[^#]' /etc/lcovrc ; echo "check_data_consistency = 0") > $DATA/lcovrc
    software_dependencies() { : ; }
    capture() {
        mode="$1"
        shift
        cap="$1"
        shift
        find $build_tree -name '*conftest*' | xargs -r rm -f
        rm -f $cap
        lcov -q $mode      \
            --directory $build_tree     \
            --exclude $build_tree       \
            --exclude /usr/     \
            --ignore-errors unused,mismatch,inconsistent      \
            --external  \
            -o $cap "$@"
    }
    postprocess_capture() {
        cap0="$1"
        cap1="$2"
        cap="$3"
        rm -f $cap1
        lcov -a $cap0   \
            --exclude "$src_tree/utils/embedded/*" \
            --exclude "$src_tree/linalg/bwc/flint-fft/*" \
            --exclude "$src_tree/gf2x/*" \
            --exclude "$(realpath $build_tree)"    \
            --ignore-errors unused,mismatch,inconsistent \
            -o $cap1
        rm -f $cap
        lcov -a $cap1 --substitute "s#^$src_tree/##" -o $cap
    }
    before_run() {
        # remove coverage results of the previous run
        (cd $build_tree ; find . -name '*.gcda' | xargs -r rm -f)
        capture -z $DATA/pre0.info -i
        capture -c $DATA/pre0.info -i
        postprocess_capture $DATA/pre0.info $DATA/pre1.info $DATA/base.info
    }
    after_run() {
        capture -c $DATA/post0.info
        postprocess_capture $DATA/post0.info $DATA/post1.info $DATA/app.info
    }
    final_processing() {
        rm -rf "$RENDER" || :
        mkdir $RENDER
        genhtml -q --config-file=$DATA/lcovrc -o $RENDER $DATA/base.info $DATA/app.info
    }
    open_url() { gio open $RENDER/index.html ; }
else
    (grep '^[^#]' /etc/lcovrc ; echo "check_data_consistency = 0") > $DATA/lcovrc
    software_dependencies() {
        type -p llvm-profdata
        type -p llvm-cov
    }
    before_run() { : ; }
    after_run() {
        llvm-profdata merge -sparse $DATA/*.profraw -o $DATA/profdata
        # all_binaries_with_coverage=()
        # for a in `list_cado_binaries` ; do
        #     if llvm-cov export --format=lcov --instr-profile $DATA/profdata $a > /dev/null ; then
        #         all_binaries_with_coverage+=("$a")
        #     fi
        # done
        all_binaries_with_coverage=(`list_cado_binaries`)
        # FIXME: does this even work?
        llvm-cov export --format=lcov --instr-profile $DATA/profdata "${all_binaries_with_coverage[@]}" > $DATA/app.info
    }
    final_processing() {
        rm -rf "$RENDER" || :
        mkdir "$RENDER" || :
        genhtml -q --config-file=$DATA/lcovrc -o $RENDER $DATA/app.info
    }
    open_url() { gio open $RENDER/index.html ; }
fi

software_dependencies
before_run
run "$@"
after_run
final_processing
open_url
