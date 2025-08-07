#!/bin/bash

set -e
set -x

# Usage: (see also dev_docs/README.coverage.md)
#
# COV=1 make -j8 all all_test_dependencies && COV=1 scripts/cov.sh -R numbertheory
#
# (it's not absolutely necessary to compile everything)
#
# The COV=1 magic assumes that you have this code inside `local.sh`:
#
# if [ "$COV" ] ; then
#     build_tree="${build_tree}.cov"
#     DEBUG=1
#     CFLAGS="-O0 -g --coverage -fprofile-update=atomic"
#     CXXFLAGS="-O0 -g --coverage -fprofile-update=atomic"
#     LDFLAGS="--coverage"
# fi
 
export COV=1
eval $(make show)

# This must be an absolute path
C=/tmp/cov

commit=$(./scripts/version.sh)
commit_ref="https://gitlab.inria.fr/cado-nfs/cado-nfs/-/commit/$commit"

export src_tree=$PWD

run() {
    OMP_DYNAMIC=true make check ARGS="-E builddep $*"
    # COV=1 ./cado-nfs.py 90377629292003121684002147101760858109247336549001090677693
}

case "$1" in
    --lcov) lcov=1; shift;;
    --gcovr) gcovr=1; shift;;
    *) lcov=1;;
esac

if [ "$gcovr" ] ; then
    software_dependencies() {
        if ! [ -e "cov/bin/activate" ] ; then
            python3 -m venv cov
        fi
        . cov/bin/activate
        if ! [ -x "cov/bin/gcovr" ] ; then
            pip3 install gcovr
        fi
    }

    before_run() {
        # remove coverage results of the previous run
        (cd $build_tree ; find . -name '*.gcda' | xargs -r rm -f)
        (cd $build_tree ; time gcovr --merge-mode-functions=separate -r $src_tree --json ${C}-base.json)
    }

    after_run() {
        (cd $build_tree ; time gcovr --merge-mode-functions=separate -r $src_tree --json ${C}-app.json)
    }

    final_processing() {
        rm -rf "$C" || :
        mkdir $C
        gcovr --merge-mode-functions=separate -a ${C}-base.json -a ${C}-app.json --html-title "Coverage for commit <a href=\"$commit_ref\">$commit</a>" --html-nested $C/coverage.html --print-summary
    }

    open_url() { gio open $C/coverage.html ; }

else
    software_dependencies() { : ; }
    capture() {
        cap="$1"
        shift
        find $build_tree -name '*conftest*' | xargs -r rm -f
        lcov -q -c      \
            --directory $build_tree     \
            --exclude $build_tree       \
            --exclude /usr/     \
            --ignore-errors unused      \
            --external  \
            -o $cap "$@"
    }
    postprocess_capture() {
        cap0="$1"
        cap1="$2"
        cap="$3"
        lcov -a $cap0   \
            --exclude "$src_tree/utils/embedded/*" \
            --exclude "$src_tree/linalg/bwc/flint-fft/*" \
            --exclude "$src_tree/gf2x/*" \
            --exclude "$(realpath $build_tree)"    \
            --ignore-errors unused \
            -o $cap1
        lcov -a $cap1 --substitute "s#^$src_tree/##" -o $cap
    }
    before_run() {
        # remove coverage results of the previous run
        (cd $build_tree ; find . -name '*.gcda' | xargs -r rm -f)
        capture ${C}-pre0.info -i
        postprocess_capture ${C}-pre0.info ${C}-pre1.info ${C}-base.info
    }
    after_run() {
        capture ${C}-post0.info
        postprocess_capture ${C}-post0.info ${C}-post1.info ${C}-app.info
    }
    final_processing() {
        rm -rf "$C" || :
        mkdir $C
        genhtml -q -o $C ${C}-base.info ${C}-app.info
    }
    open_url() { gio open $C/index.html ; }
fi

software_dependencies
before_run
run "$@"
after_run
final_processing
open_url
