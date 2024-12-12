#!/usr/bin/env bash

eval $(make show)

changed_python_files() {
    if [ "$1" = "--all" ] ; then
        git ls-files '*.py'
    elif [ "$*" ] ; then
        echo "$@"
    else
        git diff --cached --name-only --diff-filter=ACM '*.py'
    fi
}

files=(`changed_python_files "$@"`)

if [ "${#files[@]}" = 0 ] ; then
    exit 0
fi

if ! [ -d "$build_tree/venv" ] ; then
    python3 -m venv "$build_tree/venv"
fi
if ! [ -x "$build_tree/venv/bin/flake8" ] ; then
    "$build_tree/venv/bin/pip" install flake8
fi

"$build_tree/venv/bin/flake8" "${files[@]}"

# mostly in good shape, but I'm modifying this in another branch too, so
# don't duplicate work here just yet!
# "$build_tree/venv/bin/flake8" --ignore=E741 tests/sagemath/cado_sage
