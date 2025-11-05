#!/usr/bin/env bash

if git diff --cached --name-only --diff-filter=ACM | grep -q 'ci/ci$' ; then
    w() { echo -e "### \e[01;31m$1\e[00m" ; }
    w "Your commit includes changes to ci/ci."
    w "If this is not intentional, please rework your commit before pushing."
    w "You might want to update your submodule checkout to its latest remote"
    w "status with the following command:"
    echo "###    git submodule update --remote"
    echo "### To make this a warning only: git config set policies.submodule-commits warn"
    exit 1
fi
