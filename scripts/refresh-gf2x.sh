#!/usr/bin/env bash

set -e
set -x

: ${gf2x_url:=https://gitlab.inria.fr/gf2x/gf2x}
: ${gf2x_rev:=master}
# : ${gf2x_rev:=gf2x-1.2.x}
rm -rf gf2x

checkout_that() {
    url=$1
    rev=$2
    path=$3
    mkdir $path
    OPWD="$PWD"
    cd "$path"
    git init
    git remote add origin ${url}
    git fetch origin ${rev}
    git reset --hard FETCH_HEAD
    cd "$OPWD"
    rm -rf "$path/.git"
    rm -rf "$path/apps"
    sed -e s/apps// -i $path/Makefile.am
    sed -e s/apps.Makefile// -i $path/configure.ac
}

checkout_that ${gf2x_url} ${gf2x_rev} gf2x

cp -f gf2x/toom-gpl-placeholder.c gf2x/toom-gpl.c

sed -e "/^AM_MAINTAINER_MODE/ s/enable/disable/" -i gf2x/configure.ac

(cd gf2x/ ; autoreconf -i)
(cd gf2x/ ; xargs -r rm -f < no-distribute.txt)
# find gf2x/ -type f | xargs -r git add
rm -rf gf2x/jenkins
rm -f gf2x/.mailmap || :
git add gf2x/
