#!/bin/bash
set -ex

if [ $# -gt 1 ] && [ "$1" = "-g" ] ; then
    git="$2"
    URL="https://github.com/fmtlib/fmt/archive/$git.tar.gz"
elif [ $# -gt 1 ] && [ "$1" = "-h" ] ; then
    git="refs/heads/$2"
    URL="https://github.com/fmtlib/fmt/archive/$git.tar.gz"
else
    release=latest
    if [ $# -gt 1 ] && [ "$1" = "-r" ] ; then
        release="tags/$2"
    fi
    GH_CATALOG="https://api.github.com/repos/fmtlib/fmt/releases/$release"
    PY_PARSER='import sys,json; print(json.load(sys.stdin)["tarball_url"]);'
    URL=$(curl -sL "$GH_CATALOG" | python3 -c "$PY_PARSER")
fi

tmp=$(mktemp -d /tmp/XXXXXXX)
trap "rm -rf $tmp" EXIT

curl -sL "$URL" | (cd $tmp ; tar xzf -)
rm -rf utils/embedded/fmt/*
# find $tmp/*/*.rst
cp $tmp/*/include/fmt/* $tmp/*/src/* utils/embedded/fmt
cp $tmp/*/ChangeLog.md $tmp/*/README.md $tmp/*/LICENSE utils/embedded/fmt
git add utils/embedded/fmt/*
