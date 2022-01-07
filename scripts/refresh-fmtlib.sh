#!/bin/bash
set -ex
GH_CATALOG=https://api.github.com/repos/fmtlib/fmt/releases/latest
PY_PARSER='import sys,json; print(json.load(sys.stdin)["tarball_url"]);'
URL=$(curl -sL "$GH_CATALOG" | python3 -c "$PY_PARSER")
tmp=$(mktemp -d /tmp/XXXXXXX)
trap "rm -rf $tmp" EXIT
curl -sL "$URL" | (cd $tmp ; tar xzf -)
rm -rf utils/embedded/fmt/*
cp $tmp/*/*.rst $tmp/*/include/fmt/* $tmp/*/src/* utils/embedded/fmt
git add utils/embedded/fmt/*
