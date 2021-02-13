#!/usr/bin/env bash

tmp=$(mktemp -d /tmp/XXXXXXXXXX)
# trap "rm -rf $tmp" EXIT
curl -sL -o $tmp/archive.zip "https://gitlab.inria.fr/cado-nfs/cado-nfs/-/jobs/artifacts/master/download?job=merge%20coverage%20tests"
cd $tmp
unzip archive.zip
xdg-open coverage/index.html
