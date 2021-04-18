#!/bin/bash

set -x

# As such, this test is unfortunately a bit long.

dash_args=(
    -ell 6039154416103261588907535248367667
    -dlp
)

params_args=(
    slaves.nrclients=1
    slaves.hostnames=localhost
    tasks.sieve.las.threads=auto
    tasks.polyselect.import="$CADO_NFS_SOURCE_DIR/tests/estimate_matsize/p80.poly"
    tasks.workdir=${wdir?missing}
)

paramfile_args=(
    "$CADO_NFS_SOURCE_DIR/tests/estimate_matsize/p80_compsq.params"
)

freeform_args=(
    33214702831050289432509431979913881660169847863978111477101554095620630619037
)

args=("${dash_args[@]}" "${params_args[@]}" "${paramfile_args[@]}" "${freeform_args[@]}")

"$CADO_NFS_SOURCE_DIR/cado-nfs.py" "${args[@]}"

mkdir "$wdir/p80.split"
find $wdir/p80.upload -name '*.gz' | xargs "$CADO_NFS_SOURCE_DIR/scripts/estimate_matsize/build_relation_cache.pl"  -o "$wdir/p80.split" -s 1000,10000

export CADO_BUILD="$PROJECT_BINARY_DIR"

# can also go in the .params file if we want.
export relation_cache="$wdir/p80.split"

"$CADO_NFS_SOURCE_DIR/scripts/estimate_matsize/estimate_matsize.sh" -params "$CADO_NFS_SOURCE_DIR/tests/estimate_matsize/p80_compsq.ems.params" "$CADO_NFS_SOURCE_DIR/tests/estimate_matsize/p80.poly"
