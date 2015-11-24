#!/bin/sh

CADO_NFS_SOURCE_DIR=$1
CADO_NFS_BINARY_DIR=$2

NCPUS=$("`dirname $0`"/ncpus.sh)

t=`mktemp -d /tmp/cado-check.XXXXXXX`

POLYFILE=$t/p3dd7-f4g3-GJL-1.poly
PARAMFILE=$t/p3dd7-f4g3-GJL-1.params
WDIR=$t/p3dd7

mkdir $WDIR

cat > $POLYFILE <<EOF
n: 8005493
# ell: 64087926178543
skew: 1.000
poly0: 7066873,2189883,5866272,1
poly1: 1,3,4,-2,1
# f := x^4 - 2*x^3 + 4*x^2 + 3*x + 1;
# sgn(f) = (0, 2) => rk = 1
# g := x^3 + 5866272*x^2 + 2189883*x + 7066873
# sgn(g) = (1, 1) => rk = 1
EOF


cat > $PARAMFILE <<EOF
name = p3dd7-f4g3-GJL-1
dlp = true
N = 8005493
ell = 64087926178543

slaves.nrclients = $(((1+NCPUS)/2))
tasks.threads = 2
tasks.linalg.bwc.threads = $NCPUS
tasks.execpath = $CADO_NFS_BINARY_DIR
slaves.scriptpath = $CADO_NFS_SOURCE_DIR
tasks.workdir = $WDIR
slaves.basepath= $WDIR/client
slaves.hostnames = localhost

tasks.polyselect.import = $POLYFILE

# for the record, the computations involve 1 unit on side 1, 1 SM on side 0
tasks.lcideals = true

tasks.I = 11
tasks.polyselect.degree = 4
tasks.polyselect.admax = 0
tasks.polyselect.incr = 60
tasks.polyselect.adrange = 500
tasks.polyselect.P = 420
tasks.polyselect.nq = 1000

rlim = 20000
alim = 20000
lpbr = 15
lpba = 15
tasks.sieve.mfbr = 32
tasks.sieve.mfba = 32
tasks.sieve.qrange = 1000
tasks.sieve.rels_wanted = 10000

tasks.linalg.allow_zero_on_rhs = 1
tasks.reconstructlog.partial = true
checkdlp = false
EOF

cleanup() {
    if ! [ "$CADO_DEBUG" ] ; then
        rm -rf $t
    else
        echo "(debug mode, temporary files are kept in $t)"
    fi
}

${CADO_NFS_SOURCE_DIR}/cado-nfs.py $PARAMFILE && cleanup
