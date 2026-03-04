#!/bin/sh

build_tree="$CADO_NFS_BINARY_DIR"
source_tree="$CADO_NFS_BINARY_DIR"

while [ $# -gt 0 ] ; do
    if [ "$1" = "-b" ] ; then
        shift
        build_tree="$1"
        shift
    elif [ "$1" = "-s" ] ; then
        shift
        source_tree="$1"
        shift
    else
        echo "bad arg: $1" >&2
        exit 1
    fi
done


NCPUS=$("`dirname $0`"/ncpus.sh)

# We expect to be called with provide-wdir.sh, so that we don't need to
# create a directory by ourselves.
: ${wdir=${TMPDIR-/tmp}}

t="${wdir}"

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
computation = DLP
N = 8005493
ell = 64087926178543
gfpext = 3

slaves.nrclients = $(((1+NCPUS)/2))
tasks.threads = 2
tasks.linalg.bwc.threads = $NCPUS
tasks.execpath = $build_tree
slaves.scriptpath = $source_tree
tasks.workdir = $WDIR
slaves.basepath= $WDIR/client
slaves.hostnames = localhost

tasks.polyselect.import = $POLYFILE

# for the record, the computations involve 1 unit on side 1, 1 SM on side 0

tasks.I = 11
tasks.polyselect.degree = 4
tasks.polyselect.admax = 0
tasks.polyselect.incr = 60
tasks.polyselect.adrange = 500
tasks.polyselect.P = 420
tasks.polyselect.nq = 1000

lim0 = 20000
lim1 = 20000
qmin = 20000
lpb0 = 16
lpb1 = 16
tasks.sieve.mfb0 = 32
tasks.sieve.mfb1 = 32
tasks.sieve.qrange = 1000
tasks.sieve.rels_wanted = 20000

tasks.linalg.allow_zero_on_rhs = 1
tasks.reconstructlog.partial = true
checkdlp = false
EOF

${build_tree}/cado-nfs.py $PARAMFILE
