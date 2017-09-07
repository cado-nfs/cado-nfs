#!/usr/bin/env bash
#
# The purpose of this script is to help predict the size of the matrix
# that one would obtain with a given set of parameters.
# It takes as input a polynomial file, and sieving / filtering parameters
# given as env variables.
# The CADO_BUILD variable must point to the directory where cado-nfs was
# compiled.
#
# Typical usage will be something like:
#   wdir=/tmp/test I=12 lim0=... ./estimate_matsize.sh toto.poly
# The root files and the renumber file are cached, so if you change lim or
# lpb without changing wdir, you'll have to remove those cache-files
# first in wdir.
#
# If two-sided sieving is wanted, put comma-separated sides and qmax in
# the corresponding variables, for instance:
#   sqside=0,1
#   qmax=500000000,600000000

if [ $# != 1 ]; then
    echo "Usage: $0 <polyfile>"
    exit 1
fi

## default parameters: can be overriden using env variables
## these correspond more or less to a DLP-512.
: ${I=15}
: ${lim0=30000000}
: ${lim1=30000000}
: ${lpb0=27}
: ${lpb1=27}
: ${mfb0=54}
: ${mfb1=54}
: ${qmax=100000000}
: ${sqside=1}
: ${maxlevel=25}
: ${target_density=150}

## read poly file on command line
polyfile=$1
if [ ! -e $1 ]; then
    echo "Error: file $1 does not exist?"
    exit 1
fi

## if wdir is not set, build one in /tmp
if [ -z ${wdir+x} ]; then
    wdir=`mktemp -d ${TMPDIR-/tmp}/est_mat.XXXXXXX`;
else
    mkdir -p $wdir || (echo "mkdir -p $wdir failed"; false) || exit 1
fi
echo "Working directory is $wdir"

## if wdir does not contain a rootfile, build it
rootfile0="$wdir/roots0.gz"
rootfile1="$wdir/roots1.gz"
if [ ! -e $rootfile0 ]; then
    echo "Root file for side 0 not in wdir: building it..."
    $CADO_BUILD/sieve/makefb -poly $polyfile -lim $lim0 -side 0 -t 2 -out $rootfile0
fi
if [ ! -e $rootfile1 ]; then
    echo "Root file for side 1 not in wdir: building it..."
    $CADO_BUILD/sieve/makefb -poly $polyfile -lim $lim1 -side 1 -t 2 -out $rootfile1
fi

## if wdir does not contain a renumber table, build it
renumberfile=$wdir/renumber.gz
if [ ! -e $renumberfile ]; then
    echo "Renumber file not in wdir: building it..."
    $CADO_BUILD/sieve/freerel -poly $polyfile -renumber $renumberfile \
      -out /dev/null -pmax 1 -lpb0 $lpb0 -lpb1 $lpb1 -t 2
fi

## How many sides to sieve?
# convert sqside and qmax to arrays
IFS=',' read -ra array <<< "$sqside"
sqside=(${array[@]})
IFS=',' read -ra array <<< "$qmax"
qmax=(${array[@]})

nsides=${#sqside[@]}
qmin=()
echo "We sieve on $nsides sides."
if [ ${#sqside[@]} != ${#qmax[@]} ]; then
    echo "For multi-side sieving, length of sqside and qmax arrays must agree"
    exit 1;
fi
for i in `seq 0 $((nsides-1))`; do
    side=${sqside[$i]}
    qm=${qmax[$i]}
    lpb=lpb$side
    qqmax=`echo "2 ^ ${!lpb}" | bc`
    if [ $qqmax -le $qm ]; then
        echo "Error on side $side: qmax should be less then lpb"
        exit 1
    fi
    lim=lim$side
    qmin+=(${!lim})
    echo "  Side $side: qmin=${!lim} qmax=$qm"
done

## Sampling / faking on each side
NCHUNKS=5
NBSAMPLE=50
fakefiles=()
if [ $nsides == 1 ]; then
    dupqmax="0,0"
elif [ $nsides == 2 ]; then
    dupqmax="${qmax[0]},${qmax[1]}"
else
    echo "Can't deal with more than 2 sides yet."
    exit 1
fi
for i in `seq 0 $((nsides-1))`; do
    side=${sqside[$i]}
    echo "******************************"
    echo "Dealing with side $side..."
    qmax=${qmax[$i]}
    qmin=${qmin[$i]}
    qrange=$(((qmax-qmin)/NCHUNKS))

    ## sample with real sieving and build fake rels
    for i in `seq 1 $NCHUNKS`; do
        q0=$((qmin + (i-1)*qrange))
        q1=$((qmin + i*qrange))
        echo "Dealing with qrange=[$q0,$q1]"
        echo "  Sampling..."
        $CADO_BUILD/sieve/las -I $I -poly $polyfile -q0 $q0 -q1 $q1 \
          -lim0 $lim0 -lim1 $lim1 -lpb0 $lpb0 -lpb1 $lpb1 -sqside $side \
          -mfb0 $mfb0 -mfb1 $mfb1 \
          -fb0 $rootfile0 -fb1 $rootfile1 -random-sample $NBSAMPLE \
          -t 1 -dup -dup-qmax $dupqmax > $wdir/sample.side${side}.${q0}-${q1}
        echo "  Building fake relations..."
        $CADO_BUILD/sieve/fake_rels -poly $polyfile -lpb0 $lpb0 -lpb1 $lpb1 \
          -q0 $q0 -q1 $q1 -sqside $side \
          -sample $wdir/sample.side${side}.${q0}-${q1} \
          -renumber $renumberfile > $wdir/fakerels.side${side}.${q0}-${q1}
        fakefiles+=("$wdir/fakerels.side${side}.${q0}-${q1}")
    done
done
nrels=`cat ${fakefiles[@]} | grep -v "^#" | wc -l`
echo "We have $nrels fake relations"

## Filtering
# take a huge upper bound for the number of primes...
nprimes=`echo "2*2^$lpb0/l(2^$lpb0) + 2*2^$lpb1/l(2^$lpb1)" | bc -l | cut -d "." -f 1`

# purge
$CADO_BUILD/filter/purge -out $wdir/purged.gz -nrels $nrels -keep 3 \
    -col-min-index 10000 -col-max-index $nprimes -t 2 ${fakefiles[@]} \
    2>&1 | tee $wdir/purge.log

# Did we get a positive excess?
if (grep "number of rows < number of columns + keep" $wdir/purge.log > /dev/null); then
    echo "Negative excess: no way to build a matrix!"
    echo "You'll have to change the parameters"
    exit 1
fi

# merge
$CADO_BUILD/filter/merge-dl -mat $wdir/purged.gz -out $wdir/history.gz \
    -maxlevel $maxlevel -keep 3 -skip 0 -target_density $target_density
