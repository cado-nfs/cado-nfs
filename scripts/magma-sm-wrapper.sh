#!/usr/bin/env bash
#
# OUTPUT: SM's for side 0..1 in that order, whatever SM's or units are used,
#         are put in file $OUT, each line corresponding to a relation-set
#

unset POLY
unset RENUMBER
unset BADIDEALINFO
unset PURGED
unset INDEX
unset OUT
unset ELL
unset SMEXP0
unset SMEXP1
unset SMEXP2
unset NMAPS0
unset NMAPS1
unset NMAPS2
unset NMAPS
unset ABUNITS
unset mtopts

EXPLICIT0="no"
EXPLICIT1="no"
EXPLICIT2="no"

NSIDES=0

while [ -n "$1" ]
do
  if [ "$1" = "-poly" ]
  then
    POLY="$2"
    shift 2
  elif [ "$1" = "-renumber" ]
  then
    RENUMBER="$2"
    shift 2
  elif [ "$1" = "-badidealinfo" ]
  then
    BADIDEALINFO="$2"
    shift 2
  elif [ "$1" = "-purged" ]
  then
    PURGED="$2"
    shift 2
  elif [ "$1" = "-index" ]
  then
    INDEX="$2"
    shift 2
  elif [ "$1" = "-out" ]
  then
    OUT="$2"
    shift 2
  elif [ "$1" = "-gorder" ]
  then
    ELL="$2"
    shift 2
  elif [ "$1" = "-nsm0" ]
  then
    NMAPS0="$2"
    NSIDES=`expr $NSIDES '+' 1`
    shift 2
  elif [ "$1" = "-nsm1" ]
  then
    NMAPS1="$2"
    NSIDES=`expr $NSIDES '+' 1`
    shift 2
  elif [ "$1" = "-nsm2" ]
  then
    NMAPS2="$2"
    NSIDES=`expr $NSIDES '+' 1`
    shift 2
  elif [ "$1" = "-abunits" ]
  then
    ABUNITS="$2"
    shift 2
  elif [ "$1" = "-t" ]
  then
    mtopts="-t $2"
    shift 2
  elif [ "$1" = "-explicit_units0" ]
  then
    EXPLICIT0="yes"
    shift 1
  elif [ "$1" = "-explicit_units1" ]
  then
    EXPLICIT1="yes"
    shift 1
  elif [ "$1" = "-explicit_units2" ]
  then
    EXPLICIT2="yes"
    shift 1
  else
    echo "Unknown parameter: $1"
    exit 1
  fi
done

# where am I ?
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

## FIXME: units and NSIDES > 2
## set parameters depending on units or not: nsm0 (resp. nsm1) will be the
## number of true SM's to be computed (not the number of units)
side=-1
if [ $EXPLICIT0 = "yes" ]; then
    nsm0=0
    side=0
    usef=false
else
    nsm0=$NMAPS0
fi
nsm="-nsm $nsm0"
if [ $EXPLICIT1 = "yes" ]; then
    nsm1=0
    side=1
    usef=true
else
    nsm1=$NMAPS1
fi
nsm="$nsm,$nsm1"
## FIXME: copy-paste without thinking!
if [ $NSIDES -gt 2 ]; then
    if [ $EXPLICIT2 = "yes" ]; then
	nsm2=0
	side=2
	usef=true
    else
	nsm2=$NMAPS2
    fi
    nsm="$nsm,$nsm2"
fi

## build required SM's
if [ $NMAPS0 -gt 0 -o $NMAPS1 -gt 0 ]; then
    # at least one side requires SM's
    if [ -s $OUT ]; then
	echo "File $OUT already exists"
    else
	# out contains SM's for side=0..1 in that order, or 0 or 1 depending
	# on the unit side if any
	CMD="$DIR/../filter/sm -poly $POLY -purged $PURGED -index $INDEX -out $OUT -gorder $ELL $smopts $nsm $mtopts"
	if [ $nsm0 -ne 0 -o $nsm1 -ne 0 ]; then
	    echo $CMD; $CMD
	else
	    echo "WARNING: no SM to compute at all"
	fi
    fi
else
    echo "No map to be computed"
fi

## now, use units if needed: we suppose that we use either 0 or 1, not both
if [ $side -ne -1 ]; then
    # operates in 3 steps:
    # 1. extract (p, r) ideals from $RENUMBER and put them in ...algpr.gz
    # 2. from ...algpr.gz, compute a generator for each ideal and put it in
    #    ...generators.gz
    # 3. compute the unit contribution for each relation set from $INDEX
    #
    # algpr.$side.gz and generators.$side.gz should be computed only once and
    # can be done in the meantime for instance.
    #
    # TODO: it might be possible to skip the algpr file (together with the
    # debug_renumber prgm) by rewriting / using this more cleverly in C.

    prgm=$DIR/../misc/debug_renumber

    if [ ! -s $prgm ]; then echo "Please compile $prgm"; exit; fi

    # guess names
    algpr=`echo $RENUMBER | sed "s/freerel.renumber/algpr.$side/"`
    generators=`echo $RENUMBER | sed "s/freerel.renumber/generators.$side/"`
    relsdel=`echo $PURGED | sed 's/purged/relsdel/'`
    abunits=`echo $RENUMBER | sed "s/freerel.renumber.gz/units.abunits.$side/"`

    # extract (p, r)
    if [ ! -s $algpr ]; then
	echo "Building file $algpr for side $side using debug_renumber"
	# we have to make it work for rat/alg or side 0/1
	$prgm -poly $POLY -renumber $RENUMBER |\
        egrep "(alg side|side $side)" | sed 's/=/ /g' |\
        awk '{print $2, $6, $8}' |\
        gzip -c > $algpr
    fi
    # find generators for (p, r)
    if [ ! -s $generators ]; then
	echo "Building file $generators using Magma for side $side"
	CMD="magma -b polyfile:=$POLY renumber:=$RENUMBER algpr:=$algpr generators:=$generators badidealinfo:=$BADIDEALINFO purged:=$PURGED relsdel:=$relsdel index:=$INDEX ficunits:=$OUT abunits:=$abunits usef:=$usef ww:=true $DIR/nfsunits.mag"
	echo $CMD; $CMD
    fi

    if [ -s $OUT.$side ]; then
	echo "File $OUT.$side already exists"
    else
	echo "## Using magma to compute units for side $side"
        # results are put in $OUT.$side
	CMD="magma -b polyfile:=$POLY renumber:=$RENUMBER algpr:=$algpr generators:=$generators badidealinfo:=$BADIDEALINFO purged:=$PURGED relsdel:=$relsdel index:=$INDEX ficunits:=$OUT.$side abunits:=$abunits usef:=$usef ww:=false $DIR/nfsunits.mag"
	echo $CMD; $CMD

	# split abunits
	# typicaly ABUNITS="...../p3dd15-f4g3-GJL-1.sm.abunits.dir"
	abunitsdir=$ABUNITS".$side"
	if [ -d $abunitsdir ]; then
	    echo "Directory $abunitsdir already exists"
	else
	    echo "I create directory $abunitsdir"
	    sort -k1n -k2n $abunits > /tmp/foo$$
	    mkdir $abunitsdir
	    echo "I split /tmp/foo$$ into files + index"
	    $DIR/split.py /tmp/foo$$ $abunitsdir/ 100
	    /bin/rm -f /tmp/foo$$
	fi
    fi
    # paste files if needed: we require sides 0..1 in that order
    if [ $side -eq 0 ]; then
	f0=$OUT.$side; f1=$OUT
	if [ $nsm1 -ne 0 ]; then
	    # paste $OUT.0 and $OUT to get a new $OUT
	    dopaste=true
	else
            # nsm1 = 0 means no SM on side 1
	    head -1 $f0
	    cat $f0 | (read sr0 nsm0 p0; echo "$sr0 $nsm0 $ELL"; cat) > $OUT.$$
	    mv $OUT.$$ $OUT
	    dopaste=false
	fi
    else
	f0=$OUT; f1=$OUT.$side
	if [ $nsm0 -ne 0 ]; then
	    # paste $OUT and $OUT.1 to get a new $OUT
	    dopaste=true
	else
	    # nsm0 = 0 means no SM on side 0
	    head -1 $f1
	    cat $f1 | (read sr1 nsm1 p1; echo "$sr1 $nsm1 $ELL"; cat) > $OUT.$$
	    mv $OUT.$$ $OUT
	    dopaste=false
	fi
    fi
    if $dopaste; then
	# first line is "number_of_lines_of_sms nmaps0 1"
	# we need to output: "number_of_lines_of_sms (nmaps0+nmaps1) ell"
        head -1 $f0
        head -1 $f1
	paste -d" " $f0 $f1 | (read sr0 nsm0 p0 sr1 nsm1 p1 ; echo "$sr0 $((nsm0+nsm1)) $ELL" ; cat) > $OUT.$$
	mv $OUT.$$ $OUT
    fi
fi
