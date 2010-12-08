#!/usr/bin/env bash

# Run cadofactor.pl on the current machine
#
# All partial data are put in a temp dir that is removed at the end in
# case of success. The temp dir is left here in case of failure. If
# CADO_DEBUG is set, the temp dir is never removed.

# Mac OS X's mktemp requires a prefix.
# And we need $t to be an absolute pathname, or else ./new_cadofactor.pl
# will fail.

#CADO_DEBUG=1
usage="Usage: $0 <integer> [options] [arguments cadofactor.pl]\n
       \t <integer> \t \t - integer must be at least 60 digits [optimize for more than 85 digits] without small prime factors\n
       \t options:\n
       \t \t -t <integer> \t - numbers of cores\n
       \t \t -ssh \t \t - use ssh (doc README) for distribute the polynomial selection and
sieve on localhost. (on several machines, use the main script
cadofactor.pl)"

: ${t:=`mktemp -d /tmp/cado.XXXXXXXXXX`}
chmod 755 $t

n=$1
shift
if [ ! "$(grep "^[ [:digit:] ]*$" <<< $n)" ]; then
  echo -e $usage
  exit 1
fi

ssh=0
cores=1
for ((i=1; i<=2; i=i+1)) ; do
  if [ "$1" = "-t" ]; then
    cores=$2
    if [ ! "$(grep "^[ [:digit:] ]*$" <<< $cores)" ]; then
      echo -e $usage
      exit 1
    fi
    shift 2
  fi
  if [ "$1" = "-ssh" ]; then
    ssh=1
    shift
  fi
done
    
b=$(perl -e '{print int(sqrt(<STDIN>))."\n"}' <<< $cores)
until [ `expr $cores % $b` -eq  0 ]; do
      b=`expr $b - 1`
done
a=`expr $cores / $b`
bwmt=${a}x$b

size=${#n}
file="params.c$size"


#########################################################################
# Set paths properly.
cado_prefix="@CMAKE_INSTALL_PREFIX@"
example_subdir="@example_subdir@"
mpiexec="@MPIEXEC@"

for ((i=1; i<=4; i=i+1)) ; do
  if [ -d "$cado_prefix/$example_subdir" ] ; then
      # We're called in the install tree.
      if [ -f "$cado_prefix/$example_subdir/params.c$size" ] ; then
          file="$cado_prefix/$example_subdir/params.c$size" 
      fi
      bindir="$cado_prefix/bin"
      cadofactor="$bindir/cadofactor.pl"
  elif [ -x "@CADO_NFS_SOURCE_DIR@/cadofactor.pl" ] ; then
      # Otherwise we're called from the source tree (or we hope so)
      if [ -f "@CADO_NFS_SOURCE_DIR@/params/params.c$size" ] ; then
          file="@CADO_NFS_SOURCE_DIR@/params/params.c$size"
      fi
      cadofactor="@CADO_NFS_SOURCE_DIR@/cadofactor.pl"
      # Make the path absolute.
      bindir=$(cd "`dirname $0`" ; pwd)
  else
      # Otherwise we're called from the source tree (or we hope so)
      call_cmake="`dirname $0`/scripts/call_cmake.sh"
      if ! [ -x "$call_cmake" ] ; then
          echo "I don't know where I am !" >&2
      fi
      eval `$call_cmake show`
      if [ -f "${up_path}params/params.c$size" ] ; then
          file="${up_path}params/params.c$size"
      fi
      cadofactor="${up_path}cadofactor.pl"
      # Make the path absolute.
      build_tree=`cd "$build_tree" ; pwd`
      bindir="$build_tree"
  fi
  if [ -f $file ] ; then 
    break
  fi
  size=`expr $size + \( 2 \* \( $i % 2 \) - 1 \) \* $i`
done

if [ ! -f $file ] ; then
    echo "$file not found" >&2
    echo -e $usage
    exit 1
fi

#########################################################################
[ "$CADO_DEBUG" ] && echo "(debug mode, temporary files will be kept in $t)"

cp $file $t/param

$bindir/fixup_params_file.sh $t/param

# Sets the machine description file
if [ -z "$CADO_USEHOST" ] ; then
   host=localhost
else
   host=`hostname --short`
fi

mkdir $t/tmp

cat > $t/mach_desc <<EOF
[local]
tmpdir=$t/tmp
bindir=$bindir
$host cores=$cores
EOF

if [ $ssh -eq 0 ]; then
  $cadofactor $t/param n=$n bindir=$bindir parallel=0 \
  sieve_max_threads=$cores bwmt=$bwmt wdir=$t \
  sievenice=0 selectnice=0 logfile=$t/out "$@"
else
  $cadofactor $t/param n=$n bindir=$bindir parallel=1 \
  machines=$t/mach_desc bwmt=$bwmt wdir=$t \
  sievenice=0 selectnice=0 logfile=$t/out "$@"
fi

rc=$?

#########################################################################
# Check result, clean up the mess afterwards.
if [ "$rc" = 0 ] ; then
    echo OK
    if ! [ "$CADO_DEBUG" ] ; then
        rm -rf $t
    fi
else
    echo "FAILED ; data left in $t"
fi
exit $rc
