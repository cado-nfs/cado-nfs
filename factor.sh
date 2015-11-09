#!/usr/bin/env bash

# Run cadofactor.py on the current machine.
# This script is recommended only for small factorizations (less than 100
# digits). For larger numbers, you should use the cadofactor.py script
# (see the documentation in that script).
#
# All partial data are put in a temp dir that is removed at the end in
# case of success. The temp dir is left here in case of failure. If
# CADO_DEBUG is set, the temp dir is never removed.
# If the environment variable $t exists, this sets the directory to use.

usage() {
    cat >&2 <<EOF
Usage: $0 <integer> [options] [arguments passed to cadofactor.py]
    <integer>     - integer to be factored (must be at least 58 digits)

Options:
    -dlp          - don't factor, solve dlp in GF(integer)
    -gfpext <integer>
                  - (with DLP) solve dlp in GF(integer^gfpext)
    -ell <integer>
                  - (with DLP) solve dlp modulo ell
    -t <integer>  - numbers of cores to be used, or 'auto' (see below)
    -s <integer>  - numbers of slave scripts to be used, or 'auto' (see below)
    -h <host1>,<host2>,...,<hostn>
                  - comma-separated list of hosts to use for factorization.
                    With -s <n>, run <n> many scripts per host.
                    Default: localhost
    -dlpnokeep    - disable the feature that CADO_DEBUG is set by default
                    in dlp mode

    If the shell environment variable CADO_DEBUG is set to any non-empty
    string, the working directory is not deleted even if the factorization
    succeeds.

    If CADO_USEHOST is set to any non-empty string, the script will
    listen on "*" even if all clients given with -h are on localhost.

Use of 'auto' for -t and -s:
    With -t auto, the various multithreaded programs in cado-nfs are run
    with as many threads as we have available cores on the running
    machine (the value is guessed by shell commands within this script).
    This affects the tasks.threads parameter of cadofactor. The meaning
    of the -t switch may be changed by -s auto (see below).

    With -s auto (which *must* come with -t <some integer>, not -t auto),
    the las and polyselect programs are capped to the number of threads
    given by -t, however several slaves are run, up to the number which
    keeps all the cores busy. For other multithreaded programs, we use as
    many threads as we have available cores.
EOF
}

# Set "t" to a temp directory, if not already set.
# The ":" shell built-in causes variable expansion but nothing else to happen
: ${t:=`mktemp -d /tmp/cado.XXXXXXXXXX`}
if ! [ -d "$t" ] ; then
    echo "No directory $t ???" >&2
    exit 1
fi
# Make temp dir writable by us and readable by everyone
chmod 755 $t

n=
cores=1
slaves=1
verbose=false
dlp=false
dlpkeep=true
gfpext=1
ell=1
cadofactor_args=()
is_number() { [ "$(grep "^[ [:digit:] ]*$" <<< "$1")" ] ; }
want_number() {
    if ! is_number "$3" ; then
        echo "Error, argument to $1 ($2) must be an integer" >&2
        exit 1
    fi
}


while [ "$#" -gt 0 ] ; do
    if ! [ "$n" ] && is_number $1 ; then
        n=$1
        shift
    elif [ "$1" = "-t" ]; then
        cores=$2
        [ "$2" = "auto" ] || want_number -t "number or cores" $2
        shift 2
    elif [ "$1" = "-s" ]; then
        slaves=$2
        [ "$2" = "auto" ] || want_number -s "number or slaves" $2
        shift 2
    elif [ "$1" = "-h" ]; then
        hostnames="$2"
        shift 2
    elif [ "$1" = "-dlp" ]; then
        dlp=true
        shift
    elif [ "$1" = "-gfpext" ]; then
        gfpext="$2"
        if [ $gfpext != "2" ]; then 
            echo "Error, only extensions of degree 2 are supported for now" >&2
            usage
            exit 1
        fi
        shift 2
    elif [ "$1" = "-ell" ]; then
        ell=$2
        want_number -ell "DLP subgroup order" $2
        shift 2
    elif [ "$1" = "-dlpnokeep" ]; then
        dlpkeep=false
        shift
    elif [ "$1" = "-v" ]; then
        verbose=true
        shift
    else
        cadofactor_args=("${cadofactor_args[@]}" "$1")
        shift
    fi
done

if ! [ "$n" ] ; then
  echo "Error, must give number to factor (or p modulo which dlp is to be solved)" >&2
  usage
  exit 1
fi

# In DLP mode, keep data, unless -dlpnokeep has been passed.
if $dlp; then
  if $dlpkeep; then
    CADO_DEBUG="Yes, please"
  fi
fi

# If gfpext is 2, then ell must be given
if [ $gfpext = "2" ]; then
  if [ $ell = "1" ]; then
    echo "Error, with gfpext=2, ell (a factor of p+1) must be given" >&2
    usage
    exit 1
  fi
fi

# Set default value
: ${hostnames:=localhost}

#########################################################################
# Set paths properly.

# Try to find out from which location this script is running:
# the install directory (i.e., where 'make install' put it),
# the build directory (where 'cmake' put it, substituting the 
# CMake variables below), or in the source directory tree.
cado_prefix="@CMAKE_INSTALL_PREFIX@"
cado_source_dir="@CADO_NFS_SOURCE_DIR@"
cado_build_dir="@CADO_NFS_BINARY_DIR@"
mpiexec="@MPIEXEC@"
archive_name="@CADO_DIST_ARCHIVE_NAME@"


if [ "$0" -ef "$cado_prefix/lib/$archive_name/factor.sh" ] ; then
    # We're called in the install tree.
    bindir="$cado_prefix/lib/$archive_name"
    scriptpath="$bindir"
    cadofactor="$scriptpath/cadofactor.py"
    paramdir="$cado_prefix/share/$archive_name"
    cpubinding_file="$cado_prefix/share/$archive_name/misc/cpubinding.conf"
    if ! [ -f "$cpubinding_file" ] ; then cpubinding_file= ; fi
elif [ "$0" -ef "$cado_prefix/bin/factor.sh" ] ; then
    # We're called in the install tree.
    bindir="$cado_prefix/bin"
    scriptpath="$bindir"
    cadofactor="$scriptpath/cadofactor.py"
    paramdir="$cado_prefix/share/$archive_name"
    cpubinding_file="$cado_prefix/share/$archive_name/misc/cpubinding.conf"
    if ! [ -f "$cpubinding_file" ] ; then cpubinding_file= ; fi
elif [ "$0" -ef "$cado_build_dir/factor.sh" ] ; then
    # We're called in the build tree.
    scriptpath="$cado_source_dir/scripts/cadofactor"
    paramdir="$cado_source_dir/params"
    cadofactor="$scriptpath/cadofactor.py"
    # Make the path absolute.
    bindir="$cado_build_dir"
    cpubinding_file="$cado_source_dir/linalg/bwc/cpubinding.conf"
elif [ -f "`dirname $0`/cado_config_h.in" ] ; then
    # Otherwise we're called from the source tree (or we hope so)
    srcdir=$(cd "`dirname $0`" ; pwd)
    call_cmake="`dirname $0`/scripts/call_cmake.sh"
    unset build_tree
    if [ -x "$call_cmake" ] ; then
      # This should set $build_tree:
      eval `cd $srcdir ; $call_cmake show`
    fi
    if ! [ -z "$build_tree" ] ; then
      scriptpath="${srcdir}/scripts/cadofactor"
      cadofactor="${scriptpath}/cadofactor.py"
      paramdir="${srcdir}/params"
      # Make the path absolute.
      bindir=`cd "$build_tree" ; pwd`
      cpubinding_file="$srcdir/linalg/bwc/cpubinding.conf"
    fi
fi
if [ -z "$cadofactor" ] ; then
    echo "I don't know where I am !" >&2
    # but I do care!
    exit 1
fi

if $dlp ; then
    paramdir="${paramdir}_dl"
fi

if ! [ -d "$paramdir" ] ; then
    echo "Parameter dir $paramdir not found." >&2 ; exit 1
elif ! [ -d "$bindir" ] ; then
    echo "Binary dir $bindir not found." >&2 ; exit 1
elif ! [ -x "$cadofactor" ] ; then
    echo "Script $cadofactor not found." >&2 ; exit 1
else
    # Ok, everything looks good.
    :
fi

size=${#n}

# we round to the nearest multiple of 5
size=`expr \( \( $size + 2 \) \/ 5 \) \* 5`

if $dlp ; then
    if [ $gfpext = "1" ]; then
      file="$paramdir/params.p$size"
      export CADO_HINTFILE_DIR=$paramdir
    else
      file="$paramdir/params.p${gfpext}dd$size"
    fi
else
    file="$paramdir/params.c$size"
fi
echo $file

if [ ! -f $file ] ; then
    echo "no parameter file found for c${#n} (tried $file)" >&2
    exit 1
fi

ncpus() {
    if [ "$NCPUS_FAKE" ] ; then
        echo $NCPUS_FAKE
        exit 0
    fi

    if [ -f /proc/cpuinfo ] ; then
        # grep -c fails with no match, while |wc -l is happy
        nphysical=$(sort -u /proc/cpuinfo  | grep '^physical' | wc -l)
        if [ "$nphysical" -eq 0 ] ; then
            grep -c ^processor /proc/cpuinfo
            exit 0
        fi
        variants="$(grep '^cpu cores' /proc/cpuinfo | uniq | wc -l)"
        if [ "$variants" = 0 ] ; then
            grep -c ^processor /proc/cpuinfo
            exit 0
        fi
        if [ "$variants" != 1 ] ; then
            echo "inhomogeneous platform ?" >&2
            exit 1
        fi
        cores_per_cpu=$(grep '^cpu cores' /proc/cpuinfo | head -1 | cut -d: -f2)
        echo $((nphysical*cores_per_cpu))
    elif [ "$(uname -s)" = Darwin ] ; then
        # does this count hyperthreading or not ?
        sysctl -n hw.ncpu
    elif [ "$(uname -s)" = OpenBSD ] ; then
        # does this count hyperthreading or not ?
        sysctl -n hw.ncpu
    elif [ "$(uname -s)" = MINGW32_NT-6.1 ] ; then
        # not clear whether it's physical or logical.
        wmic cpu get Caption | tail -n +2 | grep -c .
        # we don't believe mingw will support multithreading well.
    else
        # this would work as well on linux and darwin, but pretty surely does
        # not count hyperthreading. Does not work on openbsd5.3
        getconf _NPROCESSORS_ONLN || (echo "Selecting 1 core only" >&2 ; echo 1)
    fi
}

cores_args=()
if [ "$cores" = auto ] && [ "$slaves" = auto ] ; then
        echo "Having both -s auto and -t auto is forbidden" >&2
        exit 1
elif [ "$cores" = auto ] ; then
    cores=`ncpus`
    cores_args=("tasks.threads=$cores" slaves.nrclients=$slaves)
elif [ "$slaves" = auto ] ; then
    ncpus=`ncpus`
    if [ "$cores" -gt "$ncpus" ] ; then
        cores=$ncpus
    fi
    slaves=$((($cores-1+$ncpus)/$cores))
    cores_args=("tasks.threads=$ncpus" "tasks.polyselect.threads=$cores" "tasks.sieve.threads=$cores" "slaves.nrclients=$slaves")
else
    cores_args=("tasks.threads=$cores" "slaves.nrclients=$slaves")
fi

#########################################################################
[ "$CADO_DEBUG" ] && echo "(debug mode, temporary files will be kept in $t)"

# This copy is not strictly necessary (it used to be when we used the
# fixup_params_file.sh script), but it might not be a bad idea to have the
# parameter file together with the rest in the working directory, which also
# allows for easy modification of parameters for the specific factorization.
cp $file $t/param

# Sets the machine description file
if [ -z "$CADO_USEHOST" ] && echo "$hostnames" | grep -E "^(localhost, *)*localhost$" > /dev/null
then
   server_address="server.address=localhost" # For the Python script
else
   server_address=""
fi

mkdir $t/tmp

args=(
    tasks.execpath="$bindir"
    "${cores_args[@]}"
    tasks.workdir="$t"
    slaves.hostnames="$hostnames"
    slaves.scriptpath="$scriptpath"
    "$server_address"
    slaves.basepath="$t/client/"
    "${cadofactor_args[@]}"
)

if [ "$cpubinding_file" ] ; then
    args=("${args[@]}" tasks.linalg.bwc.cpubinding="$cpubinding_file")
fi

if [ $ell != "1" ]; then
    args=("${args[@]}" ell=$ell)
fi

# $PYTHON is there to expand a shell variable having that name, if
# provided, in the case there is no python3 script in the path, or if
# one which is named otherwise, or placed in a non-prority location
# is the path, is preferred. If $PYTHON is empty, this is a no-op

$PYTHON $cadofactor "$t/param" N=$n "${args[@]}"

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
