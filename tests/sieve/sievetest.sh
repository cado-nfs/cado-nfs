#!/usr/bin/env bash

# This is a generic test harness for las. We recognize several
# environment variables, passed via cmake, that are used as arguments to
# las.
#
# $LAS_BINARY is the actual las binary that gets used.
#
# data is put in $WORKDIR.
#
# Any error fails the whole sript.
#
# After relations are created, post-mortem checks are done according to
# the REFERENCE_SHA1 and REGEX values. See the end of this script.

if [ "$CADO_DEBUG" ] ; then
    set -x
fi
set -e

# Print and run a command
run() {
  echo "Running: $*"
  "$@"
}

if ! [ -x "${LAS_BINARY:?missing}" ] ; then
  echo "Las binary \$LAS_BINARY = ${LAS_BINARY} is not an executable file" >&2
  exit 1
fi

if ! [ -d "${WORKDIR:?missing}" ] ; then
    echo "\$WORKDIR = $WORKDIR is not a directory" >&2
    exit
fi

BASENAME="`basename "${poly:?missing}"`"
BASENAME="${BASENAME%%.*}"

# Make temp directory world-readable for easier debugging
chmod a+rx "${WORKDIR}"

# We may put RELS in a special place.
: ${RELS="${WORKDIR}/${BASENAME}.rels"}

args=()

forwardable=(
    fbc
    lambda{0,1}
    ncurves{0,1}
    descent_hint_table
    bkmult
    bkthresh{,1}
    hint_table
    adjust_strategy
    trace{ab,Nx,ij}
    A B I
)
for var in "${forwardable[@]}" ; do
    # Those are optional
    if [ "${#var}" -gt 1 ] && [ "${!var}" ] ; then args+=(-$var "${!var}") ; fi
    # also recognize via las_$var, because it's a bit ugly to depend on
    # what the external shell variables A, B, and I might be set to.
    # Well, it's ugly to rely on shell variables altogether, of course.
    xvar="las_$var"
    if [ "${!xvar}" ] ; then
        args+=(-$var "${!xvar}")
        eval "$var=${!xvar}"
    fi
    if [[ $var =~ ^trace ]] ; then : ${TRACE="${WORKDIR}/${BASENAME}.trace"} ; fi
done
if [ "$TRACE" ] && [[ $LAS_BINARY =~ tracek ]] ; then
    args+=(-traceout "$TRACE")
fi

# create las command line from environment variables, moan if any is
# missing.
for var in poly lim{0,1} ; do
    args+=(-$var "${!var:?missing}")
done

if ! [ "$hint_table" ] ; then
    # The ones below are not be needed at all if we have a
    # hint_table (unless we use -t auto).
    for var in lpb{0,1} mfb{0,1} ; do
        args+=(-$var "${!var:?missing}")
    done
    if ! [ "$A$I" ] ; then
        echo "missing A and I" >&2
        exit 1
    fi
fi

if [ "$fb0" ] ; then args+=(-fb0 "$fb0") ; fi
if [ "$fb1" ] ; then args+=(-fb1 "$fb1") ; fi
if ! [ "$fb0" ] && ! [ "$fb1" ] ; then
    echo "neither fb nor fb0/fb1 provided" >&2 ; exit 1
fi

for var in fbc batchfile{0,1} ; do
    # Those are optional too. Being filenames, we allow that they be
    # passed as just ".", which means that we expect to have them in the
    # work directory.
    value="${!var}"
    if [ "$value" ] ; then
        if [ "$value" = "." ] ; then
            value="${WORKDIR}/${BASENAME}.$var"
            eval "$var=\"\$value\""
        fi
        args+=(-$var "$value")
    fi
done

if [ "$todo" ] ; then
    end=(-todo "${todo}")
    zero_qs=(-todo /dev/null)
elif [ "$rho" ] ; then
    end=(-q0 "${q0:?missing}" -rho "$rho")
    zero_qs=(-q0 "${q0:?missing}" -q1 "$q0")
elif [ "$nq" ] ; then
    end=(-q0 "${q0:?missing}" -nq "$nq")
    zero_qs=(-q0 "${q0:?missing}" -q1 "$q0")
elif [ "$q1" ] ; then
    end=(-q0 "${q0:?missing}" -q1 "$q1")
    zero_qs=(-q0 "${q0:?missing}" -q1 "$q0")
else
    echo "q1 or rho: missing" >&2
    exit 1
fi

# Warm up the cache files if needed
if [ "$fbc" ] || [ "$batchfile0" ] || [ "$batchfile1" ] ; then
    run "$LAS_BINARY" "${args[@]}" "${zero_qs[@]}" -out "${RELS}" "$@"
fi

if [ "$fbc" ] ; then
    # We should have created a cache file now. Use our companion script
    # to parse the file headers. This can serve as an automated check
    # that the companion script and the source code are kept in sync.
    if [ "$fbc" = "." ] ; then
        real_fbc="${WORKDIR}/${BASENAME}.fbc"
    else
        real_fbc="$fbc"
    fi
    # Note that fbc might not be supported by the current platform.
    # We don't want the test to fail in that case.
    if [ -f "$real_fbc" ] ; then
        "${CADO_NFS_SOURCE_DIR}/sieve/inspect-fbc-file.pl" -fbc "$real_fbc" > "$real_fbc.txt"
    fi
fi

if [ "$batchfile0" ] ; then
    # We should have created a cache file now. Use our companion script
    # to parse the file headers. This can serve as an automated check
    # that the companion script and the source code are kept in sync.
    if [ "$file" = "." ] ; then
        real_file="${WORKDIR}/${BASENAME}.batch0"
    else
        real_file="$batchfile0"
    fi
    "${CADO_NFS_SOURCE_DIR}/sieve/inspect-batch-file.pl" -batch "$real_file" > "$real_file.txt"
fi

if [ "$batchfile1" ] ; then
    # We should have created a cache file now. Use our companion script
    # to parse the file headers. This can serve as an automated check
    # that the companion script and the source code are kept in sync.
    if [ "$file" = "." ] ; then
        real_file="${WORKDIR}/${BASENAME}.batch1"
    else
        real_file="$batchfile1"
    fi
    "${CADO_NFS_SOURCE_DIR}/sieve/inspect-batch-file.pl" -batch "$real_file" > "$real_file.txt"
fi

# then use the cache file created above
run "$LAS_BINARY" "${args[@]}" "${end[@]}" -out "${RELS}" "$@"


### Now do the final checks.

# We depend on several environment variables, presumably defined in the
# caller scripts.
#
# RELS : a relation file to check
# REFERENCE_SHA1 : the sha1sum we should get for the relations
# REFERENCE_REVISION : the git rev that produces REFERENCE_SHA1
#
# environment variables that trigger optional checks:
#
# CHECKSUM_FILE : file that should hold a memory of the sieve regions
#   checksums. The corresponding check is not done if this variable isn't
#   defined.
# REGEX : regular expression that should match in the relation files.

checks_passed=0

if [ "$REL_COUNT" ] ; then
    if ! tail -n 100 "$RELS" | grep "Total $REL_COUNT reports" ; then
        echo "Expected $REL_COUNT reports, got: `tail -n 100 $RELS | grep 'Total.*reports'`" >&2
        exit 1
    fi
    let checks_passed+=1
fi

if [ "$REFERENCE_SHA1" ] ; then
    SHA1BIN=sha1sum
    if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=sha1 ; fi
    if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=shasum ; fi
    if ! type -p "$SHA1BIN" > /dev/null ; then
        echo "Could not find a SHA-1 checksumming binary !" >&2
        exit 1
    fi

    # This was used so that the unix sort -n produce some well-defined
    # ordering on the integers.
    # Now that we sort primes as well, we're doing it in perl anyway.
    # Setting the locale to C can't hurt, but it's perhaps less
    # important.
    export LC_ALL=C
    export LANG=C
    export LANGUAGE=C

    sort_rels() {
        read -s -r -d '' perl_code <<- 'EOF'
            /^[^#]/ or next;
            chomp($_);
            my ($ab,@sides) = split(":", $_, 3);
            for (@sides) {
                $_ = join(",", sort({ hex($a) <=> hex($b) } split(",",$_)));
            }
            print join(":", ($ab, @sides)), "\n";
            
EOF
        perl -ne "$perl_code" "$@" | sort -n
    }

    SHA1=`grep "^[^#]" "${RELS}" | sort_rels | ${SHA1BIN}` || exit 1
    SHA1="${SHA1%% *}"
    echo "$0: Got SHA1 of ${SHA1}"
    echo "$0: expected ${REFERENCE_SHA1}"
    if [ "${SHA1}" != "${REFERENCE_SHA1}" ] ; then
      if [ -n "${REFERENCE_REVISION}" ] ; then
        REFMSG=", as created by Git revision ${REFERENCE_REVISION}"
      fi
      if [ "$CADO_DEBUG" ] ; then
          REFMSG=". Files remain in ${WORKDIR}"
      else
          REFMSG=". Set CADO_DEBUG=1 to examine log output"
      fi
      echo "$0: Got SHA1(sort(${RELS}))=${SHA1} but expected ${REFERENCE_SHA1}${REFMSG}"
      exit 1
    fi
    let checks_passed+=1
fi

if [ -n "${CHECKSUM_FILE}" ] ; then
  MYCHECKSUM_FILE="${WORKDIR}/$(basename "$CHECKSUM_FILE")"
  grep "# Checksums over sieve region:" "${RELS}" > "${MYCHECKSUM_FILE}"
  # copy first, so that it's atomic, and if that happens to fail, then
  # compare the contents with what we have. If we do otherwise, we risk
  # concurrency problems with the three checks
  # F9_quick_{tracek,lambda0,lambda1} which are possibly run in parallel.
  if cp "${MYCHECKSUM_FILE}" "${CHECKSUM_FILE}" 2>/dev/null ; then
    # File with checksums does not exists, create it
    echo "Created checksum file"
  elif [ -f "$CHECKSUM_FILE" ] ; then
    # File with checksums already exists, compare
    if diff -b "${CHECKSUM_FILE}" "${MYCHECKSUM_FILE}" > /dev/null ; then
      echo "Checksums agree"
    else
      echo "Error, reference checksums in ${CHECKSUM_FILE} differ from mine in ${MYCHECKSUM_FILE}" >&2
      exit 1
    fi
    let checks_passed+=1
  else
    echo "Cannot create $CHECKSUM_FILE, this is probably an error" >&2
    exit 1
  fi
fi

if [ -n "${REGEX}" ] ; then
  echo "Searching for regex \"${REGEX}\"" >&2
  if ! grep "${REGEX}" "${RELS}" >&2 ; then
    echo "Error, regular expression \"${REGEX}\" does not match output file"
    exit 1
  fi
  let checks_passed+=1
fi

if [ -n "${TRACE}" ] && [ -n "${TRACE_REGEX}" ] ; then
  echo "Searching for regex \"${TRACE_REGEX}\" on trace file" >&2
  if ! grep "${TRACE_REGEX}" "${TRACE}" >&2 ; then
    echo "Error, regular expression \"${TRACE_REGEX}\" does not match trace file"
    exit 1
  fi
  let checks_passed+=1
fi

if [ "$PHONY_CHECK" ] ; then
  let checks_passed+=1
fi

if [ "$checks_passed" = 0 ] ; then
    echo "Error, zero checks done !" >&2
    exit 1
fi
