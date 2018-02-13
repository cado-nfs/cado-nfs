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

# create las command line from environment variables, moan if any is
# missing.
args=()
for var in poly fb I lim{0,1} lpb{0,1} mfb{0,1} ; do
    args=("${args[@]}" -$var $(eval "echo \${$var:?missing}"))
done

for var in fbc lambda{0,1} ncurves{0,1} descent_hint bkmult bkthresh{,1} ; do
    # Those are optional
    value=$(eval "echo \${$var}")
    if [ "$value" ] ; then
        args=("${args[@]}" -$var "$value")
    fi
done

for var in fbc batch{0,1} ; do
    # Those are optional too. Being filenames, we allow that they be
    # passed as just ".", which means that we expect to have them in the
    # work directory.
    value=$(eval "echo \${$var}")
    if [ "$value" ] ; then
        if [ "$value" = "." ] ; then
            value="${WORKDIR}/${BASENAME}.$var"
            eval "$var=\"\$value\""
        fi
        args=("${args[@]}" -$var "$value")
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
if [ "$fbc" ] || [ "$batch0" ] || [ "$batch1" ] ; then
    run "$LAS_BINARY" "${args[@]}" "${zero_qs[@]}" -out "${RELS}" "$@"
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
    if ! tail "$RELS" | grep "Total $REL_COUNT reports" ; then
        echo "Expected $REL_COUNT reports, got: `tail -1 $RELS`" >&2
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
  if [ -f "${CHECKSUM_FILE}" ] ; then
    # File with checksums already exists, compare
    if diff -b "${CHECKSUM_FILE}" "${MYCHECKSUM_FILE}" > /dev/null ; then
      echo "Checksums agree"
    else
      echo "Error, reference checksums in ${CHECKSUM_FILE} differ from mine in ${MYCHECKSUM_FILE}" >&2
      exit 1
    fi
    let checks_passed+=1
  else
    # File with checksums does not exists, create it
    cp "${MYCHECKSUM_FILE}" "${CHECKSUM_FILE}"
    echo "Created checksum file"
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

if [ "$PHONY_CHECK" ] ; then
  let checks_passed+=1
fi

if [ "$checks_passed" = 0 ] ; then
    echo "Error, zero checks done !" >&2
    exit 1
fi
