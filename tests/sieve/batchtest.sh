#!/usr/bin/env bash

# Print and run a command
function run {
  echo "Running: $*"
  "$@"
}

# Input files. First parameter is the polynomial file to use, second is root of the build directory, third is the expected SHA1 value
FB="$1"
LAS="$2"
POLY="$3"
REFERENCE_SHA1="$4"
REFERENCE_REVISION="$5"
CHECKSUM_FILE="$6"
shift 6

if [[ -z "${FB}" || ! -f "${FB}" ]]
then
  echo "Factor base ${FB} is not a file" >&2
  exit 1
fi

if [[ -z "${LAS}" || ! -x "${LAS}" ]]
then
  echo "Las binary ${LAS} is not an executable file" >&2
  exit 1
fi

if [ -z "${REFERENCE_SHA1}" ]
then
  echo "No reference SHA1 value specified" >&2
  exit 1
fi


if [[ -z "$rlim" || -z "$alim" || -z "$lpbr" || -z "$lpba" || -z "$maxbits" \
      || -z "$mfbr" || -z "$mfba" || -z "$rlambda" || -z "$alambda" || -z "$I" \
      || -z "$q0" || ( -z "$q1" && -z "$rho" ) ]]
then
  echo "Required shell environment variable not set" >&2
  exit 1
fi

if [ "$1" = "-regex" ]
then
  REGEX="$2"
  shift 2
fi

BASENAME="`basename "${POLY}"`"
BASENAME="${BASENAME%%.*}"

# Output files
TMPDIR=`mktemp -d /tmp/cadotest.XXXXXXXXXX`
# Make temp direcotry world-readable for easier debugging
chmod a+rx "${TMPDIR}"
RELS="${TMPDIR}/${BASENAME}.rels"

if [ -n "$q1" ]
then
  end=("-q1" "$q1")
fi
if [ -n "$rho" ]
then
  end=("-rho" "$rho")
fi

bailout() {
    if [ "$EXPECTED_FAIL" ] && ! [ "$KEEP_SIEVETEST" ] ; then
        rm -rf "$TMPDIR"
    fi
    exit 1
}

run "$LAS" -poly "$POLY" -fb "${FB}" -I "$I" -rlim "$rlim" -lpbr "$lpbr" -mfbr "$mfbr" -rlambda "$rlambda" -alim "$alim" -lpba "$lpba" -mfba "$mfba" -alambda "$alambda" -q0 "$q0" "${end[@]}" -out "${RELS}" -batch "$@" || bailout


SHA1BIN=sha1sum
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=sha1 ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=shasum ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then
    echo "Could not find a SHA-1 checksumming binary !" >&2
    exit 1
fi

# Try to make sort produce some well-defined ordering on the integers
export LC_ALL=C
export LANG=C
export LANGUAGE=C

SHA1=`grep -v "^#" "${RELS}" | sort -n | ${SHA1BIN}` || exit 1
echo "$0: Got SHA1 of ${SHA1}"
echo "$0: expected ${REFERENCE_SHA1}"
SHA1="${SHA1%% *}"
if [ "${SHA1}" != "${REFERENCE_SHA1}" ]
then
  if [ -n "${REFERENCE_REVISION}" ]
  then
    REFMSG=", as created by Git revision ${REFERENCE_REVISION}"
  fi
  echo "$0: Got SHA1 of ${SHA1} but expected ${REFERENCE_SHA1}${REFMSG}. Files remain in ${TMPDIR}"
  exit 1
fi

if [ -n "${REGEX}" ]
then
  echo "Searching for regex \"${REGEX}\"" >&2
  if ! grep "${REGEX}" "${RELS}" >&2
  then
    echo "Error, regular expression \"${REGEX}\" does not match output file"
    exit 1
  fi
fi

if [ -z "$KEEP_SIEVETEST" ]
then
  rm -f "${RELS}"
  if [ -n "${MYCHECKSUM_FILE}" ]
  then
    rm -f "${MYCHECKSUM_FILE}"
  fi
  rmdir "${TMPDIR}"
else
  echo "Keeping files in ${TMPDIR}"
fi
