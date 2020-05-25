#!/usr/bin/env bash

OUTPUT_FILE="$1"

if ! [ -d .git ] ; then
    if ! [ -f "$OUTPUT_FILE" ] ; then
        cat > "$OUTPUT_FILE" <<'EOF'
#include "cado.h" // IWYU pragma: keep
#include "version_info.h"       // IWYU pragma: keep
const char * cado_modified_files = "# (tarball extracted)\n";
EOF
    fi
    exit 0
fi

SHA1BIN=sha1sum
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=sha1 ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=shasum ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then
    echo "Could not find a SHA-1 checksumming binary !" >&2
    exit 1
fi

function list_modified() {
  cat << 'EOF'
#include "cado.h" // IWYU pragma: keep
#include "version_info.h"       // IWYU pragma: keep
const char * cado_modified_files = 
EOF
  git status --porcelain -uno | while read STATUS FILE
  do
    if [[ "$STATUS" = M  ||  "$STATUS" = A ]]
    then
      SHA1="`$SHA1BIN "$FILE" | cut -d " " -f 1`"
      echo \""# $STATUS $FILE $SHA1"'\n'\"
    else
      echo \""# $STATUS $FILE"'\n'\"
    fi
  done
  echo '"";'
}

TEMPFILE="`mktemp ${TMPDIR-/tmp}/cado-nfs.XXXXXXXXX`"
list_modified >> "$TEMPFILE"
if ! [ -f "$OUTPUT_FILE" ] || ! cmp -s "$TEMPFILE" "$OUTPUT_FILE"
then
  rm -f "$OUTPUT_FILE"
  cp "$TEMPFILE" "$OUTPUT_FILE"
fi
rm -f "$TEMPFILE"
