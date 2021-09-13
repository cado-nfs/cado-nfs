#!/bin/bash

set -ex

usage() {
    echo "Usage: $0 -o <directory> [info files...]" >&2
    exit 1
}

while [ $# -gt 1 ] ; do
    if [ "$1" = "-o" ] ; then
        shift
        OUTPUTDIR=$1
        shift
    else
        break
    fi
done

: ${OUTPUTDIR?missing}

rev=

for f in "$@" ; do
    f=$f
    f=${f##coverage-}
    f=${f%%-*}
    if ! [ "$rev" ] || [ "$rev" = "$f" ] ; then
        rev="$f"
    else
        echo "Inconsistent git commit across coverage files: $@" >&2
        exit 1
    fi
done

TMPDIR=`mktemp -d /tmp/cado-cov.XXXXXXXXXXXX`
cleanup() {
    rm -rf $TMPDIR
}
trap cleanup EXIT

cp "$@" $TMPDIR
files=($(find $TMPDIR -type f))
ls "${files[@]}" | xargs -n 1 sed -e s,^SF:,SF:$PWD/,g -i

/bin/cp -pf $(which genhtml) $TMPDIR/genhtml
commit=$rev
ex $TMPDIR/genhtml <<EOF
/sub get_date_string()\$
/^{
mark a
/^}
mark b
'a,'b c
{
       return scalar localtime;
}
.
/headerValue.*date
a
        push(@row_left, [[undef, "headerItem", "Commit:"],
                         [undef, "headerValue", "<a href=\"https://gitlab.inria.fr/cado-nfs/cado-nfs/-/commit/$commit\">$commit</a>"]]);
.
wq
EOF
# right now we don't recover the test results.
$TMPDIR/genhtml -o $OUTPUTDIR -p $PWD "${files[@]}"

# now=$(date +%Y%m%d%H%M%S)
# today=$(date +%Y%m%d)
# rm -f $OUTPUTDIR/coverage-${now} || :
# rm -f $OUTPUTDIR/coverage-${today} || :
# rm -f $OUTPUTDIR/coverage-previous || :
# ln -s `readlink $OUTPUTDIR/coverage-latest` $OUTPUTDIR/coverage-previous
# rm -f $OUTPUTDIR/coverage-latest || :
# ln -s coverage-${rev} $OUTPUTDIR/coverage-${now}
# ln -s coverage-${rev} $OUTPUTDIR/coverage-${today}
# ln -s coverage-${rev} $OUTPUTDIR/coverage-latest
# cp $tarfile $OUTPUTDIR/coverage-${rev}.tar.gz
# 
# rsync -av --filter '+ /coverage*' --filter '- /*' --delete  coverage/ $HOME/.webdir/coverage/
