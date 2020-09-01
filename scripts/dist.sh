#!/bin/sh

# Usage: scripts/dist.sh <package basename>
# This creates <package basename>.tar.gz

set -e

pkg=$1
shift

t=`mktemp -d ${TMPDIR-/tmp}/cado-nfs.XXXXXXXXXX`

if [ "$CHECK_SCM_FILES" ] ; then
    # check whether we have everything the SCM system knows about in our
    # files.* lists
    if ! scripts/check_file_lists.pl ; then
        echo "Refusing to make the distribution. Sorry" >&2
        exit 1
    fi
fi

read_files_dist() {
    # This command reads the files.dist list. All nodes ending in / imply
    # that the whole subdirectory is matched.
    grep '^[^#]' files.dist | while read x ; do
        case $x in
            */) find $x ;;
            *) echo $x;;
        esac
    done | sed -e 's+//*+/+g' -e 's+^+/+'
    # On some machines, find foo/ will produce foo//file, which does not
    # play well with rsync.
}

read_files_dist | rsync -a     \
                    --include-from=-    \
                    --include='*/'      \
                    --exclude='*'       \
                    --delete-excluded   \
                    --prune-empty-dirs  \
                    ./ $t/$pkg/

scripts/version.sh > $t/$pkg/.git_version

if tar --version 2>/dev/null | grep -q "GNU tar" ; then
	if tar --version | perl -ne '/^tar.*(\d+\.\d+)/ && do { $ver=$1; }; END { exit 1 unless defined($ver) && $ver >= 1.13; exit 0; }' ; then
		taropts="--format=oldgnu"
	fi
fi
(cd $t ; tar -c $taropts -f - $pkg) | gzip -9 > $pkg.tar.gz

echo "Built $pkg.tar.gz"
rm -rf $t
