#!/bin/bash

# This script downloads, builds, and runs cado-nfs on several machines.
# The machines from the gcc compile farm are taken as an example. Tests
# are done in parallel on the different nodes. Once all jobs are
# finished, a digest is sent by mail, with a tar.gz attachment containing
# all the results. The output can probably be improved (e.g. give
# warnings in a text attachments, not buried inside a tar.gz).  The
# requirement is that cmake is present on the target machines (as per the
# companion script, here gcc-script.sh).

# The companion script gcc-script.sh is of course required.

# In order to access the machines through ssh, it is mandatory to have
# proper support in  ~/.ssh/config. For the gcc compile farm, e.g.:
# GET http://gcc.gnu.org/wiki/CompileFarm | perl ~thome/.ssh/parse-cfarm.pl

gccs=(\
# (/home/zimmerma not found) # gcc11
# (/home/zimmerma not found) # gcc12
# (/home/zimmerma not found) # gcc13
# (/home/zimmerma not found) # gcc14
# (/home/zimmerma not found) # gcc15
# (/home/zimmerma not found) # gcc16
# (/home/zimmerma not found) # gcc17
gcc40
# gcc42 # disabled 20101202 because the C++ compiler is non-functional
# (ssh timeout) # gcc43
# (ssh timeout) # gcc50
gcc51
# (ssh timeout) # gcc53
gcc54
# (ssh timeout) # gcc55
# (ssh timeout) # gcc56
# (/home/zimmerma not found) # gcc57
gcc60
gcc61
gcc62
gcc63
gcc64
# (ssh timeout) # gcc200
# (ssh timeout) # gcc201
# (ssh timeout) # gcc01
# (ssh timeout) # gcc02
# (ssh timeout) # gcc03
# (ssh timeout) # gcc05
# (ssh timeout) # gcc06
# (ssh timeout) # gcc07
# (ssh timeout) # gcc09
# (/home/zimmerma not found) # gcc08
# (dns failure) gcc30
# (dns failure) gcc31
# (ssh timeout) # gcc41
# (/home/zimmerma not found) # gcc10
gcc33
gcc34
# (ssh timeout) # gcc35
# (ssh timeout) # gcc36
gcc37
gcc38
# (/home/zimmerma not found) # gcc52
# (/home/zimmerma not found) # gcc100
# (/home/zimmerma not found) # gcc101
)

tmpdir=`mktemp -d /tmp/XXXXXXXX`

check_paul_on_one_node() {
    host=$1
    ssh -x -oStrictHostKeyChecking=no $host ls /home/zimmerma/bin 2>&1 | sed -e "s/^/$host:/"
}

run_on_one_node() {
    host=$1
    ssh -R1234:localhost:22 $host < "`dirname $0`/gcc-script.sh" > $tmpdir/$host.output 2>&1
}

for host in "${gccs[@]}" ; do
    run_on_one_node $host & :
done
wait
for host in "${gccs[@]}" ; do
    tail -2 $tmpdir/$host.output | sed -e s,^,$host:,
    echo
done > $tmpdir/digest.txt

# TODO: Include a log of warnings and/or compilation errors.

TODAY=`date +%Y%m%d`
TAR=$tmpdir/nightlytest-$TODAY.tar.gz

(cd $tmpdir ; tar cz --transform s,^,nightlytest-$TODAY/, -f - *.output digest.txt > $TAR)

mutt -s "nightly test $TODAY" -i $tmpdir/digest.txt -a $TAR -- `id -u -n` < /dev/null >/dev/null 2>&1

rm -rf $tmpdir
