#!/usr/bin/env bash

: ${wdir:?missing}
: ${CADO_NFS_BINARY_DIR:?missing}
: ${CADO_NFS_SOURCE_DIR:?missing}

set -mex

db_backend=sqlite3

while [ $# -gt 0 ] ; do
    if [ "$1" = "--db-backend" ] ; then
        shift
        db_backend="$1"
    else
        echo "unexpected: $1" >&2
        exit 1
    fi
    shift
done

# I love sed -i, but it's not there on os x
sed -e '/tasks.sieve.rels_wanted/ d'                            \
    -e '/tasks.sieve.qrange/ d'                                 \
    < ${CADO_NFS_SOURCE_DIR}/parameters/factor/parameters.F7    \
    > $wdir/parameters.F7

server="${CADO_NFS_BINARY_DIR}/cado-nfs.py"
server_args=(340282366920938463463374607431768211457
    $wdir/parameters.F7
    tasks.polyselect.import=${CADO_NFS_SOURCE_DIR}/tests/misc/F7.poly
    tasks.maxwu=200
    tasks.filter.run=False
    server.whitelist=localhost
    server.threaded=True
    # the default parameter set has 20000, but in fact with qrange=10 and
    # tasks.filter.run=false, we're happy with less.
    tasks.sieve.rels_wanted=15000
)


db_directory=$wdir
if [ $db_backend = mysql ] ; then
    db_directory=/var/lib/mysql
    mysql -e "DROP DATABASE cado_nfs;" 2>/dev/null || :
    # This requires cooperation from the early setup. For example, on CI
    # runners, this setup gets done in the shell function
    # after_package_install in ci/ci.sh
    server_args+=(database=db:mysql://$(id -nu)@_unix_socket//run/mysqld/mysqld.sock/cado_nfs)
fi


if df $db_directory | grep -q tmpfs ; then
    server_args+=(tasks.sieve.qrange=10)
else
    echo "$db_directory is not on a tmpfs filesystem: we reduce the number of server roundtrips because committing to the physical block device is likely to be a huge source of contention"
    server_args+=(tasks.sieve.qrange=100)
fi

server_args+=(
    --wdir "${wdir}/cado-nfs"
    --filelog transaction
    --screenlog warning
    --server
)

logfile="${wdir}/cado-nfs/F7.log"

"$server" "${server_args[@]}" > "$wdir/server.log" 2>&1 &

server_pid=$!
url=
i=0

# As odd as it may seem, we do come across situations where merely
# starting the server takes a very long time
while ! [ "$url" ] && [ $i -lt 45 ] ; do
    if ! [ -f "$logfile" ] ; then
        echo "Waiting for server to create $logfile" >&2
    elif [[ $(grep 'additional.*client.*server' "$logfile") =~ --server=([^ ]*) ]] ; then
        url="${BASH_REMATCH[1]}"
        echo "Server started at $url"
        break
    fi
    let i+=1
    sleep 1
done

if ! [ "$url" ] ; then
    set +e
    echo "server did not start correctly" >&2
    cat "$wdir/server.log" >&2
    cat "$logfile" >&2
    exit 1
fi

client_pids=()
for i in {1..16} ; do
    # drop --bindir intentionally, after all it is another way to exert
    # pressure on the server.
    client_args=(
        --nocncheck
        --server "$url"
        --exit-on-server-gone
        --downloadretry 1
        --basepath "$wdir"
        # --daemon
    )
    client="${CADO_NFS_BINARY_DIR}/cado-nfs-client.py"
    "${client}" "${client_args[@]}" --clientid "client+$i" > "$wdir/client+$i.log" 2>&1 &
    client_pids+=($!)
done

tail -f "$wdir/server.log" < /dev/null &
tail_pid=$!

echo "server pid $server_pid"
echo "tail -f pid $tail_pid"
echo "client pids ${client_pids[@]}"

echo "Returning to server control"
set +e
fg %1
rc=$?
echo "$server finished with rc=$rc"
grep "Finishing early" "$logfile" || :

kill $tail_pid

# Give clients some leeway to finish by themselves
sleep 2
echo "Now killing remaining clients. Errors are not unexpected here"
kill -9 "${client_pids[@]}" 2>/dev/null || :

if [ $rc != 0 ] ; then
    echo "Errors encountered!"
    echo "Errors encountered!" >&2
    exit $rc
elif grep -q Error "$wdir/server.log" ; then
    echo "Errors encountered!"
    echo "Errors encountered!" >&2
    exit 1
else
    echo "All good!"
fi
