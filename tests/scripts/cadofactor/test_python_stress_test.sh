#!/usr/bin/env bash

: ${wdir:?missing}
: ${CADO_NFS_BINARY_DIR:?missing}
: ${CADO_NFS_SOURCE_DIR:?missing}

set -m
set -e

server="${CADO_NFS_BINARY_DIR}/cado-nfs.py"
server_args=(340282366920938463463374607431768211457
    ${CADO_NFS_SOURCE_DIR}/parameters/factor/parameters.F7
    tasks.polyselect.import=${CADO_NFS_SOURCE_DIR}/tests/misc/F7.poly
    tasks.maxwu=200
    tasks.sieve.qrange=10
    tasks.filter.run=False
    server.whitelist=localhost
    server.threaded=True
    # the default parameter set has 20000, but in fact with qrange=10 and
    # tasks.filter.run=false, we're happy with less.
    tasks.sieve.rels_wanted=15000
    --wdir "${wdir}/cado-nfs"
    --server
)

logfile="${wdir}/cado-nfs/F7.log"

"$server" "${server_args[@]}" > "$wdir/server.log" 2>&1 &

server_pid=$!
url=
i=0

while ! [ "$url" ] && [ $i -lt 4 ] ; do
    if [[ $(grep 'additional.*client.*server' "$logfile") =~ --server=([^ ]*) ]] ; then
        url="${BASH_REMATCH[1]}"
        echo "Server started at $url"
        break
    fi
    let i+=1
    sleep 1
done

if ! [ "$url" ] ; then
    echo "server did not start correctly" >&2
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
        # --daemon
    )
    client="${CADO_NFS_BINARY_DIR}/cado-nfs-client.py"
    "${client}" "${client_args[@]}" --clientid "client+$i" > "$wdir/client+$i.log" 2>&1 &
    client_pids+=($!)
done

echo "server pid $server_pid"
echo "client pids ${client_pids[@]}"

echo "Returning to server control"
set +e
fg %1
rc=$?
echo "$server finished with rc=$rc"
grep "Finishing early" "$wdir/server.log" || :

# Give clients some leeway to finish by themselves
sleep 2
echo "Now killing remaining clients. Errors are not unexpected here"
kill -9 "${client_pids[@]}"

if [ $rc = 0 ] ; then
    echo "All good!"
else
    echo "Errors encountered!"
    echo "Errors encountered!" >&2
    exit $rc
fi
