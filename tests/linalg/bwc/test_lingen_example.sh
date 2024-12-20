#!/usr/bin/env bash

set -ex

: ${wdir?missing}
: ${CADO_NFS_BINARY_DIR?missing}

: ${tuning_thresholds=recursive:10,flint:10,notiming:0}

files=()

while [ $# -gt 0 ] ; do
    if [[ $1 =~ ^tuning_thresholds=[a-z0-9:,-]*$ ]] ; then
        eval "$1"
    elif [[ $1 =~ ^layer=[a-z0-9]*$ ]] ; then
        eval "$1"
    else
        f="$1"
        # os x bash has inconsistent regexps.
        case "$f" in
            *.E.aux|*.E.single.data|*.E) : ;;
            *)  echo "Unexpected arg: $f" >&2
                exit 1
        esac
        f=${f%%.aux}
        f=${f%%.single.data}
        f=${f%%.E}
        if ! [ -f "$f.E.aux" ] ; then
            echo "Missing companion file: $f.E.aux" >&2
            exit 1
        fi
        if ! [ -f "$f.E.single.data" ] ; then
            echo "Missing companion file: $f.E.aux" >&2
            exit 1
        fi
        files+=("$f")
    fi
    shift
done

: ${layer:?missing}

prg="$CADO_NFS_BINARY_DIR/tests/linalg/bwc/test_lingen_mini_${layer}"
check="$CADO_NFS_BINARY_DIR/linalg/bwc/lingen_verify_checkpoints_${layer}"
args=(tuning_thresholds="$tuning_thresholds")

# for f in bug/cp*/5*.E.aux ; do f0=`basename "$f" .E.aux` ; D=`dirname "$f"` ; g="${f////_}"; g=`basename "$g" .E.aux`; mkdir /tmp/$g > /dev/null 2>&1  || : ; /bin/cp -f "$D/$f0.E.aux" "$D/$f0.E.single.data" "$D/$f0.pi.aux" /tmp/$g ; /bin/cp -f "$D/$f0.pi.single.data" "/tmp/$g/pi_c.data" ;  ./cado/build/coffee/tests/linalg/bwc/test_lingen_mini_p1 tuning_thresholds=recursive:10,flint:10,notiming:0  /tmp/$g/$f0.E.aux > /dev/null 2>&1 ; diff /tmp/$g/pi_c.data /tmp/$g/$f0.E.aux.pi ; for  v in $f0.E.aux.pi pi_c.data ; do ln -sf $v /tmp/$g/$f0.pi.single.data ;  ./cado/build/coffee.debug/linalg/bwc/lingen_verify_checkpoints_p1 m=4 n=4 prime=67 -- /tmp/$g/$f0.{E,pi} && echo "$v ok" || echo "$v NOK" ; done; done

i=0
for f in "${files[@]}" ; do
    FD="$wdir/test$i"
    mkdir "$FD"
    read ff fn m n level t0 t1 t nc p < <(head -n 5 "$f.E.aux" | xargs echo)
    level=0
    F="$FD/$level.$t0.$t1"
    cp "$f.E.aux" "$F.E.aux"
    cp "$f.E.single.data" "$F.E.single.data"
    "$prg" "${args[@]}" checkpoint-directory="$FD" "$F.E.aux"
    # cp "$F.E.aux.pi" "$F.pi"
    "$check" m=$m n=$n prime=$p -- "$F.E" "$F.pi"
    let i+=1
done
