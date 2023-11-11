#!/usr/bin/env bash

# just to be sure
unset DISPLAY

set -e
if [ "$CADO_DEBUG" ] ; then set -x ; fi

if ! [ "$WDIR" ] ; then
    echo "Want \$WDIR" >&2
    exit 1
fi

: ${bindir:=$PROJECT_BINARY_DIR}
: ${bindir?missing variable}

# inject the variables that were provided by guess_mpi_configs
if [ "$mpi" ] ; then
    eval "$exporter_mpirun"
    eval "$exporter_mpi_extra_args"
    set -- "$@" mpi="$mpi"
fi

# Create a fake sequence

wordsize=$(awk '/ULONG_BITS/ { print $3 }' $PROJECT_BINARY_DIR/cado_config.h)

SHA1BIN=sha1sum
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=sha1 ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=shasum ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then
    echo "Could not find a SHA-1 checksumming binary !" >&2
    exit 1
fi

# This is not accurate, unfortunately. There seems to be no way to do
# long integer arithmetic in pure bash.
sizeinbase2() {
    a=$1
    if [ "${#a}" -le 10 ] ; then 
        x=$(printf '%x' $1)
        len=$((4*${#x}))
        case $x in
            0) echo $len; return;;
            1*) echo $((len-3)); return;;
            [23]*) echo $((len-2)); return;;
            [4567]*) echo $((len-1)); return;;
            [89a-fA-F]*) echo $len; return;;
        esac
        echo "sizeinbase2 error $1 -> $x -> uh ?" >&2; exit 1
    else
        # a is M * 10^e, so log2(a) = log2(M) + e*log2(10), and too bad
        # for the inaccuracy we get...
        logM=$(sizeinbase2 "${a:0:10}")
        elog10=$((332*(${#a}-10)/100))
        loga=$((logM+elog10))
        echo $loga
        return
    fi
}


dotest() {
    REFERENCE_SHA1="$1"
    shift

    m="$1"; shift
    n="$1"; shift
    kmax="$1"; shift
    p="$1"; shift
    seed="$1"; shift

    unset mpi
    unset ascii

    args=()
    mpi_specific_args=()
    mpi_extra_args=()
    ONLY_TUNE=
    for x in "$@" ; do
        case "$x" in
            lingen_program=*) eval "$x";;
            skip_single=*) eval "$x";;
            lingen_mpi_threshold*) mpi_specific_args+=("$x");;
            # the mpi_magic= argument is interpreted by guess_mpi_configs.sh
            mpi_magic=*) mpi_magic="${x#mpi_magic=}";;
            # mpi_extra_args[@] is used by guess_mpi_configs.sh to form
            # the mpirun[@] precommand
            mpi_extra_args=*)
                eval mpi_extra_args+=(${x#mpi_extra_args=});;
            *) args+=("$x");
                if [[ "$x" =~ ascii ]] ; then
                    if [ "$p" -gt 1048576 ] ; then
                        echo "This test code support ascii tests only for small p" >&2
                        exit 1
                    fi
                    ascii=1
                elif [[ "$x" =~ --tune ]] ; then
                    ONLY_TUNE=1
                fi
                ;;
        esac
    done

    nbits_prime=$(sizeinbase2 $p)
    nwords=$((1+nbits_prime/wordsize))

    : ${lingen_program:=lingen_p$nwords}

    # The perl code for generating random data is using the seed argument
    # in a weird way. seeds that are less than 1000 actually lead do
    # random matrices that are likely to have some very non-random
    # behaviour. It would be best to change that eventually, e.g. with
    # the commented-out choices below. However that would entail changing
    # all sha1sums in the tests, and I'm lazy.
    F="$WDIR/base"
    if [ "$ascii" ] ; then
        # The perl code below generates ascii test cases which are good
        # provided that p is small. Otherwise, the smallish coefficients
        # we generate are inappropriate and lead to failure, since
        # lingen guesses the length of the ascii input file.
        read -s -r -d '' code <<-'EOF'
            my ($m, $n, $kmax, $p, $seed) = @ARGV;
            my $u = int($seed / 1000);
            my $v = $seed % 1000;
            # my $u = 1 + ($seed & 0x55555555);
            # my $v = 1 + ($seed & 0xaaaaaaaa);
            sub newx {
                    $u = ($u * 2) % 1048573;
                    $v = ($v * 3) % 1048573;
                    return $u + $v;
                    # + int(1000 * $u / 7);
                }
            for my $kmax (1..$kmax) {
                for my $i (1..$m) {
                    print join(" ", map { newx; } (1..$n)), "\n";
                }
                print "\n";
            }
            
EOF
        perl -e "$code" $m $n $((kmax/3)) $p $seed > $F
        args+=(--input-length $((3*(kmax/3))))
    else
        # generate $F with exactly ($kmax/3)*$m*$n*$nwords_per_gfp_elt
        # machine words of random data.
        read -s -r -d '' code <<-'EOF'
            my ($wordsize, $m, $n, $kmax, $nwords, $seed) = @ARGV;
            my $pack = ($wordsize == 64) ? "Q" : "L";
            my $u = int($seed / 1000);
            my $v = $seed % 1000;
            # my $u = 1 + ($seed & 0x55555555);
            # my $v = 1 + ($seed & 0xaaaaaaaa);
            sub newx {
                    $u = ($u * 2) % 1048573;
                    $v = ($v * 3) % 1048573;
                    return $u + $v;
                    # + int(1000 * $u / 7);
                }
            for my $kmax (1..$kmax) {
                for my $i (1..$m) {
                    for my $j (1..$n) {
                        for my $s (1..$nwords) {
                            print pack($pack, newx);
                        }
                    }
                }
            }
            
EOF
        perl -e "$code" $wordsize $m $n $((kmax/3)) $nwords $seed > $F
    fi

    G="$WDIR/seq.txt"
    cat $F $F $F > $G
    rm -f $F

    # For mpi uses of this script, we expect to be called from
    # do_with_mpi.sh. In this case, $mpi and $mpirun[@] are set.
    if [ "$mpi" ] ; then
        args+=("${mpi_specific_args[@]}")
        if [ "$ONLY_TUNE" ] ; then
            # push --tune at the very end of the argument list, otherwise
            # openmpi gobbles it...
            nargs=()
            for x in "${args[@]}" ; do
                if [ "$x" = "--tune" ] ; then : ; else nargs+=("$x") ; fi
            done
            args=("${nargs[@]}" tuning_mpi="$mpi" --tune)
            set -- "${mpirun[@]}"
            mpirun=()
            while [ $# -gt 0 ] ; do
                mpirun+=("$1")
                if [ "$1" = "-n" ] ; then
                    shift
                    mpirun+=(1)
                fi
                shift
            done
        else
            args+=(mpi="$mpi")
        fi
    fi

    run=("${mpirun[@]}" $bindir/linalg/bwc/$lingen_program m=$m n=$n prime=$p --afile $G "${args[@]}")
    echo "${run[@]}"
    "${run[@]}"

    if [ "$ONLY_TUNE" ] ; then exit 0 ; fi

    [ -f "$G.gen" ]

    if [ "$wordsize" = 64 ] || [ $((nwords % 2)) = 0 ] || [ "$ascii" ] ; then
        feed() { cat "$@" ; }
    else
        # We need to expand each n-(32-bit)-word integer to n+1
        feed() {
        read -s -r -d '' code <<-'EOF'
        my ($nwords) = @ARGV;
        my $ewords = $nwords + 1;
        while (sysread(STDIN, my $x, $nwords * 4)) {
            $x = pack("L$ewords", unpack("L$nwords", $x), 0);
            syswrite(STDOUT, $x);
        }

EOF
        cat "$@" | perl -e "$code" $nwords
    }
    fi
    SHA1=$(feed $G.gen | $SHA1BIN)
    SHA1="${SHA1%% *}"

    if [ "$REFERENCE_SHA1" ] ; then
        ok_sha=(${REFERENCE_SHA1//,/ })
        ok=
        for a in "${ok_sha[@]}" ; do
            if [ "${SHA1}" = "$a" ] ; then
                ok=1
                break
            fi
        done
        if ! [ "$ok" ] ; then
            one_of=
            if [ ${#ok_sha[@]} -gt 1 ] ; then one_of="one of " ; fi
            echo "$0: Got SHA1 of ${SHA1} but expected ${one_of}${REFERENCE_SHA1}${REFMSG}. Files remain in ${WDIR}" >&2
            exit 1
        fi
        echo "$SHA1 (as expected)"
    else
        echo "========= $SHA1 ========"
    fi
}

dotest "$@"

