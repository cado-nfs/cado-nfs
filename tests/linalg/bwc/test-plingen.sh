#!/usr/bin/env bash

set -e
set -x
# Create a fake sequence

# Note that if we arrive here, we are 64-bit only, since the GFP backends
# for bwc are explicitly disabled on i386 (for now -- most probably
# forever too).

wordsize=64

SHA1BIN=sha1sum
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=sha1 ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then SHA1BIN=shasum ; fi
if ! type -p "$SHA1BIN" > /dev/null ; then
    echo "Could not find a SHA-1 checksumming binary !" >&2
    exit 1
fi

bindir="$1"
shift

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

    : ${TMPDIR:=/tmp}
    TMPDIR=`mktemp -d $TMPDIR/plingen-test.XXXXXXXXXX`

    m="$1"; shift
    n="$1"; shift
    kmax="$1"; shift
    p="$1"; shift
    seed="$1"; shift

    unset thr
    unset ascii

    args=()
    mt_args=()
    for x in "$@" ; do
        case "$x" in
            plingen_program=*) eval "$x";;
            lingen_mpi_threshold*) mt_args=("${mt_args[@]}" "$x");;
            thr*) mt_args=("${mt_args[@]}" "$x"); thr="${x#thr=}";;
            *) args=("${args[@]}" "$x");
                if [[ "$x" =~ ascii ]] ; then
                    if [ "$p" -gt 1048576 ] ; then
                        echo "This test code support ascii tests only for small p" >&2
                        exit 1
                    fi
                    ascii=1
                fi
                ;;
        esac
    done

    nbits_prime=$(sizeinbase2 $p)
    nwords=$((1+nbits_prime/wordsize))

    : ${plingen_program:=plingen_p_$nwords}

    F="$TMPDIR/base"
    if [ "$ascii" ] ; then
        # The perl code below generates ascii test cases which are good
        # provided that p is small. Otherwise, the smallish coefficients
        # we generate are inappropriate and lead to failure, since
        # plingen guesses the length of the ascii input file.
        read -s -r -d '' code <<-'EOF'
            my ($m, $n, $kmax, $p, $seed) = @ARGV;
            my $u = int($seed / 1000);
            my $v = $seed % 1000;
            sub newx {
                    $u = ($u * 2) % 1048573;
                    $v = ($v * 3) % 1048573;
                    return ($u + $v) % $p;
                }

            for my $kmax (1..$kmax) {
                for my $i (1..$m) {
                    print join(" ", map { newx; } (1..$n)), "\n";
                }
                print "\n";
            }
            
EOF
        perl -e "$code" $m $n $((kmax/3)) $p $seed > $F
    else
        # generate $F with exactly ($kmax/3)*$m*$n*$nwords_per_gfp_elt
        # machine words of random data.
        read -s -r -d '' code <<-'EOF'
            my ($m, $n, $kmax, $nwords, $seed) = @ARGV;
            my $u = int($seed / 1000);
            my $v = $seed % 1000;
            sub newx {
                    $u = ($u * 2) % 1048573;
                    $v = ($v * 3) % 1048573;
                    return $u + $v;
                }
            for my $kmax (1..$kmax) {
                for my $i (1..$m) {
                    for my $j (1..$n) {
                        for my $s (1..$nwords) {
                            print pack("Q", newx);
                        }
                    }
                }
            }
            
EOF
        perl -e "$code" $m $n $((kmax/3)) $nwords $seed > $F
    fi

    G="$TMPDIR/seq.txt"
    cat $F $F $F > $G
    rm -f $F

    $bindir/linalg/bwc/$plingen_program m=$m n=$n prime=$p --afile $G "${args[@]}"
    [ -f "$G.gen" ]
    SHA1=$($SHA1BIN < $G.gen)
    SHA1="${SHA1%% *}"

    if [ "$REFERENCE_SHA1" ] ; then
        if [ "${SHA1}" != "${REFERENCE_SHA1}" ] ; then
            echo "$0: Got SHA1 of ${SHA1} but expected ${REFERENCE_SHA1}${REFMSG}. Files remain in ${TMPDIR}" >&2
            exit 1
        fi
    else
        echo "========= $SHA1 ========"
    fi
    rm -f $G.gen

    mpi_bindir=$(perl -ne '/HAVE_MPI\s*"(.*)"\s*$/ && print "$1\n";' $bindir/cado_mpi_config.h)

    if [ "$mpi_bindir" ] && [ "${mt_args[*]}" ] && [ "$thr" ] ; then
        set `echo $thr | tr 'x' ' '`
        if ! [ "$1" ] || ! [ "$2" ] ; then
            echo "Bad test configuration, thr should be of the form \d+x\d+ for MPI test" >&2
            exit 1
        fi
        njobs=$(($1*$2))
        $mpi_bindir/mpiexec -n $njobs $bindir/linalg/bwc/$plingen_program m=$m n=$n prime=$p --afile $G "${args[@]}" "${mt_args[@]}"

        [ -f "$G.gen" ]
        SHA1=$($SHA1BIN < $G.gen)
        SHA1="${SHA1%% *}"

        if [ "$REFERENCE_SHA1" ] ; then
            if [ "${SHA1}" != "${REFERENCE_SHA1}" ] ; then
                echo "$0: Got SHA1 of ${SHA1} but expected ${REFERENCE_SHA1}${REFMSG}. Files remain in ${TMPDIR}" >&2
                exit 1
            fi
            if ! [ "$CADO_DEBUG" ] ; then
                rm -f $G.gen
            fi
        else
            echo "========= $SHA1 ========"
        fi
    fi

    if ! [ "$CADO_DEBUG" ] ; then
        rm -rf "$TMPDIR"
    fi
}

dotest "$@"

