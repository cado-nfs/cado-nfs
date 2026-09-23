#!/usr/bin/env bash

# Performance harness for the cofactorization code, based on testbench.
#
# Usage: bench_facul.sh [options] <testbench A> [<testbench B>]
#
# Each case is a testbench call. All physical cores of the machine run a
# copy of the case at the same time, so that the timings are those of a
# fully busy machine. When two binaries are given, A runs on even cores
# and B on odd cores, so that both see the same conditions, and the
# output of both is compared.
#
# Every process runs its own copy of the binary. Where the text of a
# binary file lands in memory changes its speed by up to 2%: two copies
# of the same file consistently differ by that much. With one copy per
# process, the median averages over as many placements as there are
# processes.
#
# Options:
#   -w <dir>       work directory for the input files (default: /tmp/bench_facul)
#   -t <seconds>   target duration of each case (default: 2)
#   -n <ncores>    number of physical cores to use (default: all)
#   -smt           also load the SMT siblings (two copies per core)
#   -only <regex>  only run the cases whose name matches <regex>
#   -list          list the cases and exit
#
# The output has one line per case: the median time per call in
# microseconds for A (and B, and the ratio B/A), and the interquartile
# range relative to the median.

set -e

wdir=/tmp/bench_facul
target=2
ncores=
smt=
only=.
list=
while [ $# -gt 0 ] ; do
    case "$1" in
        -w) wdir="$2"; shift 2;;
        -t) target="$2"; shift 2;;
        -n) ncores="$2"; shift 2;;
        -smt) smt=1; shift;;
        -only) only="$2"; shift 2;;
        -list) list=1; shift;;
        -*) echo "Unknown option $1" >&2; exit 1;;
        *) break;;
    esac
done

# The cofactors below are primes just above 2^41, 2^73, 2^103, 2^128.
# Multiplied by the primes just above 2^21 that are in primes.txt, they
# give numbers that are handled by the modredc_ul, modredc_15ul,
# modredc_2ul2 and mpz layers respectively, close to the maximal size of
# each layer.
declare -A cof=(
    [ul]=2199023255579
    [15ul]=9444732965739290427421
    [2ul2]=10141204801825835211973625643089
    [mpz]=340282366920938463463374607431768211507
)

# The default strategy has P-1 with B1=315, P+1 with B1=525, and ECM
# curves with B1 from 105 upwards.
methods=(
    "pm1:-pm1 315 2205"
    "pp1_27:-pp1_27 525 3255"
    "pp1_65:-pp1_65 525 3255"
    "ecm_b12:-ecm 315 5355 11"
    "ecm_m12:-ecmm12 105 3255 2"
    "ecm_m16:-ecmm16 105 3255 1"
    "ecm_em12:-ecmem12 315 5355 1"
    "ecm_m12_big:-ecmm12 1000 50000 5"
)

cases=()
for w in ul 15ul 2ul2 mpz ; do
    for m in "${methods[@]}" ; do
        cases+=("${m%%:*}.$w:primes.txt:-cof ${cof[$w]} ${m#*:}")
    done
done
# The full default strategy, on numbers without prime factors below the
# factor base bound, as las gives them to facul. Most of them are not
# smooth, and go through the whole chain.
for b in 60 90 120 ; do
    cases+=("strat.$b:surv$b.txt:-strat -fbb 1048576 -lpb 31")
done

if [ "$list" ] ; then
    for c in "${cases[@]}" ; do
        name="${c%%:*}" ; rest="${c#*:}"
        if ! [[ $name =~ $only ]] ; then continue ; fi
        echo "$name: testbench -inp ${rest%%:*} ${rest#*:}"
    done
    exit 0
fi

A="$1"
B="$2"
if ! [ -x "$A" ] || { [ "$B" ] && ! [ -x "$B" ] ; } ; then
    echo "Usage: $0 [options] <testbench A> [<testbench B>]" >&2
    exit 1
fi

mkdir -p "$wdir"

if ! [ -f "$wdir/surv120.txt" ] ; then
    python3 - "$wdir" <<'EOF'
import sys, random, math
wdir = sys.argv[1]
random.seed(1)
N = 1 << 21
# primes in [2^21, 2^21 + 4*10^6]
L = 4 * 10**6
small = [p for p in range(2, int(math.isqrt(N + L)) + 1)
         if all(p % q for q in range(2, int(math.isqrt(p)) + 1))]
seg = bytearray([1]) * L
for p in small:
    for k in range((-N) % p, L, p):
        seg[k] = 0
with open(f"{wdir}/primes.txt", "w") as f:
    for i in range(L):
        if seg[i]:
            print(N + i, file=f)
# numbers of b bits without prime factors below fbb: sieve 25 random
# intervals of length 100000 by the primes below fbb.
fbb = 1 << 20
isp = bytearray([1]) * fbb
isp[0] = isp[1] = 0
for i in range(2, math.isqrt(fbb) + 1):
    if isp[i]:
        isp[i*i::i] = bytearray(len(isp[i*i::i]))
fbprimes = [i for i in range(fbb) if isp[i]]
for b in (60, 90, 120):
    with open(f"{wdir}/surv{b}.txt", "w") as f:
        for j in range(25):
            X = random.getrandbits(b - 1) | (1 << (b - 1))
            L = 100000
            seg = bytearray([1]) * L
            for p in fbprimes:
                seg[(-X) % p::p] = bytearray(len(seg[(-X) % p::p]))
            for i in range(L):
                if seg[i]:
                    print(X + i, file=f)
EOF
fi

# physical cores: the first logical cpu of each core.
mapfile -t cpus < <(lscpu -p=CPU,CORE | grep -v '^#' | awk -F, '!seen[$2]++ {print $1}')
mapfile -t siblings < <(lscpu -p=CPU,CORE | grep -v '^#' | awk -F, 'seen[$2]++ {print $1}')
: ${ncores:=${#cpus[@]}}

tmp=$(mktemp -d "$wdir/run.XXXXXX")
trap 'rm -rf "$tmp"' EXIT

mkdir "$tmp/bin"
for k in $(seq 0 $((ncores - 1))) ; do
    bin="$A"
    if [ "$B" ] && [ $((k % 2)) = 1 ] ; then bin="$B" ; fi
    cp "$bin" "$tmp/bin/$k"
    if [ "$smt" ] ; then cp "$bin" "$tmp/bin/$k.smt" ; fi
done

stats() {
    sort -g | awk '{a[NR]=$1} END {
        q1=a[int(NR/4)+1]; med=a[int(NR/2)+1]; q3=a[int(3*NR/4)+1];
        printf "%.4f %.2f\n", med, 100*(q3-q1)/med }'
}

# run <binary> <args...> ; prints the time per call in microseconds
percall() {
    "$@" | awk '/^Total time/ { print $7 }'
}

if [ "$B" ] ; then
    printf "%-18s %10s %10s %8s %7s %7s %s\n" case "A(us)" "B(us)" "B/A" "iqrA%" "iqrB%" output
else
    printf "%-18s %10s %7s\n" case "A(us)" "iqr%"
fi

for c in "${cases[@]}" ; do
    name="${c%%:*}" ; rest="${c#*:}"
    input="$wdir/${rest%%:*}" ; args="${rest#*:}"
    if ! [[ $name =~ $only ]] ; then continue ; fi

    # calibrate so that each copy runs for about $target seconds. The
    # calibration run must itself last long enough to be measured.
    n=200
    while true ; do
        t=$(percall taskset -c ${cpus[0]} "$A" -inp "$input" -inpstop $n $args)
        if [ "$t" ] && awk -v t="$t" -v n=$n 'BEGIN { exit !(t * n >= 50000) }' ; then
            break
        fi
        if [ $n -gt 10000000 ] ; then
            echo "$name: cannot calibrate" >&2 ; exit 1
        fi
        n=$((n * 4))
    done
    count=$(awk -v t="$t" -v T="$target" 'BEGIN { n = int(T * 1e6 / t); print n < 200 ? 200 : n }')

    rm -f "$tmp"/a.* "$tmp"/b.*
    for k in $(seq 0 $((ncores - 1))) ; do
        out="$tmp/a.$k"
        if [ "$B" ] && [ $((k % 2)) = 1 ] ; then out="$tmp/b.$k" ; fi
        percall taskset -c ${cpus[$k]} "$tmp/bin/$k" -inp "$input" -inpstop $count $args > "$out" &
        if [ "$smt" ] ; then
            percall taskset -c ${siblings[$k]} "$tmp/bin/$k.smt" -inp "$input" -inpstop $count $args > "$out.smt" &
        fi
    done
    wait
    read medA iqrA < <(cat "$tmp"/a.* | stats)

    if [ "$B" ] ; then
        read medB iqrB < <(cat "$tmp"/b.* | stats)
        # the factors that A and B find must be the same.
        dA=$("$A" -inp "$input" -inpstop 2000 -q -vf -vnf $args | md5sum)
        dB=$("$B" -inp "$input" -inpstop 2000 -q -vf -vnf $args | md5sum)
        same=same ; [ "$dA" = "$dB" ] || same=DIFFERENT
        ratio=$(awk -v a="$medA" -v b="$medB" 'BEGIN { printf "%.4f", b/a }')
        printf "%-18s %10s %10s %8s %7s %7s %s\n" "$name" "$medA" "$medB" "$ratio" "$iqrA" "$iqrB" "$same"
    else
        printf "%-18s %10s %7s\n" "$name" "$medA" "$iqrA"
    fi
done
