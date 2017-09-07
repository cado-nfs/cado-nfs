#!/usr/bin/env bash

# This shell script may be used to generate the proper crc code with
# magma, to be filled eventually as a magic test name in Makefile.am

# These constants are used in our random generation process.
cstart=0x32104567
cratio=0x76540123


while [ $# -gt 0 ] ; do
    if [ "$1" = "--n1" ] ; then n1=$2; shift; shift;
    elif [ "$1" = "--n2" ] ; then n2=$2; shift; shift;
    elif [ "$1" = "--n" ] ; then n1=$2; n2=$2; shift; shift;
    elif [ "$1" = "--k" ] ; then k=$2; shift; shift;
    elif [ "$1" = "--start" ] ; then cstart=$2; shift; shift;
    elif [ "$1" = "--ratio" ] ; then cratio=$2; shift; shift;
    else
        echo "Usage: $0 [--n1 <n1>] [--n2 <n2>] [--n <n>] [--k <k>] [--start <start>] [--ratio <ratio>]" >&2
        exit 1
    fi
done

read -s -r -d '' code <<-EOF
    start:=$cstart;
    ratio:=$cratio;
    k:=$k;
    n1:=$n1;
    n2:=$n2;
    KP<x>:=PolynomialRing(GF(2));
    p:=x^32 + x^7 + x^6 + x^2 + 1;
    E:=ext<GF(2)|p>;
    v:=Integers(2^32)!start;
    function gen(v,n1)
        fs:=[E|];
        project:=Vector([E|(E!x^32)^i:i in [0..n1-1]]);
        for i in [1..k^2] do
            L:=[E|];
            for j in [1..n1] do
                Append(~L, E!KP!Intseq(Integers()!v,2));
                v *:= ratio;
            end for;
            Append(~fs, (project, Vector(L)));
        end for;
        return v,fs;
    end function;
    v,fs:=gen(v,n1);
    v,gs:=gen(v,n2);

    f:=Matrix(E,k,k,fs);
    g:=Matrix(E,k,k,gs);

    h:=f*g;
    h:=Evaluate(Polynomial(Eltseq(h)),E!x^(32*(n1+n2)));
    check_code:=IntegerToString(Seqint(ChangeUniverse(Eltseq(h),Integers()),2),16);
    print check_code;
    quit;
EOF

echo $code > /tmp/a.m
expect=$(echo $code|magma -b)
expect=$(printf '%08x\n' 0x$expect)

echo check_${n1}_${n2}_${k}_${expect}


