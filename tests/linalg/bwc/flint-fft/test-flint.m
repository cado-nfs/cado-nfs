n:=2^fti_depth;
Z:=Integers();
R:=Integers(2^(n*fti_w) + 1);
RP<T>:=PolynomialRing(R);
/* sqrt(2) in R */
rho:=fti_w mod 2 eq 0 select
        R!2^(fti_w div 2)
    else
        R!(2^(3*u)-2^u)^fti_w
    where u is (n*fti_w div 4);
assert rho^(4*n) eq 1;
assert rho^2 eq 2^fti_w;

/* These make sense only in the context of the matrix algorithm (MFA) */
depth1 := fti_depth div 2;
depth2 := fti_depth + 1 - depth1;
n1 := 2^depth1; // for MFA

/* compute the truncation point */
tr:=fti_trunc0;
if tr le 2*n then tr:=2*n+1; end if;
if fti_alg eq 0 then
// trunc must be even and greater than 2n
tr:=fti_trunc0 + (fti_trunc0 mod 2);
else
// trunc must be greater than 2n and multiple of 2*n1
tr:= 2 * n1 * Ceiling(tr / (2 * n1));
end if;

/* transform coefficients are created in bitrev order, at least for the
 * non-MFA algorithm. This script is presently incapable of checking the
 * validity of transformed data for the MFA, but fixing it should not be
 * terribly hard. */
bitrev:=func<x,n|Seqint(Reverse(Intseq(x,2,n)),2)>;
bitrevseq:=func<n|[bitrev(i,n):i in [0..2^n-1]]>;

// mfaorder:=func<x|Seqint([b[1+i]:i in mfabitorder],2) where b is Intseq(x,2,fti_depth+2)>;


/* precompute powers */
seqmatch:=func<S,T,n|S[1..n] eq T[1..n]>;
pows:=[R|1];
for i in [1..2^(fti_depth+2)-1] do Append(~pows, pows[#pows]*rho); end for;


load "P0_before_dft.m"; Q0:=Polynomial([R!Seqint(x,2^64):x in data]);
load "P1_before_dft.m"; Q1:=Polynomial([R!Seqint(x,2^64):x in data]);
load "P0_after_dft.m";  tQ0:=([R!Seqint(x,2^64):x in data]);
load "P1_after_dft.m";  tQ1:=([R!Seqint(x,2^64):x in data]);
load "P2_before_ift.m"; tQ2:=([R!Seqint(x,2^64):x in data]);
load "P2_after_ift.m";  cQ2:=Polynomial([R!Seqint(x,2^64):x in data]);
cQ2 mod:= T^tr;

if fti_alg eq 0 then
ch:=func<Q,tQ|seqmatch([Evaluate(Q,pows[i+1]):i in bitrevseq(fti_depth+2)], tQ, tr)>;
else
ch:=func<Q,tQ|"skipped (MFA)">;
/* I don't understand the twiddles in the MFA code.
 *
function ch(Q, tQ)
    // weird...
    // this is only checking half of the matrix, the "top part" if we're
    // making it a (somewhat weird) matrix of size e.g. 128*32. I think that
    // there's more to it...
    assert #tQ mod 2^depth1 eq 0;
    nrows:=#tQ div 2^depth1;
    tQmat:=Matrix(nrows,2^depth1,tQ);
    for j in [0..2^depth1-1] do
        eval_code:=Transpose(tQmat)[1+j];
        eval_magma:=Vector([Evaluate(Q,pows[1+2*i+128*bitrev(j,depth1)]):i in [0..2^depth2-1]]);
        if eval_magma ne Vector(Eltseq(eval_code)[1..2^depth2]) then
            return false;
        end if;
        for i in [0..2^depth2-1] do
            if tQmat[1+i][1+j] ne Evaluate(Q,pows[1+2*i+128*bitrev(j,depth1)]) then
                return false;
            end if;
        end for;
        for i in [0..3] do
            if tQmat[1+i+2^depth2][1+j] ne Evaluate(Q,pows[1+1+2*i+128*bitrev(j,depth1)]) then
                print i, false;
            end if;
        end for;

tQmat[1+2^depth2][1+j] eq Evaluate(Q,pows[1+1+128*bitrev(j,depth1)]);
tQmat[1+2^depth2+2][1+j] eq Evaluate(Q,pows[1+1+4+128*bitrev(j,depth1)]);


// bits of integers from 0 to 4095 (most significant bit first), and the index
// in tQ, starting from 0, of the evaluation ot primitive root to that power.

 abcde fghijk l    l fghijk edcba   
 01010 000000 0    0 000000 01010   
 01110 000000 0    0 000000 01110   
 01110 000000 1    1 000000 01110   
 01110 000001 0    0 000001 01110   
 01110 000010 0    0 000010 01110   
 01110 000100 0    0 000100 01110   
 01110 001000 0    0 001000 01110   
 01110 010000 0    0 010000 01110   
 01110 100000 0    0 100000 01110   
 01101 110000 0    0 110000 10110   

(conversely, A BCDEFG HIJKL goes to LKIJH BCDEFG A)

function matrixposition_to_evalindex(ABCDEFG,HIJKL)
    A,BCDEFG:=Quotrem(ABCDEFG, 2^depth2);
    LKIJH:=bitrev(HIJKL, depth1);
    return LKIJH*2^(depth2+1)+2*BCDEFG+A;
end function;
function matrixindex_to_evalindex(ABCDEFGHIJKL)
    ABCDEFG,HIJKL := Quotrem(ABCDEFGHIJKL, 2^depth1);
    return matrixposition_to_evalindex(ABCDEFG,HIJKL);
end function;
function evalindex_to_matrixposition(abcdefghijkl)
    abcde,fghijkl:=Quotrem(abcdefghijkl,2^(depth2+1));
    edcba:=bitrev(abcde,depth1);
    fghijk,l:=Quotrem(fghijkl,2);
    return l*2^depth2+fghijk, edcba;
end function;
function evalindex_to_matrixindex(abcdefghijkl)
    lfghijk,edcba:=evalindex_to_matrixposition(abcdefghijkl);
    return lfghijk*2^depth1+edcba;
end function;






            eval_magma:=Vector([Evaluate(Q,pows[1+1+2*i+128*bitrev(j,depth1)]):i in [0..nrows-2^depth2-1]]);
            if eval_magma ne Vector(Eltseq(eval_code)[2^depth2+1..nrows]) then
                return false;
            end if;
        end for;
    end for;
    return true;
    // return &and[Evaluate(Q,pows[1+2*i+128*bitrev(j,depth1)]) eq tQ[1+j+2^depth1*i]:i in [0..2^depth2-1],j in [0..2^depth1-1]];
end function;
*/


end if;

/* work around magma bugs */
if not assigned P0 then P0:=0; end if;
if not assigned P1 then P1:=0; end if;
if not assigned P2 then P2:=0; end if;
if not assigned A0 then A0:=0; end if;
if not assigned A1 then A1:=0; end if;
if not assigned A2 then A2:=0; end if;

print "Active test: ", check;

if check eq "test_mul0" then
    assert Q0 eq Polynomial(R,Intseq(A0,2^fti_bits));
    assert Q1 eq Polynomial(R,Intseq(A1,2^fti_bits));
    Q2:=Q0*Q1 mod (T^(4*n)-1);
    assert Q2 eq cQ2;
    assert A2 eq Evaluate(ChangeRing(Q2,Z),2^fti_bits);
elif check eq "test_mul" then
    assert Q0 eq Polynomial(R,Intseq(A0,2^fti_bits));
    assert Q1 eq Polynomial(R,Intseq(A1,2^fti_bits));
    Q2:=(2*Q0*Q1+Q0+Q1) mod (T^(4*n)-1);
    assert Q2 eq cQ2;
    assert A2 eq Evaluate(ChangeRing(Q2,Z),2^fti_bits);
elif check eq "test_mulmod" then
    assert Q0 eq Polynomial(R,Intseq(A0,2^fti_bits));
    assert Q1 eq Polynomial(R,Intseq(A1,2^fti_bits));
    Q2:=(2*Q0*Q1+Q0+Q1) mod (T^(4*n)-1);
    assert Q2 eq cQ2;
    assert A2 eq Evaluate(ChangeRing(Q2,Z),2^fti_bits) mod (2^(4*n*fti_bits)-1);
elif check eq "test_mul_fppol" then
    zP0:=Evaluate(ChangeRing(P0,Z),2^fti_ks_coeff_bits);
    zP1:=Evaluate(ChangeRing(P1,Z),2^fti_ks_coeff_bits);
    assert Q0 eq Polynomial(R,Intseq(zP0,2^fti_bits));
    assert Q1 eq Polynomial(R,Intseq(zP1,2^fti_bits));
    Q2:=(2*Q0*Q1+Q0^2+Q0+Q1) mod (T^(4*n)-1);
    assert Q2 eq cQ2;
    assert P2 eq 2*P0*P1+P0^2+P0+P1;
elif check eq "test_mp_fppol" then
    zP0:=Evaluate(ChangeRing(P0,Z),2^fti_ks_coeff_bits);
    zP1:=Evaluate(ChangeRing(P1,Z),2^fti_ks_coeff_bits);
    assert Q0 eq Polynomial(R,Intseq(zP0,2^fti_bits));
    assert Q1 eq Polynomial(R,Intseq(zP1,2^fti_bits));
    // assert fti_alg eq 0;
    // if fti_alg eq 0 then
        // it's not satisfactory, and should be investigated. Apparently
        // things are a bit different with MFA. Or maybe it's the wraparound
        // thing.
        Q2:=(Q0*Q1) mod (T^(4*n)-1);
        assert Q2 eq cQ2;
    // end if;
    zP2:=Evaluate(ChangeRing(Q2,Z),2^fti_bits);
    P2:=Intseq(zP2,2^fti_ks_coeff_bits)[Degree(P0)+1..Degree(P1)+1];
    MP:=func<P0,P1|[Coefficient(P,i) :i in [Min(Degree(P0),Degree(P1))..Max(Degree(P0),Degree(P1))]] where P is P0*P1>;
    assert MP(P0,P1) eq P2;
else
    print "check not understood";
end if;

print "checking transform T0: ", ch(Q0, tQ0);
print "checking transform T1: ", ch(Q1, tQ1);
print "checking transform T2: ", ch(Q2, tQ2);


