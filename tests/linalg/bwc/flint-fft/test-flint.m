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
assert rho^(4*n) eq 1;

load "P0_before_dft.m"; Q0:=Polynomial([R!Seqint(x,2^64):x in data]);
load "P1_before_dft.m"; Q1:=Polynomial([R!Seqint(x,2^64):x in data]);
load "P0_after_dft.m";  tQ0:=([R!Seqint(x,2^64):x in data]);
load "P1_after_dft.m";  tQ1:=([R!Seqint(x,2^64):x in data]);
load "P2_before_ift.m"; tQ2:=([R!Seqint(x,2^64):x in data]);
load "P2_after_ift.m";  cQ2:=Polynomial([R!Seqint(x,2^64):x in data]);
cQ2 mod:= T^tr;

ch:=func<Q,tQ|fti_alg eq 0 select  seqmatch([Evaluate(Q,pows[i+1]):i in bitrevseq(fti_depth+2)], tQ, tr) else "skipped (MFA)">;

/* work around magma bugs */
if not assigned P0 then P0:=0; end if;
if not assigned P1 then P1:=0; end if;
if not assigned P2 then P2:=0; end if;
if not assigned A0 then A0:=0; end if;
if not assigned A1 then A1:=0; end if;
if not assigned A2 then A2:=0; end if;

print "Active test: ", check;

if check eq "test_mul" then
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
    if fti_alg eq 0 then
        // it's not satisfactory, and should be investigated. Apparently
        // things are a bit different with MFA. Or maybe it's the wraparound
        // thing.
        Q2:=(Q0*Q1) mod (T^(4*n)-1);
        assert Q2 eq cQ2;
    end if;
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


