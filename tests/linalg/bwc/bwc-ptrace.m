
/* TODO: report this to Allan:
 *
 * > Vector(GF(11),Matrix(2,2,[GF(11)|1,2,3,4]));
 * ( 1  2)
 * ( 3  4)
 *
 */

/* This magma script must be run from the directory containing the data files
 * from a bwc run (converted to magma format by bwc-ptrace.sh */


load "mn.m";
load "t.m";
M:=var;
M_orig:=M;

/* We'll pad the matrix to something which has nice dimensions, multiples of
 * nh*nv. Otherwise I'll get crazy very soon */
load "placemats.m";
/* Now M is M_orig zero-padded in rows & cols */
/* while Mt is the matrix as dispatched on the different threads/jobs */

KP<X>:=PolynomialRing(GF(p));

if p eq 2 then
    g:=func<var|Transpose(Matrix([Vector(GF(2),Intseq(x,2,64)):x in var]))>;
    ung:=func<M|[Seqint(Eltseq(x),2):x in Rows(Transpose(ChangeRing(M,Integers())))]>;
else
    function group(nlimbs)
        return func<var|Vector([GF(p)|Seqint([x mod 2^64:x in var[i*k+1..i*k+k]],2^64):i in [0..#var div k - 1]] where k is nlimbs)>;
    end function;
    plimbs:=Ceiling(Log(2^64,p));
    g:=group(plimbs);
    ung:=func<S|&cat[Intseq(x,2^64,plimbs):x in S]>;
end if;

splitwidth:=p eq 2 select 64 else 1;

// assert nc le nr;
// nc:=nr;
// x:=Matrix(GF(p),nr,nc,[]);InsertBlock(~x,M,1,1);M:=x;

load "rw.m"; rw:=var;
load "cw.m"; cw:=var;

load "b.m"; 
nr_orig:=nr;
nc_orig:=nc;
nr:=tr;
nc:=tc;
load "vectorspace.m";

/*
//L:=[]; for i in [0..nh-1] do
//t:=nr div nh + (i lt nr mod nh select 1 else 0); 
//s:=i*nrp;
//L cat:=[s+1..s+t];
//end for;Lr:=L;
//for j in [1..#Lr] do Cr[Lr[j]][j]:=1; end for;
Cr:=Matrix(GF(p),nh*nrp,nr,[]);
InsertBlock(~Cr,IdentityMatrix(GF(p), Minimum(nh*nrp,nr)),1,1);
//L:=[]; for j in [0..nv-1] do
//t:=nc  div nv + (j lt nc mod nv select 1 else 0);
//s:=j*ncp;
//L cat:=[s+1..s+t];
//end for;Lc:=L;
// for j in [1..#Lc] do Cc[Lc[j]][j]:=1; end for;
Cc:=Matrix(GF(p),nv*ncp,nc,[]);
InsertBlock(~Cc,IdentityMatrix(GF(p),Minimum(nv*ncp,nc)),1,1);

function is_almost_one(M, n)
    if not IsDiagonal(M) then return false; end if;
    for i in [1..n] do if M[i,i] ne 1 then return false; end if; end for;
    return true;
end function;

assert is_almost_one(Transpose(Cr)*Cr, Minimum(nh*nrp,nr));
assert is_almost_one(Transpose(Cc)*Cc, Minimum(nv*ncp,nc));
*/

ncx := nv*ncp;
nrx := nh*nrp;

nz:=nrp div nv;
assert nz eq ncp div nh;
pr:=func<x|(qr*nh+qq)*nz+r+1 where qq,qr is Quotrem(q, nv) where q,r is Quotrem(x-1,nz)>;
prinv:=func<x|(qr*nv+qq)*nz+r+1 where qq,qr is Quotrem(q, nh) where q,r is Quotrem(x-1,nz)>;
prp:=SymmetricGroup(nrx)![pr(x):x in [1..nrx]];
Pr:=PermutationMatrix(GF(p),[pr(x):x in [1..nrx]]);
assert Pr eq PermutationMatrix(GF(p),prp);



sc:=SymmetricGroup(nv*ncp)!colperm;
if assigned rowperm then
    sr:=SymmetricGroup(nh*nrp)!rowperm;
else
    sr:=SymmetricGroup(nh*nrp)!colperm;
end if;
sr_inv:=sr^-1;
sc_inv:=sc^-1;
Sr:=PermutationMatrix(GF(p),sr);
Sc:=PermutationMatrix(GF(p),sc);
if assigned rowperm then
    Sr:=Pr^-1*Sr;
    sr:=prp^-1*sr;
end if;

// Mx:=Cr*M*Transpose(Cc);
Mx:=M;

Mp:=Matrix(GF(p),Matrix(M_orig));
assert Mp eq Submatrix(M,1,1,Nrows(Mp),Ncols(Mp));
assert Mp eq Submatrix(Mx,1,1,Nrows(Mp),Ncols(Mp));


assert Mx[2] eq (Sr*Mx)[2^(sr^-1)];
assert Transpose(Mx)[2] eq (Sc*Transpose(Mx))[2^(sc^-1)];

S:=Sr;
P:=Pr;

// de-correlation permutation:
q:=SymmetricGroup(ncx)![preshuf(x):x in [1..ncx]];
Q:=PermutationMatrix(GF(p),q);
q0:=SymmetricGroup(nc)![preshuf(x):x in [1..nc]];
Q0:=PermutationMatrix(GF(p),q0);

/* XXX Does this hold, really ? */
assert Transpose(Pr*Sr)*Mt*Sc eq Mx*Q;
// , " (this boolean may sometimes be wrong, I think)";
// t means twisted.

/* Note that Mt being equal to Pr times Sr*Cr*M*Transpose(Cc)*Sc^-1 is a
 * direct consequence of the fact that balancing_workhorse.c permutes the
 * rows which are being sent.
 */
if Mt eq P*S*Mx*Q*S^-1 then
    print "Shuffled product detected";
elif  Mt eq S*Mx*Q*S^-1 then
    print "non-shuffled product detected";
    /* This has been deprecated */
    assert false;
    Pr:=Identity(Parent(Pr));
else
    assert false;
end if;
assert Transpose(Pr*Sr)*Mt*Sc eq Mx*Q;
assert Mt eq Pr*Sr*Mx*Q*Sc^-1;

/****************************/

colweights:=func<M|[#[i:i in [1..Nrows(M)]|M[i,j] ne 0]: j in [1..Ncols(M)]]>;
rowweights:=func<M|[#[j:j in [1..Ncols(M)]|M[i,j] ne 0]: i in [1..Nrows(M)]]>;
// colsums:=func<M|[&+[Integers()!M[i,j]:i in [1..Nrows(M)]]: j in [1..Ncols(M)]]>;
// rowsums:=func<M|[&+[Integers()!M[i,j]:j in [1..Ncols(M)]]: i in [1..Nrows(M)]]>;

assert colweights(Mt*Sc) eq colweights(Mx*Q);
// assert rowweights(Transpose(Pr*Sr*Cr)*Mt) eq rowweights(M);
assert rowweights(Transpose(Pr*Sr)*Mt) eq rowweights(M);


assert IsZero((Mx*Q)[nr_orig+1..nr]);
assert IsZero(Transpose(Mx*Q)[nc_orig+1..nc]);



/* The matrix by which we multiply, really */

// assert Submatrix(Transpose(Mx*Q),1,1,nr_orig,nr_orig) eq Submatrix(Transpose(Q),1,1,nr_orig,nr_orig)*Transpose(Matrix(M_orig));

assert Submatrix(Transpose(Mx*Q),1,1,nr_orig,nr_orig) eq Submatrix(Transpose(Q),1,1,nr_orig,nr_orig)*Transpose(Matrix(GF(p),nr_orig,nr_orig,Eltseq(M_orig)));


Qsmall:=Submatrix(Q,1,1,nr_orig,nr_orig);
// Msmall:=Matrix(GF(p),Matrix(M_orig));
Msmall:=Matrix(GF(p),nr_orig,nr_orig,Eltseq(M_orig));

MQsmall:=Msmall*Qsmall;

/* MM is really the matrix with which we are going to work. It's MQsmall in
 * the case nullspace=right, which was the first-class citizen being thought
 * of as this script was written at first.
 */
if nullspace eq "left" then
    MM:=Transpose(MQsmall);
else
    MM:=MQsmall;
end if;



// load "R.m";

// truncvec:=func<x|Vector(Eltseq(x)[1..Maximum(nr_orig, nc_orig)])>;
// expandvec:=func<x|Vector(xx cat [0:i in [#xx+1..nr]]) where xx is Eltseq(x)>;

load "VV.m";
RHS:=VV[1][1..nrhs];

load "x.m";
nchecks:=splitwidth;
Xt:=Matrix(GF(p),m,nr_orig,[]);
for i in [1..m] do for u in var[i] do Xt[i,u] +:=1; end for; end for;

function canvec(V,i) x:=V!0; x[i]:=1; return x; end function;

Xfull:=Matrix(GF(p),#var,nr_orig,[]);
for i in [1..#var] do for u in var[i] do Xfull[i,u] +:=1; end for; end for;

x0:=Rows(Xfull);
print "Dimension of the span of the vectors of X is", Dimension(sub<Universe(x0)|x0>);
xi:=[v*MM : v in x0];


load "Ct.m";
if Category(var) ne RngIntElt then
    Ct:=g(var);
elif assigned Ct then
    delete Ct;
end if;

// if System("test -f Cv0.m") eq 0 then
// magma can't do if (...) then load "..."; (sigh).
load "Cv0.m";
if Category(var) ne RngIntElt then
    Cv0:=VS!g(var);
elif assigned Cv0 then
    delete Cv0;
end if;
//load "C1.m";
//if Category(var) ne RngIntElt then
//    C1:=Vector(GF(p),nr_orig,g(var));
//elif assigned C1 then
//    delete C1;
//end if;
// load "C50.m"; C50:=Vector(GF(p),nr_orig,g(var));

/*
M2:=MM*MM;
M4:=M2*M2;
M8:=M4*M4;
xi:=x0;
xblock:=[Universe(x0)|];
for i in [1..8] do
    xblock cat:= xi;
    xi:=[v*MM : v in xi];
end for;
MX:=Matrix(xblock);


MXN:=MX;
i:=0;
rr:=[Universe(Rows(MX))|];
while i * m lt Nrows(MM) do
    rr cat:=Rows(MXN);
    MXN:=MXN*M8;
    i+:=8;
    printf ".";
    if i le 32 or i mod 128 eq 0 or (i ge 240 and i mod 32 eq 0) then
        printf "\n";
        time U:=sub<Universe(x0)|rr>;
        time d:=Dimension(U);
        print i, d, i*m-d;
    end if;
end while;


// one by one

X0:=Matrix(x0);
Xi:=X0;
xblock:=[Universe(x0)|];
for i in [1..64] do
    xblock cat:= Rows(Xi);
    Xi:=Xi*MM;
    printf ".";
    if #xblock le 32 or #xblock mod 128 eq 0 or (#xblock ge 3600 and #xblock mod 32 eq 0) then
        printf "\n";
        time d:=Rank(Matrix(xblock));
        print #xblock, d, #xblock-d;
    end if;
end for;
for i in [1..64] do
    xblock cat:= Rows(Xi);
    Xi:=Xi*MM;
end for;


Y0:=VV[1];
Yi:=Matrix(Y0);
Yblock:=[Universe(Y0)|];
for i in [1..256] do
    Yblock cat:= Rows(Yi);
    Yi:=Yi*Transpose(MM);
    printf ".";
    if #Yblock le 32 or #Yblock mod 128 eq 0 or (#Yblock ge 2800 and #Yblock mod 32 eq 0) then
        printf "\n";
        time d:=Rank(Matrix(Yblock));
        print #Yblock, d, #Yblock-d;
    end if;
end for;

time for i in [1..32] do
    Yblock cat:= Rows(Yi);
    Yi:=Yi*Transpose(MM);
    printf ".";
end for;


*/


// interval:=50;
// interval:=1;

//if assigned C1 then
//    print "Checking consistency of check vector C1";
//    assert C1 eq xi[1];
//    // C1 eq Vector(Xt)*Transpose(MM);
//    print "Checking consistency of check vector C1: done";
//end if;
//// end if;

if assigned Cv0 then
    printf "Checking consistency of check vector Cv.0-%o.0\n", splitwidth;
    /* Recall that what we think of as rows in the matrix, or also rows
     * (consecutive bits) in the data file correspond to columns in the magma
     * variable. Therefore Ct is rightly put non-transposed in the check below
     *
     * The "file" Ct represents a matrix with m rows and splitwidth columns.
     * Which means that in magma, we're making it a matrix with splitwidth
     * rows and m columns. This fits well with Xt which has m rows.
     */
    assert Cv0 eq Ct*Xt;
    printf "Checking consistency of check vector Cv.0-%o.0: done\n", splitwidth;
end if;

load "Cvi.m";
if Category(var) ne RngIntElt then
    Cvi:=VS!g(var);
elif assigned Cvi then
    delete Cvi;
end if;

load "Cr.m";
if Category(var) ne RngIntElt then
    if p ne 2 and assigned plimbs then
        Cr:=[g(var[j..j+nchecks*plimbs-1]) : j in [1..#var by nchecks*plimbs]];
        assert #Eltseq(Cr[1]) eq 1;
        Cr:=[x[1]:x in Cr];
    else
        Cr:=[g(var[j..j+nchecks*1-1]) : j in [1..#var by nchecks*1]];
        assert Nrows(Cr[1]) eq nchecks;
        assert Ncols(Cr[1]) eq nchecks;
    end if;
elif assigned Cr then
    delete Cr;
end if;

load "Cdi.m";
if Category(var) ne RngIntElt then
    Cdi:=VS!g(var);
elif assigned Cdi then
    delete Cdi;
end if;

check_together:=[];
if assigned Cvi then
    Append(~check_together, Sprintf("Cv.0-%o.%o", splitwidth, interval));
end if;
if assigned Cdi and assigned Cr then
    Append(~check_together, Sprintf("Cd.0-%o.%o", splitwidth, interval));
    Append(~check_together, Sprintf("Cr.0-%o.%o", splitwidth, splitwidth));
end if;
if #check_together gt 0 then
    printf "Checking consistency of check vectors %o\n", check_together;
    foo:=Cv0;
    bar:=VS!0;
    for i in [1..interval] do
        if assigned Cr then
            bar +:= Cr[i]*foo;
        end if;
        foo *:= MM;
    end for;
    if assigned Cvi then assert Cvi eq foo; end if;
    if assigned Cdi and assigned Cr then assert Cdi eq bar; end if;
    printf "Checking consistency of check vectors %o: done\n", check_together;
end if;


/* Note that (below, "where i is 2" always means that the result holds for all
 * i>=0 -- we just want something that is easy to copy-paste.
 *
    Cv0 eq Ct*Xt;
    Cv0*MM^i eq Ct*Xt*MM^i                    where i is 2;

    Xt*MM^i*Transpose(VV[1][1]) eq mcoeff(A,i)   where i is 2;

    &+[Cr[1+i]*Ct*Xt*MM^i:i in [0..interval-1]] eq Cdi;

    &+[Cr[1+i]*Ct*mcoeff(A,i):i in [0..interval-1]] eq Cdi * Transpose(VV[1][1]);
*/

print "Checking consistency of all vectors V_i";
for i in [1..#VV-1] do
    vx:=VerticalJoin(VV[i]);
    for j in [1..interval] do vx := vx*Transpose(MM); end for;
    assert vx eq VerticalJoin(VV[i+1]);
end for;
print "Checking consistency of all vectors V_i: done";





mcoeff:=func<M,k|
ChangeRing(Parent(M)![Coefficient(x,k):x in Eltseq(M)],BaseRing(BaseRing(M)))>;
mdiv:=func<M,k|Parent(M)![k ge 0 select x div BaseRing(M).1^k else x: x in Eltseq(M)]>;
mmulx:=func<M,k|Parent(M)![x * BaseRing(M).1^k: x in Eltseq(M)]>;
mmod:=func<M,k|Parent(M)![x mod BaseRing(M).1^k: x in Eltseq(M)]>;
mval:=func<M|Minimum([Valuation(x):x in Eltseq(M)])>;
mdeg:=func<M|Maximum([Degree(x):x in Eltseq(M)])>;
mcol:=func<M,j|Transpose(M)[j]>;
mrow:=func<M,i|M[i]>;
mrows:=func<M,L|Matrix(BaseRing(M),#L,Ncols(M),[mrow(M,i):i in L])>;
mcols:=func<M,L|Transpose(mrows(Transpose(M),L))>;
function reciprocal(P,k)
    assert k ge Degree(P);
    return Parent(P)![Coefficient(P, k-i):i in [0..k]];
end function;
mrev:=func<M,k|Parent(M)![reciprocal(x,k):x in Eltseq(M)]>;

function mpol_eval(y,M,P)
    /* compute y * P(M) */
    assert CoefficientRing(P) cmpeq CoefficientRing(M);
    s:=Parent(y)!0;
    if Degree(P) lt 0 then return s; end if;
    i:=Degree(P);
    s:=Coefficient(P, i) * y;
    while i gt 0 do
        s:=s*M;
        i:=i-1;
        s+:=Coefficient(P, i) * y;
    end while;
    return s;
end function;

/* Same, but for y a sequence of starting vectors (in the form of a matrix)
 * and P a matrix of polynomials.
 *
 * We want the different vectors in y be matched with the rows of P, so that
 * we'll do a Transpose() behind the scenes.
 * */
function mpol_eval_polmat(y,M,P)
    if splitwidth eq 1 then
        return mpol_eval(y,M,P[1,1]);
    end if;
    /* compute y * P(M) */
    assert BaseRing(CoefficientRing(P)) cmpeq CoefficientRing(M);
    s:=Parent(y)!0;
    if mdeg(P) lt 0 then return s; end if;
    i:=mdeg(P);
    s:=mcoeff(Transpose(P), i) * y;
    while i gt 0 do
        s:=s*M;
        i:=i-1;
        s+:=mcoeff(Transpose(P), i) * y;
    end while;
    return s;
end function;



load "A.m";
function matpol_from_sequence(seq, m, n)
    assert Nrows(seq)*Ncols(seq) mod m*n eq 0;
    deg:=Nrows(seq)*Ncols(seq) div (m*n)-1;
    if p gt 2 then
        enum0:=func<s|[<i-1,s[i]>:i in [1..#s]]>;
        smat:=Matrix(GF(p),deg+1,m*n,Eltseq(seq));
        F:=&+[KP.1^i*Matrix(KP,m,n,Eltseq(r)) where i,r is Explode(ir) : ir in enum0(Rows(smat))];

        // Fr:= &+[KP.1^(deg-i)*Matrix(KP,m,n,seq[i*m*n+1..(i+1)*m*n]):i in [0..deg]];
    else
        // it's quite messy...
        tmp:=Matrix((deg+1)*m, n, Eltseq(Transpose(seq)));
        F:=&+[KP.1^i*Matrix(KP,Submatrix(tmp,i*m+1,1,m,n)):i in [0..deg]];
    end if;
    return F;
end function;
A:=matpol_from_sequence(g(var),m,n);

print "Checking consistency of A";
for k in [0..9] do
    if k*interval gt mdeg(A) then break; end if;
    vx:=VerticalJoin(VV[k+1]);
    assert mcoeff(A, k*interval) eq Matrix(m,n,[(x0[i],vx[j]):j in [1..n], i in[1..m]]);
end for;
print "Checking consistency of A: done";


/* Note that bwc-ptrace already uses g() to parse data which goes to Fchunks
 * and Rchunks */
load "Fchunks.m";
print "Reading generator saved in F.*.sols*";
F:=BlockMatrix(
        n div splitwidth,
        n div splitwidth,
        [matpol_from_sequence(Fchunks[i,j],splitwidth,splitwidth)
            :i,j in [1..n div splitwidth]
        ]);
// F:=Matrix(n,n,[Polynomial(GF(p),Fchunks[i][j]): i,j in [1..n]]);
degsF:={Ncols(Fchunks[i][j]) div splitwidth -1:i,j in [1..n div splitwidth]};
assert #degsF eq 1;
degF:=Setseq(degsF)[1];
Fr:=Parent(F)![reciprocal(x,degF):x in Eltseq(F)];
rhscoeffs:=Matrix(GF(p),#RHS,n,[Polynomial(Eltseq(Rchunks[i][j])):
j in [1..n],
i in [1..#RHS]
]);
print "Reading generator saved in F.*.sols*: done";

degA:=mdeg(A);
coldegs:=[mdeg(mcol(Fr,j)):j in [1..n]];

/* At this point, in the case where no "nrhs" argument has been passed to
 * lingen, the computed matrix F should be such that the following predicates
 * hold.
 */


// degF:=mdeg(F);

if #RHS eq 0 then
    print "Checking generator computed by (p)lingen (homogeneous case)";
    AF:=mmod(mdiv(A,1)*Fr,degA);
    assert exists(e) { e : e in [0..4] | &and[IsZero(mdiv(mcol(AF,j),e+coldegs[j])):j in [1..n]]};
    print "Degree offset for generator is", e;
else
    print "Checking generator computed by (p)lingen (inhomogeneous case)";
    /* In this case, we have operated on a matrix A where only the last
     * columns have been chopped off in the computation.
     */
    r:=#RHS;
    At:=Transpose(A);
    Ax:=Transpose(VerticalJoin(Matrix(At[1..r]), mdiv(Matrix(At[#RHS+1..n]),1)));
    /* And Fr is also a bit peculiar, as some rows of it have been shifted out
     * to form rhscoeffs. The rhscoeffs were in the constant coefficient in F,
     * and therefore in the head coefficients in Fr.
     */

    Ax0:=Submatrix(Ax, 1, 1, n, r);
    Ax1:=Submatrix(Ax, 1, r+1, n, n-r);
    assert Ax0 eq Submatrix(A,1,1,n,r);
    assert Ax1 eq mdiv(Submatrix(A,1,r+1, n, n-r),1);

    Fx:=VerticalJoin(rhscoeffs + mmulx(Matrix(F[1..#RHS]),1),Matrix(F[#RHS+1..n]));
    Fx0:=Submatrix(Fx, 1, 1, r, n);
    Fx1:=Submatrix(Fx, r+1, 1, n-r, n);
    assert Fx0 eq rhscoeffs + X*Submatrix(F, 1, 1, r, n);
    assert Fx1 eq Submatrix(F, r+1, 1, n-r, n);

    // we're a bit lazy here. The top coefficient in Fx0 (and then in Fx)
    // actually only affects the RHS vectors. Therefore we should separate it
    // in the equation AxFx.
    Frx :=Parent(Fx)![reciprocal(x,1+degF):x in Eltseq(Fx)];
    Frx0:=Parent(Fx0)![reciprocal(x,1+degF):x in Eltseq(Fx0)];
    Frx1:=Parent(Fx1)![reciprocal(x,degF):x in Eltseq(Fx1)];
    /* A has an inherent O(X^(degA+1)). Since we've shifted a few columns,
     * it's even O(X^degA).
     */
    AxFx:=mmod(Ax*Frx,degA);
    assert exists(e) { e : e in [0..4] | &and[IsZero(mdiv(mcol(AxFx,j),e+coldegs[j])):j in [1..n]]};
    // when #RHS gt 0, we should always have e==1
    print "Degree offset for generator is", e;
end if;
print "Checking generator computed by lingen: done";

/* ok. Now think about this a bit more. This means we have many zeroes. So
 * something ought to be the null vector, right ?
 */
print "Checking that we have a solution to our system";
if #RHS gt 0 then
    images:=Transpose(rhscoeffs)*Matrix(RHS);
else
    images:=Matrix(GF(p),n, nr_orig, []);
end if;
solutions:=Parent(images)!0;
z := VerticalJoin(VV[1]);
for k in [0..degF] do
    for j in [1..n] do
        solutions[j]+:=mcoeff(mcol(F,j),k)*z;
    end for;
    z:=z*Transpose(MM);
end for;
zz:=images + solutions * Transpose(MM);
MM_nilpotency:=[];
for j in [0..n div splitwidth - 1] do
    j0:=j*splitwidth;
    j1:=j0+splitwidth;
    printf "Checking that sols%o-%o yields a solution... ", j0, j1;
    if IsZero(images) then
        /* Begin with raw solutions first, and check for trivial zeros each
         * time */
        tzz:=Matrix(solutions[j0+1..j1]);
        for c in [0..10] do
            if nullspace eq "right" and IsZero(Eltseq(tzz*Transpose(Qsmall))[1..nc_orig]) then
                printf " ### WARNING: This choice leads to a zero solution vector because of the padding: MM^%o*v is non-zero only on the padded coordinates", c;
                assert false;
            end if;
            tzz:=tzz*Transpose(MM);
            if IsZero(tzz) then
                printf "(for MM^%o) ", c+1;
                Append(~MM_nilpotency, c);
                break;
            end if;
        end for;
    else
        /* Then we can't make sense of this nilpotency thing */
        tzz:=Matrix(zz[j0+1..j1]);
        Append(~MM_nilpotency, 0);
    end if;
    assert IsZero(tzz);
    print " ok";
end for;

load "S.m";
nsols:=#vars * splitwidth;
printf "Computed %o solution vectors with mksol\n", nsols;

// we no longer impose that we compute #RHS solutions. One should be enough
// for everybody.
// assert nsols ge #RHS;

// SSblocks:=[[Vector(GF(p),g(v)):v in w]:w in vars];
SSblocks:=[[VS!g(v):v in w]:w in vars];
SS:=[&+[x:x in y]:y in SSblocks];

function sub_vblock(vv,j)
    if splitwidth eq 1 then
        return vv[1+j];
    else
        return VS!vv[1+j*splitwidth..(j+1)*splitwidth];
    end if;
end function;

function sub_polmat(F,i,j)
    return Submatrix(F,1+i*splitwidth,1+j*splitwidth,splitwidth,splitwidth);
end function;

print "Checking that mksol has computed what we expect";
for i in [0..nsols div splitwidth - 1] do
    vv:=VerticalJoin(V_0);
    for k in [0..#SSblocks[i+1]-1] do
        contribs:=VS!0;
        for j in [0..n div splitwidth - 1] do
            vb:=sub_vblock(vv,j);
            pol:=mmod(mdiv(sub_polmat(F,j,i),k*interval),interval);
            contribs+:=mpol_eval_polmat(vb, Transpose(MM), pol);
        end for;
        assert SSblocks[i+1][k+1] eq contribs;
        for s in [1..interval] do
            vv:=vv*Transpose(MM);
        end for;
    end for;
    assert SS[i+1] eq &+[mpol_eval_polmat(sub_vblock(VerticalJoin(V_0),j),Transpose(MM),sub_polmat(F,j,i)):j in [0..n div splitwidth-1]];
end for;
print "Checking that mksol has computed what we expect: done";

/* all_rhs is the *contribution* of the RHS columns to the sum */
all_rhs:=Transpose(rhscoeffs)*VerticalJoin(RHS);
/* XXX so I think that all_rhs unconditionally has n rows... */
assert Nrows(all_rhs) eq n;

/* here, we're checking only nsols, since we have *computed* only nsols. But
 * had we computed more, we could extend the check somewhat.
 */
for s in [0..nsols div splitwidth - 1] do
    assert IsZero(sub_vblock(all_rhs,s) + (SS[s+1])*Transpose(MM)^(MM_nilpotency[s+1]+1));
end for;

// RHS has the RHS in rows.
Mpad:=HorizontalJoin(Msmall, Transpose(VerticalJoin(RHS)));
ww:=[];
for c in [0..nsols div splitwidth - 1] do
    // v:=&+[mpol_eval(V_0[j],Transpose(MM),F[j,c]): j in [1..n]];
    v:=SS[c+1];
    assert IsZero(sub_vblock(all_rhs,c) + v * Transpose(MM)^(MM_nilpotency[c+1]+1));
    if nullspace eq "right" then
        w0:=v * Transpose(MM)^(MM_nilpotency[c+1]) * Transpose(Qsmall);
        if splitwidth eq 1 then
            w1:=Transpose(rhscoeffs)[c+1];
        else
            w1:=Submatrix(Transpose(rhscoeffs),1+c*splitwidth,1,splitwidth,Nrows(rhscoeffs));
        end if;
        assert IsZero(w0 * Transpose(Msmall) + w1*VerticalJoin(RHS));
        w:=Vector(HorizontalJoin(w0,w1));
    else
        // as we've defined things, the inhomogeneous systems are only though
        // of as M*X=B, not X*M=B. 
        assert #RHS eq 0;
        w:=v * Transpose(MM)^(MM_nilpotency[c+1]);
        assert IsZero(w * Msmall);
    end if;
    Append(~ww, w);
end for;

tr:=nullspace eq "left" select func<x|x> else Transpose;
if nullspace eq "left" then
    equation:="w * M";
elif #RHS gt 0 then
    equation:="(M||RHS) * w";
else
    equation:="M * w";
end if;
printf "Check that we have w such that %o = 0: %o\n",
    equation,
    IsZero(VerticalJoin(ww) * tr(Mpad));


print "Checking that gather has computed what we expect";
load "K.m";
ker:=[Vector(GF(p),g(x)):x in vars];

// Data for all solutions come concatenated, but we have 1+nilpotency order
// each time. Note that when the gather step fails because of a zero solution,
// we miss the last zero-but-non-zero-padding solution.
assert #ker eq nsols div splitwidth + &+MM_nilpotency[1..nsols div splitwidth];
ker_bits:=[];
stowed:=0;
for s in [0..nsols div splitwidth - 1] do
    ker_chunk:=ker[stowed + 1..stowed + 1 + MM_nilpotency[s+1]];
    early:=ker_chunk[1];
    late:=ker_chunk[#ker_chunk];

    if nullspace eq "right" then
        /* There are two sets of asserts here, and it's admittedly somewhat
         * confusing.  Depending on whether it succeeded or not, gather may
         * put varying stuff in the K.sols file.
         */
        if Dimension(Parent(early)) eq Ncols(Msmall) + #RHS then
            /* In contrast, when gather was happy, the K file has nrhs extra
             * coordinates, so that the following holds...
             */
            assert IsZero(Matrix(late)*Transpose(HorizontalJoin(Msmall, Transpose(Matrix(RHS)))));
        elif Dimension(Parent(early)) eq Ncols(Msmall) then
            /* It may occur that gather wasn't happy with the solution it got
             * for the inhomogenous linear system. In that case, the binary K
             * file misses the rhs coefficients, which is weird. So then, we
             * just have the following. And most probably this is all going to
             * fail.
             */
            assert early eq Matrix(SS[s+1])*Transpose(Qsmall);
            assert false; // don't know how to make sense of the next assert. Anyway the gather program has failed at this point.
            /*
            assert IsZero(
                HorizontalJoin(Matrix(late),Matrix(Transpose(rhscoeffs)[1..#ker])) *
                VerticalJoin(Transpose(Msmall), Matrix(RHS)));
                */
        else
            /* uh oh */
            error "I don't understand the dimension of the sequence computed by gather in K.m";
        end if;
        /* For completeness, note that the ascii K file is almost identical to
         * the above except that we acknowledge there the presence of zero
         * columns in M, so that a handful of coefficients in ker are useless,
         * and hence chopped out.
         */
    else
        /* Then it's simpler, because we have no RHS to mess things up */
        // assert IsZero(VerticalJoin(late)*Msmall);
        assert IsZero(late*Msmall);
    end if;
end for;
print "Checking that gather has computed what we expect: done";

