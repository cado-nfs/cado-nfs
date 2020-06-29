SetColumns(0);
ps:=func<x|Join([Sprintf("%o",y):y in x], " ")>;
for i in [1..40] do
    d:=Random(10)+2;
    repeat
        f:=Polynomial([Random(-100,100):i in [0..d]]);
    until IsIrreducible(f) and Degree(f) ge 2;
    disc:=LeadingCoefficient(f)*Discriminant(f);
    primes:={ pi[1]:pi in TrialDivision(disc) |
        pi[1] le 2^63 and Valuation(disc,pi[1]) ge 2 };
    primes join:={2,3,5,7};
    for p in primes do
        K:=NumberField(f);
        E:=Order([LeadingCoefficient(f)*K.1]);
        O:=pMaximalOrder(E, p);
        B:=BasisMatrix(O);
        Bd:=Denominator(B);
        Bz:=Matrix(Integers(),B*Bd);
        O`Maximal:=true;
        Fp:=Factorization(ideal<O|p>);
        print ps(Eltseq(f));
        iprint:=function(Im)
            I,m:=Explode(Im);
            // MI:=BasisMatrix(I);
            // return Seqlist(Eltseq(BasisMatrix(I))) cat [* "", m *];
            g0,g1 := TwoElement(I);
            assert g0 eq p;
            g1:=K!g1;
            d1:=Denominator(g1);
            return Seqlist([d1] cat Eltseq(d1*g1)) cat [* "", m *];
        end function;
        samples:=[];
        for k in [1..20] do
            gens:=[K![Random(-5,5):i in [1..Degree(f)]]:i in [1..2]];
            // magma bug.
            I:=&+[ideal<O|x>:x in gens];
            // this would be nice, but unfortunately we can't do that since we
            // don't really have a maximal order.
            // I*:=Random(Fp)[1]^Random(-2,2);
            // gens:=[a,b] where a,b is TwoElement(I);
            vals:=[*Valuation(I,fkp[1]):fkp in Fp*];
            if k le 2 or not &and [x eq 0:x in vals] then
                Append(~samples, [*"composite"*] cat Seqlist([#gens] cat &cat [Eltseq(g):g in gens]) cat [*"valuations"*] cat [*Valuation(I,fkp[1]):fkp in Fp*]);
                if #samples ge 5 then break; end if;
            end if;
        end for;
        Append(~samples, [*"composite", 0*]);
        print ps([p]),
            ps(["order"]),
            ps([Bd] cat Eltseq(Bz)),
            ps(["ideals"]),
            ps([*#Fp*] cat &cat [[*""*] cat iprint(Im):Im in Fp]),
            Join([ps(x):x in samples], " ");
    end for;
end for;

