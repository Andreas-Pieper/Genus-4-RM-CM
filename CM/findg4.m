declare verbose CM, 1;

intrinsic FullEnumerationG4(f::RngUPolElt : prec := 300, exp := Infinity(), FixCMType := false) -> .
  {Finds the equivalence classes of CM curves given a polynomial.}

/* Global structures */
F := RationalsExtra(prec); CC := F`CC; R<x> := PolynomialRing(F);
K := NumberFieldExtra(f);
Phis := AllCMTypesUpToEquivalence(K : Primitive := true);
//Phis := AllCMTypes(K : Primitive := true);
precsmall := 50; CCSmall := ComplexFieldExtra(precsmall);

/* Determine totally real subfield and corresponding involution */
test, tup, inv := IsCMField(K);
K0, hK0K := Explode(tup);
assert inv(inv(K.1)) eq K.1;
K0sub := sub< K | hK0K(K0.1) >;

/* Determine class group */
ZZK := Integers(K);
Diff := Different(ZZK);
//aas := ClassGroupSmallRepresentatives(ZZK, exp);
//vprint CMExp, 1 : "Done determining small representatives of class group!", #aas;

/* Determine unit group quotients */
U, phiU := UnitGroup(ZZK : GRH := true);
//for i in [1..NumberOfGenerators(U)] do print [ EmbedExtra(phiU(U.1) : iota := inf) : inf in InfinitePlacesExtra(K) ]; end for;
gensU := Generators(U);
gensV := [ (phiU(gen)*inv(phiU(gen))) @@ phiU : gen in gensU ];
V := sub< U | gensV >;
Q, pQ := quo< U | V >;
vprint CMExp, 1 : "Done determining quotient of unit group!", #Q;

/* Now follow Streng's method */
taus := [ ];
counter := 0;
aas := [ideal< ZZK | 1>]; // FIXME WARNING! assuming trivial class group
for aa in aas do
    vprint CMExp, 1 : "Norm of ideal:", Norm(aa);
    aabar := ideal< ZZK | [ inv(gen) : gen in Generators(aa) ] >;

    /* Need generator; can also exclude fields with i here */
    test1, xi0 := IsPrincipal((aa*aabar*Diff)^(-1));
    if not test1 then
        vprint CMExp, 1 : "Stopping at test1";
        continue;
    end if;
    //print [ EmbedExtra(xi0 : iota := inf) : inf in InfinitePlacesExtra(K) ];

    for q in Q do
        u := q @@ pQ; xi := K ! (phiU(u)*xi0);

        /* Insist that xi be totally imaginary */
        test2 := xi^2 in K0sub and not xi in K0sub;
        if not test2 then
            vprint CMExp, 1 : "Stopping at test2";
            continue;
        end if;
        assert ideal< ZZK | xi > eq (aa*aabar*Diff)^(-1);

        /* Find corresponding CM type if possible */
        test3 := false;
        if FixCMType then NPhis := 1; else NPhis := #Phis; end if;
        for iPhi in [1..NPhis] do
            Phi := Phis[iPhi];
            embs := [ EmbedExtra(xi : iota := inf) : inf in Phi ];
            if &and[ Im(emb) gt CC`epscomp : emb in embs ] then
                test3 := true;
                break;
            end if;
        end for;
        if not test3 then
            vprint CMExp, 1 : "Stopping at test3";
            continue;
        end if;

        /* NOTE that the tests above costs time and are stupid */
        /* If all tests pass, recover period matrix */
        B := Basis(aa);
        B := [ K ! b : b in B ];
        P := Matrix(CC, [ [ EmbedExtra(b : iota := iota) : b in B ] : iota in Phi ]);
        P := ChangeRing(P, CC);
        //P0, U := MakeSmallPeriodMatrix(P);
        //print U;
        //PrintFileMagma("Ps1.m", P);

        /* Recover principal polarization */
        E := Matrix([ [ Trace(xi*inv(x)*y) : x in B ] : y in B ]);
        E := ChangeRing(E, Rationals());
        E := -E;
        /*
        test4 := false;
        if IsPolarization(E, P) then
            test4 := true;
        elif IsPolarization(-E, P) then
            test4 := true;
            E := -E;
        end if;
        */

        /* Declare success */
        counter +:= 1;
        vprint CMExp, 1 : "New abelian variety number", counter;

        /* Convert to big period matrix */
        E := ChangeRing(E, Integers());
        E0, T := FrobeniusFormAlternating(E);
        P := P*ChangeRing(Transpose(T), CC);
        //assert IsBigPeriodMatrix(P);
        //PrintFileMagma("Ps2.m", P);

        /* Convert to small period matrix */
        tau := SmallPeriodMatrix(P);
        /*
        print BaseRing(tau);
        print ChangeRing(tau, CCSmall);
        print Maximum([ Abs(c) : c in Eltseq(tau - Transpose(tau)) ]);
        A := Matrix([ [ Im(c) : c in Eltseq(row) ] : row in Rows(tau) ]);
        print Eigenvalues(A);
        */

        /* Reduce small period matrix */
        vprint CMExp, 1 : "Reducing period matrix...";
        tau := ReduceSmallPeriodMatrix(tau);
        vprint CMExp, 1 : "done with reduction.";
        tausmall := ChangeRing(tau, CCSmall);
        Append(~taus, <tau, Phi>);
    end for;
end for;
return taus;
end intrinsic;
