declare verbose CM, 1;

declare verbose CMExp, 2;

// does full enumeration, not taking Galois action into account
intrinsic FullEnumerationG4(f::RngUPolElt : prec := 300, exp := Infinity(), FixCMType := false) -> .
  {Finds the equivalence classes of CM curves given a polynomial.}

/* Global structures */
F := RationalsExtra(prec); CC := F`CC; R<x> := PolynomialRing(F);
K := NumberFieldExtra(f : prec := prec);
//printf "K done";
Phis := AllCMTypesUpToEquivalence(K : Primitive := true);
//printf "Phis done";
//Phis := AllCMTypes(K : Primitive := true);
precsmall := 50; CCSmall := ComplexFieldExtra(precsmall);

/* Determine totally real subfield and corresponding involution */
test, tup, inv := IsCMField(K);
//printf "IsCMField done";
K0, hK0K := Explode(tup);
assert inv(inv(K.1)) eq K.1;
K0sub := sub< K | hK0K(K0.1) >;

/* Determine class group */
ZZK := Integers(K);
Diff := Different(ZZK);
//printf "Diff computed";
//aas := ClassGroupSmallRepresentatives(ZZK, exp);
//vprint CMExp, 1 : "Done determining small representatives of class group!", #aas;

/* Determine unit group quotients */
U, phiU := UnitGroup(ZZK : GRH := true);
//printf "Unit group done";
//for i in [1..NumberOfGenerators(U)] do print [ EmbedExtra(phiU(U.1) : iota := inf) : inf in InfinitePlacesExtra(K) ]; end for;
gensU := Generators(U);
//printf "generators U done";
gensV := [ (phiU(gen)*inv(phiU(gen))) @@ phiU : gen in gensU ];
V := sub< U | gensV >;
Q, pQ := quo< U | V >;
//printf "quotients done";
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
    //printf "is principal done";
    if not test1 then
        vprint CMExp, 1 : "Stopping at test1";
        continue;
    end if;
    //print [ EmbedExtra(xi0 : iota := inf) : inf in InfinitePlacesExtra(K) ];

    comp := 0;
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
            //"test3";
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

// copy-pasta'd from findg3.m
import "reflex.m": DescendIdeal;
import "smallcl.m": MakeSmallSequence;


/* First functions copied from curve_reconstruction */

/* The first two functions are for compatibility with Labrande's conventions */
/* The vectors returned are column vectors, not because I think that this is so
 * great, but because this respects Magma's conventions */
function VectorFromIndex(g, i)
j := i; v := [ ];
for k:=1 to 2*g do
    Append(~v, (j mod 2)/2);
    j := j div 2;
end for;
return Transpose(Matrix([Reverse(v)]));
end function;


/* Second return value is sign transition from standard value */
function IndexFromVector(g, v);
w := Reverse(Eltseq(v));
i := &+[ ((Integers() ! (2*w[j])) mod 2)*2^(j - 1) : j in [1..(2*g)] ];
if i eq 0 then
    return 2^(2*g);
end if;
v0 := VectorFromIndex(g, i); n := v - v0;
v0 := Eltseq(v0); n := Eltseq(n);
s := &+[ v0[i]*n[g + i] : i in [1..g] ];
if IsIntegral(s) then
    s := 1;
else
    s := -1;
end if;
return i, s;
end function;


/* Checks if a vector is even */
function IsEvenVector(v)
w := Eltseq(v); g := #w div 2;
w1 := w[1..g]; w2 := w[(g + 1)..(2*g)];
return IsZero((Integers() ! (4 * &+[ w1[i]*w2[i] : i in [1..g] ])) mod 2);
end function;


/* Finds even vanishing thetas */
function FindDelta(thetas : prectest := 50)
CP := CartesianPower([0,1/2], 8);
vs := [ Transpose(Matrix([[ c : c in tup ]])) : tup in CP ];
v0s := [ ];
for v in vs do
    theta := thetas[IndexFromVector(3, v)];
    test := (Abs(theta) lt 10^(-prectest)) and IsEvenVector(v);
    if test then
        Append(~v0s, v);
    end if;
end for;
return v0s;
end function;


function MakeSmallElement(alpha, K0info : precred := 200)
/* We assume the number field involved to be totally real */

K0, ZZK0, DiffK0, infsK0, ClK0info, UK0info, hK0K := Explode(K0info); ClK0, phiClK0 := Explode(ClK0info); UK0, phiUK0 := Explode(UK0info);
assert IsTotallyReal(K0);

/* Build lattice and find closest vectors, ignoring -1 */
rows := [ ];
for i:=2 to NumberOfGenerators(UK0) do
    u := phiUK0(UK0.i);
    row := [ Log(Abs(Evaluate(u, infK0 : Precision := precred))) : infK0 in infsK0 ];
    Append(~rows, row);
end for;
M := Matrix(rows);
Lat := Lattice(M);
v := Vector([ Log(Abs(Evaluate(alpha, infK0 : Precision := precred))) : infK0 in infsK0 ]);
w := ClosestVectors(Lat, v)[1];

/* Find integral combination and correct */
//M0 := Submatrix(M, 1,1, n-1,n-1);
//w0 := Vector(Eltseq(w)[1..(n-1)]);
//x := NumericalSolution(M0, w0);
x := NumericalSolution(M, w);
x := [ Round(c) : c in Eltseq(x) ];
/* Once more we ignore -1 */
u := (phiUK0(UK0 ! ([0] cat x)))^(-1);
alpha0 := alpha*u^(-1);

/* Final result is a small alpha0 as well as a unit u such that alpha = u alpha0 */
return K0 ! alpha0, K0 ! u;

end function;


function MakeSmallUnit(alpha, B, K0info : precred := 200)
/* B is a basis of a sublattice of the unit group */
/* We assume the number field involved to be totally real */

K0, ZZK0, DiffK0, infsK0, ClK0info, UK0info, hK0K := Explode(K0info); ClK0, phiClK0 := Explode(ClK0info); UK0, phiUK0 := Explode(UK0info);
assert IsTotallyReal(K0);

/* Build lattice and find closest vectors */
rows := [ ];
for i:=1 to #B do
    u := phiUK0(B[i]);
    row := [ Log(Abs(Evaluate(u, infK0 : Precision := precred))) : infK0 in infsK0 ];
    Append(~rows, row);
end for;
M := Matrix(rows);
Lat := Lattice(M);
v := Vector([ Log(Abs(Evaluate(alpha, infK0 : Precision := precred))) : infK0 in infsK0 ]);
w := ClosestVectors(Lat, v)[1];

/* Find integral combination and correct */
//M0 := Submatrix(M, 1,1, n-1,n-1);
//w0 := Vector(Eltseq(w)[1..(n-1)]);
//x := NumericalSolution(M0, w0);
x := NumericalSolution(M, w);
x := [ Round(c) : c in Eltseq(x) ];
u := phiUK0(&+[ x[i]*B[i] : i in [1..#B] ]);
alpha0 := alpha*u^(-1);

/* Final result is a small alpha0 as well as a unit u such that alpha = u alpha0 */
return K0 ! alpha0, K0 ! u;

end function;


function MakeSmallXi(xi, B, Kinfo, K0info : precred := 200)
/* B is a basis of a sublattice of the unit group */
/* We assume the number field involved to be totally real */

K, ZZK, DiffK, infsK, ClKinfo, UKinfo, invK := Explode(Kinfo); ClK, phiClK := Explode(ClKinfo); UK, phiUK := Explode(UKinfo);
K0, ZZK0, DiffK0, infsK0, ClK0info, UK0info, hK0K := Explode(K0info); ClK0, phiClK0 := Explode(ClK0info); UK0, phiUK0 := Explode(UK0info);
assert IsTotallyReal(K0);

/* Build lattice and find closest vectors */
rows := [ ];
for i:=1 to #B do
    u := hK0K(phiUK0(B[i]));
    row := [ Log(Abs(Evaluate(u, infK : Precision := precred))) : infK in infsK ];
    Append(~rows, row);
end for;
M := Matrix(rows);
Lat := Lattice(M);
v := Vector([ Log(Abs(Evaluate(xi, infK : Precision := precred))) : infK in infsK ]);
w := ClosestVectors(Lat, v)[1];

/* Find integral combination and correct */
//M0 := Submatrix(M, 1,1, n-1,n-1);
//w0 := Vector(Eltseq(w)[1..(n-1)]);
//x := NumericalSolution(M0, w0);
x := NumericalSolution(M, w);
x := [ Round(c) : c in Eltseq(x) ];
u := hK0K(phiUK0(&+[ x[i]*B[i] : i in [1..#B] ]));
xi0 := xi*u^(-1);

/* Final result is a small xi0 as well as a unit u such that xi = u xi0 */
return K ! xi0, K ! u;

end function;


function IsPrincipalSmall(aa, K0info : precred := 200)
/* We assume the number field involved to be totally real */

K0, ZZK0, DiffK0, infsK0, ClK0info, UK0info, hK0K := Explode(K0info); ClK0, phiClK0 := Explode(ClK0info); UK0, phiUK0 := Explode(UK0info);
assert IsTotallyReal(K0);
test, gen := IsPrincipal(aa);
if not test then return test, 0; end if;
return test, MakeSmallElement(gen, K0info : precred := precred);

end function;


function SmallRepresentativesIdeals(Kinfo, K0info : exp := 4, precred := 200)
/* Determines generators of subquotient by exp */

/* Recover information */
K, ZZK, DiffK, infsK, ClKinfo, UKinfo, invK := Explode(Kinfo); ClK, phiClK := Explode(ClKinfo); UK, phiUK := Explode(UKinfo);
K0, ZZK0, DiffK0, infsK0, ClK0info, UK0info, hK0K := Explode(K0info); ClK0, phiClK0 := Explode(ClK0info); UK0, phiUK0 := Explode(UK0info);
if (#ClK eq 1) then return [ ideal< ZZK | 1 > ], [ ideal< ZZK | 1 > ]; end if;
nClK := NumberOfGenerators(ClK);
gensClK := [ MakeSmallClassRepresentative(phiClK(ClK.i), infsK : precred := precred) : i in [1..nClK] ];

/* First restrict to classes of ideals for which aa aabar is principal in K0 */
gensh1 := [ ];
for i in [1..nClK] do
    aa := gensClK[i];
    aabar := ideal< ZZK | [ invK(genaa) : genaa in Generators(aa) ] >;
    aaaabar0 := DescendIdeal(aa*aabar, ZZK0, hK0K);
    genh1 := aaaabar0 @@ phiClK0;
    Append(~gensh1, genh1);
end for;
h1 := hom< ClK -> ClK0 | gensh1 >;
ClK1 := Kernel(h1);
if (#ClK1 eq 1) then return [ ideal< ZZK | 1 > ], [ ideal< ZZK | 1 > ]; end if;

/* TODO: Sanity check */
for i in [1..NumberOfGenerators(ClK1)] do
    aa := MakeSmallSequence(Eltseq(ClK ! ClK1.i), gensClK, infsK : precred := precred);
    aabar := ideal< ZZK | [ invK(genaa) : genaa in Generators(aa) ] >;
    aaaabar0 := DescendIdeal(aa*aabar, ZZK0, hK0K);
    assert (aa @@ phiClK) eq (ClK ! ClK1.i);
    assert IsPrincipal(aaaabar0);
end for;

/* Then restrict classes of ideals for which aa aabar is principal in K0 and generated by a totally positive element */
/* To this end we first determine the signs to which units give rise and take the relevant quotient */
C2 := CyclicGroup(2); C2prod, incs := DirectProduct([ C2 : i in [1..#infsK0] ]);
gensC2prodsub := [ ]; Msign := [ ];
for i in [1..NumberOfGenerators(UK0)] do
    u0 := phiUK0(UK0.i);
    seq := [ ]; Mrow := [ ];
    for inf0 in infsK0 do
        sg := Sign(Evaluate(u0, inf0 : Precision := precred));
        if sg eq 1 then Append(~seq, C2.1^0); else Append(~seq, C2.1^1); end if;
        if sg eq 1 then Append(~Mrow, 0); else Append(~Mrow, 1); end if;
    end for;
    gen := &*[ incs[i](seq[i]) : i in [1..#infsK0] ];
    Append(~gensC2prodsub, gen); Append(~Msign, Mrow);
end for;
Q, phiQ := quo< C2prod | gensC2prodsub >;
Msign := Matrix(FiniteField(2), Msign);

/* Now take images of generators and form kernel */
gensh2 := [ ];
for i in [1..NumberOfGenerators(ClK1)] do
    aa := MakeSmallSequence(Eltseq(ClK ! ClK1.i), gensClK, infsK : precred := precred);
    aabar := ideal< ZZK | [ invK(genaa) : genaa in Generators(aa) ] >;
    aaaabar0 := DescendIdeal(aa*aabar, ZZK0, hK0K);

    test, gen := IsPrincipalSmall(aaaabar0, K0info : precred := precred);
    assert test;
    seq := [ ];
    for inf0 in infsK0 do
        sg := Sign(Evaluate(gen, inf0 : Precision := precred));
        if sg eq 1 then Append(~seq, C2.1^0); else Append(~seq, C2.1^1); end if;
    end for;
    genh2 := phiQ(&*[ incs[i](seq[i]) : i in [1..#infsK0] ]);
    Append(~gensh2, genh2);
end for;
h2 := hom< ClK1 -> Q | gensh2 >;
ClK2 := Kernel(h2);

/* Form representatives of ClK1 / ClK2 */
Q, phiQ := quo< ClK1 | ClK2 >;
expQ := Exponent(Q);
aas := [ ];
for q in Q do
    aa := q @@ phiQ;
    aa := ClK1 ! [ c mod expQ : c in Eltseq(aa) ];
    aa := MakeSmallSequence(Eltseq(ClK ! aa), gensClK, infsK : precred := precred);
    Append(~aas, aa);
end for;

/* Simple case */
if IsOdd(#ClK2) then return aas, [ ideal< ZZK | 1 > ]; end if;

/* TODO: Sanity check */
for i in [1..NumberOfGenerators(ClK2)] do
    bb := MakeSmallSequence(Eltseq(ClK ! ClK1 ! ClK2.i), gensClK, infsK : precred := precred);
    bbbar := ideal< ZZK | [ invK(genbb) : genbb in Generators(bb) ] >;
    bbbbbar0 := DescendIdeal(bb*bbbar, ZZK0, hK0K);
    test, eta := IsPrincipalSmall(bbbbbar0, K0info : precred := precred);
    assert test; assert ideal< ZZK0 | eta > eq bbbbbar0;
    v := [ ];
    for inf0 in infsK0 do
        sg := Sign(Evaluate(eta, inf0 : Precision := precred));
        if sg eq 1 then Append(~v, 0); else Append(~v, 1); end if;
    end for;
    v := Vector(FiniteField(2), v);
    /* The fact that this does not give an error means that the sanity check passed */
    u := Solution(Msign, v);
end for;

/* Choose relevant generators */
gens2 := [ ClK ! ClK1 ! ClK2.i : i in [1..NumberOfGenerators(ClK2)] ];
gens2 := [ gen2 : gen2 in gens2 | GCD(Order(gen2), exp) ne 1 ];

/* Construct small generators */
gcds2 := [ GCD(Order(gen2), exp) : gen2 in gens2 ];
CP := CartesianProduct([ [0..(gcd2 - 1)] : gcd2 in gcds2 ]);
gens2 := [ MakeSmallSequence(Eltseq(gen2), gensClK, infsK : precred := precred): gen2 in gens2 ];

/* Construct small ideals */
tupold := 0; runningtupold := 0; bbold := 0; bbs := [ ];
counter := 0;
for tup in CP do
    if &and[ c eq 0 : c in tup ] then
        runningtup := [ ideal< ZZK | 1 > : gen2 in gens2 ];
        Append(~bbs, ideal< ZZK | 1 >);
    else
        i := Minimum([ j : j in [1..#tup] | tup[j] ne tupold[j] ]);
        runningtup[i] := MakeSmallClassRepresentative(runningtup[i]*gens2[i], infsK : precred := precred);
        for j in [(i + 1)..#tup] do runningtup[j] := ideal< ZZK | 1 >; end for;
        Append(~bbs, MakeSmallClassRepresentative(&*runningtup, infsK : precred := precred));
    end if;
    tupold := tup; runningtupold := runningtup; bbold := bbs[#bbs];
end for;
return aas, bbs;

end function;


function SmallRepresentativesUnits(Kinfo, K0info : precred := 200)
/* Determines generators of subquotient by norms */

/* Recover information */
K, ZZK, DiffK, infsK, ClKinfo, UKinfo, invK := Explode(Kinfo); ClK, phiClK := Explode(ClKinfo); UK, phiUK := Explode(UKinfo);
K0, ZZK0, DiffK0, infsK0, ClK0info, UK0info, hK0K := Explode(K0info); ClK0, phiClK0 := Explode(ClK0info); UK0, phiUK0 := Explode(UK0info);

/* Create codomain of totally positive units */
C2 := CyclicGroup(2); C2prod, incs := DirectProduct([ C2 : i in [1..#infsK0] ]);
gensh1 := [ ];
for i in [1..NumberOfGenerators(UK0)] do
    u0 := phiUK0(UK0.i);
    seq := [ ];
    for inf0 in infsK0 do
        sg := Sign(Evaluate(u0, inf0 : Precision := precred));
        if sg eq 1 then Append(~seq, C2.1^0); else Append(~seq, C2.1^1); end if;
    end for;
    genh1 := &*[ incs[i](seq[i]) : i in [1..#infsK0] ];
    Append(~gensh1, genh1);
end for;
h1 := hom< UK0 -> C2prod | gensh1 >;
UK1 := Kernel(h1);

/* Build norm map and its image */
gensNmKK0 := [ ];
for i in [1..NumberOfGenerators(UK)] do
    u := phiUK(UK.i);
    u0 := (u*invK(u)) @@ hK0K;
    Append(~gensNmKK0, UK1 ! (u0 @@ phiUK0));
end for;
NmKK0 := hom< UK -> UK1 | gensNmKK0 >;
ImKK0 := Image(NmKK0);

/* We need this for minimization */
Q, phiQ := quo< UK0 | UK1 >;
expQ := Exponent(Q);
B1 := [ UK0 ! UK1.i : i in [1..NumberOfGenerators(UK1)] ];

/* Take preimages and ensure that these are small */
u0s := [ ];
for q in Q do
    u0 := q @@ phiQ;
    u0 := UK0 ! [ c mod expQ : c in Eltseq(u0) ];
    u0 := phiUK0(u0);
    u0 := MakeSmallUnit(u0, B1, K0info : precred := precred);
    Append(~u0s, u0);
end for;

/* We need this for minimization */
Q, phiQ := quo< UK1 | ImKK0 >;
expQ := Exponent(Q);
BImKK0 := [ UK0 ! UK1 ! ImKK0.i : i in [1..NumberOfGenerators(ImKK0)] ];

/* Take preimages and ensure that these are small */
v0s := [ ];
for q in Q do
    v0 := q @@ phiQ;
    v0 := UK1 ! [ c mod expQ : c in Eltseq(v0) ];
    v0 := phiUK0(UK0 ! v0);
    v0 := MakeSmallUnit(v0, BImKK0, K0info : precred := precred);
    Append(~v0s, v0);
end for;
assert &and[ &and[ Evaluate(v0, inf : Precision := precred) gt 0 : inf in infsK0 ] : v0 in v0s ];

return u0s, v0s;

end function;


function SingleClassFromPolynomial(f : exp := 4, prec := 300, precred := 200)
/* We follow Streng's argument in Proposition 4.4 of "Computing Igusa class polynomials" and optimize along the way  */

/* Create CM field and totally real subfield */
F := RationalsExtra(prec); R := PolynomialRing(F);
K := NumberFieldExtra(R ! f); CC := K`CC;
test, tup, invK := IsCMField(K);
K0, hK0K := Explode(tup);

vprint CMExp, 2 : "Determining field information...";
/* Create rings of integers */
ZZK := Integers(K); ZZK0 := Integers(K0);
/* Create class groups */
ClK, phiClK := ClassGroup(ZZK); ClK0, phiClK0 := ClassGroup(ZZK0);
vprint CMExp, 2 : "Structure of Cl (K):", ClK;
vprint CMExp, 2 : "Structure of Cl (K0):", ClK0;
/* Create unit groups */
UK, phiUK := UnitGroup(ZZK : GRH := true); UK0, phiUK0 := UnitGroup(ZZK0 : GRH := true);
/* Create differents */
DiffK := Different(ZZK); DiffK0 := Different(ZZK0);
/* Create infinite places */
infsK := InfinitePlaces(K); infsK0 := InfinitePlaces(K0);
/* Store */
Kinfo := [* K, ZZK, DiffK, infsK, [* ClK, phiClK *], [* UK, phiUK *], invK *];
K0info := [* K0, ZZK0, DiffK0, infsK0, [* ClK0, phiClK0 *], [* UK0, phiUK0 *], hK0K *];
vprint CMExp, 2 : "done.";

/* Calculate small representatives */
vprint CMExp, 2 : "Determining small representatives...";
nClK := NumberOfGenerators(ClK);
gensClK := [ MakeSmallClassRepresentative(phiClK(ClK.i), infsK : precred := precred) : i in [1..nClK] ];
aas, bbs := SmallRepresentativesIdeals(Kinfo, K0info : exp := exp, precred := precred);
vprint CMExp, 2 : "Ideals:", #bbs;
u0s, v0s := SmallRepresentativesUnits(Kinfo, K0info : precred := precred);
assert &and[ Abs(Norm(u0)) eq 1 : u0 in u0s ];
assert &and[ Norm(v0) eq 1 : v0 in v0s ];
vprint CMExp, 2 : "Units:", #v0s;
vprint CMExp, 2 : "done.";

/* Creation of z */
z := K.1 - invK(K.1);
assert z ne 0;
/* Creation of ideal hh in K0 */
B := Basis(ZZK);
hh0 := ideal< ZZK0 | [ (K !(z*(b - invK(b)))) @@ hK0K : b in B ] >;

/* Special case of trivial class group */
if #ClK eq 1 then
    assert #ClK0 eq 1;
    aa00 := ideal< ZZK | [ ZZK ! 1 ] >;
    test, y := IsPrincipalSmall(hh0^(-1)*DiffK0^(-1), K0info : precred := precred);
    assert test;
    xi000 := hK0K(y)*z;

    /* Determine suitable CM type */
    Phis := AllCMTypes(K : Primitive := true);
    for aa in aas do
        aabar := ideal< ZZK | [ invK(genaa) : genaa in Generators(aa) ] >;
        aaaabar := DescendIdeal(aa*aabar, ZZK0, hK0K);
        test, gen := IsPrincipalSmall(aaaabar, K0info : precred := precred);
        aa0 := aa00*aa;
        aa0bar := ideal< ZZK | [ invK(genaa0) : genaa0 in Generators(aa0) ] >;
        xi00 := xi000*hK0K(gen^(-1));

        for u0 in u0s do
            xi0 := xi00*hK0K(u0);
            for Phi0 in Phis do
                if &and[ Im(EmbedExtra(xi0 : iota := iota)) gt CC`epscomp : iota in Phi0 ] then
                    /* Check and return */
                    assert IsTotallyPositive((-xi0^2) @@ hK0K);
                    assert ideal< ZZK | xi0 > eq (aa0*aa0bar*DiffK)^(-1);
                    return true, [* Phi0, aa0, xi0 *], Kinfo, K0info, bbs, v0s;
                end if;
            end for;
        end for;
    end for;
    return false, 0, 0, 0, 0, 0;
end if;

/* Create norm map if the class group is non-trivial */
vprint CMExp, 2 : "Creating norm map...";
gensNmKK0 := [ ];
for i in [1..nClK] do
    aa := gensClK[i];
    aabar := ideal< ZZK | [ invK(genaa) : genaa in Generators(aa) ] >;
    aaaabar0 := DescendIdeal(aa*aabar, ZZK0, hK0K);
    Append(~gensNmKK0, aaaabar0 @@ phiClK0);
end for;
NmKK0 := hom< ClK -> ClK0 | gensNmKK0 >;
assert Image(NmKK0) eq ClK0;
vprint CMExp, 2 : "done.";

/* Take a small preimage of hh0 */
vprint CMExp, 2 : "Determining small class representative...";
Q := quo< ClK | Kernel(NmKK0) >;
expQ := Exponent(Q);
/* Note that indeed an inverse is taken, in the form of a minus sign */
seq00 := [ -c mod expQ : c in Eltseq(((hh0*DiffK0) @@ phiClK0) @@ NmKK0) ];

/* Make representative small using geometry of numbers */
aa00 := MakeSmallSequence(seq00, gensClK, infsK : precred := precred);
aa00bar := ideal< ZZK | [ invK(genaa00) : genaa00 in Generators(aa00) ] >;
aa00aa00bar0 := DescendIdeal(aa00*aa00bar, ZZK0, hK0K);
vprint CMExp, 2 : "done.";

/* Determine x0 (minimize?) */
vprint CMExp, 2 : "Determining small generator...";
test, y := IsPrincipalSmall(hh0^(-1)*DiffK0^(-1)*aa00aa00bar0^(-1), K0info : precred := precred);
assert test;
xi000 := hK0K(y)*z;
vprint CMExp, 2 : "done.";

/* Determine suitable CM type */
Phis := AllCMTypes(K : Primitive := true);
for aa in aas do
    aabar := ideal< ZZK | [ invK(genaa) : genaa in Generators(aa) ] >;
    aaaabar := DescendIdeal(aa*aabar, ZZK0, hK0K);
    test, gen := IsPrincipalSmall(aaaabar, K0info : precred := precred);
    aa0 := aa00*aa;
    aa0bar := ideal< ZZK | [ invK(genaa0) : genaa0 in Generators(aa0) ] >;
    xi00 := xi000*hK0K(gen^(-1));

    for u0 in u0s do
        xi0 := xi00*hK0K(u0);
        for Phi0 in Phis do
            if &and[ Im(EmbedExtra(xi0 : iota := iota)) gt CC`epscomp : iota in Phi0 ] then
                /* Check and return */
                assert IsTotallyPositive((-xi0^2) @@ hK0K);
                assert ideal< ZZK | xi0 > eq (aa0*aa0bar*DiffK)^(-1);
                return true, [* Phi0, aa0, xi0 *], Kinfo, K0info, bbs, v0s;
            end if;
        end for;
    end for;
end for;
return false, 0, 0, 0, 0, 0;

end function;


intrinsic EnumerationUpToGalois(f::RngUPolElt : exp := 4, prec := 50, precred := 50, prectheta := 50) -> .
{Finds all abelian varieties with primitive CM by the maximal order of the number field defined by f up to Galois conjugation.}

// FIXME: determine bound on exponenent of quotient of group
/* Recover information from calculation for a single case */
vprint CMExp, 1 : "Calculating single solution and field information...";
test0, AV0, Kinfo, K0info, bbs, v0s := SingleClassFromPolynomial(f : exp := exp, prec := prec, precred := precred);
vprint CMExp, 1 : "done calculating single solution and field information.";
if not test0 then
    vprint CMExp, 1 : "Only imprimitive CM types!";
    return [], [];
end if;
Phi, aa0, xi00 := Explode(AV0);
K, ZZK, DiffK, infsK, ClKinfo, UKinfo, invK := Explode(Kinfo); ClK, phiClK := Explode(ClKinfo); UK, phiUK := Explode(UKinfo);
K0, ZZK0, DiffK0, infsK0, ClK0info, UK0info, hK0K := Explode(K0info); ClK0, phiClK0 := Explode(ClK0info); UK0, phiUK0 := Explode(UK0info);
CC := K`CC; CCtheta := ComplexFieldExtra(prectheta);

/* Create totally positive units */
vprint CMExp, 1 : "Determining totally positive units...";
C2 := CyclicGroup(2); C2prod, incs := DirectProduct([ C2 : i in [1..#infsK0] ]);
gensh1 := [ ]; Msign := [ ];
for i in [1..NumberOfGenerators(UK0)] do
    u0 := phiUK0(UK0.i);
    seq := [ ]; Mrow := [ ];
    for inf0 in infsK0 do
        sg := Sign(Evaluate(u0, inf0 : Precision := precred));
        if sg eq 1 then Append(~seq, C2.1^0); else Append(~seq, C2.1^1); end if;
        if sg eq 1 then Append(~Mrow, 0); else Append(~Mrow, 1); end if;
    end for;
    genh1 := &*[ incs[i](seq[i]) : i in [1..#infsK0] ];
    Append(~gensh1, genh1); Append(~Msign, Mrow);
end for;
h1 := hom< UK0 -> C2prod | gensh1 >;
UK1 := Kernel(h1);
B1 := [ UK0 ! UK1.i : i in [1..NumberOfGenerators(UK1)] ];
Msign := Matrix(FiniteField(2), Msign);
vprint CMExp, 1 : "done.";

/* Sanity check */
for u0 in B1 do
    assert &and[ Sign(Evaluate(phiUK0(u0), inf0 : Precision := precred)) eq 1 : inf0 in infsK0 ];
end for;

vprint CMExp, 1 : "Determining abelian varieties...";
AVs := [ ];
for bb in bbs do
    bbbar := ideal< ZZK | [ invK(genbb) : genbb in Generators(bb) ] >;

    /* Descend aa*aabar and find small totally positive generator */
    bbbbbar0 := DescendIdeal(bb*bbbar, ZZK0, hK0K);
    test, eta := IsPrincipalSmall(bbbbbar0, K0info : precred := precred);
    assert test; assert ideal< ZZK0 | eta > eq bbbbbar0;
    v := [ ];
    for inf0 in infsK0 do
        sg := Sign(Evaluate(eta, inf0 : Precision := precred));
        if sg eq 1 then Append(~v, 0); else Append(~v, 1); end if;
    end for;
    v := Vector(FiniteField(2), v);
    u := Solution(Msign, v);
    u0 := UK0 ! [ Integers() ! c : c in Eltseq(u) ];
    u0 := phiUK0(u0);
    u0 := MakeSmallUnit(u0, B1, K0info : precred := precred);
    eta := u0*eta;

    /* At this point eta is a totally positive generator of bb*bbbar */
    /* This implies that the following lines are justified */
    aa := aa0*bb;
    xi0 := xi00*hK0K(eta^(-1));

    /* We can still make things a bit smaller */
    aa, alpha := MakeSmallClassRepresentative(aa, infsK : precred := precred);
    xi0 := xi0*(alpha*invK(alpha))^(-1);
    xi0 := MakeSmallXi(xi0, B1, Kinfo, K0info : precred := precred);

    /* Check and create abelian variety */
    aabar := ideal< ZZK | [ invK(genaa) : genaa in Generators(aa) ] >;
    for v0 in v0s do
        xi := xi0*hK0K(v0);
        assert ideal< ZZK | xi > eq (aa*aabar*DiffK)^(-1);
        assert IsTotallyPositive((-xi^2) @@ hK0K);
        assert &and[ Im(EmbedExtra(xi : iota := iota)) gt CC`epscomp : iota in Phi ];
        Append(~AVs, [* Phi, aa, xi *]);
    end for;
end for;
vprint CMExp, 1 : "done.";

/* Creating and period matrices and counting vanishing even theta constants */
vprint CMExp, 1 : "Calculating thetas...";
schottky_vals := [];
for AV in AVs do
    Phi, aa, xi := Explode(AV);

    // Extract and make into function?
    /* Recover period matrix */
    B := Basis(aa); B := [ K ! b : b in B ];
    P := Matrix(CC, [ [ EmbedExtra(b : iota := iota) : b in B ] : iota in Phi ]);

    /* Recover principal polarization */
    E := Matrix([ [ Trace(xi*invK(x)*y) : x in B ] : y in B ]);
    E := ChangeRing(E, Rationals());
    E := -E;

    /* Convert to big period matrix */
    E := ChangeRing(E, Integers());
    E0, T := FrobeniusFormAlternating(E);
    P := P*ChangeRing(Transpose(T), CC);

    /* Convert to small period matrix and reduce */
    tau := SmallPeriodMatrix(P);
    /*
    vprint CMExp, 1 : "Reducing period matrix...";
    tau := ReduceSmallPeriodMatrix(tau);
    vprint CMExp, 1 : "done with reduction.";
    tausmall := ChangeRing(tau, CCtheta);
    */

    /* Calculate Schottky modular form and return */
    vprint CMExp, 1 : "Reducing period matrix...";
    tau_red := SiegelReduction(tau);
    tau_red := (tau_red + Transpose(tau_red))/2; // hack
    s := SchottkyModularFormFlint(tau_red : prec := prectheta);
    Append(~schottky_vals, [* K, Phi, aa, xi, Abs(s) *]); // field, CM-type, ideal class, polarization, s;
end for;
vprint CMExp, 1 : "done.";
return schottky_vals;

end intrinsic;


intrinsic EnumerationUpToGalois(datum::List : exp := 4, prec := 300, precred := 200, prectheta := 100, Labrande := false, ClassBound := Infinity()) -> .
{Finds the equivalence classes of CM curves for the given datum. ClassBound bounds the class number.}

/* Recover polynomial */
R<x> := PolynomialRing(Rationals());
f := R ! datum[2];

/* Calculate if passes class test */
K := NumberField(f);
if #ClassGroup(K) gt ClassBound then return [ ]; end if;
return EnumerationUpToGalois(f : exp := exp, prec := prec, precred := precred, prectheta := prectheta);

end intrinsic;

// Now compute invariants for abvars that pass the Schottky check

/*
Input: K, Phi, aa, xi, Abs(s)
Output: Invariants of the associated genus 4 curve

1) Make period matrix for abvar with high precision (~1000+)
2) Reconstruct over CC quadric and cubic using reconstructing-g4
4) Compute invariants over CC
4) Recognize invariants as algebraic numbers (bound on the field? Shimura reciprocity?)
    a) First find likely candidate for field K by making minpoly of early invariants
    b) Then try to recognize the rest of the elements over K
5) Reconstruct curve from recognized invariants

*/

function ClassGroupDataf(f : exp := 4, prec := 300, precred := 200)
/* We follow Streng's argument in Proposition 4.4 of "Computing Igusa class polynomials" and optimize along the way  */

/* Create CM field and totally real subfield */
F := RationalsExtra(prec); R := PolynomialRing(F);
K := NumberFieldExtra(R ! f); CC := K`CC;
test, tup, invK := IsCMField(K);
K0, hK0K := Explode(tup);

/* Create rings of integers */
ZZK := Integers(K); ZZK0 := Integers(K0);
/* Create class groups */
ClK, phiClK := ClassGroup(ZZK); ClK0, phiClK0 := NarrowClassGroup(ZZK0);
if (#ClK eq 1) then return 1, 1, 1; end if;
n1 := #ClK / #ClK0;

nClK := NumberOfGenerators(ClK);
infsK := InfinitePlaces(K);
gensClK := [ MakeSmallClassRepresentative(phiClK(ClK.i), infsK : precred := precred) : i in [1..nClK] ];

/* First restrict to classes of ideals for which aa aabar is principal in K0 */
gensh1 := [ ];
for i in [1..nClK] do
    aa := gensClK[i];
    aabar := ideal< ZZK | [ invK(genaa) : genaa in Generators(aa) ] >;
    aaaabar0 := DescendIdeal(aa*aabar, ZZK0, hK0K);
    genh1 := aaaabar0 @@ phiClK0;
    Append(~gensh1, genh1);
end for;
h1 := hom< ClK -> ClK0 | gensh1 >;
G2 := Kernel(h1);
if (#G2 eq 1) then return n1, 1, 1; end if;

gens := Generators(G2);
H2 := sub< G2 | [ 4*gen : gen in gens ] >;
Q := quo< G2 | H2 >;
n2 := #Q; n3 := Exponent(Q);
return n1, n2, n3;

end function;


intrinsic ClassGroupData(datum::List : exp := 4, prec := 300, precred := 200) -> .
{bla}

/* Recover polynomial */
R<x> := PolynomialRing(Rationals());
f := R ! datum[2];

/* Calculate if passes class test */
return ClassGroupDataf(f : exp := exp, prec := prec, precred := precred);

end intrinsic;
