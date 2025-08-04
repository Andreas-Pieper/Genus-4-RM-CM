/***
 *  Reflex field functionality
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */

SetClassGroupBounds("GRH");
declare attributes FldNum : Cl;


intrinsic NumericalReflexField(Phi::., K::.) -> .
{Recovers the reflex field numerically (not used for now).}
prec := Precision(K`CC);
B := Basis(K);
/* Use exact definition of the reflex field as a subfield of CC */
BCC := [ &+[ EmbedExtra(b : iota := iota) : iota in Phi ] : b in B ];
return NumberFieldExtra(BCC, RationalsExtra(prec));
end intrinsic;


intrinsic LiftedCMType(PhiK::., K::.) -> .
{Lifts the CM type PhiK to the Galois closure of K.}

/* Recover Galois closure from PhiK */
assert assigned K`embL; L, hKL, invL, phiG := Explode(K`embL);
auts := [ phiG(g) : g in Domain(phiG) ];
rts := [ aut(L.1) : aut in auts ];
PhiL := [ ];
/* Take all automorphisms that act appropriately */
for rt in rts do
    h := hom< L -> L | rt >;
    if h(hKL(K.1)) in PhiK then
        Append(~PhiL, rt);
    end if;
end for;
return Set(PhiL);

end intrinsic;


intrinsic ReflexField(PhiK::., K::.) -> .
{Returns the reflex field of CM type PhiK along with an embedding in the Galois closure and the reflex type.}

/* Lift CM type */
assert assigned K`embL; L, hKL, invL, phiG := Explode(K`embL);
PhiL := LiftedCMType(PhiK, K);
G := Domain(phiG);
S := [ g : g in G | phiG(g)(L.1) in PhiL ];

/* Define reflex field as fixed field */
Hr := [ g : g in G | &and[ s*g in S : s in S ] ];
Hr := sub< G | Hr >;
Kr := FixedField(L, [ phiG(h) : h in Hr ]);
hKrL := hom< Kr -> L | L ! Kr.1 >;

/* Improve representation */
Kr0, hKrKr0 := Polredabs(Kr);
hKr0Kr := Inverse(hKrKr0);
hKr0L := hKr0Kr*hKrL;
embL := [* L, hKr0L, invL, phiG *];
Kr0`embL := embL;

/* Determine reflex type */
Sinv := [ s^(-1) : s in S ];
PhiKr0 := { phiG(g)(hKr0L(Kr0.1)) : g in Sinv };
assert #PhiKr0 eq (Degree(Kr0) div 2);
return Kr0, PhiKr0;

end intrinsic;


intrinsic TypeNorm(PhiK::., K::., Kr::.) -> .
{Returns the type norm as a function.}

assert assigned K`embL; L, hKL, invL, phiG := Explode(K`embL);
assert assigned Kr`embL; L, hKrL, invL, phiG := Explode(Kr`embL);
/* Interpret CM type as maps */
hs := [ hom< K -> L | rt > : rt in PhiK ];

/* Take function that returns product */
function Nm(a)
b := &*[ h(a) : h in hs ];
return (b @@ hKrL);
end function;
return Nm;

end intrinsic;


intrinsic ReflexTypeNorm(PhiKr::., Kr::., K::.) -> .
{Returns the reflex type norm as a function.}

assert assigned K`embL; L, hKL, invL, phiG := Explode(K`embL);
assert assigned Kr`embL; L, hKrL, invL, phiG := Explode(Kr`embL);
/* Interpret CM type as maps */
hs := [ hom< Kr -> L | rt > : rt in PhiKr ];

/* Take function that returns product */
function RefNm(b)
a := &*[ h(b) : h in hs ];
return (a @@ hKL);
end function;
return RefNm;

end intrinsic;


function DescendIdeal(J, ZZK, h)
/* Descends ideal J in ZZL to ZZK under inclusion given by h */

ZZL := Order(J); BK := Basis(ZZK); BJ := Basis(J);
/* Take integral bases of ZZK as well as J */
MK := Matrix(Integers(), [ Eltseq(ZZL ! h(b)) : b in BK ]);
MJ := Matrix(Integers(), [ Eltseq(b) : b in BJ ]);
LatK := Lattice(MK);
LatJ := Lattice(MJ);

/* Intersect and solve relevant equations */
LatI := LatK meet LatJ;
MI := Matrix(Basis(LatI));
sol := Solution(MK, MI);
rows := [ Eltseq(row) : row in Rows(sol) ];
I := ideal< ZZK | [ &+[ row[i]*BK[i] : i in [1..#BK] ] : row in rows ] >;

/* Check and return */
IL := ideal< ZZL | [ h(gen) : gen in Generators(I) ] >;
assert J eq IL;
return I;

end function;


intrinsic TypeNormIdeals(PhiK::., K::., Kr::.) -> .
{Returns the type norm on ideals as a function.}

assert assigned K`embL; L, hKL, invL, phiG := Explode(K`embL);
assert assigned Kr`embL; L, hKrL, invL, phiG := Explode(Kr`embL);
ZZK := Integers(K); ZZKr := Integers(Kr); ZZL := Integers(L);
/* Interpret CM type as maps */
hs := [ hom< K -> L | rt > : rt in PhiK ];

function Nm(aa)
/* Take images of the ideal and descend the product */
gens := Generators(aa);
aas := [ ideal< ZZL | [ h(gen) : gen in gens ] > : h in hs ];
bb := &*aas;
return DescendIdeal(bb, ZZKr, hKrL);
end function;
return Nm, ZZK, ZZKr, ZZL;

end intrinsic;


intrinsic ReflexTypeNormIdeals(PhiKr::., Kr::., K::.) -> .
{Returns the type norm on ideals as a function.}

assert assigned K`embL; L, hKL, invL, phiG := Explode(K`embL);
assert assigned Kr`embL; L, hKrL, invL, phiG := Explode(Kr`embL);
ZZK := Integers(K); ZZKr := Integers(Kr); ZZL := Integers(L);
/* Interpret CM type as maps */
hs := [ hom< Kr -> L | rt > : rt in PhiKr ];

function RefNm(bb)
/* Take images of the ideal and descend the product */
gens := Generators(bb);
bbs := [ ideal< ZZL | [ h(gen) : gen in gens ] > : h in hs ];
aa := &*bbs;
return DescendIdeal(aa, ZZK, hKL);
end function;
return RefNm, ZZK, ZZKr, ZZL;

end intrinsic;


intrinsic TypeNormIdealClasses(PhiK::., K::., Kr::.) -> .
{Returns the type norm on ideal classes as a function.}

assert assigned K`embL; L, hKL, invL, phiG := Explode(K`embL);
assert assigned Kr`embL; L, hKrL, invL, phiG := Explode(Kr`embL);
ZZK := Integers(K); ZZKr := Integers(Kr); ZZL := Integers(L);
/* Interpret CM type as maps */
hs := [ hom< K -> L | rt > : rt in PhiK ];

/* Take images of ideals in the class group and descend the product */
ClK, phiClK := ClassGroup(K); ClKr, phiClKr := ClassGroup(Kr);
function Nm(cl)
aa := phiClK(cl);
gens := Generators(aa);
aas := [ ideal< ZZL | [ h(gen) : gen in gens ] > : h in hs ];
bb := &*aas;
bb := DescendIdeal(bb, ZZKr, hKrL);
return bb @@ phiClKr;
end function;
return Nm, [* ClK, phiClK *], [* ClKr, phiClKr *];

end intrinsic;


intrinsic ReflexTypeNormIdealClasses(PhiKr::., Kr::., K::.) -> .
{Returns the reflex type norm on ideal classes as a function.}

assert assigned K`embL; L, hKL, invL, phiG := Explode(K`embL);
assert assigned Kr`embL; L, hKrL, invL, phiG := Explode(Kr`embL);
ZZK := Integers(K); ZZKr := Integers(Kr); ZZL := Integers(L);
/* Interpret CM type as maps */
hs := [ hom< Kr -> L | rt > : rt in PhiKr ];

/* Take images of ideals in the class group and descend the product */
ClKr, phiClKr := ClassGroup(Kr); ClK, phiClK := ClassGroup(K);
function RefNm(cl)
bb := phiClKr(cl);
gens := Generators(bb);
bbs := [ ideal< ZZL | [ h(gen) : gen in gens ] > : h in hs ];
aa := &*bbs;
aa := DescendIdeal(aa, ZZK, hKL);
return aa @@ phiClK;
end function;
return RefNm, [* ClKr, phiClKr *], [* ClK, phiClK *];

end intrinsic;


intrinsic ClassGroupQuotient(PhiK::., K::.) -> .
{Returns the quotient of the totally positive part of the class group by the image of the reflex type norm.}

/* Construct all relevant objects */
assert assigned K`embL; L, hKL, invL, phiG := Explode(K`embL);
Kr, PhiKr := ReflexField(PhiK, K);
assert assigned Kr`embL; L, hKrL, invL, phiG := Explode(Kr`embL);
ZZK := Integers(K); ZZKr := Integers(Kr); ZZL := Integers(L);
ClKr, phiClKr := ClassGroup(Kr); ClK, phiClK := ClassGroup(K);
K0, tup, invK := IsCMField(K); K0, hK0K := Explode(tup);
ZZK0 := Integers(K0); ClK0, phiClK0 := ClassGroup(K0);

/* First restrict to classes of ideals for which aa aabar is principal in K */
gensh := [ ]; gensClK := Generators(ClK);
for i in [1..#gensClK] do
    genClK := ClK.i;
    aa := phiClK(genClK);
    aabar := ideal< ZZK | [ invK(genaa) : genaa in Generators(aa) ] >;
    genh := (aa*aabar) @@ phiClK;
    Append(~gensh, genh);
end for;
h := hom< ClK -> ClK | gensh >;
ClK1 := Kernel(h);

/* Then restrict to classes of ideals for which aa aabar is principal in K0 */
gensh := [ ]; gensClK1 := Generators(ClK1);
for i in [1..#gensClK1] do
    genClK := ClK1.i;
    aa := phiClK(genClK);
    ZZK := Order(aa);
    aabar := ideal< ZZK | [ invK(genaa) : genaa in Generators(aa) ] >;
    bb := aa*aabar;
    bb0 := DescendIdeal(bb, ZZK0, hK0K);
    genh := bb0 @@ phiClK0;
    Append(~gensh, genh);
end for;
h := hom< ClK1 -> ClK0 | gensh >;
ClK2 := Kernel(h);

/* Finally restrict classes of ideals for which aa aabar is principal in K0 and generated by a totally positive element */
/* To this end we first determine the signs to which units give rise and take the relevant quotient */
infs := InfinitePlaces(K0);
C2 := CyclicGroup(2); C2prod, incs := DirectProduct([ C2 : i in [1..#infs] ]);
U, phiU := UnitGroup(K0 : GRH := true);
gensC2prodsub := [ ];
for genU in Generators(U) do
    u := phiU(genU);
    seq := [ ];
    for inf in infs do
        sg := Sign(Evaluate(u, inf));
        if sg eq 1 then
            Append(~seq, C2.1^0);
        else
            Append(~seq, C2.1^1);
        end if;
    end for;
    gen := &*[ incs[i](seq[i]) : i in [1..#infs] ];
    Append(~gensC2prodsub, gen);
end for;
Q, pQ := quo< C2prod | gensC2prodsub >;

/* Now take images of generators and form kernel */
gensh := [ ]; gensClK := Generators(ClK2);
for i in [1..#gensClK] do
    genClK := ClK2.i;
    aa := phiClK(genClK);
    ZZK := Order(aa);
    aabar := ideal< ZZK | [ invK(genaa) : genaa in Generators(aa) ] >;
    bb := aa*aabar;
    bb0 := DescendIdeal(bb, ZZK0, hK0K);

    test, gen := IsPrincipal(bb0);
    assert test;
    seq := [ ];
    for inf in infs do
        sg := Sign(Evaluate(gen, inf));
        if sg eq 1 then
            Append(~seq, C2.1^0);
        else
            Append(~seq, C2.1^1);
        end if;
    end for;
    genh := &*[ incs[i](seq[i]) : i in [1..#infs] ];
    Append(~gensh, genh);
end for;
h := hom< ClK2 -> Q | [ pQ(genh) : genh in gensh ] >;
ClK3 := Kernel(h);

/* Interpret CM type as maps */
hs := [ hom< Kr -> L | rt > : rt in PhiKr ];
function RefNm(cl)
bb := phiClKr(cl);
gens := Generators(bb);
bbs := [ ideal< ZZL | [ h(gen) : gen in gens ] > : h in hs ];
aa := &*bbs;
aa := DescendIdeal(aa, ZZK, hKL);
return aa @@ phiClK;
end function;

/* Reflex norm has image in smaller group; take appropriate quotient */
ims := [ ClK3 ! RefNm(clr) : clr in Generators(ClKr) ];
ImRef := sub< ClK3 | ims >;
Q := quo< ClK3 | ImRef >;
return Q, ClK3, ImRef, [* ClKr, phiClKr *], [* ClK, phiClK *];

end intrinsic;
