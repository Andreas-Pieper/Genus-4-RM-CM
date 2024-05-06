/*
  Here's some information about the apparent CM point on the g=4 Schottky locus.

  The Jacobian of the curve  y^3 = x^5 + 1 has complex multiplication by the ring of integers in the 15th cyclotomic field, which is a CM field with real subfield K.  This seems to match the modular form spaces 1125.2.a.g and 1125.2.a.m, at least over C.(*)

  To get tau_1, ..., tau_4, start from the four complex numbers
    ```code
    -r^3 - 3*r^2 + 2*r + 13/2 - (4*r^3 + 2*r^2 - 13*r - 7/2)*sqrt(-3)
    ```
    with r ranging over the roots of `x^4 - x^3 - 4*x^2 + 4*x + 1 = 0`; all but one have positive imaginary part, and taking the complex conjugate of the fourth gives the desired point.

    This point is of course on the Schottky locus (indeed both E4 and F8 vanish at that point so certainly E4^2 - F8 does).  However, if we replace each $$\tau_k$$ by $$2\tau_k$$, we get an isogenous 4-fold that looks like it's also on the Schottky locus: the sums $$2^{-4} \sum \theta^{8n}$$ for n=1,2,...,6 scale to
          ```
        50625, 2562890625, 229673318429025, 25483874882915700225,
        3095290713834020200752225, 392715761148530200969393042305
        ```
        and the first two are 15^4 and 15^8, so in particular the F_8/E_4^2 ratio is 1.

        Can you find a model over Q for this genus-4 curve?

*/

// copied from 17T7 project

AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");
Attach("~/github/Genus-4-RM-CM/maximal-isotropic.m");
AttachSpec("~/github/Genus-4/magma/spec");
AttachSpec("~/github/Reconstruction/magma/spec");
// TODO: Fix Flint theta port
//Attach("~/github/Genus-4-RM-CM/FlintWrapper.m");
//Attach("~/github/Genus-4-RM-CM/schottky-fast-theta.m");
// Given a totally real field F of degree g and a g-tuple of points in the upper half-plane, return the corresponding big period matrix
function ModuliToBigPeriodMatrixNoam(F, points)
//intrinsic ModuliToBigPeriodMatrixNoam(F, points) -> AlgMatElt
//{ Modified version of the ModuliToBigPeriodMatrix. }
    prec := Min([Precision(Parent(elt)) : elt in points]);
    CC := ComplexFieldExtra(prec);
    //assert &and[Abs(Re(p)) lt CC`epscomp : p in points];
    OF := Integers(F);
	B := Basis(OF);
    g := Degree(F);
    betas := [[CC | Evaluate(B[i], pl : Precision := prec+10) : pl in InfinitePlaces(F)] : i in [1..g]];
    Pi1 := Transpose(Matrix(CC, betas));
    Pi2 := DiagonalMatrix(points)*(Transpose(Pi1)^-1);
    return HorizontalJoin(Pi2, Pi1);
end function;

// construct original curve
// TODO add reconstruction-g4 as a git submodule
SetVerbose("Reconstruction",true);
SetVerbose("Theta",true);
SetDebugOnError(true);
prec := 300;
g := 4;
CC<I> := ComplexFieldExtra(prec);
R<x,y> := PolynomialRing(QQ,2);
/*
  f := y^3-(x^5+1);
  //C := Curve(Spec(R), y^3-(x^5+1));
  S := RiemannSurface(f : Precision := prec);
  Pi := BigPeriodMatrix(S);
  Pi1 := Submatrix(Pi,1,1,g,g);
  Pi2 := Submatrix(Pi,1,g+1,g,g);
  _<t> := PolynomialRing(QQ);
  roots := [el[1] : el in Roots(t^4 - t^3 - 4*t^2 + 4*t + 1, CC)];
  taus := [-r^3 - 3*r^2 + 2*r + 13/2 - (4*r^3 + 2*r^2 - 13*r - 7/2)*Sqrt(CC!-3) : r in roots];
  taus[4] := ComplexConjugate(taus[4]);
  //taus := [2*el : el in taus]; // 2-isogeny
  F<nu> := NumberFieldExtra(t^4 - t^3 - 4*t^2 + 4*t + 1);
  Pi_taus := ModuliToBigPeriodMatrixNoam(F,taus);
  Pi_taus := Submatrix(Pi_taus,1,g+1,g,g);
  //Pi_taus_small := -SmallPeriodMatrix(Pi_taus);
  Pi_taus_small := Pi_taus2^-1*Pi_taus1;
  //Pi_taus_small := SiegelReduction(Pi_taus_small);
  //printf "Schottky modular form = %o\n", SchottkyModularForm(Pi_taus_small);
  thetas := ComputeThetas(Pi_taus_small);
  pt := [2^-g*&+[t^(8*n) : t in thetas] : n in [1..6]];
  pt := WPSMultiply([1..6], WPSNormalize([1..6], pt), 50625);
  //print "attempting to reconstruct curve";
  //ReconstructCurveG4(Pi_taus_small);
  //
*/
f := y^3 - (x^5+1);
S := RiemannSurface(f : Precision := prec);
Pi_small := SmallPeriodMatrix(S);
Pi_big := BigPeriodMatrix(S);
Pi1 := Submatrix(Pi_big,1,1,g,g);
Pi2 := Submatrix(Pi_big,1,g+1,g,g);
//Pi_big_swap := HorizontalJoin(Pi2,Pi1);
P := HorizontalJoin(Pi1,Pi2);
torsion := maximal_isotropic(4,2);
jacobians := [];
for i->V in torsion do
  Q := QFromPVFor4(P, V);
  Q1 := Submatrix(Q,1,1,g,g);
  Q2 := Submatrix(Q,1,g+1,g,g);
  Q_sm := Q1^-1*Q2;
  s := SchottkyModularForm(Q_sm : prec := 30);
  printf "subgroup %o had size %o\n", i, Abs(s);
  if Abs(s) lt 10^-20 then
    Append(~jacobians, [* i, V *]);
    break;
  end if;
end for;

curves_CC := [];
invs_CC := [];
for pair in jacobians do
  i, V := Explode(pair);
  assert torsion[i] eq V;
  Q := QFromPVFor4(P, V);
  Q1 := Submatrix(Q,1,1,g,g);
  Q2 := Submatrix(Q,1,g+1,g,g);
  Q_sm := Q1^-1*Q2;
  Q_red, M := SiegelReduction(Q_sm);
  C_CC := ReconstructCurveG4(Q_red);
  Append(~curves_CC, C_CC);
  invs, wts := InvariantsGenus4Curves(C_CC[1], C_CC[2] : normalize := true);
  Append(~invs_CC, [* invs, wts *]);
end for;

quadric, cubic := Explode(curves_CC[1]);
invs_old, wts := InvariantsGenus4Curves(quadric, cubic);
lambda := (1/invs_old[1])^(1/wts[1]);
invs := WPSMultiply(wts, invs_old, lambda);
m, err :=  MinimalPolynomial(invs[2], 2);
Polredabs(m);
K<nu> := NumberField(Polredabs(m));
OK := Integers(K);
basis := Basis(OK);
v := InfinitePlaces(K)[1];
basis_CC := [CC!Evaluate(el, v : Precision := prec) : el in basis];

ps := [];
es := [];
invs_K := [];
for i->inv in invs do
  printf "ps = %o, es = %o\n", ps, es;
  scale := (&*[Integers() | ps[k]^(Round(es[k] * wts[i])) : k in [1..#ps]]);
  inv_scaled := inv * scale;
  rel := IntegerRelation(basis_CC cat [inv_scaled], 10^20);
  printf "Integer relation found %o\n", rel;
  d := rel[#rel];
  print "Factoring";
  fact := Factorization(d);
  psd := [[p, j]: j->p in ps| not IsDivisibleBy(d, p)];
  for pind in psd do
      p, j := Explode(pind);
      exp := Min([Valuation(rel[1], p), Valuation(rel[2], p)]);
      e_new := (Round(es[j] * wts[i]) - exp)/ wts[i];
      es[j] := Max([es[j], e_new]);
  end for;
  for pair in fact do
    p, e := Explode(pair);
    e_new := e/wts[i];
    if p in ps then
      j := Index(ps, p);
      es[j] := Round(es[j] * wts[i]) / wts[i] + e_new;
    else
      Append(~ps,p);
      Append(~es, e_new);
    end if;
  end for;
  inv_K := (-1/d)*&+[rel[i]*basis[i] : i in [1..#basis]];
  inv_K /:= scale;
  printf "Invariant recognized: %o\n", inv_K;
  Append(~invs_K, inv_K);
end for;

/*
inv_K;
Evaluate(inv_K, v);
invs_new[2];
*/

invs := [];
for pair in invs_CC do
  inv_CC, wt := Explode(pair);
  inv := [b where _, b := RationalReconstruction(el) : el in inv_CC];
  inv := WPSMinimize(wt, inv);
  Append(~invs, [* inv, wt *]);
end for;

//ReconstructionGenus4(invs[1][1]);
// Not implemented yet :(

curves := [];

import "~/github/Reconstruction/magma/reconstruction_genus4.m": ReconstructionGenus4Rank3;


/*
V := torsion[15];
Q := QFromPVFor4(P, V);
Q1 := Submatrix(Q,1,1,g,g);
Q2 := Submatrix(Q,1,g+1,g,g);
Q_sm := Q1^-1*Q2;
s := SchottkyModularForm(Q_sm : prec := 80);
s;
*/



//InvariantsGenus4Curves(quadric,cubic : normalize := true);

AttachSpec("~/github/CHIMP/CHIMP.spec");
P3 := ProjectiveSpace(RationalsExtra(), 3);
_<X,Y,Z,T> := CoordinateRing(P3);
Q1 := X*Z - Y^2;
E1 := -3/20*T^3 + T*(X^2 - 9*X*Z + 9*Y*Z + 27*Z^2) + (X^3 - 15*X^2*Z + 45*X*Z^2 + 18*Y*Z^2 + 135*Z^3);
C := Curve(P3, [Q1, E1]);

// computing endomorphism ring
// non-maximal order of QQ(zeta15)
/*
E := GeometricEndomorphismRepresentation(C);
ends := [E[j][2] : j in [1..#E]];
ends;
for el in ends do
  print MinimalPolynomial(el);
end for;
[ChangeRing(el, ZZ) : el in ends];
ends_ZZ := $1;
M := Matrix([Eltseq(el) : el in ends_ZZ]);
D, P, Q := SmithForm(M);
D;
P;
Q;
for el in ends do
for el in ends_ZZ do
for el in ends_ZZ do
MinimalPolynomial(el);
end for;
for el in ends_ZZ do
Polredabs(MinimalPolynomial(el));
end for;
$1[1];
MinimalPolynomial(ends_ZZ[#ends_ZZ]);
m := $1;
K<nu> := NumberField(m);
K;
Factorization(ChangeRing(m,K));

MinimalPolynomial(ends_ZZ[#ends_ZZ]);
m := $1;
K<nu> := NumberField(m);
K;
Factorization(ChangeRing(m,K));
MinimalPolynomial(ends_ZZ[2]);
Roots($1,K);
Discriminant(m);
Factorization($1);
MinimalPolynomial(ends_ZZ[#ends_ZZ]);
m := $1;
Discriminant(m);
Factorization($1);
O := EquationOrder(m);
UnitGroup(O);
ClassGroup(O);
PicardGroup(O);
GorensteinIndex;
IsGorenstein;
D := Different(O);
Conductor(O);
c := $1;
Norm(c);
IsPrincipal(c);
*/
