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
prec := 100;
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
  end if;
end for;
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
