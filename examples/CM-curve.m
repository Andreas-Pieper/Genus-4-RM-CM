CC<I> := ComplexField(100);
CC4<X,Y,Z,T> := PolynomialRing(CC, 4);

AttachSpec("Genus-4/magma/spec");
import "Genus-4/magma/Invariants_genus_4.m" : ChangeOfBasis, CubicNewBasis, InvariantsGenus4CurvesRank3;
import "Genus-4-RM-CM/PotentialCMCurves.m" : Eqs;

// We take an example for which the invariants seem rational 
Q, E := Explode(Eqs[5]);
D, P := DiagonalForm(Q);
Q := CC4!(D-MonomialCoefficient(D, (Parent(D).4)^2)*(Parent(D).4)^2);
E := CC4!(ChangeOfBasis(E, P));

R1 := RealField(30);
inv, wgt := InvariantsGenus4Curves(Q, E);
inv := WPSMultiply(wgt, inv, (1/inv[1])^(1/6));
inv_r := [MinimalPolynomial(R1!Real(i), 2) : i in inv];
// not precise enough to recognize the fifth invariant as rational, we try to compute directly the Igusa invariants. 

// construction of the degree 6 form 
f0 := CubicNewBasis(Q, E);
R0 := BaseRing(Parent(f0));
R<s, t, w> := PolynomialRing(R0, [1,1,2]);
f_weighted := Evaluate(f0, [s^2, s*t, t^2, w]);

alpha := MonomialCoefficient(f_weighted, w^3);
f_weighted /:= alpha;  

S<[x]> := PolynomialRing(R0, 2);
f_weighted := Evaluate(f_weighted, [s, t, w-ExactQuotient(Terms(f_weighted, w)[3], 3*w^2)]);

f := S!Evaluate(f_weighted, [x[1], x[2], 0]); // degree 6 form
Inv, wgt2 := IgusaInvariants(f);
Inv := WPSMultiply(wgt2, Inv, (1/Inv[1])^(1/2));
inv_r := [BestApproximation(R1!Real(i), 10^20) : i in Inv]; // we find that the invariant of the binary sextic are rational!

HyperellipticPolynomials(HyperellipticCurveFromIgusaInvariants(inv_r)); // a rational model of that sextic form

R0<[a]> := PolynomialRing(Rationals(), 5);
S<[x]> := PolynomialRing(R0, 2);
f_rec := (x[1]^6 - 15*x[1]^4*x[2]^2 + 45*x[1]^2*x[2]^4 + 18*x[1]*x[2]^5 + 135*x[2]^6);
v_rec := &+[a[i+1]*x[1]^i*x[2]^(4-i) : i in [0..4]]; // indeterminate binary quartic form
Inv2 := InvariantsGenus4CurvesRank3(f_rec, v_rec);

I := Ideal([Inv2[7]-3375/2, Inv2[11]-23625/8, Inv2[16]+2541375/8, Inv2[19]+50625/4, Inv2[20]+2641500/11, Inv2[27]+102254400/11, Inv2[30]+830250, Inv2[33]-42342750]); // we recognize some invariants as rationals
GroebnerBasis(I); // we know the coefficients of v also

f_rec := -20/3*(x[1]^6 - 15*x[1]^4*x[2]^2 + 45*x[1]^2*x[2]^4 + 18*x[1]*x[2]^5 + 135*x[2]^6);
v_rec := -20/3*(x[1]^4 - 9*x[1]^2*x[2]^2 + 9*x[1]*x[2]^3 + 27*x[2]^4);
Inv2 := InvariantsGenus4CurvesRank3(f_rec, v_rec);
Inv2 := WPSNormalize(wgt, Inv2);
inv1 := WPSNormalize(wgt, inv);

/* check that we have the same invariants
for i in [1..60] do
        i;
        CC!Inv2[i];
        inv1[i];
        "";
end for;
*/

P3 := ProjectiveSpace(Rationals(), 3);
_<X,Y,Z,T> := CoordinateRing(P3);
Q1 := X*Z - Y^2;
E1 := -3/20*T^3 + T*(X^2 - 9*X*Z + 9*Y*Z + 27*Z^2) + (X^3 - 15*X^2*Z + 45*X*Z^2 + 18*Y*Z^2 + 135*Z^3);

inv_rec := InvariantsGenus4Curves(Q1, E1 : normalize := true);

/* check that we have the same invariants
for i in [1..60] do
        i;
        CC!inv_rec[i];
        inv1[i];
        "";
end for;
*/

// the original curve 
P2 := ProjectiveSpace(Rationals(), 2);
_<x,y,z> := CoordinateRing(P2);
f := y^3*z^2 - x^5 - z^5;
D := Curve(P2, f);
D := Image(CanonicalEmbedding(D));
E2, Q2 := Explode(Basis(Ideal(D)));

inv_orig := InvariantsGenus4Curves(Q2, E2 : normalize := true);
inv_orig eq inv_rec; // not isomorphic to y^2 = x^5 + 1


// now we check that the reconstructed curve is irreducible, reduced, smooth
C := Curve(P3, [Q1, E1]);
"Smooth:", IsNonSingular(C);
"Reduced:", IsReduced(C);
"Irreducible:", IsIrreducible(C);
