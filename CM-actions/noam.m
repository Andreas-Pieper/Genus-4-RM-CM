P3 := ProjectiveSpace(Rationals(), 3);
_<X,Y,Z,T> := CoordinateRing(P3);
Q1 := X*Z - Y^2;
E1 := -3/20*T^3 + T*(X^2 - 9*X*Z + 9*Y*Z + 27*Z^2) + (X^3 - 15*X^2*Z + 45*X*Z^2 
+ 18*Y*Z^2 + 135*Z^3);
C:= Curve(Scheme(P3, [Q1, E1]));
proj := Projection(C, ProjectiveSpace(Rationals(),2));
R := PolynomialRing(Rationals(),2);
f := Evaluate(DefiningEquation(proj), [R.1, R.2, 1]);
S:= RiemannSurface(f: Precision:= 300);
M, act := EndomorphismAlgebra(BigPeriodMatrix(S));
eigens := NumericalEigenvalues(act[1]);
eigenvecs := [NumericalEigenvectors(act[1], lambda): lambda in eigens];
eigenvecs := [vec/vec[2]: vec in eigenvecs];
diag := Matrix([[BestApproximation(Real(v),100): v in Eltseq(vec)]: vec in eigenvecs]);
diag := Transpose(diag);


