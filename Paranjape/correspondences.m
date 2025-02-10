R<X> := PolynomialRing(Rationals());   
f1 := X^4 + 36/7*X^3 + 6*X^2 - 4*X + 1;
F := Evaluate(f1, X^2)*X;
C := HyperellipticCurve(F);
endos := GeometricEndomorphismRepresentation(C);
K:= NumberFieldExtra(X^2-7);
CK:= ChangeRing(C, K);
P:= CK![1, 8/K.1,1];
b, corr := Correspondence(CK, CK, endo[2]: P:= P, Q:= P, Al := "Cantor");

