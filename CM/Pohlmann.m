R<x> := PolynomialRing(Rationals());
f:= x^6 - 3*x^5 + 9*x^4 - 13*x^3 + 14*x^2 - 8*x + 2;
K := ext<Rationals()| f>;
L := SplittingField(K);
G := AutomorphismGroup(Rationals(), K);

