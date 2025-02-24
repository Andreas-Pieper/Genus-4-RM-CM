AttachSpec("~/github/CHIMP/CHIMP.spec");
SetVerbose("EndoCheck",true);
R<X> := PolynomialRing(Rationals());   
f1 := X^4 + 36/7*X^3 + 6*X^2 - 4*X + 1;
F := Evaluate(f1, X^2)*X;
C := HyperellipticCurve(F);
//endos := GeometricEndomorphismRepresentation(C);
K:= NumberFieldExtra(X^2-7);
CK:= ChangeRing(C, K);
P:= CK![1, 8/K.1,1];
//b, corr := Correspondence(CK, CK, endo[2]: P:= P, Q:= P, Al := "Cantor");

E := EllipticCurve([1,0]);
mats, maps := GeometricHomomorphismRepresentation(PeriodMatrix(C), PeriodMatrix(E), QQ);
_, corr := Correspondence(CK, HyperellipticCurve(E), mats[1] : P := P, Al := "Cantor");

H := HyperellipticCurve(X^7 + 7*X^5 + 14*X^3 + 7*X);
L<nu> := NumberFieldExtra(X^2 - 29);
HL := ChangeRing(H,L);
Q := HL![1,L.1,1];
//matsH, mapsH := GeometricHomomorphismRepresentation(PeriodMatrix(C), PeriodMatrix(H), QQ);
matsH, mapsH := GeometricHomomorphismRepresentation(PeriodMatrix(H), PeriodMatrix(C), QQ);
//_, corrH := Correspondence(CK, H, matsH[1] : P := P, Al := "Cantor");
//_, corrH := Correspondence(CK, H, matsH[2] : P := P, Al := "Cantor");
_, corrH := Correspondence(HL, CK, matsH[1] : P := Q, Q := P, Al := "Cantor");
