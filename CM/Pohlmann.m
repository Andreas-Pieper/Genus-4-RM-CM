/*R<x> := PolynomialRing(Rationals());
f:= x^6 - 3*x^5 + 9*x^4 - 13*x^3 + 14*x^2 - 8*x + 2;
K := ext<Rationals()| f>;
L := SplittingField(K);
G := AutomorphismGroup(Rationals(), K);*/
file := Open("Noot-cm-fields-jacobians.txt", "r");
string := Read(file);
split := Split(string, "\n");
coeffs := [eval(Split(str, "|")[2]): str in split];
R<x> := PolynomialRing(Rationals());
polys := ChangeUniverse(coeffs, R);
fields := [NumberField(f): f in polys];
K3 := CyclotomicField(3);
K5 := CyclotomicField(5);
fieldsfil := [K: K in fields| not (IsSubfield(K3, K) or IsSubfield(K5, K))];

//8.0.391881616.1|[16,-24,8,8,-10,4,2,-3,1]|1|9.8001694745000484506503377348537074841863176962905910006835581192306551218425209E-80

