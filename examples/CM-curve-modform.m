AttachSpec("~/github/CHIMP/CHIMP.spec");
P3 := ProjectiveSpace(RationalsExtra(), 3);
_<X,Y,Z,T> := CoordinateRing(P3);
Q1 := X*Z - Y^2;
E1 := -3/20*T^3 + T*(X^2 - 9*X*Z + 9*Y*Z + 27*Z^2) + (X^3 - 15*X^2*Z + 45*X*Z^2 + 18*Y*Z^2 + 135*Z^3);
C := Curve(P3, [Q1, E1]);
// compute aps
aps := AssociativeArray();
for p in PrimesInInterval(7,100) do
  Cp := Curve(Reduction(C,p));
  //print LPolynomial(Cp);
  aps[p] := p + 1 - #Points(Cp, GF(p));
  print aps[p];
end for;
for p in Sort(Setseq(Keys(aps))) do
  print <p, aps[p]>;
end for;
