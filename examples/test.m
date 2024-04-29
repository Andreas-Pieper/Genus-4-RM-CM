R<x,y>:= PolynomialRing(Rationals(),2);
F := -3/20*y^3 + y* (x^4 - 9*x^2 + 9*x + 27) + (x^6 - 15*x^4 + 45*x^2 + 18*x + 
135);
C:= Curve(AffineSpace(R), F);
candidates:= [];
for nu1, nu2, nu3 in [1..14]  do
F2 := y^3+x^5+2^nu1*3^nu2*5^nu3;
C2:= Curve(AffineSpace(R), F2);

for p in PrimesInInterval(6,1000) do

   Cp := Curve(Reduction(C,p));
   C2p := Curve(Reduction(C2,p));
   S<t> := PolynomialRing(GF(p));
   inf := -3/20*t^3 + t + 1;
   num := #Points(Cp) + #Roots(inf);
   num2 := #Points(C2p) + 1;
   if num ne num2 then
	   continue nu3;
   end if;
end for;
print "hurray", nu1, nu2, nu3;
Append(~candidates, [nu1, nu2, nu3]);
end for;
//;
F2 := y^3+x^5+58320000000000;
C2:= Curve(AffineSpace(R), F2);

for p in PrimesInInterval(6,1000) do
   Cp := Curve(Reduction(C,p));
   C2p := Curve(Reduction(C2,p));
   S<t> := PolynomialRing(GF(p));
   inf := -3/20*t^3 + t + 1;
   num := #Points(Cp) + #Roots(inf);
   num2 := #Points(C2p) + 1;
   print p, num eq num2;
end for;


