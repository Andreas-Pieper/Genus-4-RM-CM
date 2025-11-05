QQ := Rationals();
t := 3;
  S<X,Y,Z,T> := PolynomialRing(QQ, 4);
  E := 3*T^3 + t*(t - 1)*(3*T*(5*X^2 + 6*X*Y + 2*t*Y*Z + 3*t*Z^2) + (-2*t + 9)*X^3 + 22*t*X^2*Y + 21*t*X^2*Z + (-14*t^2 + 18*t)*X*Y*Z + t^2*X*Z^2 + 6*t^2*Y*Z^2 + (-3*t^3 + 6*t^2)*Z^3);
R<x,y> := PolynomialRing(QQ, 2);
F := Evaluate(E, [y^2,y,1, x]);
Fy := Derivative(F, y);
Fx := Derivative(F, x);
wahl := [Fy,Fx,2*y*Fx,-y*Fy-x*Fx+3*F,-y^2*Fy-2*x*y*Fx+108/18*y*F,-y^2*Fx];
mons := Sort(Setseq(Seqset(&cat[Monomials(c): c in wahl])));
mat := Matrix([[MonomialCoefficient(w, m): w in wahl]: m in mons]);
R<v0,v1,v2,v3> := PolynomialRing(Rationals(), 4);
mm := Matrix([[v0,v1,v2,0,0,0,0,v3,0,0,0,0,0,0,0], [0,v0,v1,v2,0,0,0,0,v3,0,0,0,0,0,0], [0,0,v0,v1,v2,0,0,0,0,v3,0,0,0,0,0], [0,0,0,v0,v1,v2,0,0,0,0,v3,0,0,0,0], [0,0,0,0,v0,v1,v2,0,0,0,0,v3,0,0,0], [0,0,0,0,0,0,0,v0,v1,v2,0,0,v3,0,0], [0,0,0,0,0,0,0,0,v0,v1,v2,0,0,v3,0], 
[0,0,0,0,0,0,0,0,0,v0,v1,v2,0,0,v3],[0,0,0,0,0,0,0,0,0,0,0,0,v0,v1,v2]]);
mm[9] := Vector([mm[9,i]- v3/3*MonomialCoefficient(F, mons[i]): i in [1..15]]);
RR := PolynomialRing(R);
eva := Evaluate(F, [-v0-v1*RR.1-v2*RR.1^2 , RR.1]);
igusa := IgusaInvariants(eva);
clebsch  := ClebschInvariants(eva);

