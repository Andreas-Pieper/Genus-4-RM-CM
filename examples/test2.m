K<nu> := QuadraticField(-15);
P3 := ProjectiveSpace(K,3);
R4<X,Y,Z,T> := CoordinateRing(P3);
Q1 := (64202067720000*nu + 1777004013960000)*X^2 + (39914835912000*nu +
    951312813192000)*X*Y + (12224250181500*nu + 254418115868700)*X*Z +
    (6112125090750*nu + 127209057934350)*Y^2 + (3682347994080*nu +
    68024929079520)*Y*Z + (547806864823*nu + 9087980262431)*Z^2;
Q1 := (64202067720000*nu + 1777004013960000)*X^2 + (39914835912000*nu +
    951312813192000)*X*Y + (12224250181500*nu + 254418115868700)*X*Z +
    (6112125090750*nu + 127209057934350)*Y^2 + (3682347994080*nu +
    68024929079520)*Y*Z + (547806864823*nu + 9087980262431)*Z^2;
E1 := (95892897844943325000000*nu + 1670256537361323093000000)*X^3 +
    (85248803193107948100000*nu + 1338410452870427220900000)*X^2*Y +
    (24985375916836894020000*nu + 357365666476727973540000)*X^2*Z +
    (24985375916836894020000*nu + 357365666476727973540000)*X*Y^2 +
    (14541507254544649596000*nu + 190702024553013747228000)*X*Y*Z +
    (2102221513957388090400*nu + 25428898291230896666400)*X*Z^2 +
    (2418983748512432766000*nu + 31795515984750073038000)*Y^3 +
    (2099173784599308740400*nu + 25432426477428439316400)*Y^2*Z +
    (603807950218076831675*nu + 6777537793364764438915)*Y*Z^2 +
    (57610935549584088512*nu + 601748568833795438528)*Z^3 +
    157624875243999450000000000*T^3; 
res := Resultant(Q1, E1, X);
R2K:= PolynomialRing(K,2);
R2<x,y>:= PolynomialRing(Rationals(),2);

F2 := y^3+x^5+1;
F:= Evaluate(res, [0, R2K.1, R2K.2, 1]);
I := ideal<RingOfIntegers(K)|Coefficients(F)>;
OK := RingOfIntegers(K);
v:= Decomposition(K,2)[1];
J:= I/ Ideal(v[1]);
_, a := IsPrincipal(J);
F:= F/a;



coeffs, mons := CoefficientsAndMonomials(F);
C2:= Curve(AffineSpace(R2), F2);
C := Curve(AffineSpace(R2K), F);


for p in PrimesInInterval(11,11) do
  print p;
  v := Decomposition(K, p)[1][1];
  exp := InertiaDegree(v );
  k := ResidueClassField(v);
  R2k := PolynomialRing(k,2);
  Fk := &+[Evaluate(c, v) * R2k!mons[i]: i-> c in coeffs];
  Cp := Curve(AffineSpace(R2k), Fk);
  C2p := Curve(Reduction(C2,p));
  Lp := LPolynomial(Cp);
  L2p := LPolynomial(C2p);
  KL := NumberField(Lp);
  KL2 := NumberField(L2p);
  r:= Roots(Lp, KL);
  r2:= Roots(L2p, KL2);
  &or[MinimalPolynomial(r[1][1]^i) eq MinimalPolynomial(r2[1][1]^(exp*i)): i in [1..15]];
end for;

