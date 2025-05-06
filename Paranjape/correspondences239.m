AttachSpec("~/github/CHIMP/CHIMP.spec");
SetVerbose("EndoCheck",true);

QQ := Rationals();
R<x> := PolynomialRing(Rationals());
/*
S<t, X,Y,Z,T> := PolynomialRing(QQ, 5);
  Q := X*Z - Y^2;
  E := 3*T^3 + t*(t - 1)*(3*T*(5*X^2 + 6*X*Y + 2*t*Y*Z + 3*t*Z^2) + (-2*t + 9)*X^3 + 22*t*X^2*Y + 21*t*X^2*Z + (-14*t^2 + 18*t)*X*Y*Z + t^2*X*Z^2 + 6*t^2*Y*Z^2 + (-3*t^3 + 6*t^2)*Z^3);
f:=Evaluate(E, [t, Y^2, Y, 1, T]);
*/


// Fiber over t=1
/*
f1 := Evaluate(f, [t+1,0,Y,0,T]);
f2 := Evaluate(f1, [t^3, 0, Y, 0, t*T]) div t^3;
fev := Evaluate(f2, [0,0,Y,0,T]);
*/
//Factorization(fev /3-T^3);
R3<u,v> := PolynomialRing(QQ,2);
f := 7*u^6 + 22*u^5 + 21*u^4 + 4*u^3 + u^2 + 6*u + 3*v^3 + 3;
C4 := Curve(Spec(R3),f);
E := EllipticCurve([0,-1]);
//Get the elliptic curve y^2=x^3-1 as a component
/*
f3 := Evaluate(f2, [t,0,Y-1, 0, T]);
f4 := Evaluate(f3, [t^3,0,t^3*Y,0,t^4*T]) div t^12;
fev := Evaluate(f4, [0,0,Y,0,T]);
fev;
*/
//Get the Picard curve y^3=x^4+x as another component. It has CM by Z[zeta_9]
F3 := v^3 - (u^4+u);
C3 := Curve(Spec(R3),F3);

// t = oo
QQ := RationalsExtra();
CC := QQ`CC;
R<x> := PolynomialRing(QQ);
S<t, X,Y,Z,T> := PolynomialRing(QQ, 5);
  Q := X*Z - Y^2;
  E := 3*T^3 + t*(t - 1)*(3*T*(5*X^2 + 6*X*Y + 2*t*Y*Z + 3*t*Z^2) + (-2*t + 9)*X^3 + 22*t*X^2*Y + 21*t*X^2*Z + (-14*t^2 + 18*t)*X*Y*Z + t^2*X*Z^2 + 6*t^2*Y*Z^2 + (-3*t^3 + 6*t^2)*Z^3);
f:=Evaluate(E, [t, Y^2, Y, 1, T]);
disc:= Discriminant(f, T);
Del := Discriminant(disc, Y);
fac := Factorization(Del);
R2:=PolynomialRing(QQ, 2);
f1 := S!(t^5*Evaluate(f, [1/t, 0,Y, 0, T]));
f2 := S!Evaluate(f1, [t^3, 0, Y, 0, T/t^5]);
f3 := S!Evaluate(f2, [t, 0, Y/t, 0, T]);
fev := Evaluate(f3, [0,0,Y,0,T]);
f4 := Evaluate(fev, [0,0,R2.1, 0, R2.2]);
RS := RiemannSurface(f4 : Precision := Precision(CC));
Pi := BigPeriodMatrix(RS);
Pi := ChangeRing(Pi, CC);
GeometricEndomorphismRepresentationCC(Pi);

// genus 1 factor

E1 := EllipticCurve([0,-1]);
matsE, mapsE := GeometricHomomorphismRepresentation(PeriodMatrix(E1), Pi, QQ);

// genus 1 factor
P:= E![1,0];
Q := HK![1,nu,1];
matsH, mapsH := GeometricHomomorphismRepresentation(PeriodMatrix(H), PeriodMatrix(C), QQ);
i0 := Position([Max([Degree(MinimalPolynomial(a)): a in Eltseq(m[1])]): m in matsH],1);
homo := [*ChangeRing(matsH[i0][1], RationalsExtra()), matsH[i0][2]*];
_, corrH := Correspondence(HK, C, homo : P := Q, Q := P, Al := "Cantor");
corrH;




//b, corr := Correspondence(CK, CK, endo[2]: P:= P, Q:= P, Al := "Cantor");

E := EllipticCurve([1,0]);
mats, maps := GeometricHomomorphismRepresentation(PeriodMatrix(C), PeriodMatrix(E), QQ);
_, corr := Correspondence(CK, HyperellipticCurve(E), mats[1] : P := P, Al := "Cantor");

//matsH, mapsH := GeometricHomomorphismRepresentation(PeriodMatrix(C), PeriodMatrix(H), QQ);
//matsH, mapsH := GeometricHomomorphismRepresentation(PeriodMatrix(H), PeriodMatrix(C), QQ);
//_, corrH := Correspondence(CK, H, matsH[1] : P := P, Al := "Cantor");
//_, corrH := Correspondence(CK, H, matsH[2] : P := P, Al := "Cantor");
//
//
//
/*lhs := Matrix([[(&cat[Eltseq(m): m in Eltseq(M[1])])[i]: i in [1..144]]: M in 
matsH]);
rhs := Matrix([[(&cat[Eltseq(L!m): m in Eltseq(M)])[i]: i in [1..144]]: M in 
Basis(KM)]);
X := VerticalJoin(rhs, lhs);
Kernel(X);
subs := Subfields(L);
K:=subs[2][1];
map:=subs[2][2];
rhs := Matrix(&cat[[[(&cat[Eltseq(L!m): m in 
Eltseq(map(a)*ChangeRing(M,K))])[i]: i in [1..144]]: M in Basis(KM)]: a in 
Basis(K)]);
X := VerticalJoin(lhs, rhs);
Ke := Kernel(X);
lin := Vector(Eltseq(Ke.1)[1..6]);

homo := [*ChangeRing(matsH[2][1], K), matsH[2][2]*];



_, corrH := Correspondence(HL, CK, homo : P := Q, Q := P, Al := "Cantor");*/
