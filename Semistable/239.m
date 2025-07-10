QQ := RationalsExtra();
R<x> := PolynomialRing(Rationals());
S<t, X,Y,Z,T> := PolynomialRing(QQ, 5);
  Q := X*Z - Y^2;
  E := 3*T^3 + t*(t - 1)*(3*T*(5*X^2 + 6*X*Y + 2*t*Y*Z + 3*t*Z^2) + (-2*t + 9)*X^3 + 22*t*X^2*Y + 21*t*X^2*Z + (-14*t^2 + 18*t)*X*Y*Z + t^2*X*Z^2 + 6*t^2*Y*Z^2 + (-3*t^3 + 6*t^2)*Z^3);
f:=Evaluate(E, [t, Y^2, Y, 1, T]);
disc:= Discriminant(f, T);
Del := Discriminant(disc, Y);
fac := Factorization(Del);
R2:=PolynomialRing(QQ, 2);


//Start with t=0
f1 := Evaluate(f, [t^2, 0, t*Y, 0, T])/3;
f2 := Evaluate(f1, [t, 0, Y, 0, t^2*T]) div t^6;
f3 := Evaluate(f2, [t,0, Y, 0, T-1]);
f4alt := Evaluate(f3, [t, 0, Y, 0, t*T]) div t^2;
//f4alt is a generalized toggle model

CMat := Matrix(Parent(t), [[0,0,0], [1,0,0], [0,1,0]]);
for i in [1..3] do CMat[i,3] := -Coefficients(f3, T)[i]; end for;
q := &+[T^(i-1)*c: i->c in Coefficients(f3, T)[2..4]];
f4 := Determinant(T-Evaluate(q, [t, 0, Y, 0, CMat]));
f5 := Evaluate(f4, [t,0, Y, 0, t*T]) div t^3;
fev := Evaluate(f5, [0,0,Y,0,T]);
g := 3*Factorization(fev)[2][1];
q:=0;
r:= f5;
for i in [4..2 by -1] do;
	q := q+Coefficients(r, T)[i]*T^(i-2);
	r := f5-q*g;
end for;
CMat := Matrix(Parent(t), [[0,0,0], [1,0,0], [0,1,0]]);
for i in [1..3] do CMat[i,3] := -Coefficients(f5, T)[i]; end for;
f6 := Determinant(T-Evaluate(q, [t, 0, Y, 0, CMat]));
f7 := Evaluate(f6, [t^2,0, Y, 0, t*T]) div t^3;
fev := Evaluate(f7, [0,0,Y,0,T]);
fac := Factorization(fev);
Factorization(fac[2][1]-1/27*T^2);
//Get the genus 4 curve y^2 = x*(x^8 - 4*x^6 + 6*x^4 - 44/27*x^2 + 1)
//Does it have CM by Q(i)+Q(zeta_9+zeta_9^(-1), i) ?


//Now t=1
f1 := Evaluate(f, [t+1,0,Y,0,T]);
f2 := Evaluate(f1, [t^3, 0, Y, 0, t*T]) div t^3;
fev := Evaluate(f2, [0,0,Y,0,T]);
Factorization(fev /3-T^3);
//Get the elliptic curve y^2=x^3-1 as a component
f3 := Evaluate(f2, [t,0,Y-1, 0, T]);
f4 := Evaluate(f3, [t^3,0,t^3*Y,0,t^4*T]) div t^12;
fev := Evaluate(f4, [0,0,Y,0,T]);
fev;
//Get the Picard curve y^3=x^4+x as another component. It has CM by Z[zeta_9]



//Now t=infty
f1 := S!(t^5*Evaluate(f, [1/t, 0,Y, 0, T]));
f2 := S!Evaluate(f1, [t^3, 0, Y, 0, T/t^5]);
f3 := S!Evaluate(f2, [t, 0, Y/t, 0, T]);
fev := Evaluate(f3, [0,0,Y,0,T]);
C := Curve(AffineSpace(R2), Evaluate(fev, [0,0,R2.1, 0, R2.2]));
//Get the genus 4 curve with affine equation -2/3*Y^6 + 5*Y^4*T - 14/3*Y^3 + 2*Y*T + T^3 - 1
//Does it have CM by Q(zeta_3)+Q(zeta_9) ?
