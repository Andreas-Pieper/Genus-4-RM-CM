// relating Shioda's example y^2 = x^9 - 1 to our family
// Cf., Zhu's paper Constructing a CM Mumford fourfold from Shioda's fourfold Prop 4.1

AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/Genus-4/magma/spec");
AttachSpec("../CM/spec");

// copy pasta

QQ := RationalsExtra();
R<x> := PolynomialRing(QQ);
K<z> := NumberFieldExtra(x^6 + x^3 + 1);
S<t, X,Y,Z,T> := PolynomialRing(QQ, 5);
  Q := X*Z - Y^2;
  E := 3*T^3 + t*(t - 1)*(3*T*(5*X^2 + 6*X*Y + 2*t*Y*Z + 3*t*Z^2) + (-2*t + 9)*X^3 + 22*t*X^2*Y + 21*t*X^2*Z + (-14*t^2 + 18*t)*X*Y*Z + t^2*X*Z^2 + 6*t^2*Y*Z^2 + (-3*t^3 + 6*t^2)*Z^3);
f:=Evaluate(E, [t, Y^2, Y, 1, T]);

disc:= Discriminant(f, T);
Del := Discriminant(disc, Y);
fac := Factorization(Del);
R2:=PolynomialRing(QQ, 2);

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
fev; // genus 3 factor

// computing possible polarizations; there is only 1
// as shown by number of totally positive elements in 
// quotient below
K;
ZZK := Integers(K);
U, phiU := UnitGroup(ZZK : GRH := true);

inv := [el : el in Automorphisms(K) | el(z) eq z^-1][1];
U, phiU := UnitGroup(ZZK : GRH := true);

gensU := Generators(U);
gensV := [ (phiU(gen)*inv(phiU(gen))) @@ phiU : gen in gensU ];
V := sub< U | gensV >;
Q, pQ := quo< U | V >;
Q;
K0 := sub<K | z + z^-1>;
real := [el : el in Q | phiU(el @@ pQ) in K0];
[IsTotallyPositive(K0!phiU(el @@ pQ)) : el in real];
