//With a rational curve
 R<X> := PolynomialRing(Rationals());
 f1 := X^4 + 36/7*X^3 + 6*X^2 - 4*X + 1;
 K := SplittingField(f1);
 ro := [r[1]: r in Roots(f1, K)];
 Po<X0, X1>:= PolynomialRing(K, 2);
 FS := FieldOfFractions(Po);
 S<x0, x2, a1, a2, a3, a4>:= PolynomialRing(Rationals(), 6);
 T := quo<S|a1*a2*a3*a4-1>;
 A:= Matrix([[-(x0/a1+x2*a1), 1, -(x2/a1+x0*a1)*(x0/a1+x2*a1), (x2/a1+x0*a1)] ,
 [-(x0/a2+x2*a2), 1, -(x2/a2+x0*a2)*(x0/a2+x2*a2), (x2/a2+x0*a2)] ,
 [-(x0/a3+x2*a3), 1, -(x2/a3+x0*a3)*(x0/a3+x2*a3), (x2/a3+x0*a3)] ,
 [-(x0/a4+x2*a4), 1, -(x2/a4+x0*a4)*(x0/a4+x2*a4), (x2/a4+x0*a4)]]);
 B:= ChangeRing(DiagonalMatrix([a1^2,a2^2,a3^2,a4^2])*A, T);
 K:= Kernel(Transpose(B));
 polys := [Evaluate(po, [Po.1, 1] cat ro): po in Eltseq(K.4)];
 Cre := [1/X0, (polys[1]*X1+polys[2])/(X0*(polys[3]*X1+polys[4]))];
 f := f1^2*X;
  R2<ff, x2, x1, xx0, xx1> := PolynomialRing(Rationals(), 5);
  I2 := ideal<R2|[ff-Evaluate(f, x1)*Evaluate(f, x2), x1+x2+xx1, x1*x2-xx0]>; 
  IE := EliminationIdeal(I2, {ff, xx0, xx1});
  fm := [g: g in Generators(IE)|Degree(g, ff) eq 1][1];
  fm, fden := Explode([FS|Evaluate(c, [0,0,0,X0, X1]): c in Coefficients(fm, ff)]);
  f2 := fm / fden; 
  //Pullback with respect to the involution
  f2pri := Evaluate(Numerator(f2), Cre)/Evaluate(Denominator(f2), Cre);
  //Does the Involution lift to W?
  print "Does the involution act by u -> F u?";
  print IsPower(f2/f2pri, 4); //False
  print "\n Does the involution act by u -> F u^-1?";
  print IsPower(f2*f2pri, 4);	//True
  

/* R3 := PolynomialRing(Rationals(),3);
 R7 := PolynomialRing(Rationals(),7);
 List := [X0+Cre[1], X1+Cre[2], X1*Cre2]];
 I6 := ideal<R7| [Evaluate(Denominator(List[i]), [R7.5, R7.6])*R7.i - Evaluate(Numerator(List[i]), [R7.5, R7.6]) :i in [1..3]] cat [(R7.4*R7.7 - R7.7^2)*Evaluate(Denominator(F2), [R7.5, R7.6]) - Evaluate(Numerator(F2), [R7.5, R7.6]), R7.7^4-Evaluate(Numerator(f2), [R7.5, R7.6])]>;
 */

Crel := [Cre[1], Evaluate(Cre[2], [Po.1^2, Po.2]) ];
f2l := Evaluate(Sqrt(-f2/Po.1), [Po.1^2, Po.2])*Po.1;
R2 := PolynomialRing(Rationals(),2);
f2l := R2!(Numerator(f2l))/R2!(Denominator(f2l));

I := Matrix([[Coefficients(l, Po.2)[i]: i in [2,1]]: l in [Numerator(Crel[2]), Denominator(Crel[2])]]);
Isig := Matrix(2,2, [Evaluate(c, [1/X0, 0]): c in Eltseq(I)]);
brau := (I*Isig)[1,1];
L := OptimizedRepresentation( SplittingField(Evaluate(Numerator(brau), [X, 0])));
G, Aut, tau := AutomorphismGroup(L);
KK:= FixedField(L, [tau(g): g in Subgroups(G)[3]`subgroup]);
fac := Factorization(Evaluate(Numerator(brau), [X, 0]), KK);
per := [[i,j]: i, j in [1..4]|fac[i][1] eq -Evaluate(fac[j][1], -Parent(fac[j][1]).1)];
spli := fac[per[1][1]][1] * fac[per[1][2]][1];
RKK := PolynomialRing(KK);
FKK := FieldOfFractions(RKK);
IKK:= 1/spli * Matrix(2,2, [FKK|Evaluate(c, [RKK.1,0]): c in Eltseq(I)]);
IKKsig := Matrix(2,2, [Evaluate(c, 1/RKK.1): c in Eltseq(IKK)]);
bool, sqr := IsSquare((IKK* IKKsig)[1,1]);
assert bool;
spli := sqr*spli;
IKK := 1/sqr * IKK;
IKKsig := 1/sqr *IKKsig;

function rec(foo)
	eva := Evaluate(foo, (RKK.1-1)/(RKK.1+1));
        num := Numerator(eva);
        den := Denominator(eva);
end function;

function trah(foo)
       res :=  1/2*(foo + Evaluate(foo, 1/RKK.1));
       return res;	
end function;
a := RKK.1-1/RKK.1;
IKKoubl := BlockMatrix(2,2, [Matrix(FKK, [[trah(c), a^2*trah(c/a)],[trah(c/a), trah(c)]]): c in Eltseq(IKK)]);
D:= DiagonalMatrix([FKK|1,-1,1, -1]);
Ke := Kernel(Transpose(IKKoubl- D));
SS := Transpose(Matrix([[Eltseq(Ke.i)[1]+ Eltseq(Ke.i)[2]*a, Eltseq(Ke.i)[3]+ a*Eltseq(Ke.i)[4]]: i in [1..2]]));
SSsig := Matrix(2,2, [Evaluate(c, 1/RKK.1): c in Eltseq(SS)]);
assert IKK*SS eq SSsig;
R2KK := PolynomialRing(KK,2);
SSc := [Evaluate(c, R2KK.1): c in Eltseq(SS)];
subs := [R2KK.1, (R2KK.2*SSc[1]+SSc[2])/(R2KK.2*SSc[3]+SSc[4])];
f2leva := Numerator(Evaluate(f2l, subs));
f2ltra := Evaluate(f2leva, [(R2KK.1+1)/(R2KK.1-1), R2KK.2]);
fmin := Evaluate(f2ltra, [-R2KK.1, R2KK.2]);
_, g := IsSquare(fmin*f2ltra);
res := f2ltra + 2*g+ fmin;
assert res eq Evaluate(res, [-R2KK.1, R2KK.2]);
coeffnum := Coefficients(Numerator(res), R2KK.1);
coeffden := Coefficients(Denominator(res), R2KK.1);
resnum := &+[R2KK.1^i*coeffnum[2*i+1]: i in [0..(#coeffnum div 2)]];
resden := &+[R2KK.1^i*coeffden[2*i+1]: i in [0..(#coeffden div 2)]];
sqfnum := SquareFreeFactorisation(resnum);
sqfden := SquareFreeFactorisation(resden);
Ram := &*[fa[1]^(fa[2] mod 2): fa in sqfnum cat sqfden];
RFKK := PolynomialRing(FKK);
RamRFKK := Evaluate(Ram, [FKK.1+1, RFKK.1]);
C := GenusOneModel(RamRFKK);
E:= Jacobian(C);
MinE:= MinimalDegreeModel(E);


E := MinimalModel(WeierstrassModel(E));
L := SplittingField(Numerator(Discriminant(E)));
FL := RationalFunctionField(L);
EFL := ChangeRing(E, FL);
alpha := [Numerator(r[1]): r in Roots(DivisionPolynomial(EFL, 2))];
rr1 := Roots(Numerator((alpha[2]-alpha[1])/Universe(alpha).1));
rr3 := Roots(Numerator((alpha[2]-alpha[3])/Universe(alpha).1));
pts1 := [[Evaluate(alpha[2], r[1]), r[1]] : r in rr1];
RL<x,t> := PolynomialRing(L,2);
f0 := (x-Evaluate(Derivative(alpha[3]), 0)*t)^3;
f1 := (x-LeadingCoefficient(alpha[2])*t^3)^2*(x-LeadingCoefficient(alpha[1])*t^3);
f1 -:= Evaluate(Derivative(alpha[3]), 0)*MonomialCoefficient(f1, x^2*t^3)*x*t^4;
f01 := f0+f1-x^3;
pts3 := [[Evaluate(alpha[2], r[1]), r[1]] : r in rr3];


f1p := (x-LeadingCoefficient(Numerator(alpha[2]))*t^3)*t^5;
f1p1 := (x-LeadingCoefficient(Numerator(alpha[2]))*t^3)^2*t^2;
f0p := (x-Evaluate(Derivative(alpha[3]), 0)*t)^2*t^2;
f01p := f0p+f1p1-x^2*t^2;
f0pp := (x-Evaluate(Derivative(alpha[3]), 0)*t)*t^4;
list := [f01, f01p, f1p, f0pp, t^7, t^6];
Mat := Matrix([[Evaluate(l, pt): pt in pts1 cat pts3]: l in list]);
Ker := Kernel(Mat);
polys := [&+[list[i]*Eltseq(Ker.j)[i]: i in [1..6]]: j in [1..2]];
ram := [x-Evaluate(a, t) : a in alpha ];

eva := [Evaluate(pol, [x+ConstantCoefficient(Derivative(alpha[3]))*t+LeadingCoefficient(alpha[2])*t^3, t]): pol in polys];
evaram := [Evaluate(pol, [x+ConstantCoefficient(Derivative(alpha[3]))*t+LeadingCoefficient(alpha[2])*t^3, t]): pol in ram ];
eva1 := [Evaluate(pol, [x*t^2, t]) div t^6: pol in eva];
eva1ram := [Evaluate(pol, [x*t^2, t]): i->pol in evaram];

assert &and[Degree(pol, t) eq 1 : pol in eva1];
RFL := PolynomialRing(FL);
den := [Term(pol ,t,1) div t: pol in eva1];
num := [Term(pol ,t,0): pol in eva1];
param := -(Evaluate(num[1], [RFL.1,0])-FL.1*Evaluate(num[2], [RFL.1,0]))/(Evaluate(den[1], [RFL.1,0])-FL.1*Evaluate(den[2], [RFL.1,0]));
rami := [Evaluate(pol, [RFL.1, param]): pol in eva1ram];
rami1 := [Numerator(ra)*Denominator(ra): ra in rami];

sqf := [SquareFreeFactorization(ra): ra in rami1];
rami1sqf := [rami1[i] div &*[sq[1]^(sq[2] div 2): sq in s]^2: i->s in sqf];
ramm12 := rami1sqf[1]*rami1sqf[2] div GCD(rami1sqf[1], rami1sqf[2])^2;
ramm := ramm12*rami1sqf[3] div GCD(ramm12, rami1sqf[3])^2;

 C := GenusOneModel(ramm);
 E:= Jacobian(C);
 fac := Factorization(Numerator(Discriminant(E)));
 r0 := [Roots(f[1])[1,1]: f in fac| f[2] eq 1];
 a := 4/(r0[2]-r0[1]);
 b := 2-a*r0[2];
 ainvs := aInvariants(E);
ainvspr := [Evaluate(ainv, (FL.1-b)/a): ainv in ainvs];
Epr:= EllipticCurve(ainvspr);
E0 := EllipticCurve(RFL.1^3-3*RFL.1+Evaluate(DicksonFirst(7,1), FL.1));
jInvariant(Epr) eq jInvariant(E0);

