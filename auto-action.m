function AutoAction()
  R<x> := PolynomialRing(Rationals());
  S := RiemannSurface(x^3-1, 5);
  Chains, IntMat, Sympl := HomologyBasis(S);
  Rams := [Coordinates(ram)[1]: ram in RamificationPoints(S)];
  EdgesList := [edge`EP: edge in Chains`Edges];
  i0:= &meet [Seqset(S): S in EdgesList];
  assert #i0 eq 1;
  i0 := Random(i0);
  assert &and [edge[1] eq i0: edge in EdgesList];
  i1 := EdgesList[1][2];
  i2 := EdgesList[2][2];
  angle := (Rams[i1]-Rams[i0])/(Rams[i2]-Rams[i0]);
  orientation := Sign(Imaginary(angle));
  idx := [i: i in [5..8]|IntMat[1, i] eq orientation];
  assert #idx eq 1;
  idx := idx[1];
  path1 := Matrix(8,1, [-1,0,0,0,0,0,0,0]); 
  path1[idx,1] := 1;
  path2 := Matrix(8,1, [-1,0,0,0,0,0,0,0]);
  X0 := Matrix([[0,0,0,-1], [1,0,0,-1], [0,1,0,-1], [0,0,1,-1]]);
  X := Matrix(BlockDiagMat(<X0, X0>));
  for i in [0..4] do
    Y := HorizontalJoin([X^(i+j)*path1: j in [0..3]] cat [X^j*path2: j in [0..3]]);
    if Transpose(Y) * IntMat * Y eq IntMat then
       ret :=  Y;
    end if;
  end for;
  return Transpose(Sympl)^-1 *  ret * Transpose(Sympl), S;
end function;

// copy pasta
/*
load "auto-action.m";
Pi, S := AutoAction();
Pi;
M, S := AutoAction();
Pi := BigPeriodMatrix(S);
Nrows(M);
M^15;
M^5;
M^3;
K<z> := CyclotomicField(15);
Codifferent;
Dinv := Codifferent(K);
D := Different(K);
OK := Integers(K);
Dinv := Codifferent(K);
Dinv := Codifferent(OK);
D := Different(OK);
D^-1;
Dinv := $1;
IsPrincipal(Dinv);
_, d := IsPrincipal(Dinv);
d;
IsTotallyReal(d*I);
[Evaluate(d,v) : v in InfinitePlaces(K)];
[Re(Evaluate(d,v)) : v in InfinitePlaces(K)];
IsNarrowlyPrincipal(Dinv);
_, d:= $1;
[Im(Evaluate(d,v)) : v in InfinitePlaces(K)];
_, d := IsPrincipal(Dinv);
d;
ComplexConjugate(d);
d*ComplexConjugate(d);
IsTotallyReal($1);
[Im(Evaluate(d*ComplexConjugate(d),v)) : v in InfinitePlaces(K)
];
IsTotallyPositive(d*ComplexConjugate(d));
TotallyRealSubAlgebra;
sub< K | d*ComplexConjugate(d) >;
Kre, res := $1;
Degree(K);
Degree(Kre);
K.1;
sub< K | K.1*ComplexConjugate(K.1) >;
Kre, res := $1;
Degree(Kre);
sub< K | K.1+ComplexConjugate(K.1) >;
Kre, res := $1;
Degree(Kre);
OKre := Integers(Kre);
Dreinv := Codifferent(OKre);
Dre := Different(OKre);
Dreinv := Dre^-1;
Dreinv;
IsNarrowlyPrincipal(Dreiv);
IsNarrowlyPrincipal(Dreinv);
_, dre := $1;
MinimalPolynomial(K.1,Kre);
Discriminant($1);
disc := $1;
dre/disc;
R<x> := PolynomialRing(K);
Roots(K, x^2 - disc);
Roots(x^2 - disc,K);
Roots(x^2 - disc,K)[2][1];
r := $1;
dre/r;
d_new := $1;
[Re(Evaluate(d_new,v)) : v in InfinitePlaces(K)];
ideal< OK | d_new > eq Dinv;
M := [];
for i := 1 to Rank(OK) do
for j := 1 to Rank(OK) do
for i := 1 to Rank(OK) do
M := [Trace(b[i]*ComplexConjugate(b[j])*d_new) : i, j in [1..Rank(OK)]] where b := [OK.i : i in [1..Rank(OK)]];
M := [Trace(b[i]*ComplexConjugate(b[j])*d_new) : i, j in [1..Degree(K)]] where b := [OK.i : i in [1..Degree(K)]];
M := Matrix(8,8,M);
*/
