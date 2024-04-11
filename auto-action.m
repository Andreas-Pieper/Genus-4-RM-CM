function RealPart(A)
  return Matrix(Nrows(A), Ncols(A), [[Real(A[i,j]) : j in [1..Ncols(A)]] : i in [1..Nrows(A)]]);
end function;




function AutoAction()
  R<x> := PolynomialRing(Rationals());
  S := RiemannSurface(x^3-1, 5);
  holo := HolomorphicDifferentials(S);
  CMtype := [-5-5*om[1]+3*om[2]: om in holo];
  X := DiagonalMatrix([Evaluate(K.1^i, plK[1]): i in CMtype ]);
  XR := BlockMatrix(2,2, [RealPart(X), RealPart((-CC.1)*X), RealPart((CC.1)*X), RealPart(X)]);
  XZappr := PiR^-1*ChangeRing(XR, CoefficientRing(PiR))*PiR;
  XZ := Matrix(8,8, [Round(el): el in Eltseq(XZappr)]);
  return XZ, S;



/*  Chains, IntMat, Sympl := HomologyBasis(S);
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
  return Transpose(Sympl)^-1 *  ret * Transpose(Sympl), S;*/
end function;

// copy pasta
/*
load "auto-action.m";

M, S := AutoAction();
Pi := BigPeriodMatrix(S);
K<z> := CyclotomicField(15);
OK := Integers(K);
D := Different(OK);
Dinv := D^-1;
Kre, res := sub< K | K.1+ComplexConjugate(K.1) >;
function foo(x,y)
	return Sign(Real(Evaluate(Kre.1,x) - Evaluate(Kre.1, y)));
end function;
plK := InfinitePlaces(K);
assert [Sign(Imaginary(Evaluate(K.1, p))): p in plK] eq [1,1,1,1];
plKre := InfinitePlaces(Kre);
Sort(~plK, foo);
Sort(~plKre, foo);



OKre := Integers(Kre);
Dre := Different(OKre);
Dreinv := Dre^-1;
_, dre := IsNarrowlyPrincipal(Dreinv);
mini := MinimalPolynomial(K.1,Kre);
disc := Discriminant(mini);
d_new := dre/Sqrt(K!disc);
assert ideal< OK | d_new > eq Dinv;
Pair := [Integers()!Trace(b[i]*ComplexConjugate(b[j])*d_new) : i, j in [1..Degree(K)]] where b := [OK.i : i in [1..Degree(K)]];
Pair := Matrix(8,8,Pair);
id4 := IdentityMatrix(Integers(), 4);
zer := ZeroMatrix(Integers(), 4,4);
J := BlockMatrix([[zer, id4], [-id4, zer]]);
v:= Vector([(J*M^i)[2,2]: i in [0..7]]);
X:= Matrix([[Integers()!Trace(K.1^(i)*OK.j*d_new): i in [0..7]]: j in [1..8] ]);
sol := v*X^-1;
nrm := OK!Eltseq(sol);

//Tr(zeta^i *a* conjugate(a) * d_new)



*/
