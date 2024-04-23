function NormalizedCMBigPeriod( : prec:= 100)
  K<z> := CyclotomicField(15);
  OK := Integers(K);
  D := Different(OK);
  Dinv := D^-1;
  Kre, res := sub< K | K.1+ComplexConjugate(K.1) >;
  function foo(x,y)
        return Sign(Real(Evaluate(Kre.1,y) - Evaluate(Kre.1, x)));
  end function;
  plK := InfinitePlaces(K);
  assert [Sign(Imaginary(Evaluate(K.1, p))): p in plK] eq [1,1,1,1];
  plKre := InfinitePlaces(Kre);
  Sort(~plK, foo);
  Sort(~plKre, foo);
  OKre := Integers(Kre);
  /*Dre := Different(OKre);
  Dreinv := Dre^-1;
  _, dre := IsNarrowlyPrincipal(Dreinv);
  mini := MinimalPolynomial(K.1,Kre);
  disc := Discriminant(mini);
  d_new := dre/Sqrt(K!disc);*/
  d_new := K![ -1/5, 1/3, -1/3, -1/15, 1/3, -2/3, 3/5, -1/3 ];
  assert ideal< OK | d_new > eq Dinv;
  assert [Sign(Imaginary(Evaluate(d_new, plK[i]))): i in [1..4]] eq [-1,-1,-1,-1];
  Pair := Matrix(8,8, [Integers()!Trace(OK.i*ComplexConjugate(OK.j)*d_new) : i, j in [1..Degree(K)]]);
 // _, Sy := FrobeniusFormAlternating(Pair);
  Sy := Matrix(8, 8, [ 0, 1, 0, 0, 0, 0, 0, 0,
		       1, 1, 1, 0, 0, 0, 0, 0,
		       2, 1, -1, 1, 2, 0, 0, 0,
		       1, 1, 0, 0, -1, 0, 2, 1,
		       1, 0, 0, 0, 0, 0, 0, 0,
		       2, 2, 0, 0, 1, 0, 0, 0,
		       0, 0, 1, 0, -2, 0, 1, 0,
		       3, 1, -1, 1, -1, 1, 3, 0 ]);
  id4 := IdentityMatrix(Integers(), 4);
  zer := ZeroMatrix(Integers(), 4,4);
  J := BlockMatrix([[zer, id4], [-id4, zer]]);
  function RealPart(A)
    return Matrix(Nrows(A), Ncols(A), [[Real(A[i,j]) : j in [1..Ncols(A)]] : i in [1..Nrows(A)]]);
  end function;
  

  CC := ComplexField(prec);
  R<x> := PolynomialRing(Rationals());
  S := RiemannSurface(x^3-1, 5: Precision:=prec);
  Pi := BigPeriodMatrix(S);
  holo := HolomorphicDifferentials(S);
  CMtype := [-5-5*om[1]+3*om[2]: om in holo];
  X := DiagonalMatrix([Evaluate(K.1^i, plK[1]): i in CMtype ]);
  XR := BlockMatrix(2,2, [RealPart(X), RealPart((-CC.1)*X), RealPart((CC.1)*X), RealPart(X)]);
  PiR := VerticalJoin([RealPart(Pi), RealPart((-CC.1)*Pi)]);
  XZappr := PiR^-1*ChangeRing(XR, CoefficientRing(PiR))*PiR;
  XZ := Matrix(8,8, [Round(el): el in Eltseq(XZappr)]);
  v:= Vector([(J*XZ^(-i))[1,1]: i in [0..7]]);
  X:= Matrix([[Integers()!Trace(K.1^(-i)*OK.j*d_new): i in [0..7]]: j in [1..8] ]);
  sol := v*X^-1;
  nrm := -K!Eltseq(sol);
  assert Norm(nrm) eq 1;
  mini := MinimalPolynomial(K.1,Kre);
  F:= ext<Kre|mini>;
  OF := RingOfIntegers(F);
  b, u := NormEquation(RingOfIntegers(F), OKre!nrm);
  assert b;
  u := Kre!u[1][1]+K!u[1][2]*K.1;
  M := &+[XZ^(i-1)*n:  i->n in Eltseq(u)];
  Sy2 := Matrix([Transpose(M^-1*XZ^i)[1]: i in [0..7]]);
  Sytot := Transpose(Sy2) * Transpose(Sy);
  assert Transpose(Sytot) *J*Sytot eq J;
  return Pi * ChangeRing(Sytot, CC), Sytot^-1* XZ *Sytot;
end function;



