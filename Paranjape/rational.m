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
  

 R3 := PolynomialRing(Rationals(),3);
 R5 := PolynomialRing(Rationals(),5);
 List := [X0+Cre[1], X1+Cre[2], X1*Cre[2]];
 I5 := ideal<R5| [R5.i * Evaluate(Denominator( List[i]), [R5.4, R5.5])-Evaluate(Numerator(List[i]), [R5.4, R5.5]): i in [1..3]]>;

