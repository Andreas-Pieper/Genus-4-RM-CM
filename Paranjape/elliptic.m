//With a rational curve
/* R<X> := PolynomialRing(Rationals());
 f1 := X^4 + 36/7*X^3 + 6*X^2 - 4*X + 1;
 K := SplittingField(f1);
 ro := [r[1]: r in Roots(f1, K)];
 Po<X0, X1, X2>:= PolynomialRing(K, 3);
 S<x0, x2, a1, a2, a3, a4>:= PolynomialRing(Rationals(), 6);
 T := quo<S|a1*a2*a3*a4-1>;
 A:= Matrix([[-(x0/a1+x2*a1), 1, -(x2/a1+x0*a1)*(x0/a1+x2*a1), (x2/a1+x0*a1)] ,
 [-(x0/a2+x2*a2), 1, -(x2/a2+x0*a2)*(x0/a2+x2*a2), (x2/a2+x0*a2)] ,
 [-(x0/a3+x2*a3), 1, -(x2/a3+x0*a3)*(x0/a3+x2*a3), (x2/a3+x0*a3)] ,
 [-(x0/a4+x2*a4), 1, -(x2/a4+x0*a4)*(x0/a4+x2*a4), (x2/a4+x0*a4)]]);
 B:= ChangeRing(DiagonalMatrix([a1^2,a2^2,a3^2,a4^2])*A, T);
 K:= Kernel(Transpose(B));
 polys := [Evaluate(po, [Po.1, Po.3] cat ro): po in Eltseq(K.4)];
 Cre := [X2* (polys[3]*X1+polys[4]), (polys[1]*X1+polys[2]), X0*(polys[3]*X1+polys[4])];

[
    X0^2*X1*X2 - 36/7*X0^2*X2^2 - 4*X0*X2^3 - X1*X2^3,
    X0^4 - 6*X0^3*X2 - 4*X0^2*X1*X2 - 36/7*X0*X1*X2^2 + 6*X0*X2^3 - X2^4,
    X0^3*X1 - 36/7*X0^3*X2 - 4*X0^2*X2^2 - X0*X1*X2^2
]*/

//With an elliptic curve
  E := EllipticCurve([0, 0, 0, -1889771, 961189146]);
  FE := FunctionField(E); 
  T := E![119, 27166];
  P1 := E![629, 4624];
  P2 := E![965, 6016];
  P3 := E![-774, 44274];
  P4 := -P1-P2-P3;
  Pi := [P1, P2, P3, P4];
  R<y,x,t> := PolynomialRing(Rationals(), 3);
  S := quo<R|-y^2+x^3-1889771*x+961189146>;
  FS := FieldOfFractions(S);
  A:= Matrix([[(Pi[i][2]+y)/(Pi[i][1]-x), 1, -(Pi[i][2]+y)/(Pi[i][1]-x)*(Pi[i][2]-y)/(Pi[i][1]-x), -(Pi[i][2]-y)/(Pi[i][1]-x)] : i in [1..4]]);
  A:= ChangeRing(LCM([Denominator(c): c in Eltseq(A)])*A, S);
  K:=Kernel(Transpose(A));
  polys := ChangeUniverse(Eltseq(K.1), R);
  Cre := [FS|-y, x, (polys[1]*t+polys[2])/(polys[3]*t+polys[4])];
  D := &+[2*Divisor(P): P in Pi]-4*Divisor(E!0)-4*Divisor(T);
  bool, f := IsPrincipal(D);
  assert bool;
  //Write f(x1, y1)*f(x2, y2) as a rational function in (x3, y3, t). Here (x3,y3)=(x1,y1) + (x_2, y_2) where the the sum is taken in E.
  num := Numerator(f);
  den := Denominator(f);
  R2<ff, y2, x2, y1, x1, y3, x3, tt> := PolynomialRing(Rationals(), 8);
  I2 := ideal<R2|[ff- (Evaluate(den, [x1,y1])* Evaluate(den, [x2,y2])), -y1^2+x1^3-1889771*x1+961189146,-y2^2+x2^3-1889771*x2+961189146,-y3^2+x3^3-1889771*x3+961189146, -x3+tt^2-x1-x2,-y3+tt*(x1-x3)-y1, -(y1-y2)+tt*(x1-x2)]>;
  IE := EliminationIdeal(I2, {ff, x3,y3,tt});
  denm := [g: g in Generators(IE)|Degree(g, ff) eq 1][1];
  dennum, denden := Explode([FS|Evaluate(c, [0,0,0,0,0,y,x,t]): c in Coefficients(denm, ff)]);
  I2 := ideal<R2|[ff- (Evaluate(num, [x1,y1])* Evaluate(num, [x2,y2])), -y1^2+x1^3-1889771*x1+961189146,-y2^2+x2^3-1889771*x2+961189146,-y3^2+x3^3-1889771*x3+961189146, -x3+tt^2-x1-x2,-y3+tt*(x1-x3)-y1, -(y1-y2)+tt*(x1-x2)]>;
  IE := EliminationIdeal(I2, {ff, x3,y3,tt});
  numm := [g: g in Generators(IE)|Degree(g, ff) eq 1][1];
  numnum, numden := Explode([FS| Evaluate(c, [0,0,0,0,0,y,x,t]): c in Coefficients(numm, ff)]);
  f2 := (numnum*denden)/(numden*dennum);
  
  //Pullback with respect to the involution
  f2pri := Evaluate(Numerator(f2), Cre)/Evaluate(Denominator(f2), Cre);
  //Does the Involution lift to W?
  print "Does the involution act by u -> F u?";
  print IsPower(f2/f2pri, 4); //False
  print "\n Does the involution act by u -> F u^-1?";
  print IsPower(f2*f2pri, 4);	//True
