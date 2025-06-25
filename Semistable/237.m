function EvaluateCoeff(f, x0 : mult := 1)
  R := Parent(f);
  Coeff := Coefficients(f);
  d := #Coeff-1;
  return &+[R!(mult*Evaluate(Coeff[i+1], x0))*(R.1)^i : i in [0..d]];
end function;

function DivisionCoeff(f, x0)
  R := Parent(f);
  Coeff := Coefficients(f);
  d := #Coeff-1;
  return &+[(Coeff[i+1] div x0)*(R.1)^i : i in [0..d]];
end function;


f := GenericShimuraCurve237();
S<x> := Parent(f);
R<t> := BaseRing(S);
S0<y> := PolynomialRing(QQ);

// fiber at 0
g := EvaluateCoeff(f, t^2);
g := Evaluate(g, t*x);
g := DivisionCoeff(g, t^9);
g := EvaluateCoeff(g, 0);
g := Evaluate(g, y);

g_red := S0!((1*y+1)^10*Evaluate(g, (2*y+1)/(1*y+1))); // make it a degree 10 polynomial by choosing a random transformation

MinRedBinaryForm(g_red);
// -25*x^10 - 126*x^9 - 252*x^8 - 168*x^7 + 252*x^6 + 504*x^5 + 168*x^4 - 288*x^3 - 252*x^2 - 56*x
// gives a smooth curve

// fiber at 1
g := EvaluateCoeff(f, t+1);
g := EvaluateCoeff(g, 0);
g := MinRedBinaryForm(&*[l[1] : l in Factorization(Evaluate(g, y))]);
// 2*y^4 + 3*y^3 - 3*y^2 + 4

S_1<X,Y> := PolynomialRing(QQ, 2);

g := S_1!(Y^4*Evaluate(g, X/Y));
g := Evaluate(g, [X, Y-X]); // put a root at infinity
g := Evaluate(g, [X+Y, Y]); // remove coeff in x^2
g := Evaluate(g, [X, -1/7*Y]);
// X^3*Y - 5/343*X*Y^3 + 2/2401*Y^4
// j-invariant: -3375 -> potential CM by Sqrt(-7)



/////
g := EvaluateCoeff(f, t+1);
Factorization(EvaluateCoeff(g, 0));

g := Parent(g)!((-x+1)^10*Evaluate(g, 2*x/(-x+1)));

Coeff := Coefficients(g);
S_1<X,Y> := PolynomialRing(QQ, 2);
Coeff1 := [Evaluate(c, Y) : c in Coeff];
g := S_1!(&+[Coeff1[i]*X^(i-1) : i in [1..#Coeff]]);
newt := NewtonPolygon(g);
LowerSlopes(newt);

// ListSignatures(Type(newt) : Isa := false);
g := Evaluate(S_1!(Y^6*Evaluate(g, [X/Y^2, Y^7])), [X,0]);
g := x^10*Evaluate(g, [1/x, 0]);
// y^2 = -1152*x^7 + 3456

Terms(S_1!(Y^6*Evaluate(g, [X/Y^2, Y^7])), X);


g := Parent(g)!((-x+1)^10*Evaluate(g, 2*x/(-x+1)));

[Factorization(c) : c in Coefficients(g) | c ne 0];

g := EvaluateCoeff(g, t^7);
g := Parent(g)!(t^12*Evaluate(g, x/t));


EvaluateCoeff(g, 0);
Factorization(EvaluateCoeff(g, 0));




// fiber at infinity
g := EvaluateCoeff(f, 1/t : mult := t^4);

Coeff := Coefficients(g);
S_1<X,Y> := PolynomialRing(QQ, 2);
Coeff1 := [Evaluate(c, Y) : c in Coeff];
g := S_1!(&+[Coeff1[i]*X^(i-1) : i in [1..#Coeff]]);
newt := NewtonPolygon(g);
LowerSlopes(newt);

g := EvaluateCoeff(g, t^3);
g := Parent(g)!(t*Evaluate(g, x/t));
EvaluateCoeff(g, 0);
// y^2 = x^10 - 84*x^7 + 84*x^4 - 28*x
