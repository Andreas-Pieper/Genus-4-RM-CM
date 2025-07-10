// trying to express the periods of the 2,3,7 Mumford-Shimura family in terms of hypergeometric functions

AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/Genus-4/magma/spec");
AttachSpec("spec");
load "Shimura_models.m";
prec := 300;
RR := RealField(prec);

//A := 2; B := 3; C := 7;
A := 2; B := 7; C := 3;
// double check with Wolfart equation (3)
a := (1/2)*(1 + 1/A - 1/B - 1/C);
b := (1/2)*(1 + 1/A - 1/B + 1/C);
c := 1 + 1/A;

lambda := 1-c;
mu := b-a;
nu := c-b-a;

R<x1,x2> := PolynomialRing(QQ,2);
mons := MonomialsOfDegree(R,3);
t0s := [el/1000 : el in [1..50]]; // put parameters here
M := [];
for t0 in t0s do
  // hypergeometric side
  w1 := BetaFunction(RR!b,RR!c-b)*HypergeometricSeries2F1(a,b,c,t0);
  w2 := (-1)^(1-nu)*t0^lambda*BetaFunction(RR!1+a-c,RR!1-a)*HypergeometricSeries2F1(1+a-c,1+b-c,2-c,t0);
  // period matrix side
  f := GetShimuraCurve237(t0);
  S := RiemannSurface(f,2);
  Pi := BigPeriodMatrix(S);
  //Append(~M, [Evaluate(m,[w1,w2]) : m in mons] cat Eltseq(Pi));
  Append(~M, [Evaluate(m,[w1,w2]) : m in mons] cat [Pi[1,1]]);
end for;
X := Matrix(M);
printf "# rows = %o, #cols = %o\n", Nrows(X), Ncols(X);
Nrows(NumericalKernel(X));
Nrows(NumericalKernel(Transpose(X)));

//Solution(Transpose(X),Transpose(Y));
