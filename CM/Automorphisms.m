K<a,b,c,d,t> := PolynomialRing(QQ, 6);
S<X,Y> := PolynomialRing(K, 2);

f := 1/3*t*(t-1)*((-2*t + 9)*X^6 + 22*t*X^5*Y + 21*t*X^4*Y^2 + (-14*t^2 + 18*t)*X^3*Y^3 + t^2*X^2*Y^4 + 6*t^2*X*Y^5 + (-3*t^3 + 6*t^2)*Y^6);
v := t*(t - 1)*(5*X^4 + 6*X^3*Y + 2*t*X*Y^3 + 3*t*Y^4);

q1 := Transvectant(f, Transvectant(f, f, 4), 4); // deg 3 in f
q2 := Transvectant(q1, Transvectant(f, f, 4), 2); // deg 5 in f
q3 := Transvectant(q2, Transvectant(f, f, 4), 2); // deg 7 in f

Q := [q1, q2, q3];
Q1 := [Evaluate(Q[i], [a*X+b*Y, c*X+d*Y]) : i in [1..3]];

Q0 := [Q[i]-Q1[i] : i in [1..3]];

// 10th invariant is of weight 10, 11th of weight 9, and they must remain the same, thus det^9 = det^10 = 1, which implies det = 1.
// so we're only interested in the elements in SL_2 which stabilize the quadratic forms q1, q2, q3 simultaneously

I := Ideal(&cat[Coefficients(Q0[i]) : i in [1..3]] cat [(a*d-b*c)-1]);
RadicalDecomposition(I);
 
/*
[
Ideal of Polynomial ring of rank 6 over Rational Field
Order: Lexicographical
Variables: a, b, c, d, t, $.6
Inhomogeneous, Dimension 4, Radical, Prime
Groebner basis:
[
a*d - b*c - 1,
t - 1
],
Ideal of Polynomial ring of rank 6 over Rational Field
Order: Lexicographical
Variables: a, b, c, d, t, $.6
Inhomogeneous, Dimension 4, Radical, Prime
Groebner basis:
[
a*d - b*c - 1,
t
],
Ideal of Polynomial ring of rank 6 over Rational Field
Order: Lexicographical
Variables: a, b, c, d, t, $.6
Inhomogeneous, Dimension 2, Radical, Prime
Groebner basis:
[
a + 1,
b,
c,
d + 1
],
Ideal of Polynomial ring of rank 6 over Rational Field
Order: Lexicographical
Variables: a, b, c, d, t, $.6
Inhomogeneous, Dimension 2, Radical, Prime
Groebner basis:
[
a - 1,
b,
c,
d - 1
]
]
*/

// so when t is not 0 or 1, the automorphism group of the corresponding fibre is trivial.
