//AttachSpec("~/github/CHIMP/CHIMP.spec");
SetVerbose("EndoCheck",true);
R<X> := PolynomialRing(Rationals());   
f1 := X^4 + 36/7*X^3 + 6*X^2 - 4*X + 1;
F := Evaluate(f1, X^2)*X;
F:= 7*F;
C := HyperellipticCurve(F);
P:= C![1, 8,1];
H := HyperellipticCurve(2*(X^7 + 7*X^5 + 14*X^3 + 7*X));
K<nu>:= NumberFieldExtra(X^2-58);
assert nu^2 eq 58;
HK := ChangeRing(H,K);
Q := HK![1,nu,1];
matsH, mapsH := GeometricHomomorphismRepresentation(PeriodMatrix(H), PeriodMatrix(C), QQ);
i0 := Position([Max([Degree(MinimalPolynomial(a)): a in Eltseq(m[1])]): m in matsH],1);
homo := [*ChangeRing(matsH[i0][1], RationalsExtra()), matsH[i0][2]*];
_, corrH := Correspondence(HK, C, homo : P := Q, Q := P, Al := "Cantor");
corrH;




//b, corr := Correspondence(CK, CK, endo[2]: P:= P, Q:= P, Al := "Cantor");

E := EllipticCurve([1,0]);
mats, maps := GeometricHomomorphismRepresentation(PeriodMatrix(C), PeriodMatrix(E), QQ);
_, corr := Correspondence(CK, HyperellipticCurve(E), mats[1] : P := P, Al := "Cantor");

//matsH, mapsH := GeometricHomomorphismRepresentation(PeriodMatrix(C), PeriodMatrix(H), QQ);
//matsH, mapsH := GeometricHomomorphismRepresentation(PeriodMatrix(H), PeriodMatrix(C), QQ);
//_, corrH := Correspondence(CK, H, matsH[1] : P := P, Al := "Cantor");
//_, corrH := Correspondence(CK, H, matsH[2] : P := P, Al := "Cantor");
//
//
//
/*lhs := Matrix([[(&cat[Eltseq(m): m in Eltseq(M[1])])[i]: i in [1..144]]: M in 
matsH]);
rhs := Matrix([[(&cat[Eltseq(L!m): m in Eltseq(M)])[i]: i in [1..144]]: M in 
Basis(KM)]);
X := VerticalJoin(rhs, lhs);
Kernel(X);
subs := Subfields(L);
K:=subs[2][1];
map:=subs[2][2];
rhs := Matrix(&cat[[[(&cat[Eltseq(L!m): m in 
Eltseq(map(a)*ChangeRing(M,K))])[i]: i in [1..144]]: M in Basis(KM)]: a in 
Basis(K)]);
X := VerticalJoin(lhs, rhs);
Ke := Kernel(X);
lin := Vector(Eltseq(Ke.1)[1..6]);

homo := [*ChangeRing(matsH[2][1], K), matsH[2][2]*];



_, corrH := Correspondence(HL, CK, homo : P := Q, Q := P, Al := "Cantor");*/
