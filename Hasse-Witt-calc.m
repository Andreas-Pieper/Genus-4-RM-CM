AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/Genus-4/magma/spec");
load "Shimura/Shimura_models.m";
prec := 500;
QQ := RationalsExtra(prec);
t0 := QQ!1/2;
Q, E := GetShimuraCurve239(t0);
C := Curve(Proj(Parent(Q)), [Q,E]);

prime_data := [* *];
for p in PrimesInInterval(5, 75) do
  k := GF(p);
  Cp := ChangeRing(C,k);
  M := CartierRepresentation(Cp);
  //lpoly := LPolynomial(Cp);
  //printf "L-poly: %o\n", lpoly;
  ranks := [Rank(M^i) : i in [1..4]];
  //for i := 1 to 4 do
    //print M^i;
    //printf "rank: %o\n", Rank(M^i);
    //print "-------";
  //end for;
  charpoly := CharacteristicPolynomial(M);
  printf "charpoly: %o\n", charpoly;
  //Append(~prime_data, <p, lpoly, ranks, charpoly>);
  Append(~prime_data, <p, ranks, charpoly>);
end for;
