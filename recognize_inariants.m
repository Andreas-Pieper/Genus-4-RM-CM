function recognize_invariants(invs, wts, K)
OK := Integers(K);
dK := Degree(K);
basis := Basis(OK);
v := InfinitePlaces(K)[1];
basis_CC := [CC!Evaluate(el, v : Precision := prec) : el in basis];
Cl, m:= ClassGroup(K);
eles :=  ElementaryDivisors(Cl);

ps := [];
es := [];
invs_K := [];
Lat := ZeroMatrix(Integers(), 0, #eles);
for i->inv in invs do
  printf "ps = %o, es = %o\n", ps, es;
  opt := Vector([es[k] * wts[i] : k in [1..#ps]]);
  optLat := opt*Lat^-1;
  L := LatticeWithGram(GLat);
  closLat := ClosestVector(L, optLat);
  clos := Eltseq(opt*Lat);
  scaleId := &*[ps[j]^clos[j]: j in [1..#ps]];
  _,   scale := IsPrincipal(scaleId);
  inv_scaled := inv * Evaluate(scale, v);
  rel := IntegerRelation(basis_CC cat [inv_scaled], 10^20);
  printf "Integer relation found %o\n", rel;
  d := rel[#rel];
  print "Factoring";
  fact := Factorization(d);
  psd := [[p, j]: j->p in ps| not IsDivisibleBy(d, p)];
  for pind in psd do
      p, j := Explode(pind);
      exp := Valuation(&+[rel[i], p);
      e_new := (Round(es[j] * wts[i]) - exp)/ wts[i];
      es[j] := Max([es[j], e_new]);
  end for;
  has_new := false;
  for pair in fact do
    p, e := Explode(pair);
    e_new := e/wts[i];
    if p in ps then
      j := Index(ps, p);
      es[j] := Round(es[j] * wts[i]) / wts[i] + e_new;
    else
      has_new := true;
      Append(~ps,p);
      Append(~es, e_new);
    end if;
  end for;
  inv_K := (-1/d)*&+[rel[i]*basis[i] : i in [1..#basis]];
  inv_K /:= scale;
  printf "Invariant recognized: %o\n", inv_K;
  Append(~invs_K, inv_K);
  if has_new then
      M := VerticalJoin(Matrix([Eltseq(p@@m): p in ps]), DiagonalMatrix( eles));
      KerM := KernelMatrix(M);
      Lat := Submatrix(KerM, 1,1, #ps, #ps);
      G := DiagonalMatrix([Log(Norm(p)): p in ps]);
      GLat := Transpose(Lat)*G*Lat;
  end if;
end for;



end function;
