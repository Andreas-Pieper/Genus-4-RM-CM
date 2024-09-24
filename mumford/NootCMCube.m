// see Noot - On Mumford's families of abelian varieties, section 3.6

intrinsic CheckNootCondition(K::FldNum, CMType::SetEnum) -> BoolElt
  {Given a (CM) number field of degree 8 and a CM type, check if Noot s condition on the action of the Galois group on the CM type is satisfied.}

  prec := 30; //FIXME
  Gal, roots, data := GaloisGroup(K : Type := "Complex");
  G := TransitiveGroup(8,24);
  b, g := IsConjugate(Sym(8),G,Gal);
  if not b then return false; end if;

  // match roots with CMType
  inds := [];
  for r in CMType do
    m, i := Min([Abs(r-r2) : r2 in roots]);
    assert m lt 10^(-prec/3);
    Append(~inds,i);
  end for;
  inds_act := {i^(g^-1) : i in inds}; // g or g^-1???
  faces := [
    {1,7,2,4}, {1,5,8,4}, {1,5,3,7},
    {6,2,7,3}, {6,8,4,2}, {6,8,5,3}
  ];
  return inds_act in faces;
end intrinsic;