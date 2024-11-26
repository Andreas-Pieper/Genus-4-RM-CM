intrinsic Corestriction(K::Fld, a::RngElt, b::RngElt) -> AlgBas
{Given elements a, b of a cyclic cubic number field defining the quaternion algebra B=(a,b|F), return the Corestriction of B to QQ}

  assert Degree(K) eq 3;
  assert IsCyclic(K);

  auts := Automorphisms(K);
  aconjs := [m(a) : m in auts];
  bconjs := [m(b) : m in auts];
end intrinsic;
