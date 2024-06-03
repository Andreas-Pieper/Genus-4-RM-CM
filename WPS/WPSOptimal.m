intrinsic WPSOptimal(wts :: SeqEnum, V::SeqEnum : Precision := 100) -> SeqEnum
{Given a point in weighted projective space with over the real or complex numbers returns a normalized representative of optimal size}
  K := Universe(V);
  N := #V;
  assert N eq #wts;
  nonz := [i: i in [1..N]| Abs(V[i]) ge 10^(-Precision/2)];
  no := #nonz;
  wtsnonz := [wts[nonz[i]]: i in [1..no]];
  Vnonz := [V[nonz[i]]: i in [1..no]];
  gc, si := XGCD(wtsnonz);
  sivec := Vector(Integers(), si);
  La := KernelMatrix(Matrix(Integers(), no, 1, wtsnonz));
  LaK := ChangeRing(La, K);
  D := DiagonalMatrix([Abs(Log(Abs(v))): v in Vnonz]);
  G := LaK * D * Transpose(LaK);
  G := G+Transpose(G);
  vec := Vector(K, si) * Transpose(LaK) * (LaK * Transpose(LaK))^-1;
  L := LatticeWithGram(G);
  cl := ClosestVector(L, vec);
  ch := Eltseq(sivec - Vector(ChangeRing(cl, Integers()) * La));
  mul := [K|1: i in [1..N]];
  for i in [1..no] do
       mul[nonz[i]] := V[nonz[i]]^ch[i];
  end for;
  lambdagc := &*mul;
  if Type(K) eq FldCom then
       lambda := lambdagc^(1/gc);
  elif Type(K) eq FldRe then
       lambda := Abs(lambdagc)^(1/gc);
  end if;
  return WPSMultiply(wts, V, 1/lambda);
end intrinsic;
