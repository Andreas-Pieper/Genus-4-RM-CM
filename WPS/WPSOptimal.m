intrinsic WPSOptimal(wts::SeqEnum, V::SeqEnum : Precision := 100) -> SeqEnum
{Given a point in weighted projective space with over the real or
complex numbers returns a normalized representative of optimal size} 
     K := Universe(V);
     N := #V;
     assert N eq #wts;
     if Type (K) eq FldCom or Type(K) eq FldRe then 
          _, i0 := Max([Abs((V[i]) ^ (1 / wts[i])) : i in[1..N]]); // find the biggest normalizing factor possible
          V0 := WPSMultiply(wts, V, V[i0] ^ (-1 / wts[i0])); // renormlize by that factor, to see which invariants are actually 0
          V := [Abs(V0[i]) ge 10 ^ (-Precision / 2) select V[i] else 0 : i in[1..N]];
          nonz := [i:i in[1..N] | Abs(V0[i]) ge 10 ^ (-Precision / 2)];
     else 
          nonz := [i:i in[1..N] | V[i] ne 0];
     end if;
     nonz;
     no := #nonz;
     wtsnonz := [wts[nonz[i]]:i in[1..no]];
     Vnonz := [V[nonz[i]]:i in[1..no]];
     gc, si := XGCD(wtsnonz);
     sivec := Vector(Integers(), si);
     La := KernelMatrix(Matrix(Integers(), no, 1, wtsnonz));
     LaK := ChangeRing(La, RealField(Precision));
     if Type (K) eq FldNum then 
          D := DiagonalMatrix([Abs(Log(Abs(RealField(Precision)!Norm(v)))) : v in Vnonz]);
     else
          D := DiagonalMatrix([Abs(Log(Abs(v))) : v in Vnonz]);
     end if;
     G := LaK * D * Transpose(LaK);
     G := G + Transpose(G);
     vec := Vector(K, si) * Transpose(LaK) * (LaK * Transpose(LaK)) ^ -1;
     L := LatticeWithGram(G);
     cl := ClosestVector(L, vec);
     ch := Eltseq(sivec - Vector(ChangeRing(cl, Integers()) * La));
     mul := [K | 1 : i in[1..N]];
     for i in [1..no] do 
          mul[nonz[i]] := V[nonz[i]] ^ ch[i];
     end for;
     lambda := &*mul;
     wts := [i in nonz select wts[i] div gc else 0 : i in[1..N]]; // so that it works also for number fields
     return WPSMultiply(wts, V, 1 / lambda);
end intrinsic;
