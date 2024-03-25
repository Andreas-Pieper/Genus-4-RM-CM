intrinsic maximal_isotropic(n::RngIntElt, p::RngIntElt) -> SeqEnum
  {}
	set := [1..n] cat [2*n..n by -1];
	ret := [];
	for v in VectorSpace(GF(2), n) do
		pivots := Sort([i+Integers()!v[i]*n: i in [1..n]]);
		nonpivots := Sort([i+Integers()!(1+v[i])*n: i in [1..n]]);
		entries := &cat[[[i,j]: j  in nonpivots|set[pivots[i]] lt set[j]]: i in [1..n]];
		pairs := &cat[[[i1, i2]: i2 in [i1+1..n]]: i1 in [1..n]];
		A := ZeroMatrix(GF(p), #entries, #pairs);
		for nu -> pair in pairs do
			i1 := pair[1];
			i2 := pair[2];
			j1 := (n+pivots[i2]-1) mod (2*n) + 1;
			j2 := (n+pivots[i1]-1) mod (2*n) + 1;
			if set[pivots[i1]] lt set[j1] then
				entry1 := Position(entries, [i1, j1]);
				A[entry1, nu] := 1;
			end if;
                        if set[pivots[i2]] lt set[j2] then
                        	entry2 := Position(entries, [i2, j2]);
				// TODO: Check correctness of signs here
                        	A[entry2, nu] := -1;
			end if;
		end for;
		K := Kernel(A);
		d := #entries;
		for k in K do
			X := ZeroMatrix(GF(p), n, 2*n);
			for i in [1..n] do
				X[i, pivots[i]] := 1;
			end for;
			for nu in [1..d] do
				X[entries[nu][1], entries[nu][2]] := k[nu];
			end for;
			Append(~ret, Image(X));
			end for;
	end for;
	return ret;
end intrinsic;

function LiftFF(c, n)
  return (Integers() ! c)/n;
end function;

function StandardSymplecticMatrix(g)
  id:=IdentityMatrix(Integers(), g);
  zer:=ZeroMatrix(Integers(), g,g);
  J:=BlockMatrix(2,2, [[zer, id], [-id, zer]]);
  return J;
end function;

intrinsic QFromPVFor4(P::ModMatFldElt, V::ModTupFld)-> Any
  {Creates quotient of abelian variety corresponding to P by symplectic subgroup corresponding to V}

  p := Characteristic(BaseRing(V));
  L1 := Lattice(IdentityMatrix(Rationals(), 8));
  M2 := Matrix(Basis(V));
  M2 := Matrix(Rationals(), #Rows(M2), #Rows(Transpose(M2)), [ LiftFF(c, p) : c in Eltseq(M2) ]);
  L2 := Lattice(M2);
  L := L1 + L2;
  T := Matrix(Basis(L));

  E3 := StandardSymplecticMatrix(4);
  T1 := Transpose(T);
  E := Transpose(T1)*E3*T1;

  E0, T2 := FrobeniusFormAlternating(ChangeRing(p*E, Integers()));
  BT := T1*Transpose(T2);

  Q := P*ChangeRing(BT, BaseRing(P));
  //assert IsBigPeriodMatrix(Q);
  return Q;
end intrinsic;

