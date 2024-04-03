function AutoAction()
  R<x> := PolynomialRing(Rationals());
  S := RiemannSurface(x^3-1, 5);
  Chains, IntMat, Sympl := HomologyBasis(S);
  Rams := [Coordinates(ram)[1]: ram in RamificationPoints(S)];
  EdgesList := [edge`EP: edge in Chains`Edges];
  i0:= &meet [Seqset(S): S in EdgesList];
  assert #i0 eq 1;
  i0 := Random(i0);
  assert &and [edge[1] eq i0: edge in EdgesList];
  i1 := EdgesList[1][2];
  i2 := EdgesList[2][2];
  angle := (Rams[i1]-Rams[i0])/(Rams[i2]-Rams[i0]);
  orientation := Sign(Imaginary(angle));
  idx := [i: i in [5..8]|IntMat[1, i] eq orientation];
  assert #idx eq 1;
  idx := idx[1];
  path1 := Matrix(8,1, [-1,0,0,0,0,0,0,0]); 
  path1[idx,1] := 1;
  path2 := Matrix(8,1, [-1,0,0,0,0,0,0,0]);
  X0 := Matrix([[0,0,0,-1], [1,0,0,-1], [0,1,0,-1], [0,0,1,-1]]);
  X := Matrix(BlockDiagMat(<X0, X0>));
  for i in [0..4] do
    Y := HorizontalJoin([X^(i+j)*path1: j in [0..3]] cat [X^j*path2: j in [0..3]]);
    if Transpose(Y) * IntMat * Y eq IntMat then
       ret :=  Y;
    end if;
  end for;
  return Transpose(Sympl)^-1 *  ret * Transpose(Sympl);
end function;
