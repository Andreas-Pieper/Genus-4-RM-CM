//AttachSpec("~/github/CHIMP/CHIMP.spec"); // change to CHIMP spec path
// requires CHIMP

intrinsic PolynomialToString(F::RngMPolElt, i::RngIntElt) -> MonStgElt
  {}
  s := Sprint(F);
  if "$.1" in s then
    s := ReplaceString(s, "$.1", "K.1");
  end if;
  return s;
end intrinsic;

intrinsic WriteArrayToFile(A::Assoc : filename:="CM-curves-class-1.txt") -> Any
  {}

  Write(filename, MakeHeader()*"\n");
  R<x> := PolynomialRing(Rationals());
  for k in Keys(A) do
    rec := A[k];
    assert #rec[1] eq #rec[2];
    assert #rec[1] eq #rec[3];
    if #rec[1] eq 0 then
      Write(filename, Sprintf("%o|||", k));
    else
      for i := 1 to #rec[1] do
        Ki := rec[1][i];
        Qi := rec[2][i];
        Ei := rec[3][i];
        line := strip(Join([k, Sprint(R!DefiningPolynomial(Ki)), PolynomialToString(Qi,i), PolynomialToString(Ei,i)], "|"));
        Write(filename, line);
      end for;
    end if;
  end for;
  return Sprintf("Data written to file %o", filename);
end intrinsic;

intrinsic MakeHeader() -> MonStgElt
  {}
  return Join(["CM_field_label", "base_field", "quadric", "cubic"], "|");
end intrinsic;

intrinsic FileToArray( : filename:="CM-curves-class-1.txt") -> Assoc
  {Load the contents of filename ("CM-curves-class-1.txt" by default) into an associative array. File should be of the form output by WriteArrayToFile.}

  A := AssociativeArray();
  R<x> := PolynomialRing(Rationals());
  file := Open(filename, "r");
  line := Gets(file); // skip header
  line := Gets(file);
  eof := false;
  while not eof do
    line := Gets(file);
    if IsEof(line) then
      eof := true;
      break;
    end if;
    spl := Split(line, "|");
    if #spl eq 1 then // no curves
      label := spl[1];
      A[label] := [* [], [*  *], [*  *] *];
    else
      label, Kpolystr, Qstr, Estr := Explode(spl);
      Kpoly := eval Kpolystr;
      K := NumberField(Kpoly);
      S<X,Y,Z,T> := PolynomialRing(K,4);
      Q := eval Qstr;
      E := eval Estr;
      if IsDefined(A, label) then // if defined, then append
        Append(~A[label][1], K);
        Append(~A[label][2], Q);
        Append(~A[label][3], E);
      else
        A[label] := [* [K], [* Q *], [* E *] *];
      end if;
    end if;
  end while;
  return A;
end intrinsic;
