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

function ProcessNumberString(s,prec)
  if "*" in s then
    parts := Split(s, "*");
    head, tail := Explode(parts);
    assert (tail eq "i") or (tail eq "I");
    head *:= Sprintf("p%o", prec);
    //val := head*"*"*tail;
    val := head*"*"*"I";
  else
    val := s*Sprintf("p%o", prec);
  end if;
  return val;
end function;

function ProcessComplexNumberString(Phistr,prec)
  CC<I> := ComplexField(prec);
  if Phistr[1] eq "-" then
    re_sign := CC!-1;
    t := Phistr[2..#Phistr];
  else
    re_sign := CC!1;
    t := Phistr;
  end if;
  if (not "-" in t) and (not "+" in t) then // pure real or imag
    return re_sign*(eval ProcessNumberString(t,prec));
  end if;
  if "-" in t then
    im_sign := CC!-1;
    spl_char := "-";
  else
    im_sign := CC!1;
    spl_char := "+";
  end if;
  parts := Split(t, spl_char);
  re := ProcessNumberString(parts[1],prec);
  im := ProcessNumberString(parts[2],prec);
  return re_sign*(eval re) + im_sign*(eval im);
end function;

// remember to make all complex numbers use I (not i)
intrinsic JacobianFileToArray( : filename:="../cm-fields-jacobians.txt") -> Assoc
  {Load the contents of filename ("cm-fields-jacobians.txt" by default) into an associative array. File should be of the form output by WriteArrayToFile.}

  A := AssociativeArray();
  R<x> := PolynomialRing(Rationals());
  file := Open(filename, "r");
  line := Gets(file); // skip header
  line := Gets(file);
  eof := false;
  cnt := 0;
  while not eof do
  if cnt mod 100 eq 0 then
    printf "read in %o lines\n", cnt;
  end if;
    line := Gets(file);
    if IsEof(line) then
      eof := true;
      break;
    end if;
    spl := Split(line, "|");
    if #spl ne 7 then
      print "Error! Wrong number of data fields.";
    end if;
    // K, Phi, aa, xi, invK, sch := Explode(val);
    label, Kpolystr, Phistr, aastr, xistr, invKstr, schstr := Explode(spl);
    // process CM field
    Kpoly := eval Kpolystr;
    K<nu> := NumberField(Polynomial(Kpoly));
    // process CM type
    prec := 300;
    Phistr := Phistr[2..#Phistr-1]; // get rid of brackets
    Phiparts := Split(Phistr,",");
    Phi := {ProcessComplexNumberString(part,prec) : part in Phiparts};
    // process ideal
    aagens := eval aastr;
    // gens written wrt to monomial or integral basis?
    OK := Integers(K);
    aa := ideal< OK | [OK!el : el in aagens] >;
    // process polarizing element xi
    xi := K!(eval xistr);
    // process involution
    invKelt := K!(eval invKstr);
    // process Schottky value
    sch := eval schstr;
    invK := hom< K -> K | invKelt >;
    data := <K, Phi, aa, xi, invK, sch>;
    if IsDefined(A, label) then // if defined, then append
      Append(~A[label], data);
    else
      A[label] := [* data *];
    end if;
    cnt +:= 1;
  end while;
  return A;
end intrinsic;

/*
intrinsic CompleteErrorFile(ErrorFile, CoeffsFile)
  {}

  while not eof do
  if cnt mod 100 eq 0 then
    printf "read in %o lines\n", cnt;
  end if;
    line := Gets(ErrorFile);
    if IsEof(line) then
      eof := true;
      break;
    end if;
    spl := Split(line, "|");
    if #spl ne 7 then
      print "Error! Wrong number of data fields.";
    end if;
    // K, Phi, aa, xi, invK, sch := Explode(val);
    label, Kpolystr, Phistr, aastr, xistr, invKstr, schstr := Explode(spl);
    // process CM field

end intrinsic;
*/

intrinsic GetMissingFields( : AllFieldsPath:="../CM-fields-deg-8.txt", ComputedFieldsPath:="../cm-fields-schottky-values.txt", MissingPath:="../missing_fields.txt") -> Any
  {}

  all_labels := [];
  all_data := [];
  file := Open(AllFieldsPath, "r");
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
    Append(~all_labels, <spl[1], spl[2]>);
    Append(~all_data, <spl[1], spl[2], spl[3]>);
  end while;
  all_labels := Set(all_labels);

  computed_labels := [];
  file := Open(ComputedFieldsPath, "r");
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
    Append(~computed_labels, <spl[1], spl[2]>);
  end while;
  computed_labels := Set(computed_labels);
  missing_labels := all_labels diff computed_labels;
  missing_data := [el : el in all_data | <el[1], el[2]> in missing_labels];

  printf "%o missing labels found\n", #missing_data;

  // write to file
  for el in missing_data do
    s := Join([el[1], el[2], el[3]], "|");
    Write(MissingPath, s);
  end for;
  return Sprintf("Data written to %o", MissingPath);
end intrinsic;
