function IsLabel(line)
  LABEL_RE := "^8.0.[0-9]+.[0-9]+";
  return Regexp(LABEL_RE, line);
end function;

function FieldPolynomial(line)
  R<x> := PolynomialRing(Rationals());
  if "K := Rational Field" in line then // QQ
    return x-1;
  elif "K := Number Field with" in line then
    i1 := Index(line, "polynomial") + 10; // cuz I don't know how to regex in Magma
    i2 := Index(line, "over") - 2;
    return eval line[i1..i2];
  else
    error "Unknown field type";
  end if;
end function;

function CurvePolynomial(line, K)
  _<X,Y,Z,T> := PolynomialRing(K,4);
  s := Split(line, "=")[2];
  s := Split(s, ";")[1];
  //s := s[1..#s-1]; // get rid of semicolon, for some reason
  return eval s;
end function;

procedure ProcessLine(line, ~Ks, ~Qs, ~Es)
  if IsLabel(line) then // label
    printf "Found label %o\n", line;
  elif "K := " in line then // field
    K := NumberField(FieldPolynomial(line));
    Append(~Ks, K);
  elif "Q := " in line then // quadric 
    K := Ks[#Ks]; // should be the most recently added field
    Q := CurvePolynomial(line, K);
    Append(~Qs, Q);
  elif "E := " in line then // cubic 
    K := Ks[#Ks]; // should be the most recently added field
    E := CurvePolynomial(line, K);
    Append(~Es, E);
  elif "Not able to recognize" in line then
    print line;
  else
    error "Unknown line type";
  end if;
end procedure;

// top-level function
function ProcessFile(: path := "curves_class1.m");
  file := Open(path, "r");
  A := AssociativeArray();
  eof := false;
  Ks := false; Qs := false; Es := false;
  while not eof do
    line := Gets(file);
    if IsEof(line) then
      eof := true;
      break;
    end if;
    if #line eq 1 then // skip blank lines
      continue;
    end if;
    label_bool, new_label := IsLabel(line);
    if label_bool then
      if assigned label then
        if not Ks cmpeq false then
          A[label] := [* Ks, Qs, Es *];
        end if;
      end if;
      label := new_label;
      Ks := [];
      Qs := [* *];
      Es := [* *];
    end if;
    ProcessLine(line, ~Ks, ~Qs, ~Es);
  end while;
  if not IsDefined(A,label) then // think I have to do the last one separately
    A[label] := [* Ks, Qs, Es *];
  end if;
  return A;
end function;

/*
function ProcessFile(: path := "curves_class1.m");
  s := Read(path);
  lines := Split(s,"\n");
  lines := [el : el in lines | #el gt 1]; // no blank lines
  A := AssociativeArray();
  for line in lines do
    label_bool, label := IsLabel(line);
    if #line eq 1 then // blank line, skip
      continue;
    end if;

  end for;
end function;
*/


