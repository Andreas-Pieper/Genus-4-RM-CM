import "CMTypes.m" : Exponents;

intrinsic InitializeDataFiles() -> Any 
  {}

  all_filename := "cm-fields-schottky-values.txt";
  jac_filename := "cm-fields-jacobians.txt";
  errors_filename := "cm-fields-errors.txt";
  names := [all_filename, jac_filename, errors_filename];
  for name in names do
    System(Sprintf("touch %o", name));
    Write(name, "Format: label, coeffs, CM type, ideal class rep, polarizing element, complex conjugation, absolute value of Schottky form\n");
  end for;
  return Sprintf("Files %o initialized", Join(names, ", "));
end intrinsic;

intrinsic CheckForJacobians(label::MonStgElt, coeffs_str::MonStgElt, gal_label::MonStgElt : precred := 300, prectheta := 50, preccheck := 150, all_filename := "cm-fields-schottky-values.txt", jac_filename := "cm-fields-jacobians.txt", done_filename := "cm-fields-done.txt", special_filename:="special-fields.txt") -> Any
  
  {}

  exps := Exponents();
  coeffs := eval coeffs_str;
  QQ := RationalsExtra(preccheck);
  CC<I> := QQ`CC;
  eps := CC`epscomp;
  R<x> := PolynomialRing(QQ);
  f := R!coeffs;
  d,k := Explode(Split(gal_label,"T"));
  k := eval k;
  if k in [13, 24] then
    Write(done_filename, Join([label,coeffs_str,gal_label],"|"));
    return Sprintf("Galois group 8T13 or 8T24; skipping");
  end if;
  vals := EnumerationUpToGalois(f : exp := exps[k], prec:=precred,precred:=precred, prectheta:=prectheta);
  Write(done_filename, Join([label,coeffs_str,gal_label],"|"));
  for val in vals do
    K, Phi, aa, xi, invK, sch := Explode(val);
    output_list := [ ];
    Phi_str := Sprint(Phi);
    Phi_str := ReplaceString(Phi_str,"$.1", "I");
    Append(~output_list, Phi_str);
    Append(~output_list, Sprint(Generators(aa)));
    Append(~output_list, Sprint(Eltseq(xi)));
    Append(~output_list, Sprint(Eltseq(invK(K.1))));
    Append(~output_list, Sprint(sch));
    output := StripWhiteSpace(Sprintf(Join([ label, Sprint(coeffs) ] cat output_list, "|")));
    // write all data to all file
    Write(all_filename, output);
    // write Jacobians to Jacobian file
    if val[#val] lt 10^(-preccheck) then
      Write(jac_filename, output);
    end if;
  end for;
  return Sprintf("Data written to %o and %o\n", all_filename, jac_filename);
end intrinsic;
