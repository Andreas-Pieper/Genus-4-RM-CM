import "CMTypes.m" : Exponents;

intrinsic InitializeDataFiles() -> Any 
  {}

  all_filename := "cm-fields-schottky-values.txt";
  jac_filename := "cm-fields-jacobians.txt";
  for name in [all_filename, jac_filename] do
    System(Sprintf("touch %o", name));
    Write(name, "Format: label, coeffs, CM type, ideal class rep, polarizing element, complex conjugation, absolute value of Schottky form\n");
  end for;
  return Sprintf("Files %o and %o initialized", all_filename, jac_filename);
end intrinsic;

intrinsic CheckForJacobians(label::MonStgElt, coeffs_str::MonStgElt, gal_label::MonStgElt : precred := 50, prectheta := 50, preccheck := 500) -> Any
  {}

  exps := Exponents();
  coeffs := eval coeffs_str;
  QQ := RationalsExtra(prec);
  CC<I> := QQ`CC;
  eps := CC`epscomp;
  R<x> := PolynomialRing(QQ);
  f := R!coeffs;
  d,k := Split(gal_label,"T");
  vals := EnumerationUpToGalois(f : exp := exps[k]);
  all_filename := "cm-fields-schottky-values.txt";
  jac_filename := "cm-fields-jacobians.txt";
  for val in vals do
    K, Phi, aa, xi, invK, sch := Explode(datum);
    output_list := [* *];
    Phi_str := Sprint(Phi);
    Phi_str := ReplaceString(Phi_str,"$.1", "I");
    Append(~output_list, Phi_str);
    Append(~output_list, Sprint(Generators(aa)));
    Append(~output_list, Sprint(Eltseq(xi)));
    Append(~output_list, Sprint(Eltseq(invK(K.1))));
    Append(~output_list, Sprint(sch));
    output := StripWhiteSpace(Sprintf(Join("|", [* label, coeffs *] cat output_list)));
    // write all data to all file
    Write(all_filename, output);
    // write Jacobians to Jacobian file
    if val[#val] lt 10^(-preccheck/2) then
      Write(jac_filename, output);
    end if;
  end for;
  return Sprintf("Data written to %o and %o\n", all_filename, jac_filename);
end intrinsic;
