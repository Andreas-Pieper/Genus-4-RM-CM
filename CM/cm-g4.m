intrinsic InitializeDataFiles() -> Any 
  {}

  all_filename := "cm-fields-schottky-values.txt";
  jac_filename := "cm-fields-jacobians.txt";
  for name in [all_filename, jac_filename] do
    System(Sprintf("touch %o", name));
    Write(name, "Format: label, coeffs, index in list of taus, absolute value of Schottky form\n");
  end for;
  return Sprintf("Files %o and %o initialized", all_filename, jac_filename);
end intrinsic;

intrinsic CheckForJacobians(label::MonStgElt, coeffs_str::MonStgElt : prec := 40) -> Any
  {}

  coeffs := eval coeffs_str;
  QQ := RationalsExtra(prec);
  CC<I> := QQ`CC;
  eps := CC`epscomp;
  R<x> := PolynomialRing(QQ);
  f := R!coeffs;
  taus := FullEnumerationG4(f);
  all_filename := "cm-fields-schottky-values.txt";
  jac_filename := "cm-fields-jacobians.txt";
  for i->tau in taus do
    tau0 := tau;
    if not IsSymmetric(tau0) then
      print "tau not symmetric: replacing by (tau + tau^T)/2";
      tau0 := (tau0 + Transpose(tau0))/2;
    end if;
    tau_red := SiegelReduction(tau0);
    if not IsSymmetric(tau_red) then
      print "tau_red not symmetric: replacing by (tau_red + tau_red^T)/2";
      tau_red := (tau_red + Transpose(tau_red))/2;
    end if;
    sch := SchottkyModularForm(tau_red : prec := prec);
    abs := Abs(sch);
    output := StripWhiteSpace(Sprintf("%o|%o|%o|%o", label, coeffs, i, abs));
    Write(all_filename, output);
    if abs lt 10^(-prec*(3/8)) then
      Write(jac_filename, output);
    end if;
  end for;
  return Sprintf("Data written to %o and %o\n", all_filename, jac_filename);
end intrinsic;
