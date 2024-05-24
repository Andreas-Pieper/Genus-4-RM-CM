AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/cm-calculations-g3/magma/spec");
AttachSpec("~/github/Genus-4/magma/spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");
AttachSpec("~/github/Reconstruction/magma/spec");
AttachSpec("spec");

intrinsic InitializeDataFiles() -> .
  {}

  all_filename := "cm-fields-schottky-values.txt";
  jac_filename := "cm-fields-jacobians.txt";
  for name in [all_filename, jac_filename] do
    System(Sprintf("touch %o", name));
    Write(name, "Format: label, coeffs, index in list of taus, absolute value of Schottky form");
  end for;
  return Sprintf("Files %o and %o\n initialized\n", all_filename, jac_filename);
end intrinsic;

intrinsic CheckForJacobians(label::MonStgElt, coeffs:SeqEnum, filename::MonStgElt : prec := 40) -> .
  {}

  QQ := RationalsExtra(prec);
  CC<I> := QQ`CC;
  eps := CC`epscomp;
  R<x> := PolynomialRing(QQ);
  f := R!coeffs;
  taus := FullEnumerationG4(f);
  for i->tau in taus do
    tau_red := SiegelReduction(tau);
    sch := SchottkyModularForm(tau_red : prec := prec);
    abs := Abs(sch);
    output := Sprintf("%o|%o|%o|%o\n", label, coeffs, i, abs);
    Write(all_filename, output);
    if abs lt 10^(-prec*(3/8)) then
      Write(jac_filename, output);
    end if;
  end for;
  return Sprintf("Data written to %o and %o\n", all_filename, jac_filename);
end intrinsic;

prec := 300;
//8.0.1265625.1
//coeffs := [1, -1, 0, 1, -1, 1, 0, -1, 1];
//8.0.13140625.1
coeffs := [1, 3, 5, 3, 4, -3, 5, -3, 1];
QQ := RationalsExtra(prec);
CC<I> := QQ`CC;
eps := CC`epscomp;
R<x> := PolynomialRing(QQ);
f := R!coeffs;
K<nu> := NumberFieldExtra(f);
_, K0incl, inv := IsCMField(K);
K0, incl := Explode(K0incl);
taus := FullEnumerationG4(f);
tau := taus[1];
// reducing precision
tau0 := tau;
tau := ChangeRing(tau, CC);
print "checking if Jacobian";
tau_red := SiegelReduction(tau);
sch := SchottkyModularForm(tau_red : prec := (prec div 3));
print sch;
assert Abs(sch) lt 10^(-prec/4);

