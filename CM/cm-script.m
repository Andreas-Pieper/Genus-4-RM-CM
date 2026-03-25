/*
  Usage: Make an input file CM-fields-deg-8.txt each of whose rows consists of an LMFDB field label, the coefficients of the defining polynomial, and the label of its Galois group. E.g., 

  8.0.1265625.1|[1,-1,0,1,-1,1,0,-1,1]|8T2
  8.0.4000000.1|[1,0,-1,0,1,0,-1,0,1]|8T2
  8.0.5308416.1|[1,0,0,0,-1,0,0,0,1]|8T3

  Then run the following command in the directory ~/github/Genus-4-RM-CM/CM

parallel -j 25 --joblog joblog --eta --colsep '\|' -a CM-fields-deg-8.txt magma -b label:={1} coeffs:={2} gal_label:={3} CM/cm-script.m

parallel -j 25 --joblog joblog --eta --colsep '\|' -a CM-fields-deg-8.txt --dry-run magma -b label:={1} coeffs:={2} gal_label:={3} CM/cm-script.m

Use --dry-run to see what the commands would be run without actually executing
*/

AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");
AttachSpec("~/github/Genus-4-RM-CM/CM/spec");
//AttachSpec("~/github/Genus-4/magma/spec");
//AttachSpec("~/github/Reconstruction/magma/spec");

try 
  print coeffs;
  CheckForJacobians(label,coeffs,gal_label);
  exit 0;
catch e
  //WriteStderr(e);
  E := Open("cm-fields-errors.txt", "a");
  s := Join([Sprint(el) : el in [* label, coeffs, gal_label *]], "|");
  s := ReplaceAll(" ", "", s);
  Write(E, s);
  //Write(E, e);
  exit 1;
end try;
