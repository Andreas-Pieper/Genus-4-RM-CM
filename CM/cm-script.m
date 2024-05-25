/*
  Usage: Make an input file cm-fields-deg-8-h-1.txt each of whose rows consists of an LMFDB field label and an LMFDB ideal label separated by a space using GenerateFieldsAndLevels. E.g., 

  8.0.1265625.1|[1, -1, 0, 1, -1, 1, 0, -1, 1]
  8.0.4000000.1|[1, 0, -1, 0, 1, 0, -1, 0, 1]
  8.0.5308416.1|[1, 0, 0, 0, -1, 0, 0, 0, 1]

  Then run the following command in the directory ~/github/Genus-4-RM-CM/CM

parallel -j 16 --joblog joblog --eta --colsep '|' -a cm-fields-deg-8-h-1.txt magma -b label:={1} coeffs:={2} cm-script.m
*/

AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/cm-calculations-g3/magma/spec");
AttachSpec("~/github/Genus-4/magma/spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");
AttachSpec("~/github/Reconstruction/magma/spec");
AttachSpec("spec");

try 
  CheckForJacobians(label, coeffs);
  exit 0;
catch e
  WriteStderr(e);
  exit 1;
end try;
