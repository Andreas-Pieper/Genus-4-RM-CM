/*
  Usage: Make an input file input.txt each of whose rows consists of an LMFDB field label and an LMFDB ideal label separated by a space using GenerateFieldsAndLevels. E.g., 

  8.0.1265625.1|[1, -1, 0, 1, -1, 1, 0, -1, 1]
  8.0.4000000.1|[1, 0, -1, 0, 1, 0, -1, 0, 1]
  8.0.5308416.1|[1, 0, 0, 0, -1, 0, 0, 0, 1]

  Then run the following command in the directory ~/github/Genus-4-RM-CM/CM

  parallel -j 16 --joblog joblog --eta --colsep '|' -a input.txt magma -b label:={1} coeffs:={2} ModFrmHilD/CanonicalRing/CanonicalRingsScript.m
*/

AttachSpec("spec");

try 
  CheckForJacobians(label, coeffs);
  exit 0;
catch e
  WriteStderr(e);
  exit 1;
end try;
