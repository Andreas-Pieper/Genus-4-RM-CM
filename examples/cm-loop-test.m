AttachSpec("spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");
AttachSpec("~/github/CHIMP/CHIMP.spec");
SetVerbose("CMExp",2);
SetDebugOnError(true);

QQ := RationalsExtra();
R<x> := PolynomialRing(QQ);
f := x^8 - x^7 + x^5 - x^4 + x^3 - x + 1;
// f := x^8 + 18*x^6 + 107*x^4 + 246*x^2 + 169; // class number 4
//f := x^8 + 25*x^6 + 215*x^4 + 725*x^2 + 725; // class number 20
time FullEnumerationG4(f);
time EnumerationUpToGalois(f);
