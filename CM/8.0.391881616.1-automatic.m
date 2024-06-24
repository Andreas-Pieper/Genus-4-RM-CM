load "full_proc.m";
label := "8.0.391881616.1";
coeffs := [16,-24,8,8,-10,4,2,-3,1];
R<x> := PolynomialRing(QQ);
f := R!coeffs;
K := NumberField(f);
taus := FullEnumerationG4(f);
//InitializeEmbedding(K);
//cm_types := [];
//for tau in taus do

//end for;
