AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");
AttachSpec("~/github/Genus-4-RM-CM/CM/spec");
//AttachSpec("~/github/Genus-4/magma/spec");
//AttachSpec("~/github/Reconstruction/magma/spec");

import "CMTypes.m" : Exponents;

SetDebugOnError(true);
SetClassGroupBounds("GRH");
SetVerbose("CMExp",true);

res_djn3 := [];
res_djn5 := [];
res_hyp := [];
res_others := [];
R<x> := PolynomialRing(Rationals());

//filename := "../cm-fields-jacobians.txt";
filename := "../missing-fields-jacobians.txt";
file := Open(filename, "r");
line := Gets(file); // skip header
line := Gets(file);
eof := false;
cnt := 0;
while not eof do
if cnt mod 100 eq 0 then
  printf "read in %o lines\n", cnt;
end if;
  line := Gets(file);
  if IsEof(line) then
    eof := true;
    break;
  end if;
  spl := Split(line, "|");
  if #spl ne 7 then
    print "Error! Wrong number of data fields.";
  end if;
  // K, Phi, aa, xi, invK, sch := Explode(val);
  label, Kpolystr, Phistr, aastr, xistr, invKstr, schstr := Explode(spl);
  // process CM field
  Kpoly := eval Kpolystr;
  L<nu> := NumberField(Polynomial(Kpoly));
    try
        _ := Sqrt(L!-1);
        Append(~res_hyp, spl);
        "hyp", #res_hyp;
    catch e
        try 
            _ := RootOfUnity(3, L);
            Append(~res_djn3, spl);
            "djn3", #res_djn3;
        catch e
            try 
                _ := RootOfUnity(5, L);
                Append(~res_djn5, spl);
                "djn5", #res_djn5;
            catch e 
                Append(~res_others, spl);
                "other", #res_others;
            end try;
        end try;
    end try; 
    cnt +:= 1;
end while;

discs := [];
for el in res_others do
  label := el[1];
  spl := Split(label, ".");
  disc := eval spl[3];
  Append(~discs, disc);
end for;
ParallelSort(~discs, ~res_others);

print "Checking possible interesting examples at higher precision";
printf "%o examples to check", #res_others;
interesting := [];
//for j->spl in [res_djn3[1]] do 
for i->spl in res_others do 
    label, Kpolystr, Phistr, aastr, xistr, invKstr, schstr := Explode(spl);
    if j mod 100 eq 0 then
      printf "%o examples checked", j;
    end if;
    print label;
    // process CM field
    prec := 5000;
    Kpoly := eval Kpolystr;
    K<nu> := NumberFieldExtra(Polynomial(Kpoly) : prec := prec);
    inf_places := InfinitePlacesExtra(K); 
    CC<i> := K`CC;
    // process CM type
    Phistr := Phistr[2..#Phistr-1]; // get rid of brackets
    Phiparts := Split(Phistr,",");
    Phi := {eval part : part in Phiparts};
    //Phi_new := {phi : phi in inf_places | #[phi0 : phi0 in Phi | Abs(phi-phi0) le 10^(-25)] eq 1};
    //assert #Phi_new eq 4;
    // process ideal
    aagens := eval aastr;
    // gens written wrt to monomial or integral basis?
    OK := Integers(K);
    aa := ideal< OK | [OK!el : el in aagens] >;
    // process polarizing element xi
    xi := K!(eval xistr);
    // process involution
    invKelt := K!(eval invKstr);
    // process Schottky value
    sch := eval schstr;
    invK := hom< K -> K | invKelt >;
    data := [* K, Kpoly, aa, xi, invK, sch *];
    RealField(20)!ComputeSchottky(K, Phi, aa, xi, invK, sch);
    sch_high := ComputeSchottky(K, Phi, aa, xi, invK, sch);
    if Abs(sch_high) lt 10^(-prec/2) then
      print "---------------------------";
      print "Interesting example found!";
      print label;
      print "---------------------------";
      Append(~interesting,spl);
    end if;
end for;
