Attach("~/github/Genus-4-RM-CM/CM/findg4.m");
AttachSpec("~/github/CHIMP/CHIMP.spec");
Attach("~/github/Genus-4-RM-CM/flint-stuff/FlintWrapper.m");
Attach("~/github/Genus-4-RM-CM/flint-stuff/schottky-fast-theta.m");
Attach("~/github/Genus-4-RM-CM/mumford/NootCMCube.m");
AttachSpec("~/github/cm-calculations/magma/spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");
load "~/lmfdb_8t13.m";
load "~/github/Genus-4-RM-CM/CM/minimize.m";

SetDebugOnError(true);
load "~/github/Genus-4-RM-CM/CM/full_proc.m";

/*
d := data[10];
taus := FullEnumerationG4(R!coeffs : prec := 2000);
tau := [t[1] : t in taus];
tau_red := SiegelReduction(tau[1]);
inv := FindCurve(tau_red);

AttachSpec("~/github/Reconstruction/magma/spec");

QQ1 := RationalsExtra(2000);
inv1 := ChangeUniverse(WPSNormalize([1,2,3,5],inv[1..4]), ComplexField(2000));
inv2 := ChangeUniverse(WPSNormalize([2,3], inv[6..7]), ComplexField(2000));
wgt := [6,12,18,30,45,4,6,8,10,10,9,13,14,12,12,15,14,14,15,15,17,16,16,20,19,19,18,16,18,18,17,17,21,19,22,22,23,21,20,20,21,25,26,24,21,23,23,24,25,29,25,28,27,32,29,31,35,33,37,41];

inv0 := ChangeUniverse(inv, ComplexField(2000));
inv0 := WPSMultiply(wgt, inv0, (inv0[6]/inv0[1])^(1/2));
L, inv_alg := NumberFieldExtra([(wgt[i] mod 2) eq 0 select inv else 0 : i->inv in inv0], QQ1);

SetVerbose("Reconstruction", 1);

Q, E := ReconstructionGenus4(inv_alg);
Q1, E1 := MinimizeQuadric(Q,E);

inv1 := ChangeUniverse(WPSNormalize([1,2,3,5],inv[1..4]), ComplexField(2000));
CC<I> := Universe(inv0);

inv2 := inv0[6..7];
inv2 := WPSMultiply(wgt[6..7], inv2, Exp(4*I*Pi(CC)/3));
L, inv_alg := NumberFieldExtra(inv1,QQ1);
L2, inv_alg2 := NumberFieldExtra(inv2,QQ1);

f := HyperellipticPolynomials(HyperellipticCurveFromIgusaInvariants(OliveToIgusa(inv_alg) : minimize := true));

ReconstructionGenus4(inv_alg cat [0 : i in [1..56]]);

Polredabs(MinimalPolynomial(inv1[4], 4));


_<x> := PolynomialRing(RationalsExtra(2000));
L<nu> := NumberFieldExtra(x^3-3*x-1 : prec := 2000);
AlgebraizeElementsExtra(inv1,L);


*/

//Write("curves_class.txt", "" : Overwrite := true);
for d in data[2686..#data] do
  if d[5] ne [] then
    R<x> := PolynomialRing(Rationals());
    d[1];
    coeffs := d[2];
    R!coeffs;
    taus := FullEnumerationG4(R!coeffs : prec := 2000);
    tau := [t[1] : t in taus];
    for i in [1..#tau] do 
      tau_red := SiegelReduction(tau[i]);
      Q, E, djn := FindCurve(tau_red);
      if E eq 0 then
        if Q ne 0 then
          S := Parent(Q);
          try 
            K := BaseRing(S);
          catch e
            K := Rationals();
          end try;
          Write("curves_class_8t13.txt", Sprint(d[1]) cat "|" cat Sprint(djn) cat "|" cat Sprint(R!DefiningPolynomial(K)) cat "|" cat Sprint(Q) cat "|");
        else
          Write("curves_class_8t13.txt", Sprint(d[1]) cat "|Not able to recognize invariants");
        end if;
      elif E ne 1 then
        "yay!";
        S := Parent(Q);
        try 
          K := BaseRing(S);
        catch e
          K := Rationals();
        end try;
        //if Type(K) eq FldRat then
          //Q, E := MinimizeG4(Q, E);
        //else
          //"Minimizing...";
          //Q, E := MinimizeQuadric(Q, E);
        //end if;
        //Q;
        //E;
        Write("curves_class_8t13.txt", Sprint(d[1]) cat "|" cat Sprint(djn) cat "|" cat Sprint(R!DefiningPolynomial(K)) cat "|" cat Sprint(Q) cat "|" cat Sprint(E) cat "");
      end if;
    end for;
  end if;
end for;


Eqs := ChangeUniverse(Eqs, ComplexFieldExtra(20000));
L, ros := NumberFieldExtra(Eqs, RationalsExtra(20000));


S<x> := PolynomialRing(Universe(Eqs));
f := x*(x-1)*&*[x-r : r in Eqs];
L, ros := NumberFieldExtra(Coefficients(f), RationalsExtra(20000));
f := Evaluate(f, x-1/9*MonomialCoefficient(f, x^8));
u := MonomialCoefficient(f, x^7)^(1/2);
f1 := 1/u^9*Evaluate(f, u*x);
L, ros := NumberFieldExtra(Coefficients(f1), RationalsExtra(20000));


R := ProductProjectiveSpace(RationalsExtra(60), [1, 1]); 
S<x,y,u,v> := CoordinateRing(R);
M := [m : m in MonomialsOfDegree(S, 6) | Degree(m, x)+Degree(m, y) eq Degree(m, u)+Degree(m, v)];
f := &+[Random(5)*m : m in M];

C := Curve(R, f);
time P := PeriodMatrix(C);
