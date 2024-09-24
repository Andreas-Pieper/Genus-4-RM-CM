Attach("~/github/Genus-4-RM-CM/CM/findg4.m");
AttachSpec("~/github/CHIMP/CHIMP.spec");
Attach("~/github/Genus-4-RM-CM/flint-stuff/FlintWrapper.m");
Attach("~/github/Genus-4-RM-CM/flint-stuff/schottky-fast-theta.m");
Attach("~/github/Genus-4-RM-CM/mumford/NootCMCube.m");
AttachSpec("~/github/cm-calculations/magma/spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");

load "~/CM/new_fields.m";

SetDebugOnError(true);
load "~/github/Genus-4-RM-CM/CM/full_proc.m";


function OliveToIgusa(I)
    J1 := -15*I[1];
    J2 := 45/8*(3*I[1]^2 - 25/2*I[2]);
    J3 := 15/8*(9/2*I[1]^3 - 25*I[1]*I[2] - 375/4*I[3]);
    J4 := (J1*J3-J2^2)/4;
    J5 := 81/16*(-3*I[1]^5 + 7625/256*I[1]^3*I[2] + 13125/256*I[1]^2*I[3] - 9375/128*I[1]*I[2]^2 - 28125/128*I[2]*I[3] - 28125/256*I[4]);
    return [J1,J2,J3,J4,J5];
end function;

for d in data[#data-25..#data] do
  d[1];
  //Write(d[1] cat ".m", "" : Overwrite := true);
  //if d[5] eq [] then
    coeffs := d[2];
    coeffs;
    taus := FullEnumerationG4(PolynomialRing(Rationals())!coeffs : prec := 2000);
    tau := [t[1] : t in taus];
    phi := [t[2] : t in taus];
    for i in [1..#taus] do
      if CheckNootCondition(NumberField(PolynomialRing(Rationals())!coeffs), phi[i]) then
        "nay!", i;
        tau_i := (tau[i]+Transpose(tau[i]))/2;
        tau_red := SiegelReduction(tau_i);
        tau_red := (tau_red+Transpose(tau_red))/2;

        abs := Abs(SchottkyModularFormFlint(ChangeRing(tau_red, ComplexField(100))));
        abs;
       if abs le 10^(-90) then

          Q, E := FindCurve(tau_red);
        
          if E eq 0 then
            if Q ne 0 then
              K := BaseRing(Parent(Q));
              //Write(d[1] cat ".m", "K :=" cat Sprint(K) cat "\n" cat Sprint(Q) cat ";\n");
            end if;
          else
            "yay!";
            Q;
            E;
            K := BaseRing(Parent(Q));
            if Type(K) eq FldRat then
              Q, E := MinimizeG4(Q, E);
            end if;
            //Write(d[1] cat ".m", "K := " cat Sprint(K) cat "\nQ := " cat Sprint(Q) cat ";\nE := " cat Sprint(E) cat ";\n");
          end if;
        end if;
      else "no :(";
      end if;
    end for;
  //end if;
end for;
