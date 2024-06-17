/*AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/cm-calculations-g3/magma/spec");
AttachSpec("~/github/Genus-4/magma/spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");
AttachSpec("~/github/Reconstruction/magma/spec");

SetVerbose("Reconstruction",true);
SetVerbose("Theta",true);
SetDebugOnError(true);*/
Attach("~/Genus-4-RM-CM/WPS/WPSOptimal");

function FindCurveHyperelliptic(taus, F)
  for tau in taus do
    rosens := ReconstructCurveG4(tau);
    CC := Parent(rosens[1]);
    // recognize Rosenhain invariants

    _<x,y> := PolynomialRing(F, 2);
    f := y*x*(x-y)*(&*[x-r*y : r in rosens]);
  end for;
  inv, wgt := InvariantsGenus4Curves(f);// : normalize := true);
  vprint CM: "Computing an optimized normalized model...";
  inv_opti := WPSOptimal(wgt, inv);// can try WPSNormalize if this takes too long
//Andreas: Still need to figure out the best way to normalize in a WPS here.
  Append(~invs, inv_opti);
  vprint CM: "Computing better model...";        
  t, lambda := IsSquare(inv[1]/inv_opti[1]);
  if not t then
      // TODO: introduce quadratic extension
  end if;

  C_24 := Transvectant(f, f, 8);
  C_28 := Transvectant(f, f, 6);
  c_3 := 1/lambda^3*Transvectant(C_28, f, 8);
  c_5_1 := 1/lambda^2*Transvectant(C_24, c_3, 2);
  c_5_2 := 1/lambda^5*Transvectant(C_28, Transvectant(C_24, f, 4), 6);

  C := [c_3, c_5_1, c_5_2];

  _<[X]> := PolynomialRing(Parent(inv_opti[1]), 3);
  Q := &+[Rationals()!(Evaluate(Transvectant(C[i], C[j], 2), [0,0]))*X[i]*X[j] : i in [1..3], j in [1..3]];
  E := &+[Rationals()!(1/lambda*Evaluate(Transvectant(C[i1]*C[i2]*C[i3]*C[i4]*C[i5], f, 10), [0,0]))*X[i1]*X[i2]*X[i3]*X[i4]*X[i5] : i1 in [1..3], i2 in [1..3], i3 in [1..3], i4 in [1..3], i5 in [1..3]];

  Con := Conic(ProjectiveSpace(Parent(Q)), Q);

  try 
      p := Parametrization(Con);
      f := E @@ p;

      t, f_bis := IsCoercible(PolynomialRing(Rationals(), 4), f);
      if not t then
          return f;
      else
          return MinRedBinaryForm(f_bis);
      end if;

  catch e 
      print "No parametrization found, try parametrizing Q using ConicParametrization (https://github.com/JRSijsling/hyperelliptic/blob/main/magma/toolbox/diophantine.m)";
      print "Mestre conic Q and f_5 returned";
      return Q, E;
  end try;

end function;

// Input: taus, a Galois orbit of period matrices
//        F, intersection of field of moduli with the reflex field
function FindCurve(taus, F)
  invs := [];
  for tau in taus do
    tau0 := (tau + Transpose(tau))/2;
    tau_red := SiegelReduction(tau0);
    tau_red := (tau_red + Transpose(tau_red))/2;
    
    vprint CM: "Checking if Jacobian...";
    sch := SchottkyModularForm(tau_red : prec := (prec div 3));
    if sch gt 10^(-prec/4) then
        print "Not a jacobian.";
        return Abs(sch);
    end if;
    print "It is probably a jacobian.";
    printf "Norm of the Schottky modular form: %o\n", Abs(sch);
    Eqs := ReconstructCurveG4(tau_red);

    if #Eqs eq 7 then // hyperelliptic case
        vprint CM: "Hyperelliptic case";
        return FindCurveHyperelliptic(); // TODO

    elif #Eqs eq 2 then
        vprint CM: "Non-hyperelliptic curve";
        Q, E := Explode(Eqs);
        CC4<X,Y,Z,T> := Parent(Q);
        CC<I> := BaseRing(CC4);
        inv, wgt := InvariantsGenus4Curves(Q, E : normalize := true);
	//Andreas: When the rank is 3, then normalizing by putting the first coefficient to 1 is dangerous. In the rank 4 case you do not have this issue.
        if #inv eq 65 then // checks if it is actually rank 3 even if the precision is wrong
            inv_rk3 := InvariantsGenus4Curves(X*T-Y*Z, (Y-Z)^3 : normalize := true);
            if Max([Abs(CC!inv_rk3[i]-inv[i]) : i in [1..65]]) le 10^(-2) then // it is actually of rank 3
                while #inv eq 65 do
                    CC<I> := ComplexField(Floor(3/4*Precision(CC)));
                    inv, wgt := InvariantsGenus4Curves(ChangeRing(Q, CC), ChangeRing(E, CC));
                end while;
            end if;
        end if;
        Append(~invs, inv);
  end for;
  assert #Seqset([#el : el in invs]) eq 1;
  invs_al := RecognizeInvariantsOrbits(invs, wgt, F); // see WPS
  curves := [ReconstructionGenus4(el) : el in invs_al];
  return curves;
end function;
