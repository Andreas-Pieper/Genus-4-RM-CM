AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/cm-calculations/magma/spec");
AttachSpec("~/github/Genus-4/magma/spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");
AttachSpec("~/github/Reconstruction/magma/spec");
Attach("~/github/Genus-4-RM-CM/flint-stuff/schottky-fast-theta.m");

SetVerbose("Reconstruction",true);
SetVerbose("Theta",true);
SetVerbose("Genus4",true);
SetDebugOnError(true);
//Attach("~/Genus-4-RM-CM/WPS/WPSOptimal");

import "~/Genus-4-RM-CM/WPS/recognize_invariants.m" : recognize_invariants;

function Normalize(inv, wgt)
    _, i0 := Max([Abs((inv[i]) ^ (1 / wgt[i])) : i in [1..#inv]]); // find the biggest normalizing factor possible
    inv0 := WPSMultiply(wgt, inv, inv[i0] ^ (-1 / wgt[i0])); // renormlize by that factor, to see which invariants are actually 0
    inv := [Abs(inv0[i]) ge 10 ^ (-Precision / 2) select inv[i] else 0 : i in [1..#inv]];
    return inv;
end function;

function FindCurveHyperelliptic(rosens, F)
    CC := Parent(rosens[1]);
    // recognize Rosenhain invariants

    _<x,y> := PolynomialRing(F, 2);
    f := y*x*(x-y)*(&*[x-r*y : r in rosens]);
    inv, wgt := InvariantsGenus4Curves(f);// : normalize := true);
    vprint CM: "Computing an optimized normalized model...";
    inv_opti := WPSNormalize(wgt, inv);

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

function FindCurveNonHyperelliptic(Q, E, F)
    CC4<X,Y,Z,T> := Parent(Q);
    CC<I> := BaseRing(CC4);
    inv, wgt := InvariantsGenus4Curves(Q, E);
    inv := Normalize(inv, wgt);
        
    if #inv eq 65 then // checks if it is actually rank 3 even if the precision is wrong
        inv_rk3 := InvariantsGenus4Curves(X*T-Y*Z, (Y-Z)^3 : normalize := true);
        if Max([Abs(CC!inv_rk3[i]-inv[i]) : i in [1..65]]) le 10^(-2) then // it is actually of rank 3
            while #inv eq 65 do
                CC<I> := ComplexField(Floor(3/4*Precision(CC)));
                inv, wgt := InvariantsGenus4Curves(ChangeRing(Q, CC), ChangeRing(E, CC));
                inv := Normalize(inv, wgt);
            end while;
        end if;
    end if;

    inv_rec := recognize_invariants(inv, wgt, F);

    return ReconstructionGenus4(inv_rec);
end function;

// Input: taus, a Galois orbit of period matrices
//        F, intersection of field of moduli with the reflex field
function FindCurve(tau, F : prec := Precision(BaseRing(tau)))
    invs := [];
    //for tau in taus do
    tau0 := (tau + Transpose(tau))/2;
    tau_red := SiegelReduction(tau0);
    tau_red := (tau_red + Transpose(tau_red))/2;
    
    vprint CM: "Checking if Jacobian...";
    sch := SchottkyModularFormMagma(tau_red : prec := (prec div 3));
    if Abs(sch) gt 10^(-prec/4) then
        print "Not a jacobian.";
        return Abs(sch);
    end if;
    vprint CM: "It is probably a jacobian.";
    vprint CM: "Norm of the Schottky modular form: %o\n", Abs(sch);

    Eqs := ReconstructCurveG4(tau_red);

    if #Eqs eq 7 then // hyperelliptic case
        vprint CM: "Hyperelliptic case";
        return FindCurveHyperelliptic(Eqs, F);
    elif #Eqs eq 2 then
        vprint CM: "Non-hyperelliptic curve";
        Q, E := Explode(Eqs);
        return FindCurveNonHyperelliptic(Q, E, F);
    end if;

    vprint CM: "Not the right number of equations, error in reconstructing-g4 package";
    return 0,0;
end function;
