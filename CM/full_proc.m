AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/cm-calculations/magma/spec");
AttachSpec("~/github/Genus-4/magma/spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");
AttachSpec("~/github/Reconstruction/magma/spec");
Attach("~/github/Genus-4-RM-CM/flint-stuff/schottky-fast-theta.m");
Attach("~/github/Genus-4-RM-CM/CM/findg4.m");

function Normalize(inv, wgt : prec := Precision(Parent(inv[1])), j := "not def")
    _, i0 := Max([Abs((inv[i]) ^ (1 / wgt[i])) : i in [1..#inv]]); // find the biggest normalizing factor possible

    inv0 := WPSMultiply(wgt, inv, inv[i0] ^ (-1 / wgt[i0])); // renormalize by that factor, to see which invariants are actually 0
    
    if Type(j) eq MonStgElt then
        inv := [Abs(inv0[i]) ge 10 ^ (-prec/2) select inv[i] else 0 : i in [1..#inv]];
    else
        inv := [Abs(inv0[i]) gt Abs(inv0[1]) select inv[i] else 0 : i in [1..#inv]];
    end if;

    res := WPSNormalize(wgt, inv);

    return res;
end function;

function IsDeJongNoot(inv)
    if #inv eq 60 then
        non_zero := [i : i in [1..#inv] | inv[i] ne 0];
        if Max(non_zero) le 5 then 
            return true, 1;
        elif &and[i in [ 4, 5, 9, 20, 24, 49, 57 ] : i in non_zero] then
            return true, 2;
        end if;
    end if;
    return false, 0;
end function;

function FindCurveHyperelliptic(rosensCC) // Fix me warning
    K<x,y> := PolynomialRing(Universe(rosensCC),2);
    return y*x*(x-y)*(&*[x-r*y : r in rosensCC]);

    prec := Precision(Parent(rosensCC[1]));
    K := RationalsExtra(prec);
    L, rosens := NumberFieldExtra(rosensCC, K);
    _<x,y> := PolynomialRing(Parent(rosens[1]), 2);
    
    f := y*x*(x-y)*(&*[x-r*y : r in rosens]);
    inv, wgt := InvariantsGenus4Curves(f);// : normalize := true);
    vprint CM: "Computing an optimized normalized model...";
    inv_opti := WPSNormalize(wgt, inv);

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

        t, f_bis := IsCoercible(PolynomialRing(Rationals(), 2), f);
        printf t;
        if not t then
            return f;
        else
            return MinRedBinaryForm(f_bis), 0, true;
        end if;

    catch e 
        print "No parametrization found, try parametrizing Q using ConicParametrization (https://github.com/JRSijsling/hyperelliptic/blob/main/magma/toolbox/diophantine.m)";
        print "Mestre conic Q and f_5 returned";
        return Q, E, true;
    end try;

end function;

function OliveToIgusa(I)
    J1 := -15*I[1];
    J2 := 45/8*(3*I[1]^2 - 25/2*I[2]);
    J3 := 15/8*(9/2*I[1]^3 - 25*I[1]*I[2] - 375/4*I[3]);
    J4 := (J1*J3-J2^2)/4;
    J5 := 81/16*(-3*I[1]^5 + 7625/256*I[1]^3*I[2] + 13125/256*I[1]^2*I[3] - 9375/128*I[1]*I[2]^2 - 28125/128*I[2]*I[3] - 28125/256*I[4]);
    return [J1,J2,J3,J4,J5];
end function;

function FindCurveNonHyperelliptic(Q, E : K := [Rationals()])
    CC4<X,Y,Z,T> := Parent(Q);
    CC<I> := BaseRing(CC4);
    inv, wgt := InvariantsGenus4Curves(Q, E);
    inv := Normalize(inv, wgt);
    //ChangeUniverse(inv, ComplexField(5));

    prec := Precision(Parent(inv[1]));
    F := RationalsExtra(prec);

    if #inv eq 60 then
        "Rank 3 case";
        if Max([i : i->el in inv | el ne 0]) le 5 then 
            "Only sextic is non-zero";
            return 2, 2, true;
            invs := WPSNormalize([2, 4, 6, 10], inv[1..4]);
            try 
                "Trying LLL...";
                L, invs_alg := NumberFieldExtra(invs, RationalsExtra(prec));
                "LLL worked";
                L;
                printf "De Jong-Noot : %o\n", true;
                Q, E := ReconstructionGenus4(invs_alg cat [0 : i in [1..56]]);
                return Q, E, true;
            catch e 
                "LLL didn't work";
                return 0, 0, 0;
            end try;
        else
            try 
                if IsDeJongNoot(inv) then
                    return 2, 2, true;
                else
                    return 2, 2, false;
                    "Trying LLL sextic...";
                    invs := WPSNormalize([2, 4, 6, 10], inv[1..4]);
                    L, invs_alg := NumberFieldExtra(invs, RationalsExtra(prec));
                    "LLL sextic worked";
                    try 
                        L1, invs_alg := NumberFieldExtra(inv, L);
                        "LLL all worked";
                        printf "De Jong-Noot : %o\n", IsDeJongNoot(invs_alg);
                        Q, E := ReconstructionGenus4(invs_alg);
                        return Q, E, IsDeJongNoot(invs_alg);
                    catch e
                        "LLL all didn't work";
                        return 0,0,0;
                    end try;
                end if;
            catch e 
                "LLL didn't work, should try something more precise";
                return 0,0,0;
            end try;
        end if;
/* need to finish this to cover all cases
        igu := WPSNormalize([2, 4, 6, 8, 10], OliveToIgusa(inv[1..4]));
        ChangeUniverse(igu, ComplexField(30));
        try 
            "Trying LLL...";
            L, igu_alg := NumberFieldExtra(igu, RationalsExtra(prec));
            L;
            "LLL worked";
        catch e 
            "LLL didn't work...";
            return 0,0;
        end try;

        if igu_alg ne [0,0,0,0,0] then  
            igu_alg;
            f_rec := HyperellipticPolynomials(HyperellipticCurveFromIgusaInvariants(igu_alg : minimize := true));
            lis := &cat[Flat(c) : c in Coefficients(f_rec)];
            f_rec *:= GCD([Denominator(c) : c in lis])/GCD([Numerator(c) : c in lis]);
            nu := InfinitePlaces(L);
            i0 := Min([i : i in [1..#nu] | Min([Abs(Evaluate(igu_alg[j], nu[i])-igu[j]) : j in [1..#igu_alg]]) le 10^(-prec/2)]);
            i0; // now we know which embedding it is
        else
            f_rec := 0;
        end if;
        // recognize all invariants which could work
        // 
        if inv[8] ne 0 then
            sec_inv := [1, inv[13]/(inv[8]*inv[1]), inv[24]/(inv[8]*inv[1]^2), inv[43]/(inv[8]*inv[1]^3), inv[54]/(inv[8]*inv[1]^4)];
            L1, sec_alg := NumberFieldExtra(sec_inv, RationalsExtra(prec));
            sec_alg;
        else
            // todo
        end if;*/
   else
        "Rank 4 case";
        return 2, 2, "rank 4";    
        try
            "Trying LLL...";
            L, inv_rec := NumberFieldExtra(inv, F);
            Q, E := ReconstructionGenus4(inv_rec);
            return Q, E, IsDeJongNoot(inv_rec);
        catch e
            "LLL did not work...";
        end try;
    end if;
    return 0,0,0;
end function;

// Input: taus, a Galois orbit of period matrices
//        F, intersection of field of moduli with the reflex field
function FindCurve(tau : prec := Precision(BaseRing(tau)), K := Rationals())
    invs := [];
    //for tau in taus do
    tau0 := (tau + Transpose(tau))/2;
    tau_red := SiegelReduction(tau0);
    tau_red := (tau_red + Transpose(tau_red))/2;
    
    vprint CM: "Checking if Jacobian...";
    sch := SchottkyModularForm(tau_red : prec := (prec div 3), flint := true);
    RealField(5)!Abs(sch);
    if Abs(sch) gt 10^(-prec/4) then
        //Abs(SchottkyModularForm(SiegelReduction(tau_red) : prec := 50, flint := false));
        vprint CM: "Not a jacobian.";
        return 0, 1, 0;
    end if;
    vprint CM: "It is probably a jacobian.";
    vprint CM: "Norm of the Schottky modular form: %o\n", Abs(sch);

    Eqs := ReconstructCurveG4(tau_red : flint := true);

    if #Eqs eq 7 then // hyperelliptic case
        vprint CM: "Hyperelliptic case";
        //return 1, 0, "hyp";
        return FindCurveHyperelliptic(Eqs), 0;
    elif #Eqs eq 2 then
        vprint CM: "Non-hyperelliptic curve";
        Q, E := Explode(Eqs);
        return Q, E;
        //K;
        Q, E, djn := FindCurveNonHyperelliptic(Q, E : K := K);
        return Q, E, djn;
    end if;

    vprint CM: "Not the right number of equations, error in reconstructing-g4 package";
    return 0, 0, 0;
end function;