AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/cm-calculations/magma/spec");
AttachSpec("~/github/Genus-4/magma/spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");
AttachSpec("~/github/Reconstruction/magma/spec");
Attach("~/FlintWrapper.m");
Attach("~/github/Genus-4-RM-CM/flint-stuff/schottky-fast-theta.m");
Attach("~/github/Genus-4-RM-CM/CM/findg4.m");
/*
SetVerbose("CM",true);
SetVerbose("Reconstruction",true);
SetVerbose("Theta",true);
SetVerbose("Genus4",true);
SetDebugOnError(true);
*/

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

function FindCurveHyperelliptic(rosens)
    prec := Precision(Parent(rosens[1]));
    _<x,y> := PolynomialRing(Parent(rosens[1]), 2);
    // recognize Rosenhain invariants

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
    ChangeUniverse(inv, ComplexField(5));
    /*
    if #inv eq 65 then // checks if it is actually rank 3 even if the precision is wrong
        inv_rk3 := InvariantsGenus4Curves(X*T-Y*Z, (Y-Z)^3 : normalize := true);
        if Max([Abs(CC!inv_rk3[i]-inv[i]) : i in [1..65]]) le 10^(-2) then // it is actually of rank 3
            while #inv eq 65 do
                CC<I> := ComplexField(Floor(3/4*Precision(CC)));
                "New prec : ", Precision(CC);
                Q := ChangeRing(Q, CC);
                E := Parent(Q)!E;
                inv, wgt := InvariantsGenus4Curves(Q,E);
            end while;
        end if;
    end if;*/
    prec := Precision(Parent(inv[1]));
    F := RationalsExtra(prec);

    if #inv eq 60 then
        "Rank 3 case";
        if Max([i : i->el in inv | el ne 0]) le 5 then 
            "Only sextic is non-zero";
            invs := WPSNormalize([2, 4, 6, 10], inv[1..4]);
            try 
                "Trying LLL...";
                L, invs_alg := NumberFieldExtra(invs, RationalsExtra(prec));
                "LLL worked";
                L;
                return ReconstructionGenus4(invs_alg cat [0 : i in [1..56]]);
            catch e 
                "LLL didn't work";
                return 0,0;
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
        try
            "Trying LLL...";
            L, inv_rec := NumberFieldExtra(inv, F);
            return ReconstructionGenus4(inv_rec);
        catch e
            "LLL did not work...";
        end try;
    end if;
    return 0, 0;
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
    RealField(10)!Abs(sch);
    prec;
    if Abs(sch) gt 10^(-prec/4) then
        vprint CM: "Not a jacobian.";
        return 0, 0;
    end if;
    vprint CM: "It is probably a jacobian.";
    vprint CM: "Norm of the Schottky modular form: %o\n", Abs(sch);

    Eqs := ReconstructCurveG4(tau_red : flint := true);

    if #Eqs eq 7 then // hyperelliptic case
        vprint CM: "Hyperelliptic case";
        return FindCurveHyperelliptic(Eqs), 0;
    elif #Eqs eq 2 then
        vprint CM: "Non-hyperelliptic curve";
        Q, E := Explode(Eqs);
        //K;
        return FindCurveNonHyperelliptic(Q, E : K := K);
    end if;

    vprint CM: "Not the right number of equations, error in reconstructing-g4 package";
    return 0, 0;
end function;

 
/*
load "fields_pot.m";

R<x> := PolynomialRing(Rationals());

for l in L_fields do
    l;
    f := R!l;
    taus := FullEnumerationG4(f : prec := 300); 
    sch := [SchottkyModularForm(t) : t in taus];
    for i in [1..#taus] do
        Write("schottky-values.txt", Sprint(l) cat "|" cat Sprint(i) cat "|" cat Sprint(RealField(100)!Abs(sch[i])));
        if Abs(sch[i]) le 10^(-100) then
            Write("jacobians.txt", Sprint(l) cat "|" cat Sprint(i) cat "|" cat Abs(sch[i]));
        end if;
    end for;
end for;
*/