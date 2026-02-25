AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/cm-calculations/magma/spec");
AttachSpec("~/github/Genus-4/magma/spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");
AttachSpec("~/github/Reconstruction/magma/spec");
Attach("~/github/Genus-4-RM-CM/flint-stuff/schottky-fast-theta.m");
Attach("~/github/Genus-4-RM-CM/CM/findg4.m");

import "CMTypes.m" : Exponents;
exps := Exponents();

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

function InvsRank3(invs)
  invs := [invs[i] : i in [1, 2, 3, 4, 8, 13, 24, 43, 54]];
  invs := Normalize(invs, [6, 12, 18, 30, 8, 14, 20, 26, 32]);
  try 
    prec := Precision(Parent(rosensCC[1]));
    K := RationalsExtra(prec);
    L, _ := NumberFieldExtra(invs[1..2], K);// Alternative when class group not trivial: store them and wait to know the full Galois orbit 
    I := AlgebraizeElementsExtra(invs, L); 
  catch e
    error("Not able to recognize field of moduli");
  end try;
  a11 := I[3];
  a12 := 1/3*(-1/3*I[1]^2*I[2] + I[1]*I[3] + 2*I[2]^2);
  a13 := 1/54*I[1]^3*I[2] - 1/18*I[1]^2*I[3] - 1/9*I[1]*I[2]^2 + 1/3*I[2]*I[3] + 1/2*I[4];
  a22 := 1/54*I[1]^3*I[2] - 1/18*I[1]^2*I[3] - 1/9*I[1]*I[2]^2 + 1/3*I[2]*I[3] + 1/2*I[4];
  a23 := 1/3*(-1/6*I[1]^2*I[2]^2 + 1/3*I[1]*I[2]*I[3] + I[2]^3 + 1/2*I[3]^2);
  a33 := 1/2*(5/162*I[1]^3*I[2]^2 - 7/54*I[1]^2*I[2]*I[3] - 5/27*I[1]*I[2]^3 + 1/9*I[1]*I[3]^2 + 5/9*I[2]^2*I[3] + 1/2*I[2]*I[4]);

  v11 := I[7];
  v12 := I[8];
  v13 := -1/18*I[1]^2*I[2]*I[5] + 1/12*I[1]*I[3]*I[5] + 1/3*I[2]^2*I[5] + 1/4*I[2]*I[7] + 1/4*I[3]*I[6] - 1/2*I[9];
  v22 := I[9];
  v23 := 1/54*I[1]^3*I[2]*I[5] - 1/36*I[1]^2*I[2]*I[6] - 1/18*I[1]^2*I[3]*I[5] -
          1/9*I[1]*I[2]^2*I[5] + 1/36*I[1]*I[2]*I[7] + 1/12*I[1]*I[3]*I[6] +
          1/6*I[2]^2*I[6] + 1/6*I[2]*I[3]*I[5] - 1/12*I[3]*I[7] + 1/4*I[4]*I[5];
  v33 := -1/324*I[1]^4*I[2]*I[5] + 1/108*I[1]^3*I[2]*I[6] + 1/108*I[1]^3*I[3]*I[5] -
          1/27*I[1]^2*I[2]^2*I[5] - 1/36*I[1]^2*I[3]*I[6] - 1/18*I[1]*I[2]^2*I[6] +
          1/18*I[1]*I[2]*I[3]*I[5] + 1/9*I[1]*I[2]*I[8] - 1/12*I[1]*I[4]*I[5] +
          1/3*I[2]^3*I[5] + 1/6*I[2]*I[3]*I[6] - 1/2*I[2]*I[9] + 1/6*I[3]^2*I[5] -
          1/3*I[3]*I[8] + 1/4*I[4]*I[6];

  f111 := I[4];
  f112 := 1/3*(1/54*I[1]^4*I[2] - 1/18*I[1]^3*I[3] - 2/9*I[1]^2*I[2]^2 + 1/3*I[1]*I[2]*I[3] + 1/2*I[1]*I[4] + 2/3*I[2]^3 + I[3]^2);
  f113 := 1/3*(-1/6*I[1]^2*I[2]*I[3] + 1/2*I[1]*I[3]^2 + I[2]^2*I[3] + 1/2*I[2]*I[4]);
  f122 := 1/3*(-1/6*I[1]^2*I[2]*I[3] + 1/2*I[1]*I[3]^2 + I[2]^2*I[3] + 1/2*I[2]*I[4]);
  f123 := 1/3*(1/108*I[1]^4*I[2]^2 - 1/36*I[1]^3*I[2]*I[3] - 1/9*I[1]^2*I[2]^3 + 1/6*I[1]*I[2]^2*I[3] + 1/12*I[1]*I[2]*I[4] +
          1/3*I[2]^4 + 1/2*I[2]*I[3]^2 + 1/2*I[3]*I[4]);
  f133 := 1/3*(-1/972*I[1]^5*I[2]^2 + 1/162*I[1]^4*I[2]*I[3] + 1/81*I[1]^3*I[2]^3 - 1/108*I[1]^3*I[3]^2 - 5/36*I[1]^2*I[2]^2*I[3]
            - 1/36*I[1]^2*I[2]*I[4] - 1/27*I[1]*I[2]^4 + 1/4*I[1]*I[2]*I[3]^2 + 1/12*I[1]*I[3]*I[4] + 11/18*I[2]^3*I[3] +
          1/4*I[2]^2*I[4] + 1/6*I[3]^3);
  f222 := 1/3*(1/36*I[1]^4*I[2]^2 - 1/6*I[1]^3*I[2]*I[3] - 1/3*I[1]^2*I[2]^3 + 1/4*I[1]^2*I[3]^2 + I[1]*I[2]^2*I[3] +
            1/12*I[1]*I[2]*I[4] + I[2]^4 - 1/4*I[3]*I[4]);
  f223 := -1/729*I[1]^5*I[2]^2 + 2/243*I[1]^4*I[2]*I[3] + 4/243*I[1]^3*I[2]^3 - 1/81*I[1]^3*I[3]^2 - 2/27*I[1]^2*I[2]^2*I[3]
          - 1/27*I[1]^2*I[2]*I[4] - 4/81*I[1]*I[2]^4 + 1/12*I[1]*I[2]*I[3]^2 + 1/9*I[1]*I[3]*I[4] + 4/27*I[2]^3*I[3] +
          1/4*I[2]^2*I[4] - 1/36*I[3]^3;
  f233 := 1/5832*I[1]^6*I[2]^2 - 1/972*I[1]^5*I[2]*I[3] - 1/1944*I[1]^4*I[2]^3 + 1/648*I[1]^4*I[3]^2 +
          1/324*I[1]^3*I[2]^2*I[3] + 1/108*I[1]^3*I[2]*I[4] - 1/81*I[1]^2*I[2]^4 - 1/216*I[1]^2*I[2]*I[3]^2 -
          1/36*I[1]^2*I[3]*I[4] + 1/54*I[1]*I[2]^3*I[3] - 11/216*I[1]*I[2]^2*I[4] + 1/18*I[2]^5 + 1/18*I[2]^2*I[3]^2 +
          11/72*I[2]*I[3]*I[4] + 1/8*I[4]^2;
  f333 := -1/1944*I[1]^5*I[2]^3 + 1/648*I[1]^4*I[2]^2*I[3] + 1/162*I[1]^3*I[2]^4 + 1/216*I[1]^3*I[2]*I[3]^2 -
          1/54*I[1]^2*I[2]^3*I[3] - 13/648*I[1]^2*I[2]^2*I[4] - 1/72*I[1]^2*I[3]^3 - 1/54*I[1]*I[2]^5 -
          1/72*I[1]*I[2]^2*I[3]^2 + 1/27*I[1]*I[2]*I[3]*I[4] + 1/18*I[2]^4*I[3] + 1/8*I[2]^3*I[4] + 1/24*I[2]*I[3]^3 +
          5/72*I[3]^2*I[4];
  res := [a11,a12,a13,a22,a23,a33,v11,v12,v13,v22,v23,v33,f111,f112,f113,f122,f123,f133,f222,f223,f233,f333];
  wgt := [18,24,30,30,36,42, 20,26,32,32,38,44, 30,36,42,42,48,54,48,54,60,66];
  
  return res, wgt;
end function;

function FindCurveHyperelliptic(rosensCC) // Fix me warning
    K<x,y> := PolynomialRing(Universe(rosensCC),2);
    f :=  y*x*(x-y)*(&*[x-r*y : r in rosensCC]);

    C_24 := Transvectant(f, f, 8);
    C_28 := Transvectant(f, f, 6);
    c_3 := Transvectant(C_28, f, 8);
    c_5_1 := Transvectant(C_24, c_3, 2);
    c_5_2 := Transvectant(C_28, Transvectant(C_24, f, 4), 6);

    C := [c_3, c_5_1, c_5_2];
    wgt := [3, 5, 5];

    invsQ := [Evaluate(Transvectant(C[i], C[j], 2), [0,0]) : j in [i..3], i in [1..3]];
    invsE := [Evaluate(Transvectant(C[i1]*C[i2]*C[i3]*C[i4]*C[i5], f, 10), [0,0]) : i5 in [i4..3], i4 in [i3..3], i3 in [i2..3], i2 in [i1..3], i1 in [1..3]];
    
    WgtQ := [wgt[i]+wgt[j] : j in [i..3], i in [1..3]];
    WgtE := [wgt[i1]+wgt[i2]+wgt[i3]+wgt[i4]+wgt[i5] : i5 in [i4..3], i4 in [i3..3], i3 in [i2..3], i2 in [i1..3], i1 in [1..3]];

    invs := Normalize(invsQ cat invsE, WgtQ cat WgtE);
    invsQNorm := invs[1..6];
    invsENorm := invs[7..#invs];
    try 
        prec := Precision(Universe(invsQNorm));
        K := RationalsExtra(prec);
        L, _ := NumberFieldExtra(invsQNorm[1..2], K);// Alternative when class group not trivial: store them and wait to know the full Galois orbit 
        _, InvsQInt := AlgebraizeElementsExtra(invsQNorm, L); // linear algebra
        _, InvsEInt := AlgebraizeElementsExtra(invsENorm, L); // linear algebra
    catch e
        error("Not able to recognize field of moduli");
    end try;
    
    LQ := [[i, j] : j in [i..3], i in [1..3]];
    LE := [ [i1, i2, i3, i4, i5]: i5 in [i4..3], i4 in [i3..3], i3 in [i2..3], i2 in [i1..3], i1 in [1..3]];
    
    _<[X]> := PolynomialRing(L, 3);

    Q := &+[InvsQNorm[Index(Sort([i, j]), LQ)]*X[i]*X[j] : i in [1..3], j in [1..3]];
    E := &+[InvsENorm[Index(Sort([i1, i2, i3, i4, i5]), LE)]*X[i1]*X[i2]*X[i3]*X[i4]*X[i5] : i1 in [1..3], i2 in [1..3], i3 in [1..3], i4 in [1..3], i5 in [1..3]];

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

function FindCurveNonHyperelliptic(Q, E : K := [Rationals()])
    CC4<X,Y,Z,T> := Parent(Q);
    CC<I> := BaseRing(CC4);
    inv, wgt := InvariantsGenus4Curves(Q, E);
    if #inv eq 60 then
        invs_rec := InvsRank3(inv);
        Q, E := ReconstructionGenus4Rank3(invs_rec);
    elif #inv eq 65 then
        try
            invsNorm := Normalize(inv, wgt);
            "Trying LLL...";
            prec := Precision(Universe(invsQNorm));
            K := RationalsExtra(prec);
            L := NumberFieldExtra(invsNorm[1..2], K);
            _, inv_rec := AlgebraizeElementsExtra(invsNorm, L);
            Q, E := ReconstructionGenus4(inv_rec);
            return Q, E;
        catch e
            "LLL did not work...";
        end try;
    end if;
end function;

function FindCurveFromCMPoly(datum, expos)
    if expos[eval Split(datum[4], "T")[2]] in {13, 24} then // FIXME later, edge cases C2xA4, C2xS4
        return [* *]; 
    end if;
    R<x> := PolynomialRing(Rationals());
    coeffs := datum[2];
    f := R!coeffs;
    taus := EnumerationUpToGalois(f : exp := expos[eval Split(datum[4], "T")[2]], prec := 5000);
    taus := [t[1] : t in taus];
    curves_f := [* *];
    for tau in taus do 
        L := FindCurve(tau);
        if L ne [] then // actual curve
            Append(~curves_f, L);
        end if;
    end for;
    return curves_f;
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
    sch := SchottkyModularForm(tau_red : prec := (prec div 3));
    if Abs(sch) gt 10^(-prec/4) then
        vprint CM: "Not a jacobian.";
        return [];
    end if;
    vprint CM: "It is probably a jacobian.";
    vprint CM: "Norm of the Schottky modular form: %o\n", Abs(sch);

    Eqs := ReconstructCurveG4(tau_red : flint := true);

    if #Eqs eq 7 then // hyperelliptic case
        vprint CM: "Hyperelliptic case";
        return [FindCurveHyperelliptic(Eqs)];
    elif #Eqs eq 2 then
        vprint CM: "Non-hyperelliptic curve";
        Q, E := Explode(Eqs);
        Q_rec, E_rec := FindCurveNonHyperelliptic(Q, E);
        return [Q_rec, E_rec];
    end if;

    error("Not the right number of equations, error in reconstructing-g4 package"); 
end function;