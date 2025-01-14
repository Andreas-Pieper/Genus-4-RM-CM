AttachSpec("~/github/CHIMP/CHIMP.spec");
AttachSpec("~/github/cm-calculations/magma/spec");
AttachSpec("~/github/Genus-4/magma/spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");
AttachSpec("~/github/Reconstruction/magma/spec");
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

/*
R<lambda, mu> := PolynomialRing(Rationals(), 2);
S<x,y,z,t> := PolynomialRing(R, 4);

Q := x*z-y^2;
E := mu*x^2*y-t*(t-z)*(t-lambda*z);
InvariantsGenus4Curves(Q,E);
*/

// Shimura curves
function IsOnShimuraHyp(f : epsilon := Precision(BaseRing(Parent(f))) div 2)
    invs_f, wgt := InvariantsGenus4Curves(f);
    wgt_even := [w : w in wgt | w mod 2 eq 0];
    invs_even := [inv : i->inv in invs_f | wgt[i] mod 2 eq 0]; 
    invs_norm := Normalize(invs_even, wgt_even);
    mu := invs_norm[4];
    invs_shim := 
    [0,
    0,
    8/225*mu,
    mu,
    9555/1552661*mu,
    -136/7425*mu,
    -136/7425*mu,
    22321/1111320*mu,
    41905/373527*mu,
    -1088/121275*mu,
    32453/271656*mu,
    9378917/16180819200*mu,
    -80580799/82771113600*mu,
    3965075087/373776923520*mu,
    -10106041/343429632000*mu,
    6679079/135136512000*mu,
    8415/33575584*mu,
    6432239/92989769600*mu,
    2964307/10679340672*mu,
    799098587/649689064113375*mu^2,
    -35403551/965268798243750*mu^2,
    -737872919/3467287952689275*mu^2,
    104367704/1449535097278125*mu^2 - 27781896571/2628148399027200*mu,
    202166944/543307229053875*mu^2 - 47535787/533045971584*mu,
    24971192/1353269001459375*mu^2 - 1527999253/419581586511360*mu,
    -13955952/19843977993125*mu^2 + 82387485713/827109665907840*mu,
    -33747224/55202338780875*mu^2 - 6182478247/708951142206720*mu,
    45834047/222225660421875*mu^2 - 38197116529/1618258041993600*mu,
    88783492/83293469758125*mu^2 - 91289913/875248323712*mu,
    -112232344/6243925288535*mu^2 + 8403525/51353855728*mu,
    59920516/113310063813375*mu^2 - 146795/4506427296*mu,
    -61538022538343/300181687280798400000*mu^2 + 21026458789611847/4307290703943769350144*mu,
    419176166767579/550333093348130400000*mu^2 - 2817979013040763/89735222998828528128*mu,
    17360603883629/83274086493467100000*mu^2 + 1125455259401629/27156712223329686144*mu,
    92179526694239/150090843640399200000*mu^2 - 52350137286370495/2153645351971884675072*mu,
    1397787163426537/1284110551145637600000*mu^2 - 8124306944780585/209382186997266565632*mu,
    1918564676253251/22623693164729506080000*mu^2 - 57104093245637555/11934784658844194241024*mu,
    36255126482249/431511175466147700000*mu^2 + 88287259051679129/24766921547676673763328*mu,
    1096637001599/35086171240612800000*mu^2 + 402739710504703/1926241805490381324288*mu,
    -5117181688447/64324647274456800000*mu^2 - 53975425119187/40130037614382944256*mu,
    -20211972884261/126533352204359100000*mu^2 + 21556912169221/12144616646457996288*mu,
    -15266198370251/228060113063983200000*mu^2 - 1002711837811255/963120902745190662144*mu,
    -19578822922741/150090843640399200000*mu^2 - 155612557479665/93636754433560203264*mu,
    -9230069858983/327836412529475850000*mu^2 + 1691049620266721/11075890381569692614656*mu,
    169805623/20391829387776000*mu^2,
    -89401/105931581235200*mu^2,
    237610238567/1567601492355892224000*mu^2,
    -1524484215210463/6342919314726500352000*mu^2,
    -7661196259085981/74138017964335718400000*mu^2,
    -7753465235717563/14684835500764910910000000*mu^3 + 960201643633802107/1494684528918996602880000*mu^2];
    RealField(10)!Max([Abs(invs_shim[i]-invs_norm[i]) : i in [1..#invs_shim]]);
    return Max([Abs(invs_shim[i]-invs_norm[i]) : i in [1..#invs_shim]]) lt epsilon;
end function;


function GetMu(f)
  c3 := Transvectant(f, f, 8);
  c4 := Transvectant(f, f, 6);
  c5 := Transvectant(f, f, 4);
  c6 := Transvectant(f, f, 2);
  c7 := Transvectant(c4, f, 8);
  c9 := Transvectant(c5, f, 8);
  c12 := Transvectant(c6, f, 7);
  c22 := Transvectant(c12, f, 9);

  i4 := Transvectant(c7, c7, 2);
  i6 := Transvectant(c3*c9, f, 10);
  i7 := Transvectant(c9*c22, f, 10);
  inv := Evaluate(i4*(i6/i7)^3, [0,0]);
  
  _, L := NumberFieldExtra([inv], NumberFieldExtra(PolynomialRing(Rationals())![1,-2,-1,1] : prec := Precision(Parent(inv))));
  
  return 1/46814919504949789567799425*(L[1]-359165724775/3794886144);
end function;

function GetShimuraCurve(mu)
  S<X,Y> := PolynomialRing(Parent(mu), 2);
  return (-1432353454804332119596592722603696979968*mu - 4886589248453472)*X^10 - 309519*X^9*Y + 9851363924882199552*mu*X^8*Y^2 + (-104516331448118333562141904207872*mu^2 + 118855296*mu)*X^7*Y^3 - 2837192810366073470976*mu^2*X^6*Y^4 - 17115162624*mu^2*X^5*Y^5 - 90790169931714351071232*mu^3*X^4*Y^6 - 1408333381632*mu^3*X^3*Y^7 - 26288889790464*mu^4*X*Y^9;
end function;

function GetShimuraCurve(f)
  try 
    mu := GetMu(f);
  catch e
    "Not enough precision, try again.";
    return 0;
  end try;
  S<X,Y> := PolynomialRing(Parent(mu), 2);
  return (-1432353454804332119596592722603696979968*mu - 4886589248453472)*X^10 - 309519*X^9*Y + 9851363924882199552*mu*X^8*Y^2 + (-104516331448118333562141904207872*mu^2 + 118855296*mu)*X^7*Y^3 - 2837192810366073470976*mu^2*X^6*Y^4 - 17115162624*mu^2*X^5*Y^5 - 90790169931714351071232*mu^3*X^4*Y^6 - 1408333381632*mu^3*X^3*Y^7 - 26288889790464*mu^4*X*Y^9;
end function;
