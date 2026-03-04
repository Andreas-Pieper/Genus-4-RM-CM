Attach("~/github/Genus-4-RM-CM/CM/findg4.m");
AttachSpec("~/github/CHIMP/CHIMP.spec");
Attach("~/github/Genus-4-RM-CM/flint-stuff/FlintWrapper.m");
Attach("~/github/Genus-4-RM-CM/flint-stuff/schottky-fast-theta.m");
Attach("~/github/Genus-4-RM-CM/mumford/NootCMCube.m");
AttachSpec("~/github/cm-calculations/magma/spec");
AttachSpec("~/github/reconstructing-g4/magma/spec");
load "~/lmfdb_8t13.m";
load "~/github/Genus-4-RM-CM/CM/minimize.m";
load "~/github/Genus-4/magma/gordan-10.dat";
load "~/github/Genus-4/magma/invS10.m";
SetDebugOnError(true);
//SetVerbose("CM", 1);

load "~/github/Genus-4-RM-CM/CM/full_proc.m";
load "~/Decomposition/reconstruction_hyp_special2.m";
load "~/Decomposition/Evaluation_interpolation.m";
load "~/invs1.m";

str := Read("~/github/Genus-4-RM-CM/CM/curves_class_8t13.txt");
str := Split(str, "\n");
list_str := [Split(l, "|") : l in str];
list_labels := Multiset([l[1] : l in list_str | l[2] eq "false"]);//, "Not able to recognize invariants"]]); // we keep only the ones we're interested in
list_labels := Setseq(Set([l : l in list_labels | Multiplicity(list_labels, l) eq 3]));// | Multiplicity(list_labels, l) eq 3])); // we keep only the ones with 3 period matrices

function Normalize(inv, wgt : prec := Precision(Universe(inv)))
  _, i0 := Max([Abs((inv[i]) ^ (1 / wgt[i])) : i in [1..#inv]]); // find the biggest normalizing factor possible
  inv0 := WPSMultiply(wgt, inv, inv[i0] ^ (-1 / wgt[i0])); // renormalize by that factor, to see which invariants are actually 0
  inv := [Abs(inv0[i]) ge 10 ^ (-prec/2) select inv[i] else 0 : i in [1..#inv]];
  res := WPSNormalize(wgt, inv);
  return res;
end function;

function Relations(invs, basis_weights)
    prec := Precision(Parent(invs[1][1]));
    m, n := Explode(<#basis_weights, #invs>); 
    mat := Matrix([[Real(Power(inv, b)) : b in basis_weights] : i->inv in invs]);
    //mat := Matrix([[mat[i][j]/mat[i][1] : j in [1..m]] : i in [1..n]]);

    rank := NumericalRank(mat : Epsilon := RR!10^(-prec/2));
    rank;
    if rank ne #basis_weights-1 then
      return [];
    end if;
    M := NumericalKernel(Transpose(mat) : Epsilon := RR!10^(-prec/2));
    if Nrows(M) eq 1 then
        M := Eltseq(M[1]);
        if #M gt 1 then
            if Abs(M[#M]) lt RR!10^(-prec/2) then
              "Invariant not linearly dependent to the others";
              return [];
            else
              return [BestApproximation(-M[i]/M[#M], 10^(3*prec div 4)) : i in [1..#M-1]];
            end if;
        end if;
    end if;
    return [];
end function;

list_curves := [];

for l in list_labels do
  d := data[Index([el[1] : el in data], l)];
  R<x> := PolynomialRing(Rationals());
  d[1];
  coeffs := d[2];
  R!coeffs;
  taus := FullEnumerationG4(R!coeffs : prec := 3000);
  taus := [t[1] : t in taus];
  S<X,Y,Z,T> := PolynomialRing(BaseRing(Parent(taus[1])), 4);
  for i in [1..#taus] do
    tau_red := SiegelReduction(taus[i]);
    Q, E := FindCurve(tau_red);
    if Q ne 0 then
      Append(~list_curves, [S!Q,S!E]);
    end if;
  end for;
end for;

invs := [];

for i in [1..#list_curves] do
  I_all, wgt_all := InvariantsGenus4Curves(list_curves[i][1], list_curves[i][2]);
  I := [wgt_all[i] mod 2 eq 0 select I_all[i] else 0 : i in [1..#I_all]];
  I := WPSNormalize(wgt_all, I);
  Append(~invs, I);
end for;

wgt := wgt_all;

ChangeUniverse([inv[6]/inv[1] : inv in invs], CC);

invs0 := invs;
wgt0 := [6,12,18,30,45,4,6,8,10,10,9,13,14,12,12,15,14,14,15,15,17,16,16,20,19,19,18,16,18,18,17,17,21,19,22,22,23,21,20,20,21,25,26,24,21,23,23,24,25,29,25,28,27,32,29,31,35,33,37,41];//wgt_all;

//l := 1;
ind := [1..60];
//invs := [[inv[i] : i in ind] cat [invs_rec(inv)[l]] : inv in invs0];
invs := [[inv[i] : i in ind] : inv in invs0];

//wgtrec := [ 18, 24, 30, 30, 12, 14, 16, 18, 20, 20, 22, 24, 22, 24, 26, 28, 24, 26, 28, 30, 32, 26, 28, 30, 32, 34, 36 ];
//wgt := [wgt0[i] : i in ind] cat [wgtrec[l]];
wgt := [wgt0[i] : i in ind];

//d := wgtrec[l];

function i3(invs)
  return (-24989524896/817400375*invs[1]^4+283164527296/36783016875*invs[1]^3*invs[7]-201946018576/2979424366875*invs[1]^2*invs[6]^3+10643356864/2150921334375*invs[1]*invs[6]^3*invs[7]-46425728/2150921334375*invs[6]^6)/(invs[1]-28/135*invs[7]);
end function;

function i7_sq(invs)
  return 91125/748*invs[1]^2-4995/187*invs[1]*invs[7]+81/374*invs[6]^3;
end function;

function i3(invs)
  return (-24989524896/817400375*invs[1]^4+283164527296/36783016875*invs[1]^3*invs[7]-201946018576/2979424366875*invs[1]^2*invs[6]^3+10643356864/2150921334375*invs[1]*invs[6]^3*invs[7]-46425728/2150921334375*invs[6]^6)/(invs[1]-28/135*invs[7]);
end function;

function i7mod_sq(invs)
  return 10497600/34969*invs[1]^2+81/374*invs[6]^3;
end function;

/*
S<x,y> := PolynomialRing(QQ, 2);
R<x,y,z,t> := PolynomialRing(QQ, 4);
i1 := t*(x-28/135*z)-(-24989524896/817400375*x^4+283164527296/36783016875*x^3*z-201946018576/2979424366875*x^2*y^3+10643356864/2150921334375*x*y^3*z-46425728/2150921334375*y^6);
i2 := z^2-(91125/748*x^2-4995/187*x*z+81/374*y^3);
i4 := y^6*z-(-1316158875415125/1085201392*x^5+164523172629375/542600696*x^4*z+8497801809045/1381165408*x^3*y^3-377331413765625/1299920384*x^2*t-552080810475/271300348*x^2*y^3*z+8602435278/474775609*x*y^6+48309463125/23212864*t*y^3);
I := Ideal([i1,i2,i4]);

*/

invs0 := invs;
wgt0 := [6,12,18,30,45,4,6,8,10,10,9,13,14,12,12,15,14,14,15,15,17,16,16,20,19,19,18,16,18,18,17,17,21,19,22,22,23,21,20,20,21,25,26,24,21,23,23,24,25,29,25,28,27,32,29,31,35,33,37,41];//wgt_all;

[ChangeUniverse([invs0[i][6]^3/(invs0[i][6]^3-27*invs0[i][7]^2)], ComplexField(8)) : i in [1..60]];


l := 6;
ind := [1,3,6,7];
invs := [[inv[i] : i in ind] cat [invs_rec(inv)[l]]: inv in invs0];
wgtrec := [18,24,30,30,36,42, 20,26,32,32,38,44, 30,36,42,42,48,54,48,54,60,66];
wgt := [wgt0[i] : i in ind] cat [wgtrec[l]];

d := wgtrec[l];

basis_weights := DegreeDBasisWeights(wgt, d);
basis := [NormalForm(R!Power([x,t,y,z], b[1..#b-1]), I) : b in basis_weights[1..#basis_weights-1]];
basis_weights := [basis_weights[i] : i in [1..#basis_weights-1] | #Monomials(basis[i]) eq 1] cat [basis_weights[#basis_weights]];

//basis_weights := [b : b in basis_weights | b[Index(ind, 3)] le 1 and b[Index(ind, 7)] le 1 and b[Index(ind, 7)]*b[Index(ind, 3)] eq 0 and (b[Index(ind, 7)] ne 1 or (b[Index(ind, 6)] mod 6 ne 0 or b[Index(ind, 6)] eq 0))];
//basis_weights := 

#basis_weights;

rel := Relations(invs, basis_weights);
rel := [(r ne 0) select (Sprint(r)[1] eq "-" select Sprint(r) else "+" cat Sprint(r)) else "" : r in rel];
pol := [&cat["*X[" cat Sprint(ind[i]) cat "]^" cat Sprint(basis_weights[d][i]) : i in [1..#ind] | basis_weights[d][i] ne 0] : d in [1..#basis_weights-1]];
&cat[(rel[i] ne "") select rel[i] cat pol[i] else "" : i in [1..#pol]];
//basis_weights;

//6^3, 
//4^5,6*4^3,6*4^3,6^3*(6/4)
d := 12;
ind := [6,7];
invs := [[inv[i] : i in ind] : inv in invs0];
//for i in [1..#invs] do
//  invs[i][Index(ind, 7)] +:= 4995/374*invs[i][1]; 
//end for;
wgt := [wgt0[i] : i in ind];
basis_weights := DegreeDBasisWeights(wgt, d);
basis_weights := [b : b in basis_weights | b[2] le 1];
#basis_weights;

rel := Relations(invs, basis_weights);
rel := [Sprint(r)[1] eq "-" select Sprint(r) else "+" cat Sprint(r) : r in rel];
pol := [&cat["*X[" cat Sprint(ind[i]) cat "]^" cat Sprint(basis_weights[d][i]) : i in [1..#ind] | basis_weights[d][i] ne 0] : d in [1..#basis_weights-1]];
&cat[rel[i] cat pol[i] : i in [1..#pol]];


f := x^2*z^3-54675/392*y^2*x^2-165350706912/32696015*y^5*x+23614191424/780704780625*y*x*z^6+11546288733332/360473565375*y^3*x*z^3;
I := Ideal([z^3-54675/392*y^2]);
NormalForm(f, I);
//18 : dim 5 4 5
//20 : dim 4 3 3
//22 : dim 5 4 5
//24 : dim 6 5 6 (relation)
//26 : dim 5 4 5
//28 : dim 6 5 6 (relation)
//30 : 

//[basis_weights[1]] cat [basis_weights[4]] cat [basis_weights[10]]
// qd renormalisé, inv1 = inv6, et sinon on peut essayer d'en trouver un autre qui a une relation avec eux de plus haut degré

Relations(invs, [basis_weights[1]] cat [basis_weights[2]] cat [basis_weights[4]] cat [basis_weights[6]] cat [basis_weights[10]]);
Relations(invs,  [[-1,0,0,0,0,3] cat [0 : i in [7..60]]] cat [[-4,0,0,0,0,6] cat [0 : i in [7..60]]] cat [[0,0,0,0,0,0,1] cat [0 : i in [8..60]]]);

    //Write("invariants_hyp.m", Sprint(d[1]) cat "|" cat Sprint(R!DefiningPolynomial(Rationals())) cat "|" cat Sprint(res) cat "|");
    
   // L0, invs := NumberFieldExtra(invsCC0[1..7], NumberFieldExtra(R![1,-2,-1,1] : prec := 5000));
    //L, invs := NumberFieldExtra(invsCC0, L0);
   // invs1 := WPSNormalize(wgteven[1..7], invs);
   // invs1;

   // Write("invariants_hyp.m", Sprint(d[1]) cat "|" cat Sprint(R!DefiningPolynomial(L0)) cat "|" cat Sprint(invs1) cat "|");
//  end if;
//end for;


//for d in data do
  if d[5] ne [] and d[1] in L then
    R<x> := PolynomialRing(Rationals());
    d[1];
    coeffs := d[2];
    R!coeffs;
    taus := FullEnumerationG4(R!coeffs : prec := 5000);
    tau := [t[1] : t in taus];
    for i in [1..#tau] do 
      tau_red := SiegelReduction(tau[i]);
      f := FindCurve(tau_red);//Q, E, djn := FindCurve(tau_red);
      if E eq 0 then
        if Q ne 0 then
          S := Parent(Q);
          try 
            K := BaseRing(S);
          catch e
            K := Rationals();
          end try;
          //Write("curves_class_8t13.txt", Sprint(d[1]) cat "|" cat Sprint(djn) cat "|" cat Sprint(R!DefiningPolynomial(K)) cat "|" cat Sprint(Q) cat "|");
        else
          //Write("curves_class_8t13.txt", Sprint(d[1]) cat "|Not able to recognize invariants");
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
//end for;


R<c0,c1,c2,c3,c4,c5,c6,d0,d1,d2,d3,d4,a,b> := PolynomialRing(QQ, 14); 
S<x,y,t> := PolynomialRing(R, 3);

//f := x^3+x+t^7+b;
f := y^2-(x^3+x+t^7-7*1*t^5+14*1^2*t^3-7*1^3*t+1);

//f1 := Evaluate(f, [d0+d1*t+d2*t^2+d3*t^3+1*t^4,(3/4*d1^2+3/2*d0*d2-1/8*d2^3-3/4*d1*d2*d3-3/8*d0*d3^2+9/32*d2^2*d3^2+3/16*d1*d3^3-15/128*d2*d3^4+7/512*d3^6-3/2*d3)/2+(3/2*d1*d2+3/2*d0*d3-3/8*d2^2*d3-3/8*d1*d3^2+3/16*d2*d3^3-3/128*d3^5+1)/2*t+(3*d0+3/4*d2^2+3/2*d1*d3-3/8*d2*d3^2+3/64*d3^4)/2*t^2+(3*d1+3/2*d2*d3-1/8*d3^3)/2*t^3+(3*d2+3/4*d3^2)/2*t^4+3/2*d3*t^5+1*t^6,t]); 
f1 := Evaluate(f, [d0+d1*t+d2*t^2+d3*t^3+d4*t^4,c0+c1*t+c2*t^2+c3*t^3+c4*t^4+c5*t^5+c6*t^6,t]); 
f1;

I := Ideal(Coefficients(f1));
EliminationIdeal(I, 13);

GroebnerBasis(Coefficients(f1), 5);