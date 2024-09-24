function QuadraticFormToMatrix(Q)
        R := Parent(Q);
        K := CoefficientRing(R);
        Q_mat := Matrix(K, [[MonomialCoefficient(Q, R.i*R.j)/2 : j in [1..4]] : i in [1..4]]);

        for i in [1..4] do
                Q_mat[i][i] := 2*Q_mat[i][i];
        end for;

        return Q_mat;
end function;

function NewBasis(Q)
    K := BaseRing(Parent(Q));
    D, P := DiagonalForm(Q);
    D := QuadraticFormToMatrix(D);
    if IsExact(K) then
        t := Rank(D);
    else
        vprint Genus4 : "Base ring is precision field; using numerical algorithms";
        prec := Precision(K);
        RR := RealField(prec);
        eps := RR!10^(-prec*7/8);
        t := NumericalRank(D : Epsilon:=eps);
    end if;
    vprint Genus4 : "rank =", t;

        if t lt 3 then
                "The quadric is not of rank 3 or 4";
                return D, t, 1;

        elif t eq 4 then
                L := [-D[4][4]/D[1][1], -D[3][3]/D[2][2]];
                bool1 := IsPower(L[1], 2);
                bool2 := IsPower(L[2], 2);

                if bool1 and bool2 then
                        S := K;
                        _, sq1 := IsPower(L[1], 2);
                        _, sq2 := IsPower(L[2], 2);
                        Sq := [sq1, sq2];
                else
                        _<x> := PolynomialRing(K);
                        S := SplittingField([x^2-L[1], x^2-L[2]]);
                        Sq := [Sqrt(S!L[1]), Sqrt(S!L[2])];
                end if;

                M2 := KMatrixSpace(S,4,4);
                P := ChangeRing(P, S);
                P_fin := (M2![S!1/(2*D[1][1]),0,0,S!1/(2*D[1][1]*Sq[1]),0,S!-1/(2*D[2][2]),S!-1/(2*D[2][2]*Sq[2]),0,0,S!1/2,S!-1/(2*Sq[2]),0,S!1/2,0,0,S!-1/(2*Sq[1])])*P;
                //Determinant(P_fin)^2;
                return P_fin, 4, 1;

        else
                i := 1;
                if IsExact(K) then
                    while (D[i][i] ne 0) and (i lt 4) do
                            i := i+1;
                    end while;
                else
                    _, i := Min([Abs(D[j][j]) : j in [1..4]]);
                    //"min", i;
                end if;

                L_swap := [1,2,3,4];
                L_swap[i] := 4;
                L_swap[4] := i;
                P_swap := PermutationMatrix(K, L_swap);
                //P_swap;

                D := P_swap*D*P_swap;
                P := P_swap*P;
                L := [-D[3][3]/D[1][1], -D[2][2]];
                bool1 := IsPower(L[1], 2);
                bool2 := IsPower(L[2], 2);

                if bool1 and bool2 then
                        S := K;
                        _, sq1 := IsPower(L[1], 2);
                        _, sq2 := IsPower(L[2], 2);
                        //sq1^2, L[1];
                        //sq2^2, L[2];
                        Sq := [sq1, sq2];
                else
                        _<x> := PolynomialRing(K);
                        S := SplittingField([x^2-L[1], x^2-L[2]]);
                        Sq := [Sqrt(S!L[1]), Sqrt(S!L[2])];
                end if;

                M2 := KMatrixSpace(S,4,4);
                P := ChangeRing(P, S);
                P_fin := (M2![S!1/(2*D[1][1]),0,S!1/(2*D[1][1]*Sq[1]),0,0,S!1/Sq[2],0,0,S!1/2,0,S!-1/(2*Sq[1]),0,0,0,0,S!1])*P;
                //Determinant(P_fin)^2;
                return P_fin, 3, Sq;
        end if;
end function;


// given a form and a change of variables (given as a matrix), returns the form after the change of variables
function ChangeOfBasis(C, P)
        R := ChangeRing(Parent(C), Parent(P[1][1]));
        C := R!C;
        P := Transpose(Matrix([[R!P[i][j] : j in [1..NumberOfColumns(P)]] : i in [1..NumberOfRows(P)]]));
        Y := ElementToSequence(P*Matrix([[R.i] : i in [1..Rank(R)]]));
        return Evaluate(C, Y);
end function;


/*
for d in Dec do
    cont := true;
    while cont do
        _, p := IsPrincipal(d[1]);
        Fp, m := ResidueClassField(OK, d[1]);//^d[2]);
        Fpx<x,y,z> := PolynomialRing(Fp, 3);
        if IsPrime(#Fp) then
            inc := hom<Fp -> Fpx | >;
        else
            inc := hom<Fp -> Fpx | Fp.1>;
        end if;
        h := hom<S -> Fpx | m*inc,[x, y, z]>;
        if Characteristic(Fp) ne 2 then
            Diag, I := DiagonalForm(h(S!Q1));
        
        /*Fac := Factorization(h(S!Q1));
        Fac;
        cont := false;
/*        if #Fac eq 1 then // weight [0,0,1]
            l := Fac[1][1];
            if Degree(l) eq 1 then
                ind := Min([i : i in [1..Rank(Fpx)] | MonomialCoefficient(l, Fpx.i) ne 0]);
                I := IdentityMatrix(Fp, Rank(Fpx));
                for j in [1..Rank(Fpx)] do
                    if j ne ind then
                        I[ind][j] := -MonomialCoefficient(l, Fpx.j)/MonomialCoefficient(l, Fpx.ind); ;
                    end if;
                end for;
                ind := #Coefficients(Diag);
                I := Transpose(Matrix(OK, [[I[i,j] @@ m : j in [1..Ncols(I)]] : i in [1..Ncols(I)]]))*DiagonalMatrix(OK, [i gt ind select 1 else p : i in [1..Rank(S)]]);
                Q2 := Q1^I;
                e := Min([Valuation(c, d[1]) : c in Coefficients(Q2)]);
                if e gt 2*ind/3 then
                    Q2 := Q2 div p^e;
                    T := T*I;
                    Q1 := Q2;
                else 
                    cont := false;
                end if;
        else
            cont := false;
        end if;
            //else
            //    cont := false;
            //end if;
        //else
        //    cont := false;
        //elif #Fac eq 2 then // weight [0,1,1]
        //end if;
    end while;
end for;

D := &*Coefficients(Q1);
Dec := Decomposition(OK!D);

Q2 := Q^T0;

D1 := &*Coefficients(Q1);
D;
D1;
*/
/*K := NumberFieldExtra(PolynomialRing(QQ)![1,-1,-4,0,1] : prec := 300); 
OK := RingOfIntegers(K);
S<X,Y,Z,T> := PolynomialRing(K, 4);
Q := 1/128625*(575000*K.1^3 + 221900*K.1^2 - 2208810*K.1 - 1419318)*X^2 + 1/4501875*(-153115500*K.1^3 - 60826800*K.1^2 + 588431340*K.1 + 387024964)*X*Y + 1/157565625*(9159221750*K.1^3 + 3628138250*K.1^2 - 35197618575*K.1 - 23100730811)*X*Z + 1/315131250*(9159221750*K.1^3 + 3628138250*K.1^2 - 35197618575*K.1 - 23100730811)*Y^2 + 1/5514796875*(859371821500*K.1^3 + 340574405200*K.1^2 - 3302487234240*K.1 - 2168159942836)*Y*Z + 1/386035781250*(-42873025928375*K.1^3 - 16992558908450*K.1^2 + 164757508535595*K.1 + 108174078005139)*Z^2;
E := 1/157565625*(13058247000*K.1^3 + 5172827250*K.1^2 - 50181188625*K.1 - 32935574353)*X^3 + 1/11029593750*(1590636346250*K.1^3 + 630207030350*K.1^2 - 6112630599645*K.1 - 4012401111193)*X^2*Y + 1/64339296875*(830526951625*K.1^3 + 328717613650*K.1^2 - 3191541106215*K.1 - 2093665398753)*X^2*Z + 1/64339296875*(830526951625*K.1^3 + 328717613650*K.1^2 - 3191541106215*K.1 - 2093665398753)*X*Y^2 + 1/6755626171875*(-17827862551428875*K.1^3 - 7065902545205600*K.1^2 + 68510982690997815*K.1 + 44981533906551157)*X*Y*Z + 1/472893832031250*(-188877559323955375*K.1^3 - 74859896879675275*K.1^2 + 725840666494777545*K.1 + 476557844374624832)*X*Z^2 + 1/13511252343750*(-10227561207178375*K.1^3 - 4053586423233325*K.1^2 + 39303658442463855*K.1 + 25805154188013094)*Y^3 + 1/236446916015625*(1071473607310420750*K.1^3 + 424666091100292675*K.1^2 - 4117582772521909665*K.1 - 2703429618689275709)*Y^2*Z + 1/22068378828125000*(-23183327866866549000*K.1^3 - 9188452638183550125*K.1^2 + 89091578053401622950*K.1 + 58493777103891816551)*Y*Z^2 + 1/579294944238281250*(-4363276308376373695625*K.1^3 - 1729334545308948675800*K.1^2 + 16767703335710978781960*K.1 + 11008963299922388808644)*Z^3 + T^3;
*/
K := NumberFieldExtra(PolynomialRing(QQ)![-2,-11,-13,-1,1] : prec := 300); 
OK := RingOfIntegers(K);
S<X,Y,Z,T> := PolynomialRing(K, 4);
Q := 1/1895803875*(-113192379240*K.1^3 - 393317496940*K.1^2 - 289176226020*K.1 - 50465818102)*X^2 + 1/199059406875*(27909354977480*K.1^3 + 97016599239280*K.1^2 + 71399792929240*K.1 + 12486656012488)*X*Y + 1/9593668114340625*(63564159650145615200*K.1^3 + 220940077675254027200*K.1^2 + 162562718109700509600*K.1 + 28402911703969489776)*X*Z + 1/9593668114340625*(31782079825072807600*K.1^3 + 110470038837627013600*K.1^2 + 81281359054850254800*K.1 + 14201455851984744888)*Y^2 + 1/10782216997395046875*(2318375735664822965564640*K.1^3 + 8058349296616251492213440*K.1^2 + 5929157904454592826599520*K.1 + 1035946612557618344216096)*Y*Z + 1/10189195062538319296875*(1344547221932520565005812480*K.1^3 + 4673457773144208550537008480*K.1^2 + 3438628374280220932880979840*K.1 + 600799493520328936501787296)*Z^2;
E := 1/9593668114340625*(121863315039620964500*K.1^3 + 423579740878817711000*K.1^2 + 311660512091604493500*K.1 + 54453420629511691548)*X^3 + 1/97039952976555421875*(60847943026249576628246240*K.1^3 + 211498925271436643186755040*K.1^2 + 155616299759256930918736320*K.1 + 27189386374042105683629016)*X^2*Y + 1/10189195062538319296875*(-1106280940364585110919801360*K.1^3 - 3845277558233094218003667360*K.1^2 - 2829271437963319159976180880*K.1 - 494332251981999769089831952)*X^2*Z + 1/10189195062538319296875*(-1106280940364585110919801360*K.1^3 - 3845277558233094218003667360*K.1^2 - 2829271437963319159976180880*K.1 - 494332251981999769089831952)*X*Y^2 + 1/18187713186630899944921875*(-637928338556078404881716107588480*K.1^3 - 2217349533528498274674952144434880*K.1^2 - 1631477468763782701802479767431040*K.1 - 285052862357244240830245737885504)*X*Y*Z + 1/551906156648314658828654296875*(-596909276718814918907335327825035328640*K.1^3 - 2074773021194512542636924219366595507040*K.1^2 - 1526572777849809630228108740565380872320*K.1 - 266723842809660003442969394065530110784)*X*Z^2 + 1/54563139559892699834765625*(627599312497303610636041802445760*K.1^3 + 2181447286028790938403344806562560*K.1^2 + 1605061377442051416913233111924480*K.1 + 280437424793363654627775196657248)*Y^3 + 1/551906156648314658828654296875*(345549311080612538902173117593800111360*K.1^3 + 1201080995194575522607110146908985272960*K.1^2 + 883729223643035265672911062807608967680*K.1 + 154405776097231187297102292544687089216)*Y^2*Z + 1/30679489296038667799592841796875*(113299791831126545439035406007565690355200*K.1^3 + 393814203542554807994838122654869906491200*K.1^2 + 289759909405693717155500644765548234185600*K.1 + 50627050115307238388599847225698887706496)*Y*Z^2 + 1/310323034229431124792881594775390625*(17245553711074802706516180557963752986578304320*K.1^3 + 59943128664349286718267114167330793651451550720*K.1^2 + 44104847857266964863100849611965491468806413760*K.1 + 7706029268774172986616484993717003332487659648)*Z^3 + T^3;


function NumberFieldLCM(L)
  if #L eq 0 then
    return 1;
  end if;
  K := Parent(L[1]);
  OK := Integers(K);
  if K eq Rationals() then
    return LCM(ChangeUniverse(L, Integers()));
  end if;
  I := [OK*x : x in L];
  _, a := IsPrincipal(&meet I);
  return K!a;
end function;

function NumberFieldGCD(L)
  if #L eq 0 then
    return 1;
  end if;
  K := Parent(L[1]);
  OK := Integers(K);
  if Type(K) eq FldRat then
    return GCD(ChangeUniverse(L, Integers()));
  elif Type(K) eq FldNum then
    I := [OK*x : x in L];
    _, a := IsPrincipal(&+I);
    return K!a;
  end if;
  I := [OK*x : x in L];
  _, a := IsPrincipal(&+I);
  return K!a;
end function;

function ElimGCD(F)
    return F div NumberFieldGCD(Coefficients(F));
end function;

function MinimizeQuadric(Q, E)
    OK := RingOfIntegers(BaseRing(Parent(Q)));
    S := PolynomialRing(OK, Rank(Parent(Q)));
    
    Q, T0 := DiagonalForm(Q);
    T0 := Transpose(T0);
    E := E^T0;
    
    Q1 := Q*LCM([Denominator(c) : c in Coefficients(Q)]);
    E1 := E*LCM([Denominator(c) : c in Coefficients(E)]);
    Q1 := ElimGCD(S!Q1);
    E1 := ElimGCD(S!E1);

    Dec := Decomposition(&*Coefficients(Q1));
    T := IdentityMatrix(OK, Rank(S));

    for d in Dec do
        cont := true;
        while cont do
            _, p := IsPrincipal(d[1]);
            Fp, m := ResidueClassField(OK, d[1]);
            Fpx := PolynomialRing(Fp, Rank(S));
            if IsPrime(#Fp) then
                inc := hom<Fp -> Fpx | >;
            else
                inc := hom<Fp -> Fpx | Fp.1>;
            end if;
            h := hom<S -> Fpx | m*inc,[Fpx.i : i in [1..Rank(S)]]>;
            if Characteristic(Fp) ne 2 then
                Diag, I := DiagonalForm(h(S!Q1));
                ind := #Coefficients(Diag);
                I := Transpose(Matrix(OK, [[I[i,j] @@ m : j in [1..Ncols(I)]] : i in [1..Ncols(I)]]))*DiagonalMatrix(OK, [i gt ind select 1 else p : i in [1..Rank(S)]]);
                Q2 := Q1^I;
                e := Min([Valuation(c, d[1]) : c in Coefficients(Q2)]);
                if e gt 2*ind/3 then
                    Q2 := Q2 div p^e;
                    T := T*I;
                    Q1 := Q2;
                else 
                    cont := false;
                end if;
            else
                cont := false;
            end if;
        end while;
    end for;
    return Q1, ElimGCD(E1^T);
end function;


//////////////////////////////////////////////////////////////////
SubsForLFac := function(fac, rk, p, m, h)
 local lf, tran,j,i;

 Coef, Mon := CoefficientsAndMonomials(fac);
 Coef_lift := [c @@ m : c in Coef];
 S := Domain(h);
 S1 := Parent(fac);
 Mon_lift := [&*[(S.j)^(Degree(mon, S1.j)) : j in [1..Rank(S)]] : mon in Mon];
 lf := rk!Polynomial(Coef_lift, Mon_lift);
 j := 1;
 while MonomialCoefficient(lf, rk.j) ne 1 do j := j+1; end while;
 tran := [rk.i : i in [1..4]];
 tran[j] := Evaluate(2 * Parent(lf).j - lf,j,p*rk.j);
 
 return tran;
end function;

TestRedMinimization := function(f,p,m,h)
 local tran, res, tmp;

/* Weight vector 0,0,0,1 -- meaning reducible reduction -- test it */
 fac := Factorization(h(f));
 for tmp in fac do 
  if TotalDegree(tmp[1]) eq 1 then
   tran := SubsForLFac(tmp[1], Parent(f), p, m, h);
   res := Evaluate(f,tran);
   gcd := NumberFieldGCD(Coefficients(res));
   try  /* Equation gets divisible by the transformation */
    res := res div gcd;
   catch e
    "p :=", p;
    "gcd :=", gcd;
    "p should divide gcd";
     return false, 0,0;
    end try;
   return true,res,tran;
  end if;
 end for;
 return false, 0,0;
end function;

/* The case of singular lines */
ListSingularLines := function(f)
 local sing,lines,pd,i,gb,tmp;  

 lines := [];
 sing := ideal<Parent(f) | f, Derivative(f,1), Derivative(f,2), Derivative(f,3), Derivative(f,4)>; 
 pd := PrimaryDecomposition(Radical(sing));
 for i := 1 to #pd do
  gb := GroebnerBasis(pd[i]);
  if #gb eq 2 and &and [Degree(tmp) eq 1 : tmp in gb] then Append(~lines,gb); end if;
 end for;
 return lines;
end function;

/* Moves singlar line to rk.1 = rk.2 = 0 

   The lines is given by 2 linear equations over GF(p). 
   The resulting transformations is defined in rk.       */
NormalizeLine := function(line, rk, m)
 local i,j,tran, k,perm,c;

  i := 1;  
  while MonomialCoefficient(line[1],Parent(line[1]).i) ne 1 do i := i+1; end while;
  j := i+1;
  while MonomialCoefficient(line[2],Parent(line[2]).j) ne 1 do j := j+1; end while;
 /* We replace the i-th variable using the first and the j-th variable using the second equation. */

  tran := [rk.i : i in [1..4]];
  for k := 1 to 4 do 
   if i ne k then
    tran := [Evaluate(akt,i,rk.i - rk.k*(MonomialCoefficient(line[1],Parent(line[1]).k)) @@ m) : akt in tran];
   end if;
  end for; 
  for k := 1 to 4 do
   if j ne k then
    tran := [Evaluate(akt,j,rk.j - rk.k*(MonomialCoefficient(line[2],Parent(line[2]).k)) @@ m) : akt in tran];
   end if;
  end for;

  perm := [rk!0 : k in [1..4] ];
  c := 3;
  for k := 1 to 4 do 
   if k eq i then perm[k] := rk.1; else
    if k eq j then perm[k] := rk.2; else
     perm[k] := rk.c; c := c + 1;
    end if;
   end if;
  end for;
  tran := [Evaluate(akt,perm) : akt in tran];
 return(tran);
end function;

PointsOnLine := function(line, m)
 local tran,res,i,j,rj, akt;
 
 tran := NormalizeLine(line, Parent(line[1]), m);
 res := []; 
 for i := 0 to 1 do
  if i eq 0 then rj := [1]; else rj := BaseRing(Parent(line[1])); end if;
  for j in rj do
   Append(~res,[BaseRing(Parent(line[1]))!Evaluate(akt,[0,0,i,j]) : akt in tran]);
  end for;
 end for;
 return res; 
end function;

CritPointsOnLine := function(line,f,p,m,h)
 local r2, akt, tran, rf, fac, res, r2p, gcd, tmp, i;

 
 r2 := PolynomialRing(BaseRing(Parent(f)), 2);
 tran := NormalizeLine(line, PolynomialRing(BaseRing(Parent(f)),4), m);
 tran := [Evaluate(akt,[0,0,r2.1,r2.2]) : akt in tran]; 
 rf := Evaluate(f,tran);
 gcd := NumberFieldGCD(Coefficients(rf));
 try
  rf := rf div p;
 catch e
  "p :=", p;
  "gcd :=", gcd;
  "p should divide gcd";
  return false, 0,0;
 end try;  /* Equation is zero mod p but not mod p^2 on the line */
 try
  rf := rf div p;
  "p^2 should not divide gcd";
  return false, 0,0;
 catch e
  vprintf MinRedCubSurf,2:"Success\n";
 end try;

Fp := Codomain(m);
Fpx := PolynomialRing(Fp, 2);
if IsPrime(#Fp) then
    inc := hom<Fp -> Fpx | >;
else
    inc := hom<Fp -> Fpx | Fp.1>;
end if;
h := hom<r2 -> Fpx | m*inc,[Fpx.i : i in [1..2]]>;

r2p := Codomain(h);

 fac := Factorization(h(rf)); 
 res := [];
 for akt in fac do
  if Degree(akt[1]) eq 1 then
   tmp := [MonomialCoefficient(akt[1],r2p.2) @@ m, -MonomialCoefficient(akt[1],r2p.1) @@ m];
   Append(~res,[Evaluate(tran[i],tmp) : i in [1..4]]);
  end if;
 end for;
 return res;
end function;

/* The homogeneous case */
GroebnerBasisToPoint := function(gb) 
 local akt, mat, ker;

 if (#gb eq 3) and Max([TotalDegree(akt) : akt in gb]) eq 1 then
  mat := Matrix([[MonomialCoefficient(f,Parent(f).k) : k in [1..4]] : f in gb]);
  ker := BasisMatrix(Kernel(Transpose(mat)));
  assert Rank(ker) eq 1; /* Just a projective point */
  return true, [ker[1,k] : k in [1..4]];
 end if; 

 return false,[0,0,0,0];
end function;

ListIsolatedSingularPoints := function(f)
 local res,gen,i,sing,pd,j,suc,pt;
 
 res := []; 
 sing := ideal<Parent(f) | [f] cat [Derivative(f,i) : i in [1..4]]  /* gen cat [Parent(f).i - 1]*/>;
 pd := PrimaryDecomposition(sing);
 for j := 1 to #pd do 
  suc,pt := GroebnerBasisToPoint(GroebnerBasis(Radical(pd[j])));
  if suc then Append(~res,pt); end if;
 end for;  
 return(res);
end function;

/* Compute a unimodlar substitution, that evaluates on [1,0,0,0] to pt */
NormalizePoint := function(pt,rk,m)
 local subs, i, j, akt, perm;

 i := 1; 
 while pt[i] eq 0 do i := i + 1; end while;
 subs := [ ((pt[1]*pt[i]^(-1)) @@ m)*rk.i,((pt[2]*pt[i]^(-1)) @@ m)*rk.i,
           ((pt[3]*pt[i]^(-1)) @@ m)*rk.i,((pt[4]*pt[i]^(-1)) @@ m)*rk.i ];
 for j := 1 to 4 do if j ne i then subs[j] := subs[j] + rk.j; end if; end for;
 if i ne 1 then
  perm := [rk.1, rk.2, rk.3, rk.4];
  perm[i] := rk.1;
  perm[1] := rk.i;
  subs := [Evaluate(akt,perm) : akt in subs];
 end if;
 return subs;
end function;

ListCriticalPoints :=  function(f, lines, p, m, h)
 local crit_pnt,rp3, gb,pt,akt, suc;
 
 crit_pnt := ListIsolatedSingularPoints(h(f));
 crit_pnt := [c @@ m : c in crit_pnt];
 for line in lines do
  if Characteristic(Codomain(m)) le 3 then
   crit_pnt := crit_pnt cat PointsOnLine(line,m);
  else
   crit_pnt := crit_pnt cat CritPointsOnLine(line,f,p,m,h);
  end if;
 end for;
 crit_pnt := [m(c) : c in crit_pnt];

// if IsIrreducible(h(f)) and (not IsIrreducible(PolynomialRing(GF(p^3),4)!f)) then  // we don't do it here
/* We have 3 conjugate planes */
/*  rp3 := PolynomialRing(GF(p^3),4);
  fac := Factorization(rp3!f);
  gb := GroebnerBasis(ideal<rp3 | [akt[1] : akt in fac]>);
  suc, pt := GroebnerBasisToPoint(gb);
  if suc then
   if &and [akt in GF(p) : akt in pt] then
    Append(~crit_pnt, [GF(p)!akt : akt in pt]);
   end if;
  end if;
 end if;*/
 return crit_pnt;
end function;

/* f is an intermediate result:
   - treated with the wight-vector (0,1,1,1)
   - divided by p^2  

  We have to try, weather (0,0,1,1) gives us a gcd of p^2 or
                          (0,1,1,2) gives us a gcd of p^4.
*/
TryToComplete := function(f,p,m,h) 
 local fac, rp, l_fac, q_fac, tmp, ttc_2, gcd_2, sing, gb, i, tran, res, gcd, triv, tran_2, ttc, suc; 

 vprintf MinRedCubSurf,2:"Try to complete for weight vectors (0,1,2,2) and (0,2,2,3).\n";   

/* As we can go back to the beginning form ttc_1 with the weight-vector (0,0,0,1), we must be reducible */
 rp := Codomain(h);
 fac := Factorization(h(f));
 triv := [tmp : tmp in fac | tmp[1] eq rp.1];
 assert #triv gt 0; /* Using the trivial factor we get back to the beginning */
 if triv[1][2] ge 2 then
 /* printf"Trivial factor is multiple\n";  */
  return false,0,0; /* The trivial factor must have multiplicity 1, or no improvement is possible */
 end if;  
 
 l_fac := [tmp[1] : tmp in fac | (TotalDegree(tmp[1]) eq 1) and (tmp[1] ne rp.1)];
 q_fac := [tmp[1] : tmp in fac | TotalDegree(tmp[1]) eq 2]; 

 l_fac := [tmp : tmp in l_fac | tmp ne rp.1]; /* Throw it away */
 if (#q_fac gt 0) or (#l_fac eq 2) then /* Only the (0,1,2,2)-case is possible */
  q_fac := &*(q_fac cat l_fac);
  assert TotalDegree(q_fac) eq 2; /* The remaining degree 2 factor */
  sing := ideal<rp | [q_fac] cat [Derivative(q_fac,i) : i in [1..4]]>;
  gb := GroebnerBasis(Radical(sing));
 
  if (#gb eq 2) and (Max([Degree(i) : i in gb]) eq 1) then
   tran := NormalizeLine(gb,Parent(f),m);
   tran := [Evaluate(Evaluate(tmp,1,p*Parent(f).1),2,p*Parent(f).2) : tmp in tran];
   res := Evaluate(f,tran);
   gcd := NumberFieldGCD(Coefficients(res));
   
   try  /* Equation gets divisible by the transformation */
    _ := gcd div p;
   catch e
    "p :=", p;
    "gcd :=", gcd;
    "p should divide gcd";
     return false, 0,0;
    end try; /* Singular line is on the surface */
/*   printf "gcd %o found\n",gcd; */
   try  /* Equation gets divisible by the transformation */
    _ := gcd div p^2;
    vprintf MinRedCubSurf,2:"Success\n";
    return true, res div gcd, tran;
   catch e
    vprintf MinRedCubSurf,2:"No success\n";
   end try;
  end if;
  return false,0,0; /* No improvement possible */ 
 end if;

/* printf"Multiple linear factor found %o\n",l_fac;  */
  
/* Now we have a remaining linear factor of multiplicity 2 */ 
 tran := SubsForLFac(l_fac[1], Parent(f), p, m, h);
  
 ttc := Evaluate(f,tran);
 gcd := NumberFieldGCD(Coefficients(ttc));
 
 try  /* Equation gets divisible by the transformation */
  _ := gcd div p;
 catch e
  "p :=", p;
  "gcd :=", gcd;
  "p should divide gcd";
  return false, 0,0;
 end try; /* Singular line is on the surface */
/*   printf "gcd %o found\n",gcd; */
 try  /* Equation gets divisible by the transformation */
  _ := gcd div p^2;
  vprintf MinRedCubSurf,2:"Success\n";
  return true, ttc div gcd, tran;
 catch e
  vprintf MinRedCubSurf,2:"No success\n";
 end try;
 ttc := ttc div p;
 suc, res, tran_2 := TestRedMinimization(ttc,p,m,h);
 if suc then
  vprintf MinRedCubSurf,2:"Success\n";
  return true,res, [Evaluate(tmp,tran_2) : tmp in tran];
 end if; 
 
/* Now we are finished with the wight vector (0,1,2,2). 
   (0,2,2,3) is only possible if ttc is modulo p the union of 3 conjugate planes that have a common line. 
   If one plane is definied over GF(p) then we are already done.   */

/* printf "Intermediate equation %o\n",ttc; */
 sing := ideal<rp | [ttc] cat [Derivative(ttc,i) : i in [1..4]]>;
 gb := GroebnerBasis(Radical(sing));

 if (#gb eq 2) and (Max([Degree(i) : i in gb]) eq 1) then
  tran_2 := NormalizeLine(gb,Parent(f),m);
  tran_2 := [Evaluate(Evaluate(tmp,1,p*Parent(f).1),2,p*Parent(f).2) : tmp in tran_2];
  
  ttc_2 := Evaluate(ttc,tran_2);
  gcd := NumberFieldGCD(Coefficients(ttc_2));
  try  /* Equation gets divisible by the transformation */
   _ := gcd div p;
  catch e
   "p :=", p;
   "gcd :=", gcd;
   "p should divide gcd";
   return false, 0,0;
  end try; /* Singular line is on the surface */
/*   printf "gcd %o found\n",gcd; */
  try  /* Equation gets divisible by the transformation */
   _ := gcd div p^2;
   vprintf MinRedCubSurf,2:"Success\n";
   return true, res div gcd, [Evaluate(tmp,tran_2) : tmp in tran];
  catch e
   vprintf MinRedCubSurf,2:"No success\n";
  end try;
  ttc_2 := ttc_2 div p;
  fac := Factorization(rp!ttc_2);
  for tmp in fac do
   if TotalDegree(tmp[1]) eq 1 then
    tran_3 := SubsForLFac(tmp[1], Parent(f), p, m, h); 
    ttc_3 := Evaluate(ttc_2,tran_3);
    gcd := NumberFieldGCD(Coefficients(ttc_3));
    try  /* Equation gets divisible by the transformation */
     _ := gcd div p;
    catch e
     "p :=", p;
     "gcd :=", gcd;
     "p should divide gcd";
     return false, 0,0;
    end try; /* Singular line is on the surface */
/*   printf "gcd %o found\n",gcd; */
    try  /* Equation gets divisible by the transformation */
     _ := gcd div p^2;
     tran := [Evaluate(tmp,tran_2) : tmp in tran];
     vprintf MinRedCubSurf,2:"Success\n";
     return true, ttc_3 div gcd, [Evaluate(tmp,tran_3) : tmp in tran];
    catch e
     vprintf MinRedCubSurf,2:"No success\n";
    end try;
   end if; 
  end for; 
 end if;   
 return false,0,0;
end function;

/* Returns a better equation and the applied substitution */
LocalMinimization := function(f,p,m,h) 
 local suc,res,tran, lines, gcd, crit_pnt, ttc, tran_2, i;

 vprintf MinRedCubSurf,2:"Trying reducible reduction\n";
 suc,res,tran := TestRedMinimization(f,p,m,h);
 if suc then 
  vprintf MinRedCubSurf,2:"Success using weight vector (0,0,0,1). \n";
  return res, tran; 
 end if;
/* Now we assume that the reduction of the surface is irreducible over GF(p)
   only singular points and lines are possible */

/* Weight vector 0,0,1,1 -- meaning reduction with singular line -- test it */
/* printf "Trying (0,0,1,1)\n"; */
 vprintf MinRedCubSurf,2:"Searching for singluar lines\n";
 lines := ListSingularLines(h(f));
 vprintf MinRedCubSurf,2:"%o singular lines found.\n",#lines;
 for line in lines do 
  tran := NormalizeLine(line, Parent(f),m);    
  tran := [Evaluate(akt,1,p*Parent(akt).1) : akt in tran];
  tran := [Evaluate(akt,2,p*Parent(akt).2) : akt in tran];
  res := Evaluate(f,tran);
  gcd := NumberFieldGCD(Coefficients(res));
  try 
   _ := gcd div (p^2);
   vprintf MinRedCubSurf,2:"Success using weight vector (0,0,1,1). \n";   
   res := res div gcd;
   return res,tran;
  catch e
   vprintf MinRedCubSurf,2:"No success using weight vector (0,0,1,1). \n";   
  end try;
 end for;
/* when we arrive here, then the surface is irreducible and no singular line is mod p^2 contained in f = 0 */

 vprintf MinRedCubSurf,2:"Searching for critical points. \n";
 crit_pnt := ListCriticalPoints(f,lines,p, m, h);
 vprintf MinRedCubSurf,2:"%o critical points found. \n",#crit_pnt;
/* crit_pnt; */
/* printf"Trying remaining weight-vectors\n"; */
 for i := 1 to #crit_pnt do
  time tran := NormalizePoint(crit_pnt[i], Parent(f), m);
/*  printf "Normalize singular point by %o\n",tran;  */
  tran := [Evaluate(akt,[Parent(f).1, Parent(f).2*p, Parent(f).3*p, Parent(f).4*p]) : akt in tran];
  time res := Evaluate(f, tran);
  time gcd := NumberFieldGCD(Coefficients(res));
/*  printf "gcd %o found\n",gcd; */
  try 
   time _ := gcd div (p^3);    
   vprintf MinRedCubSurf,2:"Success using weight vector (0,1,1,1). \n";   
   return (res div gcd),tran;
  catch e
   vprintf MinRedCubSurf,2:"No success using weight vector (0,1,1,1). \n";         
  end try;
  try 
   _ := gcd div (p^2);
   suc,ttc,tran_2 := TryToComplete(res div gcd,p);
   if suc then 
    return ttc, [Evaluate(akt,tran_2) : akt in tran];
   end if;
  catch e
   vprintf MinRedCubSurf,2:"No success. \n"; 
  end try;
 end for; 
 return f, [Parent(f).i : i in [1..4]];
end function;



function MinimizeCubicSurfaceNbFld(f, d)
    OK := BaseRing(Parent(f));
    S := Parent(f);
    _, p := IsPrincipal(d[1]);
    Fp, m := ResidueClassField(OK, d[1]);
    Fpx := PolynomialRing(Fp, 4);
    if IsPrime(#Fp) then
        inc := hom<Fp -> Fpx | >;
    else
        inc := hom<Fp -> Fpx | Fp.1>;
    end if;
    h := hom<S -> Fpx | m*inc,[Fpx.i : i in [1..4]]>;
    res := ElimGCD(f);
    tran := [Parent(f).i : i in [1..4]];

    if Characteristic(Fp) ge 5 then 
        repeat
        res,t_neu := LocalMinimization(res,p,m,h);
        tran := [Evaluate(akt,t_neu) : akt in tran];
        until (t_neu eq [Parent(f).i : i in [1..4]]);
    end if;
    /* Simplify transformation by using LLL */
    //Lat := LLL(Matrix([[MonomialCoefficient(tran[j],Parent(f).i) : j in [1..#tran]] : i in [1..Rank(Parent(f))]]));
    //Lat := Transpose(Lat);

    //res := f^Lat;
    res := ElimGCD(res);
    
    return res, tran;//Lat;
end function;

Q1, E1 := MinimizeQuadric(Q,E);


P, t0, r := NewBasis(Parent(Q)!Q1);
f0 := ChangeOfBasis(Parent(Q)!E1,P);

if t0 eq 3 then // compute an invariant faster than the discriminant, to find the places of potentially bad reduction
	R<s, t, w> := PolynomialRing(BaseRing(Parent(f0)), [1,1,2]);
	f_weighted := Evaluate(f0, [s^2, s*t, t^2, w]);

	alpha := MonomialCoefficient(f_weighted, w^3);
	f_weighted /:= alpha;        
	f_weighted := Evaluate(f_weighted, [s, t, w-ExactQuotient(Terms(f_weighted, w)[3], 3*w^2)]);

	S<[x]> := PolynomialRing(BaseRing(Parent(f_weighted)), 2);

	f := S!Evaluate(f_weighted, [x[1], x[2], 0]);
	v := S!Evaluate(ExactQuotient(Terms(f_weighted, w)[2], w), [x[1], x[2], 0]);

	h24 := Transvectant(f, f, 4);
	inv := Evaluate(Transvectant(h24, v, 4), [0,0])*alpha^8/Determinant(P)^(24);
    if inv eq 0 then
        inv := InvariantsGenus4Curves(Parent(Q)!Q1, Parent(Q)!E1)[1];//Evaluate(Transvectant(f, f, 6), [0,0])*(alpha^6/(Determinant(P)^(18)));
    end if;
else 
	R<x, y, u, v> := PolynomialRing(BaseRing(Parent(f0)), 4);
	f_bic := 1/Determinant(P)^3*Evaluate(f0, [x*u, y*u, x*v, y*v]);
		
	inv := Transvectant(f_bic, f_bic, 3, 3 : invariant := true);
end if;
if inv ne 0 then
	l_p := Decomposition(OK!inv);
end if;

SetDebugOnError(true);
SetVerbose("MinRedCubSurf", 2);
for i->d in l_p do
    E1, T := MinimizeCubicSurfaceNbFld(E1, d);
    Q1 := ElimGCD(Evaluate(Q1, T));
end for;

/* the curves are still isomorphic, no problem of transposition with the transformations
inv := InvariantsGenus4Curves(Parent(Q)!Q1,Parent(Q)!E1 : normalize := true);
inv1 := InvariantsGenus4Curves(Parent(Q)!Q,Parent(Q)!E : normalize := true);*/