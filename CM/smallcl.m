/***
 *  Better class group representatives
 *
 *  Written by Jeroen Sijsling (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
 */


function MakeSmallPeriodMatrix(P)
/* LLL on columns of a period matrix (does not seem to be very useful) */

PRR := VerticalSplitMatrix(P);
Lat := Lattice(Transpose(PRR));
Lat0, M := LLL(Lat);
P0RR := Matrix([ Eltseq(b) : b in Basis(Lat0) ]);
P0 := CombineVerticallySplitMatrix(Transpose(P0RR), BaseRing(P));
return P0, M;

end function;


intrinsic MakeSmallClassRepresentative(aa::., infs::. : precred := 300) -> .
{Makes an ideal class in the sense of Minkowski.}

/* Create corresponding lattice */
B := Basis(aa^(-1));
M := [ ];
for b in B do
    row := [ ];
    for inf in infs do
        ev := Evaluate(b, inf : Precision := precred);
        pair := [ Re(ev), Im(ev) ];
        row := row cat pair;
    end for;
    Append(~M, row);
end for;
M := Matrix(M);
Lat := Lattice(M);

/* Create small element in aa */
vs := ShortestVectors(Lat);
v := Vector(vs[1]);
c := NumericalSolution(M, v);
c := [ Round(x) : x in Eltseq(c) ];
alpha := &+[ c[i]*B[i] : i in [1..#B] ];
return alpha*aa, alpha;

end intrinsic;


function MakeSmallSequence(seq, gens, infs : precred := 300)
/* Shortcut use of previous function in the presence of generators */

assert #gens gt 0;
aa := gens[1]^0;
for i in [1..#seq] do
    /* Use binary representation */
    q := seq[i]; power := gens[i];
    while true do
        r := q mod 2;
        if r ne 0 then
            aa := MakeSmallClassRepresentative(aa*power, infs : precred := precred);
        end if;
        q := (q - r) div 2;
        if q eq 0 then
            break;
        end if;
        power := MakeSmallClassRepresentative(power^2, infs : precred := precred);
        B := Basis(power);
        power := ideal< Parent(B[1]) | B >;
    end while;
end for;
return aa;

end function;


intrinsic ClassGroupSmallRepresentatives(ZZK, exp : precred := 300) -> .
{Finds representatives of small norm of all elements of the class group. If an exponent is specified, only representatives for Cl / exp Cl are returned.}

/* Basic data */
Cl, phiCl := ClassGroup(ZZK);
vprint CMExp, 1 : "Done determining class group!", Cl;
infs := InfinitePlaces(NumberField(ZZK));
if #Cl eq 1 then return [ ideal< ZZK | 1 > ]; end if;

/* Throw away superfluous generators */
n := NumberOfGenerators(Cl); gens := [ Cl.i : i in [1..n] ];
if exp ne Infinity() then
    gcds := [ GCD(exp, Order(gens[i])) : i in [1..n] ];
    gens := [ gens[i] : i in [1..n] | gcds[i] ne 1 ];
    n := #gens;
    /* (Stupid recalculation of gcd that costs no time) */
    gcds := [ GCD(exp, Order(gens[i])) : i in [1..n] ];
    if #gcds eq 0 then return [ ideal< ZZK | 1 > ]; end if;
    CP := CartesianProduct([ [0..(gcds[i] - 1)] : i in [1..n] ]);
else
    CP := CartesianProduct([ [0..(Order(gens[i]) - 1)] : i in [1..n] ]);
end if;

/* Start with small generators and cycle */
gens := [ MakeSmallClassRepresentative(phiCl(gens[i]), infs : precred := precred) : i in [1..n] ];
tupold := 0; runningtupold := 0; aaold := 0; aas := [ ];
counter := 0;
for tup in CP do
    if &and[ c eq 0 : c in tup ] then
        runningtup := [ ideal< ZZK | 1 > : gen in gens ];
        Append(~aas, ideal< ZZK | 1 >);
    else
        i := Minimum([ j : j in [1..n] | tup[j] ne tupold[j] ]);
        runningtup[i] := MakeSmallClassRepresentative(runningtup[i]*gens[i], infs : precred := precred);
        for j in [(i + 1)..n] do runningtup[j] := ideal< ZZK | 1 >; end for;
        Append(~aas, MakeSmallClassRepresentative(&*runningtup, infs : precred := precred));
    end if;
    tupold := tup; runningtupold := runningtup; aaold := aas[#aas];
end for;
return aas;

end intrinsic;
