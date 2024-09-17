G := TransitiveGroup(9,11);
Order(G);
G := TransitiveGroup(6,11);
Order(G);
TransitiveGroupDescription(G);
subgroups := Subgroups(G);
Hs := [el : el in subgroups | Order(el`group) eq 8];
subgroups[1];
Hs := [el : el in subgroups | Order(el`subgroup) eq 8];
G;
g := G!(1,5)(2,4)(3,6);
g;
[el : el in H | not g in el];
[el : el in Hs | not g in el];
[el : el in Hs | not g in el`subgroup];
[el`subgroup : el in $1];
Hs := $1;
#Hs;
Subgroups;
Hs0 := [el : el in subgroups | Order(el`subgroup) eq 8];
#Hs0;
#Hs;
g in Center(G);
g;
#Center(G);
Center(G);
g := Generators($1)[1];
Generators(Center(G));
Random($1);
g := $1;
Hs := [el : el in Hs0 | no g in el`subgroup];
Hs := [el : el in Hs0 | not g in el`subgroup];
#Hs;
Hs;
Q, mpQ := quo< G | Center(G) >;
#Q;
Q;
Q eq Sym(4);
G2 := TransitiveGroup(8,24);
#G2;
G2;
Center(G2);
g := Random(Generators($1));
Hs := [el : el in Subgroups(G2) | Order(el`subgroup) eq 8 and (not g in el`subgroup)];
#Hs;
Hs;
Hs_old := Hs;
Hs := [el`subgroup | el in Hs_old];
Hs := [el`subgroup : el in Hs_old];
Hs;
Orbit;
H := Hs[1];
Orbit(H,1);
Orbit(H,2);
H;
Center(H);
g;
[Orbit(H, el) : el in [1..8]];
Seqset($1);
for H in Hs do
for H in Hs do
print H;
{Orbit(H, el) : el in [1..8]};
print "-------------------"
;
end for;
g
;
Stabilizer(H,1);
IsMaximal($1);
MaximalSubgroups;
H;
MaximalPartition(H);
BlocksAction;
G2;
Stabilizer(G2,1);
Order($1);
Conjugates;
G2;
S := Stabilizer(G2,1);
#S;
N := Normalizer(G2,S);
T := Transversal(G2,N);
#T;
T[1];
T[2]*S*(T[2])^-1;
T[2];
S;
gens := Generators(S);
gens[1]^T[2];
gens[1]^(T[2]);
gens;
Setseq(gens);
gens := $1;
gens[1]^(T[2]);
conjs := [];
for t in T do
C := sub< G2 | [g^t : g in gens] >;
Append(~conjs, C);
end for;
#conjs;
subgroups_containing := []l
subgroups_containing := [];
for Hrec in Subgroups(G2) do
H := Hrec`subgroup;
con_bool := false;
for c in C do
if c subset H then
con_bool := true;
end if;
end for;
if con_bool and (not g in H) then
Append(~subgroups_containing,H);
end if;
end for;
C := [sub< G2 | [g^t : g in gens] > : t in T];
#C;
C;
subgroups_containing := [];
for Hrec in Subgroups(G2) do
subgroups_containing := [];
for Hrec in Subgroups(G2) do
  H := Hrec`subgroup;
  con_bool := false;
  for c in C do
    if c subset H then
      con_bool := true;
    end if;
  end for;
  if con_bool and (not g in H) then
    Append(~subgroups_containing,H);
  end if;
end for;
#subgroups_containing;
subgroups_containing;
S;
GaloisGroup;
R<x> := PolynomialRing(QQ);
f := x^8 + 8*x^6 + 18*x^4 + 11*x^2 + 1;