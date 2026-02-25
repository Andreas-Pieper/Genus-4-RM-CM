d:= 4;
// test
function IsGaloisGroupOfCMField(G)
	if IsPrimitive(G) then return false; end if;
	X := MinimalPartition(G);
	if #X ne d then return false; end if;
	if Sym(2*d)![[Random(x diff {y}): x in X | y in x][1]: y in [1..2*d]] notin G then return false; end if;	
	return true;
end function;

function IsPrimitiveType(G, Phi)
	S := {g: g in G| 1^g in Phi};
	H := [h: h in G| {h*s: s in S} eq S];
	return #G/#H eq 2*d;
end function;


function ReflexType(G, Phi)
	S := {g: g in G| 1^g in Phi};
        Hr := [h: h in G| {s*h: s in S} eq S];
	So := {{s*h: h in Hr}: s in S};
	Phir := {Setseq(o)[1]^(-1): o in So};
	return sub<G|Hr>, Phir;
end function;

// Given CM type Phi, identify which element of CM types
// it is conjugate to
function IdentifyCMType(G, Phi, CMtypes)
  orbit := Setseq({{x^g: x in Phi}: g in G});
  for Phi0 in CMtypes do
    if Phi0 in orbit then
      return Phi0;
    end if;
  end for;
end function;

// Storing exponent bounds in nested associative array
// each key is G; each value is another associative array whose keys are CM types Phi; each value of that is the exponent bound
Gs := TransitiveGroups(2*d, IsGaloisGroupOfCMField);
exps := AssociativeArray();
for G in Gs do
    expsG := AssociativeArray();
	A := GroupAlgebra(Rationals(), G);
	X := MinimalPartition(G);
	H := sub<G|[h:h in G| 1^h eq 1]>;
	T,m := Transversal(G, H);
	sig := G![[Random(x diff {y}): x in X | y in x][1]: y in [1..2*d]];
	CMtypes := CartesianProduct([X[i]: i in [1..d]]);
	cmtypeorbits := {};
	for Phi in CMtypes do
		orbit := {{x^g: x in Phi}: g in G};
	        Include(~cmtypeorbits, orbit);
		if #CMtypes eq &+[#o: o in cmtypeorbits] then break Phi; end if;
	end for;
	CMtypespri := [Setseq(o)[1]: o in cmtypeorbits| IsPrimitiveType(G, Setseq(o)[1])];
	for Phi in CMtypespri do
		Hr, Phir := ReflexType(G, Phi);
		refl := &+[A!phi: phi in Phir];
		refln := &+[A![g: g in G|1^g eq phi][1]: phi in Phi];
		Mr  := PermutationModule(G, Hr, Rationals());
		M := PermutationModule(G, Rationals());
		T := [t: t in T];
		T0 := [[i: i in [1..#T]| M.i eq M.1*T[j]][1]: j in [1..#T]];
		ParallelSort(~T0 ,~T);
		assert Basis(M) eq [M.1*T[i]: i in [1..#T]];
		Mat := Matrix([Eltseq(&+[Mr.1*phi*T[i]: phi in Phir] ): i in [1..#T]]);
		f := hom<M->Mr|Mat>;
		IsModuleHomomorphism(Hom(M,Mr)!f);
        	print G, Phi;
        // Want {phi: Mr -> M such that phi*f*iota = iota}
        // where iota: Msig -> M
        V := GHom(Mr,M);
        B := Basis(V);
        Bint := [];
        for b in B do
          bint := Denominator(b)*b;
          Append(~Bint, bint);
        end for;
        Msig := sub< M | [b*sig - b : b in Basis(M)] >;
        Mrsig := sub< Mr | [b*sig - b : b in Basis(Mr)] >;
        // create inclusion map Mrsig -> Mr
        iota := Matrix([Coordinates(M,b) : b in Basis(Msig)]);
        rows := [Eltseq(iota*Mat*b) : b in Bint];
        Append(~rows, Eltseq(iota));
        Rel := Matrix(rows);
        // make matrix integral
        Rel *:= Denominator(Rel);
        Rel := ChangeRing(Rel,Integers());
        k := KernelMatrix(Rel);
        exp := GCD(Eltseq(Rows(Transpose(k))[Ncols(k)]));
        expsG[Phi] := exp;
	end for;
    exps[G] := expsG;
end for;
