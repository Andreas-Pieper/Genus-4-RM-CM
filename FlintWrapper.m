// wrapper that ports Jean Kieffer's Flint code for computing theta functions to Magma
// written by Edgar Costa and Raymond van Bommel
// requires CHIMP

function FlintToMagma(s)
	sSplit := Split(s, "/");
	sNumber := sSplit[1][[1..#sSplit[1]-1]];
	sError := sSplit[2][[2..#sSplit[2]]];
	sNumber := ReplaceCharacter(sNumber, "j", "*j");
	sNumber := ReplaceCharacter(sNumber, "e", "E");
	sError := ReplaceCharacter(sError, "j", "*1");
	sError := ReplaceCharacter(sError, "e", "E");
	mError := Max([ Abs(eval s) : s in Split(sError, ",() ") ]);
	mNumber := Max([1] cat [ Abs(eval s) : s in  Split(sNumber, " *,()j") | s ne "+" and s ne "-"]);
	prec := 2*Ceiling(Log(10, mNumber/mError)) + 3;
	CC<j> := ComplexField(prec);
	return eval sNumber;
end function;

function ThetaFlint(char, z, tau)
	assert IsZero(char);
	arb_print_real := func<elt | ReplaceCharacter(Sprintf("%o +/- %.*o\n", elt, 3, Abs(elt)*10^(1-Precision(Parent(elt)))), "E", "e")>;
	arb_print_complex := func<elt | arb_print_real(Real(elt)) cat arb_print_real(Imaginary(elt))>;
	arb_print_matrix := func<elt | &cat[ arb_print_complex(elt[i,j]) : i in [1..NumberOfRows(elt)], j in [1..NumberOfColumns(elt)] ] >;
	i := Random(10^20);
	//return arb_print_matrix(z) cat arb_print_matrix(tau);
	Write("input" cat Sprint(i), arb_print_matrix(z) cat arb_print_matrix(tau));
	cmd := "~/flint/build/examples/acb_theta " cat Sprint(NumberOfRows(tau)) cat " " cat Sprint(Precision(BaseRing(Parent(tau)))) cat " ~/CODE/input" cat Sprint(i) cat " /dev/null 1";
	output := Pipe(cmd, "");
	System("rm ~/CODE/input" cat Sprint(i));
	output2 := Split(output, ":");
	output3 := output2[#output2];
	return [FlintToMagma(s) : s in Split(output3, "\n")];
end function;
