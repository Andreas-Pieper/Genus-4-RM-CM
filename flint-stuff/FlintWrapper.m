// You first need to compile acb_theta.c

function mktemp()
    return Split(Pipe("mktemp", ""), "\n")[1];
end function;

function FlintToMagma(s)
	sSplit := Split(s, "/");
	sNumber := sSplit[1][[1..#sSplit[1]-1]];
	sError := sSplit[2][[2..#sSplit[2]]];
	sNumber := ReplaceCharacter(sNumber, "j", "*j");
	sNumber := ReplaceCharacter(sNumber, "e", "E");
    sNumber := StripWhiteSpace(sNumber);
	sError := ReplaceCharacter(sError, "j", "*1");
	sError := ReplaceCharacter(sError, "e", "E");
	mError := Max([ Abs(eval s) : s in Split(sError, ",() ") ]);
	mNumber := Max([1] cat [ Abs(eval s) : s in  Split(sNumber, " *,()j") | s ne "+" and s ne "-"]);
	prec := Ceiling(Log(10, mNumber/mError)) + 3;
	if prec le 0 then
		return "precision error";
	end if;
    sNumber cat:=Sprintf("p%o", prec);
	if sNumber eq "nanp3" then
		return "infinity error";
	end if;
	return eval sNumber;
end function;

intrinsic ThetaFlint(char::Mtrx, z::Mtrx, tau::Mtrx[FldCom]) -> SeqEnum
    {}
	assert IsZero(char);
	arb_print_real := func<elt | ReplaceCharacter(Sprintf("%o +/- %.*o\n", elt, 3, Max(Abs(elt)*10^(1-Precision(Parent(elt))), 10^(1-Precision(Parent(elt)))) ), "E", "e")>;
	arb_print_complex := func<elt | arb_print_real(Real(elt)) cat arb_print_real(Imaginary(elt))>;
	arb_print_matrix := func<elt | &cat[ arb_print_complex(elt[i,j]) : i in [1..NumberOfRows(elt)], j in [1..NumberOfColumns(elt)] ] >;
	i := Random(10^20);
    input_filename := mktemp();
    output_filename := mktemp();
	Write(input_filename, arb_print_matrix(z) cat arb_print_matrix(tau));
    digits := Precision(BaseRing(Parent(tau)));
    bits_precision := Ceiling(Log(10)*digits/Log(2));
    cmd := Sprintf("~/github/Genus-4-RM-CM/flint-stuff/a.out %o %o %o %o 0", NumberOfRows(tau), bits_precision+100, input_filename, output_filename);
	
    _ := System(cmd);
	output := Read(output_filename);
	output2 := Split(output, "\n");
	reals_list := [* FlintToMagma(s) : s in output2 *];
	err := [el : el in reals_list | Type(el) eq MonStgElt];
	if #err ge 1 then
		Set(err);
		return "error";
	end if;
	CC<j> := ComplexField(Min([Precision(r) : r in reals_list]));
	return [ reals_list[i] + reals_list[i+1]*j : i in [1..#reals_list by 2] ];
end intrinsic;