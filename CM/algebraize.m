function AlgebraizedInvariants(ICC, K : Base := false)
    if Base then
        "Trying to algebraize invariants...";
        test, I := AlgebraizeElementsExtra(ICC, K);    
        if not test then
            print "Failed to algebraize invariants.";
            return 0, 0, false;
        end if;
        L := K; 
        hKL := CanonicalInclusionMap(K, L);
    else
        L, I, hKL := NumberFieldExtra(ICC, K);
    end if;
    return I, hKL, true;
end function;

function AlgRecG4(inv::SeqEnum)
    prec := Precision(Parent(inv[1]));
    if #inv eq 60 then 
        wgt := [ 6, 12, 18, 30, 45, 4, 6, 8, 10, 10, 9, 13, 14, 12, 12, 15, 14, 14, 15, 15, 17, 16, 16, 20, 19, 19, 18, 16, 18, 18, 17, 17, 21, 19, 22, 22, 23, 21, 20, 20, 21, 25, 26, 24, 21, 23, 23, 24, 25, 29, 25, 28, 27, 32, 29, 31, 35, 33, 37, 41 ];
        inv := WPSMultiply(wgt, inv, (1/(inv[1]^(1/6)))); // doable iff inv[1] is not 0, add later
        inv1 := AlgebraizedInvariants(inv, RationalsExtra(Floor(1/3*prec)));
    elif #inv eq 65 then
        wgt := [ 2, 4, 4, 6, 6, 8, 8, 10, 12, 14, 6, 8, 8, 10, 10, 10, 10, 10, 10, 12, 12, 12, 12, 12, 12, 12, 12, 12, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18 ];
        inv := WPSMultiply(wgt, inv, (1/inv[1])^(1/2)); // doable iff inv[1] is not 0, do later
        inv1 := AlgebraizedInvariants(inv, RationalsExtra(prec));
    end if;
    Q1, E1 := ReconstructionGenus4(inv1);
    return Q1, E1;
end function;
