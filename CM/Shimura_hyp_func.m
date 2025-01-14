function GetCurveShimura(param : prec := prec)
    K<mu> := Parent(param);

    M<mu0> := PolynomialRing(K);
    S<X,Y> := PolynomialRing(FieldOfFractions(M), 2);

    f := (-1432353454804332119596592722603696979968*mu0 - 4886589248453472)*X^10 - 309519*X^9*Y + 9851363924882199552*mu0*X^8*Y^2 + (-104516331448118333562141904207872*mu0^2 + 118855296*mu0)*X^7*Y^3 - 2837192810366073470976*mu0^2*X^6*Y^4 - 17115162624*mu0^2*X^5*Y^5 - 90790169931714351071232*mu0^3*X^4*Y^6 - 1408333381632*mu0^3*X^3*Y^7 - 26288889790464*mu0^4*X*Y^9;
    Coeff, Mon := CoefficientsAndMonomials(f);
    return Polynomial([Evaluate(c, param) : c in Coeff], Mon);
end function;