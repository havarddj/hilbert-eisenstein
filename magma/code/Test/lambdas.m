/* load "Discriminants.m"; */
/* load "DiagonalRestrictions.m"; */
load "pAdicL.m";
/* import "/home/havard/Projects/HilbertEisenstein/magma/code/Test/pAdicL.m": pAdicLfunction; */

function getLambda(F,p : m :=8)
    P,Q := pAdicLfunction(F,p,m);
    lambda := 0;
    L := Coefficients(Q);
    n := 1;		
    while((ZZ!L[n] mod p) eq 0) do
	lambda := lambda + 1;
	n := n + 1;
    end while;
    return P, Q, lambda;
end function;

p := 5;
d := 7;				/* Q(sqrt(d)) */
Discs := [**];
for f := 2 to 20 do
    /* if (not IsSquare(d*D)) and (d*D eq FundamentalDiscriminant(d*D)) then */
    /* 	Append(~Discs, d*D); */
    /* end if; */
    D := f*f*d;
    if (((D mod 4) eq 0) or ((D mod 4) eq 0 )) and not (f eq p) then
	Append(~Discs, D);
    end if;
end for;

Lambdas := [**];

for D in Discs do
    for F in ReducedForms(D) do
	/* F := QuadraticForms(D)!1; */
	F;
	print "D = ", D, "\n F = ", F;
	Symb := KroneckerSymbol(D,p);
	/* if not (Symb eq 0) then */
	P,Q, lambda := getLambda(F,p);
	Append(~Lambdas,  [*F, D, P, Q, ZZ!lambda, ZZ!Symb*]);
    end for;
    /* end if; */
end for;
/* Discs; */




