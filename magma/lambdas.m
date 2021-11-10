load "Discriminants.m";
load "DiagonalRestrictions.m";
load "pAdicL.m";
/* import "/home/havard/Projects/HilbertEisenstein/magma/code/Test/pAdicL.m": pAdicLfunction; */

function getLambda(F,p : m :=10)
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
d := 3;				/* Q(sqrt(d)) */
Discs := [**];
for D := 2 to 100 do
    if (not IsSquare(d*D)) and (d*D eq FundamentalDiscriminant(d*D)) then
	Append(~Discs, d*D);
    end if;
end for;
D := Discs[1];
Lambdas := [**];

for D in Discs do
    F := QuadraticForms(D)!1;
    F;
    Symb := KroneckerSymbol(D,p);
    /* if not (Symb eq 0) then */
    P,Q, lambda := getLambda(F,p);
    Append(~Lambdas,  [P,Q,lambda,Symb]);
    /* end if; */
end for;
/* Discs; */




