/*

November 1, 2018

Looking at Jan's example.


*/

D:=13;
f:=5;
K:=QuadraticField(D);
O:=Integers(K);
Fs:=Factorisation(ideal<O|f>);
ff:=Fs[1][1];
G,m:=RayClassGroup(ff);

function computeXn(D,n) // no condition on p

	K:=QuadraticField(D);
    assert ClassNumber(K) eq 1;
    O:=Integers(K);
    Qx<x>:=PolynomialRing(Rationals());
    
    Bs:={i: i in [0 .. Ceiling(n*SquareRoot(D)) - 1]}; 
    Bs:=Bs join {-B: B in Bs};
    
    delta:=Roots(x^2 - D,K)[1][1];
    
    Xn:=[**];
    for B in Bs do
    	alpha:=(-B + n*delta)/(2*delta);
    	assert (B^2 - n^2*D) lt 0;
        assert Trace(alpha) eq n;
        if (B mod 2) eq (n*D mod 2) then
    		alphadelta:=ideal<O|alpha*delta>;
            divs:=Divisors(alphadelta);
            for I in divs do
                Append(~Xn,[*alpha,I*]);
            end for;
    	end if;
    end for;
    
    return Xn,delta;

end function;


function jan(NN) // weight 2k

	k:=2;
    Z:=Integers();
	Zq<q>:=PowerSeriesRing(Z,NN);

	D:=13;
	f:=5;
	K:=QuadraticField(D);
	O:=Integers(K);
	Fs:=Factorisation(ideal<O|f>);
	ff:=Fs[1][1];
	G,m:=RayClassGroup(ff);
    
    minv:=Inverse(m);

    Gk:=Zq!0;
    O1:=ideal<O|1>;
    
    for n:=1 to NN-1 do
    	Gkn:=0;
        Xn:=computeXn(D,n);
    	for alphaI in Xn do
        	I:=alphaI[2];
            gcd:=GreatestCommonDivisor(I,ff);
            if gcd eq O1 then
            	J:=minv(I);
                // J, J eq G.1, Order(J);
            	if J eq G.1 then
        			Gkn:=Gkn  - Norm(I)^(k-1); 
            	else
                	Gkn:=Gkn  + Norm(I)^(k-1);
            	end if;
            end if;
        end for;
        Gk:=Gk + 4*Gkn*q^n;
    end for;

	return Gk;

end function;
