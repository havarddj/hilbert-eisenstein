/*

November 6, 2018

Looking for an example using wide class group order 2, with modulus
a single split prime.

*/



function splitsearch(B,pB: narrow:=false)

	out:=[];
	for b:=2 to B do
    	if IsFundamentalDiscriminant(b) eq true then
        	K:=QuadraticField(b);
            O:=Integers(K);
            for p in PrimesUpTo(pB) do
    			if b mod p ne 0 then
                	Fs:=Factorisation(ideal<O|p>);
                    if #Fs eq 2 then
                    	ff:=Fs[1][1];
                        if narrow eq false then
                        	G,m:=RayClassGroup(ff);
                        else
                        	G,m:=RayClassGroup(ff,[1,2]);
                        end if;
                        Append(~out,[*b,p,Order(G)*]);
                	end if;
                end if;
    		end for;
    	end if;
    end for;

	return out;


end function;

load "./Magma-Files/LV/janexample.m";

function split(NN) // weight 2k

	k:=2;
    Z:=Integers();
	Zq<q>:=PowerSeriesRing(Z,NN);

	D:=29;
	f:=5;
	K:=QuadraticField(D);
	O:=Integers(K);
	Fs:=Factorisation(ideal<O|f>);
	ff:=Fs[1][1];
    Norm(ff);
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


f:=split(12);
f;
M:=ModularForms(Gamma1(5),4);
B:=Basis(M);
4*B[1] + 4*B[2] - 28*B[3] - 104*B[4] + 228*B[5];

// Now look at example with narrow class group

function splitnarrow(NN) // weight 2k

	k:=2;
    Z:=Integers();
	Zq<q>:=PowerSeriesRing(Z,NN);

	D:=89;
	f:=5;
	K:=QuadraticField(D);
	O:=Integers(K);
	Fs:=Factorisation(ideal<O|f>);
	ff:=Fs[1][1];
    Norm(ff);
	G,m:=RayClassGroup(ff,[1,2]); // order 4
    
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
            	if (J eq 0*G.1) or (J eq 2*G.1) then
        			Gkn:=Gkn  + Norm(I)^(k-1); 
            	else
                	Gkn:=Gkn  - Norm(I)^(k-1);
            	end if;
            end if;
        end for;
        Gk:=Gk + 4*Gkn*q^n;
    end for;

	return Gk;

end function;


f:=splitnarrow(12);
M:=ModularForms(Gamma1(5),4);
B:=Basis(M);
24*B[1] + 24*B[2] - 168*B[3] - 624*B[4] + 1368*B[5];
f;



