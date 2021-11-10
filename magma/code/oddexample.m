/*

November 6, 2018

Looking for an example using wide class group order 2, with modulus
a single split prime.

November 8, 2018

Now looking for totally odd character.

*/



function splitsearchodd(B,pB)

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
                        G1:=RayClassGroup(ff,[1]);
                        G2:=RayClassGroup(ff,[2]);
                        G12:=RayClassGroup(ff,[1,2]);
                        Append(~out,[*b,p,G1,G2,G12*]);
                	end if;
                end if;
    		end for;
    	end if;
    end for;

	return out;


end function;

load "./Magma-Files/LV/janexample.m";

function splitodd(NN) // weight 2k

	k:=1; // maybe 3?
    Z:=Integers();
	Zq<q>:=PowerSeriesRing(Z,NN);

	D:=44;
	f:=7;
	K:=QuadraticField(D);
	O:=Integers(K);
	Fs:=Factorisation(ideal<O|f>);
	ff:=Fs[1][1];
    Norm(ff);
	G,m:=RayClassGroup(ff,[1,2]);
    
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
            	if (J eq G.1) or (J eq G.2)then
        			Gkn:=Gkn  - Norm(I)^(k-1); 
            	elif (J eq G!0) or (J eq (G.1 + G.2)) then
                	Gkn:=Gkn  + Norm(I)^(k-1);
            	end if;
            end if;
        end for;
        Gk:=Gk + 4*Gkn*q^n;
    end for;

	return Gk;

end function;


f:=splitodd(12);
M:=ModularForms(Gamma1(7),2);
B:=Basis(M);
B[1] + 4*B[2] + 12*B[3] + 16*B[4] + 28*B[5];
// Also checked for k = 3 and M_6(Gamma1(7)).






