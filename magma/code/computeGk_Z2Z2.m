/*

December 4, 2018

Fastest code so far for computing the modular form G_k.

Actually perhaps slightly faster with computeXn_newb as a procedure.

Function below is for the Z/2 x Z/2 narrow class group case.

*/

load "./Magma-Files/LV-2/signs.m";

function checklist(P,normP,Pnormlist,Pbothlist)

	if #Pnormlist eq 0 then
    	return false,1;
    end if;

	found:=false;
    normtoolarge:=false;
	// normP:=Norm(P);
    i:=1;
    while (found eq false) and (normtoolarge eq false) do
    	if Pnormlist[i] eq normP then
    		if Pbothlist[i][1] eq P then 
    			found:=true;
                Pinv:=Pbothlist[i][2];
    		end if;
    	end if;
        if Pnormlist[i] gt normP then
        	normtoolarge:=true;
            Pinv:=1;
        end if;
        i:=i+1;
        if i gt #Pnormlist then
        	return false,1;
        end if;
    end while;
    
    return found,Pinv;

end function;

function computeXn_newb(D,n,minv,O1,ff,Pnormlist, Pbothlist: BOUND:=100) 

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
            facts:=Factorisation(alphadelta);
            coprimedivnorms:=[1];
            coprimedivminv:=[Codomain(minv)!0];
            //
            for fs in facts do
            	P:=fs[1];
                if GreatestCommonDivisor(P,ff) eq O1 then
                	normP:=Norm(P);
                    if normP lt BOUND then
                    	found,Pinv:=checklist(P,normP,Pnormlist, Pbothlist);
                        if found eq true then // Pinv already computed
                			Ppowsnorm:=[normP^i: i in [0 .. fs[2]]];
                    		Ppowsinv:=[i*Pinv: i in [0 .. fs[2]]];
                			coprimedivnorms:=[A*Norm(B): A in coprimedivnorms, B in Ppowsnorm];
                			coprimedivminv:=[A + B: A in coprimedivminv, B in Ppowsinv];
                    	else
                			Ppowsnorm:=[normP^i: i in [0 .. fs[2]]];
                    		Pinv:=minv(P);
                    		Ppowsinv:=[i*Pinv: i in [0 .. fs[2]]];
                			coprimedivnorms:=[A*Norm(B): A in coprimedivnorms, B in Ppowsnorm];
                			coprimedivminv:=[A + B: A in coprimedivminv, B in Ppowsinv];
                        	// now add new one
                        	Append(~Pnormlist,normP);
                        	Append(~Pbothlist,[*P,Pinv*]);
                        	ParallelSort(~Pnormlist,~Pbothlist);
                    	end if;
                    else // don't bother and just compute it, possibly not for first time
                    	Ppowsnorm:=[normP^i: i in [0 .. fs[2]]];
                    	Pinv:=minv(P);
                    	Ppowsinv:=[i*Pinv: i in [0 .. fs[2]]];
                		coprimedivnorms:=[A*Norm(B): A in coprimedivnorms, B in Ppowsnorm];
                		coprimedivminv:=[A + B: A in coprimedivminv, B in Ppowsinv];
                    end if;
            	end if;
            end for;
            //            
            for i:=1 to #coprimedivnorms do
                Append(~Xn,[*alpha,coprimedivnorms[i],coprimedivminv[i]*]);
            end for;
    	end if;
    end for;
    
    return Xn,delta,Pnormlist,Pbothlist;

end function;

/*

We need the narrow ray class group to be Z/2 x Z/2 and each of the less-narrow
ones to be Z/2. We take then a character on the full group of order 2 which has
order 2 when restrict to each kernel from the full to narrow groups.

*/


function compute_Gk_Z2Z2(k,D,f,NN) // weight 2k

	assert IsOdd(k);
    Z:=Integers();
	Zq<q>:=PowerSeriesRing(Z,NN);

	K:=QuadraticField(D);
	O:=Integers(K);
	Fs:=Factorisation(ideal<O|f>);
	ff:=Fs[1][1]; // note here f may be split or inert.
    // Norm(ff);
	G,m:=RayClassGroup(ff,[1,2]);
    assert (Order(G) eq 4) and (Exponent(G) eq 2);
    
    minv:=Inverse(m);

    Gk:=Zq!0;
    O1:=ideal<O|1>;
    
    Pnormlist:=[];
    Pbothlist:=[]; // Ps,Pinvs in pairs
    
    k1,k2:=from_narrow_kernels(ff);
   	assert #k1 eq 2;
    assert #k2 eq 2;
    // character we use should have order two on each of these kernels.
    assert k1[2] ne k2[2]; // perhaps this is always true.
      
    
    for n:=1 to NN-1 do
    	Gkn:=0;
        t0:=Cputime();
        Xn,_,Pnormlist,Pbothlist:=computeXn_newb(D,n,minv,O1,ff,Pnormlist,Pbothlist);
        // already checked if coprime
    	for alphaI in Xn do
        	NormI:=alphaI[2];
            J:=alphaI[3];
            if (J eq k1[2]) or (J eq k2[2])then
            	Gkn:=Gkn  - NormI^(k-1); 
            else
                Gkn:=Gkn  + NormI^(k-1);
            end if;
        end for;
        Gk:=Gk + 4*Gkn*q^n;
    end for;

	return Gk;

end function;

// NOW DO THE P-DEPLETING THING


// Now with a p in it.

function pdividesI(p,I)

	gens:=Generators(I);
    O:=Parent(gens[1]);
    
    for g in gens do	
    	if ((g/p) in O) eq false then
        	return false;
        end if;
    end for;

	return true;

end function;


function compute_Gk_Z2Z2_p(D,k,f,N,NN,p) // weight 2k

	R:=pAdicRing(p,N);
	Zq<q>:=PowerSeriesRing(R,NN);

	assert IsOdd(k);
    Z:=Integers();
	Zq<q>:=PowerSeriesRing(Z,NN);

	K:=QuadraticField(D);
	O:=Integers(K);
	Fs:=Factorisation(ideal<O|f>);
	ff:=Fs[1][1]; // note here f may be split or inert.
    // Norm(ff);
	G,m:=RayClassGroup(ff,[1,2]);
    assert (Order(G) eq 4) and (Exponent(G) eq 2);
    
    minv:=Inverse(m);

    Gk:=Zq!0;
    O1:=ideal<O|1>;
    
    Pnormlist:=[];
    Pbothlist:=[]; // Ps,Pinvs in pairs
    
    k1,k2:=from_narrow_kernels(ff);
   	assert #k1 eq 2;
    assert #k2 eq 2;
    // character we use should have order two on each of these kernels.
    assert k1[2] ne k2[2]; // perhaps this is always true.
      
    
    for n:=1 to NN-1 do
    	Gkn:=0;
        t0:=Cputime();
        Xn,_,Pnormlist,Pbothlist:=computeXn_newb(D,n,minv,O1,ff,Pnormlist,Pbothlist);
        // already checked if coprime
    	for alphaI in Xn do
        	NormI:=alphaI[2];
            if NormI mod p ne 0 then // I'M NOT COMPLETELY SURE THIS IS CORRECT CONDITION!!!
            	J:=alphaI[3];
            	if (J eq k1[2]) or (J eq k2[2])then
                	Gkn:=Gkn  - NormI^(k-1); 
            	else
                	Gkn:=Gkn  + NormI^(k-1);
                end if;
            end if;
        end for;
        Gk:=Gk + 4*Gkn*q^n;
    end for;

	return Gk;

end function;




