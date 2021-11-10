/*

December 22, 2018

Code using Hecke operator, in an attempt to be competitive with the previous algorithm. 
Uses the computation of X_1, and attempts to compute the sets X_n by an application of the Hecke operators.

It remains to be seen whether this will be competitive with the existing implementation you made! 

*/


// ===================================================================
// ===                  Frequently used matrices                   ===
// ===================================================================

IMat := Matrix(2,2,[[1, 0],[0,1]]);  // Identity
SMat := Matrix(2,2,[[0,-1],[1,0]]);  // Order 4
TMat := Matrix(2,2,[[1, 1],[0,1]]);  // Translation
STMat:= Matrix(2,2,[[0,-1],[1,1]]);  // Order 6 

ZZ := Integers();

// ===================================================================
// ===                     Computation of X_n                      ===
// ===================================================================


// Find the (ordered) roots of a quadratic form, as an element in a degree 2 ext of Q.
function FormRoots(Q,K)
	
	S<x> := PolynomialRing(K);
	f := Q[1]*x^2+Q[2]*x+Q[3];
	Rts := Roots(f);
	return [Minimum(Rts)[1],Maximum(Rts)[1]];

end function;


// Given F, this function returns all SL_2(Z)-translates Q' of Q or Qinv that are separated by 0.
// This means that one of the roots of Q' is less than 0, and the other is greater than 0.
function NearlyReducedForms(F)

	DF := Discriminant(F);
	KK := QuadraticField(DF);

	SepForms := [];
	Red := ReductionOrbit(F);

	for Q in Red do
		RtsQ := FormRoots(Q,KK);
		A := Floor(RtsQ[1])+1;
		B := Floor(RtsQ[2]);
		
		for i := A to B do
			Mi := Matrix([[1,i],[0,1]]);
			QNear := Q*Mi;
			Append(~SepForms,QNear);

			if i ne A and i ne B then
				QNear := Q*Mi*SMat;
				Append(~SepForms,QNear);
			end if;
		end for;
	end for;
		
	return SepForms;

end function;


// Given F, this function returns all SL_2(Z)-translates Q' of Q or Qinv that are separated by 0.
// This means that one of the roots of Q' is less than 0, and the other is greater than 0.
function NearlyReducedForms2(F : M := 1)

	F0 := ReducedForm(F);
	As := [];
	F  := F0;
	started := false;
	while F ne F0 or started eq false do
		started := true;
		a := F[1];
		b := F[2];
		c := F[3];
		F := ReductionStep(F);
		s := Abs(ZZ!((b+F[2])/(2*c))); // The quantity 's' from section 6.4 in Buchmann-Vollmer
		//print "F is: ", F;
		//print "s is: ", s;
		if a gt 0 then
			for i := 1 to s do
				A := a+i*b+i^2*c;
				if GreatestCommonDivisor(A,M) eq 1 then
					Append(~As,A);
					//print "Added ", A,b+2*i*c,c;
					//assert (b+2*i*c)^2 - 4*A*c eq Discriminant(F);
				end if;
			end for;
		end if;
		if a lt 0 then
			if GreatestCommonDivisor(c,M) eq 1 then
				for i := 1 to s do
					Append(~As,c);
					//print "Added ", c,-b+2*i*c,a-i*b+i^2*c;
					//assert (-b+2*i*c)^2-4*c*(a-i*b+i^2*c) eq Discriminant(F);
				end for;
			end if;
		end if;		
	end while;
		
	return As;

end function;


// Given F, this function returns all SL_2(Z)-translates Q' of Q or Qinv that are separated by 0.
// This means that one of the roots of Q' is less than 0, and the other is greater than 0.
function NearlyReducedForms3(F : M := 1)

	F0 := ReducedForm(F);
	As := [];
	F  := F0;
	started := false;
	while F ne F0 or started eq false do
		started := true;
		a := F[1];
		b := F[2];
		c := F[3];
		F := ReductionStep(F);
		s := Abs(ZZ!((b+F[2])/(2*c))); // The quantity 's' from section 6.4 in Buchmann-Vollmer
		//print "F is: ", F;
		//print "s is: ", s;
		if a gt 0 then
			if IsDivisibleBy(c,M) then
				for i := 1 to s do
					A := a+i*b+i^2*c;
					if GreatestCommonDivisor(A,M) eq 1 then
						Append(~As,A);
					//print "Added ", A,b+2*i*c,c;
					//assert (b+2*i*c)^2 - 4*A*c eq Discriminant(F);
					end if;
				end for;
			end if;
		end if;
		if a lt 0 then
			if GreatestCommonDivisor(c,M) eq 1 then
				for i := 1 to s do
					C := a-i*b+i^2*c;
					if IsDivisibleBy(C,M) then
						Append(~As,c);
					//print "Added ", c,-b+2*i*c,a-i*b+i^2*c;
					//assert (-b+2*i*c)^2-4*c*(a-i*b+i^2*c) eq Discriminant(F);
					end if;
				end for;
			end if;
		end if;		
	end while;
		
	return As;

end function;


// Given an integer n>0, return the matrices that define the Hecke operator T_n on quadratic forms of disc DM.
function HeckeMatrices(M,n)

	Sd := [d : d in Divisors(n) | GreatestCommonDivisor(ZZ!(n/d),M) eq 1];
	Mn := [];
	for d in Sd do
		for j := 0 to d-1 do
			Append(~Mn,Matrix([[d,j],[0,ZZ!(n/d)]]));
		end for;
	end for;
	
	return Mn;

end function;


// As input, we give:
// D = Discriminant
// M = (positive) generator of Z cap ff
// n > 0 integer
//
// This function computes the set Xn, in terms of quadratic forms. 
function ComputeXn_Hecke(D,M,n)
	
	DM  := D*M;
	x,y := SquarefreeFactorisation(DM);
	QQ  := QuadraticForms(DM);
	Red := ReducedOrbits(QQ);
	Forms0 := [];
	//for d in Divisors(y) do
	for d in [1] do
		Red := ReducedOrbits(QuadraticForms(ZZ!(DM/d^2)));
		for orb in Red do
			Fy := orb[1];
			Append(~Forms0,QQ!<d*Fy[1],d*Fy[2],d*Fy[3]>);
		end for;
	end for;
	print Forms0;
		
	Mn  := HeckeMatrices(1,n);
	Xn := [];
	for F in Forms0 do
		a := F[1];
		b := F[2];
		c := F[3];
		// print "F is: ", F;
		RedFn := {@ @};
		for Mat in Mn do
			// Matrix coefficients;
			q  := Mat[1,1];
			r  := Mat[1,2];
			s  := Mat[2,1];
			t  := Mat[2,2];
			// Action of M2(Z);
			A :=   q^2*a +       q*s*b +   s^2*c;
			B := 2*q*r*a + (q*t+r*s)*b + 2*s*t*c;
			C :=   r^2*a +       r*t*b +   t^2*c;
			// Action of M2(Z), again (test);
			AA :=   q^2*a;
			BB := 2*q*r*a +   n*b;
			CC :=   r^2*a + r*t*b + t^2*c;
			assert AA eq A and BB eq B and CC eq C; // Test whether they are the same.
			MF  := QuadraticForms(DM*n^2)!<A,B,C>;
			RedFn := RedFn join {@ MFsep : MFsep in NearlyReducedForms(MF) @};
		end for;
		Append(~Xn,RedFn);
	end for;

	// return Xn;
	return Sort([F[1] : F in Xni, Xni in Xn | F[1] gt 0]);

end function;



// As input, we give:
// D = Discriminant
// M = (positive) generator of Z cap ff
// n > 0 integer
//
// This function computes the set Xn, in terms of quadratic forms. 
function ComputeXn_Hecke2(D,M,n)
	
	DM  := D*M;
	x,y := SquarefreeFactorisation(DM);
	QQ  := QuadraticForms(DM);
	Forms0 := [];
	// for d in Divisors(y) do
	for d in [1] do
		Red := ReducedOrbits(QuadraticForms(ZZ!(DM/d^2)));
		for orb in Red do
			Fy := orb[1];
			Append(~Forms0,QQ!<d*Fy[1],d*Fy[2],d*Fy[3]>);
		end for;
	end for;
		
	Mn  := HeckeMatrices(1,n);
	Xn  := [];
	Fn  := [];
	for F in Forms0 do
		a := F[1];
		b := F[2];
		c := F[3];
		RedFn := [];
		for Mat in Mn do
			// Matrix coefficients;
			q := Mat[1,1];
			r := Mat[1,2];
			s := Mat[2,1];
			t := Mat[2,2];
			// Action of M2(Z);
			A :=   q^2*a;
			B := 2*q*r*a +   n*b;
			C :=   r^2*a + r*t*b + t^2*c;
			MF := QuadraticForms(DM*n^2)!<A,B,C>;
			if [IsEquivalent(MF,G) : G in Fn] eq [false : G in Fn] then
				RedFn := RedFn cat NearlyReducedForms3(MF : M := M);
				Append(~Fn,MF);
			end if;
		end for;
		Xn := Xn cat RedFn;
	end for;

 	return Sort(Xn);

end function;


function ComputeXn(D,n : K:=QuadraticField(D), O:=Integers(K), M := 1) // no condition on p


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
    
    //return Xn,delta;
    return Sort([Norm(a[2]) : a in Xn]);

end function;


// Naive function for Xn, to check we get the same answer. 
function ComputeXn2(D,n : K:=QuadraticField(D), O:=Integers(K), M := 1) // no condition on p

    Qx<x>:=PolynomialRing(Rationals());
    Bs:={i: i in [0 .. Ceiling(n*SquareRoot(D)) - 1]}; 
    Bs:=Bs join {-B: B in Bs};
    
    delta:=Roots(x^2 - D,K)[1][1];
    
    Xn:=[];
    for B in Bs do
    	alpha:=(-B + n*delta)/(2*delta);
    	assert (B^2 - n^2*D) lt 0;
        assert Trace(alpha) eq n;
        if (B mod 2) eq (n*D mod 2) then
    		alphadelta := ideal<O|alpha*delta>;
            divs := Divisors(alphadelta);
            for I in divs do
            	Nm := ZZ!Norm(I);
            	if GreatestCommonDivisor(Nm,M) eq 1 then
                	Append(~Xn,Nm);
                end if;
            end for;
    	end if;
    end for;
    
    // return Xn;
    return Sort(Xn);

end function;


function TestHecke(D,M)

    K  := QuadraticField(D);
	O  := Integers(K);
	for n := 1 to 20 do
		print " ";
		print " ------------ ";
		t0 := Cputime();
		A := ComputeXn2(D,n: K:=K,O:=O,M:=M);
		print "Time to compute naively: ", Cputime() - t0;
		t0 := Cputime();
		B  := ComputeXn_Hecke2(D,M,n);
		print "Time to compute via Hecke: ", Cputime() - t0;
		print "n = ", n, ", same = ", A eq B;
		// print A;
		// print B;
	end for;
	
	return "Done.";		

end function;


// ===================================================================
// ===                     Computation of G_k                      ===
// ===================================================================


function G_k(D,k,NN) // weight 2k

    Z:=Integers();
	Zq<q>:=PowerSeriesRing(Z,NN);

    K:=QuadraticField(D);
	O:=Integers(K);
	pp:=3;
	P:=Factorisation(ideal<O|pp>)[1][1];
	ff:=P;
    // ff:=ideal<O|1>;
	G,m:=RayClassGroup(ff,[1,2]);
    minv:=Inverse(m);

    Gk:=Zq!0;
    O1:=ideal<O|1>;
    
    for n:=1 to NN-1 do
    	Gkn:=0;
        Xn:=ComputeXn(D,n);
    	for alphaI in Xn do
        	I:=alphaI[2];
            gcd:=GreatestCommonDivisor(I,ff);
            if gcd eq O1 then
            	J:=minv(I);
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



// ===================================================================
// ===                     PREVIOUS VERSION                        ===
// ===================================================================


/*

December 4, 2018

Fastest code so far for computing the modular form G_k.

Actually perhaps slightly faster with computeXn_newb as a procedure.

Function below is for the Z/2 x Z/2 narrow class group case.

*/

// load "./Magma-Files/LV-2/signs.m";
// 
// function checklist(P,normP,Pnormlist,Pbothlist)
// 
// 	if #Pnormlist eq 0 then
//     	return false,1;
//     end if;
// 
// 	found:=false;
//     normtoolarge:=false;
// 	// normP:=Norm(P);
//     i:=1;
//     while (found eq false) and (normtoolarge eq false) do
//     	if Pnormlist[i] eq normP then
//     		if Pbothlist[i][1] eq P then 
//     			found:=true;
//                 Pinv:=Pbothlist[i][2];
//     		end if;
//     	end if;
//         if Pnormlist[i] gt normP then
//         	normtoolarge:=true;
//             Pinv:=1;
//         end if;
//         i:=i+1;
//         if i gt #Pnormlist then
//         	return false,1;
//         end if;
//     end while;
//     
//     return found,Pinv;
// 
// end function;
// 
// function computeXn_newb(D,n,minv,O1,ff,Pnormlist, Pbothlist: BOUND:=100) 
// 
// 	K:=QuadraticField(D);
//     assert ClassNumber(K) eq 1;
//     O:=Integers(K);
//     Qx<x>:=PolynomialRing(Rationals());
//     
//     Bs:={i: i in [0 .. Ceiling(n*SquareRoot(D)) - 1]}; 
//     Bs:=Bs join {-B: B in Bs};
//     
//     delta:=Roots(x^2 - D,K)[1][1];
//     
//     Xn:=[**];
//     for B in Bs do
//     	alpha:=(-B + n*delta)/(2*delta);
//     	assert (B^2 - n^2*D) lt 0;
//         assert Trace(alpha) eq n;
//         if (B mod 2) eq (n*D mod 2) then
//     		alphadelta:=ideal<O|alpha*delta>;
//             facts:=Factorisation(alphadelta);
//             coprimedivnorms:=[1];
//             coprimedivminv:=[Codomain(minv)!0];
//             //
//             for fs in facts do
//             	P:=fs[1];
//                 if GreatestCommonDivisor(P,ff) eq O1 then
//                 	normP:=Norm(P);
//                     if normP lt BOUND then
//                     	found,Pinv:=checklist(P,normP,Pnormlist, Pbothlist);
//                         if found eq true then // Pinv already computed
//                 			Ppowsnorm:=[normP^i: i in [0 .. fs[2]]];
//                     		Ppowsinv:=[i*Pinv: i in [0 .. fs[2]]];
//                 			coprimedivnorms:=[A*Norm(B): A in coprimedivnorms, B in Ppowsnorm];
//                 			coprimedivminv:=[A + B: A in coprimedivminv, B in Ppowsinv];
//                     	else
//                 			Ppowsnorm:=[normP^i: i in [0 .. fs[2]]];
//                     		Pinv:=minv(P);
//                     		Ppowsinv:=[i*Pinv: i in [0 .. fs[2]]];
//                 			coprimedivnorms:=[A*Norm(B): A in coprimedivnorms, B in Ppowsnorm];
//                 			coprimedivminv:=[A + B: A in coprimedivminv, B in Ppowsinv];
//                         	// now add new one
//                         	Append(~Pnormlist,normP);
//                         	Append(~Pbothlist,[*P,Pinv*]);
//                         	ParallelSort(~Pnormlist,~Pbothlist);
//                     	end if;
//                     else // don't bother and just compute it, possibly not for first time
//                     	Ppowsnorm:=[normP^i: i in [0 .. fs[2]]];
//                     	Pinv:=minv(P);
//                     	Ppowsinv:=[i*Pinv: i in [0 .. fs[2]]];
//                 		coprimedivnorms:=[A*Norm(B): A in coprimedivnorms, B in Ppowsnorm];
//                 		coprimedivminv:=[A + B: A in coprimedivminv, B in Ppowsinv];
//                     end if;
//             	end if;
//             end for;
//             //            
//             for i:=1 to #coprimedivnorms do
//                 Append(~Xn,[*alpha,coprimedivnorms[i],coprimedivminv[i]*]);
//             end for;
//     	end if;
//     end for;
//     
//     return Xn,delta,Pnormlist,Pbothlist;
// 
// end function;
// 
// /*
// 
// We need the narrow ray class group to be Z/2 x Z/2 and each of the less-narrow
// ones to be Z/2. We take then a character on the full group of order 2 which has
// order 2 when restrict to each kernel from the full to narrow groups.
// 
// */
// 
// 
// function compute_Gk_Z2Z2(k,D,f,NN) // weight 2k
// 
// 	assert IsOdd(k);
//     Z:=Integers();
// 	Zq<q>:=PowerSeriesRing(Z,NN);
// 
// 	K:=QuadraticField(D);
// 	O:=Integers(K);
// 	Fs:=Factorisation(ideal<O|f>);
// 	ff:=Fs[1][1]; // note here f may be split or inert.
//     // Norm(ff);
// 	G,m:=RayClassGroup(ff,[1,2]);
//     assert (Order(G) eq 4) and (Exponent(G) eq 2);
//     
//     minv:=Inverse(m);
// 
//     Gk:=Zq!0;
//     O1:=ideal<O|1>;
//     
//     Pnormlist:=[];
//     Pbothlist:=[]; // Ps,Pinvs in pairs
//     
//     k1,k2:=from_narrow_kernels(ff);
//    	assert #k1 eq 2;
//     assert #k2 eq 2;
//     // character we use should have order two on each of these kernels.
//     assert k1[2] ne k2[2]; // perhaps this is always true.
//       
//     
//     for n:=1 to NN-1 do
//     	Gkn:=0;
//         t0:=Cputime();
//         Xn,_,Pnormlist,Pbothlist:=computeXn_newb(D,n,minv,O1,ff,Pnormlist,Pbothlist);
//         // already checked if coprime
//     	for alphaI in Xn do
//         	NormI:=alphaI[2];
//             J:=alphaI[3];
//             if (J eq k1[2]) or (J eq k2[2])then
//             	Gkn:=Gkn  - NormI^(k-1); 
//             else
//                 Gkn:=Gkn  + NormI^(k-1);
//             end if;
//         end for;
//         Gk:=Gk + 4*Gkn*q^n;
//     end for;
// 
// 	return Gk;
// 
// end function;
// 
// // NOW DO THE P-DEPLETING THING
// 
// 
// // Now with a p in it.
// 
// function pdividesI(p,I)
// 
// 	gens:=Generators(I);
//     O:=Parent(gens[1]);
//     
//     for g in gens do	
//     	if ((g/p) in O) eq false then
//         	return false;
//         end if;
//     end for;
// 
// 	return true;
// 
// end function;
// 
// 
// function compute_Gk_Z2Z2_p(D,k,f,N,NN,p) // weight 2k
// 
// 	R:=pAdicRing(p,N);
// 	Zq<q>:=PowerSeriesRing(R,NN);
// 
// 	assert IsOdd(k);
//     Z:=Integers();
// 	Zq<q>:=PowerSeriesRing(Z,NN);
// 
// 	K:=QuadraticField(D);
// 	O:=Integers(K);
// 	Fs:=Factorisation(ideal<O|f>);
// 	ff:=Fs[1][1]; // note here f may be split or inert.
//     // Norm(ff);
// 	G,m:=RayClassGroup(ff,[1,2]);
//     assert (Order(G) eq 4) and (Exponent(G) eq 2);
//     
//     minv:=Inverse(m);
// 
//     Gk:=Zq!0;
//     O1:=ideal<O|1>;
//     
//     Pnormlist:=[];
//     Pbothlist:=[]; // Ps,Pinvs in pairs
//     
//     k1,k2:=from_narrow_kernels(ff);
//    	assert #k1 eq 2;
//     assert #k2 eq 2;
//     // character we use should have order two on each of these kernels.
//     assert k1[2] ne k2[2]; // perhaps this is always true.
//       
//     
//     for n:=1 to NN-1 do
//     	Gkn:=0;
//         t0:=Cputime();
//         Xn,_,Pnormlist,Pbothlist:=computeXn_newb(D,n,minv,O1,ff,Pnormlist,Pbothlist);
//         // already checked if coprime
//     	for alphaI in Xn do
//         	NormI:=alphaI[2];
//             if NormI mod p ne 0 then // I'M NOT COMPLETELY SURE THIS IS CORRECT CONDITION!!!
//             	J:=alphaI[3];
//             	if (J eq k1[2]) or (J eq k2[2])then
//                 	Gkn:=Gkn  - NormI^(k-1); 
//             	else
//                 	Gkn:=Gkn  + NormI^(k-1);
//                 end if;
//             end if;
//         end for;
//         Gk:=Gk + 4*Gkn*q^n;
//     end for;
// 
// 	return Gk;
// 
// end function;






