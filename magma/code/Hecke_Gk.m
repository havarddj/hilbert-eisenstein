/*

April 10, 2019

The code was put together properly again, after some experimentation had rendered some of the routines slightly dysfunctional. 
Attempts to speed up the Hecke approach, by computing the n-th Hecke operator, using the coefficients that were previously already computed,

i.e. T_n(quadratic forms level 1) = T_p(quadratic forms level n/p)        (U_p if p divides n/p)

where we choose p as small as possible. This avoids having to check lots of the same equivalences over and over again.



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
QQ := Rationals();

// ===================================================================
// ===                     Computation of X_n                      ===
// ===================================================================

// Return the gcd of the coefficients of a quadratic form.
function GcdForm(Q)

	return Gcd(Gcd(Q[1],Q[2]),Q[3]);
	
end function; 


// Check whether a form is primitive.
function IsPrimitiveForm(Q)

	return GcdForm(Q) eq 1;

end function;

// Acts on the form [a,b,c] by matrix Mat
function ActionOnForm(Mat,a,b,c)

	// Matrix coefficients;
	q := Mat[1,1];
	r := Mat[1,2];
	s := Mat[2,1];
	t := Mat[2,2];
	// Action of M2(Z);
	A :=   q^2*a +    q*s*b    + s^2*c;
	B := 2*q*r*a + (q*t+q*s)*b + 2*s*t*c;
	C :=   r^2*a +    r*t*b    + t^2*c;
	
	return [A,B,C];

end function;


// Given F, this function returns all SL_2(Z)-translates Q' of Q or Qinv that are separated by 0.
// This means that one of the roots of Q' is less than 0, and the other is greater than 0.
function NearlyReducedForms(F)

	// assert IsPrimitiveForm(F);
	F0 := ReducedForm(F);
	As := [];
	Fs := [];
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
				Append(~As,A);
				Append(~Fs,[A,b+2*i*c,c]);
				//print "Added ", A,b+2*i*c,c;
				assert (b+2*i*c)^2 - 4*A*c eq Discriminant(F);  // REMOVE !!!
			end for;
		end if;
		if a lt 0 then
			for i := 1 to s do
				Append(~As,c);
				Append(~Fs,[c,-b+2*i*c,a-i*b+i^2*c]);
				//C := a-i*b+i^2*c;
				//print "Added ", c,-b+2*i*c,a-i*b+i^2*c;
				assert (-b+2*i*c)^2-4*c*(a-i*b+i^2*c) eq Discriminant(F); // REMOVE !!!
			end for;
		end if;		
	end while;
		
	return As, Fs;

end function;


// Given an integer n>0, return the matrices that define the Hecke operator T_n in level M.
function HeckeMatrices(n : M := 1)

	Sd := [d : d in Divisors(n) | GreatestCommonDivisor(ZZ!(n/d),M) eq 1];
	Mn := [];
	for d in Sd do
		for j := 0 to d-1 do
			Append(~Mn,Matrix([[d,j],[0,ZZ!(n/d)]]));
		end for;
	end for;
	
	return Mn;

end function;



// Given an integer d and discriminant D, find the number of ideals in O_D of norm d, coprime to a fixed ideal of norm M.
// ************
// WARNING !!!
// We assume that for any split prime p dividing M, only one of the two primes above it divides the ideal of norm M.
// This is not relevant if M is 1, which we focus on for now. 
function SplittingFactors(d,D : M := 1)

	mult := 1;
	for fac in Factorisation(d) do
		if KroneckerSymbol(D,fac[1]) eq 1 and IsDivisibleBy(M,fac[1]) eq false then
			mult := mult*(fac[2]+1);
		end if;
		if KroneckerSymbol(D,fac[1]) eq -1 and (IsOdd(fac[2]) or IsDivisibleBy(M,fac[1])) then
			mult := 0;
		end if;
	end for;
	
	return mult;

end function;



// ============================================================
// ===                  Continued fractions                 ===
// ============================================================


// Partial quotients in continued fraction expansion of rational r.
function PartialQuotients(r)

	CC := ContinuedFraction(r);
	d := #CC;
	if d eq 1 then
		return CC,[1],d;
	else
		P := [CC[1],CC[2]*CC[1]+1];
		Q := [1,CC[2]];
		for i in [3..d] do
			Ci := CC[i];
			Append(~P,Ci*P[i-1]+P[i-2]);
			Append(~Q,Ci*Q[i-1]+Q[i-2]);
		end for;	
		
		return P,Q,d;
	end if;
	
end function;

// Returns a list of matrices in SL_2(Z) whose corresponding paths connect r to infinity.
function GammasUnimodular(r)

	if Type(r) eq Type(Infinity) then
  		return [];
  	else
		P,Q,d := PartialQuotients(r);	
		List := [Matrix(2,2,[1,P[1],0,Q[1]])^(-1)];
	
		for i in [1..d-1] do
			Append(~List,Matrix(2,2,[P[i],(-1)^i*P[i+1],Q[i],(-1)^i*Q[i+1]])^(-1));
		end for;
	
		return List;
	end if;
	
end function;

// For r in Q, finds the rational number Gamma*r.
function MatrixAction(Gamma,r)

	if Type(r) eq Type(Infinity) then
		Gammar := Gamma[1,1]/Gamma[2,1];
	else
		if Gamma[2,1]*r+Gamma[2,2] eq 0 then
			Gammar := Infinity;
		else
			Gammar := (Gamma[1,1]*r+Gamma[1,2])/(Gamma[2,1]*r+Gamma[2,2]);
		end if;
	end if;
	
	return Gammar;

end function;


// ===================================================================
// ===                     Computation of X_n                      ===
// ===================================================================


// Naive function for the union of X(n,C) over all ideal classes.
// *********
// WARNING: 
// Already here we ignore completely the difficulty of evaluating the character at the various ideals.
// This is a royal pain, and I imagine it adds a lot of complexity to the naive method. The Hecke method 
// naturally splits up the contributions according to classes.  

function ComputeXn_Naive(D,n : K:=QuadraticField(D), O:=Integers(K), ff := 1*O) 

    Qx<x>:=PolynomialRing(Rationals());
    Bs:={i: i in [0 .. Ceiling(n*SquareRoot(D)) - 1]}; 
    Bs:=Bs join {-B: B in Bs};
    
    delta:=Roots(x^2 - D,K)[1][1];
    
    Xn:=[];
    for B in Bs do
    	alpha := (-B + n*delta)/(2*delta);
    	assert (B^2 - n^2*D) lt 0;
        assert Trace(alpha) eq n;
        if (B mod 2) eq (n*D mod 2) then
    		alphadelta := ideal<O|alpha*delta>;
            divs := Divisors(alphadelta);
            for I in divs do
            	// Nm := ZZ!Norm(I);
            	if GreatestCommonDivisor(I,ff) eq 1*O then
                	Append(~Xn,[* alpha, I *]);
                end if;
            end for;
    	end if;
    end for;
    
    return Xn;
    // return Sort([Norm(a[2]) : a in Xn]);

end function;



// As input, we give:
// D = Discriminant
// n > 0 integer
//
// This function computes the set X(n,C) for a fixed ideal class C, in terms of Hecke operators on quadratic forms. 
// The computation relies on a set of representatives Forms_d of the classes in X(d,C) for d a large divisor of n.

function ComputeXnC_Hecke(D,n,d,Forms_d)
	
	Mn      := HeckeMatrices(ZZ!(n/d) : M := 1); // Hecke operator in level 1!
	Forms_n := [];
	Red_nAs := [];
	for F in Forms_d do
		a := F[1];
		b := F[2];
		c := F[3];
		for Mat in Mn do
			// Matrix coefficients;
			q := Mat[1,1];
			r := Mat[1,2];
			s := Mat[2,1];
			t := Mat[2,2];
			// Action of M2(Z);
			A  :=   q^2*a;
			B  := 2*q*r*a + ZZ!(n/d)*b;
			C  :=   r^2*a + r*t*b + t^2*c;
			MF := QuadraticForms(D*n^2)!<A,B,C>;
			if [IsEquivalent(MF,G) : G in Forms_n] eq [false : G in Forms_n] then
				for ind := 1 to SplittingFactors(GcdForm(MF),D) do
					As, Fs  := NearlyReducedForms(MF);
					Red_nAs := Red_nAs cat As;
				end for;
				Append(~Forms_n,ReducedForm(MF));
			end if;
		end for;
		// Xn := Xn cat RedFn;
		// Append(~Xn,[* RednAs, Forms_n *] *]);
	end for;

//	return Sort(Xn); 
 	return Red_nAs, Forms_n;

end function;


// As input, we give:
// D = Discriminant
// M = (positive) generator of Z cap ff
// n > 0 integer
//
// This function computes the set Xn, in terms of Hecke operators on quadratic forms. 
function ComputeXn_ModularSymbol(F,n)
	
	ZZ := Integers();
	DM2:= Discriminant(F);
	D  := FundamentalDiscriminant(DM2);
	M  := ZZ!((DM2/D)^(1/2));
	
	X1, F1 := NearlyReducedForms(F);
	Xns     := [* [*X1,F1*] *];
	for i := n to n do
		fac := Factorisation(n);
		ndash := ZZ!(n/fac[1][1]);
		
		// Setting up the recursion.
		// STEP 1: Find SL_2(Z) matrices connecting -j/d -> oo to 0 -> oo.
		Gammasn := []; 
		Sd := [d : d in Divisors(n) | GreatestCommonDivisor(ZZ!(n/d),M) eq 1];
		for d in Sd do
			for j := 0 to d-1 do
				Gammasn := Gammasn cat GammasUnimodular(-j/d);
			end for;
		end for;
		
		// STEP 2: Find SL_2(Z) matrices connecting -j/d -> oo to 0 -> oo.	
		//Xndash := Xns[ndash]; // The nearly reduced forms
		Xndash := X1;
		Fndash := F1;
		Xn     := [];
		Fn     := [];
		for gamma in Gammasn do
			for Q in Fndash do
				gammaQ := ActionOnForm(gamma,Q[1],Q[2],Q[3]);
				Append(~Xn,gammaQ[1]);
				Append(~Fn,gammaQ);
			end for;
		end for;
		
		Append(~Xns,[*Xn,Fn*]);
		
	end for;
	
	return Xns;

end function;


// ===================================================================
// ===                     Computation of G_k                      ===
// ===================================================================

	
// Computes the diagonal restriction straight from the definition, i.e. by computing Xn with ideals.
function Gk_Naive(D,k,NN) // weight 2k

	// Not needed yet!
	K  := QuadraticField(D); 
	O  := Integers(K);
	ff := 1*O;

	// assert IsOdd(k) eq true;
	Zq<q> := PowerSeriesRing(ZZ,NN);
    Gk := Zq!0;    
    for n := 1 to NN-1 do
    	Gkn := 0;
        Xn  := ComputeXn_Naive(D,n); // Already coprime to ff!
    	for alphaI in Xn do
        	I   := alphaI[2];
       		Gkn := Gkn + Norm(I)^(k-1);// This is the identity character!! 
        end for;
        Gk := Gk + 4*Gkn*q^n;
    end for;

	return Gk;
	
end function;



// Computes the same series, but with Hecke operators.
// Here, we assume that the character is 
function Gk_Hecke(D,k,NN) // weight 2k

	// Not needed yet!
	K  := QuadraticField(D); 
	O  := Integers(K);
	ff := 1*O;
	M  := 1;

	// Actual computation.
	DM      := D*M;
	Forms_0 := ReducedForms(QuadraticForms(DM));
	Zq<q>   := PowerSeriesRing(ZZ,NN);
    Gk      := Zq!0;
    Forms   := [[* [F] *] : F in Forms_0]; // Save data for all previously computed coefficients.
    for n := 1 to NN-1 do
    	Gkn := 0;
    	As  := [];
    	if n eq 1 then
    		d := 1;
    		for i := 1 to #Forms do
    			As_C, Forms_1 := ComputeXnC_Hecke(D,1,1,Forms[i][1]);
    			Append(~As, As_C);
    		end for;
		else
	    	p := PrimeDivisors(n)[1]; // Smallest prime divisor of n.
	    	d := ZZ!(n/p^(Valuation(n,p)));
	    	// d := ZZ!(n/p);   // This also seems to work, without Hecke relation correction, and is marginally faster. 
    		for i := 1 to #Forms do
    			As_C, Forms_n := ComputeXnC_Hecke(D,n,d,Forms[i][d]);
    			Append(~As, As_C);
    			Append(~Forms[i] , Forms_n);
    		end for;
    	end if;
    	
    	// Now that the As are computed in every class, sum them up with the character to get n-th Fourier coefficient.
    	for AsC in As do
    		C_coeff_n := 0;
    		for a in AsC do
    			C_coeff_n := C_coeff_n + a^(k-1); 
    		end for;
			Gkn := Gkn + C_coeff_n; // This is the identity character!! Need psi(F)*Class_coeff_n here.
        end for;
        // Hecke operator relations:
//         if IsDivisibleBy(d,ZZ!(n/d)) then
//         	p := ZZ!(n/d);
//         	Gkn := Gkn - ZZ!(p^(2*k-1)*Coefficient(Gk, ZZ!(n/p^2))/4);
//         end if;
        
        Gk := Gk + 4*Gkn*q^n;
    end for;

	return Gk;

end function;

// ===================================================================
// ===                            TESTS                            ===
// ===================================================================


function TestHecke(D,k,NN)

	// Compute naively.
	t0 := Cputime();
	Gk_1 := Gk_Naive(D,k,NN);
	print "Total time naive algorithm: ", Cputime() - t0;
	
	// Compute using Hecke.
	t0 := Cputime();
	Gk_2 := Gk_Hecke(D,k,NN);
	print "Total time Hecke algorithm: ", Cputime() - t0;
	
	return "Difference of q-series (should be zero): ", Gk_1 - Gk_2;		

end function;


function TestHecke_Thorough(D,NN)

	// For naive:
    K  := QuadraticField(D);
	O  := Integers(K);
	// For Hecke:
	QQD  := QuadraticForms(D);
	Forms := [[* [F] *] : F in ReducedForms(QQD)]; // Save data for all previously computed coefficients.
	
	Totaltime_naive := 0;
	Totaltime_Hecke := 0;
	
	boo := true;
	
	for n := 1 to NN do
		print " ";
		print " ------------ ";
		
		
		t0 := Cputime();
		A := [Norm(I[2]) : I in ComputeXn_Naive(D,n)];
		time_stepn := Cputime() - t0;
		Totaltime_naive := Totaltime_naive + time_stepn;
		print "Time to compute naively: ", time_stepn;
		
		t0 := Cputime();
    	B := [];
    	if n eq 1 then
    		d := 1;
    		for i := 1 to #Forms do
    			As_C, Forms_1 := ComputeXnC_Hecke(D,1,1,Forms[i][1]);
    			B := B cat As_C;
    		end for;
		else
	    	p := PrimeDivisors(n)[1]; // Smallest prime divisor of n.
	    	d := ZZ!(n/p^(Valuation(n,p)));
    		for i := 1 to #Forms do
    			As_C, Forms_n := ComputeXnC_Hecke(D,n,d,Forms[i][d]);
    			B := B cat As_C;
    			Append(~Forms[i] , Forms_n);
    		end for;
    	end if;
		time_stepn := Cputime() - t0;
		Totaltime_Hecke := Totaltime_Hecke + time_stepn;
		print "Time to compute via Hecke: ", time_stepn;
		
		print "n = ", n, ", same = ", Sort(A) eq Sort(B);
		if Sort(A) eq Sort(B) eq false then
			boo := false;
		end if;
		// print Sort(A);
		// print Sort(B);
	end for;
	
	print "Total time naive algorithm: ", Totaltime_naive;
	print "Total time Hecke algorithm: ", Totaltime_Hecke;

	return "All test passed: ", boo;		

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






