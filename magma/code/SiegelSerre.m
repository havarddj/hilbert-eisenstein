// Code for the paper on p-adic L-functions with Alan.

// -----------------------------------------------------
// -----------------------------------------------------

//        TABLE OF CONTENTS:
// 
//        1) Cubic example
//        2) Real quadratic fields
//        3) Tests and checks


load "Discriminants.m";

// Common structures:
ZZ := Integers();
QQ := Rationals(); 

IMat := Matrix(2,2,[[1, 0],[0,1]]);  // Identity
SMat := Matrix(2,2,[[0,-1],[1,0]]);  // Order 4
TMat := Matrix(2,2,[[1, 1],[0,1]]);  // Translation
STMat:= Matrix(2,2,[[0,-1],[1,1]]);  // Order 6 

// -------------------------------------------------------------------------------
//                              1) Cubic example
// -------------------------------------------------------------------------------

// Test case for cubic field generated by zeta + zeta^{-1} for zeta = 18th root of unity. 

// INPUT:
// n = Precision O(q^n)
// k = Weight, i.e. parallel weight (k,k,k) and diagonal restriction weight 3k.
//
// OUTPUT: 
// The coefficients of q, q^2, ..., q^(n-1) of the diagonal restriction with respect to
//     psi = Unique nontrivial totally even character of F with conductor 5. 

function CubicExample(n,Ks)

	if Type(Ks) eq Type(0) then
		Ks := [Ks];
	end if;
	
	for k in Ks do
		assert IsEven(k); // The character psi chosen below (hard coded) is totally even! 
	end for;
	
	R<x> := PolynomialRing(Rationals());
	F<a> := NumberField(x^3-3*x-1);
	OF   := Integers(F);
	dd   := Different(OF);
	d    := Degree(F);
	gen  := 1/(3*a^2 - 3);
	print "Traces (Should be 0 0 1) :", Trace(gen), Trace(a*gen), Trace(a^2*gen);
	
	// Character.
	psi := HeckeCharacterGroup(5*OF).1;

	vec1 := Vector(RealEmbeddings(gen));
	vec2 := Vector(RealEmbeddings(a*gen));
	vec3 := Vector(RealEmbeddings(a^2*gen));
	A := Matrix([vec1,vec2]);
	
	R<q>  := PowerSeriesRing(ZZ,n);
	Diags := [R!0 : k in Ks];
	
	for j := 1 to n-1 do
		print "Doing coefficient ", j;
		w := -j*vec3;
		X := []; // x-coordinates of triangle points.
		Y := []; // y-coordinates of triangle points.
		for i := 1 to d do
			Ai := RemoveColumn(A,i);
			wi := RemoveColumn(w,i);
			sol := Solution(Ai,wi)[1];
			Append(~X,sol[1]);
			Append(~Y,sol[2]);
		end for;
		Nus := [];
		for x := Floor(Minimum(X)) to Ceiling(Maximum(X)) do
			for y := Floor(Minimum(Y)) to Ceiling(Maximum(Y)) do
				nu := x*gen + y*a*gen + j*a^2*gen;
				if IsTotallyPositive(nu) then
					Append(~Nus, nu);
				end if;
			end for;
		end for;
		
		coeff_js := [0 : k in Ks];
		for nu in Nus do
        	assert Trace(nu) eq j;
            divs := Divisors(nu*dd);
            for I in divs do
            	NmI := Norm(I);
            	if IsDivisibleBy(NmI,5) eq false then
            		psiI := psi(I);
            		for i in [1..#Ks] do
            			k := Ks[i];
            			coeff_js[i] := coeff_js[i] + psiI*NmI^(k-1);
            		end for;
				end if;	
			end for;
		end for;
		
		for i in [1..#Ks] do
			Diags[i] := Diags[i] + 2^d*coeff_js[i]*q^j;
		end for;
		
	end for;
	
	return Diags;

end function;




// -------------------------------------------------------------
//                  2) Real quadratic fields
// -------------------------------------------------------------


// Return the gcd of the coefficients of a quadratic form.
function GcdForm(Q)

	return Gcd(Gcd(Q[1],Q[2]),Q[3]);
	
end function; 


// Check whether a form is primitive.
function IsPrimitiveForm(Q)

	return GcdForm(Q) eq 1;

end function;


// Check whether a given 2x2 matrix 'g' is integral
function IsIntegralMatrix(g)

	d := Degree(Parent(g));
	return [IsIntegral(g[i,j]) : i,j in [1..d]] eq [true : i,j in [1..d]];

end function;


// Check whether pair1 = (w1, delta1) and pair2 = (w2, delta2) are equivalent in RM(n,C)
function EquivalenceRM(pair1,pair2)

	boo := false;
	w1  := pair1[1];
	w2  := pair2[1];
	delta1 := pair1[2];
	delta2 := pair2[2];
	equiv, gamma := IsEquivalent(w1,w2);
	if equiv then
		delta2 := ChangeRing(delta2,Rationals());
		boo := IsIntegralMatrix(delta1*gamma*delta2^(-1));
	end if;
	
	return boo;

end function;


// Given F, this function returns all SL_2(Z)-translates Q' of Q or Qinv that are separated by 0.
// This means that one of the roots of Q' is less than 0, and the other is greater than 0.
function NearlyReducedForms(F)
 
 	F0 := ReducedForm(F);
 	As_p := [];  // Positive a's
   	As_m := [];  // Negative a's
//  	Fs_p := [];  // All forms with positive a
//  	Fs_m := [];  // All forms with negative a

	F  := F0;
	started := false;
 	while F ne F0 or started eq false do
 		started := true;
 		a := F[1];
 		b := F[2];
 		c := F[3];
 		F := ReductionStep(F);
 		s := Abs(ZZ!((b+F[2])/(2*c))); // The quantity 's' from section 6.4 in Buchmann-Vollmer
 		if a gt 0 then
 			for i := 1 to s do
 				A := a+i*b+i^2*c;
 				Append(~As_p,A);
 				Append(~As_m,c);
// 				Append(~Fs,[A,b+2*i*c,c]);
// 				print "Added ", A,b+2*i*c,c;
//  				assert (b+2*i*c)^2 - 4*A*c eq Discriminant(F);  // REMOVE !!!
 			end for;
 		end if;
 		if a lt 0 then
 			for i := 1 to s do
 				A := a-i*b+i^2*c;
 				Append(~As_p,c);
 				Append(~As_m,A);
// 				Append(~Fs,[c,-b+2*i*c,a-i*b+i^2*c]);
// 				print "Added ", c,-b+2*i*c,a-i*b+i^2*c;
//  				assert (-b+2*i*c)^2-4*c*(a-i*b+i^2*c) eq Discriminant(F); // REMOVE !!!
 			end for;
 		end if;
 	end while;
 		
 	return As_p, As_m;
 	
end function;


// Given an integer n>0, return the matrices that define the Hecke operator T_n in level N.
function HeckeMatrices(n : N := 1)

	Sd := [d : d in Divisors(n) | GreatestCommonDivisor(ZZ!(n/d),N) eq 1];
	Mn := [];
	for d in Sd do
		for j := 0 to d-1 do
			Append(~Mn,Matrix([[d,j],[0,ZZ!(n/d)]]));
		end for;
	end for;
	
	return Mn;

end function;


// This function computes the set RM+-(n,C) for a fixed ideal class C, using Hecke operators on quadratic forms. 
// The computation relies on a set of representatives Forms_d of the classes in RM+-(n,C) for d a large divisor of n.
function RM_Points(n,d,Forms_d : D := ZZ!(Discriminant(Forms_d[1][1])/d^2))

 	Mn      := HeckeMatrices(ZZ!(n/d) : N := 1); // Hecke operator in level 1!
 	Forms_n := [];
 	A_pos_n := [];
 	A_neg_n := [];
 	for pair in Forms_d do
 		a := pair[1][1];
 		b := pair[1][2];
 		c := pair[1][3];
  		delta_d := pair[2];
 		for Mat in Mn do
 			// Matrix coefficients;
 			q := Mat[1,1];
 			r := Mat[1,2];
 			s := Mat[2,1];
 			t := Mat[2,2];
 			// Action of M2(Z), uses the fact that s=0 for these Hecke representatives!! 
 			A  :=   q^2*a;
 			B  := 2*q*r*a + ZZ!(n/d)*b;
 			C  :=   r^2*a + r*t*b + t^2*c;
 			MF := QuadraticForms(D*n^2)!<A,B,C>;
 			delta_n := delta_d*Mat;
 			RM_n := [* MF,delta_n *];
 			if [EquivalenceRM(RM_n,G) : G in Forms_n] eq [false : G in Forms_n] then
 				As_p, As_m  := NearlyReducedForms(MF);
				A_pos_n := A_pos_n cat As_p;
				A_neg_n := A_neg_n cat As_m;
 				Append(~Forms_n,RM_n);
 			end if;
 		end for;
 	end for;
 
  	return A_pos_n, A_neg_n, Forms_n;
 
end function;


// Computes the nearly reduced forms for the first m coefficients of the diagonal restriction.
function Diagonal_Restriction_data(F,m) // weight 2k

	// Initial values:
	A_p,A_n,Forms_0 := RM_Points(1,1,[[*F,IMat*]]);
	Forms := [Forms_0];
	As    := [[A_p,A_n]];
	// Hecke operators
	for n := 2 to m-1 do
		Diag_F_n := 0;
		if n eq 1 then
			p := 1;
			d := 1;
		else 
 	    	p := PrimeDivisors(n)[1]; // Smallest prime divisor of n.
 	    	d := ZZ!(n/p^(Valuation(n,p)));
 	    end if;
		A_pos_n, A_neg_n, Forms_n := RM_Points(n,d,Forms[d]);
		Append(~As, [A_pos_n, A_neg_n]);
   		Append(~Forms, Forms_n);
   	end for;
   	
   	return As, Forms;
   	
end function;


// Computes the first m coefficients (EXCEPT the constant term!) of the diagonal restriction of 
//            E_{k,k}^{p}(1,psi_F) 
// where F is a quadratic form, whose discriminant we assume for now to be fundamental.  
function diagonal_restriction(F,k,m,As,Forms : Stabilisation := 1)

	R<q>   := PowerSeriesRing(ZZ,m);
	Diag_F := R!0;

	f := Conductor(Parent(F));
	p := Stabilisation;
	assert IsPrime(p) or p eq 1;
	
	for n in [1..m-1] do
		A_p := As[n][1];
		A_m := As[n][2];
		coeff_n := 0;
		for i := 1 to #A_p do
			if Gcd(A_p[i],p*f) eq 1 then
				coeff_n := coeff_n + A_p[i]^(k-1);
			end if;
			if Gcd(A_m[i],p*f) eq 1 then
				coeff_n := coeff_n + (-1)^k*(-A_m[i])^(k-1);
			end if;
		end for;
		
		Diag_F := Diag_F + 4*coeff_n*q^n;
		
	end for;
	
	return Diag_F;

end function;

// Computes the first m coefficients (EXCEPT the constant term!) of the diagonal restriction of 
//            d/dk E_{k,k}^{p}(1,psi_F) | k=1
// where F is a quadratic form, whose discriminant we assume for now to be fundamental.  
function diagonal_restriction_derivative(F,p,m,As,Forms)

	Rp     := pAdicRing(p,m);
	R<q>   := PowerSeriesRing(Rp,m);
	Diag_F := R!0;
	
	f := Conductor(Parent(F));
	Logs   := AssociativeArray(); // This is to avoid having to compute logarithms of the same number over and over again 
	                              // appears to be taking a lot of time!
	for n in [1..m-1] do
		A_p := As[n][1];
		A_m := As[n][2];
		coeff_n := 0;
		for i := 1 to #A_p do
			// Positive a's
			a := A_p[i];
			if Gcd(a,p*f) eq 1 then
				boo, log := IsDefined(Logs,a);
				if boo eq false then
					log := Log(Rp!a);
					Logs[a] := log;
				end if;
				coeff_n := coeff_n + log;
			end if;
			// Negative a's
			a := A_m[i];
			if Gcd(a,p*f) eq 1 then
				boo, log := IsDefined(Logs,a);
				if boo eq false then
					log := Log(Rp!a);
					Logs[a] := log;
				end if;
				coeff_n := coeff_n - log;
			end if;
 		end for;
		Diag_F := Diag_F + 4*coeff_n*q^n;
	end for;
	
	return Diag_F;

end function;



// Computes the first m coefficients (EXCEPT the constant term!) of the diagonal restriction of 
//            E_{k,k}^{p}(1,psi_F) 
// where F is a quadratic form, whose discriminant we assume for now to be fundamental.  
function Diagonal_Restriction(F,k,m : Stabilisation := 1)

	As, Forms := Diagonal_Restriction_data(F,m);
	
	return diagonal_restriction(F,k,m,As,Forms : Stabilisation := Stabilisation);
	
end function;



// Computes the first m coefficients (EXCEPT the constant term!) of the diagonal restriction of 
//            d/dk E_{k,k}^{p}(1,psi_F) | k=1
// where F is a quadratic form, whose discriminant we assume for now to be fundamental.  
function Diagonal_Restriction_Derivative(F,p,m)

	As, Forms := Diagonal_Restriction_data(F,m);
	
	return diagonal_restriction_derivative(F,p,m,As,Forms);
	
end function;



// Computes the first m coefficients (EXCEPT the constant term!) of the first order deformation (in k) of the diagonal restriction of 
//            E_{k,k}^{p}(1,psi_F)
// where F is a quadratic form, whose discriminant we assume for now to be fundamental.  
function WeakHarmonicMaassForm(F,p,m)
	
	sgn := KroneckerSymbol(Discriminant(F),p);
	if sgn eq 1 then
		print "Coherent case";
		print "Kronecker symbol ", sgn;
	else
		print "Incoherent case";
		print "Kronecker symbol ", sgn;
	end if;
	
	t0 := Cputime();
	As, Forms := Diagonal_Restriction_data(F,m);
	print "Time for data: ", Cputime() - t0;
	
	t0 := Cputime();
	G  := diagonal_restriction(F,1,m,As,Forms : Stabilisation := p);  // Weight k=1
	print "Time for diagonal restriction: ", Cputime() - t0;
	
	t0 := Cputime();
	Geps := diagonal_restriction_derivative(F,p,m,As,Forms);
	print "Time for derivative: ", Cputime() - t0;

	return G, Geps;
	
end function;



// -------------------------------------------------------------
//                    3) Tests and checks
// -------------------------------------------------------------



// Algebraic recognition for element in Z/p^mZ.
function algdep(a,deg)

// 	N  := Modulus(Parent(a));
	p  := Prime(Parent(a));
	m  := Precision(Parent(a))-5;
	N  := p^m;
	ZZ := Integers();
	RR := RealField(500);
	M  := ZeroMatrix(RR,deg+2,deg+2);

	for i:=0 to deg do
		M[i+1,1] := ZZ!((a)^i);
	end for;
	M[deg+2][1] := N;

	for j:=1 to deg+1 do
		M[j,j+1] := 1;
	end for;

	Y,T := LLL(M);

	PolZ<t> := PolynomialRing(ZZ);
	P       := Floor(Y[2][1] - Y[2][2]);
	for j := 3 to deg+2 do
		P := P - Floor(Y[2][j])*t^(j-2);
	end for;

	return P;

end function;




// Check that a given q-expansion is (modulo the constant term) a member of the space 
// M_k(chi) of classical forms of weight k.
function IsClassicalForm(Gk,k : N:=1)

	m := Precision(Parent(Gk));
	B := Basis(ModularForms(Gamma1(N),k),m);
	print "Classical space is of dimension ", #B;
	// Matrix of (higher) coefficients
	C := Matrix([[Coefficient(B[i],j) : j in [1..m-1]] : i in [1..#B]] cat [[Coefficient(Gk,j) : j in [1..m-1]]]);
	K := Basis(Kernel(C));
	if #K eq 0 then
		print "Not a classical form!";
		CT := 0;
	else
		assert #K eq 1; // There must be a unique linear combination.
		comb := K[1];
// 		print "Classical form, linear combination is: ", K[1];
// 		print "k is: ", k;
		CT := -(&+ [Coefficient(B[i],0)*comb[i] : i in [1..#B]])/comb[#B+1];
	end if;
	
	return CT;

end function;




// Naive function for X(n,C) and X(n,C*), where C = [I].
function ComputeXn_Naive(n,F) 

	D := Discriminant(F);
	K := QuadraticField(D);
	O := Integers(K);
    Qx<x> := PolynomialRing(Rationals());
    delta := Roots(x^2 - D,K)[1][1];	    
	J := ideal<O | F[1], (-F[2] + delta)/2>;
	
    Bs := {i: i in [0 .. Ceiling(n*SquareRoot(D)) - 1]}; 
    Bs := Bs join {-B: B in Bs};
    
    Xn := [];
    for B in Bs do
    	alpha := (-B + n*delta)/(2*delta);
    	assert (B^2 - n^2*D) lt 0;
        assert Trace(alpha) eq n;
        if (B mod 2) eq (n*D mod 2) then
    		alphadelta := ideal<O|alpha*delta>;
            divs := Divisors(alphadelta);
            for I in divs do
            	if IsPrincipal(I*J) then
                	Append(~Xn,[* alpha, I *]);
                	print " a b c : ", Norm(I), B, Norm(alphadelta*I^(-1));
                end if;
            end for;
    	end if;
    end for;
    
    return Sort([Norm(a[2]) : a in Xn]);

end function;


