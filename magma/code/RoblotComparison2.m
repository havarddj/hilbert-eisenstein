
K := QuadraticField(401);
O := Integers(K);
D := 401;	
p := 13;
m := 20;
O1 := ideal<O|1>;
ff := O1;
G,mm  := RayClassGroup(ff);
mminv := Inverse(mm);
Rp := pAdicRing(p,m);

	
	
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




// Naive function for X(n,C), where C = trivial class.
function ComputeXn_Naive(n) 

    Qx<x> := PolynomialRing(Rationals());
    delta := Roots(x^2 - D,K)[1][1];	    
	
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
            	if Gcd(I,ff) eq O1 then
            		if mminv(I) eq G.0 then
                		Append(~Xn,[* alpha, I *]);
                	end if;
                end if;
            end for;
    	end if;
    end for;
    
    return Xn;

end function;



// Roblot's example.
function DiagonalRestriction(k,NN) // weight 2k

	Cq<q>:=PowerSeriesRing(Rp,NN);    

    Gk := Cq!0;
    
    for n:=1 to NN-1 do
    	Gkn := 0;
        Xn  := ComputeXn_Naive(n);
    	for alphaI in Xn do
        	I := alphaI[2];
            gcd := GreatestCommonDivisor(I,ff);
            if gcd eq O1 then
            	Gkn := Gkn + Norm(I)^(k-1); 
            end if;
        end for;
        Gk:=Gk + 4*Gkn*q^n;
    end for;

	return Gk;

end function;



// All weights at once. 
function DiagonalRestrictions(Ks,NN)

	Cq<q>:=PowerSeriesRing(Rp,NN);    
	pp:= 19;
	P := Factorisation(ideal<O|pp>)[2][1];
	ff:= P*ideal<O|2>;
	G,m  := RayClassGroup(ff);

	GGs := [Cq!0 : k in Ks];

	for n:=1 to NN-1 do
		Gkn := 0;
		Xn  := ComputeXn_Naive(n);	
		for alphaI in Xn do
			I := alphaI[2];
			gcd := GreatestCommonDivisor(I,ff);
			if gcd eq O1 then
				for i in [1..#Ks] do
					k := Ks[i];
					GGs[i] := GGs[i] + 4*Norm(I)^(k-1)*q^n; 
				end for;
			end if;
		end for;
	end for;

	return GGs;

end function;


// EXPLICIT COMPUTATION OF BOUNDS.

// computes upper bound on degree of P_k(psi) mod p^m
function Find_mdash(p,m)

	mdash:=1;
	if p ne 2 then
		mdash := Floor((p-1)*m/(p-2));
	else
		mdash := m;
	end if;
	
	return mdash; 

end function;




// Main function. 
function pAdicBounds(p,m)

	// q-parameter. 
	q := p;
	if p eq 2 then
		q := 4;
	end if;
	
	// Weights adjusted for p=2. 
	mdash := Find_mdash(p,m);
	Rp := pAdicField(p,mdash);
	Ks := [1 + j*EulerPhi(q): j in [1 .. mdash]];
	S  := PrecisionBound(ModularForms(1,2*Ks[#Ks])) + 1;
	
	print "mdash: ", mdash;
	print "Ks is: ", Ks;
	print "S is : ", S;
	
// 	Fs := Diagonal_Restrictions(F,Ks,S); // S instead of m
// 	ellp  := Precision(Parent(Fs[1]));
// 
// 	k0 := 2*Ks[1];
// //  	Ms := ComputeAllMi(k0,p,ellp,mdash,n);
// 	Ls := [* *]; 
// 	for i := 1 to #Ks do // Sept19 + 1
// 		Append(~Ls,IsClassicalForm(Fs[i],2*Ks[i] : N := N));
// // 		Append(~Ls,QFIsClassicalForm_padic(Fs[i],Ms,i,s));
// 	end for; 

	return "Done";


end function;
