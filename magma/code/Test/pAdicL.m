/*

Master file to compute the p-adic L-series L_p(psi,T). 

*/


load "./DiagonalRestrictions.m";
load "./ClassicalBases.m";


// Works only for unramified characters when h=2. 
function EulerData(p,D,OL : f := 1)
        /* print "p =", p, "d = ", D, "F = ", OL, "N =", f; */
	kron := KroneckerSymbol(D,p);
	
	if Gcd(f,p) gt 1 then
	
		return [[1,p^2]];
		
	else
	
		if kron eq -1 then
			Euler := [[1, p^2]];
		else
			pp := Factorisation(p*OL)[1][1];
			boo, g := IsPrincipal(pp);
			assert boo eq true;
			sgn := 1;
			/* "pp = ", pp; */
			/* "g = ", g; */
			if (IsTotallyPositive(g) eq false) and (IsTotallyPositive(-g) eq false) then
				sgn := -1;
			end if;
			if kron eq 1 then
				Euler := [[sgn,p],[sgn,p]];
			else
				Euler := [[sgn,p]];
			end if;
		end if;
	
		return Euler;
		
	end if;

end function;


// Computes the Euler factor. 
function EulerFac(data,n)

	fac := 1;
	for data_i in data do
		fac := fac*(1-data_i[1]*data_i[2]^n);
	end for;
		
	return fac;

end function;

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


// Interpolation algorithm
function NewtonSolve(kseq,vseq)

    d:=#kseq;

    Qp:=Parent(vseq[1]);
    aseq:=[Qp!0: i in [1 .. d]];
    
    // create matrix: just seems clearer to see things and avoids null sequences
    M:=MatrixRing(Qp,d)!0;
    for i := 1 to d do
        M[i,1] := 1;
    end for;

    for i := 2 to d do
        for j := 2 to i do
            M[i,j] := &*[kseq[i] - kseq[l]: l in [1 .. j-1]];
        end for;
    end for;

    
    aseq[1] := vseq[1];
    for i := 2 to d do
        aseq[i] := M[i,i]^-1*(vseq[i] - &+[aseq[j]*M[i,j]: j in [1 .. i-1]]);
    end for;
    
    Aseq := RSpace(Qp,d)!aseq;
    Vseq := RSpace(Qp,d)!vseq;
    // (Aseq*Transpose(M) - Vseq); // should be zero
    Qpk<k> := PolynomialRing(Qp);
    
    // construct Newton polynomials
    Nseq := [Qpk!0: i in [1 .. d]];
    Nseq[1] := 1;
    
    for i := 2 to d do
        Nseq[i] := &*[(k - kseq[j]): j in [1 .. i-1]];
    end for;
    
    return &+[aseq[i]*Nseq[i]: i in [1 .. d]];
    
end function;




// Main function. 
function pAdicLfunction(F,p,m : trivial := false)

	DF := Discriminant(F);
	L  := QuadraticField(DF);
	OL := Integers(L);
	
	// q-parameter. 
	q := p;
	if p eq 2 then
		q := 4;
	end if;
	
	// Weights adjusted for p=2. 
	mdash := Find_mdash(p,m);
	Rp := pAdicField(p,mdash);
//	n     := 2*m-1;
	Ks := [1 + j*EulerPhi(q): j in [1 .. mdash]];
	N  := Conductor(Parent(F));
	print "Conductor is: ", N;
	S  := PrecisionBound(ModularForms(N,2*Ks[#Ks])) + 1;
	print "Sturm bound is: ", S, Ks;
	Fs := Diagonal_Restriction_Ks(F,Ks,S); // S instead of m

	ellp  := Precision(Parent(Fs[1]));

	k0 := 2*Ks[1];
	print "k0 =", k0;
//  	Ms := ComputeAllMi(k0,p,ellp,mdash,n);
	Ls := [* *]; 
	for i := 1 to #Ks do // Sept19 + 1
		Append(~Ls,IsClassicalForm(Fs[i],2*Ks[i] : N := N));
// 		Append(~Ls,QFIsClassicalForm_padic(Fs[i],Ms,i,s));
	end for;
	
	data := EulerData(p,DF,OL : f := N);
	Lsp  := [EulerFac(data,Ks[i]-1)*Ls[i] : i in [1 .. #Ls]];
	/* print "[EulerFac(data,Ks[i]-1)*Ls[i] : i in [1 .. #Ls]]" */
	/* print "Fs:", Fs; */
	/* print "Base field is: ", L; */
	/* print "Weights are: ", Ks; */
	/* print "L-values are: ", Ls; */
	/* print "Lp-values are: ", Lsp; */
	/* print ""; */
 	/* print "Kronecker symbol is: ", KroneckerSymbol(DF,p); */
	print "Euler data is: ", data;
	Zp   := pAdicField(p,m);
	RR   := PolynomialRing(pAdicRing(p,m));
	RRR<s> := PolynomialRing(Zp);
	if trivial eq true then
		P := RRR!RR!NewtonSolve([1-k : k in Ks],[(-Ks[i])*Lsp[i]: i in [1..#Ks]]);
 		print "Interpolation check: ", [Valuation(Evaluate(P,1-Ks[i])+(Ks[i])*Lsp[i]) : i in [1..#Ks]];
	else
	    P := RRR!RR!NewtonSolve([1-k : k in Ks],[a: a in Lsp]);
	    print "Interpolation check: ", [Valuation(Evaluate(P,1-Ks[i])-Lsp[i]) : i in [1..#Ks]];
	end if;
	print " * * * ";
	print "s polynomial: ", P;
	print "Value at s = 0: ", Evaluate(P,0);
	algdep(Evaluate(P,0),1);
	
	kseqnew := [(1+q)^(1-k) - 1: k in Ks];
	SS<T>    := PolynomialRing(Zp);
	if trivial eq true then
		Q := SS!NewtonSolve(kseqnew,[((1+q)^(1-Ks[i]) - (1+q))*Lsp[i]: i in [1..#Ks]]);
	else
		Q := SS!NewtonSolve(kseqnew,[a: a in Lsp]);
	end if;

	print " * * * ";
	print "T polynomial: ", Q;
	print "Value at T = 0: ", Evaluate(Q,0);

	return P,Q;

end function;



function FunctionalEquationCheck(F,p,m)

	

end function;


// t0:=Cputime();
// F:=D1573F2; // 11^2 * 13
// G:=DirichletGroup(11);
// eps:=G!1;
// 
// Ks:=[1 + j*(p-1): j in [0 .. m-1]];
// 
// s:=PrecisionBound(ModularForms(eps,2*Ks[#Ks])) + 1;
// print "s = ", s;
// 
// Fs:=Diagonal_Restriction_Ks(F,Ks,s); // s instead of m
// 
// t1:=Cputime();
// 
// // load "./Magma-Files/LV-v-4/fran.m"; // for overconvergent case
// 
// load "./Magma-Files/LV-v-5/compute_all_Mi_N.m"; // for ComputeAllMi
// 
// n:=2*m-1;
// mdash:=m;
// ellp:=Precision(Parent(Fs[1]));
// 
// k0:=2;
// // Ms:=ComputeAllMi(k0,p,ellp,mdash,n);
// Ms:=ComputeAllMi(k0,p,mdash,ellp,n,eps);
// 
// t2:=Cputime();
// 
// function QFIsClassicalForm_padic(Gk,Ms,j)
// 
// 	m := Precision(Parent(Gk)); 
// 	// B := Basis(ModularForms(Gamma1(N),k),m);
//     B:=Ms[2*j - 1];
//     // B;
// 	// print "Classical space is of dimension ", #B;
// 	// Matrix of (higher) coefficients
// 	C := Matrix([[Coefficient(B[i],j) : j in [1..m-1]] : i in [1..#B]] cat [[Coefficient(Gk,j) : j in [1..m-1]]]);
// 	K := Basis(Kernel(C));
// 	if #K eq 0 then
// 		print "Not a classical form!";
// 		CT := 0;
// 	else
//     	K;
//         print "#K = ", #K;
//         // remove assertion because could be error forms in basis of high valuation.
// 		// assert #K eq 1; // There must be a unique linear combination.
// 		comb := K[1];
// // 		print "Classical form, linear combination is: ", K[1];
// // 		print "k is: ", k;
// 		CT := -(&+ [Coefficient(B[i],0)*comb[i] : i in [1..#B]])/comb[#B+1];
// 	end if;
// 	
//     // CT;
// 	return CT;
// 
// end function;
// 
// Ls:=[*0*]; // Fs[1] is zero, and the classical space is zero.
// for i:=2 to #Fs do // Sept19 + 1
// 	Append(~Ls,QFIsClassicalForm_padic(Fs[i],Ms,i));
// end for;
// 
// kseq:=[1 + j*(p-1): j in [0 .. #Ls - 1]];
// 
// Lsp:=[(1 - (p^2)^(1 + (i-1)*(p-1) - 1))*Ls[i]: i in [1 .. #Ls]];
// 
// Zp:=pAdicField(p,m);
// P:=NewtonSolve(kseq,[Zp!a: a in Lsp]);
// 
// kseqnew:=[(1+p)^k - 1: k in kseq];
// Q:=NewtonSolve(kseqnew,[Zp!a: a in Lsp]);
