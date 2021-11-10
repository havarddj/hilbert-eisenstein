/*

March 22, 2019

Building level one classical spaces using Katz expansion code.

March 25, 2019

Put all code in one file. This should compute p-adic L-series around kappa0 = 1 (and 1 + (p-1)/2).
Here the conductor is trivial and Q(sqrt(D)) should be a field of class number 1 and narrow class number 2. Also p should be inert in K.


*/

//  *** LEVEL 1 *** // Dec 18, 2015: from webpage

// Returns Eisenstein series E_k over given ring S modulo q^NN normalised
// to have constant term a_0 = 1.

function ESeries(k,NN,S)

	R<q> := PowerSeriesRing(S,NN);
	a1 := -S!(2*k/Bernoulli(k));
	Ek := 1 + a1*q;
	for n := 2 to NN-1 do
		coeffn := 0;
		divs := Divisors(n);
		for d in divs do
			coeffn := coeffn + d^(k-1);
		end for;
		Ek := Ek + a1*coeffn*q^n;
	end for;

	return Ek;

end function;

// Returns dimension of space of modular forms of weight k and level 1.

function DimensionMk(k)

	if ((k mod 2) eq 1) or (k lt 0) then
		return 0;
	end if;
	kmod12 := k mod 12;
	if (kmod12 eq 2) then
		return Floor(k/12);
	else
		return Floor(k/12) + 1;
	end if;

end function;

// Returns basis for W_i (those elements in the Miller basis of weight k which 
// are not in E_(p-1) x {weight (k - (p-1)) forms} and a power of delta
// which can be reused on the next call of function.

function ComputeWi(k,p,delta,deltaj,E4,E6)
	
	// Define a and b
	a := k mod 3;
	b := (k div 2) mod 2;
		
	// Compute dimensions required for Miller basis
	d := DimensionMk(k) - 1;
	e := DimensionMk(k - (p-1)) - 1;
	
	// Construct basis for Wi
	Wi := [];
	for j := e+1 to d do
		// compute aj = delta^j*E6^(2*(d-j) + b)*E4^a
		aj := deltaj*E6^(2*(d-j) + b)*E4^a;
		deltaj := deltaj*delta;
		Wi := Append(Wi,aj);
	end for;

	return Wi,deltaj;

end function;

// CLASSICAL SPACES

function ComputeAllWi(k0,p,ellp,mdash,n) 
	
	S := IntegerRing(p^mdash); // See Footnote (1) below
	Ep1 := ESeries(p-1,ellp,S);

	E4 := ESeries(4,ellp,S);
	E6 := ESeries(6,ellp,S);
	q := Parent(E4).1;
	delta := Delta(q);
	
	deltaj := Parent(delta)!1;
	Wis := [**];	
	for i := 0 to n do
		Wi,deltaj := ComputeWi(k0 + i*(p-1),p,delta,deltaj,E4,E6);
		Append(~Wis,Wi);
	end for;
	
	return Wis,Ep1;

end function;

// Computes bases for all spaces M_(k0 _ i*(p-1)) for i = 0,..,n
// to precision q^ellp, p^mdash

function ComputeAllMi(k0,p,ellp,mdash,n)

	Wis,Ep1:=ComputeAllWi(k0,p,ellp,mdash,n);
    
    Mis:=[**];
    for i:=0 to n do
    	Mi:=[**];
    	for j:=0 to i do
        	Ep1j:=Ep1^(i-j);
    		for e in Wis[j+1] do
            	Append(~Mi,Ep1j*e);
            end for;
    	end for;
        Append(~Mis,Mi);
    end for;

	return Mis;

end function;

// Checking code

function FormInBasis(ord,f)

	dimord:=#ord;
	g:=f;
	
	S:=Parent(Coefficient(f,0)); // IntegerRing()
	ZZ:=Integers();
	fvec:=RSpace(S,dimord)!0;
	
	for j:=1 to dimord do
		ejleadingpos:=Valuation(ord[j]);
		ejleadingcoeff:=Coefficient(ord[j],ejleadingpos);
		fvec[j]:=S!((ZZ!Coefficient(g,ejleadingpos))/(ZZ!ejleadingcoeff));
		g:=g - fvec[j]*ord[j];
	end for;
	
    g;
    /*
	if g ne 0 then
		print "Error in OrdFormInBasis: form NOT l.c. of basis!";
		print "... leading coeff of g is: ", Coefficient(g,Valuation(g)); // Apr 10, 2017
	end if;
	*/
	
	return fvec;

end function;

function ComputeAllMiNaive(k0,p,ellp,mdash,n)

	Mis:=[**];
    for i:=0 to n do
    	Mi:=Basis(ModularForms(1,k0 + i*(p-1)),ellp);
        Append(~Mis,Mi);
    end for;

	return Mis;

end function;

function compare(k0,p,ellp,mdash,n)

	time Mis:=ComputeAllMi(k0,p,ellp,mdash,n);
    time MisN:=ComputeAllMiNaive(k0,p,ellp,mdash,n);

	for i:=1 to #Mis do
    	for f in Mis[i] do
    		_:=FormInBasis(MisN[i],f);
    	end for;
    end for;
    
    for i:=1 to #MisN do
    	for f in MisN[i] do
    		_:=FormInBasis(Mis[i],f);
    	end for;
    end for;

	return 0;

end function;

// *************************************************************** //

// ***** CODE FOR COMPUTING G_K ****** //

// This computes unstabilised version, which is what we need here.
// from LV-v-3/Gk_faster

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


procedure computeXn_faster(D,n,minv,O1,~Pnormlist,~Pbothlist, ~Xn: BOUND:=100) 

	K:=QuadraticField(D);
    assert ClassNumber(K) eq 1;
    O:=Integers(K);
    Qx<x>:=PolynomialRing(Rationals());
    
    Bs:={i: i in [0 .. Ceiling(n*SquareRoot(D)) - 1]}; 
    Bs:=Bs join {-B: B in Bs};
    
    delta:=Roots(x^2 - D,K)[1][1];
    
    // Xn:=[**];
    for B in Bs do
    	alpha:=(-B + n*delta)/(2*delta);
    	assert (B^2 - n^2*D) lt 0;
        assert Trace(alpha) eq n;
        if (B mod 2) eq (n*D mod 2) then
    		alphadelta:=ideal<O|alpha*delta>;
            //divs:=Divisors(alphadelta);
            facts:=Factorisation(alphadelta);
            coprimedivnorms:=[1];
            coprimedivminv:=[Codomain(minv)!0];
            for fs in facts do
            	P:=fs[1];
                // if GreatestCommonDivisor(P,ff) eq O1 then
                	normP:=Norm(P);
                    if normP lt BOUND then
                    
                	found,Pinv:=checklist(P,normP,Pnormlist, Pbothlist);
                    
                    if found eq true then // Pinv already computed
                    	// normP:=Norm(P);
                		Ppowsnorm:=[normP^i: i in [0 .. fs[2]]];
                    	// Pinv:=minv(P);
                    	Ppowsinv:=[i*Pinv: i in [0 .. fs[2]]];
                		coprimedivnorms:=[A*Norm(B): A in coprimedivnorms, B in Ppowsnorm];
                		coprimedivminv:=[A + B: A in coprimedivminv, B in Ppowsinv];
                    else
                    	// normP:=Norm(P);
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
                    
            	// end if;
            end for;
                        
            for i:=1 to #coprimedivnorms do
                Append(~Xn,[*alpha,coprimedivnorms[i],coprimedivminv[i]*]);
            end for;
    	end if;
    end for;
    
    // return Xn,delta,Pnormlist,Pbothlist;

end procedure;

// WARNING - I THINK THIS COMPUTES THE UNSTABILISED FORM

function G_kp_seq_faster(D,ks,p,N,NN) // weight 2k

	R:=pAdicRing(p,N);
	// assert IsOdd(k) eq true;
	Zq<q>:=PowerSeriesRing(R,NN);

	K:=QuadraticField(D);
    G,m:=NarrowClassGroup(K);
    minv:=Inverse(m);
    assert Order(G) eq 2;
    assert ClassNumber(K) eq 1;

    O:=Integers(K);
    O1:=ideal<O|1>;
        
    Pnormlist:=[];
    Pbothlist:=[]; // Ps,Pinvs in pairs

	Gk:=[Zq!0: i in [1 .. #ks]]; // sequence of forms
    for n:=1 to NN-1 do
    	Gkn:=[0: i in [1 .. #ks]];
        Xn:=[**];
        computeXn_faster(D,n,minv,O1,~Pnormlist,~Pbothlist,~Xn);
    	for alphaI in Xn do
        	NormI:=alphaI[2];
            // if pdividesI(p,I) eq false then
                J:=alphaI[3];
            	if Order(J) eq 2 then
                	for i:=1 to #ks do
        				Gkn[i]:=Gkn[i]  - (R!NormI)^(ks[i]-1);
                    end for;
            	else
                	for i:=1 to #ks do
                		Gkn[i]:=Gkn[i]  + (R!NormI)^(ks[i]-1);
                    end for;
            	end if;
            // end if;
        end for;
        for i:=1 to #ks do
        	Gk[i]:=Gk[i] + 4*Gkn[i]*q^n;
        end for;
    end for;

	return Gk; // problem is this is not stabilised.

end function;

// ********************************************************** //

// ****** CODE FOR GETTING L-SERIES ****** //

function NewtonSolve(kseq,vseq)

    d:=#kseq;

    Qp:=Parent(vseq[1]);
    aseq:=[Qp!0: i in [1 .. d]];
    
    // create matrix: just seems clearer to see things and avoids null sequences
    
    M:=MatrixRing(Qp,d)!0;
    
    for i:=1 to d do
        M[i,1]:=1;
    end for;

    
    for i:=2 to d do
        for  j:=2 to i do
            M[i,j]:=&*[kseq[i] - kseq[l]: l in [1 .. j-1]];
        end for;
    end for;

    
    aseq[1]:=vseq[1];
    for i:=2 to d do
        aseq[i]:=M[i,i]^-1*(vseq[i] - &+[aseq[j]*M[i,j]: j in [1 .. i-1]]);
    end for;
    
    Aseq:=RSpace(Qp,d)!aseq;
    Vseq:=RSpace(Qp,d)!vseq;
    
    // (Aseq*Transpose(M) - Vseq); // should be zero
    
    
    Qpk<k>:=PolynomialRing(Qp);
    
    // construct Newton polynomials
    
    Nseq:=[Qpk!0: i in [1 .. d]];
    Nseq[1]:=1;
    
    for i:=2 to d do
        Nseq[i]:=&*[(k - kseq[j]): j in [1 .. i-1]];
    end for;
    
    return &+[aseq[i]*Nseq[i]: i in [1 .. d]];
    
end function;

// Input: E is the Katz basis modulo p^mdash, and f a modular form as a POWER SERIES (to precision at
// least the "elldash" in UpGj() modulo p^mdash) which is p/(p+1)-overconvergent.
// Output: f in the Katz basis modulo p^m.
// Note: there is a loss of precision of mdash - m here, so need coefficients of f in IntegerRing(p^mdash)

function qExpToKatzNEWb(f,E,p,m,elldash)

	ell := #E;
    // print ell, elldash;
	PS := Parent(E[1]);
    elldashp:=Precision(PS);
    Qa<a>:=PolynomialRing(Rationals());
    R<q>:=PowerSeriesRing(Qa,elldashp);
    Zq<q>:=PowerSeriesRing(Integers(),elldashp);
    
    f:=R!(Zq!f);
    Z:=Integers();

	
	Rmdashelldash:=RSpace(Qa,elldash);
	w := Rmdashelldash![Z!Coefficient(f,i-1): i in [1..elldash]]; // first elldash coefficients
	Rmdashell:=RSpace(Qa,ell);
	v := Rmdashell![0: i in [1..ell]];
    v[1]:=a;

    
    E1:=Rmdashelldash![Z!Coefficient(E[1],l-1): l in [1 .. elldash]];
    w := w - a*E1;
    
    F:=a*(R!(Zq!E[1])); // the form we're ending up with
	
	for j:=2 to ell do
		Ej:=Rmdashelldash![Z!Coefficient(E[j],l-1): l in [1 .. elldash]]; 
        // Ej;
        // E[j];
        // Ejleadpos:=LeadingPosition(Ej);
        Ejleadpos:=Depth(Ej); // Feb 26, 2013
		lj:=Integers()!(Ej[Ejleadpos]);
		v[j]:=Qa!(w[Ejleadpos]/lj);
		w:=w - v[j]*Ej;
        F:=F + v[j]*(R!(Zq!E[j]));
		// f:=f - v[j]*E[j]; // for checking only
	end for;
	
    /*
	if ValPS(f,p) lt m then // check
		print "Error in qExpToKatz: form not l.c. of katz basis!", ValPS(f,p);
        // IF YOU GET THIS ERROR PROBABLY INPUT CHARACTER IS WRONG ONE.
	end if;
	*/
        	
	Rm:=RSpace(Qa,ell);
    
	
	return Rm!v,F; // loss of precision mdash - m.

end function;

// Converts (p/(p+1))-OC form in Katz Basis E to q-expansion.

function KatzToqExp(E,v)

	ell:=#E;
	qprec:=Precision(Parent(E[1]));
	S:=Parent(v[1]);
	PW<q>:=PowerSeriesRing(S,qprec);
	
	w:=PW!0;
	for i:=1 to ell do
		w:=w + v[i]*(PW!E[i]);
	end for;

	return w;

end function;


function missing_constant_e(G,E,p,m,elldash,qprec)

	dimkatz:=#E; //
	val:=Valuation(E[dimkatz]);

	v,F:=qExpToKatzNEWb(G,E,p,m,elldash);
	Zq:=PowerSeriesRing(Integers(),qprec);
    EE:=[PowerSeriesRing(Parent(v[1]),qprec)!(Zq!e): e in E];

	GG:=KatzToqExp(EE,v);
	Gval:=Coefficient(G,val + 1);
	GGval:=Coefficient(GG,val + 1); // linear form in a
    
    // Valuation(Coefficient(GGval,0),p);

	a:=(Gval - Coefficient(GGval,0))/Coefficient(GGval,1);
	// Z:=Integers();
	Ga:=G + a;
    
    return Ga;

end function;

function Lpk_faster(D,kappa0,p,m) // kappa = k/2 is index of G_kappa

	n:=2*m;
	Qp:=pAdicField(p,m);
	Lpks:=[];
    Z:=Integers();
    //kappajs:=[kappa0 + j*Z!((p-1)/2): j in [0 .. m]];
    kjs:=[2*kappa0 + j*(p-1): j in [1 .. n]];
    
    ellp:=Dimension(ModularForms(1,2*kappa0 + n*(p-1))) + 10; // +1 enough
    // this is just working q precision!
    time Mis:=ComputeAllMi(2*kappa0,p,ellp,m,n); // (k0,p,ellps,mdash,n)
    
    kappaseq:=[kappa0 + j*Z!((p-1)/2): j in [1 .. n]];
    time G_h_seq:=G_kp_seq_faster(D,kappaseq,p,m,ellp);
    // NOTE: THIS IS IN FACT UNSTABILISED FORM, WHICH IS WHAT WE WANT.
    
    t0:=Cputime();
    for j:=1 to n do
    	// time G_h:=G_k(D,kappa0 + j*Z!((p-1)/2),ellp);
        G_h:=G_h_seq[j];
        B:=Mis[j+1]; // basis for M_(2*kappa0 + j*(p-1))
    	G:=missing_constant_e(G_h,B,p,m,ellp,ellp);
	
    	kappaj:=kappa0 + Z!(j*(p-1)/2);
    	Append(~Lpks,Qp!((1 - p^(2*(kappaj - 1)))*Coefficient(G,0)));
        // I don't think the factor is correct for odd j
    end for;
    print Cputime(t0);
    
    return Lpks,kjs;
    
    // Need to interpolate the even js and odd js separately. Latter should give
    // L-series for kappa0 + (p-1)/2 I think.

end function;

// K = Q(sqrt(D)) class number 1, narrow class number 2; p inert in K.

function main_faster(D,kappa0,p,m)

	Ls,ks:=Lpk_faster(D,kappa0,p,m);
    
    t1:=Cputime();
    F:=NewtonSolve([ks[2*i]: i in [1 .. m]],[Ls[2*i]: i in [1 .. m]]);
    G:=NewtonSolve([ks[2*i - 1]: i in [1 .. m]],[Ls[2*i - 1]: i in [1 .. m]]);
    print Cputime(t1);
    
    return F,G; 
    // F is p-adic L-series around kappa0
    // G should be the one around kappa0 + (p-1)/2, I think.
    // WARNING: sometimes G doesn't look anything like an L-series, so not sure what is going on.
    
    
end function;


/*

Tested on D = 33 and p = 5.

Found here relevant 5-unit is root of x^2 - 3*x + 5. This lies in the 5-unit group of the 
narrow class field. F'(2) and log(unit) have ratio something like -3/2.

*/





