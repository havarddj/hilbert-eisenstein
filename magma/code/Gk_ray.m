/*z
June 21, 2018

Code to compute the series G_k(psi), except the constant term. Note this
is not p-stabilised ... in fact it is not in general an eigenform I think, so
have to be careful about this "stabilisation" business.

K:=QuadraticField(D); // class number 1 and narrow class number 2, e.g. D = 12
k must be ODD (or can take k even, but then need psi the trivial map on ideals)

October 16, 2018

Based on Gk in DLV folder.

Use this example:

K:=QuadraticField(5);
O:=Integers(K);
pp:=19;
P:=Factorisation(ideal<O|pp>)[1][1];
ff:=P*ideal<O|2>;
G,m:=RayClassGroup(ff);
// InfinitePlaces(K);


*/




// Oct 16, 2018: I think the same function is fine here

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

// Oct 16, 2018: Will modify this one.

function G_k(D,k,NN) // weight 2k

	// assert IsOdd(k) eq true;
    C<z3>:=CyclotomicField(3);
	Cq<q>:=PowerSeriesRing(C,NN);

	/*
	K:=QuadraticField(D);
    G,m:=NarrowClassGroup(K);
    minv:=Inverse(m);
    assert Order(G) eq 2;
    assert ClassNumber(K) eq 1;
	*/
    
    
    K:=QuadraticField(5);
	O:=Integers(K);
	pp:=19;
	P:=Factorisation(ideal<O|pp>)[2][1];
	ff:=P*ideal<O|2>;
	G,m:=RayClassGroup(ff);
    minv:=Inverse(m);

    Gk:=Cq!0;
    O1:=ideal<O|1>;
    
    for n:=1 to NN-1 do
    	Gkn:=0;
        Xn:=computeXn(D,n);
    	for alphaI in Xn do
        	I:=alphaI[2];
            gcd:=GreatestCommonDivisor(I,ff);
            if gcd eq O1 then
            	J:=minv(I);
            	if J eq G.1 then
        			Gkn:=Gkn  + z3*Norm(I)^(k-1); 
            	elif J eq 2*G.1 then
                	Gkn:=Gkn  + z3^2*Norm(I)^(k-1);
            	elif J eq 3*G.1 then // zero
            		Gkn:=Gkn  + z3^2*Norm(I)^(k-1);
            	end if;
            end if;
        end for;
        Gk:=Gk + 4*Gkn*q^n;
    end for;

	return Gk;

end function;



function G_k2(D,k,NN) // weight 2k

	// assert IsOdd(k) eq true;
    C<i>:=QuadraticField(-1);
	Cq<q>:=PowerSeriesRing(C,NN);

    K:=QuadraticField(5);
	O:=Integers(K);
	pp:=29;
	P:=Factorisation(ideal<O|pp>)[2][1];
	ff:=P;
	G,m:=RayClassGroup(ff,[1,2]);
    minv:=Inverse(m);

    Gk:=Cq!0;
    O1:=ideal<O|1>;
    
    for n:=1 to NN-1 do
    	Gkn:=0;
        Xn:=computeXn(D,n);
    	for alphaI in Xn do
        	I:=alphaI[2];
            gcd:=GreatestCommonDivisor(I,ff);
            if gcd eq O1 then
            	J:=minv(I);
            	if J eq G.1 then
        			Gkn:=Gkn  + i*Norm(I)^(k-1); 
            	elif J eq 2*G.1 then
                	Gkn:=Gkn  + i^2*Norm(I)^(k-1);
            	elif J eq 3*G.1 then // zero
            		Gkn:=Gkn  + i^3*Norm(I)^(k-1);
                elif J eq 4*G.1 then // zero
            		Gkn:=Gkn  + i^4*Norm(I)^(k-1);
            	end if;
            end if;
        end for;
        Gk:=Gk + 4*Gkn*q^n;
    end for;

	return Gk;

end function;

function MM(f,k,NN) // no good function

	C<i>:=QuadraticField(-1);
	Cq<q>:=PowerSeriesRing(C,NN);

	M:=ModularForms(145,2*k);
    
    B:=Basis(M);
    
    b:=[Cq!PowerSeries(m,NN): m in Basis(M)];
    
    bb:=[LeadingCoefficient(m)^(-1)*m: m in b];
    
    F:=f;
    for i:=1 to Dimension(M) - 1 do
    	F:=F - Coefficient(F,i)*bb[i+1];
    end for;
    
    // F:=&+[Coefficient(f,i)*LeadingCoefficient(b[i+1])^(-1)*b[i+1]: i in [1 .. Dimension(M) - 1]];
    
	return F,bb;

end function;



function G_k3(k,NN) // weight 2k

	D:=13;
	// assert IsOdd(k) eq true;
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
        Xn:=computeXn(D,n);
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

function MM3(f,k,NN) // no good function

	C<i>:=QuadraticField(-1);
	Cq<q>:=PowerSeriesRing(C,NN);

	M:=ModularForms(Gamma1(11),2*k);
    
    B:=Basis(M);
    
    b:=[Cq!PowerSeries(m,NN): m in Basis(M)];
    
    bb:=[LeadingCoefficient(m)^(-1)*m: m in b];
    
    F:=f;
    for i:=1 to Dimension(M) - 1 do
    	F:=F - Coefficient(F,i)*bb[i+1];
    end for;
    
    // F:=&+[Coefficient(f,i)*LeadingCoefficient(b[i+1])^(-1)*b[i+1]: i in [1 .. Dimension(M) - 1]];
    
	return F,bb;

end function;



