/*

SEPTEMBER 9, 2019	

Code to compute all the classical spaces for general level N.

For level 1 this is done in LV-v-4.

*/


//  *** AUXILIARY CODE FOR GENERAL LEVEL ***

// CLASSICAL MODULAR FORMS

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

// Returns dimension of space of classical modular forms for char of weight k.	

function DimensionMNk(char,k)

	if k gt 0 then
		M:=ModularForms(char,k); // space over Q, all conjugates
        deg:=Degree(Parent(char(1))); // number of conjugates
		return Integers()!(Dimension(M)/deg);
	elif ((k eq 0) and Order(char) eq 1) then // trivial character
		return 1;
	else // k < 0 or k = 0 with non-trivial character
		return 0;
	end if;

end function;

// Check kseq contains weights congruent modulo p-1.

function IsValidWeightSequence(kseq,p)

	len:=#kseq;
	if len eq 0 then
		return false;
	end if;
	
	result:=true;
	k0:=kseq[1] mod (p-1);
	for i:=2 to len do
		if (kseq[i] mod (p-1)) ne k0 then
			result:=false;
		end if;
	end for;

	return result;

end function;

// *** LEVEL N > 1 ***

// CLASSICAL MODULAR FORMS AND LINEAR ALGEBRA

// Converts a "modular form" (power series) with coefficients in a cyclotomic field
// C(zeta) where Order(zeta)|(p-1) into a power series with coefficients in
// Z/(p^m), via a fixed choice zeta_pm of embedding of zeta in Zp.
// WARNING: never give a ModFrmElt as input since it computes each coefficient
// "separately" and takes forever - convert it to a power series first.

function qexpToPSinZpm(f,NN,degC,zeta_pm)

    Zpm:=Parent(zeta_pm);
    R<q>:=PowerSeriesRing(Zpm,NN);
    
    fPS:=R!0;
    for i:=0 to NN-1 do
        fi:=Coefficient(f,i);
        fi_pm:=&+[(Zpm!fi[j])*zeta_pm^(j-1): j in [1 .. degC]];
        fPS:=fPS + fi_pm*q^i;
    end for;

    return fPS;

end function;

// Computes a saturated basis for M_k(eps) to precision q^qprec, assuming class
// number of cyclotomic field is one: Steve showed me how to do this.

function IntegralBasis(eps,k,qprec: is_cuspidal:=false) // weight k >= 2 here

    assert k ge 2;
    
    // When modular form space defined over Q
    if (2 mod Order(eps)) eq 0 then 
        // when R = Z the creation of ideals works differently
        // and it is simpler just to do the following ...
        if is_cuspidal eq false then
            M:=ModularForms(eps,k);
        else // cusp forms only
            M:=CuspForms(eps,k);
        end if;
        b:=Basis(M,qprec);
        C:=BaseRing(Parent(eps)); // cyclotomic field
        Cq:=PowerSeriesRing(C,qprec);
        return [Cq!f: f in b]; // 
    end if;
    
    // Now consider case when space not over Q

    // Construct q-expansions of cuspidal space using modular symbols
    MM:=ModularSymbols(eps,k);
    SS:=CuspidalSubspace(MM);
    SSb:=qExpansionBasis(SS,qprec);
    
    // Directly compute eisenstein series
    M:=ModularForms(eps,k);
    if Dimension(M) eq 0 then
        return [];
    end if;
    sturm:=PrecisionBound(M);
    assert qprec ge sturm; // if this fails qprec is set too low
    
    // Put coefficients in a matrix
    if is_cuspidal eq false then
        Es:=EisensteinSeries(M); // this space is defined over Q(eps)
        dim:=#SSb + #Es;
    else // cusp forms only
        dim:=#SSb;
    end if;
    C:=BaseRing(Parent(eps)); // cyclotomic field
    
    A:=ZeroMatrix(C,dim,sturm);
    for i:=1 to #SSb do
        for j:=1 to sturm do
            A[i,j]:=Coefficient(SSb[i],j-1);
        end for;
    end for;
    if is_cuspidal eq false then // full space
        for i:=#SSb + 1 to dim do
            for j:=1 to sturm do
                A[i,j]:=Coefficient(Es[i - #SSb],j-1); 
            end for;
        end for;
    end if;
    
    // Create pseudomatrices for A and R^dim      
    R:=Integers(C);
    I1:=1*R;
    ps_A:=PseudoMatrix([I1: i in [1 .. dim]],A);
    ps_R:=PseudoMatrix(Module(R,sturm)); // R^sturm as pseudo matrix

    ps_Asat:=ps_A meet ps_R; // compute the intersection of the spaces
    // I believe one can also just saturate ps_A directly.

    assert ClassNumber(C) eq 1; // to ensure all ideals principal
    // Steve says you should take second element in TwoElement(Is[i]) otherwise, and
    // this will be a local generator.
    
    Is:=CoefficientIdeals(ps_Asat);
    Is_gen:=[];
    for i:=1 to dim do // find generators for principal ideals
        _,gen:=IsPrincipal(Is[i]);
        Append(~Is_gen,gen);
    end for;
    Asat_vecs:=Matrix(ps_Asat);
    
    Asat:=ZeroMatrix(C,dim,sturm);
    for i:=1 to dim do
        for j:=1 to sturm do
            Asat[i,j]:=Is_gen[i]*Asat_vecs[i,j];
        end for;
    end for;
    
    // Find transformation matrix from A to Asat
    B:=Solution(A,Asat); // B*A eq Asat
    
    // Transform basis elements to full q-adic precision
    b:=[];
    Cq:=PowerSeriesRing(C,qprec);
    if is_cuspidal eq false then // full space
        Mb:=[f: f in SSb] cat [PowerSeries(e,qprec): e in Es];
    else // cusp forms only
        Mb:=[f: f in SSb];
    end if; 
    Mb_sat:=[**];
    for i:=1 to dim do
        f:=Cq!0;
        for j:=1 to dim do
            f:=f + B[i,j]*(Cq!Mb[j]);
        end for;
        Append(~Mb_sat,f);
    end for;
    
    return Mb_sat;
    
end function;

// Code for saturating the weight one basis output by Steve's new function.

function WeightOneSaturatedBasis(eps,qprec: is_cuspidal:=false)

    M2:=ModularForms(eps^2,2); 
    sturm:=PrecisionBound(M2);
    // this will be a sturm bound for weight one too
    qprec:=Maximum(qprec,sturm); // increase precision if too small

    Mb:=ModularFormsBasis(eps,1,qprec: Cuspidal:=is_cuspidal);
    
    if Mb eq [] then // saturation code gives segmentation fault with empty input.
        return [];
    end if;
    
    // Put coefficients in a matrix
    dim:=#Mb;
    C:=BaseRing(Parent(eps));
    
    A:=ZeroMatrix(C,dim,sturm);
    for i:=1 to #Mb do
        for j:=1 to sturm do
            A[i,j]:=Coefficient(Mb[i],j-1);
        end for;
    end for;
        
    if C eq Rationals() then
        Asat:=Saturation(A);
        B:=Solution(A,Parent(A)!Asat); // B*A eq Asat
    else
        R:=Integers(C);
        I1:=1*R;
        ps_A:=PseudoMatrix([I1: i in [1 .. dim]],A);
        ps_R:=PseudoMatrix(Module(R,sturm)); // R^sturm as pseudo matrix

        ps_Asat:=ps_A meet ps_R;

        assert ClassNumber(C) eq 1;
    
        Is:=CoefficientIdeals(ps_Asat);
        Is_gen:=[];
        for i:=1 to dim do // find generators for principal ideals
            _,gen:=IsPrincipal(Is[i]);
            Append(~Is_gen,gen);
        end for;
        Asat_vecs:=Matrix(ps_Asat);
    
        Asat:=ZeroMatrix(C,dim,sturm);
        for i:=1 to dim do
            for j:=1 to sturm do
                Asat[i,j]:=Is_gen[i]*Asat_vecs[i,j];
            end for;
        end for;
        B:=Solution(A,Asat); // B*A eq Asat
    end if;
    
    // constructed saturated code to full precision
    b:=[];
    Cq:=PowerSeriesRing(C,qprec);
    Mb_sat:=[**];
    for i:=1 to dim do
        f:=Cq!0;
        for j:=1 to dim do
            f:=f + B[i,j]*Mb[j];
        end for;
        Append(~Mb_sat,f);
    end for;
    
    return Mb_sat;

end function;

// Returns an integral (saturated) basis for the space of modular forms of weight 
// k and character eps, to q-expansion precision q^qprec. We assume here that
// the class number of Q(eps) is one.
// March 20, 2013: added Cuspidal optional input, to specify only the cusp forms.

function IntegralBasisAllWeights(eps,k,qprec: Cuspidal:=false)

    if k eq 1 then
        return WeightOneSaturatedBasis(eps,qprec: is_cuspidal:=Cuspidal);
    else    // k >= 2
        return IntegralBasis(eps,k,qprec: is_cuspidal:=Cuspidal);
    end if;
    
end function;   

// Constructs lists of integral bases in weight <= weightbound and characters powers of "char".
// These basis are saturated, but we require the class number of Q(char) to be one. 
// The output also includes a matching list of the characters, which are stored
// as their values on a set of generators for the multiplicative group modulo Modulus(char).
// (The characters will be multiplied later on, and this is a better way to store them.) A third
// output is a root of unity used to embed Q(char) into Z_p - note that all modular forms 
// are thought of as power series having coefficients in Q(char), even if the coefficients really
// lie in a smaller field.

function LowWeightBasesWithCharEmbeddedInZp(eps,p,m,NN,weightbound) // m is really mdash here

    if (p-1) mod Order(eps) ne 0 then // Nov 14
		print "ERROR: LowWeightBasesWithChar, spaces not embeddable in Z_p.";
	end if;

	print "mdash,NN = ", m,NN;
	
	generators:=[];
	characters:=[];
    
    C:=Parent(eps(1)); // cyclotomic field or rationals
    degC:=Degree(C);
    if degC gt 1 then
        BasisC:=Basis(C);
        zeta:=BasisC[2];
        assert BasisC eq [zeta^i: i in [0 .. degC-1]];
        Zpm:=pAdicField(p,m);
        IntZpm:=IntegerRing(p^m);
        ZZ:=Integers();
        PolZpm:=PolynomialRing(Zpm);
        zeta_pm:=IntZpm!(ZZ!(Roots(PolZpm!MinimalPolynomial(zeta))[1][1])); // NOTE: chosen FIRST root here
    else // in this case zeta_pm is never used
        if Order(eps) eq 2 then
            zeta_pm:=-1;
        else
            zeta_pm:=1;
        end if;
    end if; 
    
    Cq:=PowerSeriesRing(C,NN);
    
    // February 19, 2013: speeding up multiplication of characters    
    ZZ:=Integers();
    Z_N:=Integers(Modulus(eps));
    U_N,m_N:=UnitGroup(Z_N);
    gens_N:=[ZZ!m_N(u): u in Generators(U_N)]; // generators
    

    S:=IntegerRing(p^m);
    Sq<q>:=PowerSeriesRing(S,NN);
    
	for k:=1 to weightbound do 
        print "weight k = ", k;
		basisweightk:=[];
		charsweightk:=[];
        for i:=0 to Order(eps)-1 do 
            // print "eps^i, i = ", i;
            b:=IntegralBasisAllWeights(eps^i,k,NN);
            randomb:=b;
			if #b gt 0 then // randomisation to remove echelon shape of basis.
				R:=Parent(b[1]);
				dimweightkchari:=#b;
				coeffsmat:=Random(GL(dimweightkchari,p));
				randomb:=[];
				for j:=1 to dimweightkchari do
					coeffsmatj:=coeffsmat[j];
					f:=R!0; 
					for l:=1 to dimweightkchari do
						f:=f + (Integers()!coeffsmatj[l])*b[l];
					end for;
					Append(~randomb,f);
				end for;
			else
				randomb:=b;
			end if;
            if degC eq 1 then
                for f in randomb do
                    Append(~basisweightk,Sq!f);
                    Append(~charsweightk,[(eps^i)(g): g in gens_N]); // represent character by value on generators
                end for;
            else // cyclotomic field
                for f in randomb do
                    Append(~basisweightk,qexpToPSinZpm(f.base_extend(C),NN,degC,zeta_pm));
                    Append(~charsweightk,[(eps^i)(g): g in gens_N]); // represent character by value on generators
                end for;
            end if;
		end for;
		Append(~generators,basisweightk);
		Append(~characters,charsweightk);
        
	end for;
    
	return generators, characters, zeta_pm;

end function;

// The output is a list [[a1,b1],....] where [ai,bi] is the bi character of weight ai,
// and the product of the associated modular forms is a weight "weight" and character "char".
// 19/02/13: use new approach to handling characters where you store them as a list of
// their values on generators - this reduces MASSIVE overhead on multiplying characters.

function RandomSolutionWithChar(characters,weight,char_on_gens,char_on_gens_1)

    // char_on_gens: this is just the characters "char" values on the generators of Z/N^* where
    // N is the modulus. char_on_gens_1: this is just the value of the trivial character on these generators.
	B := #characters;
	found := false;
    
	while found eq false do
		K := weight;
		sol := [];
        charprod:=char_on_gens_1; // just all ones vector of right length.
		a := [];
		// Choose elements in weights B,...,2.
		for i:=B to 2 by -1 do
            if #characters[i] gt 0 then // Nov 15: i.e. there are forms in that weight
                ai := Random(0,Floor(K/i)); // pick ai elements of weight i
                for j:=1 to ai do
                    bij := Random(1,#characters[i]); // characters[i] = chars for weight i
                    charprod:=[charprod[l]*characters[i][bij][l]: l in [1 .. #charprod]]; // #charprod is 1 or 2.
                    Append(~sol,[i,bij]);
                end for;
                K := K - ai*i;
            else
                ai:=0;
            end if;
		end for;
        // Feb 18, 2013: some code which will work even when nothing in weight one
        if #characters[1] gt 0 then // pick K elements in weight one
            for j:=1 to K do
                b1j := Random(1,#characters[1]);
                charprod:=[charprod[l]*characters[1][b1j][l]: l in [1 .. #charprod]]; // Feb 19, 2013
                Append(~sol,[1,b1j]);
            end for;
            if (charprod eq char_on_gens) then
                found := true;
            end if;
        else // nothing in weight one
            if (K eq 0) and (charprod eq char_on_gens) then
                found:=true;
            end if;
        end if;
	end while;
			
	sol:=Reverse(sol);
	
	return sol;

end function;


// COMPLEMENTARY SPACES.

// Auxiliary function used in main loop of ComplementarySpacesModp

function RandomNewBasisModp(p,k,LWBModp,OldBasisModp,weightbound,characters,char,char_on_gens,char_on_gens_1)
    
    R:=Parent(LWBModp[2][1]); // this should be non-empty, since it is the space of weight two forms with
    // trivial character.
    
    // Construct TotalBasisModp
	TotalBasisModp:=OldBasisModp; // Recall E_(p-1) = 1 mod p.
    elldash:=Degree(TotalBasisModp); // Steve 19/02/13: more efficient to use vector spaces than matrices.
	
	// Case k0 + i(p-1) = 0 + 0(p-1) = 0
	if ((k eq 0) and Order(char) eq 1) then // need trivial character and weight 0 to get anything  
		// return [R!1],[[]]; // empty string <-> 1 for NewBasis.
        // TotalBasisModp[1,1]:=1; // this is what we made when using matrices rather than vectors spaces
        TotalBasisModp:=sub< V | V.1 > where V is VectorSpace(GF(p), Degree(TotalBasisModp)); // (1,0,0,...,0)
        return TotalBasisModp, [[]], 0; // March 28: added zero for number of tries.
	elif k eq 0 then // non-trivial character in weight 0.
        return TotalBasisModp, [], 0; // here [] should correspond to "nothing".
	end if;
	
	// Case k = k0 + i(p-1) > 0
	di:=DimensionMNk(char,k);
	diminus1:=DimensionMNk(char,k-(p-1));
    
    mi:=di - diminus1;
    
    // Generate mi new forms in weight k.
	NewBasisCode:=[];
	rk:=diminus1;
    numberoftries:=0; // March 24, 2012: just a counter to see how long things are taking.
	for i:=1 to mi do // extra forms
		while (rk lt diminus1 + i) do
            numberoftries:=numberoftries+1;
            sol:=RandomSolutionWithChar(characters,k,char_on_gens,char_on_gens_1); // 19/02/13
			TotalBasisi:=R!1;
			TotalBasisiCode:=sol;
			for s in sol do
				TotalBasisi:=TotalBasisi*LWBModp[s[1]][s[2]];
			end for;
            // Steve 19/02/13: more efficient way of detecting when "new" modular form is found
            Include(~TotalBasisModp, Vector([Coefficient(TotalBasisi,j): j in [0 .. elldash-1]]), ~new);
            rk:=Dimension(TotalBasisModp);
		end while;
		Append(~NewBasisCode,TotalBasisiCode); // this one increased dimension.
	end for;
		
	return TotalBasisModp,NewBasisCode,numberoftries;

end function;

// Finds complementary spaces modulo p and returns a list of "codes" describing
// what products of basis forms were chosen.

function ComplementarySpacesModp(p,k0,n,elldash,LWBModp,weightbound,characters,char)

    // 19/02/13: change to way characters are handled
    ZZ:=Integers();
    Z_N:=Integers(Modulus(char));
    U_N,m_N:=UnitGroup(Z_N);
    gens_N:=[ZZ!m_N(u): u in Generators(U_N)];
    char_on_gens:=[char(u): u in gens_N];
    char_on_gens_1:=[1: u in Generators(U_N)];  // trivial character

	CompSpacesCode:=[];
    ell:=DimensionMNk(char,k0 + n*(p-1)); 
    OldBasisModp:=sub< VectorSpace(GF(p),elldash) | >; // Steve 19/02/13
    print "n = ", n;
    totalnumberoftries:=0;
	for i:=0 to n do
        print "Computing basis for M_(k0 + i(p-1)) mod p for k0+i(p-1),i,n:", k0+i*(p-1),i,n;
        TotalBasisModp,NewBasisCodemi,numberoftriesi:=RandomNewBasisModp(p,k0 + i*(p-1),LWBModp,OldBasisModp,weightbound,characters,char,char_on_gens,char_on_gens_1);
		Append(~CompSpacesCode,NewBasisCodemi);
		OldBasisModp:=TotalBasisModp; // basis for M_(k0 + i(p-1))
        totalnumberoftries:=totalnumberoftries + numberoftriesi;
	end for;
    print "ell, totalnumberoftries: ", ell, totalnumberoftries;

	return CompSpacesCode;

end function;

// Reduces the basis of low weight forms mod (p,q^elldash).

function LowWeightBasesModp(LWB,p,elldash)

	R:=PowerSeriesRing(GF(p),elldash);
	LWBModp:=[];	
	
	for i:=1 to #LWB do // weight k = i
		LWBModpWeightk:=[];
		for j:=1 to #LWB[i] do
			Append(~LWBModpWeightk,R!LWB[i][j]);
		end for;
		Append(~LWBModp,LWBModpWeightk);
	end for;
	
	return LWBModp;

end function;

// Returns complementary spaces W = [W_0,W_1,...,W_n] as list of 
// basis elements W_i modulo (p^mdash,q^elldashp).

function ComplementarySpaces_A(p,k0,n,mdash,elldash,elldashp,weightbound,char,eps)

	// Find q-expansions for k <= weightbound mod (p^mdash,q^elldashp)
	
	t0:=Cputime();
	LWB,characters,zeta_pm := LowWeightBasesWithCharEmbeddedInZp(eps,p,mdash,elldashp,weightbound);
	print "Time to compute low weight basis:", Cputime(t0);
    // LWB,characters;

	LWBModp:=LowWeightBasesModp(LWB,p,elldash);
	t1:=Cputime();
	CompSpacesCode:=ComplementarySpacesModp(p,k0,n,elldash,LWBModp,weightbound,characters,char);
    // CompSpacesCode;
	print "Time to find complementary spaces modulo p:", Cputime(t1);
		
	// Reconstruct products mod (p^mdash,q^elldashp)
	
	W:=[];
	
	Epm1:=ESeries(p-1,elldashp,IntegerRing(p^mdash));

	t2:=Cputime();
	OldBasis:=[];
	for i:=0 to n do
		CompSpacesCodemi:=CompSpacesCode[i+1];
		Wi:=[];
		for k:=1 to #CompSpacesCodemi do
			// "code" <-> [ ... (k/2,a) ...] where k = weight, a = index of element
			// of weight k chosen, and then one takes the product over the list.
            // Nov 24: "code" <-> [ .... (k-1,a) ... ].
			CompSpacesCodemik:=CompSpacesCodemi[k];	// this is a "sol"	
			Wik:=Parent(Epm1)!1;
			for j:=1 to #CompSpacesCodemik do 
                kminus1:=CompSpacesCodemik[j][1]; // 26/03/13: I think this is really "k" now.
				index:=CompSpacesCodemik[j][2];
                Wik:=Wik*LWB[kminus1,index];
			end for;
			Append(~Wi,Wik);
		end for;
		Append(~W,Wi); 
	end for;
	
	print "Constructed complementary spaces in time:", Cputime(t2);

	return W, zeta_pm, Epm1; // Nov 15, 2012: also return element used to embed cyclotomic ring in Zp.

end function;

// 26/03/13: new version "B" complementary spaces code

/*

In this new version we append the previous spaces constructed to the LWB and use them. This
is quite complicated.

*/

function RandomSolutionWithChar_B(characters,char_on_gens,char_on_gens_1,k0,j,p,weightbound)

    weight:=k0 + j*(p-1); // 08/04/03 - this is the target weight.
    // char_on_gens: this is just the characters "char" values on the generators of Z/N^* where
    // N is the modulus. char_on_gens_1: this is just the value of the trivial character on these generators.
	B := weightbound + (j-1); // start choosing from this position.
	found := false;
    
	while found eq false do
		K := weight;
		sol := [];
        charprod:=char_on_gens_1; // just all ones vector of right length.
		a := [];
		// Choose elements in weights B,...,2.
		for i:=B to 2 by -1 do
            if #characters[i] gt 0 then // Nov 15: i.e. there are forms in that weight
                if i le weightbound then
                    ai := Random(0,Floor(K/i)); // pick ai elements of weight i
                    K := K - ai*i;
                else 
                    j_i:=i - weightbound; // so form in position i is weight k0 + j_i(p-1)
                    ai := Random(0,Floor(K/(k0 + j_i*(p-1))));
                    K := K - ai*(k0 + j_i*(p-1)); 
                end if; 
                for m:=1 to ai do // surely bad use of j here since j is an input? 08/04/13 j -> m
                    bim := Random(1,#characters[i]); // characters[i] = chars for weight i
                    charprod:=[charprod[l]*characters[i][bim][l]: l in [1 .. #charprod]]; // #charprod is 1 or 2.
                    Append(~sol,[i,bim]);
                end for;
            else
                ai:=0;
            end if;
		end for;
        // Feb 18, 2013: some code which will work even when nothing in weight one
        if #characters[1] gt 0 then // pick K elements in weight one
            for m:=1 to K do // 08/04/13 j -> m
                b1m := Random(1,#characters[1]);
                charprod:=[charprod[l]*characters[1][b1m][l]: l in [1 .. #charprod]]; // Feb 19, 2013
                Append(~sol,[1,b1m]);
            end for;
            if (charprod eq char_on_gens) then
                found := true;
            end if;
        else // nothing in weight one
            if (K eq 0) and (charprod eq char_on_gens) then
                found:=true;
            end if;
        end if;
	end while;
			
	sol:=Reverse(sol);
	
	return sol;


end function;


// Auxiliary function used in main loop of ComplementarySpacesModp

function RandomNewBasisModp_B(p,k0,j,LWBModp,OldBasisModp,weightbound,characters,char,char_on_gens,char_on_gens_1)
    
    k:=k0 + j*(p-1); // 08/04/13
    
    R:=Parent(LWBModp[2][1]); // this should be non-empty, since it is the space of weight two forms with
    // trivial character.
    
    // Construct TotalBasisModp
	TotalBasisModp:=OldBasisModp; // Recall E_(p-1) = 1 mod p.
    elldash:=Degree(TotalBasisModp); // Steve 19/02/13: more efficient to use vector spaces than matrices.
	
	// Case k0 + i(p-1) = 0 + 0(p-1) = 0
	if ((k eq 0) and Order(char) eq 1) then // need trivial character and weight 0 to get anything  
		// return [R!1],[[]]; // empty string <-> 1 for NewBasis.
        // TotalBasisModp[1,1]:=1; // this is what we made when using matrices rather than vectors spaces
        TotalBasisModp:=sub< V | V.1 > where V is VectorSpace(GF(p), Degree(TotalBasisModp)); // (1,0,0,...,0)
        return TotalBasisModp, [[]], 0; // March 28: added zero for number of tries.
	elif k eq 0 then // non-trivial character in weight 0.
        return TotalBasisModp, [], 0; // here [] should correspond to "nothing".
	end if;
	
	// Case k = k0 + i(p-1) > 0
	di:=DimensionMNk(char,k);
	diminus1:=DimensionMNk(char,k-(p-1));
    
    mi:=di - diminus1;
    
    // March 27, 2013:
    
    // j:=Floor((k - k0)/(p-1)); // k:=k0 + j(p-1). // 08/04/13
    
    // Generate mi new forms in weight k.
	NewBasisCode:=[];
	rk:=diminus1;
    numberoftries:=0; // March 24, 2012: just a counter to see how long things are taking.
	for i:=1 to mi do // extra forms
		while (rk lt diminus1 + i) do
            numberoftries:=numberoftries+1;
            sol:=RandomSolutionWithChar_B(characters,char_on_gens,char_on_gens_1,k0,j,p,weightbound); // 19/02/13
			TotalBasisi:=R!1;
			TotalBasisiCode:=sol;
			for s in sol do
				TotalBasisi:=TotalBasisi*LWBModp[s[1]][s[2]];
			end for;
            // Steve 19/02/13: more efficient way of detecting when "new" modular form is found
            Include(~TotalBasisModp, Vector([Coefficient(TotalBasisi,j): j in [0 .. elldash-1]]), ~new);
            rk:=Dimension(TotalBasisModp);
		end while;
		Append(~NewBasisCode,TotalBasisiCode); // this one increased dimension.
        if j gt 0 then // this is case k0 + j(p-1) with j > 0
            Append(~LWBModp[weightbound + j],TotalBasisi); // add this new q-expansion to LWBModp
            Append(~characters[weightbound + j],char_on_gens); // add character "char" to list of characters
        end if;
	end for;
		
	return TotalBasisModp,NewBasisCode,LWBModp,characters,numberoftries;

end function;

// Finds complementary spaces modulo p and returns a list of "codes" describing
// what products of basis forms were chosen.

function ComplementarySpacesModp_B(p,k0,n,elldash,LWBModp,weightbound,characters,char)

    // 19/02/13: change to way characters are handled
    ZZ:=Integers();
    Z_N:=Integers(Modulus(char));
    U_N,m_N:=UnitGroup(Z_N);
    gens_N:=[ZZ!m_N(u): u in Generators(U_N)];
    char_on_gens:=[char(u): u in gens_N];
    char_on_gens_1:=[1: u in Generators(U_N)];  // trivial character

	CompSpacesCode:=[];
    ell:=DimensionMNk(char,k0 + n*(p-1)); 
    OldBasisModp:=sub< VectorSpace(GF(p),elldash) | >; // Steve 19/02/13
    print "n = ", n;
    totalnumberoftries:=0;
	for i:=0 to n do
        print "Computing basis for M_(k0 + i(p-1)) mod p for k0+i(p-1),i,n:", k0+i*(p-1),i,n;
        // Here we update both LWBModp and characters, appending the new q-expansions found.
        TotalBasisModp,NewBasisCodemi,LWBModp,characters,numberoftriesi:=RandomNewBasisModp_B(p,k0,i,LWBModp,OldBasisModp,weightbound,characters,char,char_on_gens,char_on_gens_1);
		Append(~CompSpacesCode,NewBasisCodemi);
		OldBasisModp:=TotalBasisModp; // basis for M_(k0 + i(p-1))
        totalnumberoftries:=totalnumberoftries + numberoftriesi;
	end for;
    print "ell, totalnumberoftries: ", ell, totalnumberoftries;

	return CompSpacesCode;

end function;



function ComplementarySpaces_B(p,k0,n,mdash,elldash,elldashp,weightbound,char,eps)

    // Find q-expansions for k <= weightbound mod (p^mdash,q^elldashp)
	
	t0:=Cputime();
	LWB,characters,zeta_pm := LowWeightBasesWithCharEmbeddedInZp(eps,p,mdash,elldashp,weightbound);
	print "Time to compute low weight basis:", Cputime(t0);
    // LWB,characters;

	LWBModp:=LowWeightBasesModp(LWB,p,elldash);
	t1:=Cputime();
    
    // March 27, 2013: add extra lists to be filled in during the algorithm
    
    LWBModp:=LWBModp cat [[]: j in [1  .. n]];
    LWB:=LWB cat [[]: j in [1 .. n]];
    characters:=characters cat [[]: j in [1 .. n]];

    
	CompSpacesCode:=ComplementarySpacesModp_B(p,k0,n,elldash,LWBModp,weightbound,characters,char);
    // CompSpacesCode;
	print "Time to find complementary spaces modulo p:", Cputime(t1);
    
        
	// Reconstruct products mod (p^mdash,q^elldashp)
	
	W:=[];
	
	Epm1:=ESeries(p-1,elldashp,IntegerRing(p^mdash));

	t2:=Cputime();
	OldBasis:=[];
	for i:=0 to n do
		CompSpacesCodemi:=CompSpacesCode[i+1];
		Wi:=[];
		for k:=1 to #CompSpacesCodemi do
            // Mar 27, 2013: pretty sure in version A "code" <-> [ ... (k,a) ...]
			CompSpacesCodemik:=CompSpacesCodemi[k];	// this is a "sol"	
			Wik:=Parent(Epm1)!1;
			for j:=1 to #CompSpacesCodemik do 
                kminus1:=CompSpacesCodemik[j][1]; // 26/03/13: I think this is really "k" now.
				index:=CompSpacesCodemik[j][2];
                Wik:=Wik*LWB[kminus1,index];
                // In B the (k,a) in "code" will correspond to a form in weight k when 1 <= k <= weightbound
                // but otherwise in weight (k - weightbound)*(p-1) + k0.
                // In any case, the product of these forms will have the right weight and character.
			end for;
			Append(~Wi,Wik);
            if i gt 0 then // add new form to low weight basis.
                Append(~LWB[weightbound+i],Wik);
            end if;
		end for;
		Append(~W,Wi); 
	end for;
	
	print "Constructed complementary spaces in time:", Cputime(t2);

	return W, zeta_pm, Epm1; // Nov 15, 2012: also return element used to embed cyclotomic ring in Zp.
    
end function;


// 26/03/13: comp:="A"(old) or "B"(new), method for constructing complementary spaces.

function ComplementarySpaces(p,k0,n,mdash,elldash,elldashp,weightbound,char,eps,comp)

    if comp eq "A" then
        return ComplementarySpaces_A(p,k0,n,mdash,elldash,elldashp,weightbound,char,eps);
    else
        // This version is *much less well tested* but appears to be faster
        return ComplementarySpaces_B(p,k0,n,mdash,elldash,elldashp,weightbound,char,eps);
    end if;

end function;

// Returns ldash, the Sturm bound.

function Computeelldash_new(p,char,k0,n) // INPUT n rather than m

	// n:=Floor(((p+1)/(p-1))*(m+1));
	N:=Modulus(char);
	// From Page 173-174 Stein: Corollary 9.19 and 9.20
	ind:=Index(Gamma0(N));
	elldash:=Floor(((k0 + n*(p-1))*ind)/12) + 1;
	// This is a sturm bound for M(Gamma0(N),k), and hence also
	// for M(char,k) by Corollary 9.20.

	return elldash;

end function;


// END CODE FROM DLR-V-8/pii_B_300513.m


// ***********************************************


// HERE k0 IS THREE TIMES WHAT YOU THINK.
// ALSO "n = 3 n - 2", i.e., nearly three times what you think
function ComputeAllMi(k0,p,m,ellp,n,eps) // ellp is just the q-precision desired

	// Wis,Ep1:=ComputeAllWi(k0,p,ellp,n);
    
    char:=eps;
    elldash:=Computeelldash_new(p,char,k0,n); // a Sturm bounded on biggest space
    elldash;
    // assert ellp gt elldash;
	// elldashp:=elldash*p;
    elldashp:=Maximum(ellp,elldash); // just the q-precision needed here
	// n:=Floor(((p+1)/(p-1))*(m+1));	 
	// mdash:=m + Ceiling(n/(p+1));
    mdash:=m;
    
    weightbound:=4; // 3 -> 4 Oct 1, 2019
    comp:="B";
    Wis,zetapm,Ep1:=ComplementarySpaces(p,k0,n,mdash,elldash,elldashp,weightbound,char,eps,comp);
    
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

// ********

function ComputeAllMi_Naive(k0,n,p,eps)

	qprec:=PrecisionBound(ModularForms(eps,k0 + n*(p-1)));
	Ms:=[**];
    
    for j:=0 to n do
    	M:=ModularForms(eps,k0 + j*(p-1));
        B:=Basis(M,qprec);
        Append(~Ms,B);
    end for;

	return Ms;


end function;

// Sept 13, 2019





// Sept 13, 2019: 

/*

Code below is for Jan's cubic field example.

G:=DirichletGroup(5);
eps:=G.1;
n:=3;
p:=7;
m:=15; // precision
fs:=main(p,n+1); // +1 just to make sure we have enough coefficients, otherwise things go wrong!!
qprec:=Precision(Parent(fs[1]));
time Ms:=ComputeAllMi(3*2,7,m,qprec,3*n - 2,eps); // note the 6, and also have to 

Now for p = 7 we have (7) is inert in the field and psi((7)) = -1 again. So
Euler factor should be

(1 + (7^3)^((6*i - 4)-1))

kseq:=[2 + (j-1)*(p-1): j in [1 .. 10]];

This works - I mean getting the L-series



*/

// Check that a given q-expansion is (modulo the constant term) a member of the space 
// M_k(chi) of classical forms of weight k.
function IsClassicalForm_padic(Gk,Ms,j)

	m := Precision(Parent(Gk)); 
	// B := Basis(ModularForms(Gamma1(N),k),m);
    B:=Ms[3*j - 2];
    // B;
	print "Classical space is of dimension ", #B;
	// Matrix of (higher) coefficients
	C := Matrix([[Coefficient(B[i],j) : j in [1..m-1]] : i in [1..#B]] cat [[Coefficient(Gk,j) : j in [1..m-1]]]);
	K := Basis(Kernel(C));
	if #K eq 0 then
		print "Not a classical form!";
		CT := 0;
	else
    	K;
        print "#K = ", #K;
        // remove assertion because could be error forms in basis of high valuation.
		// assert #K eq 1; // There must be a unique linear combination.
		comb := K[1];
// 		print "Classical form, linear combination is: ", K[1];
// 		print "k is: ", k;
		CT := -(&+ [Coefficient(B[i],0)*comb[i] : i in [1..#B]])/comb[#B+1];
	end if;
	
	return CT;

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









// ***************

// Sept 27, 2019: Code added by Alan.



// computes upper bound on degree of P_k(psi) mod p^m
function findmdash(p,m)

	mdash:=1;
    
    while mdash - Valuation(Factorial(mdash),p) lt m do
    	mdash:=mdash + 1;
    end while;

	return mdash - 1;

end function;



// function Ls_Ks(F,m,eps,p) // eps is character of the diagonal restrictions
// 
// 	mdash:=findmdash(p,m);
//     Ks:=[2 + j*(p-1): j in [0 .. mdash-1+1]]; // Sept 19 +1
//     
//     n:=PrecisionBound(ModularForms(eps,3*Ks[#Ks])) + 1; 
//     
//     fs:=Diagonal_Restriction_Ks(F,Ks,n);
//     
// 	// return fs;
//     
//     qprec:=Precision(Parent(fs[1])); // precision is 91 for m = 10
// 	time Ms:=ComputeAllMi(6,p,m,qprec,3*(mdash+1) - 2,eps);
//     
//     Ls:=[**];
// 	for i:=1 to mdash+1 do // Sept19 + 1
// 		Append(~Ls,IsClassicalForm_padic(fs[i],Ms,i));
// 	end for;
//     
//     return Ls,Ks; // need to put Euler factor at p in.
// 
// end function;

/// September 30, 2019: Above function works but just trying to speed it up.

function Diagonal_Restriction_Ks2(F,Ks,m : Stabilisation := 1)

	As, Forms := Diagonal_Restriction_data(F,m);
	
    // return diagonal_restriction(F,k,m,As,Forms : Stabilisation := Stabilisation);
	
    R<q>   := PowerSeriesRing(ZZ,m);
	Diag_Fs := [R!0: k in Ks];

	f := Conductor(Parent(F));
	p := Stabilisation;
	assert IsPrime(p) or p eq 1;
	
	for n in [1..m-1] do
		A_p := As[n][1];
		A_m := As[n][2];
		coeff_ns := [R!0: k in Ks];
		for i := 1 to #A_p do
			if Gcd(A_p[i],p*f) eq 1 then
            	for j:=1 to #Ks do
					coeff_ns[j] := coeff_ns[j] + A_p[i]^(Ks[j]-1);
                end for;
			end if;
			if Gcd(A_m[i],p*f) eq 1 then
            	for j:=1 to #Ks do
					coeff_ns[j] := coeff_ns[j] + (-1)^Ks[j]*(-A_m[i])^(Ks[j]-1);
                end for;
			end if;
		end for;
		
        for j:=1 to #Ks do
			Diag_Fs[j] := Diag_Fs[j] + 4*coeff_ns[j]*q^n;
        end for;
		
	end for;
	
	return Diag_Fs;

end function;



