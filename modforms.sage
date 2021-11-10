####
# Algorithm for computing M_{k_0 + j(p-1)}(char) mod p^n, q^elldash, j=0,1,...,
#
#
# Ported from Jan Vonk & Alan Lauder's magma file 'ClassicalBases.m'
# Similar to the one in sage's hecke_series.py, but with nebentypus
#
#
###

# ClassicalBases.m lines 866-898
def ComputeAllMi(k0, p, m, ellp, n, eps):  # ellp is just the q-precision desired
    # // Wis,Ep1:=ComputeAllWi(k0,p,ellp,n);

    char = eps
    # // a Sturm bounded on biggest space
    elldash = Computeelldash_new(p, char, k0, n)
    # elldash;
    # assert ellp > elldash
    # // elldashp:=elldash*p;
    elldashp = max(ellp, elldash)  # just the q-precision needed here
    # // n = floor(((p+1)/(p-1))*(m+1));
    # // mdash:=m + Ceiling(n/(p+1));
    mdash = m

    weightbound = 4  # // 3 -> 4 Oct 1, 2019
    # comp = "B"
    Wis, zetapm, Ep1 = ComplementarySpaces_B(
        p, k0, n, mdash, elldash, elldashp, weightbound, char, eps)

    Mis = []
    for i in [0..n]:
        Mi = []
        for j in [0..i]:
            Ep1j = Ep1**(i-j)
            for e in Wis[j]:
                Mi.append(Ep1j*e)
        Mis.append(Mi)

    return(Mis)

# end function;

# ClassicalBases.m lines 843-855


def Computeelldash_new(p, char, k0, n):  # INPUT n rather than m

    # n = floor(((p+1)/(p-1))*(m+1))
    N = char.modulus()
    # // From Page 173-174 Stein: Corollary 9.19 and 9.20
    ind = Gamma0(N).index()
    elldash = floor(((k0 + n*(p-1))*ind)/12) + 1
    # // This is a sturm bound for M(Gamma0(N),k), and hence also
    # // for M(char,k) by Corollary 9.20.
    return elldash
# end function;


# ClassicalBases.m lines 766-825
def ComplementarySpaces_B(p, k0, n, mdash, elldash, elldashp, weightbound, char, eps):

    # // Find q-expansions for k <= weightbound mod (p^mdash,q^elldashp)

    # t0:=Cputime();
    LWB, characters, zeta_pm = LowWeightBasesWithCharEmbeddedInZp(
        eps, p, mdash, elldashp, weightbound)

    # print "Time to compute low weight basis:", Cputime(t0);
    # // LWB,characters;

    LWBModp = LowWeightBasesModp(LWB, p, elldash)

    # // March 27, 2013: add extra lists to be filled in during the algorithm

    LWBModp = LWBModp + [[] for j in [1..n]]
    LWB = LWB + [[] for j in [1 .. n]]
    characters = characters + [[] for j in [1..n]]

    CompSpacesCode = ComplementarySpacesModp_B(p, k0, n, elldash,
                                               LWBModp, weightbound, characters, char)

    # Reconstruct products mod (p^mdash,q^elldashp)
    W = []
    Epm1 = ESeries(p-1, elldashp, Integers(p**mdash))
    OldBasis = []
    for i in [0..n]:
        CompSpacesCodemi = CompSpacesCode[i]
        Wi = []
        for k in [1..len(CompSpacesCodemi)]:
            # // Mar 27, 2013: pretty sure in version A "code" <-> [ ... (k,a) ...]
            CompSpacesCodemik = CompSpacesCodemi[k]  # // this is a "sol"
            Wik = Epm1.parent()(1)
            for j in [1..len(CompSpacesCodemik)]:
                # 26/03/13: I think this is really "k" now.
                kminus1 = CompSpacesCodemik[j][1]
                index = CompSpacesCodemik[j][2]
                Wik = Wik*LWB[kminus1, index]
                #  In B the (k,a) in "code" will correspond to a form
                # in weight k when 1 <= k <= weightbound but
                # otherwise in weight (k - weightbound)*(p-1) + k0.
                # In any case, the product of these forms will have
                # the right weight and character.
            # end for;
            Wi.append(Wik)
            if i > 0:                 # // add new form to low weight basis.
                LWB[weightbound+i].append(Wik)
                # end if;
                # end for;
        W.append(Wi)
        # end for;
        # print "Constructed complementary spaces in time:", Cputime(t2)
    # // Nov 15, 2012: also return element used to embed cyclotomic ring in Zp.
    return([W, zeta_pm, Epm1])

# end function;
# ClassicalBases.m lines 19-35


def ESeries(k, NN, S):
    R.<q> = PowerSeriesRing(S, default_prec=NN)
    a1 = -S(2*k/bernoulli(k))
    Ek = 1 + a1*q
    for n in [2..NN-1]:
        coeffn = 0
        divs = divisors(n)
        for d in divs:
            coeffn += d**(k-1)
            # end for;
        Ek += Ek + a1*coeffn*q**n
        # end for;

    return Ek


# end function;
# ClassicalBases.m lines 308-393


# m is really mdash here
def LowWeightBasesWithCharEmbeddedInZp(eps, p, m, NN, weightbound):

    if (p-1) % eps.order() != 0:  # // Nov 14
        print("ERROR: LowWeightBasesWithChar, spaces not embeddable in Z_p.")
        # end if;

    print("mdash, NN = ", m, NN)

    generators = []
    characters = []

    C = eps(1).parent()  # cyclotomic field or rationals
    degC = C.degree()
    if degC > 1:
        BasisC = C.power_basis()
        zeta = BasisC[1]
        assert BasisC == [zeta**i for i in [0 .. degC-1]]

        Zpm = Qp(p, m)          # ???
        IntZpm = Integers(p**m)
        ZZ = Integers()
        PolZpm = PolynomialRing(Zpm,'T')

        # NOTE: chosen FIRST root here:
        try:
            zeta_pm =IntZpm(ZZ(PolZpm(zeta.minpoly()).roots(multiplicities=False)[0]))
            print('zeta_pm ok')
        except: 
            zeta_pm = IntZpm(1)
            print('zeta_pm not really ok')
    else:  # in this case zeta_pm is never used
        if eps.order() == 2:
            zeta_pm = -1
        else:
            zeta_pm = 1
            # end if;
            # end if;

    Cq = PowerSeriesRing(C,'Q',default_prec=NN)

    # February 19, 2013: speeding up multiplication of characters
    # ZZ:=Integers();
    Z_N = Integers(eps.modulus())
    # U_N,m_N:=UnitGroup(Z_N);
    gens_N = [ZZ(u) for u in Z_N.unit_gens()]  # generators

    S = Integers(p^m)
    Sq = PowerSeriesRing(S,'q',default_prec=NN)

    for k in [1..weightbound]:
        print("weight k =", k)
        basisweightk = []
        charsweightk = []
        for i in [0..eps.order()-1]:
            # print("eps^i, i =", i)
            # b = IntegralBasisAllWeights(eps ^ i, k, NN)
            b = ModularForms(chi**i,k,prec=NN).basis()
            randomb = b
            if len(b) > 0:  # randomisation to remove echelon shape of basis.
                R = b[1].parent()
                dimweightkchari = len(b)
                coeffsmat = GL(dimweightkchari, p).random_element()
                randomb = []
                for j in [0..dimweightkchari-1]:
                    coeffsmatj = coeffsmat.matrix()[j]
                    f = R(0)
                    for l in [0..dimweightkchari-1]:
                        f += coeffsmatj.change_ring(ZZ)[l]*b[l]
                        # end for;
                    # make sure all are viewed as forms
                    # over appropriate cyclotomic extension:
                    randomb.append(f.base_extend(C))
                # end for;
            else:
                randomb = b
                # end if;
            if degC == 1:
                for f in randomb:
                    basisweightk.append(Sq(f))
                    charsweightk.append([(eps^i)(g) for g in gens_N])
                    # represent character by value on generators
                    # end for;
            else:  # // cyclotomic field
                for f in randomb:
                    basisweightk.append(qexpToPSinZpm(f, NN, degC, zeta_pm))
                    # represent character by value on generators
                    charsweightk.append([(eps^i)(g) for g in gens_N])
                    # end for;
                    # end if;
        # end for;
        generators.append(basisweightk)
        characters.append(charsweightk)

    # end for;

    return([generators, characters, zeta_pm])

# end function;

# ClassicalBases.m lines 289-297
# def IntegralBasisAllWeights(eps,k,qprec: Cuspidal:=false):
#     if k == 1:
#         return WeightOneSaturatedBasis(eps,qprec: is_cuspidal:=Cuspidal);
#     else:     # k >= 2
#         return IntegralBasis(eps,k,qprec: is_cuspidal:=Cuspidal);
#     # end if;
    
# # end function;   

# # ClassicalBases.m lines 103-210
# def IntegralBasis(eps, k, qprec, is_cuspidal=false) # weight k >= 2 here

#     assert k > 1;
    
#     # // When modular form space defined over Q
#     if (2 mod eps.order()) == 0: 
#         # when R = Z the creation of ideals works differently
#         #  and it is simpler just to do the following ...
#         if is_cuspidal == false:
#             M = ModularForms(eps,k);
#         else:  # cusp forms only
#             M = CuspForms(eps,k);
#         # end if;
#         b = M.q_expansion_basis(prec=qprec)
#         C = eps(1).parent() # // cyclotomic field
#         Cq = PowerSeriesRing(C,default_prec=qprec)
#         return([Cq!f: f in b])
#     # end if;
    
#     # Now consider case when space not over Q

#      # Construct q-expansions of cuspidal space using modular symbols
#     MM = ModularSymbols(eps, k)
#     SS = CuspidalSubspace(MM)
#     SSb = SS.q_expansion_basis(qprec);
    
#     # Directly compute eisenstein series
#     M = ModularForms(eps,k)
#     if M.dimension() == 0:
#         return [];
#     # end if;
#     sturm = M.sturm_bound()
#     assert qprec > sturm; # // if this fails qprec is set too low
    
#     # // Put coefficients in a matrix
#     if is_cuspidal == false:
#         Es = M.eisenstein_series() // this space is defined over Q(eps)
#         dim = len(SSb) + len(Es)
#     else:  # cusp forms only
#         dim = len(SSb)
#     # end if;
#     # C = BaseRing(Parent(eps)); // cyclotomic field
#     C = eps(1).parent() 
#     A = zero_matrix(C,dim,sturm)
#     for i in [1..len(SSb)]:
#         for j in [1 .. sturm]:
#             A[i,j] = SSb[i][j-1]
#         # end for;
#     # end for;
#     if is_cuspidal == false: # // full space
#         for i = [len(SSb) + 1 .. dim]:
#             for j in [1 .. sturm]:
#                 A[i,j] = Es[i - len(SSb)][j-1]
#             # end for;
#         # end for;
#     # end if;
    
#     # Create pseudomatrices for A and R^dim      
#     R:=Integers(C);
#     I1:=1*R;
#     ps_A:=PseudoMatrix([I1: i in [1 .. dim]],A);
#     ps_R:=PseudoMatrix(Module(R,sturm)); // R^sturm as pseudo matrix

#     ps_Asat:=ps_A meet ps_R; // compute the intersection of the spaces
#     // I believe one can also just saturate ps_A directly.

#     assert ClassNumber(C) eq 1; // to ensure all ideals principal
#     // Steve says you should take second element in TwoElement(Is[i]) otherwise, and
#     // this will be a local generator.
    
#     Is:=CoefficientIdeals(ps_Asat);
#     Is_gen:=[];
#     for i:=1 to dim do // find generators for principal ideals
#         _,gen:=IsPrincipal(Is[i]);
#         Append(~Is_gen,gen);
#     end for;
#     Asat_vecs:=Matrix(ps_Asat);
    
#     Asat:=ZeroMatrix(C,dim,sturm);
#     for i:=1 to dim do
#         for j:=1 to sturm do
#             Asat[i,j]:=Is_gen[i]*Asat_vecs[i,j];
#         end for;
#     end for;
    
#     // Find transformation matrix from A to Asat
#     B:=Solution(A,Asat); // B*A eq Asat
    
#     // Transform basis elements to full q-adic precision
#     b:=[];
#     Cq:=PowerSeriesRing(C,qprec);
#     if is_cuspidal eq false then // full space
#         Mb:=[f: f in SSb] cat [PowerSeries(e,qprec): e in Es];
#     else // cusp forms only
#         Mb:=[f: f in SSb];
#     end if; 
#     Mb_sat:=[**];
#     for i:=1 to dim do
#         f:=Cq!0;
#         for j:=1 to dim do
#             f:=f + B[i,j]*(Cq!Mb[j]);
#         end for;
#         Append(~Mb_sat,f);
#     end for;
    
#     return Mb_sat;
    
# # end function;
                                       
# ClassicalBases.m lines 735-762
def ComplementarySpacesModp_B(p, k0, n, elldash, LWBModp, weightbound, characters, char):
    Z_N = Integers(char.modulus())
    # U_N,m_N:=UnitGroup(Z_N);
    # gens_N:=[ZZ!m_N(u): u in Generators(U_N)];
    gens_N = Z_N.unit_gens()
    char_on_gens = [char(u) for u in gens_N]
    char_on_gens_1 = [1 for u in gens_N]  # // trivial character

    CompSpacesCode = []

    ell = dimension_modular_forms(char, k0+n*(p-1))
    # OldBasisModp:=sub< VectorSpace(GF(p),elldash) | >; // Steve 19/02/13
    OldBasisModp = VectorSpace(GF(p), elldash)

    print("n =", n)
    totalnumberoftries = 0
    for i in range(0, n+1):
        print("Computing basis for M_(k0 + i(p-1)) mod p for k0+i(p-1),i,n:",
              k0+i*(p-1), i, n)
        #     // Here we update both LWBModp and characters, appending the new q-expansions found.
        TotalBasisModp, NewBasisCodemi, nLWBModp, characters, numberoftriesi = RandomNewBasisModp_B(p, k0, i, LWBModp, OldBasisModp, weightbound, characters, char, char_on_gens, char_on_gens_1)
        CompSpacesCode.append(NewBasisCodemi)
        OldBasisModp = TotalBasisModp  # // basis for M_(k0 + i(p-1))
        totalnumberoftries += numberoftriesi

    print("ell, totalnumberoftries: ", ell, totalnumberoftries)

    return(CompSpacesCode)

# ClassicalBases.m lines 673-730
# function RandomNewBasisModp_B(p,k0,j,LWBModp,OldBasisModp,weightbound,characters,char,char_on_gens,char_on_gens_1)


def RandomNewBasisModp_B(p, k0, j, LWBModp, OldBasisModp, weightbound, characters, char, char_on_gens, char_on_gens_1):
    """
    """

    k = k0 + j*(p-1)

    R = LWBModp[2][1].parent()  # // this should be non-empty, since it
    # is the space of weight two forms with
    # trivial character.

    #     // Construct TotalBasisModp
    TotalBasisModp = OldBasisModp  # // Recall E_(p-1) = 1 mod p.
    V = VectorSpace(GF(p), OldBasisModp.dimension())
    elldash = TotalBasisModp.dimension()
    # Steve 19/02/13: more efficient to use vector spaces than matrices.
    # Case k0 + i(p-1) = 0 + 0(p-1) = 0
    if k == 0 and char.order() == 1:
        TotalBasisModp = V.subspace([V[1]])
        # // March 28: added zero for number of tries.
        return([TotalBasisModp, [[]], 0])

    elif k == 0:  # non-trivial character in weight 0.
        # here [] should correspond to "nothing".
        return([TotalBasisModp, [], 0])

    # // Case k = k0 + i(p-1) > 0
    di = DimensionMNk(char, k)
    diminus1 = DimensionMNk(char, k-(p-1))
    mi = di - diminus1

    NewBasisCode = []
    rk = diminus1
    # // March 24, 2012: just a counter to see how long things are taking.
    numberoftries = 0
    for i in range(1, mi+1):  # // extra forms
        while rk < diminus1 + i:
            numberoftries += 1
            sol = RandomSolutionWithChar_B(
                characters, char_on_gens, char_on_gens_1, k0, j, p, weightbound)  # // 19/02/13
            TotalBasisi = R(1)
            TotalBasisiCode = sol
            for s in sol:
                TotalBasisi *= LWBModp[s[0]][s[1]]

            Vec = V([TotalBasisi[j] for j in [0..elldash-1]])
            if not Vec in TotalBasisModp:
                TotalBasisModp += V.subspace(Vec)
                rk = TotalBasisModp.dimension()
                NewBasisCode.append(TotalBasisiCode)

        if j > 0:  # // this is case k0 + j(p-1) with j > 0
            # // add this new q-expansion to LWBModp
            LWBModp[weightbound + j-1].append(TotalBasisi)
            # // add character "char" to list of characters
            characters[weightbound+j-1].append(char_on_gens)
    return([TotalBasisModp, NewBasisCode, LWBModp, characters, numberoftries])


def DimensionMNk(char, k):
    """
    Returns dimension of space of classical modular forms for char of weight k.	
    """
    if k > 0:
        M = ModularForms(char, k)         # space over k, all conjugates(?)
        deg = char(1).parent().degree()  # number of conjugates
        return(ZZ(M.dimension()/deg))
    elif k == 0 and char.order() == 1:
        return(1)               # trivial character
    else:                       # k > 0 or k=0 with non-triv char
        return(0)


# ClassicalBases.m lines 613-668
def RandomSolutionWithChar_B(characters, char_on_gens, char_on_gens_1,
                             k0, j, p, weightbound):
    """The output is a list [[a1,b1],....] where [ai,bi] is the bi
    character of weight ai, and the product of the associated modular
    forms is a weight "weight" and character "char".  19/02/13: use new
    approach to handling characters where you store them as a list of
    their values on generators - this reduces MASSIVE overhead on
    multiplying characters.  In this new version we append the previous
    spaces constructed to the LWB and use them. This is quite
    complicated.

    """

    weight = k0 + j*(p-1)  # 08/04/03 - this is the target weight.
    # char_on_gens: this is just the characters "char" values on the
    # generators of Z/N^* where N is the modulus.
    # char_on_gens_1: this is just the value of the trivial character on these
    # generators.
    B = weightbound + (j-1)  # start choosing from this position.
    found = false
    while found == false:
        K = weight
        sol = []
        charprod = char_on_gens_1  # // just all ones vector of right length.
        a = []
        # Choose elements in weights B,...,2.
        # for i:=B to 2 by -1 do
        for i in range(B-1, 0, -1):
            # print(i)
            if len(characters[i]) > 0:  # Nov 15: i.e. there are forms in that weight
                if i <= weightbound:
                    ai = randint(0, floor(K/i))  # pick ai elements of weight i
                    K = K - ai*i
                else:
                    # so form in position i is weight k0 + j_i(p-1)
                    j_i = i - weightbound
                    ai = randint(0, floor(K/(k0 + j_i*(p-1))))
                    K = K - ai*(k0 + j_i*(p-1))
                # end if;
                for m in [0..ai-1]:
                    # characters[i] = chars for weight i
                    bim = randint(0, len(characters[i-1]))
                    charprod = [charprod[l]*characters[i][bim][l] for
                                l in [0..len(charprod)-1]]  # charprod is 1 or 2.
                    sol.append([i, bim])
            else:
                ai = 0
            # end if;
        # end for;
        # Feb 18, 2013: some code which will work even when nothing in weight one
        if len(characters[1]) > 0:  # pick K elements in weight one
            for m in [0..K-1]:  # 08/04/13 j -> m
                b1m = randint(0, len(characters[0])-1)
                charprod = [charprod[l]*characters[0][b1m][l] for l in
                            [0..len(charprod)-1]]  # Feb 19, 2013
                sol.append([1, b1m])
                # end for;
            if charprod == char_on_gens:
                found = true
                # end if;
        else:  # nothing in weight one
            if K == 0 and charprod == char_on_gens:
                found = true
                #     end if;
                # end if;
                # end while;

    sol.reverse()
    return(sol)
# end function;

# ClassicalBases.m lines 536-551


def LowWeightBasesModp(LWB, p, elldash):
    R = PowerSeriesRing(GF(p), 'z', default_prec=elldash)
    LWBModp = []
    for i in [0..len(LWB)-1]:  # k = i
        LWBModpWeightk = []
        for j in [0..len(LWB[i])-1]:
            LWBModpWeightk.append(R(LWB[i][j]))
            # end for;
        LWBModp.append(LWBModpWeightk)
            # end for;
    return LWBModp

# end function;

# ClassicalBases.m lines 84-98


def qexpToPSinZpm(f, NN, degC, zeta_pm):
    # print(f)
    # print(f.base_ring()) # for debugging
    Zpm = zeta_pm.parent()
    # print(Zpm)
    R.<q> = PowerSeriesRing(Zpm, default_prec=NN)
    fPS = R(0)

    for i in [0..NN-1]:
        fi = f[i]
        # print("fi = ", fi)
        # print()
        fi_pm = sum(Zpm(fi[j])*zeta_pm^j for j in [0..degC-1])
        fPS += fi_pm*q**i
    # end for

    return fPS

# end function;
