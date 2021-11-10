attach("~/Projects/HilbertEisenstein/worksheet.sage")

def find_chars(F,N):
    """Give a collection of totally odd ray class characters of F of
    modulus bounded above in norm by N and of order 2.

    """
    ideals = [I[0] for I in list(F.ideals_of_bdd_norm(N).values()) if I != []]

    r1 = F.signature()[0]       # nr of real embeddings = 2
    char_list = []
    for I in ideals:
        mfrak = F.modulus(F.ideal(I), range(r1))  # narrow class group
        H = HeckeCharacterGroup(mfrak)
        for chi in H:
            if is_totally_odd(chi) and chi.order() == 2 and chi.is_primitive():
                char_list.append(chi)
    return char_list


def find_lambda(Q):
    i = 0
    for coeff in Q.padded_list():
        if not coeff % 3 == 0:
            return i
        i += 1



def batch_L_values(F, k0, chars, m, verbose=False, rand_bases=True, nebentype=False):
    """
    Similar to "classical_L_values" but computing for a bunch of characters at once
    """
    if isinstance(psi,sage.quadratic_forms.binary_qf.BinaryQF):
        RM = True
        Q = psi
        D = Q.discriminant()
        M = sqrt(D/D.squarefree_part())
        Psi = "bqf,no psi"
    else:
    # elif isinstance(psi,sage.modular.hecke_character.HeckeCharacterGroup_class_with_category.element_class):
        RM = False
        for psi in chars:
            mfrak = psi.modulus()
            mfraks.append(mfrak)
            M = mfrak.finite_part().smallest_integer()
            Ms.append(M)
            mfins.append(mfrak.finite_part())
            Psi = hecke_to_dirichlet_char(psi, M)
            Psis.append(Psi)
    # else:
    #     raise ValueError("Psi has to be a character or a binary quadratic form")
    d = F.degree()    
    assert d == 2, "only quadratic fields supported atm"
    # assert is_tot_even(psi) or is_tot_odd(psi), f"psi = {psi} needs to be totally even or totally odd."
    # assert (k0 % 2 == 0 and  is_tot_even(psi)) or (k0 % 2 == 1 and is_tot_odd(psi)), "k0 has wrong parity"

    if p == 2:
        q = 4
        deltam = m
        kj = [k0+j*2 for j in [0..deltam]]
    else:
        q = p
        deltam = ceil(m*(p-1)/(p-2))
        kj = [k0+j*(p-1) for j in [0..deltam]]



    if nebentype:
        for psi in chars:
        Mtops.append(ModularForms(Psi ,weight=d*kj[deltam]))
    else:
        Mtop = ModularForms(Gamma1(M), weight=d*kj[deltam])
    Ss = [Mtop.sturm_bound() for Mtop in Mtops]
    # K = Psi.base_ring()
    Ks = [Mtop.base_ring() for Mtop in Mtops]


    # Mdkj_dims = [ModularForms(Psi, weight=(d*kj[j])).dimension() for j in [0..deltam]]
    if verbose:
        print("Ss =", Ss)
        print("Ms =", Ms)
        print("Weights:", [d*k for k in kj])
        print("delta_m =", deltam)
        if chars:
            print("Chars =", chars)
        print("Base fields:", Ks)
        # print("dims of M_{k_j}:", Mdkj_dims, "\n")
        print("Computing spaces of modular forms: \n")
        print("Computing low weight basis.")

    if nebentype:
        for r in len(chars):
            LWBs.append([[ModularForms(Gamma0(Ms[r]), weight=l,
                                 prec=Ss[r]).q_expansion_basis(prec=Ss[r]) for l in range(1,7)],
               [ModularForms(chars[r], weight=l,prec=Ss[r]).q_expansion_basis(prec=Ss[r]) for l in range(1,7)]])
    else:
        LWB.append([ModularForms(Gamma1(Ms[r]), weight=l,
                            prec=Ss[r]).q_expansion_basis(prec=Ss[r]) for l in range(1,7)])
    Mdkj = []
    if verbose:
        print("Computing higher weight bases:")
    for j in [0..deltam]:
        if Mdkj == []:
            if rand_bases:
                if nebentype:
                    Mdkj.append(MF_char_high_weight_basis(Psi, d*kj[j], prec=S,
                                                          verbose=False, LWBs = LWBs))
                else:
                    Mdkj.append(Gamma1_high_weight_basis(M, d*kj[j], prec=S,
                                                         verbose=False, LWB = LWB))
            else:
                Mdkj.append(ModularForms(Psi,d*kj[j],prec=S).q_expansion_basis(S))
        else:
            if rand_bases:
                if nebentype:
                    Mdkj.append(MF_char_high_weight_basis(Psi, k=d*kj[j], prec=S,
                                                          prev=[d*kj[j-1], Mdkj[j-1]],
                                                          verbose=False, LWBs = LWBs))
                else:
                    Mdkj.append(Gamma1_high_weight_basis(M, k=d*kj[j], prec=S,
                                                         prev=[d*kj[j-1], Mdkj[j-1]],
                                                         verbose=False, LWB = LWB))
            else:
                Mdkj.append(ModularForms(Psi,d*kj[j],prec=S).q_expansion_basis(S))
        if verbose:
            print("Done with j =",j)
    if verbose:
        print("Done computing spaces of modular forms")
        print("M_{top} has basis", Mdkj[deltam-1])
        print("M_k0 has basis", Mdkj[0])

    ##### Compute diagonal restrictions:
    # if the modulus is an integer, use bqfs
    R.<Z> = PowerSeriesRing(K)
    if verbose:
        print("\nComputing diagonal restriction data:")
    if RM:
        if verbose:
            print("Computing QF data:")
        
        Deltaj = [0 for j in [0..deltam]]
        QF_data = Delta_data(Q, S)  # as defined in RM-set.spyx
        # This will eventually be a list of non-const coeffs of Delta_k_j, j=0,..deltam
        
        for j in [0..deltam]:
            for n in [1..S-1]:                    
                coeff_n = sum(A[0]^(kj[j]-1) for A in QF_data[n-1])
                Deltaj[j] +=4*Z^n*coeff_n

                # RMforms[n-1][0][0] = 1st coeff of Q' defined above
                # (index starts at 0)     
        if verbose:
            print("Done computing diagonal restrictions")
            # print("Deltaj", Deltaj)
    else:
        if verbose:
            print("Computing nu-sets")
        nu_set = []
        Diff = F.different()
        invDiff = Diff^-1
        for n in [1..S]:
            nu_set.append(tp_elts_of_trace_n(F, n, invDiff,floor(n*a)))
        if verbose:
            print("Computing Delta_j, IN PARALLEL!!")
            start_time = time.time()
        Deltaj = sorted(list(compute_Delta_j_chars([(j,F,chars,kj,S,Diff,nu_set) for
                                       j in [0..deltam]]))) # returns (input,output)

        Deltaj = [l[1] for l in Deltaj]
        if verbose:
            end_time = time.time()
            print(f'Computing delta_j took {end_time - start_time} seconds')
    if verbose:
        print("Higher coefficients of Delta_{k_0}:", Deltaj[0])
        # print("Higher coefficients of Delta_{k_1}:", Deltaj[1])
        # print("Higher coefficients of Delta_{k_deltam}:", Deltaj[deltam])
    # Using some linear algebra, we find the constant terms, which are
    # classical L-values; viz. const term(Delta_j) = L(psi,1-k_j)
    Ls = [find_const_term(Mdkj[j], R(Deltaj[j]), field=K) for j in [0..deltam]]
    
    if verbose:
        print("Classical L-values:", Ls)
    return(Ls)




# @parallel(ncpus=8)
# def compute_Delta_j_chars(j,F,chars,kj,S,Diff,nu_set,verbose=False):
#     """
#     Computes higher coeffs of diagonal restriction, hopefully in parallel

#     RMK: this can be improved by storing pfrak, psi(prak) in a dictionary
#     """
#     R.<Z> = PowerSeriesRing(QQ)

#     Deltas = [0 for psi in chars]
#     for n in [1..S-1]:
#         coeff_n = [0 for psi in chars]
#         for nu in nu_set[n-1]:
#             factors = (nu*Diff).factor()
#             cop_factors = [[pfrak, mult] for pfrak, mult in factors]
#             ell = len(cop_factors)
#             # pfraks = [pfrak for pfrak, mult in cop_factors]

#             psi_pfrak = [[psi(pfrak) for pfrak,_ in cop_factors] for psi in chars]
#             norm_pfrak = [pfrak.norm() for pfrak,_ in cop_factors]
#             for tup in Product(*[[0..m] for _, m in cop_factors]):
#                 if tup:
#                     for r in range(len(chars)):
#                         coeff_n[r] += prod(psi_pfrak[r][i]**tup[i]*(norm_pfrak[i]**tup[i])**(kj[j]-1) for i in range(ell))
#                 else:
#                     coeff_n += 1  # corr to ideal 1
#         for r in range(len(chars)):
#             Deltas[r] += 4*Z^n*coeff_n[r]
#     return(Deltas)



def diagonal_restriction_Lp_chars(F, p, k0, chars, m, verbose=False,Lprec=512):
    """Compute p-adic L-function L_p(psi,s) as polynomial in s using
    diagonal restriction.
    Currently only works for F real quadratic

    """
    assert p > 2, "p=2 is messy"
    if p == 2:
        q = 4
        deltam = m
        kj = [k0+j*2 for j in [0..deltam]]
    else:
        q = p
        deltam = floor(m*(p-1)/(p-2))
        kj = [k0+j*(p-1) for j in [0..deltam]]

    # print('deltam')

    assert [is_totally_odd(psi) for psi in chars] == [True for psi in chars], "all chars must be totally odd"  # cf LV21 --- ensures totally odd

    R.<Z> = PowerSeriesRing(QQ)
    Zpp = Zp(p, prec=deltam+1, type='capped-rel', print_mode='val-unit')
    Qpp = Zpp.fraction_field()

    # Ls = classical_L_values(F, k0, psi, m, verbose=verbose,
    #                         rand_bases=True, nebentype=nebentype)
    Plist = []
    Qlist = []
    # print("bump")
    for psi in chars:
        print("Modulus:", psi.modulus())
        print("log-vals on gens:", psi._log_values_on_gens())
        L = psi.Lfunction(Lprec)
        Ls = []

        for j in [0..deltam]:
            L_algdep = algdep(L(1-kj[j]), 1).roots()
            # print(L_algdep)
            if L_algdep == []:

                Ls.append(0)
            else:
                Ls.append(L_algdep[0][0])
                

        if verbose:
            print("Classical L-values:", Ls)
        Lsp = [Qpp(Ls[j]*Euler_factor(F, p, psi, kj[j])) for j in [0..deltam]]
        if verbose:
            print("p-adic L-values:", Lsp)
        # R.<Z> = PolynomialRing(QQ)
        # R = Zp(p, prec=m)
        # S.<s> = PolynomialRing(R)



        # Compute the Newton interpolation polynomial from the special
        # values (1-k,L_p(psi,1-k)):

        if verbose:
            print("x-coordinates of interpolation points in s:", [1-k for k in kj], "\n \n")
        points = [[1-kj[j], Lsp[j]] for j in [0..deltam]]
        P = Newton_poln(points)
        Plist.append([psi, P])
        if verbose:
            print("P(s) = ", P, "\n \n")
        K = psi.base_ring()            
        PolnQp.<T> = PolynomialRing(Zpp.fraction_field())
        T_coords = [(1+p)^(k-1)-1 for k in kj]
        if verbose:
            print("x-coordinates of interpolation points in T:", T_coords, "\n \n")
        pointsnew = [(T_coords[j], Lsp[j]) for j in [0..deltam]]

        Q = PolnQp(Newton_poln(pointsnew))
        Qlist.append([psi,Q])

        if verbose:
            print("Q(T) = ", Q)
    print("Q(T) = ", Q)
    lambda_list = []
    for psi, Q in Qlist:
        if Q != 0:
            lambda_list.append([F, psi, find_lambda(Q)])
    return(lambda_list)
