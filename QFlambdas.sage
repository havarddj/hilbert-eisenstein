attach("~/Projects/HilbertEisenstein/cython_test/RM-set.spyx")
attach("~/Projects/HilbertEisenstein/worksheet.sage")

def p_adic_L(F, p, k0, chi, m, verbose=False):
    """Computes p-adic L-function attached to chi, a ring class
    character of a real quadratic field, by decomposing chi as a sum
    of (differences of) indicator functions and using the reduction
    theory of quadratic forms to compute L-functions of these.

    """
    d = F.degree()

    if p == 2:
        q = 4
        deltam = m
        kj = [k0+j*2 for j in [0..deltam]]
    else:
        q = p
        deltam = ceil(m*(p-1)/(p-2))
        kj = [k0+j*(p-1) for j in [0..deltam]]

    R.<Z> = PowerSeriesRing(QQ)
    Zpp = Zp(p, prec=deltam+1, type='capped-rel', print_mode='val-unit')
    Qpp = Zpp.fraction_field()
    Ls = [0 for j in [0..deltam]]  # initialise L-values
    
    for aC, C in decompose_chi(chi):
        Q = QF_of_class(C)
        D = Q.discriminant()

        D0 = fundamental_discriminant(D)
        F = QuadraticField(D0)
        f = sqrt(D/D0)

        if kronecker_symbol(D0, p) == 1:
            print(f"p = {p} is split in Q(sqrt({D0}))")
        elif kronecker_symbol(D0, p) == 0:
            print(f"p = {p} is ramified in Q(sqrt({D0}))")
        else:
            print(f"p = {p} is inert in Q(sqrt({D0}))")



        Mtop = ModularForms(Gamma1(f), weight=d*kj[deltam])
        S = Mtop.sturm_bound()
        # K = Psi.base_ring()
        K = Mtop.base_ring()
        R.<Z> = PowerSeriesRing(K)

        # Mdkj_dims = [ModularForms(Psi, weight=(d*kj[j])).dimension() for j in [0..deltam]]
        if verbose:
            print("S =", S)
            print("f =", f)
            print("Weights:", [d*k for k in kj])
            print("delta_m =", deltam)
            print("Base field:", K)
            # print("dims of M_{k_j}:", Mdkj_dims, "\n")
            print("Computing spaces of modular forms: \n")
            print("Computing low weight basis.")

        LWB = [ModularForms(Gamma1(f), weight=l,
                            prec=S).q_expansion_basis(prec=S) for l in range(1,7)]

        Mdkj = []
        if verbose:
            print("Computing higher weight bases:")

        for j in [0..deltam]:
            if Mdkj == []:
                Mdkj.append(Gamma1_high_weight_basis(f, d*kj[j], prec=S,
                                                             verbose=False, LWB=LWB))
            else:
                Mdkj.append(Gamma1_high_weight_basis(f, k=d*kj[j], prec=S,
                                                     prev=[d*kj[j-1],
                                                           Mdkj[j-1]],
                                                     verbose=False, LWB =
                                                     LWB))

            if verbose:
                print("Done with j =",j)
        if verbose:
            print("Done computing spaces of modular forms")
            # print("M_{top} has basis", Mdkj[deltam-1])
            # print("M_k0 has basis", Mdkj[0])
        QF_data = Delta_data(Q, S)  # as defined in RM-set.spyx
        # This will eventually be a list of non-const coeffs of Delta_k_j, j=0,..deltam

        # for j in [0..deltam]:
        #     for n in [1..S-1]:
        #         coeff_n = sum(A[0]^(kj[j]-1) for A in QF_data[n-1])
        #         Deltaj[j] +=4*Z^n*coeff_n

        #             # RMforms[n-1][0][0] = 1st coeff of Q' defined above
        #             # (index starts at 0)
        #     if verbose:
        #         print("Done computing diagonal restrictions")
        Deltaj = diag_restr_coeffs(Q, kj, S)
        if verbose:
            print("Deltaj", Deltaj)
        if Deltaj != [0 for j in [0..deltam]]:
            LsQ = [find_const_term(Mdkj[j], R(Deltaj[j]), field=K) for j in [0..deltam]]
        else:
            LsQ = [0 for j in [0..deltam]]

        if verbose:
            print(f"Classical L-values for {Q}:", LsQ)
        for j in [0..deltam]:
            Ls[j] += aC *LsQ[j]
    
    # Lsp = [Qpp(Ls[j]*Euler_factor(F, p, chi, kj[j])) for j in [0..deltam]]
    Lsp = [Qpp(Ls[j]*Euler_factor2(Euler_data(p, D, F,f), kj[j])) for j in [0..deltam]]
    points = [[1-kj[j], Lsp[j]] for j in [0..deltam]]
    P = Newton_poln(points)
    if verbose:
        print("P(s) = ", P, "\n \n")

    PolnQp.<T> = PolynomialRing(Zpp.fraction_field())
    T_coords = [(1+p)^(k-1)-1 for k in kj]
    if verbose:
        print("x-coordinates of interpolation points in T:", T_coords, "\n \n")
    pointsnew = [(T_coords[j], Lsp[j]) for j in [0..deltam]]

    Q = PolnQp(Newton_poln(pointsnew))
    # print(points)

    if verbose:
        print("Q(T) = ", Q)
    return(P,Q)
    

    
def decompose_chi(chi):
    """Given chi, returns pairs [a_C, C] such that chi = sum_C a_C
    (1_C - 1_{C conj})"""
    G = chi.parent().group()
    L = []
    for C in G.gens():
        L.append([chi(C),C])    # seems shady
    return(L)

def QF_of_class(C):
    I = C.ideal()
    # return(I.quadratic_form().reduced_form()) # should probably not reduce it?
    return(I.quadratic_form())
    
def find_lambda(Q):
    i = 0
    p = Q.base_ring().uniformiser()
    for coeff in Q.padded_list():
        if coeff % p != 0:
            return i
        i += 1
    return 0

def find_chars(F,N, rat=False):
    """Give a collection of totally odd ray class characters of F of
    modulus bounded above in norm by N and of order 2.

    """ 
    ideals = [I[0] for I in list(F.ideals_of_bdd_norm(N).values()) if I != []]
    if rat:
        ideals = [I for I in ideals if I in QQ]
    r1 = F.signature()[0]       # nr of real embeddings = 2
    char_list = []
    for I in ideals:
        mfrak = F.modulus(F.ideal(I), range(r1))  # narrow class group
        H = HeckeCharacterGroup(mfrak)
        for chi in H:
            if is_totally_odd(chi) and chi.order() == 2 and chi.is_primitive():
                char_list.append(chi)
    return char_list

def padicL_QF(Q, p, k0, m, verbose=False):
    D = Q.discriminant()

    D0 = fundamental_discriminant(D)
    # F = QuadraticField(D0)
    F = QuadraticField(D)
    f = sqrt(D/D0)
    N = f
    d = F.degree()

    if p == 2:
        q = 4
        deltam = m
        kj = [k0+j*2 for j in [0..deltam]]
    else:
        q = p
        deltam = ceil(m*(p-1)/(p-2))
        kj = [k0+j*(p-1) for j in [0..deltam]]

    R.<Z> = PowerSeriesRing(QQ)
    Zpp = Zp(p, prec=deltam+1, type='capped-rel', print_mode='val-unit')
    Qpp = Zpp.fraction_field()
    if kronecker_symbol(D0, p) == 1:
        print(f"p = {p} is split in Q(sqrt({D0}))")
    elif kronecker_symbol(D0, p) == 0:
        print(f"p = {p} is ramified in Q(sqrt({D0}))")
    else:
        print(f"p = {p} is inert in Q(sqrt({D0}))")



    Mtop = ModularForms(Gamma1(f), weight=d*kj[deltam])
    S = Mtop.sturm_bound()
    # K = Psi.base_ring()
    K = Mtop.base_ring()
    R.<Z> = PowerSeriesRing(K)

    # Mdkj_dims = [ModularForms(Psi, weight=(d*kj[j])).dimension() for j in [0..deltam]]
    if verbose:
        print("S =", S)
        print("f =", f)
        print("Weights:", [d*k for k in kj])
        print("delta_m =", deltam)
        print("Base field:", K)
        # print("dims of M_{k_j}:", Mdkj_dims, "\n")
        print("Computing spaces of modular forms: \n")
        print("Computing low weight basis.")

    LWB = [ModularForms(Gamma1(f), weight=l,
                        prec=S).q_expansion_basis(prec=S) for l in range(1,7)]

    Mdkj = []
    if verbose:
        print("Computing higher weight bases:")

    for j in [0..deltam]:
        if Mdkj == []:
            Mdkj.append(Gamma1_high_weight_basis(f, d*kj[j], prec=S,
                                                         verbose=False, LWB=LWB))
        else:
            Mdkj.append(Gamma1_high_weight_basis(f, k=d*kj[j], prec=S,
                                                 prev=[d*kj[j-1],
                                                       Mdkj[j-1]],
                                                 verbose=False, LWB =
                                                 LWB))

        if verbose:
            print("Done with j =",j)
    if verbose:
        print("Done computing spaces of modular forms")
        # print("M_{top} has basis", Mdkj[deltam-1])
        # print("M_k0 has basis", Mdkj[0])
    # QF_data = Delta_data(Q, S)  # as defined in RM-set.spyx
    # This will eventually be a list of non-const coeffs of Delta_k_j, j=0,..deltam

    # for j in [0..deltam]:
    #     for n in [1..S-1]:
    #         coeff_n = sum(A[0]^(kj[j]-1) for A in QF_data[n-1])
    #         Deltaj[j] +=4*Z^n*coeff_n

    #             # RMforms[n-1][0][0] = 1st coeff of Q' defined above
    #             # (index starts at 0)
    #     if verbose:
    #         print("Done computing diagonal restrictions")
    Deltaj = diag_restr_coeffs(Q, kj, S)
    if verbose:
        print("Deltaj", Deltaj)
    if Deltaj != [0 for j in [0..deltam]]:
        LsQ = [find_const_term(Mdkj[j], R(Deltaj[j]), field=K) for j in [0..deltam]]
    else:
        LsQ = [0 for j in [0..deltam]]

    if verbose:
        print(f"Classical L-values for {Q}:", LsQ)
        print("Euler data:", Euler_data(p, D, F,N))
    Lsp = [Qpp(LsQ[j]*Euler_factor2(Euler_data(p, D, F,N),kj[j])) for j in [0..deltam]]

    points = [[1-kj[j], Lsp[j]] for j in [0..deltam]]
    P = Newton_poln(points)
    if verbose:
        print("P(s) = ", P, "\n \n")

    PolnQp.<T> = PolynomialRing(Zpp.fraction_field())
    T_coords = [(1+p)^(k-1)-1 for k in kj]
    if verbose:
        print("x-coordinates of interpolation points in T:", T_coords, "\n \n")
    pointsnew = [(T_coords[j], Lsp[j]) for j in [0..deltam]]

    Q = PolnQp(Newton_poln(pointsnew))
    # print(points)

    if verbose:
        print("Q(T) = ", Q)
    return(P,Q)

def Euler_data(p, D, F, f=1):
    # print(f"p = {p}, D = {D}, F = {F}, f = {f}")
    kron = kronecker(D,p)
    if gcd(f,p) > 1:
        return [[1,p^2]]
    else:
        if kron == -1:
            return [[1, p^2]];
        else:
            pp = F.ideal(p).factor()[0][0]
            boo = pp.is_principal()
            assert boo == true
            g = pp.gens_reduced()[0]

            # print("g = ", g)
            # print("pp = ", pp)
            # print("boo", boo)
            sgn = 1
            if (not g.is_totally_positive()) and (not (-g).is_totally_positive()):
                sgn = -1
                
            if kron == 1:
                return [[sgn,p],[sgn,p]]
            else:
                return [[sgn,p]]

    
def Euler_factor2(data, n):
    return prod(1-sgn*p^n for sgn,p in data)

# table([(l.base_ring().uniformiser(),find_lambda(l) for r, l in L])
# if __name__ == "__main__":
    # # Q = BinaryQF_reduced_representatives(21)[1]
    # r1 = F.signature()[0]       # nr of real embeddings = 2
    # mfrak = F.modulus(F.ideal(4), range(r1))  # narrow class group
    # H = HeckeCharacterGroup(mfrak)
    # try:
    #     chi = H.gens()[1]
    #     # psi = H.one()
    #     # psi = psi.primitive_character()
    # except:
    #     chi = H.one()
    # # F.dirichlet_group()[0]
    # print("F = ", F,
    #       "\n H = ", H.group(),
    #       "\np = ", p,
    #       "\nk0 = ", k0,
    #       "\nm = ", m,
    #       "\npsi = ", psi)
    # N = 20
    # chi_list = find_chars(F, N, rat=True)
