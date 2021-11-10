import time
import pickle
from itertools import product as Product
from multiprocessing import Process, Manager

attach("~/Projects/HilbertEisenstein/cython_test/RM-set.spyx")

def diagonal_restriction_Lp(F, p, k0, psi, m, verbose=False, rand_bases=True,nebentype=True):
    """Compute p-adic L-function L_p(psi,s) as polynomial in s using
    diagonal restriction.
    Currently only works for F real quadratic

    """

    if p == 2:
        q = 4
        deltam = m
        kj = [k0+j*2 for j in [0..deltam]]
    else:
        q = p
        deltam = floor(m*(p-1)/(p-2))
        kj = [k0+j*(p-1) for j in [0..deltam]]

    # print('deltam')
    mfrak = psi.modulus()
    M = mfrak.finite_part().smallest_integer()
    # assert is_totally_odd(psi)  # cf LV21 --- ensures totally odd
    R.<Z> = PowerSeriesRing(QQ)
    Zpp = Zp(p, prec=deltam+1, type='capped-rel', print_mode='val-unit')
    Qpp = Zpp.fraction_field()

    Ls = classical_L_values(F,p, k0, psi, m, verbose=verbose,
                            rand_bases=True, nebentype=nebentype)
    # L = psi.Lfunction(512)
    # Ls = [algdep(L(1-kj[j]), 1).roots()[0][0] for j in [0..deltam]]
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

def diagonal_restriction_Lp_multps(F, plist, k0, psi, m, verbose=False, rand_bases=True,nebentype=True):
    """Compute p-adic L-function L_p(psi,s) as polynomial in s using
    diagonal restriction.
    Currently only works for F real quadratic
NOT CURRENTLY FUNCTIONING
    """
    p = max(plist)
    if p == 2:
        q = 4
        deltam = m
        kj = [k0+j*2 for j in [0..deltam]]
    else:
        q = p
        deltam = floor(m*(p-1)/(p-2))
        kj = [k0+j*(p-1) for j in [0..deltam]]

    # print('deltam')
    mfrak = psi.modulus()
    M = mfrak.finite_part().smallest_integer()
    # assert is_totally_odd(psi)  # cf LV21 --- ensures totally odd
    R.<Z> = PowerSeriesRing(QQ)
    Ls = classical_L_values(F,p, k0, psi, m, verbose=verbose,
                            rand_bases=True, nebentype=nebentype)
    # L = psi.Lfunction(512)
    # Ls = [algdep(L(1-kj[j]), 1).roots()[0][0] for j in [0..deltam]]
    print("Classical L-values:", Ls)
    Lplist = []
    for p in plist:
        Zpp = Zp(p, prec=deltam+1, type='capped-rel', print_mode='val-unit')
        Qpp = Zpp.fraction_field()
        Lsp = [Qpp(Ls[j]*Euler_factor(F, p, psi, kj[j])) for j in [0..deltam]]
        if verbose:
            print("p-adic L-values:", Lsp)


            # Compute the Newton interpolation polynomial from the special
            # values (1-k,L_p(psi,1-k)):

        if verbose:
            print("x-coordinates of interpolation points in s:", [1-k for k in kj], "\n \n")
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
        Lp_list.append([P,Q])
    # print(points)

    # if verbose:
    #     print("Q(T) = ", Q)
    return(Lp_list)


def Euler_factor(F, p, psi, kj):
    r"""Return $\prod_{\mathfrak p | (p)}(1-\psi(\mathfrak p)
    \mathrm{Nrm}(\mathfrak p)^{k_j-1})$

    for quadratic fields only
    """                         
    if F.discriminant() % p == 0:
        return 1 - psi(p)*p^(kj-1)     # p ramified
    else:
        return prod(1-psi(pf[0])*pf[0].norm()^(kj-1) for pf in F.ideal(p).factor())



def hecke_to_dirichlet_char(psi, M):
    """Given Hecke character psi with modulus mfrak such that mfrak
    cap Z = (M), return associated Dirichlet character mod M, Psi

    """
    # assert M == psi.modulus().finite_part().smallest_integer()
    D = DirichletGroup(M)
    gens = D.unit_gens()
    return D([psi(g) for g in gens])



def find_const_term(M, f, field=QQ):
    """given a q-expansion basis M, and f a q-expansion in M with
    missing constant term, find constant term of f.
    """
    if M == []:
        print("Basis is empty")
        return(0)
    R.<Z> = PowerSeriesRing(field)
    A = [R(M[j]).padded_list()[1:] for j in [0..len(M)-1]]
    Mat = Matrix(field, A)
    try:
        soln = Mat.transpose() \ vector(field, f.padded_list()[1:])
    except ValueError as e:
        # print(e)
        if str(e) != "matrix equation has no solutions":
            raise
        else:
            print("not in space of forms")
            return(0)
    CT = (sum(soln[n]*M[n] for n in [0..len(M)-1]))[0]
    return(CT)

def Newton_poln(points):
    """Return Newton interpolation polynomial, computed using divided differences. Input is a list of tuples points = [(x_0,y_0),(x_1,y_1),...,(x_n,y_n)]. The resulting polynomial P is of degree n and satisfies P(x_i) = y_i for every tuple in points"""
    K = parent(points[0][1]).fraction_field()
    n = len(points)
    coefs = copy(zero_matrix(K, n))
    for i in range(n):
        coefs[i,0] = points[i][1]
    for j in [1..n]:
        # print("j=",j)
        for i in [0..n-j-1]:
            # print("i=", i)
            coefs[i, j] = (coefs[i+1, j-1] - coefs[i, j-1]) / (points[i+j][0]-points[i][0])
    R.<s> = PolynomialRing(K)
    P = sum(coefs[0][i]*prod((s-points[j][0]) for j in [0..i-1]) for i in [0..n-1])
    return P


# def forward_diff(i, L):
#     """Return Delta^i (f(x_0)) where L is a list of pairs (x_j,
#     f(x_j)) with x_j - x_(j-1) = h for all j"""
#     return(sum((-1)^(i-j)*binomial(i,j)*L[j][1] for j in [0..i]))
    # S.<s> = PolynomialRing(QQ)
    # h = points[1][0] - points[0][0]
    # assert h == (points[2][0] - points[1][0])
    # return(sum(forward_diff(i, points) * prod(s - points[0][0] -k for k in [0..i-1])/factorial(i) for i in [0..len(points)-1]))

# def diagonal_restriction_Lp(F, p, k0, psi, m):

#     # F.<a> = NumberField(x^3-3*x-1)

#     if not F.is_totally_real():
#         print("F is not totally real!")

#     d = F.degree()
#     if d > 3:
#         return("degree too high for now")

#     Psi = psi
#     M = m
#     # deltam = 24
#     deltam = ceil(m*(p-1)/(p-2))
    
#     kj = [2+j*(p-1) for j in [0..deltam]]  # note k0 = 2 as above

#     S = ModularForms(Psi,d*kj[deltam]).sturm_bound()
#     print(S)
#     # Mdkj = [ModularForms(Psi, d*kj[j], base_ring=QQ).q_expansion_basis()
#     #     for j in [0..deltam]]
#     # To compute q-expansions more efficiently, we should probably use
#     # Alan's algorithm

#     # return(compute_Delta_j(psi, S, deltam, d, kj, p))
#     R.<q> = PowerSeriesRing(QQ)
#     DeltaList = []
#     mfrak = psi.conductor()
#     QFlist = BinaryQF_reduced_representatives(mfrak)
#     print(QFlist)
#     for Q in QFlist:
#         RMlist = []
#         for n in [1..S-1]:
#             # print(n)
#             RMset = get_RM_set(n, Q)
#             # print(RMset)
#             if len(RMset) !=0:
#                 RMlist.append(RMset)
#         print("Computed RM sets \n")
#         for j in [0..deltam]:
#             Delta = 0
#             for n in [1..S-1]:
#                 for l in [0..len(RMlist)-1]:
#                     RM = RMlist[l]
#                     print(RM[0])
#                     if RM != []:
#                         Q2 = RM[0]
#                         tau = stable_root(Q2[0])
#                         print(tau)
#                         Delta += q^n *2^d*psi(tau)*Q2[0][0]^(kj[j]-1)
#         DeltaList.append(Delta)
#     print("Computed Delta Q = ", Q)

#     return(DeltaList)
#     # Clean this up into consecutive computations

# Using some pari functions to compute ray class groups
# Note; we don't need to use the gp wrapper

# # load("rayclass.py") 
from sage.libs.pari.all import pari, pari_gen
def conductor(O):
    f = O.index_in(O.number_field().maximal_order())
    return(f)

def order_narrow_class_group(O):
    # if O.is_maximal():
    #     return(O.narrow_class_group())
    P = pari(O.number_field().pari_polynomial())
    f = conductor(O)
    r1 = O.number_field().signature()[0]
    # # k = O.number_field().pari_bnf(flag=1)
    # bnf = pari.bnfinit(P,flag=1)
    # bnr = pari.bnrinit(bnf,Pari(conductor(O)),flag=1)
    # Cl = pari("bnf.")
    L = pari("bnrinit(bnfinit({}),[{},{}],1).cyc".format(P, f, [1 for n in [1..r1]]))
    # The list in the end says that the modulus contains f and all the infinite places
    # cycle_structure = tuple( ZZ(c) for c in k.bnf_get_cyc() )
    if len(L) == 0:
        return(AbelianGroup([1]))
    return(AbelianGroup(L))


def mod_narrow_class_group(F,m):
    """
    Takes number field F and integral ideal m, and returns Cl_m^+
    """
    P = F.pari_polynomial()
    r1 = F.signature()[0]
    # Initialise pari field
    bnf = pari("bnfinit({})".format(P))
    # For the time being, let's assume m is principal
    # decompose m in pari
    mPari = pari("idealhnf({},{})".format(bnf, m))
    # compute Cl_m^+
    L = pari("bnrinit({bnf},[{mPari},{L}],1).cyc".format(bnf=bnf,mPari=mPari, L=[1 for n in [1..r1]]))
    if len(L) == 0:
        return(AdditiveAbelianGroup([1]))
    return(AdditiveAbelianGroup(L))

def mod_narrow_ray_char_grp(F,m):
    """
    Takes number field F and integral ideal m, and returns the
    group of characters on Cl_m^+
    """
    P = F.pari_polynomial()
    bnf = pari("bnfinit({})".format(P))
    r1 = F.signature()[0]
    mPari = pari("idealhnf({},{})".format(bnf, m))
    bnr = pari("bnrinit({}, [{mPari},{L}],1)".format(bnf, mPari=mPari,
                                                      L=[1 for n in [1..r1]]))
    g = pari("[{}.gen[1]]".format(bnr))
    return(pari("bnrchar({},{})".format(bnr, g)))

def get_RM_set2(n, Q):
    """
    Compute the set of RM-points RM(n, tau)_f using algorithm from
    [lauder-vonk?]_

    Input:
    - $n$ an integer
    - a reduced binary quadratic form with stable root tau
    Output:
    - the set of augmented RM points of discriminant n^2 D, as a list
      of lists [F,gamma_n] where F is a quadratic form and gamma_n a
      determinant n matrix with integer coeffs

    Note that the norm of the ideal corresponding to F is precisely
    the first coefficient, a.
    """
    D = Q.discriminant()
    D0 = D.squarefree_part()
    f = sqrt(D/D0)
    # First, we do a stupid computation to find possible elements of
    # stabilisers of RMpts
    Stab = []
    for a in [-1..1]:
        for b in [-1..1]:
            for c in [-1..1]:
                for d in [-1..1]:
                    if a*d-b*c == 1:
                        s = Matrix(ZZ,[[a,b],[c,d]])
                        if s^(12) == Matrix(ZZ,[[1,0],[0,1]]):
                            Stab.append(s)

    Mn = []
    for d in divisors(n):
        if gcd(d, n/d) == 1:
            for j in [0..d-1]:
                Mn.append(Matrix(ZZ,[[n/d,j], [0,d]]))

    def is_in_SL2Z(Q, gamman, gnp, Stab):
        # """Returns True if $\gamma_n'\mathrm{Stab(\tau)}\gamma_n^-1
        # \subset \Sl_2(\Z), and false otherwise"""
        for s in Stab:
            if Q.matrix_action_right(s) == Q:
                Mat = gnp*s*(gamman^-1)
                for i in [0..1]:
                    for j in [0..1]:
                        if not Mat[i][j].is_integer():
                            return False
        return True

    RMn = []
    Qlist = []
    for k in range(len(Mn)):
        gamman = Mn[k]
        flag = False                            
        for [Qnp, gnp] in RMn:
        # for gnp in Mn[:k-1]:
            if is_in_SL2Z(Q, gamman, gnp, Stab):
                flag = True
        if not flag:
            Q2 = Q.matrix_action_right(gamman)
            Qlist.append(Q2)
            # print(Q2,gcd(f,Q2[0]),f)
            # print(Q2)
            # print("Forms =", Jans_RM_algo(Q2))
            RMn += [[Qn, gamman] for Qn in magma_algo(Q2) if gcd(Qn[0],f) == 1]
    print(Qlist, "of length", len(Qlist))
    return(RMn)

# def nearly_reduced_forms(Q2):
#     """Returns forms specified in Step 4. of [Lauder-Vonk2021],
#     where we additionally test for coprimality of a(Q2) with f"""
#     Q0 = Q2.reduced_form()
#     Q2 = Q0
#     forms = []
#     has_started = False
#     while Q0 != Q2 or not has_started:  # let it run at least once
#         has_started = True
#         # print(Q2,Q0, Q2 == Q0)

#         a, b, c = list(Q2)
#         # print(a,b,c)
#         d = floor(sqrt(b^2-4*a*c))         # sqrt of discriminant
#         cabs = c.abs()
#         if cabs >= d:
#             s = c.sign() * ((cabs + b) / (2 * cabs)).floor()
#         else:
#             s = c.sign() * ((d + b) / (2 * cabs)).floor()
#             # print(s)3
#         if a >= 0:
#             for i in range(1, s+1):
#                 a2 = a+ i*b + c*i**2
#                 Q2 = BinaryQF(a2, b + 2*i*c, c)
#                 # if gcd(a2, f) == 1:  # only include bois with first coeff coprime to  
#                 forms.append(Q2)
#         else:
#             for i in range(1,s+1):
#                 Q2 = BinaryQF(c, -b + 2*i*c, a - i*b + c*i**2)
#                 # if gcd(c, f) == 1:
#                 forms.append(Q2)
#         print("Is Q2 reduced?", Q2.is_reduced())
#         print("Q2 = ", Q2)
#         print("Q0 = ", Q0)
#         # for i in range(1,s+1):
#         #     Q2 = BinaryQF(c, -b + 2*i*c, a - i*b + c*i**2)
            
#         #     if gcd(c, f) == 1:
#         #         forms.append(Q2)
#         # print("Is Q2 reduced:", Q2.is_reduced())
#         # print("Disc Q2 =", Q2.discriminant())
#         # print("Disc Q0 =", Q0.discriminant())
#         # print("Q0,Q2 is equivalent:", Q0.is_equivalent(Q2))
#     return(forms)

def Jans_RM_algo(Q):
    a,b,c = list(Q)
    D = b^2-4*a*c
    d = floor(sqrt(b^2-4*a*c))         # sqrt of discriminant
    cabs=c.abs()
    if cabs >= d:
        s = c.sign() * ((cabs + b) / (2 * cabs)).floor()
    else:
        s = c.sign() * ((d + b) / (2 * cabs)).floor()
    forms = []
    if a > 0:
        for i in range(1, s+1):
            A = a + i*b+i**2 *c
            B = b + 2*i*c
            C = c
            Q = BinaryQF([A,B,C])
            # Q2 = BinaryQF([C,-B,A])
            Qprodroots = -(-B^2+ D)/(4*A^2)
            # Q2prodroots = -(-B^2+D)/(4*C^2)
            if Qprodroots < 0 and not Q in forms:
                forms.append(Q)
            # if Q2prodroots < 0 and not Q in forms:
            #     forms.append(Q2)

    if a < 0:
        for i in range(1, s+1):
            C = a - i*b+i**2 *c
            B = -b + 2*i*c
            A = c
            Q = BinaryQF([A,B,C])
            Q2 = BinaryQF([C,-B,A])
            Qprodroots = -(-B^2+ D)/(4*A^2)
            # Q2prodroots = -(-B^2+D)/(4*C^2)
            if Qprodroots < 0 and not Q in forms:
                forms.append(Q)
            # if Q2prodroots < 0 and not Q in forms:
            #     forms.append(Q2)

    # Q is now the form with i = s
    Q0 = Q
    done = False
    while not done:
        a, b, c = list(Q)
        # d = floor(sqrt(b^2-4*a*c))         # sqrt of discriminant
        # cabs = c.abs()
        # if cabs >= d:
        #     s = c.sign() * ((cabs + b) / (2 * cabs)).floor()
        # else:
        #     s = c.sign() * ((d + b) / (2 * cabs)).floor()
        # It seems plausible that s is computed once and for all, not at every step.
        if a > 0:
            # print("a pos")
            for i in range(1,s+1):
                A = a + i*b + c*i**2 
                B = b + 2*i*c
                C = c
                Q = BinaryQF([A,B,C])
                # Q2 = BinaryQF([C,-B,A])
                Qprodroots = -(-B^2+ D)/(4*A^2)
                # Q2prodroots = -(-B^2+D)/(4*C^2)
                if Qprodroots < 0 and not Q in forms:
                    forms.append(Q)
                # if Q2prodroots < 0 and not Q in forms:
                #     forms.append(Q2)

        if a < 0:
            # print("a neg")
            for i in range(1,s+1):
                C = a - i*b + c*i**2
                B = -b + 2*i*c
                A = c
                Q = BinaryQF([A,B,C])
                Qprodroots = -(-B^2+ D)/(4*A^2)
                # Q2prodroots = -(-B^2+D)/(4*C^2)
                if Qprodroots < 0 and not Q in forms:  # actually, first cond is automatic from AC<0
                    forms.append(Q)
                # if Q2prodroots < 0 and not Q in forms:
                #     forms.append(Q2)
        # if Q0 == Q:
        if Q in forms:
            done = True

    for Q in forms:
        A,B,C = list(Q)        # add Q2 = Q|(0, -1 ; 1,0)
        Q2 = BinaryQF([C,-B,A])
        Q2prodroots = -(-B^2+D)/(4*C^2)  # prod of roots of Q2(x,1)
        if Q2prodroots < 0 and not Q2 in forms:
            forms.append(Q2)
    return(forms)


# def unred_cycle(Q):             # rho(f) is necessarily reduced by Lemma 6.5.6
#     Qdash = Q._Rho()
#     C = []
#     while Qdash != Q:
#         C.append(Qdash)
#         Qdash = Qdash._Rho()
#     return C
        
def stable_root(Q):
    """ Returns stable root of indefinite quadratic form Q
    """
    if Q[0] == 0:
        return(0)
    return((-Q[1]+sqrt(Q.discriminant()))/(2*Q[0]))

# def compute_Delta_j(psi, S, deltam, d, kj, p):
#     """Computes non-constant terms of diagonal restriction
#     Input:

#     Output: list where the j-th entry is the power series Delta_j

#     """
#     R.<q> = PowerSeriesRing(QQ)
#     DeltaList = []
#     mfrak = psi.conductor()
#     QFlist = BinaryQF_reduced_representatives(mfrak)
#     print(QFlist)
#     for Q in QFlist:
#         RMlist = []
#         for n in [1..S-1]:
#             # print(n)
#             RMset = get_RM_set(n, Q)
#             # print(RMset)
#             if len(RMset) !=0:
#                 RMlist.append(RMset)
#         print("Computed RM sets \n")
#         for j in [0..deltam]:
#             Delta = 0
#             for n in [1..S-1]:
#                 for l in [0..len(RMlist)-1]:
#                     RM = RMlist[l]
#                     print(RM[0])
#                     if RM != []:
#                         Q2 = RM[0]
#                         print("\psi(0)", psi(F(1)))
#                         Delta += q^n *2^d*psi(F(stable_root(Q2[0])))*Q2[0][0]^(kj[j]-1)
#         DeltaList.append(Delta)
#     print("Computed Delta Q = ", Q)

#     return(DeltaList)
#     # Clean this up into consecutive computations

    # if RM:
    #     R.<Z> = PowerSeriesRing(QQ)
    #     mfrak = 1
    #     QFlist = BinaryQF_reduced_representatives(mfrak)
    #     print(QFlist)
    #     for Q in QFlist:
    #         RMlist = []
    #         for n in [1..S-1]:
    #         # print(n)
    #             RMset = get_RM_set(n, Q)
    #             # print(RMset)
    #             if len(RMset) !=0:
    #                 RMlist.append(RMset)
    #     print("Computed RM sets \n")
    #     for j in [0..deltam]:
    #         Delta = 0
    #         for n in [1..S-1]:
    #             for l in [0..len(RMlist)-1]:
    #                 RM = RMlist[l]
    #                 # print(RM[0])
    #                 if RM != []:
    #                     Q2 = RM[0]
    #                     Delta += Z^n *2^d*Q2[0][0]^(kj[j]-1)
    #                     Deltaj.append(Delta)
    #     # print("Computed Deltaj")
    # else:

def compute_Xn(F,n, mod):

    Diff = F.different()
    L = []
    for nu in tp_elts_of_trace_n(F,n,floor(n*a)):
        for A in divisors(nu*Diff):
            if A.is_coprime(mod):
                L.append([nu, A])
    return(L)

# @cached_function -- doesn't make sense in case of random errors
def Gamma1_high_weight_basis(M, k, prec, prev=[], verbose=False, LWB=False):
    """Compute basis for spaces of modular forms M_{k}(Gamma1(M))
    by multiplying lower weight basis elements.

    prev = [oldk, oldBasis = [g]]
    """
    R.<q> = PowerSeriesRing(QQ,default_prec=prec)
    dim = dimension_modular_forms(Gamma1(M), k)
    if verbose:
        print("Level:", M)
        print("Weight:", k)
        print("Dim:", dim)
    if k <= 6:
        return(ModularForms(Gamma1(M),weight=k).q_expansion_basis(prec=prec))
    if verbose:
        print("Computing low weight bases")
        start_time = time.time()
    # first compute lower weight bases for k =1,2,...,6:
    if not LWB:
        LWB = [ModularForms(Gamma1(M), weight=l, prec=prec).q_expansion_basis(prec=prec)
               for l in range(1,7)]
    if verbose:
        end_time = time.time()
        print(f'Computing LWB took {end_time - start_time} seconds')
    basis = []
    if prev != []:
        # NOTE: probably faster to compute image of everything and then refine to a basis
        oldk = prev[0]
        # print(oldk)
        oldBasis = prev[1]
        if verbose:
            print("Previous basis found")
            print("Computing middle weight basis:")
        MWB = Gamma1_high_weight_basis(M, k-oldk, prec, verbose=False)
        # we first compute all the product of every pair (g,h) with g
        # in oldBasis and h in MWB.

        old_prods = [((g*h).O(prec)).padded_list() for g in oldBasis for h in MWB]
        # next we row-reduce
        basis = [sum(l[n]*q^n for n in [0..prec-1]) + O(q^prec)
                 for l in Matrix(QQ, old_prods).echelon_form().rows() if not l.is_zero()]
        
        assert (len(basis) <= dim), "too many old forms!"
        
                
                    # if verbose: 
                    #     print("Appended:", f)
                    #     print("Now len(basis) = ", len(basis))
        

    if verbose:
        print("computing random products")
        start_time = time.time()
    while len(basis) < dim:
        P = Partitions(k, min_part=1, max_part=6, max_length=ceil(k/6+1)
                       ).random_element()
        # print(P)
        f = prod(sum(QQ.random_element()*g for g in LWB[wt-1]) for wt in P)
        f = R(f).lift_to_precision(prec)
        f = f.O(prec)
        # note: QQ.random_element() doesn't seem very random
        if f.padded_list() != [0 for l in range(0,f.prec())]:
            if not is_qexp_in_space(basis, f):
                basis.append(f)
            # if verbose: 
            #     print("Appended:", f)
            #     print("Now len(basis) = ", len(basis))
    if verbose:
        end_time = time.time()
        print(f'Computing random products took {end_time - start_time} seconds')
    # print("done")
    return(basis)

# actually, why randomise basis if we can pick random vector?
# def randomise_basis(M):
#     F = M[0].base_ring()        # this ought to be Q
#     prec = M[0].prec()
#     R.<q> = PowerSeriesRing(F,default_prec=prec)
#     M_list = [f.padded_list(prec) for f in M]
#     Mdim = len(M)
#     Mat = random_matrix(F, Mdim)
#     while Mat.determinant() == 0:  # ensure that it's invertible
#         Mat = random_matrix(F,Mdim)        
#     return
# def f_in_basis(f, basis):
#     F = f.base_ring()
#     if H[0].prec() != prec:
#         raise ValueError("precision of f (= {}) must equal "
#                          " precision of H (= {})".format(prec, H[0].prec()))

#     MS = MatrixSpace(QQ, prec, len(H))
#     A = copy(MS.zero_matrix())
#     for i in range(len(H)):
#         for j in range(prec):
#             A[j, i] = H[i][j]
#     B = vector(f)
#     B = vector(B.list() + [0, 0])
#     try:
#         return(A.solve_right(B))
#     except ValueError:
#         return False 

def is_qexp_in_space(H, f):
    if H == []:
        return False
    prec = f.prec()
    if H[0].prec() != prec:
        raise ValueError("precision of f (= {}) must equal "
                         " precision of H (= {})".format(prec, H[0].prec()))

    F = H[0].base_ring()
    MS = MatrixSpace(F, prec, len(H))
    A = copy(MS.zero_matrix())
    for i in range(len(H)):
        for j in range(prec):
            A[j, i] = H[i][j]
    # print(A.nrows(),A.ncols())
    B = vector(f.padded_list(prec))
    # B = vector(B.list() + [0, 0])
    try:
        return(A.solve_right(B))
    except ValueError as e:
        if str(e) != "matrix equation has no solutions":
            raise
        else:
            return False


# @parallel(ncpus=8)
def compute_Delta_j(j,F,psi,kj,S,Diff,nu_set,mfin,return_dict,verbose=False):
    """
    Computes higher coeffs of diagonal restriction, hopefully in parallel

    RMK: this can be improved by storing pfrak, psi(prak) in a dictionary
    """
    mess = "all good in the hood"
    print(mess)
    R.<Z> = PowerSeriesRing(QQ)
    # if verbose:
    #     print("Computing diagonal restrictions:")
    # Delta = 4*sum(Z^n * sum(psi(A)*A.norm()^(kj[j]-1)
    #                         for nu in nu_set[n-1]
    #                         for A in divisors(nu*Diff) if A.is_coprime(mfin))
    #               for n in [1..S-1])
    prime_dict = {}
    Delta = 0
    for n in [1..S-1]:
        # print("n=", n)
        coeff_n = 0
        # print(mess)
        for nu in nu_set[n-1]:
            factors = (nu*Diff).factor()
            # print("nu = ", nu)
            cop_factors = [[pfrak, mult] for pfrak, mult in factors if
                           pfrak.is_coprime(mfin)]
            ell = len(cop_factors)
            # print("ell = ", ell)
            # pfraks = [pfrak for pfrak, mult in cop_factors]
            psi_pfrak = []
            norm_pfrak = []
            for pfrak,_ in  cop_factors:
                # pfrak_str = f"{pfrak}"
            #     # print(pfrak_str)
            #     # print(type(pfrak_str))
            #     # print(prime_dict)
                
            #     if not pfrak_str in prime_dict:  # if we haven't seen ideal before:
            #         uppsi = psi(pfrak)  # compute psi(p) and norm(p)
            #         upnorm = pfrak.norm()
            #         a = pickle.dumps(uppsi)  # then store (string specifying) guys in shared dict
            #         b = pickle.dumps(upnorm)
            #         prime_dict[pfrak_str] = [a, b]
            #         psi_pfrak.append(uppsi)
            #         norm_pfrak.append(upnorm)
            #         # print(type(a), type(b))
            #     else:           # otherwise we unpickle and add to list
            #         # print("prime_dict:", prime_dict)
            #         ppsi, pnorm = prime_dict[pfrak_str]
            #         psi_pfrak.append(pickle.loads(ppsi))
            #         norm_pfrak.append(pickle.loads(pnorm))
            # # psi_pfrak = [psi(pfrak) for pfrak,_ in cop_factors]
            # # norm_pfrak = [pfrak.norm() for pfrak,_ in cop_factors]
                if not pfrak in prime_dict:
                    ppsi = psi(pfrak)
                    pnorm = pfrak.norm()
                    prime_dict[pfrak] = [ppsi,pnorm]
                else:
                    ppsi,pnorm = prime_dict[pfrak]
                psi_pfrak.append(ppsi)
                norm_pfrak.append(pnorm)
                    
            assert len(psi_pfrak) == len(cop_factors)
            for tup in Product(*[[0..m] for _, m in cop_factors]):
                if tup:
                    coeff_n += prod(psi_pfrak[i]**tup[i]*(norm_pfrak[i]**tup[i])**(kj[j]-1) for i in range(ell))
                else:
                    coeff_n += 1  # corr to ideal 1
        # print(f"{n}th coefficent:", coeff_n)
        Delta += 4*Z^n*coeff_n
    return_dict[j] = pickle.dumps(Delta)
    return(Delta)
            
            

def classical_L_values(F,p, k0, psi, m, verbose=False, rand_bases=True,nebentype=True):
    """
    Same as diagonal_restriction_Lp, except we return only the classical L-values.
    """
    if isinstance(psi,sage.quadratic_forms.binary_qf.BinaryQF):
        RM = True
        nebentype = False
        Q = psi
        D = Q.discriminant()
        M = sqrt(D/D.squarefree_part())
        Psi = "bqf,no psi"
    else:
    # elif isinstance(psi,sage.modular.hecke_character.HeckeCharacterGroup_class_with_category.element_class):
        RM = False
        mfrak = psi.modulus()
        M = mfrak.finite_part().smallest_integer()
        mfin = mfrak.finite_part()
        Psi = hecke_to_dirichlet_char(psi, M)
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



    if nebentype and not RM:
        Mtop = ModularForms(Psi, weight=d*kj[deltam])
    else:
        Mtop = ModularForms(Gamma1(M), weight=d*kj[deltam])
    S = Mtop.sturm_bound()
    # K = Psi.base_ring()
    K = Mtop.base_ring()


    # Mdkj_dims = [ModularForms(Psi, weight=(d*kj[j])).dimension() for j in [0..deltam]]
    if verbose:
        print("S =", S)
        print("M =", M)
        print("Weights:", [d*k for k in kj])
        print("delta_m =", deltam)
        if Psi:
            print("Psi =", Psi)
        print("Base field:", K)
        # print("dims of M_{k_j}:", Mdkj_dims, "\n")
        print("Computing spaces of modular forms: \n")
        print("Computing low weight basis.")

    if nebentype:
        LWBs =[[ModularForms(Gamma0(M), weight=l,
                             prec=S).q_expansion_basis(prec=S) for l in range(1,7)],
               [ModularForms(Psi, weight=l,prec=S).q_expansion_basis(prec=S) for l in range(1,7)]]
    else:
        LWB = [ModularForms(Gamma1(M), weight=l,
                            prec=S).q_expansion_basis(prec=S) for l in range(1,7)]
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
    # Mdkj = [ModularForms(Gamma1,
    #                      weight=(d*kj[j])
    #                      ).q_expansion_basis(prec=S)
    #         for j in [0..deltam]]

    # # print(len(Mdkj))
    # # for basis in Mdkj:
    # #     print(len(basis))


    # We compute the non-constant coefficients of the diagonal
    # restrictions, all in one go (i.e. \Delta_j for j = 0,...,delta_m)

    # if psi.is_trivial():
        # for j in [0..deltam]:
        #     Deltaj.append(2^d*sum(Z^n * sum(A.norm()^(kj[j]-1) for nu in
        #                                 tp_elts_of_trace_n(F,n) for A in
        #                                 divisors(nu*Diff)
        #                                 ) for n in [1..S-1]))
    # else:

    ##### Compute diagonal restrictions:
    # if the modulus is an integer, use bqfs
    R.<Z> = PowerSeriesRing(K)
    if verbose:
        print("\nComputing diagonal restriction data:")

    # print()
    # if RM and [g.is_rational() for g in mfin.gens()] == [True for g in mfin.gens()] or mfin == F.ideal(1):
    #     Delta_list = []
    #     D = (M*a)^2
    #     for g in F.ray_class_group(mfrak).gens():
    #         psi_g = psi(g.ideal())
    #         Q = g.ideal().quadratic_form()
    #         A, B, C = list(Q)
    #         if A != 0:
    #             tau = (-B + sqrt(D))/(2*A)
    #         else:
    #             tau = -C/B
    #         # D = Q.discriminant()
            
    #         # define a list RMForms = [(Q', gamma_n')] for n =1,...,S-1
    #         # RMforms = [get_RM_set(n, list(Q), D, f) for n in [1..S-1]]

    #         # now add this to delta_data = [(psi(tau), RMForms(tau))]

    #         if verbose:
    #             print(f"psi({g}) = {psi_g}")
    #         if psi_g != 0: # no point in computing data if psi(tau) = 0
    #             Delta_list.append([psi_g, Delta_data(Q, S)])
    #         else:
    #             Delta_list.append([psi_g, []])
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
            nu_set.append(tp_elts_of_trace_n(F, n, invDiff, floor(n*a)))
        if verbose:
            print("Computing Delta_j, IN PARALLEL!!")
            start_time = time.time()
        manager = Manager()
        # prime_dict = manager.dict()  # We try to use internal lists instead
        return_dict = manager.dict()

        job = [Process(target=compute_Delta_j, args=(j,F,psi,kj,S,Diff,nu_set,mfin,return_dict)) for j in [0..deltam]]
        _ = [p.start() for p in job]
        _ = [p.join() for p in job]
        Deltaj = [pickle.loads(return_dict[j]) for j in [0..deltam]]
        # print(compute_Delta_j(3,F,psi,kj,S,Diff,nu_set,mfin,{},{}))
        # Deltaj = sorted(list(compute_Delta_j([(j,F,psi,kj,S,Diff,nu_set,mfin) for
        #                                j in [0..deltam]]))) # returns (input,output)

        # Deltaj = [l[1] for l in Deltaj]  # retrieve output
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
    

def is_totally_odd(psi):
    """
    Test if a Hecke character psi is totally positive  # this is probably not correct
    """
    F = psi.parent().number_field()
    a = F.gen()
    mfrak = psi.parent().modulus()
    M = mfrak.finite_part().smallest_integer()
    alpha = 1-M*a
    if not alpha.is_totally_positive():
        if psi(alpha) == -1:
            return(True)
    else:
        alpha = 1+M*a
        if psi(alpha) == -1:
            return(True)
    return(False)

def MF_char_high_weight_basis(chi, k, prec, prev=[], verbose=False, LWBs=False):
    """Compute basis for spaces of modular forms M_{k}(chi)
    by multiplying lower weight basis elements.

    prev = [oldk, oldBasis = [g]]


    test:
    chi= DirichletGroup(7)[2];
    [ModularForms(chi,find_in_space(f) for f in MF_char_high_weight_basis(chi, 12, 15, verbose=True, LWBs=False)
    """
    # print(chi)
    K = chi.base_ring()
    M = chi.modulus()
    R.<q> = PowerSeriesRing(K,default_prec=prec)
    if not chi.is_trivial():
        dim = dimension_modular_forms(chi, k)
    else:
        # This seems to be fast enough
        return ModularForms(Gamma0(M), weight=k).q_expansion_basis(prec=prec)
    if verbose:
        print("Level:", M)
        print("Weight:", k)
        print("Dim:", dim)
    if dim == 0:
        return []
    if k <= 6:
        if k == 1:
            return ModularForms(Gamma1(M), weight=k).q_expansion_basis(prec=prec)
        return ModularForms(chi, weight=k).q_expansion_basis(prec=prec)
    if verbose:
        print("Computing low weight bases")
        start_time = time.time()
    # first compute lower weight bases for k =1,2,...,6:
    if not LWBs:
        LWB = [ModularForms(Gamma0(M), weight=l, prec=prec).q_expansion_basis(prec=prec)
               for l in range(1,7)]
        LWBchar = [ModularForms(chi, weight=l, prec=prec).q_expansion_basis(prec=prec)
               for l in range(1,7)]
    else:
        LWB, LWBchar = LWBs
    if verbose:
        end_time = time.time()
        print(f'Computing LWB took {end_time - start_time} seconds')
    basis = []
    if prev != []:
        # NOTE: probably faster to compute image of everything and
        # then refine to a basis
        oldk = prev[0]
        # print(oldk)
        oldBasis = prev[1]
        if verbose:
            print("Previous basis found")
            print("Computing middle weight basis:")
        MWB = MF_char_high_weight_basis(DirichletGroup(1)[0], k-oldk, prec, verbose=False);
        # we first compute all the product of every pair (g,h) with g
        # in oldBasis and h in MWB.

        old_prods = [((g*h).O(prec)).padded_list() for g in oldBasis for h in MWB]
        # next we row-reduce and remove 0-columns 
        basis = [sum(l[n]*q^n for n in [0..prec-1]) + O(q^prec)
                 for l in Matrix(K, old_prods).echelon_form().rows() if not l.is_zero()]
        
        assert (len(basis) <= dim), "too many old forms!"
        
                
                    # if verbose: 
                    #     print("Appended:", f)
                    #     print("Now len(basis) = ", len(basis))
        

    if verbose:
        print("Computing random products")
        start_time = time.time()

    while len(basis) < dim:
        P = Partitions(k, min_part=1, max_part=6, max_length=ceil(k/6+1)
                       ).random_element()
        lenP = len(P)
        for i in range(lenP):
            fchars = LWBchar[P[i]-1]

            if not fchars == []:        
                # trivial nebentypus for all but the last
                
                f = prod(sum(QQ.random_element()*g for g in LWB[P[j]-1])
                         for j in range(lenP) if j != i)
                # note: QQ.random_element() doesn't seem very random
                f *= QQ.random_element()*sample(fchars,1)[0]
                f = R(f).lift_to_precision(prec)
                f = f.O(prec)

                if f.padded_list() != [0 for l in range(0,f.prec())]:
                    if not is_qexp_in_space(basis, f):
                        basis.append(f)
            # if verbose: 
            #     print("Appended:", f)
            #     print("Now len(basis) = ", len(basis))
    if verbose:
        end_time = time.time()
        print(f'Computing random products took {end_time - start_time} seconds')
    # print("done")
    return(basis)

def diag_restr_coeffs(Q, ks, m):
    """Compute Delta_j corresponding to the L-function L(1_[F]-1_[F
    conjugate], 1-k) for k in ks, to precision m
    Uses QF-stuff from RM-set.spyx

    """
    As, Forms = Delta_data(Q, m)
    # print("computed QF data")
    D = Q.discriminant()

    D0 = fundamental_discriminant(D)
    f = sqrt(D/D0)
    print(f)
    R.<Z> = PowerSeriesRing(ZZ, m)
    Diag_Fs = vector(R, [0 for k in ks])

    for n in [1..m-1]:
        Ap = As[n-1][0]
        Am = As[n-1][1]
        # print("n=",n)
        # print("Ap:", Ap)
        # print("Am:", Am)
        # print()
        coeff_ns = vector(R, [0 for k in ks])
        for i in range(len(Ap)):
            if gcd(Ap[i], f) == 1:
                coeff_ns += vector(R, [Ap[i]**(k-1) for k in ks])
            if gcd(Am[i], f) == 1:
                coeff_ns += vector(R, [(-1)**k*(-Am[i])**(k-1) for k in ks])
                
        Diag_Fs += 4*coeff_ns*Z^n

    return list(Diag_Fs)
   
if __name__ == "__main__":
    F.<a> = NumberField(x^2-5)
    p = 3
    k0 = 3
    m = 12
    r1 = F.signature()[0]       # nr of real embeddings = 2
    mfrak = F.modulus(F.ideal(4), range(r1))  # narrow class group
    H = HeckeCharacterGroup(mfrak)
    try:
        psi = H.gens()[1]
        # psi = H.one()
        # psi = psi.primitive_character()
    except:
        psi = H.one()
    # F.dirichlet_group()[0]
    print("F = ", F,
          "\n H = ", H.group(),
          "\np = ", p,
          "\nk0 = ", k0,
          "\nm = ", m,
          "\npsi = ", psi)
    # P = diagonal_restriction_Lp(F, p, k0, psi, m, verbose=True)
    # print("P(s) = ", P)
    # print("Q(T) = ", Q)

#     with open('pAdicLvalues.csv', newline='') as csvfile:
